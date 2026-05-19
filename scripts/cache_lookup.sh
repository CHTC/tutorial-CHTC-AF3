#!/bin/bash

#set -x #for complete debugging

# NOTE: intentionally NOT using `set -e`. The MSA cache is a best-effort
# optimization and must NEVER fail the AF3 job: every failure path here
# (network, non-200, parse error, pelican download, sha256 mismatch) is
# handled explicitly and degrades to de-novo alignment. The script always
# exits 0 except for invalid CLI arguments (exit 2).
set -uo pipefail

# ---- logging (same convention as scripts/data_pipeline.sh) ----
VERBOSE_LEVEL=1 # 0 = silent, 1 = info, 2 = verbose
function printstd() { echo "$@"; }
function printerr() { echo "ERROR: $@" 1>&2; }
function printinfo() {
  if [[ $VERBOSE_LEVEL -ge 1 ]]; then
    printstd "INFO: $@"
  fi
}
function printverbose() {
  if [[ $VERBOSE_LEVEL -ge 2 ]]; then
    printstd "DEBUG: $@"
  fi
}

# ---- defaults ----
INPUT_JSON=""
API_KEY="a66d4cd10c1cadc6cc9a6ae6e8eb4fe2f7bdc924fd1904a535dd138282dedd16"
API_URL="https://149.165.170.71.sslip.io/v1/query"
PREFERRED_SOURCES="OSG-Generated, GDM, Nvidia, Viro3D, BFVD, AllTheBacteria, Kinetoplastid, Community Contributed - Cached"
WORK_DIR="."
MAX_ATTEMPTS=4    # 1 initial attempt + 3 retries
CURL_MAX_TIME=60  # seconds, per request
PYBIN="${PYTHON_BIN:-python3}"

# ---- argument parsing ----
while [[ $# -gt 0 ]]; do
  case $1 in
    --input_json)
      INPUT_JSON="$2"; shift; shift;;
    --api_key)
      API_KEY="$2"; shift; shift;;
    --api_url)
      API_URL="$2"; shift; shift;;
    --preferred_sources)
      PREFERRED_SOURCES="$2"; shift; shift;;
    --work_dir)
      WORK_DIR="$2"; shift; shift;;
    -v|--verbose)
      VERBOSE_LEVEL=2; shift;;
    -s|--silent)
      VERBOSE_LEVEL=0; shift;;
    -*|--*)
      printerr "Unknown option $1"; exit 2;;
    *)
      printerr "Unexpected argument $1"; exit 2;;
  esac
done

# ---- validate arguments (bad args => exit 2; anything else => exit 0) ----
if [[ -z "$INPUT_JSON" ]]; then printerr "--input_json is required"; exit 2; fi
if [[ ! -f "$INPUT_JSON" ]]; then printerr "--input_json not found: $INPUT_JSON"; exit 2; fi
if [[ -z "$API_KEY" ]]; then printerr "--api_key is required"; exit 2; fi
if ! command -v "$PYBIN" >/dev/null 2>&1; then
  printerr "$PYBIN not found — skipping MSA cache lookup, proceeding de-novo"
  exit 0
fi
mkdir -p "$WORK_DIR"

SEQLIST="$(mktemp)"
MAPFILE="$(mktemp)"
PELICAN_LOG="$(mktemp)"
trap 'rm -f "${SEQLIST:-}" "${MAPFILE:-}" "${PELICAN_LOG:-}" 2>/dev/null || true' EXIT

# ---------------------------------------------------------------------------
# Step 1: extract unique protein sequences (normalized + hashed) that still
# need an MSA. Normalization MUST stay byte-identical to the cache server
# (af3-cache/app/core/hashing.py): "".join(seq.split()).upper().
# ---------------------------------------------------------------------------
if ! "$PYBIN" - "$INPUT_JSON" > "$SEQLIST" 2>/dev/null <<'PY'
import json, sys, hashlib

with open(sys.argv[1]) as fh:
    data = json.load(fh)

seen = set()
for entry in data.get("sequences") or []:
    if not isinstance(entry, dict):
        continue
    prot = entry.get("protein")
    if not isinstance(prot, dict):
        continue
    raw = prot.get("sequence")
    if not raw or not isinstance(raw, str):
        continue
    existing = prot.get("unpairedMsa")
    if isinstance(existing, str) and existing.strip():
        # respect a user-provided MSA — never overwrite it
        continue
    norm = "".join(raw.split()).upper()
    if not norm:
        continue
    h = hashlib.sha256(norm.encode()).hexdigest()
    if h in seen:
        continue
    seen.add(h)
    sys.stdout.write(h + "\t" + norm + "\n")
PY
then
  printerr "Could not read protein sequences from ${INPUT_JSON} — skipping cache lookup, proceeding de-novo"
  exit 0
fi

if [[ ! -s "$SEQLIST" ]]; then
  printinfo "No protein sequences needing an MSA in $(basename "$INPUT_JSON") — nothing to look up"
  exit 0
fi

# ---------------------------------------------------------------------------
# Step 2: per unique protein sequence — query the cache, pick a preferred-
# source hit, download + verify the .a3m. Successful hits are recorded in
# MAPFILE as "<seq_hash>\t<a3m_path>" for the single atomic JSON rewrite.
# ---------------------------------------------------------------------------
while IFS=$'\t' read -r seq_hash norm; do
  [[ -z "$seq_hash" ]] && continue
  short="${seq_hash:0:12}…"

  body="$("$PYBIN" -c 'import json,sys; print(json.dumps({"sequence": sys.argv[1]}))' "$norm")"
  resp_file="$(mktemp)"

  ok=0
  attempt=1
  while [[ $attempt -le $MAX_ATTEMPTS ]]; do
    http_code="$(curl -k -sS -X 'POST' "$API_URL" \
      -H 'accept: application/json' \
      -H "X-API-Key: $API_KEY" \
      -H 'Content-Type: application/json' \
      -d "$body" \
      -o "$resp_file" \
      -w '%{http_code}' \
      --max-time "$CURL_MAX_TIME" 2>/dev/null)"
    curl_rc=$?
    if [[ $curl_rc -eq 0 && "$http_code" == "200" ]]; then
      ok=1
      break
    fi
    printerr "Cache query attempt ${attempt}/${MAX_ATTEMPTS} failed (curl rc=${curl_rc}, http=${http_code:-none}) for ${short}"
    if [[ $attempt -lt $MAX_ATTEMPTS ]]; then
      s=$(( RANDOM % 46 + 15 ))   # 15–60s random backoff
      printinfo "Retrying cache query for ${short} in ${s}s"
      sleep "$s"
    fi
    attempt=$(( attempt + 1 ))
  done

  if [[ $ok -ne 1 ]]; then
    printerr "MSA cache not available after ${MAX_ATTEMPTS} attempts — proceeding with de-novo alignment for ${short}"
    rm -f "$resp_file"
    continue
  fi

  # Parse the response and select a hit by preferred-source order.
  # The selector also writes a human-readable hit summary to stderr; bash
  # surfaces it via printverbose (only visible at -v / VERBOSE_LEVEL >= 2).
  selstderr="$(mktemp)"
  sel="$("$PYBIN" - "$resp_file" "$PREFERRED_SOURCES" "$seq_hash" 2>"$selstderr" <<'PY'
import json, sys

resp_path, preferred_csv, exp_hash = sys.argv[1], sys.argv[2], sys.argv[3]
try:
    with open(resp_path) as fh:
        data = json.load(fh)
except Exception:
    sys.exit(1)

results = data.get("results") or []
res = next((r for r in results if r.get("seq_hash") == exp_hash), None)
if res is None and results:
    res = results[0]

hits = []
if res and res.get("found"):
    hits = [h for h in (res.get("hits") or []) if isinstance(h, dict) and h.get("osdf_uri")]

preferred = [s.strip() for s in preferred_csv.split(",") if s.strip()]
chosen = None
nopref = False
if hits:
    if preferred:
        best = None
        for h in hits:
            src = h.get("source")
            if src in preferred:
                rank = preferred.index(src)
                if best is None or rank < best:
                    best, chosen = rank, h
        if chosen is None:
            nopref = True
    else:
        # No preferred-source list configured -> accept the first hit.
        chosen = hits[0]

# Verbose hit summary on stderr (caller decides whether to display).
if not hits:
    sys.stderr.write("no hits returned\n")
else:
    sys.stderr.write("%d hit(s) returned:\n" % len(hits))
    for h in hits:
        mark = " [chosen]" if h is chosen else ""
        sys.stderr.write("  - source=%r db_version=%s osdf_uri=%s%s\n" % (
            h.get("source"), h.get("db_version"), h.get("osdf_uri"), mark))

if not hits:
    print("MISS")
elif nopref:
    avail = ", ".join(sorted({str(h.get("source")) for h in hits}))
    print("NOPREF\t" + avail)
else:
    print("HIT\t%s\t%s" % (chosen.get("osdf_uri"), chosen.get("a3m_sha256") or ""))
PY
)"
  sel_rc=$?
  # Surface the selector's per-query hit summary at verbose level.
  while IFS= read -r _l; do
    printverbose "  ${_l}"
  done < "$selstderr"
  rm -f "$selstderr" "$resp_file"
  if [[ $sel_rc -ne 0 ]]; then
    printerr "Could not parse cache response for ${short} — proceeding de-novo"
    continue
  fi

  kind="${sel%%$'\t'*}"
  case "$kind" in
    MISS)
      printinfo "No cache hit for ${short} — computing MSA de-novo"
      continue;;
    NOPREF)
      avail="${sel#*$'\t'}"
      printinfo "Cache hit(s) for ${short} but none from preferred sources (have: ${avail}) — computing MSA de-novo"
      continue;;
    HIT)
      :;;
    *)
      printerr "Unexpected cache-selection output for ${short} — proceeding de-novo"
      continue;;
  esac

  rest="${sel#HIT$'\t'}"
  osdf_uri="${rest%%$'\t'*}"
  exp_sha="${rest#*$'\t'}"
  [[ "$exp_sha" == "$rest" && "$exp_sha" == "$osdf_uri" ]] && exp_sha=""

  a3m_path="${WORK_DIR}/${seq_hash}.a3m"
  printverbose "pelican object get ${osdf_uri} -> ${a3m_path}"
  if ! pelican object get "$osdf_uri" "$a3m_path" > "$PELICAN_LOG" 2>&1; then
    printerr "pelican download failed for ${short} (${osdf_uri}) — proceeding de-novo"
    printverbose "pelican output: $(cat "$PELICAN_LOG" 2>/dev/null)"
    rm -f "$a3m_path"
    continue
  fi

  if [[ -n "$exp_sha" ]]; then
    got_sha="$("$PYBIN" -c 'import hashlib,sys; print(hashlib.sha256(open(sys.argv[1],"rb").read()).hexdigest())' "$a3m_path" 2>/dev/null)" || got_sha=""
    if [[ "$got_sha" != "$exp_sha" ]]; then
      printerr "sha256 mismatch for ${short} (expected ${exp_sha:0:12}…, got ${got_sha:0:12}…) — discarding cached MSA, proceeding de-novo"
      rm -f "$a3m_path"
      continue
    fi
  fi

  printf '%s\t%s\n' "$seq_hash" "$a3m_path" >> "$MAPFILE"
  printinfo "Cache HIT for ${short} — downloaded cached unpaired MSA"
done < "$SEQLIST"

# ---------------------------------------------------------------------------
# Step 3: single atomic rewrite of the input JSON — for every protein chain
# with a cache hit set unpairedMsa to the cached .a3m and pairedMsa to null
# (AF3 still computes the paired MSA and templates itself).
# ---------------------------------------------------------------------------
if [[ -s "$MAPFILE" ]]; then
  if n="$("$PYBIN" - "$INPUT_JSON" "$MAPFILE" 2>/dev/null <<'PY'
import json, sys, hashlib, os, tempfile

json_path, map_path = sys.argv[1], sys.argv[2]

mapping = {}
with open(map_path) as fh:
    for line in fh:
        line = line.rstrip("\n")
        if not line:
            continue
        h, _, path = line.partition("\t")
        if h and path:
            mapping[h] = path

with open(json_path) as fh:
    data = json.load(fh)

count = 0
for entry in data.get("sequences") or []:
    if not isinstance(entry, dict):
        continue
    prot = entry.get("protein")
    if not isinstance(prot, dict):
        continue
    raw = prot.get("sequence")
    if not raw or not isinstance(raw, str):
        continue
    existing = prot.get("unpairedMsa")
    if isinstance(existing, str) and existing.strip():
        continue
    norm = "".join(raw.split()).upper()
    if not norm:
        continue
    h = hashlib.sha256(norm.encode()).hexdigest()
    a3m_path = mapping.get(h)
    if not a3m_path:
        continue
    with open(a3m_path) as af:
        prot["unpairedMsa"] = af.read()
    prot["pairedMsa"] = None
    count += 1

d = os.path.dirname(os.path.abspath(json_path)) or "."
fd, tmp = tempfile.mkstemp(prefix=".cache_lookup.", dir=d)
try:
    with os.fdopen(fd, "w") as out:
        json.dump(data, out, indent=2)
        out.write("\n")
    os.replace(tmp, json_path)
except Exception:
    try:
        os.unlink(tmp)
    except OSError:
        pass
    raise

print(count)
PY
)"; then
    printinfo "Injected cached unpaired MSA into ${n} protein chain(s) in $(basename "$INPUT_JSON")"
  else
    printerr "Failed to inject cached MSAs into ${INPUT_JSON} — leaving input unchanged, proceeding de-novo"
  fi
else
  printinfo "No usable cache hits for $(basename "$INPUT_JSON") — all protein chains will be aligned de-novo"
fi

exit 0
