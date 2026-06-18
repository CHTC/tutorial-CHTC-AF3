#!/bin/bash

#set -x #for complete debugging

# NOTE: intentionally NOT using `set -e`. Contributing MSAs back to the cache
# is a best-effort, post-success side effect and must NEVER fail the AF3 job
# (the alignment output is already produced). Every failure path here degrades
# gracefully; the script always exits 0 except for invalid CLI arguments
# (exit 2).
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
AF_INPUT_DIR=""    # post-cache-lookup AF3 inputs (af_input)
AF_OUTPUT_DIR=""   # AF3 data-pipeline output dir (af_output)
API_KEY=""
API_BASE_URL="https://149.165.170.71.sslip.io"
DB_VERSION="AlphaFold Monomer v2.0 pipeline 2025-08-01T00:00:00Z"
REQUESTER=""
WORK_DIR="."
MAX_ATTEMPTS=4      # uploads/request: 1 initial attempt + 3 retries (transient)
UPLOAD_RETRIES=5    # pelican object put: retries per file AFTER the first try
CURL_MAX_TIME=120   # seconds, per request
PYBIN="${PYTHON_BIN:-python3}"

# ---- argument parsing ----
while [[ $# -gt 0 ]]; do
  case $1 in
    --af_input_dir)
      AF_INPUT_DIR="$2"; shift; shift;;
    --af_output_dir)
      AF_OUTPUT_DIR="$2"; shift; shift;;
    --api_key)
      API_KEY="$2"; shift; shift;;
    --api_base_url)
      API_BASE_URL="$2"; shift; shift;;
    --db_version)
      DB_VERSION="$2"; shift; shift;;
    --requester)
      REQUESTER="$2"; shift; shift;;
    --work_dir)
      WORK_DIR="$2"; shift; shift;;
    --upload_retries)
      UPLOAD_RETRIES="$2"; shift; shift;;
    -v|--verbose)
      VERBOSE_LEVEL=2; shift;;
    -s|--silent)
      VERBOSE_LEVEL=0; shift;;
    -*|--*)
      printerr "Unknown option $1 — skipping MSA cache upload"; exit 0;;
    *)
      printerr "Unexpected argument $1 — skipping MSA cache upload"; exit 0;;
  esac
done

# ---- validate arguments. The cache must NEVER interrupt the job, so any bad or
# ---- missing argument is logged and we exit 0 (skip upload, let the job go on).
if [[ -z "$AF_INPUT_DIR" ]]; then printerr "--af_input_dir is required — skipping MSA cache upload"; exit 0; fi
if [[ -z "$AF_OUTPUT_DIR" ]]; then printerr "--af_output_dir is required — skipping MSA cache upload"; exit 0; fi
if [[ -z "$API_KEY" ]]; then printerr "--api_key is required — skipping MSA cache upload"; exit 0; fi
if [[ -z "$DB_VERSION" ]]; then printerr "--db_version must be non-empty — skipping MSA cache upload"; exit 0; fi
[[ "$UPLOAD_RETRIES" =~ ^[0-9]+$ ]] || UPLOAD_RETRIES=5
if ! command -v "$PYBIN" >/dev/null 2>&1; then
  printerr "$PYBIN not found — skipping MSA cache upload"
  exit 0
fi
if ! command -v curl >/dev/null 2>&1; then
  printerr "curl not found — skipping MSA cache upload"
  exit 0
fi
if [[ ! -d "$AF_OUTPUT_DIR" ]]; then
  printinfo "AF3 output dir not found (${AF_OUTPUT_DIR}) — nothing to contribute"
  exit 0
fi
mkdir -p "$WORK_DIR"

PLAN="$(mktemp)"
trap 'rm -f "${PLAN:-}" "${WORK_DIR:-/nonexistent}"/.token_*.jwt 2>/dev/null || true' EXIT

# ---------------------------------------------------------------------------
# Step 1: decide what to upload. A protein chain is "recalculated de-novo" iff
# it had NO unpairedMsa in the post-cache-lookup input (neither user-provided
# nor cache-injected). For each such chain, pull the computed a3m out of the
# AF3 output *_data.json, write it to a file, and record sha256 + size.
# Normalization/hashing is byte-identical to cache_lookup.sh and the server
# (af3-cache/app/core/hashing.py): "".join(seq.split()).upper().
# Emits TSV: <seq_hash>\t<raw_sequence>\t<a3m_path>\t<sha256>\t<size_bytes>
# ---------------------------------------------------------------------------
if ! "$PYBIN" - "$AF_INPUT_DIR" "$AF_OUTPUT_DIR" "$WORK_DIR" > "$PLAN" 2>/dev/null <<'PY'
import glob, hashlib, json, os, sys

af_input_dir, af_output_dir, work_dir = sys.argv[1], sys.argv[2], sys.argv[3]


def norm_hash(seq):
    n = "".join(seq.split()).upper()
    return n, hashlib.sha256(n.encode()).hexdigest()


def proteins(path):
    try:
        with open(path) as fh:
            data = json.load(fh)
    except Exception:
        return
    for entry in data.get("sequences") or []:
        if isinstance(entry, dict) and isinstance(entry.get("protein"), dict):
            yield entry["protein"]


# Chains that already had an MSA going in => NOT recalculated => skip.
had_msa = set()
for jf in sorted(glob.glob(os.path.join(af_input_dir, "*.json"))):
    for prot in proteins(jf):
        raw = prot.get("sequence")
        msa = prot.get("unpairedMsa")
        if isinstance(raw, str) and raw.strip() and isinstance(msa, str) and msa.strip():
            had_msa.add(norm_hash(raw)[1])

emitted = set()
for df in sorted(glob.glob(os.path.join(af_output_dir, "**", "*_data.json"), recursive=True)):
    for prot in proteins(df):
        raw = prot.get("sequence")
        msa = prot.get("unpairedMsa")
        if not (isinstance(raw, str) and raw.strip()):
            continue
        if not (isinstance(msa, str) and msa.strip()):
            continue
        norm, h = norm_hash(raw)
        if h in had_msa or h in emitted:
            continue
        emitted.add(h)
        a3m_path = os.path.join(work_dir, "upload_%s.a3m" % h)
        b = msa.encode()
        with open(a3m_path, "wb") as out:
            out.write(b)
        sha = hashlib.sha256(b).hexdigest()
        sys.stdout.write("\t".join([h, norm, a3m_path, sha, str(len(b))]) + "\n")
PY
then
  printerr "Could not scan AF3 inputs/outputs for upload — skipping MSA cache contribution"
  exit 0
fi

if [[ ! -s "$PLAN" ]]; then
  printinfo "No de-novo protein alignments to contribute to the cache"
  exit 0
fi

n_total=$(wc -l < "$PLAN" | tr -d ' ')
printinfo "Found ${n_total} de-novo protein alignment(s) to contribute to the cache"

# ---------------------------------------------------------------------------
# Helper: POST a JSON body, capturing the body to $1 and printing the HTTP code.
# ---------------------------------------------------------------------------
post_json() { # $1=outfile $2=url $3=json-body
  curl -k -sS -X 'POST' "$2" \
    -H 'accept: application/json' \
    -H "X-API-Key: $API_KEY" \
    -H 'Content-Type: application/json' \
    -d "$3" \
    -o "$1" \
    -w '%{http_code}' \
    --max-time "$CURL_MAX_TIME" 2>/dev/null
}

n_uploaded=0
n_already=0
n_skipped=0

while IFS=$'\t' read -r seq_hash norm a3m_path sha size; do
  [[ -z "$seq_hash" ]] && continue
  short="${seq_hash:0:12}…"

  req_body="$("$PYBIN" -c '
import json, sys
seq, dbv, sha, size, req = sys.argv[1:6]
body = {"sequence": seq, "db_version": dbv,
        "declared_size_bytes": int(size), "declared_sha256": sha}
if req:
    body["requester"] = req
print(json.dumps(body))
' "$norm" "$DB_VERSION" "$sha" "$size" "$REQUESTER")"

  # ---- request a scoped write token (retry transient failures only) ----
  resp="$(mktemp)"
  http=""
  attempt=1
  while [[ $attempt -le $MAX_ATTEMPTS ]]; do
    http="$(post_json "$resp" "${API_BASE_URL}/v1/uploads/request" "$req_body")"
    rc=$?
    if [[ $rc -eq 0 ]]; then
      case "$http" in
        200|409|413|422|400) break;;  # terminal — no retry
      esac
    fi
    printerr "uploads/request attempt ${attempt}/${MAX_ATTEMPTS} failed (curl rc=${rc}, http=${http:-none}) for ${short}"
    if [[ $attempt -lt $MAX_ATTEMPTS ]]; then
      s=$(( RANDOM % 46 + 15 ))
      printinfo "Retrying uploads/request for ${short} in ${s}s"
      sleep "$s"
    fi
    attempt=$(( attempt + 1 ))
  done

  if [[ "$http" == "409" ]]; then
    existing="$("$PYBIN" -c 'import json,sys
try: print((json.load(open(sys.argv[1])).get("error") or {}).get("osdf_uri") or "")
except Exception: print("")' "$resp" 2>/dev/null)"
    printinfo "Already cached for ${short} (existing: ${existing:-?}) — nothing to upload"
    n_already=$(( n_already + 1 ))
    rm -f "$resp" "$a3m_path"
    continue
  fi
  if [[ "$http" == "413" ]]; then
    printerr "a3m for ${short} exceeds the cache size limit — skipping upload"
    n_skipped=$(( n_skipped + 1 ))
    rm -f "$resp" "$a3m_path"
    continue
  fi
  if [[ "$http" != "200" ]]; then
    printerr "uploads/request returned http=${http:-none} for ${short} — skipping upload"
    n_skipped=$(( n_skipped + 1 ))
    rm -f "$resp" "$a3m_path"
    continue
  fi

  # ---- parse the token + target out of the 200 response ----
  parsed="$("$PYBIN" -c 'import json,sys
d=json.load(open(sys.argv[1]))
print(d["ingestion_id"]); print(d["osdf_uri"]); print(d["write_token"])' "$resp" 2>/dev/null)"
  pp_rc=$?
  rm -f "$resp"
  if [[ $pp_rc -ne 0 ]]; then
    printerr "Could not parse uploads/request response for ${short} — skipping upload"
    n_skipped=$(( n_skipped + 1 ))
    rm -f "$a3m_path"
    continue
  fi
  ingestion_id="$(sed -n '1p' <<<"$parsed")"
  osdf_uri="$(sed -n '2p' <<<"$parsed")"
  write_token="$(sed -n '3p' <<<"$parsed")"

  if ! command -v pelican >/dev/null 2>&1; then
    printerr "pelican not found — cannot upload a3m for ${short}; skipping (not calling complete)"
    n_skipped=$(( n_skipped + 1 ))
    rm -f "$a3m_path"
    continue
  fi

  # Each upload-request returns its OWN write token, scoped to this object's
  # path (tokens are NOT interchangeable between files). Save it to a private
  # per-object file and use only that token for this object's put.
  token_file="${WORK_DIR}/.token_${ingestion_id}.jwt"
  ( umask 077; printf '%s' "$write_token" > "$token_file" )

  # ---- pelican object put, retried up to UPLOAD_RETRIES times after the first
  # ---- attempt (15-60s back-off), reusing this object's scoped token (valid
  # ---- ~1200s). Exhausting all attempts is logged but NOT fatal — contributing
  # ---- to the cache is best-effort and must never interrupt the job.
  put_ok=0
  put_attempt=1
  put_total=$(( UPLOAD_RETRIES + 1 ))
  put_log="$(mktemp)"
  while [[ $put_attempt -le $put_total ]]; do
    printverbose "pelican object put (attempt ${put_attempt}/${put_total}) ${a3m_path} -> ${osdf_uri}"
    if pelican object put "$a3m_path" "$osdf_uri" --token "$token_file" > "$put_log" 2>&1; then
      put_ok=1
      break
    fi
    printerr "pelican upload attempt ${put_attempt}/${put_total} failed for ${short}"
    printverbose "pelican output: $(cat "$put_log" 2>/dev/null)"
    if [[ $put_attempt -lt $put_total ]]; then
      s=$(( RANDOM % 46 + 15 ))
      printinfo "Retrying upload for ${short} in ${s}s"
      sleep "$s"
    fi
    put_attempt=$(( put_attempt + 1 ))
  done
  rm -f "$put_log" "$token_file"

  if [[ $put_ok -ne 1 ]]; then
    printerr "All ${put_total} upload attempts failed for ${short} — NOT calling complete; continuing (cache contribution is best-effort)"
    n_skipped=$(( n_skipped + 1 ))
    rm -f "$a3m_path"
    continue
  fi

  # ---- only after a successful put, tell the server to verify + publish ----
  cresp="$(mktemp)"
  chttp="$(curl -k -sS -X 'POST' "${API_BASE_URL}/v1/uploads/${ingestion_id}/complete" \
    -H 'accept: application/json' \
    -H "X-API-Key: $API_KEY" \
    -o "$cresp" -w '%{http_code}' --max-time "$CURL_MAX_TIME" 2>/dev/null)"
  if [[ "$chttp" == "200" ]]; then
    status="$("$PYBIN" -c 'import json,sys
d=json.load(open(sys.argv[1])); print(d.get("status","?"), d.get("reject_reason") or "")' "$cresp" 2>/dev/null)"
    printinfo "Contributed a3m for ${short} — server status: ${status}"
    n_uploaded=$(( n_uploaded + 1 ))
  else
    printerr "uploads/${ingestion_id}/complete returned http=${chttp:-none} for ${short} (object uploaded; publish unconfirmed)"
    n_skipped=$(( n_skipped + 1 ))
  fi
  rm -f "$cresp" "$a3m_path"
done < "$PLAN"

printinfo "MSA cache contribution summary: ${n_uploaded} uploaded, ${n_already} already-cached, ${n_skipped} skipped"
exit 0
