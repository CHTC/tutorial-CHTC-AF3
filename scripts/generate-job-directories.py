#!/usr/bin/env python3
"""
Generate AlphaFold3 job directories from a manifest CSV.

Manifest structure:
    job_name, mol1_type, mol1_chain, mol1_seq, mol2_type, mol2_chain, mol2_seq, ...

Where molecule types can be:
    - protein
    - rna
    - dna
    (extendable)

Each output job directory is placed under:
    output_dir/Job#_job_name/
with a companion data_inputs/fold_input.json and empty inference_inputs/ directory created beside it.

If multiple chains are desired for a molecule, use "|" to separate chain IDs in the molN_chain field:
    e.g., "A|B|C" for chains A, B, and C

Usage:
    python3 generate-job-directories.py --manifest manifest.csv --output_dir AF3_Jobs
"""

import os
import csv
import json
import argparse


# Build a minimal AF3 JSON block for a single molecule.
# Supports:
#   - standard protein, RNA, or DNA sequences
#   - single chain ("A") or multichain ("A|C|D") ➜ ["A", "C", "D"]
def build_molecule_block(molecule_type, chain_id, sequence, apply_mods=True):
    sequence = sequence.strip()

    # Convert chains: "A|B|C" → ["A", "B", "C"]
    if "|" in chain_id:
        chains = [c.strip() for c in chain_id.split("|") if c.strip()]
    else:
        chains = chain_id.strip()  # single-chain mode

    return {molecule_type.strip(): {"id": chains, "sequence": sequence}}


def parse_molecules(row_dict):
    """
    Find all molecule triplets in the row:
        molN_type, molN_chain, molN_seq
    Returns a list of molecule blocks.
    """

    molecules = []
    idx = 1

    while True:
        type_key = f"mol{idx}_type"
        chain_key = f"mol{idx}_chain"
        seq_key = f"mol{idx}_seq"

        # Stop when we no longer find molN_type
        if type_key not in row_dict:
            break

        mol_type_val = row_dict.get(type_key)
        mol_chain_val = row_dict.get(chain_key)
        mol_seq_val = row_dict.get(seq_key)

        # Allow partially empty trailing triplets; stop cleanly
        if mol_type_val is None and mol_chain_val is None and mol_seq_val is None:
            break

        # Skip this molecule if any required field is missing or empty
        if not mol_type_val or not mol_chain_val or not mol_seq_val:
            idx += 1
            continue

        mol_type = mol_type_val.strip()
        mol_chain = mol_chain_val.strip()
        mol_seq = mol_seq_val.strip()

        if mol_type and mol_chain and mol_seq:
            molecules.append(build_molecule_block(mol_type, mol_chain, mol_seq))

        idx += 1

    return molecules


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""\
        Generate AlphaFold3 job directories from a manifest CSV.

        Manifest structure:
            job_name, mol1_type, mol1_chain, mol1_seq, mol2_type, mol2_chain, mol2_seq, ...

        Where molecule types can be:
            - protein
            - rna
            - dna
            (extendable)

        Each output job directory is placed under:
            output_dir/Job#_job_name/
        with a companion data_inputs/fold_input.json and empty inference_inputs/ directory created beside it.

        If multiple chains are desired for a molecule, use "|" to separate chain IDs in the molN_chain field:
            e.g., "A|B|C" for chains A, B, and C

        Usage:
            python generate-job-directories.py [options] --manifest MANIFEST --output_dir OUTPUT_DIR
        """,
    )
    parser.add_argument("--manifest", required=True, help="Path to manifest CSV")
    parser.add_argument(
        "--output_dir", required=True, help="Where to create job directories"
    )
    parser.add_argument(
        "--jobs_list",
        default="./list_of_af3_jobs.txt",
        help="Path to write job directory list (default: ./list_of_af3_jobs.txt)",
    )
    parser.add_argument(
        "--seed",
        nargs="+",
        default=[1],
        type=int,
        help="Seed(s) to use for reproducible runs. Can specify multiple seeds for multiple runs per job (default: 1). Accepts space-separated values (e.g. --seed 1 2 3)",
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    jobs_list_path = args.jobs_list
    jobs_list_file = open(jobs_list_path, "w")

    with open(args.manifest, newline="") as csvfile:
        reader = csv.DictReader(csvfile)

        if "job_name" not in reader.fieldnames:
            raise ValueError("Manifest must include a 'job_name' column.")

        for idx, row in enumerate(reader, start=1):
            job_name = row["job_name"].strip()
            job_dir_name = f"Job{idx}_{job_name}"
            job_dir = os.path.join(args.output_dir, job_dir_name)

            data_inputs = os.path.join(job_dir, "data_inputs")
            inference_inputs = os.path.join(job_dir, "inference_inputs")

            os.makedirs(data_inputs, exist_ok=True)
            os.makedirs(inference_inputs, exist_ok=True)

            # Extract molecules (protein/RNA/etc.)
            molecules = parse_molecules(row)

            fold_json = {
                "name": job_name,
                "sequences": molecules,
                "modelSeeds": args.seed,
                "dialect": "alphafold3",
                "version": 1,
            }

            json_path = os.path.join(data_inputs, "fold_input.json")
            with open(json_path, "w") as jf:
                json.dump(fold_json, jf, indent=2)

            # Compute total tokens (sequence length across all chains)
            total_tokens = 0
            for mol in molecules:
                mol_type = list(mol.keys())[0]
                seq = mol[mol_type]["sequence"]
                chains = mol[mol_type]["id"]

                if isinstance(chains, list):
                    total_tokens += len(seq) * len(chains)
                else:
                    total_tokens += len(seq)

            print(f"[+] Created {job_dir_name} with {len(molecules)} molecules and {total_tokens} tokens.")
            jobs_list_file.write(f"{job_dir_name}\n")

    jobs_list_file.close()
    print("\nAll multi-molecule AF3 job directories generated.")
    print("\nThe list of job directories has been written to:", jobs_list_path)


if __name__ == "__main__":
    main()
