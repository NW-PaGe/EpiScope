#!/usr/bin/env bash
set -euo pipefail

# Sketch cluster_references/ into a sourmash reference DB (the sourmash
# counterpart of make_ref_db.sh).
#
# Output: $OUTDIR/ref_sourmash_db.zip   (OUTDIR defaults to the current dir)
#
# Uses `sourmash scripts manysketch` so every sketch is named after its file
# (e.g. "refA.fa"). Plain `sourmash sketch dna` leaves sketches unnamed, which
# makes the results unusable downstream.
#
# Env overrides: OUTDIR, THREADS (default 4), KSIZE (default 21),
#                SCALED (default 1000; lower = more sensitive for very close
#                genomes, but slower and larger)
#
# Requirements: sourmash, sourmash_plugin_branchwater, python3

OUTDIR="${OUTDIR:-.}"
THREADS="${THREADS:-4}"
KSIZE="${KSIZE:-21}"
SCALED="${SCALED:-1000}"
REF_DIR="cluster_references"
DB="$OUTDIR/ref_sourmash_db.zip"

if ! sourmash scripts manysketch --help >/dev/null 2>&1; then
  echo "Error: 'sourmash scripts manysketch' not found; install sourmash_plugin_branchwater." >&2
  exit 1
fi
[ -d "$REF_DIR" ] || { echo "Error: $REF_DIR/ not found." >&2; exit 1; }

mkdir -p "$OUTDIR"
TMP_CSV=$(mktemp)
trap 'rm -f "$TMP_CSV"' EXIT

# Build the manysketch input: name = file basename, FASTA files only.
python3 - "$REF_DIR" "$TMP_CSV" <<'PYEOF'
import csv, sys
from pathlib import Path

ref_dir, out = Path(sys.argv[1]), sys.argv[2]
exts = (".fa", ".fasta", ".fna", ".fas")

def is_fasta(p):
    n = p.name.lower()
    if n.endswith(".gz"):
        n = n[:-3]
    return n.endswith(exts)

files = sorted(p for p in ref_dir.iterdir() if p.is_file() and is_fasta(p))
if not files:
    sys.exit(f"Error: no FASTA files found in {ref_dir}/")
with open(out, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["name", "genome_filename", "protein_filename"])
    for p in files:
        w.writerow([p.name, str(p), ""])
print(f"Sketching {len(files)} reference file(s) from {ref_dir}/", file=sys.stderr)
PYEOF

rm -f "$DB"
sourmash scripts manysketch "$TMP_CSV" -o "$DB" -p "dna,k=${KSIZE},scaled=${SCALED}" -c "$THREADS"
echo "Wrote $DB"
