#!/usr/bin/env python3
"""
Sourmash version of get_mash_dist.sh, driven through the sourmash command line
plus the branchwater plugin for the multithreaded parts.

Pipeline:
  1. sourmash scripts manysketch   sketch every assembly (multithreaded)
  2. sourmash scripts multisearch  compare all assemblies against the reference
                                   DB (multithreaded)
  3. post-process into the same four output files as the mash version

Requirements (see environment.yml):
  sourmash, sourmash_plugin_branchwater   (provides `sourmash scripts ...`)

Build the reference DB first with make_ref_db_sourmash.sh. It must be made with
`manysketch` (which names each sketch after its file) rather than plain
`sourmash sketch dna`, which leaves the sketches unnamed.

Inputs : assemblies/            FASTA files (.fa .fasta .fna .fas, optionally .gz);
                                other files (e.g. NCBI .jsonl) are ignored
         $OUTDIR/ref_sourmash_db.zip
Outputs: written to $OUTDIR (same names as the mash version, so downstream
         steps are unchanged): unfiltered_mash.csv, passed_dist_filter.csv,
         ref_nonref_ranked.csv, ref_nearest_nonref.csv

Env overrides: OUTDIR (default .), THREADS (default 4), KSIZE (default 21),
               SCALED (default: the reference DB's scaled; must be >= it),
               THRESHOLD (default 0.005)

Distances are Mash-style distances estimated from the FracMinHash Jaccard that
multisearch reports:  D = -ln(2j / (1 + j)) / k.  Pairs with no shared hashes
are not reported by multisearch (-t 0) and are treated as distance 1.0.
"""
import csv
import math
import os
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

OUTDIR = Path(os.environ.get("OUTDIR", "."))
OUTDIR.mkdir(parents=True, exist_ok=True)

ASSEMBLY_DIR = Path("assemblies")
REF_DB = OUTDIR / "ref_sourmash_db.zip"
THRESHOLD = float(os.environ.get("THRESHOLD", 0.005))
THREADS = str(max(1, int(os.environ.get("THREADS", 4))))
K = int(os.environ.get("KSIZE", 21))

OUT_CSV = OUTDIR / "passed_dist_filter.csv"
REF_NEAREST_OUT = OUTDIR / "ref_nearest_nonref.csv"
REF_RANKED_OUT = OUTDIR / "ref_nonref_ranked.csv"
UNFILTERED_CSV = OUTDIR / "unfiltered_mash.csv"

FASTA_EXTS = (".fa", ".fasta", ".fna", ".fas")


def is_fasta(path):
    name = path.name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    return name.endswith(FASTA_EXTS)


def run(cmd):
    print("+ " + " ".join(map(str, cmd)), file=sys.stderr, flush=True)
    subprocess.run([str(c) for c in cmd], check=True)


def read_manifest(path):
    """Read a `sourmash sig manifest` CSV, which starts with a '#' comment line."""
    with open(path, newline="") as fh:
        lines = [line for line in fh if not line.startswith("#")]
    return list(csv.DictReader(lines))


def clean_ref(raw):
    """Tidy a reference name; None if unusable."""
    ref = raw
    while True:  # remove bracketed content [ ... ]
        i = ref.find("[")
        if i == -1:
            break
        j = ref.find("]", i + 1)
        if j == -1:
            break
        ref = ref[:i] + ref[j + 1:]
    p = ref.find("origname=")
    if p != -1:
        parts = ref[p + len("origname="):].split()
        ref = parts[0] if parts else ""
    ref = ref.split("/")[-1]
    ref = ref.translate(str.maketrans("", "", "\"',`\\")).strip(" \t\r\n")
    if not ref or any(c in ref for c in ":()"):
        return None
    c = ref[0]
    if not (c.isascii() and (c.isalnum() or c in "_.-")):
        return None
    return ref


def mash_distance(j, k):
    if j <= 0:
        return 1.0
    return max(0.0, -math.log(2.0 * j / (1.0 + j)) / k)


def q(s):
    """CSV-quote a string the same way get_mash_dist.sh does."""
    return '"' + s.replace('"', '""') + '"'


def main():
    probe = subprocess.run(["sourmash", "scripts", "multisearch", "--help"],
                           capture_output=True)
    if probe.returncode != 0:
        sys.exit("`sourmash scripts multisearch` not found: install "
                 "sourmash_plugin_branchwater (conda install sourmash_plugin_branchwater)")
    if not REF_DB.exists():
        sys.exit(f"Reference DB not found: {REF_DB} (run make_ref_db_sourmash.sh first)")

    all_files = sorted(p for p in ASSEMBLY_DIR.iterdir() if p.is_file())
    asm_paths = [p for p in all_files if is_fasta(p)]
    if len(all_files) > len(asm_paths):
        print(f"Ignoring {len(all_files) - len(asm_paths)} non-FASTA file(s) in {ASSEMBLY_DIR}",
              file=sys.stderr)
    if not asm_paths:
        print(f"No FASTA files found in {ASSEMBLY_DIR}")
        return 0

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)

        # 0. read the reference DB's manifest: names, k and scaled
        ref_manifest = tmp / "ref_manifest.csv"
        run(["sourmash", "sig", "manifest", REF_DB, "-o", ref_manifest])
        ref_rows = [r for r in read_manifest(ref_manifest) if r.get("ksize") == str(K)]
        if not ref_rows:
            sys.exit(f"No k={K} sketches in {REF_DB}; rebuild it with KSIZE={K} "
                     f"or set KSIZE to match the DB.")
        db_scaled = max(int(r["scaled"]) for r in ref_rows)
        scaled = int(os.environ.get("SCALED", db_scaled))
        if scaled < db_scaled:
            sys.exit(f"SCALED={scaled} is lower than the reference DB's scaled={db_scaled}; "
                     f"rebuild the DB at a lower scaled or raise SCALED.")

        ref_names = set()
        for r in ref_rows:
            raw = r.get("name") or Path(r.get("filename", "")).name
            if not raw:
                sys.exit(f"{REF_DB} contains unnamed sketches; rebuild it with "
                         f"make_ref_db_sourmash.sh (plain `sourmash sketch dna` leaves them unnamed).")
            ref_names.add(clean_ref(raw) or raw)

        # 1. sketch all assemblies, named by file basename
        sketch_csv = tmp / "manysketch.csv"
        with open(sketch_csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["name", "genome_filename", "protein_filename"])
            for p in asm_paths:
                w.writerow([p.name, str(p), ""])
        asm_zip = tmp / "assemblies.zip"
        run(["sourmash", "scripts", "manysketch", sketch_csv, "-o", asm_zip,
             "-p", f"dna,k={K},scaled={scaled}", "-c", THREADS])

        # Only assemblies that were actually sketched (manysketch exits 0 even
        # if a file fails to load, so check the manifest instead of trusting it).
        asm_manifest = tmp / "asm_manifest.csv"
        run(["sourmash", "sig", "manifest", asm_zip, "-o", asm_manifest])
        sketched = {r["name"] for r in read_manifest(asm_manifest)}
        asm_names = [p.name for p in asm_paths if p.name in sketched]
        failed = len(asm_paths) - len(asm_names)
        if failed:
            print(f"WARNING: {failed} assembly file(s) could not be sketched and are "
                  f"excluded from the results", file=sys.stderr)
        if not asm_names:
            print("No assemblies could be sketched; exiting.", file=sys.stderr)
            return 1

        # 2. all assemblies vs reference DB; -t 0 reports any overlap
        results = tmp / "multisearch.csv"
        run(["sourmash", "scripts", "multisearch", asm_zip, REF_DB, "-o", results,
             "-k", K, "-s", scaled, "-t", "0", "-c", THREADS])

        # 3. read results: dist[asm][ref] = best (smallest) distance
        dist = defaultdict(dict)
        with open(results, newline="") as fh:
            for row in csv.DictReader(fh):
                asm = row.get("query_name") or row.get("query")
                raw = row.get("match_name") or row.get("match")
                ref = clean_ref(raw) or raw
                ref_names.add(ref)
                d = mash_distance(float(row["jaccard"]), K)
                if d < dist[asm].get(ref, 2.0):
                    dist[asm][ref] = d

    refs = sorted(ref_names)
    if not refs:
        print("No reference sketches found; exiting.", file=sys.stderr)
        return 1

    # unfiltered long form: assembly_file,reference,mash_distance
    with open(UNFILTERED_CSV, "w") as fh:
        fh.write("assembly_file,reference,mash_distance\n")
        for asm in asm_names:
            row = dist.get(asm, {})
            for ref in refs:
                fh.write(f"{q(asm)},{q(ref)},{row.get(ref, 1.0):.6f}\n")

    # assemblies whose closest reference is within the threshold
    with open(OUT_CSV, "w") as fh:
        fh.write("assembly_file,mash_distance\n")
        for asm in asm_names:
            row = dist.get(asm)
            if row and min(row.values()) <= THRESHOLD:
                fh.write(f"{asm},{min(row.values()):.6f}\n")

    # per reference: rank all non-reference assemblies by distance (ascending)
    with open(REF_RANKED_OUT, "w", newline="") as ranked, \
         open(REF_NEAREST_OUT, "w", newline="") as nearest:
        ranked.write("reference,rank,assembly,mash_distance\n")
        nearest.write("reference,nearest_nonref_assembly,mash_distance\n")
        rw = csv.writer(ranked, lineterminator="\n")
        nw = csv.writer(nearest, lineterminator="\n")
        for ref in refs:
            cands = sorted(
                (dist.get(a, {}).get(ref, 1.0), a) for a in asm_names
                if a != ref and Path(a).stem != ref)
            for rank, (d, a) in enumerate(cands, 1):
                rw.writerow([ref, rank, a, f"{d:.6f}"])
            if cands:
                nw.writerow([ref, cands[0][1], f"{cands[0][0]:.6f}"])
            else:
                nw.writerow([ref, "", ""])

    print(f"Filtered results saved to {OUT_CSV}")
    print(f"Nearest non-reference per reference saved to {REF_NEAREST_OUT}")
    print(f"Ranked non-reference candidates per reference saved to {REF_RANKED_OUT}")
    print(f"Unfiltered long-form distances saved to {UNFILTERED_CSV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
