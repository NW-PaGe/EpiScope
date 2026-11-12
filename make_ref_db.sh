#!/bin/bash

# Sketch cluster_references/ into a Mash DB.
# Outputs (cluster_ref_list.txt, ref_mash_db.msh) go to $OUTDIR (default: current dir).
OUTDIR="${OUTDIR:-.}"
mkdir -p "$OUTDIR"

rm -rf "$OUTDIR/cluster_ref_list.txt"
rm -f "$OUTDIR"/ref_mash_db*
ls cluster_references/ | sed 's|^|cluster_references/|' > "$OUTDIR/cluster_ref_list.txt"
mash sketch -o "$OUTDIR/ref_mash_db" -l "$OUTDIR/cluster_ref_list.txt"
