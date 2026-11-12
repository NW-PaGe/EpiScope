#!/usr/bin/env bash
set -euo pipefail

# Compute mash distances between assemblies and a reference mash DB,
# record assemblies that meet the distance threshold, and for each
# reference determine the closest non-reference assembly.
#
# Requirements: mash, flock, awk, sort, find, xargs, tr, grep
#
# Output location: all CSVs are written to $OUTDIR (default: current dir), and
# the reference DB is read from $OUTDIR/ref_mash_db.msh.

OUTDIR="${OUTDIR:-.}"
mkdir -p "$OUTDIR"

ASSEMBLY_DIR="assemblies"
REF_DB="$OUTDIR/ref_mash_db.msh"
THRESHOLD=0.005
OUT_CSV="$OUTDIR/passed_dist_filter.csv"
REF_NEAREST_OUT="$OUTDIR/ref_nearest_nonref.csv"
REF_RANKED_OUT="$OUTDIR/ref_nonref_ranked.csv"
UNFILTERED_CSV="$OUTDIR/unfiltered_mash.csv"
TMP_DIR=$(mktemp -d)
PARALLEL="${THREADS:-4}"   # parallel assemblies; set THREADS to change
LOCKFILE="$TMP_DIR/prog.lock"
COUNTER_FILE="$TMP_DIR/progress.count"
ALL_DISTS="$TMP_DIR/all.dist.tsv"

cleanup() { rm -rf "$TMP_DIR"; }
trap cleanup EXIT

rm -f "$OUT_CSV" "$REF_NEAREST_OUT" "$REF_RANKED_OUT" "$ALL_DISTS" "$UNFILTERED_CSV"
echo "assembly_file,mash_distance" > "$OUT_CSV"
echo "reference,nearest_nonref_assembly,mash_distance" > "$REF_NEAREST_OUT"
echo "reference,rank,assembly,mash_distance" > "$REF_RANKED_OUT"
# unfiltered long form CSV: assembly_file,reference,mash_distance
echo "assembly_file,reference,mash_distance" > "$UNFILTERED_CSV"
echo 0 > "$COUNTER_FILE"

# Gather assemblies safely into an array
mapfile -d '' -t ASMS < <(find "$ASSEMBLY_DIR" -maxdepth 1 -type f -print0)
TOTAL=${#ASMS[@]}
if [ "$TOTAL" -eq 0 ]; then
  echo "No files found in $ASSEMBLY_DIR"
  exit 0
fi

# Progress display background job
show_progress() {
  while :; do
    [ ! -f "$COUNTER_FILE" ] && break
    count=$(cat "$COUNTER_FILE" 2>/dev/null || echo 0)
    count=${count:-0}
    if [ "$TOTAL" -gt 0 ]; then percent=$((100 * count / TOTAL)); else percent=0; fi
    filled=$(( percent * 50 / 100 )); empty=$((50 - filled))
    if [ "$filled" -gt 0 ]; then filled_seg=$(printf '#%.0s' $(seq 1 $filled)); else filled_seg=""; fi
    if [ "$empty" -gt 0 ]; then empty_seg=$(printf ' %.0s' $(seq 1 $empty)); else empty_seg=""; fi
    printf '\rProgress: [%s%s] %d/%d (%d%%)' "$filled_seg" "$empty_seg" "$count" "$TOTAL" "$percent"
    if [ "$count" -ge "$TOTAL" ]; then printf '\n'; break; fi
    sleep 0.5
  done
}
show_progress & progress_pid=$!

# Worker: sketch assembly, compute distances to REF_DB, append to ALL_DISTS, UNFILTERED_CSV and OUT_CSV
process_one() {
  asm="$1"
  asm_base=$(basename "$asm")
  sketch_prefix="$TMP_DIR/${asm_base%.*}"
  sketch_file="${sketch_prefix}.msh"
  dist_out="$TMP_DIR/${asm_base}.dist"
  norm_out="$TMP_DIR/${asm_base}.norm"

  mash sketch -o "$sketch_prefix" "$asm" >/dev/null 2>&1
  mash dist "$REF_DB" "$sketch_file" > "$dist_out" 2>/dev/null || true

  # Normalize mash's raw dist output into exactly 3 tab-separated fields:
  # ref<TAB>assembly_basename<TAB>distance.
  #
  # Mash's own output is normally 5 tab-separated fields (ref, query,
  # distance, pvalue, shared-hashes), but some reference IDs in this DB
  # carry messy FASTA header metadata that can include an embedded literal
  # tab or space, which shifts every field after it. Reconstruct robustly
  # by treating the LAST 4 whitespace/tab tokens as
  # query/distance/pvalue/hashes (order is fixed by mash) and joining
  # everything before that back into a single ref token, so nothing
  # downstream can ever misalign on the distance field.
  awk -F'\t' -v asm="$asm_base" '
    NF >= 5 {
      n = NF
      dist = $(n-2)
      ref = $1
      for (i = 2; i <= n-4; i++) { ref = ref "_" $i }
      printf "%s\t%s\t%f\n", ref, asm, dist
    }
  ' "$dist_out" > "$norm_out" 2>/dev/null || true

  cat "$norm_out" >> "$ALL_DISTS"

  # also append long-form CSV unfiltered (assembly_file,reference,mash_distance)
  awk -F'\t' 'BEGIN{OFS=","} { ref=$1; gsub(/"/, "\"\"", ref); printf "\"%s\",\"%s\",%.6f\n", $2, ref, $3 }' "$norm_out" >> "$UNFILTERED_CSV" || true

  # if any distance <= threshold, record first match to OUT_CSV
  awk -F'\t' -v t="$THRESHOLD" '$3 <= t { printf "%s,%.6f\n", $2, $3; exit }' "$norm_out" >> "$OUT_CSV"

  # increment counter atomically
  (
    flock 200
    [ ! -f "$COUNTER_FILE" ] && echo 0 > "$COUNTER_FILE"
    cnt=$(cat "$COUNTER_FILE"); cnt=$((cnt + 1)); echo "$cnt" > "$COUNTER_FILE"
  ) 200>"$LOCKFILE"
}
export -f process_one
export TMP_DIR REF_DB THRESHOLD OUT_CSV ALL_DISTS LOCKFILE COUNTER_FILE UNFILTERED_CSV

# Run in parallel
printf '%s\0' "${ASMS[@]}" | xargs -0 -n1 -P "$PARALLEL" -I{} bash -c 'process_one "$@"' _ {}

wait "$progress_pid"

# Ensure ALL_DISTS exists
if [ ! -s "$ALL_DISTS" ]; then
  echo "No distance records produced; exiting."
  exit 1
fi

# Create an awk script file to clean reference names without regex bracket literals
CLEAN_AWK="$TMP_DIR/clean_refs.awk"
cat > "$CLEAN_AWK" <<'AWK'
BEGIN { FS = "\t"; OFS = "\n" }
{
  ref = $1

  # remove bracketed content [ ... ] by locating '[' and matching next ']'
  while (1) {
    i = index(ref, "[")
    if (i == 0) break
    j = index(substr(ref, i+1), "]")
    if (j == 0) break
    # j is offset inside substr; convert to absolute index of closing bracket
    j = i + j
    ref = (i > 1 ? substr(ref, 1, i-1) : "") substr(ref, j+1)
  }

  # if origname= exists, take substring after it up to first whitespace
  p = index(ref, "origname=")
  if (p > 0) {
    ref = substr(ref, p + 8)  # 8 is length of "origname="
    # cut at first whitespace
    sp = index(ref, " ")
    if (sp > 0) ref = substr(ref, 1, sp-1)
  }

  # extract basename if path-like
  n = split(ref, parts, "/")
  ref = parts[n]

  # remove a small set of unwanted single chars by simple loop
  g = ""
  for (k = 1; k <= length(ref); k++) {
    c = substr(ref, k, 1)
    if (c == "\"" || c == "'" || c == "," || c == "`" || c == "\\" ) continue
    g = g c
  }
  ref = g

  # trim leading/trailing whitespace characters (space, tab, CR, LF)
  # leading
  while (length(ref) > 0 && (substr(ref,1,1) == " " || substr(ref,1,1) == "\t" || substr(ref,1,1) == "\r" || substr(ref,1,1) == "\n")) {
    ref = substr(ref, 2)
  }
  # trailing
  while (length(ref) > 0) {
    last = substr(ref, length(ref), 1)
    if (last == " " || last == "\t" || last == "\r" || last == "\n") {
      ref = substr(ref, 1, length(ref)-1)
    } else break
  }

  if (ref == "") next
  # skip header-like tokens containing ':' or '(' or ')'
  if (index(ref, ":") || index(ref, "(") || index(ref, ")")) next

  # accept names starting with alnum or _ . -
  firstc = substr(ref,1,1)
  if ( (firstc >= "0" && firstc <= "9") || (firstc >= "A" && firstc <= "Z") || (firstc >= "a" && firstc <= "z") || firstc == "_" || firstc == "." || firstc == "-" ) {
    print ref
  }
}
AWK

REF_NAMES="$TMP_DIR/ref_names.txt"
awk -f "$CLEAN_AWK" "$ALL_DISTS" | sort -u > "$REF_NAMES"

# Fallback if no cleaned refs found
if [ ! -s "$REF_NAMES" ]; then
  awk -F'\t' '{print $1}' "$ALL_DISTS" | sort -u > "$REF_NAMES"
fi

# Create another awk script for matching and selecting nearest non-ref lines (again no bracket regex)
MATCH_AWK="$TMP_DIR/match_ref.awk"
cat > "$MATCH_AWK" <<'AWK'
BEGIN { FS = "\t" }
{
  refraw = $1
  refclean = refraw

  # remove bracketed content
  while (1) {
    i = index(refclean, "[")
    if (i == 0) break
    j = index(substr(refclean, i+1), "]")
    if (j == 0) break
    j = i + j
    refclean = (i > 1 ? substr(refclean, 1, i-1) : "") substr(refclean, j+1)
  }

  # origname= handling
  p = index(refclean, "origname=")
  if (p > 0) {
    refclean = substr(refclean, p + 8)
    sp = index(refclean, " ")
    if (sp > 0) refclean = substr(refclean, 1, sp-1)
  }

  n = split(refclean, a, "/"); refclean = a[n]
  # trim leading/trailing whitespace
  while (length(refclean) > 0 && (substr(refclean,1,1) == " " || substr(refclean,1,1) == "\t" || substr(refclean,1,1) == "\r" || substr(refclean,1,1) == "\n")) refclean = substr(refclean,2)
  while (length(refclean) > 0) {
    last = substr(refclean, length(refclean), 1)
    if (last == " " || last == "\t" || last == "\r" || last == "\n") refclean = substr(refclean, 1, length(refclean)-1)
    else break
  }

  if (refclean == target) print $0
}
AWK

# For each reference, rank all non-reference assemblies by distance (ascending).
# Write full ranked list to REF_RANKED_OUT, and mirror rank 1 into REF_NEAREST_OUT
# (kept for backward compatibility with anything consuming the old single-nearest file).
#
# Formatting is done entirely inside awk (using $2/$3 directly) rather than
# via bash `read`, so it can't be corrupted even if a reference name still
# ends up containing an unexpected delimiter somewhere upstream.
while IFS= read -r ref; do
  [ -z "$ref" ] && continue

  awk -v target="$ref" -f "$MATCH_AWK" "$ALL_DISTS" \
    | sort -t $'\t' -k3,3n \
    | awk -F'\t' -v ref="$ref" -v ranked_out="$REF_RANKED_OUT" -v nearest_out="$REF_NEAREST_OUT" '
        $2 != ref {
          rank++
          printf "%s,%d,%s,%.6f\n", ref, rank, $2, $3 >> ranked_out
          if (rank == 1) {
            printf "%s,%s,%.6f\n", ref, $2, $3 >> nearest_out
          }
        }
        END {
          if (rank == 0) {
            printf "%s,%s,%s\n", ref, "", "" >> nearest_out
          }
        }
      '
done < "$REF_NAMES"

echo "Filtered results saved to $OUT_CSV"
echo "Nearest non-reference per reference saved to $REF_NEAREST_OUT"
echo "Ranked non-reference candidates per reference saved to $REF_RANKED_OUT"
echo "Unfiltered long-form mash distances saved to $UNFILTERED_CSV"