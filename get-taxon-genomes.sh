#!/bin/bash
set -euo pipefail

# Download all genomes matching a taxon/filter query, batching the download
# so no single `datasets download` call exceeds NCBI's ~50k genome ceiling
# (which otherwise causes NCBI's dehydrated download to fail/crash).
#
# Before downloading anything, this also estimates the total disk space the
# full matched genome set will require and compares it against space
# available on the target filesystem. If it won't fit, the script refuses
# to start -- there's no point downloading and rehydrating a dataset you
# don't have room to finish. It also re-checks available space before each
# batch, in case something else on the machine eats disk space mid-run.
#
# Usage: same as before -- pass whatever you'd normally pass to
#   `datasets download genome taxon <taxon> [filters...]`
# e.g.
#   ./get-taxon-genomes.sh "Staphylococcus aureus" --reference
#   ./get-taxon-genomes.sh "Escherichia coli" --assembly-level complete
#
# Extra options handled by this script (everything else is passed to `datasets`):
#   --geo-location TEXT   only download genomes collected in a matching place.
#   (or --geolocation)    Matches TEXT against each genome's BioSample
#                         "geo_loc_name": case-insensitive substring match, so
#                         "USA" also matches "USA: Washington"; a "*" inside TEXT
#                         matches any characters ("USA*Seattle"); spacing around
#                         the colon is ignored. Repeat the flag to accept several
#                         locations (OR). Genomes with no location are excluded.
#                         The filter is applied to the `datasets summary` results
#                         BEFORE anything is downloaded, so only matching genomes
#                         are fetched (and only they count toward the disk check).
#   --accessions-out FILE write the final accession list (one per line) to FILE.
#   --dry-run             resolve the genome list, apply any filters and run the
#                         disk space check, then stop without downloading.
# e.g.
#   ./get-taxon-genomes.sh "Haemophilus influenzae" --geo-location "USA*" --dry-run
#   ./get-taxon-genomes.sh "Escherichia coli" --assembly-level complete \
#       --geo-location Canada --geo-location "USA: Oregon" --accessions-out accessions.txt
#
# To filter genomes you have ALREADY downloaded, see filter_assemblies_by_geo.sh.
#
# Optional environment variable overrides:
#   BATCH_SIZE            genomes per download batch (default 40000, kept
#                          comfortably under NCBI's ~50k limit)
#   SPACE_SAFETY_MARGIN    multiplier applied to the estimated space
#                          requirement as a buffer (default 1.3, i.e. +30%)
#   SKIP_SPACE_CHECK       set to 1 to bypass the disk space check entirely
#
# e.g.
#   BATCH_SIZE=25000 SPACE_SAFETY_MARGIN=1.5 ./get-taxon-genomes.sh "Escherichia coli"
#
# Requirements: datasets, dataformat, unzip, jq, find, mv, split, df, awk

BATCH_SIZE="${BATCH_SIZE:-40000}"
SPACE_SAFETY_MARGIN="${SPACE_SAFETY_MARGIN:-1.2}"
SKIP_SPACE_CHECK="${SKIP_SPACE_CHECK:-0}"

if [ "$#" -eq 0 ]; then
  echo "Usage: $0 <taxon> [datasets filter flags...]" >&2
  exit 1
fi

for tool in jq dataformat df awk; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "Error: $tool is required but not found on PATH." >&2
    exit 1
  fi
done

# Determine target directory the same way as before, based on the raw args.
if [[ "$@" == *"--reference"* ]]; then
  directory="reference_genomes"
else
  directory="assemblies"
fi

# Split "$@" into:
#   FILTER_ARGS - taxon + filter flags (--reference, --assembly-level, etc.),
#                 passed to `datasets summary genome taxon` to build the
#                 accession list.
#   DL_ARGS     - flags that control *what data* is fetched per genome
#                 (currently just --include), which `summary` doesn't
#                 understand but `download genome accession` does, and
#                 must be forwarded to the batched download step.
FILTER_ARGS=()
DL_ARGS=()
GEO_FILTERS=()
ACCESSIONS_OUT=""
DRY_RUN=0
args=("$@")
n=${#args[@]}
i=0
while [ "$i" -lt "$n" ]; do
  arg="${args[$i]}"
  case "$arg" in
    --include)
      DL_ARGS+=("$arg" "${args[$((i + 1))]}")
      i=$((i + 2))
      continue
      ;;
    --include=*)
      DL_ARGS+=("$arg")
      i=$((i + 1))
      continue
      ;;
    --geo-location|--geolocation)
      if [ "$((i + 1))" -ge "$n" ]; then
        echo "Error: $arg requires a value." >&2
        exit 1
      fi
      GEO_FILTERS+=("${args[$((i + 1))]}")
      i=$((i + 2))
      continue
      ;;
    --geo-location=*|--geolocation=*)
      GEO_FILTERS+=("${arg#*=}")
      i=$((i + 1))
      continue
      ;;
    --accessions-out)
      if [ "$((i + 1))" -ge "$n" ]; then
        echo "Error: --accessions-out requires a file name." >&2
        exit 1
      fi
      ACCESSIONS_OUT="${args[$((i + 1))]}"
      i=$((i + 2))
      continue
      ;;
    --accessions-out=*)
      ACCESSIONS_OUT="${arg#--accessions-out=}"
      i=$((i + 1))
      continue
      ;;
    --dry-run|--dryrun)
      DRY_RUN=1
      i=$((i + 1))
      continue
      ;;
    *)
      FILTER_ARGS+=("$arg")
      i=$((i + 1))
      ;;
  esac
done

for g in ${GEO_FILTERS[@]+"${GEO_FILTERS[@]}"}; do
  if [ -z "${g//[[:space:]]/}" ]; then
    echo "Error: --geo-location cannot be empty." >&2
    exit 1
  fi
done

# A dry run must not create anything; measure free space on the nearest
# directory that exists instead.
SPACE_DIR="$directory"
if [ "$DRY_RUN" -eq 1 ]; then
  [ -d "$directory" ] || SPACE_DIR="."
else
  mkdir -p "$directory"
fi

echo "Resolving genome list for: ${FILTER_ARGS[*]}"
RAW_JSON=$(mktemp)
datasets summary genome taxon "${FILTER_ARGS[@]}" --as-json-lines > "$RAW_JSON"

if [ "${#GEO_FILTERS[@]}" -gt 0 ]; then
  # The filter uses jq's --args / $ARGS, which need jq 1.6 or newer.
  JQ_VERSION=$(jq --version 2>/dev/null || true)
  JQ_MAJOR_MINOR=$(printf '%s' "$JQ_VERSION" | sed -nE 's/^jq-([0-9]+)\.([0-9]+).*/\1 \2/p')
  read -r JQ_MAJOR JQ_MINOR <<< "${JQ_MAJOR_MINOR:-0 0}"
  if [ "$JQ_MAJOR" -lt 1 ] || { [ "$JQ_MAJOR" -eq 1 ] && [ "$JQ_MINOR" -lt 6 ]; }; then
    echo "Error: --geo-location needs jq 1.6 or newer (found: ${JQ_VERSION:-unknown})." >&2
    rm -f "$RAW_JSON"
    exit 1
  fi

  # Where the location sits in a record depends on the datasets version and
  # report type, so look for it anywhere in the record: a "geo_loc_name"
  # attribute ({"name":"geo_loc_name","value":"USA"}), a "geo_loc_name" key
  # ("geo_loc_name":"USA"), or a "geoLocName" key. Examples of values:
  # "USA: Washington, Seattle", "Colombia:Gambita", "France". Case and spacing
  # around ':' are normalised before matching.
  GEO_JQ_DEFS='
    def norm: ascii_downcase | gsub("\\s*:\\s*"; ":") | gsub("^\\s+|\\s+$"; "");
    # Escape regex metacharacters (except *, which is handled separately).
    def rx_escape:
      explode
      | map(. as $c
            | if any(46,94,36,124,63,43,40,41,123,125,91,93,92; . == $c)
              then [92, $c] else [$c] end)
      | add // [] | implode;
    # "usa" matches anywhere in the value; a "*" inside the pattern matches any
    # run of characters, so "usa*seattle" matches "usa:washington,seattle".
    def geo_match($pat): test($pat | split("*") | map(rx_escape) | join(".*"));
    def geo_raw:
      [ ( .. | objects | select(.name? == "geo_loc_name") | .value? ),
        ( .. | objects | (.geo_loc_name?, .geoLocName?) )
      ] | map(select(type == "string" and . != "")) | unique;
    def geo_values: geo_raw | map(norm);
  '

  N_ALL=$(wc -l < "$RAW_JSON")
  GEO_JSON=$(mktemp)
  jq -c "$GEO_JQ_DEFS"'
    ($ARGS.positional | map(norm)) as $p
    | select( any(geo_values[]; . as $v | any($p[]; . as $q | $v | geo_match($q))) )
  ' --args "${GEO_FILTERS[@]}" < "$RAW_JSON" > "$GEO_JSON"
  N_NO_GEO=$(jq -c "$GEO_JQ_DEFS"' select(geo_values | length == 0)' < "$RAW_JSON" | wc -l)
  N_KEPT=$(wc -l < "$GEO_JSON")

  echo "Geo filter (location contains any of: $(printf '"%s" ' "${GEO_FILTERS[@]}"))"
  echo "  kept $N_KEPT of $N_ALL genome(s); $N_NO_GEO had no location and cannot match."

  if [ "$N_ALL" -gt 0 ] && [ "$N_NO_GEO" -eq "$N_ALL" ]; then
    echo "" >&2
    echo "Error: none of the $N_ALL record(s) from 'datasets summary' contain a location, so" >&2
    echo "--geo-location can't be applied before downloading. Inspect one record with:" >&2
    echo "  datasets summary genome taxon <taxon> --as-json-lines | head -n 1 | jq ." >&2
    echo "Or download without the filter and then run filter_assemblies_by_geo.sh, which" >&2
    echo "uses the report that comes with the download." >&2
    rm -f "$RAW_JSON" "$GEO_JSON"
    exit 1
  fi

  if [ "$N_KEPT" -eq 0 ] && [ "$N_ALL" -gt 0 ]; then
    echo "  Nothing matched. Locations present in the results (count, most common first):" >&2
    jq -r "$GEO_JQ_DEFS"' geo_raw[]' < "$RAW_JSON" | sort | uniq -c | sort -rn | head -15 | sed 's/^/    /' >&2
    echo "  Check the spelling in your --geo-location values against these." >&2
  fi
  mv "$GEO_JSON" "$RAW_JSON"
fi

ACC_FILE=$(mktemp)
jq -r '.accession' "$RAW_JSON" > "$ACC_FILE"

TOTAL=$(wc -l < "$ACC_FILE")
echo "Found $TOTAL genome(s) matching query${GEO_FILTERS[0]:+ and geo filter}."

if [ "$TOTAL" -eq 0 ]; then
  echo "No genomes found; nothing to download."
  if [ -n "$ACCESSIONS_OUT" ]; then
    mkdir -p "$(dirname "$ACCESSIONS_OUT")"
    : > "$ACCESSIONS_OUT"
    echo "Wrote an empty accession list to $ACCESSIONS_OUT"
  fi
  rm -f "$RAW_JSON" "$ACC_FILE"
  exit 0
fi

if [ -n "$ACCESSIONS_OUT" ]; then
  mkdir -p "$(dirname "$ACCESSIONS_OUT")"
  cp "$ACC_FILE" "$ACCESSIONS_OUT"
  echo "Wrote $TOTAL accession(s) to $ACCESSIONS_OUT"
fi

if [ "$SKIP_SPACE_CHECK" != "1" ]; then
  echo "Estimating disk space required for $TOTAL genome(s)..."

  # Total bases across all matched assemblies, from NCBI's own metadata.
  # This is a size proxy for the uncompressed FASTA data. jq extraction is
  # unreliable across schema versions, so use dataformat (the officially
  # supported way to pull fields from a summary report).
  TOTAL_BASES=$(
    dataformat tsv genome --inputfile "$RAW_JSON" \
      --fields accession,assmstats-total-sequence-len --elide-header \
      | awk -F'\t' '{ if ($2 ~ /^[0-9]+$/) sum += $2 } END { print sum+0 }'
  )

  if [ "$TOTAL_BASES" -eq 0 ]; then
    echo "Warning: could not determine assembly sizes from NCBI metadata," >&2
    echo "so the disk space check can't run. Set SKIP_SPACE_CHECK=1 to" >&2
    echo "proceed anyway, or investigate why assmstats-total-sequence-len" >&2
    echo "came back empty for this query." >&2
    rm -f "$RAW_JSON" "$ACC_FILE"
    exit 1
  fi

  # FASTA overhead (headers/newlines) is small; ~5% is a safe padding.
  FINAL_BYTES=$(awk -v b="$TOTAL_BASES" 'BEGIN { printf "%.0f", b * 1.05 }')

  # Per-batch bytes, used to estimate peak *transient* extra space needed
  # while a batch's zip and its rehydrated/extracted copy briefly coexist
  # on disk before files are moved into $directory and the temp files removed.
  N_BATCHES=$(( (TOTAL + BATCH_SIZE - 1) / BATCH_SIZE ))
  BATCH_BYTES=$(awk -v b="$FINAL_BYTES" -v n="$N_BATCHES" 'BEGIN { printf "%.0f", b / n }')
  PEAK_EXTRA_BYTES=$(awk -v b="$BATCH_BYTES" 'BEGIN { printf "%.0f", b * 2 }')

  REQUIRED_BYTES=$(awk -v f="$FINAL_BYTES" -v e="$PEAK_EXTRA_BYTES" -v m="$SPACE_SAFETY_MARGIN" \
    'BEGIN { printf "%.0f", (f + e) * m }')

  AVAIL_BYTES=$(df -P -B1 "$SPACE_DIR" | awk 'NR==2 { print $4 }')

  REQUIRED_GB=$(awk -v b="$REQUIRED_BYTES" 'BEGIN { printf "%.1f", b / 1e9 }')
  AVAIL_GB=$(awk -v b="$AVAIL_BYTES" 'BEGIN { printf "%.1f", b / 1e9 }')

  echo "Estimated space needed (with ${SPACE_SAFETY_MARGIN}x safety margin): ${REQUIRED_GB} GB"
  echo "Available space on the filesystem holding '$directory': ${AVAIL_GB} GB"

  if [ "$AVAIL_BYTES" -lt "$REQUIRED_BYTES" ]; then
    echo "" >&2
    echo "Error: not enough disk space to complete this download." >&2
    echo "  Estimated requirement: ${REQUIRED_GB} GB" >&2
    echo "  Available:             ${AVAIL_GB} GB" >&2
    echo "" >&2
    echo "Downloading and rehydrating a dataset you don't have room to" >&2
    echo "finish just wastes time and NCBI bandwidth. Options:" >&2
    echo "  - Free up space or point $directory at a larger volume" >&2
    echo "  - Narrow the query (e.g. --assembly-level, --annotated," >&2
    echo "    --released-after) to match fewer/smaller genomes" >&2
    echo "  - Set SKIP_SPACE_CHECK=1 to bypass this check if you're sure" >&2
    rm -f "$RAW_JSON" "$ACC_FILE"
    exit 1
  fi
fi

if [ "$DRY_RUN" -eq 1 ]; then
  N_BATCHES_DRY=$(( (TOTAL + BATCH_SIZE - 1) / BATCH_SIZE ))
  echo ""
  echo "Dry run: would download $TOTAL genome(s) in $N_BATCHES_DRY batch(es) into $directory/. Nothing was downloaded."
  if [ -z "$ACCESSIONS_OUT" ]; then
    echo "(Add --accessions-out FILE to save the accession list.)"
  fi
  rm -f "$RAW_JSON" "$ACC_FILE"
  exit 0
fi

rm -f "$RAW_JSON"

# Split accessions into batches of BATCH_SIZE, one accession per line.
BATCH_PREFIX="${ACC_FILE}.batch_"
split -l "$BATCH_SIZE" -d -a 4 "$ACC_FILE" "$BATCH_PREFIX"

work_dir=""
trap 'rm -rf "${work_dir:-}"' EXIT

batch_num=0
for batch_file in "${BATCH_PREFIX}"*; do
  batch_num=$((batch_num + 1))
  n_in_batch=$(wc -l < "$batch_file")
  echo "=== Batch $batch_num: $n_in_batch genome(s) ==="

  if [ "$SKIP_SPACE_CHECK" != "1" ]; then
    avail_now=$(df -P -B1 "$directory" | awk 'NR==2 { print $4 }')
    avail_now_gb=$(awk -v b="$avail_now" 'BEGIN { printf "%.1f", b / 1e9 }')
    # Rough per-batch floor: this batch's share of bytes plus its own
    # transient zip+extract overhead. If we're already below that, stop
    # now rather than leave a half-downloaded, half-rehydrated batch.
    batch_floor=$(awk -v b="$BATCH_BYTES" 'BEGIN { printf "%.0f", b * 1.5 }')
    if [ -n "${BATCH_BYTES:-}" ] && [ "$avail_now" -lt "$batch_floor" ]; then
      echo "" >&2
      echo "Error: only ${avail_now_gb} GB free before starting batch $batch_num;" >&2
      echo "that's not enough headroom to safely process this batch. Stopping" >&2
      echo "here rather than risk a corrupted partial batch. $((batch_num - 1))" >&2
      echo "batch(es) already completed and saved in $directory/." >&2
      echo "Free up space and re-run to continue (already-downloaded" >&2
      echo "accessions can be skipped by trimming $ACC_FILE if desired)." >&2
      exit 1
    fi
  fi

  # Download and extract inside an isolated temp directory so NCBI's bundled
  # README.md (and anything else in the zip) can never overwrite or delete
  # files in the current working directory. It is created under ./ so it sits
  # on the same filesystem as $directory, keeping the final mv cheap.
  work_dir=$(mktemp -d "./.ncbi_batch_XXXXXX")
  datasets download genome accession --inputfile "$batch_file" --dehydrated \
    --filename "$work_dir/ncbi_dataset.zip" "${DL_ARGS[@]}"
  unzip -o -q "$work_dir/ncbi_dataset.zip" -d "$work_dir"
  datasets rehydrate --directory "$work_dir"
  find "$work_dir/ncbi_dataset/data/" -type f -exec mv -t "$directory" {} +
  rm -rf "$work_dir"
  work_dir=""

  echo "Batch $batch_num complete. $directory now has $(find "$directory" -type f | wc -l) file(s)."
done

rm -f "$ACC_FILE" "${BATCH_PREFIX}"*

echo "All $batch_num batch(es) complete. Genomes saved to $directory/"