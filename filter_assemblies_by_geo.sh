#!/usr/bin/env bash
set -euo pipefail

# Optional step: keep only the downloaded assemblies that were collected in the
# location(s) you want. Assemblies that don't fit are moved out of assemblies/
# into a tar file (default: filtered_out.tar.gz), so nothing is lost.
#
# It reads NCBI's assembly_data_report.jsonl (written to assemblies/ by
# get-taxon-genomes.sh) and compares each genome's BioSample "geo_loc_name"
# against your --geo-location value(s).
#
# Usage:
#   filter_assemblies_by_geo.sh --geo-location TEXT [--geo-location TEXT ...] [options]
#
# Always preview first:
#   filter_assemblies_by_geo.sh --geo-location "USA*" --dry-run
#
# Matching (same rules as before):
#   * case-insensitive substring match: "USA" also matches "USA: Washington"
#   * a "*" inside TEXT matches any characters ("USA*Seattle")
#   * spacing around the colon is ignored ("USA: Washington" = "USA:Washington")
#   * repeat the flag to accept several locations (a genome matching any is kept)
#   * genomes with no geo_loc_name in the report do NOT match, so they are moved
#
# Which files are touched:
#   Only files whose names start with a GenBank/RefSeq accession (GCA_/GCF_...)
#   and that have a record in the report. Everything else -- the report itself,
#   dataset_catalog.json, and your own sequences with other names -- is never
#   moved. Files that look like NCBI genomes but are missing from the report are
#   also left alone (and counted), since their location can't be checked.
#
# Requirements: bash, jq (>= 1.6), tar, awk

ASSEMBLY_DIR="assemblies"
REPORT=""
TAR_OUT="filtered_out.tar.gz"
DRY_RUN=0
ALLOW_EMPTY=0
GEO_FILTERS=()

usage() {
  cat <<'EOF'
Usage: filter_assemblies_by_geo.sh --geo-location TEXT [--geo-location TEXT ...] [options]

Required:
  --geo-location TEXT   Keep genomes whose BioSample geo_loc_name contains TEXT
                        (case-insensitive; "*" matches any characters; spacing
                        around ":" ignored). Repeat to accept several locations.

Options:
  --dry-run             Show what would be moved and kept; change nothing
  --assembly-dir DIR    Directory of downloaded assemblies (default: assemblies)
  --report FILE         NCBI report (default: <assembly-dir>/assembly_data_report.jsonl)
  --tar FILE            Tar file for the removed assemblies
                        (default: filtered_out.tar.gz). Must not already exist.
  --allow-empty         Allow a filter that would remove every assembly
  -h, --help            Show this help

Genomes that don't match are moved into the tar and deleted from the assembly
directory only after the tar has been written and verified. The report file is
left unchanged.
EOF
}

need_val() { [ "$#" -ge 2 ] && [ -n "$2" ] || { echo "Error: $1 requires a value." >&2; exit 1; }; }

while [ "$#" -gt 0 ]; do
  case "$1" in
    --geo-location)      need_val "$@"; GEO_FILTERS+=("$2"); shift 2 ;;
    --geo-location=*)    GEO_FILTERS+=("${1#--geo-location=}"); shift ;;
    --assembly-dir)      need_val "$@"; ASSEMBLY_DIR="${2%/}"; shift 2 ;;
    --report)            need_val "$@"; REPORT="$2"; shift 2 ;;
    --tar)               need_val "$@"; TAR_OUT="$2"; shift 2 ;;
    --dry-run|--dryrun)  DRY_RUN=1; shift ;;
    --allow-empty)       ALLOW_EMPTY=1; shift ;;
    -h|--help)           usage; exit 0 ;;
    *) echo "Error: unknown argument '$1'. Run with --help for usage." >&2; exit 1 ;;
  esac
done

[ "${#GEO_FILTERS[@]}" -gt 0 ] || { echo "Error: --geo-location is required (at least once)." >&2; exit 1; }
for g in "${GEO_FILTERS[@]}"; do
  [ -n "${g//[[:space:]]/}" ] || { echo "Error: --geo-location cannot be empty." >&2; exit 1; }
done
[ -n "$REPORT" ] || REPORT="$ASSEMBLY_DIR/assembly_data_report.jsonl"

for tool in jq tar awk; do
  command -v "$tool" >/dev/null 2>&1 || { echo "Error: $tool is required but not found on PATH." >&2; exit 1; }
done
[ -d "$ASSEMBLY_DIR" ] || { echo "Error: assembly directory '$ASSEMBLY_DIR' not found." >&2; exit 1; }
[ -f "$REPORT" ] || {
  echo "Error: report file '$REPORT' not found." >&2
  echo "It is written to the assembly directory by get-taxon-genomes.sh; pass --report if it is elsewhere." >&2
  exit 1
}

# The filter uses jq's --args / $ARGS, which need jq 1.6 or newer.
JQ_VERSION=$(jq --version 2>/dev/null || true)
JQ_MAJOR_MINOR=$(printf '%s' "$JQ_VERSION" | sed -nE 's/^jq-([0-9]+)\.([0-9]+).*/\1 \2/p')
read -r JQ_MAJOR JQ_MINOR <<< "${JQ_MAJOR_MINOR:-0 0}"
if [ "$JQ_MAJOR" -lt 1 ] || { [ "$JQ_MAJOR" -eq 1 ] && [ "$JQ_MINOR" -lt 6 ]; }; then
  echo "Error: this script needs jq 1.6 or newer (found: ${JQ_VERSION:-unknown})." >&2
  exit 1
fi

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"; [ -z "${PARTIAL:-}" ] || rm -f "$PARTIAL"' EXIT
PARTIAL=""

# BioSample location is an attribute named geo_loc_name, e.g.
# "USA: Washington, Seattle", "Colombia:Gambita" or just "France". The value is
# read from that attribute and, as a fallback, from the BioSample's geoLocName
# field. Case and spacing around ':' are normalised before matching.
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
    [ ( (.assemblyInfo.biosample.attributes // .biosample.attributes // [])[]?
        | select(.name == "geo_loc_name") | .value ),
      (.assemblyInfo.biosample.geoLocName // .biosample.geoLocName)
    ] | map(select(. != null and . != "")) | unique;
  def geo_values: geo_raw | map(norm);
'

# 1. One line per report record: accession, paired accession, 1/0 match, location
jq -r "$GEO_JQ_DEFS"'
  select(.accession != null)
  | ($ARGS.positional | map(norm)) as $p
  | geo_values as $g
  | [ .accession,
      (.pairedAccession // ""),
      (if any($g[]; . as $v | any($p[]; . as $q | $v | geo_match($q))) then "1" else "0" end),
      ((geo_raw | join("; ")) | if . == "" then "(no location)" else . end) ]
  | @tsv' --args "${GEO_FILTERS[@]}" < "$REPORT" > "$TMP/report.tsv"

N_RECORDS=$(wc -l < "$TMP/report.tsv")
[ "$N_RECORDS" -gt 0 ] || { echo "Error: no usable records found in '$REPORT'." >&2; exit 1; }

# 2. Files in the assembly directory
for f in "$ASSEMBLY_DIR"/*; do
  if [ -f "$f" ]; then printf '%s\n' "${f##*/}"; fi
done > "$TMP/files.txt"

# 3. Decide what to do with each file: keep | move | norecord | skip
#    (matched by the GCA_/GCF_ accession at the start of the file name, against
#    either the record's accession or its paired accession)
awk -F'\t' -v OFS='\t' '
  FNR == NR {
    st[$1] = $3; loc[$1] = $4
    if ($2 != "" && !($2 in st)) { st[$2] = $3; loc[$2] = $4 }
    next
  }
  {
    name = $0
    if (match(name, /^GC[AF]_[0-9]+\.[0-9]+/)) {
      acc = substr(name, RSTART, RLENGTH)
      if (acc in st) print name, (st[acc] == "1" ? "keep" : "move"), loc[acc]
      else           print name, "norecord", ""
    } else {
      print name, "skip", ""
    }
  }
' "$TMP/report.tsv" "$TMP/files.txt" > "$TMP/decisions.tsv"

count() { awk -F'\t' -v a="$1" '$2 == a { n++ } END { print n + 0 }' "$TMP/decisions.tsv"; }
N_MOVE=$(count move); N_KEEP=$(count keep); N_NOREC=$(count norecord); N_SKIP=$(count skip)

echo "Filter: geo_loc_name contains any of: $(printf '"%s" ' "${GEO_FILTERS[@]}")"
echo "Report: $REPORT ($N_RECORDS record(s))   Assembly dir: $ASSEMBLY_DIR/"
echo
echo "  keep (match the filter):            $N_KEEP"
echo "  move to tar (do not match):         $N_MOVE"
echo "  left alone, not in the report:      $N_NOREC"
echo "  left alone, not an NCBI genome file: $N_SKIP"

if [ "$N_NOREC" -gt 0 ]; then
  echo >&2
  echo "Warning: $N_NOREC genome file(s) have no record in the report, so their location can't be" >&2
  echo "checked and they are being left alone. This happens if the report only covers part of" >&2
  echo "the download (e.g. a download split into several batches overwrote it)." >&2
fi

show_locations() {  # $1 = action, $2 = heading
  echo "$2"
  awk -F'\t' -v a="$1" '$2 == a { c[$3]++ } END { for (k in c) printf "%7d  %s\n", c[k], k }' "$TMP/decisions.tsv" \
    | sort -rn | head -10 | sed 's/^/  /'
}

EMPTY_FILTER=0
if [ "$N_KEEP" -eq 0 ] && [ "$N_MOVE" -gt 0 ]; then EMPTY_FILTER=1; fi

# ---------- dry run ----------
if [ "$DRY_RUN" -eq 1 ]; then
  echo
  echo "DRY RUN: nothing will be moved, written or deleted."
  if [ "$N_MOVE" -gt 0 ]; then
    echo
    echo "Would move to $TAR_OUT (file, location):"
    awk -F'\t' '$2 == "move" { printf "  %s\t%s\n", $1, $3 }' "$TMP/decisions.tsv"
  fi
  if [ "$N_KEEP" -gt 0 ]; then
    echo
    show_locations keep "Locations of the $N_KEEP genome(s) that would be kept (top 10):"
  fi
  if [ "$N_MOVE" -gt 0 ]; then
    echo
    show_locations move "Locations of the $N_MOVE genome(s) that would be moved (top 10):"
  fi
  if [ "$EMPTY_FILTER" -eq 1 ] && [ "$ALLOW_EMPTY" -eq 0 ]; then
    echo
    echo "Warning: this filter keeps NONE of the genomes. A real run would refuse unless you add"
    echo "--allow-empty. Check the spelling against the locations listed above."
  fi
  if [ "$N_MOVE" -gt 0 ] && [ -e "$TAR_OUT" ]; then
    echo
    echo "Warning: $TAR_OUT already exists. A real run would refuse to overwrite it; choose another"
    echo "name with --tar or move the existing file away."
  fi
  exit 0
fi

# ---------- real run ----------
if [ "$N_MOVE" -eq 0 ]; then
  echo
  echo "Nothing to move: every genome in the report that is in $ASSEMBLY_DIR/ matches the filter."
  exit 0
fi

if [ "$EMPTY_FILTER" -eq 1 ] && [ "$ALLOW_EMPTY" -eq 0 ]; then
  echo >&2
  echo "Error: this filter keeps none of the $N_MOVE genome(s); refusing to move them all." >&2
  show_locations move "Locations present (top 10):" >&2
  echo "Check the spelling of --geo-location (re-run with --dry-run), or pass --allow-empty" >&2
  echo "if removing everything is really what you want." >&2
  exit 1
fi

if [ -e "$TAR_OUT" ]; then
  echo >&2
  echo "Error: $TAR_OUT already exists; refusing to overwrite it (it may hold genomes that were" >&2
  echo "previously moved out of $ASSEMBLY_DIR/). Choose another name with --tar, or move it away." >&2
  exit 1
fi

mkdir -p "$(dirname "$TAR_OUT")"
awk -F'\t' '$2 == "move" { print $1 }' "$TMP/decisions.tsv" | tr '\n' '\0' > "$TMP/move.lst"

echo
echo "Writing $TAR_OUT ..."
PARTIAL="$TAR_OUT.partial"
tar -czf "$PARTIAL" -C "$ASSEMBLY_DIR" --null -T "$TMP/move.lst"

# Verify the archive is readable and holds every file before deleting anything.
N_IN_TAR=$(tar -tzf "$PARTIAL" | wc -l)
if [ "$N_IN_TAR" -ne "$N_MOVE" ]; then
  echo "Error: archive check failed (expected $N_MOVE file(s), found $N_IN_TAR). Nothing was deleted." >&2
  exit 1
fi
mv "$PARTIAL" "$TAR_OUT"
PARTIAL=""

( cd "$ASSEMBLY_DIR" && xargs -0 rm -f -- < "$TMP/move.lst" )

echo "Moved $N_MOVE genome file(s) from $ASSEMBLY_DIR/ into $TAR_OUT (verified before deleting)."
echo "$N_KEEP genome file(s) remain. The report $REPORT was not changed."
echo "To restore them: tar -xzf $TAR_OUT -C $ASSEMBLY_DIR"
