#!/usr/bin/env bash
set -euo pipefail

# Populate cluster_references/ with reference FASTA files from BigBacter
# databases, stored locally or on S3.
#
# Usage:
#   populate_bb_references.sh --config CONFIG.csv --species SPECIES [options]
#
# The config CSV has two columns, "db" and "path" (header row required):
#
#   db,path
#   bbv1,/data/bigbacter_v1_db
#   bbv2,s3://bucket/bigbacter_db
#
#   db     bbv1 or bbv2 -- the BigBacter version of that database, which
#          decides where in the species folder the reference FASTAs are
#          found. Each value may appear only once.
#   path   the BigBacter database root, either a local path (absolute,
#          relative, or starting with ~) or an s3:// path -- bbv1 and bbv2
#          can each be local or S3, in any combination. This is the folder
#          that CONTAINS the species folders. A trailing "/" or "/*" is
#          tolerated and ignored. Symlinks are followed for local paths.
#
# --species is required and names the species folder inside each BigBacter
# database, exactly as it appears there (case-sensitive), e.g.
#   --species Staphylococcus_aureus
# It applies to every row of the config.
#
# Blank lines and lines starting with '#' are skipped. Columns may be in
# either order, and extra columns are ignored.
#
# Requirements: bash, find, cp; the AWS CLI (`aws`) only if any path is s3://.
# S3 access uses your normal AWS credentials/profile (AWS_PROFILE etc.).

# ============================================================================
# BigBacter layout settings -- where reference FASTAs live in each version:
#
#   bbv1: <path>/<species>/clusters/<cluster_number>/assembly/*.fa.gz
#   bbv2: <path>/<species>/clusters/<cluster_number>/asm/*.fa.gz
#
# Every .fa.gz in those folders is copied (all of them, not one per cluster).
#
# *_REF_GLOB   shell glob a file's full path must match to be kept. For s3://
#              paths it is also used as the aws --include filter, so only the
#              matching files are downloaded.
# *_MAX_DEPTH  how many directory levels below the species folder to search
#              (clusters/<n>/assembly/<file> is 4 levels).
# ============================================================================
BBV1_REF_GLOB='*clusters/*/assembly/*.fa.gz'
BBV1_MAX_DEPTH=4

BBV2_REF_GLOB='*clusters/*/asm/*.fa.gz'
BBV2_MAX_DEPTH=4
# ============================================================================

DEST="cluster_references"
CONFIG=""
SPECIES=""
ONLY_V1=0
ONLY_V2=0
DRY_RUN=0
CLEAN=0
OVERWRITE=0

usage() {
  cat <<'EOF'
Usage: populate_bb_references.sh --config CONFIG.csv --species SPECIES [options]

Required:
  --config FILE       CSV with columns "db" and "path":
                        db    bbv1 or bbv2 (each value at most once)
                        path  BigBacter DB root, local or s3:// (the folder
                              that contains the species folders)
  --species NAME      Species folder inside the BigBacter database(s),
                      exactly as named there (case-sensitive),
                      e.g. Staphylococcus_aureus

Options:
  --bbv1              Only process the config row whose db is bbv1
  --bbv2              Only process the config row whose db is bbv2
                      (give neither, or both, to process every row)
  --dest DIR          Destination directory (default: cluster_references)
  --clean             Delete existing FASTA files in DEST before copying
  --overwrite         Overwrite files already in DEST that share a name
  --dry-run           Show what would be copied without copying anything
  -h, --help          Show this help

Example:
  populate_bb_references.sh --config bb_dbs.csv --species Staphylococcus_aureus

  bb_dbs.csv (each path can be local or s3://, in any combination):
    db,path
    bbv1,/data/bigbacter_v1_db
    bbv2,s3://bucket/bigbacter_db

Files are copied flat into DEST under their original file names. If two
sources have the same file name, the first is kept and the rest are reported
as collisions (they are not silently overwritten).

Reference locations searched (set in the settings block at the top of the script):
  bbv1  every .fa.gz in <path>/<species>/clusters/<n>/assembly/
  bbv2  every .fa.gz in <path>/<species>/clusters/<n>/asm/
EOF
}

need_val() { [ "$#" -ge 2 ] && [ -n "$2" ] || { echo "Error: $1 requires a value." >&2; exit 1; }; }

while [ "$#" -gt 0 ]; do
  case "$1" in
    --bbv1)      ONLY_V1=1; shift ;;
    --bbv2)      ONLY_V2=1; shift ;;
    --config)    need_val "$@"; CONFIG="$2"; shift 2 ;;
    --species)   need_val "$@"; SPECIES="$2"; shift 2 ;;
    --dest)      need_val "$@"; DEST="$2"; shift 2 ;;
    --clean)     CLEAN=1; shift ;;
    --overwrite) OVERWRITE=1; shift ;;
    --dry-run|--dryrun) DRY_RUN=1; shift ;;
    -h|--help)   usage; exit 0 ;;
    *) echo "Error: unknown argument '$1'. Run with --help for usage." >&2; exit 1 ;;
  esac
done

[ -n "$CONFIG" ]  || { echo "Error: --config is required." >&2; exit 1; }
[ -n "$SPECIES" ] || { echo "Error: --species is required (e.g. --species Staphylococcus_aureus)." >&2; exit 1; }
case "$SPECIES" in
  */*) echo "Error: --species must be a folder name, not a path (got '$SPECIES')." >&2; exit 1 ;;
esac
[ -f "$CONFIG" ] || { echo "Error: config file '$CONFIG' not found." >&2; exit 1; }

# If both (or neither) version flags are given, process every row.
if [ "$ONLY_V1" -eq 1 ] && [ "$ONLY_V2" -eq 1 ]; then ONLY_V1=0; ONLY_V2=0; fi

# ---------- read config ----------
trim() {
  local v="$1"
  v="${v#"${v%%[![:space:]]*}"}"
  v="${v%"${v##*[![:space:]]}"}"
  v="${v#\"}"; v="${v%\"}"
  printf '%s' "$v"
}
lower() { printf '%s' "$1" | tr 'A-Z' 'a-z'; }

ENTRY_DBS=()
ENTRY_ROOTS=()
n_rows_skipped=0
db_col=0; path_col=0; header_seen=0; lineno=0
bbv1_line=0; bbv2_line=0

while IFS= read -r line || [ -n "$line" ]; do
  lineno=$((lineno + 1))
  line="${line%$'\r'}"
  [ -z "${line//[[:space:]]/}" ] && continue
  case "$line" in \#*) continue ;; esac
  IFS=',' read -r -a cols <<< "$line"

  if [ "$header_seen" -eq 0 ]; then
    for i in "${!cols[@]}"; do
      name=$(lower "$(trim "${cols[$i]}")")
      if [ "$name" = "db" ]; then db_col=$((i + 1)); fi
      if [ "$name" = "path" ]; then path_col=$((i + 1)); fi
    done
    if [ "$db_col" -eq 0 ] || [ "$path_col" -eq 0 ]; then
      echo "Error: config header (line $lineno) must contain 'db' and 'path' columns." >&2
      echo "Expected something like:  db,path   (species is given with --species, not in the config)" >&2
      exit 1
    fi
    header_seen=1
    continue
  fi

  db=$(lower "$(trim "${cols[$((db_col - 1))]:-}")")
  root=$(trim "${cols[$((path_col - 1))]:-}")

  case "$db" in
    bbv1|bbv2) ;;
    *) echo "Error: config line $lineno: db must be 'bbv1' or 'bbv2' (got '$db')." >&2; exit 1 ;;
  esac
  [ -n "$root" ] || { echo "Error: config line $lineno: path is empty." >&2; exit 1; }

  # Expand a leading ~ (the shell does not do this for values read from a file).
  case "$root" in
    "~")    root="$HOME" ;;
    "~/"*)  root="$HOME/${root#\~/}" ;;
  esac

  # Tolerate a trailing "/*" or "/" on the database root.
  root="${root%/\*}"
  while [[ "$root" == */ ]] && [ "${#root}" -gt 1 ]; do root="${root%/}"; done

  # Each db value may appear only once in the config.
  if [ "$db" = "bbv1" ]; then
    [ "$bbv1_line" -eq 0 ] || { echo "Error: config line $lineno: db 'bbv1' is already used on line $bbv1_line; each db value must be unique." >&2; exit 1; }
    bbv1_line=$lineno
  else
    [ "$bbv2_line" -eq 0 ] || { echo "Error: config line $lineno: db 'bbv2' is already used on line $bbv2_line; each db value must be unique." >&2; exit 1; }
    bbv2_line=$lineno
  fi

  if { [ "$ONLY_V1" -eq 1 ] && [ "$db" != "bbv1" ]; } || { [ "$ONLY_V2" -eq 1 ] && [ "$db" != "bbv2" ]; }; then
    n_rows_skipped=$((n_rows_skipped + 1))
    continue
  fi

  ENTRY_DBS+=("$db")
  ENTRY_ROOTS+=("$root")
done < "$CONFIG"

[ "$header_seen" -eq 1 ] || { echo "Error: config '$CONFIG' is empty." >&2; exit 1; }
[ "${#ENTRY_ROOTS[@]}" -gt 0 ] || {
  echo "Error: no usable rows in '$CONFIG' (after --bbv1/--bbv2 filtering)." >&2
  exit 1
}

NEED_AWS=0
for e in "${ENTRY_ROOTS[@]}"; do
  case "$e" in s3://*) NEED_AWS=1 ;; esac
done
if [ "$NEED_AWS" -eq 1 ] && ! command -v aws >/dev/null 2>&1; then
  echo "Error: config contains s3:// paths but the AWS CLI ('aws') is not on PATH." >&2
  exit 1
fi

# ---------- helpers ----------
is_fasta() {
  local lower_name
  lower_name=$(lower "$1")
  case "$lower_name" in
    *.fa|*.fasta|*.fna|*.fas|*.fa.gz|*.fasta.gz|*.fna.gz|*.fas.gz) return 0 ;;
  esac
  return 1
}

STAGE_ROOT=$(mktemp -d)
SEEN_FILE="$STAGE_ROOT/seen.txt"
: > "$SEEN_FILE"
trap 'rm -rf "$STAGE_ROOT"' EXIT

n_copied=0
n_collisions=0
n_existing=0
n_filtered=0
row_matched=0

# consider_file SRC_PATH DISPLAY_PATH
#   SRC_PATH is the file to copy; DISPLAY_PATH is its original (local or s3)
#   location, used for the path filter and messages. Uses the REF_GLOB of the
#   config row currently being processed.
consider_file() {
  local src="$1" shown="$2" base
  base=$(basename "$src")

  if ! is_fasta "$base"; then
    n_filtered=$((n_filtered + 1)); return 0
  fi
  # shellcheck disable=SC2053  # REF_GLOB is intentionally a glob
  if [[ "$shown" != $REF_GLOB ]]; then
    n_filtered=$((n_filtered + 1)); return 0
  fi
  row_matched=$((row_matched + 1))

  if grep -Fxq -- "$base" "$SEEN_FILE"; then
    echo "  collision: '$base' already taken from another source; skipping $shown" >&2
    n_collisions=$((n_collisions + 1)); return 0
  fi
  if [ -e "$DEST/$base" ] && [ "$OVERWRITE" -eq 0 ] && [ "$CLEAN" -eq 0 ]; then
    echo "  exists: $DEST/$base (use --overwrite to replace); skipping $shown" >&2
    n_existing=$((n_existing + 1)); return 0
  fi

  printf '%s\n' "$base" >> "$SEEN_FILE"
  if [ "$DRY_RUN" -eq 1 ]; then
    echo "  would copy: $shown -> $DEST/$base"
  else
    cp -- "$src" "$DEST/$base"
    echo "  copied: $shown"
  fi
  n_copied=$((n_copied + 1))
}

# list_available LOCATION  -- print the folders under a database root, to help
# the user spot a mistyped species name.
list_available() {
  local root="$1"
  echo "  Folders found in $root/:" >&2
  case "$root" in
    s3://*) aws s3 ls "$root/" 2>/dev/null | sed 's/^ *//; s/^/    /' >&2 || true ;;
    *)      ls -1 "$root" 2>/dev/null | sed 's/^/    /' >&2 || true ;;
  esac
}

# ---------- prepare destination ----------
if [ "$DRY_RUN" -eq 0 ]; then
  mkdir -p "$DEST"
fi
if [ -d "$DEST" ]; then
  existing=$(find "$DEST" -maxdepth 1 -type f | wc -l)
  if [ "$CLEAN" -eq 1 ]; then
    if [ "$DRY_RUN" -eq 1 ]; then
      echo "Dry run: would remove existing FASTA files in $DEST/"
    else
      find "$DEST" -maxdepth 1 -type f \( -iname '*.fa' -o -iname '*.fasta' -o -iname '*.fna' -o -iname '*.fas' \
        -o -iname '*.fa.gz' -o -iname '*.fasta.gz' -o -iname '*.fna.gz' -o -iname '*.fas.gz' \) -delete
      echo "Cleaned existing FASTA files from $DEST/"
    fi
  elif [ "$existing" -gt 0 ]; then
    echo "Note: $DEST/ already contains $existing file(s); new files will be added alongside them." >&2
    echo "      make_ref_db.sh sketches everything in that directory. Use --clean for a fresh set." >&2
  fi
fi

echo "Config: $CONFIG | species: $SPECIES | ${#ENTRY_ROOTS[@]} row(s) to process, $n_rows_skipped skipped by version filter"
if [ "$DRY_RUN" -eq 1 ]; then echo "(dry run: nothing will be copied)"; fi

# ---------- process each config row ----------
for ((i = 0; i < ${#ENTRY_ROOTS[@]}; i++)); do
  idx=$((i + 1))
  db="${ENTRY_DBS[$i]}"
  root="${ENTRY_ROOTS[$i]}"
  row_matched=0

  if [ "$db" = "bbv1" ]; then
    MAX_DEPTH="$BBV1_MAX_DEPTH"; REF_GLOB="$BBV1_REF_GLOB"
  else
    MAX_DEPTH="$BBV2_MAX_DEPTH"; REF_GLOB="$BBV2_REF_GLOB"
  fi
  echo "[$idx/${#ENTRY_ROOTS[@]}] $db  species=$SPECIES  $root"

  case "$root" in
    s3://*)
      if [[ "$root" == *[\*\?\[]* ]]; then
        echo "  error: wildcards are not supported in s3:// paths ('$root')" >&2
        exit 1
      fi
      s3_base="$root/$SPECIES/"

      if [ "$DRY_RUN" -eq 1 ]; then
        listing=$(aws s3 cp "$s3_base" "$STAGE_ROOT/unused" --recursive --exclude "*" --include "$REF_GLOB" --dryrun || true)
        if [ -z "$listing" ]; then
          echo "  warning: no matching objects under $s3_base" >&2
          list_available "$root"
        else
          printf '%s\n' "$listing" | sed 's/^/  /'
        fi
        continue
      fi

      stage="$STAGE_ROOT/s3_$idx"
      mkdir -p "$stage"
      aws s3 cp "$s3_base" "$stage" --recursive --exclude "*" --include "$REF_GLOB" --only-show-errors
      while IFS= read -r -d '' f; do
        consider_file "$f" "${s3_base}${f#"$stage"/}"
      done < <(find "$stage" -maxdepth "$MAX_DEPTH" -type f -print0 | sort -z)
      rm -rf "$stage"
      ;;

    *)
      species_dir="$root/$SPECIES"
      if [ ! -d "$root" ]; then
        echo "  error: database path not found: $root" >&2
        exit 1
      fi
      if [ ! -d "$species_dir" ]; then
        echo "  error: species folder not found: $species_dir" >&2
        list_available "$root"
        exit 1
      fi
      while IFS= read -r -d '' f; do
        consider_file "$f" "$f"
      done < <(find -L "$species_dir" -maxdepth "$MAX_DEPTH" -type f -print0 | sort -z)
      ;;
  esac

  if [ "$row_matched" -eq 0 ]; then
    echo "  error: no reference files found for $db species '$SPECIES' under $root/" >&2
    echo "  Check the species name (case-sensitive), the db value, and the layout settings at the top of this script." >&2
    list_available "$root"
    exit 1
  fi
done

# ---------- summary ----------
echo
if [ "$DRY_RUN" -eq 1 ]; then
  echo "Dry run complete. Would copy $n_copied local file(s); S3 entries are listed above by aws."
else
  echo "Copied $n_copied reference file(s) into $DEST/"
fi
echo "Skipped: $n_filtered non-FASTA/filtered, $n_existing already present, $n_collisions name collision(s)."
