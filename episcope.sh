#!/usr/bin/env bash
set -euo pipefail

# Wrapper for the contextual-sequence pipeline (everything after the download):
#
#   1. make_ref_db.sh      sketch cluster_references/ -> ref_mash_db.msh
#   2. get_mash_dist.sh    assemblies vs. reference DB   -> ranked CSVs
#   3. long_to_matrix.py   (optional, --matrix)          -> distance matrix CSV
#   4. get_sra.py          closest assembly with SRA     -> info + summary CSVs
#
# With --tool sourmash, steps 1-2 use make_ref_db_sourmash.sh and
# get_sourmash_dist.py instead (ref_sourmash_db.zip). The output files and
# everything downstream are the same.
#
# Everything the pipeline generates is written to --outdir (created if it
# does not exist). Inputs (cluster_references/, assemblies/) are read from
# the current directory and are never modified.
#
# get-taxon-genomes.sh is NOT run here; run it first to populate assemblies/.
#
# Run this from the directory that holds cluster_references/ and
# assemblies/. The pipeline scripts are looked up next to this wrapper, so
# they don't need to be in the working directory.
#
# Requirements: everything in environment.yml (mash, sourmash, biopython, pandas, ...)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
  cat <<'EOF'
Usage: episcope.sh --outdir DIR --email EMAIL --db Genus_species [tech flags] [options]

Required:
  --outdir DIR           Directory for ALL generated output (created if it
                         does not exist)
  --email EMAIL          Email address required by NCBI Entrez
  --db Genus_species     Taxon label, e.g. Staphylococcus_aureus
  --illumina             Accept Illumina SRA runs  } at least one
  --ont                  Accept ONT SRA runs       } technology flag
  --pacbio               Accept PacBio SRA runs    } is required

Options:
  --tool mash|sourmash   Distance backend (default: mash). sourmash uses the
                         branchwater plugin and is multithreaded
  --threads N            Parallel jobs/threads for the distance step
                         (default: 4)
  --exclude IDS_OR_FILE  Never select these sequences. A comma/space-separated
                         list or a text file of IDs; repeat the flag to give
                         several. Accepts GenBank/RefSeq accessions (GCA_/GCF_,
                         with or without .version), SRA runs/experiments/
                         samples/studies (SRR/ERR/DRR, SRX.., SRS.., SRP..) and
                         BioSamples (SAMN/SAMEA/SAMD). Excluded genomes are
                         skipped and the next-closest one is used instead.
  --max-candidates N     Give up on a reference after N ranked candidates
  --report-file PATH     NCBI assembly report
                         (default: assemblies/assembly_data_report.jsonl)
  --matrix               Also build a distance matrix from unfiltered_mash.csv
  --matrix-out NAME      Matrix file name, written inside --outdir
                         (default: distance_matrix.csv)
  --skip-ref-db          Reuse the existing reference DB in OUTDIR
                         (ref_mash_db.msh or ref_sourmash_db.zip)
  --skip-mash            Reuse existing OUTDIR/ref_nonref_ranked.csv
  --skip-dist            (skips steps 1-2; works with either --tool)
  -h, --help             Show this help

The distance threshold (THRESHOLD, default 0.005) is set at the top of
get_mash_dist.sh, or via the THRESHOLD environment variable for sourmash.

Example:
  episcope.sh --outdir results --email me@example.org --db Staphylococcus_aureus --illumina --matrix
EOF
}

OUTDIR=""
TOOL="mash"
THREADS=4
EMAIL=""
DB=""
TECH_FLAGS=()
MAX_CANDIDATES=""
EXCLUDES=()
REPORT_FILE="assemblies/assembly_data_report.jsonl"
DO_MATRIX=0
MATRIX_OUT="distance_matrix.csv"
SKIP_REF_DB=0
SKIP_MASH=0

need_val() { [ "$#" -ge 2 ] && [ -n "$2" ] || { echo "Error: $1 requires a value." >&2; exit 1; }; }

while [ "$#" -gt 0 ]; do
  case "$1" in
    --outdir)         need_val "$@"; OUTDIR="$2"; shift 2 ;;
    --tool)           need_val "$@"; TOOL="$2"; shift 2 ;;
    --threads)        need_val "$@"; THREADS="$2"; shift 2 ;;
    --email)          need_val "$@"; EMAIL="$2"; shift 2 ;;
    --db)             need_val "$@"; DB="$2"; shift 2 ;;
    --illumina|--ont|--pacbio) TECH_FLAGS+=("$1"); shift ;;
    --max-candidates) need_val "$@"; MAX_CANDIDATES="$2"; shift 2 ;;
    --exclude)        need_val "$@"; EXCLUDES+=("$2"); shift 2 ;;
    --exclude=*)      EXCLUDES+=("${1#--exclude=}"); shift ;;
    --report-file)    need_val "$@"; REPORT_FILE="$2"; shift 2 ;;
    --matrix)         DO_MATRIX=1; shift ;;
    --matrix-out)     need_val "$@"; MATRIX_OUT="$2"; shift 2 ;;
    --skip-ref-db)    SKIP_REF_DB=1; shift ;;
    --skip-mash|--skip-dist) SKIP_MASH=1; shift ;;
    -h|--help)        usage; exit 0 ;;
    *) echo "Error: unknown argument '$1'" >&2; echo >&2; usage >&2; exit 1 ;;
  esac
done

# ---------- validate arguments ----------
[ -n "$OUTDIR" ] || { echo "Error: --outdir is required." >&2; exit 1; }
[ -n "$EMAIL" ] || { echo "Error: --email is required." >&2; exit 1; }
[ -n "$DB" ]    || { echo "Error: --db is required." >&2; exit 1; }
if [ "${#TECH_FLAGS[@]}" -eq 0 ]; then
  echo "Error: specify at least one of --illumina, --ont, --pacbio." >&2
  exit 1
fi

case "$TOOL" in
  mash|sourmash) ;;
  *) echo "Error: --tool must be 'mash' or 'sourmash' (got '$TOOL')." >&2; exit 1 ;;
esac
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "Error: --threads must be a positive integer (got '$THREADS')." >&2; exit 1; }

# Backend-specific scripts and files
if [ "$TOOL" = "sourmash" ]; then
  REF_DB_FILE="ref_sourmash_db.zip"
  DB_SCRIPT="make_ref_db_sourmash.sh"
  DIST_SCRIPT="get_sourmash_dist.py"
  DIST_CMD=(python3 "$SCRIPT_DIR/get_sourmash_dist.py")
else
  REF_DB_FILE="ref_mash_db.msh"
  DB_SCRIPT="make_ref_db.sh"
  DIST_SCRIPT="get_mash_dist.sh"
  DIST_CMD=(bash "$SCRIPT_DIR/get_mash_dist.sh")
fi

case "$MATRIX_OUT" in
  */*) echo "Error: --matrix-out must be a file name (it is written inside --outdir), got '$MATRIX_OUT'." >&2; exit 1 ;;
esac
if [ -e "$OUTDIR" ] && [ ! -d "$OUTDIR" ]; then
  echo "Error: --outdir '$OUTDIR' exists but is not a directory." >&2
  exit 1
fi

# ---------- validate tools ----------
required_tools=(python3)
if [ "$SKIP_MASH" -eq 0 ]; then
  if [ "$TOOL" = "sourmash" ]; then
    required_tools+=(sourmash)
  else
    required_tools+=(mash flock xargs find awk sort)
  fi
fi
for tool in "${required_tools[@]}"; do
  command -v "$tool" >/dev/null 2>&1 || {
    echo "Error: $tool not found on PATH (is the conda env active?)." >&2
    exit 1
  }
done
if [ "$SKIP_MASH" -eq 0 ] && [ "$TOOL" = "sourmash" ]; then
  sourmash scripts multisearch --help >/dev/null 2>&1 || {
    echo "Error: the sourmash branchwater plugin is not installed ('sourmash scripts' not found)." >&2
    echo "Install it with: conda install -c conda-forge sourmash_plugin_branchwater" >&2
    exit 1
  }
fi

for s in "$DB_SCRIPT" "$DIST_SCRIPT" long_to_matrix.py get_sra.py; do
  [ -f "$SCRIPT_DIR/$s" ] || { echo "Error: $SCRIPT_DIR/$s not found." >&2; exit 1; }
done

# The helper scripts must be the versions that write to $OUTDIR / --outdir.
# Older copies silently write to the current directory instead, which would
# leave results outside --outdir (and make later steps fail).
stale=()
if [ "$SKIP_MASH" -eq 0 ]; then
  if [ "$SKIP_REF_DB" -eq 0 ] && ! grep -q 'OUTDIR' "$SCRIPT_DIR/$DB_SCRIPT"; then
    stale+=("$DB_SCRIPT")
  fi
  grep -q 'THREADS' "$SCRIPT_DIR/$DIST_SCRIPT" || stale+=("$DIST_SCRIPT")
fi
grep -q -e '--outdir' "$SCRIPT_DIR/get_sra.py" || stale+=(get_sra.py)
if [ "${#EXCLUDES[@]}" -gt 0 ] && ! grep -q -e '--validate-exclude-only' "$SCRIPT_DIR/get_sra.py"; then
  [[ " ${stale[*]-} " == *" get_sra.py "* ]] || stale+=(get_sra.py)
fi
if [ "${#stale[@]}" -gt 0 ]; then
  echo "Error: these scripts in $SCRIPT_DIR are out of date (missing --outdir/--threads/--exclude support):" >&2
  for f in "${stale[@]}"; do echo "  $f" >&2; done
  echo "Replace them with the updated versions that came with this episcope.sh" >&2
  echo "(make_ref_db.sh, get_mash_dist.sh and get_sra.py all changed), then re-run." >&2
  exit 1
fi

# Check the --exclude IDs now, so a typo fails in seconds instead of after the
# long distance step.
if [ "${#EXCLUDES[@]}" -gt 0 ]; then
  EXCLUDE_ARGS=()
  for e in "${EXCLUDES[@]}"; do EXCLUDE_ARGS+=(--exclude "$e"); done
  python3 "$SCRIPT_DIR/get_sra.py" --validate-exclude-only "${EXCLUDE_ARGS[@]}" || exit 1
fi

# ---------- validate inputs ----------
if [ "$SKIP_MASH" -eq 0 ]; then
  if [ "$SKIP_REF_DB" -eq 0 ]; then
    [ -d cluster_references ] && [ -n "$(ls -A cluster_references 2>/dev/null)" ] || {
      echo "Error: cluster_references/ is missing or empty." >&2
      if [ -d bb_cluster_references ]; then
        echo "Found bb_cluster_references/ (the old default name). Rename it with:" >&2
        echo "  mv bb_cluster_references cluster_references" >&2
      fi
      exit 1
    }
  else
    [ -f "$OUTDIR/$REF_DB_FILE" ] || {
      echo "Error: --skip-ref-db given but $OUTDIR/$REF_DB_FILE does not exist (for --tool $TOOL)." >&2
      exit 1
    }
  fi
  [ -d assemblies ] && [ -n "$(ls -A assemblies 2>/dev/null)" ] || {
    echo "Error: assemblies/ is missing or empty. Run get-taxon-genomes.sh first." >&2
    exit 1
  }
else
  [ -f "$OUTDIR/ref_nonref_ranked.csv" ] || {
    echo "Error: --skip-mash given but $OUTDIR/ref_nonref_ranked.csv does not exist." >&2
    exit 1
  }
fi

[ -f "$REPORT_FILE" ] || {
  echo "Error: assembly report not found at '$REPORT_FILE'." >&2
  echo "It is created by get-taxon-genomes.sh, or generate it with:" >&2
  echo "  datasets summary genome taxon \"<taxon>\" --as-json-lines > $REPORT_FILE" >&2
  exit 1
}

if [ "$DO_MATRIX" -eq 1 ] && [ "$SKIP_MASH" -eq 1 ] && [ ! -f "$OUTDIR/unfiltered_mash.csv" ]; then
  echo "Error: --matrix with --skip-mash needs an existing $OUTDIR/unfiltered_mash.csv." >&2
  exit 1
fi

# ---------- run ----------
step() { printf '\n=== [%s] %s ===\n' "$(date +%H:%M:%S)" "$1"; }
START=$(date +%s)

# Create the output directory and tell the pipeline scripts where to write.
mkdir -p "$OUTDIR"
export OUTDIR THREADS
echo "Output directory: $OUTDIR | distance tool: $TOOL | threads: $THREADS"

if [ "$SKIP_MASH" -eq 0 ]; then
  if [ "$SKIP_REF_DB" -eq 0 ]; then
    step "1/4 Building reference DB ($TOOL)"
    bash "$SCRIPT_DIR/$DB_SCRIPT"
    if [ ! -f "$OUTDIR/$REF_DB_FILE" ]; then
      echo "Error: $OUTDIR/$REF_DB_FILE was not created." >&2
      if [ -f "$REF_DB_FILE" ]; then
        echo "A $REF_DB_FILE exists in the current directory instead, so $DB_SCRIPT" >&2
        echo "ignored OUTDIR. Make sure you are using the updated $DB_SCRIPT." >&2
      fi
      exit 1
    fi
  else
    step "1/4 Skipping reference DB build (reusing $OUTDIR/$REF_DB_FILE)"
  fi

  step "2/4 Computing distances and ranking candidates ($TOOL)"
  "${DIST_CMD[@]}"
  [ -f "$OUTDIR/ref_nonref_ranked.csv" ] || { echo "Error: $OUTDIR/ref_nonref_ranked.csv was not created." >&2; exit 1; }
else
  step "1-2/4 Skipping reference DB and distance steps (reusing $OUTDIR/ref_nonref_ranked.csv)"
fi

if [ "$DO_MATRIX" -eq 1 ]; then
  step "3/4 Building distance matrix -> $OUTDIR/$MATRIX_OUT"
  python3 "$SCRIPT_DIR/long_to_matrix.py" "$OUTDIR/unfiltered_mash.csv" "$OUTDIR/$MATRIX_OUT"
else
  step "3/4 Skipping distance matrix (use --matrix to build it)"
fi

step "4/4 Finding closest assemblies with SRA data"
SRA_ARGS=(--email "$EMAIL" --db "$DB" --report-file "$REPORT_FILE"
          --infile "$OUTDIR/ref_nonref_ranked.csv" --outdir "$OUTDIR" "${TECH_FLAGS[@]}")
[ -n "$MAX_CANDIDATES" ] && SRA_ARGS+=(--max-candidates "$MAX_CANDIDATES")
for e in ${EXCLUDES[@]+"${EXCLUDES[@]}"}; do SRA_ARGS+=(--exclude "$e"); done
python3 "$SCRIPT_DIR/get_sra.py" "${SRA_ARGS[@]}"

ELAPSED=$(( $(date +%s) - START ))
printf '\nPipeline finished in %dm%02ds. Output is in %s/\n' $((ELAPSED / 60)) $((ELAPSED % 60)) "$OUTDIR"
