# EpiScope: Contextual Sequence Selection for Epi

Scripts for finding publicly available, SRA-backed genomes that are closely related to a set of cluster reference genomes. The goal is to pull in **contextual sequences** (genomes near your outbreak/cluster references) to put local isolates in epidemiological context.

For each reference genome, the pipeline ranks every public assembly of the same taxon by Mash distance, then walks that ranked list and keeps the closest assembly that has a usable SRA run (matching the sequencing technology you ask for).

## Quickstart

Three commands take you from nothing to a list of SRA runs for contextual sequences. Run them all from one working directory; steps 1 and 2 create `assemblies/` and `cluster_references/` there, and step 3 reads them.

**0. Set up once**

```bash
conda env create -f environment.yml
conda activate episcope
chmod +x *.sh
```

**1. Download public assemblies for your taxon** (`get-taxon-genomes.sh`)

```bash
./get-taxon-genomes.sh "Staphylococcus aureus" --assembly-level complete
```

To download only genomes collected in certain places, add `--geo-location` (the filter is applied *before* downloading, so only matching genomes are fetched), and preview with `--dry-run` first:

```bash
./get-taxon-genomes.sh "Staphylococcus aureus" --assembly-level complete --geo-location "USA*" --dry-run
./get-taxon-genomes.sh "Staphylococcus aureus" --assembly-level complete --geo-location "USA*"
```

This writes the genomes, plus NCBI's `assembly_data_report.jsonl` (needed later for the SRA lookup), to `assemblies/`. It checks that you have enough disk space before downloading anything. Anything you'd normally pass to `datasets download genome taxon` works as a filter.

**1b. (Optional) Filter genomes you have already downloaded** (`filter_assemblies_by_geo.sh`)

If you downloaded without a location filter, you can still remove the genomes you don't want afterwards:

```bash
./filter_assemblies_by_geo.sh --geo-location "USA*" --dry-run   # preview what would be removed
./filter_assemblies_by_geo.sh --geo-location "USA*"
```

Genomes whose BioSample location doesn't match are moved out of `assemblies/` into `filtered_out.tar.gz` (nothing is deleted; see [1c](#1c-optional-filter-already-downloaded-assemblies-by-geographic-location)). If you used `--geo-location` in step 1 you don't need this.

**2. Collect your cluster reference genomes** (`populate_bb_references.sh`)

If your references are in a BigBacter database (local or on S3), list the database(s) in a small CSV:

```csv
db,path
bbv1,/data/bigbacter_v1_db
bbv2,s3://my-bucket/bigbacter_db
```

then copy the references for one species into `cluster_references/`:

```bash
./populate_bb_references.sh --config bb_dbs.csv --species Staphylococcus_aureus --dry-run   # preview
./populate_bb_references.sh --config bb_dbs.csv --species Staphylococcus_aureus
```

`--species` is the species folder name exactly as it appears inside the BigBacter database. If your references aren't in a BigBacter database, skip this step and put your reference FASTA files straight into `cluster_references/`.

**3. Find the closest assembly with SRA data for each reference** (`episcope.sh`)

```bash
./episcope.sh --outdir results --email you@example.org --db Staphylococcus_aureus --illumina
```

`--outdir` is where everything is written (created if needed), `--email` is required by NCBI, `--db` is a `Genus_species` label used in output file names, and at least one of `--illumina`, `--ont` or `--pacbio` selects the sequencing technology. Useful extras: `--exclude IDS_OR_FILE` to keep specific genomes or runs out of the selection (see [Excluding specific sequences](#excluding-specific-sequences)), `--tool sourmash --threads 16` for the multithreaded backend (see [Mash or sourmash](#mash-or-sourmash)) and `--matrix` for a distance matrix.

**Result:** `results/BigBacter.<db>.closest_to_SRA_<date>.csv` lists the selected assemblies (`sample, taxa, genbank, sra`), ready for downloading reads. `results/BigBacter.<db>.info.csv` adds which reference each came from, its distance and its rank. The sections below explain every step, option and output in detail.

## Pipeline overview

```
get-taxon-genomes.sh      download all public assemblies for a taxon      -> assemblies/
                          (optionally only those from given locations)
filter_assemblies_by_geo.sh (optional) filter ALREADY-downloaded genomes    -> filtered_out.tar.gz
populate_bb_references.sh (optional) copy references from a BigBacter DB  -> cluster_references/
make_ref_db.sh            sketch cluster references into a Mash DB        -> ref_mash_db.msh
get_mash_dist.sh          assemblies vs. reference DB, rank per reference -> CSVs (see below)
  (or, with --tool sourmash:)
make_ref_db_sourmash.sh   sketch cluster references into a sourmash DB    -> ref_sourmash_db.zip
get_sourmash_dist.py      same outputs as get_mash_dist.sh, multithreaded -> CSVs (see below)

long_to_matrix.py         (optional) pivot distances to a matrix          -> matrix CSV
get_sra.py                closest assembly per reference that has SRA     -> info + summary CSVs

episcope.sh           wrapper that runs the last four steps in one go -> --outdir
```

You can run each step by hand (see [Step-by-step](#step-by-step)) or use the wrapper (see [Running everything with the wrapper](#running-everything-with-the-wrapper)). Either way, run `get-taxon-genomes.sh` first.

## Requirements

The easiest setup is the provided conda environment, which installs everything listed below:

```bash
conda env create -f environment.yml
conda activate contextual-seqs
```

Or install the tools yourself:

| Tool | Used by |
|------|---------|
| NCBI `datasets`, `dataformat` | `get-taxon-genomes.sh` |
| `jq` (1.6 or newer for `--geo-location`), `unzip`, `split`, `df`, `awk` | `get-taxon-genomes.sh` |
| `jq` (1.6 or newer), `tar`, `awk` | `filter_assemblies_by_geo.sh` |
| `mash` | `make_ref_db.sh`, `get_mash_dist.sh` (the default backend) |
| `sourmash`, `sourmash_plugin_branchwater` | `make_ref_db_sourmash.sh`, `get_sourmash_dist.py` (optional backend, see [Mash or sourmash](#mash-or-sourmash)) |
| `flock`, `xargs`, `find`, `sort` | `get_mash_dist.sh` |
| Python 3 with `pandas` | `long_to_matrix.py` |
| Python 3 with `biopython` | `get_sra.py` |

> **Renamed in this version:** the reference directory is now `cluster_references/` (previously `bb_cluster_references/`), and the reference list written next to the Mash DB is `cluster_ref_list.txt` (previously `bb_ref_list.txt`). If you have an existing directory, rename it with `mv bb_cluster_references cluster_references`; the wrapper reminds you if it finds the old name. `long_to_matrix.py` still accepts the old prefix in older distance tables.

## Expected directory layout

```
.
├── cluster_references/     # your reference genomes (FASTA), input to make_ref_db.sh
├── assemblies/                # created by get-taxon-genomes.sh
│   ├── *.fna                  # genomes
│   └── assembly_data_report.jsonl   # NCBI metadata, needed by get_sra.py
└── *.sh / *.py
```

## Step-by-step

### 1. Download public assemblies

```bash
./get-taxon-genomes.sh "Staphylococcus aureus" --assembly-level complete
```

Pass whatever you would normally pass to `datasets download genome taxon`. Behavior:

- Without `--reference`, files go to `assemblies/`. With `--reference`, they go to `reference_genomes/` instead.
- Downloads are batched (default 40,000 genomes per batch) to stay under NCBI's ~50k per-download ceiling.
- Before downloading, it estimates required disk space from NCBI metadata and refuses to start if the filesystem can't hold it. It re-checks free space before each batch.
- `--include ...` flags are forwarded to the download step only.
- Each batch is downloaded and extracted in a temporary `.ncbi_batch_*` directory (removed afterwards, including on failure), so NCBI's bundled `README.md` can't overwrite or delete files in your working directory.

#### Downloading only genomes from certain locations

```bash
./get-taxon-genomes.sh "Haemophilus influenzae" --geo-location "USA*" --dry-run
./get-taxon-genomes.sh "Haemophilus influenzae" --assembly-level complete --geo-location "USA*"
./get-taxon-genomes.sh "Escherichia coli" --geo-location Canada --geo-location "USA: Oregon" --accessions-out accessions.txt
```

| Option | Meaning |
|--------|---------|
| `--geo-location TEXT` (or `--geolocation`) | Only download genomes whose BioSample location matches `TEXT`. Repeat the flag to accept several locations (a genome matching any is kept) |
| `--dry-run` | Resolve the genome list, apply the filters and run the disk space check, then stop. Nothing is downloaded and no directory is created |
| `--accessions-out FILE` | Write the final accession list (one per line) to `FILE`, e.g. to reuse with `datasets download genome accession --inputfile FILE` |

Everything else on the command line is passed to `datasets` exactly as before (filters such as `--assembly-level`, `--annotated`, `--reference` go to `datasets summary`, and `--include` goes to the download).

**How it works.** The script asks `datasets summary` for every genome matching your query, keeps the accessions whose location matches, and then downloads only those, in batches. This is the same idea as

```bash
datasets summary genome taxon "Haemophilus influenzae" --as-json-lines | grep -i '"geo_loc_name":"USA' | jq -r '.accession' > accessions.txt
```

followed by `datasets download genome accession --inputfile accessions.txt`, but with more robust matching. That `grep` only finds records that happen to store the location as a flat `"geo_loc_name":"USA…"` key at the start of the value. The script finds the location wherever it is stored in the record (a `geo_loc_name` attribute, a `geo_loc_name` key, or a `geoLocName` key), and it supports the wildcard, case and spacing rules below. Only matching genomes are downloaded and only they count toward the disk space check.

The matching rules, and what happens to genomes with no location, are the same as for [`filter_assemblies_by_geo.sh`](#1c-optional-filter-already-downloaded-assemblies-by-geographic-location) (described there). In short: a case-insensitive substring match (`USA` also matches `USA: Washington`), `*` matches any characters, spacing around the colon is ignored, and genomes with no location are excluded.

If `--geo-location` selects nothing, the script prints the locations that are present, so you can see how they are spelled. If **no** record in the summary contains a location at all, it stops with an error (exit 1) and suggests downloading without the filter and then using `filter_assemblies_by_geo.sh`, which works from the report that comes with the download. It needs **jq 1.6 or newer** for `--geo-location`.

Environment overrides:

| Variable | Default | Meaning |
|----------|---------|---------|
| `BATCH_SIZE` | `40000` | Genomes per download batch |
| `SPACE_SAFETY_MARGIN` | `1.2` | Multiplier on estimated space needed |
| `SKIP_SPACE_CHECK` | `0` | Set to `1` to bypass the disk check |

```bash
BATCH_SIZE=25000 SPACE_SAFETY_MARGIN=1.5 ./get-taxon-genomes.sh "Escherichia coli"
```

### 1b. (Optional) Populate `cluster_references/` from BigBacter databases

If your references live in BigBacter v1 or v2 databases, `populate_bb_references.sh` copies them into `cluster_references/` for you. Each database can be a **local path or on S3**, and the two versions can be mixed (e.g. `bbv1` local and `bbv2` on S3).

```bash
./populate_bb_references.sh --config bb_dbs.csv --species Staphylococcus_aureus
```

The config is a CSV with two columns, `db` and `path` (the header row is required):

```csv
db,path
bbv1,/data/bigbacter_v1_db
bbv2,s3://bucket/bigbacter_db
```

Either row could equally be local or `s3://`; for example, both local:

```csv
db,path
bbv1,/data/bigbacter_v1_db
bbv2,/data/bigbacter_v2_db
```

| Column | Meaning |
|--------|---------|
| `db` | `bbv1` or `bbv2`, the BigBacter version. Each value may appear **only once**, so there is at most one row per version (a duplicate stops the run with an error) |
| `path` | The BigBacter database root, i.e. the folder that *contains* the species folders. Either a local path (absolute, relative to where you run the script, or starting with `~`) or an `s3://` path. A trailing `/` or `/*` is ignored |

Columns can be in any order, extra columns are ignored, and blank lines and `#` comment lines are skipped. Paths containing a comma aren't supported, and `s3://` paths can't contain wildcards. For local paths, spaces in paths are fine and symbolic links are followed (a symlinked database root, species folder, cluster folder or reference file all work).

The species is **not** part of the config. It is a required argument and applies to every row.

| Option | Meaning |
|--------|---------|
| `--config FILE` | Required; the CSV described above |
| `--species NAME` | **Required.** The species folder inside the BigBacter database(s), written exactly as it appears there (case-sensitive), e.g. `Staphylococcus_aureus`. Must be a folder name, not a path |
| `--bbv1` / `--bbv2` | Optional; only process the row with that `db` value (give neither to process every row) |
| `--dest DIR` | Destination directory (default `cluster_references`) |
| `--clean` | Delete existing FASTA files in the destination first, so stale references don't end up in the Mash DB |
| `--overwrite` | Replace files in the destination that have the same name |
| `--dry-run` | Show what would be copied without copying |

Behavior:

- The `db` value decides where inside the species folder the references are found. **Every** `.fa.gz` file in these folders is copied (not just one per cluster):

  | `db` | Files copied |
  |------|--------------|
  | `bbv1` | `<path>/<species>/clusters/<cluster_number>/assembly/*.fa.gz` |
  | `bbv2` | `<path>/<species>/clusters/<cluster_number>/asm/*.fa.gz` |

  Files are flattened into the destination under their original names (so v1 files keep any `_T1` suffix). Other FASTAs in the database, and files in those folders that aren't `.fa.gz`, are ignored.
- Only the named species is read. On S3, only that species' reference files are downloaded (into a temporary directory that is removed afterwards), not the whole database. `--dry-run` lists the matching S3 objects without downloading.
- If a row finds no reference files (a mistyped species, the wrong `db` value for that database, or a different layout), the script stops with an error and lists the folders that exist under the database root so you can spot the right species name.
- If two sources contain a file with the same name, the first is kept and the others are reported as collisions, never silently overwritten.
- Files already in the destination are left alone (and a note is printed) unless you use `--clean` or `--overwrite`. `make_ref_db.sh` sketches **everything** in `cluster_references/`, so use `--clean` when switching reference sets.
- S3 access uses the AWS CLI and your normal credentials (`AWS_PROFILE` etc.).

**Layout settings:** the reference locations above are set in the settings block at the top of the script (`BBV1_REF_GLOB`, `BBV2_REF_GLOB`, and `*_MAX_DEPTH`). Edit them if your databases are organized differently. Check with `--dry-run` first.

### 1c. (Optional) Filter already-downloaded assemblies by geographic location

```bash
./filter_assemblies_by_geo.sh --geo-location "USA*" --dry-run   # preview first
./filter_assemblies_by_geo.sh --geo-location "USA*"
./filter_assemblies_by_geo.sh --geo-location Canada --geo-location "USA: Oregon"
```

Use this when you have **already downloaded** genomes without a location filter. (If you know the location up front, `--geo-location` in step 1 is more efficient, since it avoids downloading the genomes you don't want.) Run it before the pipeline, since the pipeline compares everything left in `assemblies/`. The script reads `assemblies/assembly_data_report.jsonl`, keeps the genomes whose BioSample `geo_loc_name` matches, and **moves the rest into `filtered_out.tar.gz`**. To get them back: `tar -xzf filtered_out.tar.gz -C assemblies`.

| Option | Meaning |
|--------|---------|
| `--geo-location TEXT` (or `--geolocation`) | **Required.** Keep genomes whose location matches `TEXT`. Repeat the flag to accept several locations (a genome matching any is kept) |
| `--dry-run` | Show what would be moved and kept, with each file's location, and change nothing |
| `--assembly-dir DIR` | Directory of downloaded genomes (default `assemblies`) |
| `--report FILE` | NCBI report (default `<assembly-dir>/assembly_data_report.jsonl`) |
| `--tar FILE` | Archive for the removed genomes (default `filtered_out.tar.gz`; parent directories are created) |
| `--allow-empty` | Allow a filter that would remove every genome |

How matching works:

- It is a **case-insensitive substring** match, so `USA` also keeps `USA: Washington`, `USA:California` and plain `USA`. Spacing around the colon is ignored, so `usa:washington`, `USA: Washington` and `USA:Washington` behave the same.
- A `*` inside `TEXT` matches any run of characters (`"USA*Seattle"` matches `USA: Washington, Seattle`). A leading or trailing `*` changes nothing, and every other character, including `.`, `(`, `)` and `[`, is matched literally. Quote patterns containing `*` so your shell doesn't expand them.
- Submitters enter locations as free text, usually `Country: Region, Locality`, so a broad term can keep more than you intend: `"USA: Washington"` also matches `USA: Washington DC`. Make the text as specific as the data allows, and check with `--dry-run`, which also shows how many genomes sit at each location.
- Genomes whose BioSample has **no** `geo_loc_name` don't match, so they are moved to the archive.
- The location is found wherever the record stores it: a `geo_loc_name` attribute, a `geo_loc_name` key, or a `geoLocName` key.

Which files are touched:

- Only files whose names start with a GenBank/RefSeq accession (`GCA_…` / `GCF_…`) **and** that have a record in the report, matched by either the record's accession or its paired accession. Any other file with that accession in its name (e.g. a gzipped copy) moves with it.
- Everything else is left alone: the report itself, `dataset_catalog.json`, and **your own sequences** with other names, since they have no location to check.
- Genome files that have no record in the report are left alone too, and the script counts them and warns. That typically means the report covers only part of the download (for example a download split into several batches overwrote it; see [the quirks below](#known-quirks-and-things-to-watch)).
- The report file is not modified.

Safety:

- The archive is written and verified (readable, correct number of files) **before** anything is deleted.
- It never overwrites an existing archive: if `filtered_out.tar.gz` already exists and there is something to move, it stops and asks you to pick another name with `--tar` or move the old one away. Re-running when nothing needs moving does nothing.
- If the filter would remove every genome (usually a mistyped location), it refuses and lists the locations that are present, unless you pass `--allow-empty`.
- It needs **jq 1.6 or newer** and stops with a clear message otherwise.

### 2. Build the reference Mash DB

Put your cluster reference FASTAs in `cluster_references/` (by hand or with step 1b), then:

```bash
./make_ref_db.sh
```

This writes `cluster_ref_list.txt` and `ref_mash_db.msh`, and deletes any existing `ref_mash_db*` files first. By default they go in the current directory; set `OUTDIR` to write them elsewhere (created if needed), e.g. `OUTDIR=results ./make_ref_db.sh`.

> **Using sourmash instead?** Steps 2 and 3 have sourmash equivalents (`make_ref_db_sourmash.sh` and `get_sourmash_dist.py`); see [Mash or sourmash](#mash-or-sourmash). Everything from step 4 on is identical.

### 3. Compute Mash distances and rank candidates

```bash
./get_mash_dist.sh
```

Every file in `assemblies/` is sketched and compared to `ref_mash_db.msh`, 4 jobs in parallel. Settings are variables at the top of the script (`ASSEMBLY_DIR`, `THRESHOLD`, `PARALLEL`).

Run standalone, it reads the DB from and writes the CSVs below to the current directory. Set `OUTDIR` to use another directory for both (the DB is expected at `$OUTDIR/ref_mash_db.msh`), e.g. `OUTDIR=results ./get_mash_dist.sh`.

| Output | Contents |
|--------|----------|
| `unfiltered_mash.csv` | `assembly_file,reference,mash_distance` for every assembly/reference pair |
| `passed_dist_filter.csv` | Assemblies within `THRESHOLD` (default 0.005) of at least one reference |
| `ref_nonref_ranked.csv` | `reference,rank,assembly,mash_distance`, all non-reference assemblies per reference, closest first |
| `ref_nearest_nonref.csv` | Rank-1 row per reference (kept for backward compatibility) |

Reference names in the ranked/nearest files are cleaned: bracketed text removed, `origname=` unwrapped, path stripped to basename, and stray quotes/commas removed.

### 4. (Optional) Distance matrix

```bash
python long_to_matrix.py unfiltered_mash.csv distance_matrix.csv
```

Pivots the long-form table into an assembly × reference matrix. Input and output paths are whatever you pass on the command line. The `cluster_references/` prefix is stripped from reference names. If a pair appears more than once, the first value is used. If your reference directory has a different name, edit `PREFIX` in the script.

### 5. Find the closest assembly with SRA data

```bash
python get_sra.py --email you@example.org --db Staphylococcus_aureus --illumina
```

For each reference, candidates are tried in rank order. A candidate is accepted if it is in the NCBI assembly report with a matching sequencing technology, has a BioSample, and that BioSample has SRA runs of the requested technology. References where no candidate qualifies are skipped and reported.

| Option | Meaning |
|--------|---------|
| `--email` | Required by Entrez |
| `--db` | Taxon label in `<Genus_species>` format, used in output file names and the `taxa` column |
| `--illumina`, `--ont`, `--pacbio` | Sequencing technologies to accept. **At least one is required** |
| `--infile` | Ranked CSV (default `ref_nonref_ranked.csv`) |
| `--max-candidates` | Give up on a reference after this many ranked candidates |
| `--exclude IDS_OR_FILE` | Never select these sequences (see [Excluding specific sequences](#excluding-specific-sequences)). A comma/space-separated list or a text file of IDs; repeat the flag to give several |
| `--validate-exclude-only` | Check the `--exclude` values for unrecognised IDs, print a summary and exit without reading any input |
| `--outdir` | Directory for the output CSVs, created if needed (default: current directory). Applies to the default output names below |
| `--outfile` | Full info CSV path (default `<outdir>/BigBacter.<db>.info.csv`) |
| `--summary-outfile` | Summary CSV path (default `<outdir>/BigBacter.<db>.closest_to_SRA_<YYYYMMDD>.csv`) |
| `--report-file` | NCBI report (default `assemblies/assembly_data_report.jsonl`) |

Outputs:

- `BigBacter.<db>.info.csv`: `reference,sample,taxa,genbank,sra,mash_distance,candidate_rank`, one row per SRA run.
- `BigBacter.<db>.closest_to_SRA_<date>.csv`: `sample,taxa,genbank,sra`, de-duplicated, ready to feed into downstream SRA download/analysis.

## Running everything with the wrapper

`episcope.sh` runs steps 2-5 (reference DB, Mash distances, optional matrix, SRA lookup) in sequence. It does **not** run `get-taxon-genomes.sh`; do step 1 first.

```bash
conda activate contextual-seqs
./episcope.sh --outdir results --email you@example.org --db Staphylococcus_aureus --illumina
```

Run it from the directory containing `cluster_references/` and `assemblies/`. The pipeline scripts are found next to the wrapper, so they don't have to be in the working directory. If the wrapper isn't executable, run `chmod +x episcope.sh` or call it with `bash episcope.sh`.

| Option | Meaning |
|--------|---------|
| `--outdir DIR` | **Required.** Directory for everything the pipeline generates. Created (including parent directories) if it doesn't exist |
| `--tool mash\|sourmash` | Distance backend (default `mash`); see [Mash or sourmash](#mash-or-sourmash) |
| `--threads N` | Parallel assemblies (mash) or threads (sourmash) for the distance step (default 4) |
| `--email`, `--db` | Required; passed to `get_sra.py` |
| `--illumina`, `--ont`, `--pacbio` | At least one required; passed to `get_sra.py` |
| `--max-candidates N` | Passed to `get_sra.py` |
| `--exclude IDS_OR_FILE` | Never select these sequences: GenBank/RefSeq accessions, SRA IDs or BioSamples, as a list or a file; repeat to give several. Checked before the pipeline starts. See [Excluding specific sequences](#excluding-specific-sequences) |
| `--report-file PATH` | Passed to `get_sra.py` (default `assemblies/assembly_data_report.jsonl`) |
| `--matrix` | Also run `long_to_matrix.py` (off by default, since the matrix can be very large for big taxa) |
| `--matrix-out NAME` | Matrix file name, written inside `--outdir` (default `distance_matrix.csv`; must not contain a `/`) |
| `--skip-ref-db` | Reuse the existing reference DB in the outdir (`ref_mash_db.msh` or `ref_sourmash_db.zip`, matching `--tool`) instead of rebuilding it |
| `--skip-dist` (alias `--skip-mash`) | Skip the DB and distance steps and reuse an existing `OUTDIR/ref_nonref_ranked.csv`, e.g. to rerun only the SRA lookup with a different technology or `--max-candidates` |

**Outputs.** Everything below is written to `--outdir`; nothing is written to the working directory, and the inputs (`cluster_references/`, `assemblies/`) are never modified.

| File | Step |
|------|------|
| `cluster_ref_list.txt`, `ref_mash_db.msh` (or `ref_sourmash_db.zip` with `--tool sourmash`) | reference DB |
| `unfiltered_mash.csv`, `passed_dist_filter.csv`, `ref_nonref_ranked.csv`, `ref_nearest_nonref.csv` | Mash distances |
| `distance_matrix.csv` (or `--matrix-out`) | only with `--matrix` |
| `BigBacter.<db>.info.csv`, `BigBacter.<db>.closest_to_SRA_<date>.csv` | SRA lookup |

Use a different `--outdir` per run (e.g. per species) to keep results separate. Re-using an outdir overwrites the files above.

Before running anything, the wrapper checks that the required tools are on `PATH`, that the input directories and report file exist, that `--outdir` is usable, and that at least one technology flag was given, so it fails fast rather than partway through a long Mash run. It stops at the first failing step.

The distance threshold is not a wrapper option. Edit `THRESHOLD` at the top of `get_mash_dist.sh`, or set the `THRESHOLD` environment variable when using `--tool sourmash`.

## Excluding specific sequences

Use `--exclude` to make sure certain genomes or sequencing runs are never chosen as contextual sequences (for example known outliers, contaminated samples, or sequences you already have).

```bash
./episcope.sh --outdir results --email you@example.org --db Staphylococcus_aureus --illumina \
    --exclude GCF_000783955.2,SRR1234567
./episcope.sh ... --exclude exclude_list.txt                       # IDs from a file
./episcope.sh ... --exclude exclude_list.txt --exclude SAMN0000001   # both
```

The value is either a comma/space-separated list of IDs or the path to a text file (one or more IDs per line, commas allowed, `#` starts a comment). Repeat the flag to combine lists. These IDs are accepted, in any mix and in upper or lower case:

| Kind of ID | Examples | What is excluded |
|------------|----------|------------------|
| GenBank / RefSeq assembly | `GCA_000783955.2`, `GCF_000783955.2`, `GCF_000783955` | That genome. A GCA_/GCF_ ID also excludes its paired accession (excluding `GCA_…` removes the `GCF_…` genome and vice versa, using the report). With a `.version` only that version is excluded; without one, every version is |
| SRA run | `SRR…`, `ERR…`, `DRR…` | That run only. If the genome has other runs it is still used with the remaining ones |
| SRA experiment / sample / study | `SRX…`/`ERX…`/`DRX…`, `SRS…`/`ERS…`/`DRS…`, `SRP…`/`ERP…`/`DRP…` | Every run belonging to it |
| BioSample | `SAMN…`, `SAMEA…`, `SAMD…` | The genome with that BioSample |

How it works: for each reference, the SRA step walks the ranked candidates from closest to farthest and takes the first one with usable SRA data. An excluded genome is simply skipped, and a genome whose runs are all excluded is rejected, so the **next-closest** genome is considered instead. If every candidate for a reference is excluded, that reference is skipped and the message shows the `excluded=` count.

Good to know:

- **IDs are checked before anything runs.** An unrecognised ID (a typo, or a file name that doesn't exist) stops the wrapper straight away, listing the offending value and its file and line number, rather than failing at the end of a long run. A valid list is summarised first (for example `Exclusions OK: 2 assembly accession(s), 1 SRA ID(s), 0 BioSample(s)`).
- **Exclusion happens at selection time.** Excluded genomes still appear in the distance tables and the matrix; they just can't be selected. This means you can change your exclusions without recomputing distances: re-run with `--skip-dist` and a new `--exclude`.
- **The run ends with a report:** how many candidates were skipped and runs dropped, and any excluded IDs that were never matched. An ID is only matched if its genome was reached while walking a reference's ranked candidates, so IDs for genomes that rank behind a selected one (or that aren't in your data) show as "not encountered". That is normal and not an error.
- It is implemented in `get_sra.py` (`--exclude`), so you can use it directly too.

## Mash or sourmash

The distance step (reference DB + assemblies-vs-references) can run on either backend. Everything else, including the output file names, is the same, so the rest of the pipeline doesn't care which one you used.

```bash
./episcope.sh --outdir results --tool sourmash --threads 16 \
    --email you@example.org --db Staphylococcus_aureus --illumina
```

| | `--tool mash` (default) | `--tool sourmash` |
|---|---|---|
| Scripts | `make_ref_db.sh`, `get_mash_dist.sh` | `make_ref_db_sourmash.sh`, `get_sourmash_dist.py` |
| Reference DB | `ref_mash_db.msh` | `ref_sourmash_db.zip` |
| Parallelism | `--threads N` assemblies at a time, one thread each | `--threads N` threads inside sketching and comparison |
| Ranking step | serial `awk` scan per reference | in memory (Python) |
| Extra requirements | `mash` | `sourmash` and the branchwater plugin |

Things to know when using sourmash:

- **Install the plugin.** `conda install -c conda-forge sourmash_plugin_branchwater` (included in `environment.yml`). The wrapper stops with this message if `sourmash scripts` isn't available.
- **Distances are estimates, not Mash's numbers.** They are Mash-style distances computed from the FracMinHash Jaccard that sourmash reports (`D = -ln(2j/(1+j))/k`, k = 21). They rank genomes the same way but won't match Mash's values exactly, so the `0.005` threshold may select a slightly different set. If you move an existing analysis to sourmash, re-check the threshold.
- **Resolution is set by `scaled`.** Both scripts default to `SCALED=1000`, comparable to Mash's default sketch size. For very close genomes (outbreak-level distances) a lower value, such as `SCALED=100`, gives finer resolution at the cost of time and memory. Set it when building the DB (`SCALED=100 ./make_ref_db_sourmash.sh`); the distance step uses the DB's value unless you set a higher one. K is set with `KSIZE` (default 21) and must match between DB and comparison.
- **Build the DB with the provided script.** It uses `manysketch` so each reference sketch is named after its file (e.g. `refA.fa`). A DB made with plain `sourmash sketch dna` has unnamed sketches, which the distance script rejects with an explanation.
- **Non-FASTA files are ignored.** Only `.fa`, `.fasta`, `.fna`, `.fas` (optionally `.gz`) in `assemblies/` are sketched, so NCBI's `assembly_data_report.jsonl` and `dataset_catalog.json` don't cause errors. Assemblies that fail to sketch are reported and excluded.
- **Names are the same as with Mash.** References are labelled by file name and the unfiltered table is still called `unfiltered_mash.csv`, so `long_to_matrix.py` and `get_sra.py` work unchanged.
- **Reuse flags work for either backend.** `--skip-ref-db` looks for the DB that matches `--tool`; `--skip-dist` (alias `--skip-mash`) reuses `ref_nonref_ranked.csv` whichever backend made it.
- **Memory.** `multisearch` holds the sketches in memory, and the unfiltered table has one row per assembly × reference pair, as it does with Mash.

I checked the two backends against each other on synthetic genomes: they produced the same reference ranking and the same selected assemblies, with distances within a few ten-thousandths of each other. That's a small test, so compare them on a subset of your own data before relying on sourmash for results.

`--threads` also applies to `--tool mash`, where it sets how many assemblies are processed in parallel (default 4).

## Known quirks and things to watch

- **`--exclude` doesn't remove genomes from the distance tables.** It only prevents them from being selected, so excluded genomes still appear in `unfiltered_mash.csv`, the ranked lists and the matrix. To remove genomes from everything, take them out of `assemblies/` before running (for example with `filter_assemblies_by_geo.sh`, or by deleting the files).
- **Multi-batch downloads overwrite the metadata report.** `get-taxon-genomes.sh` moves every file from each batch into one directory, so `assembly_data_report.jsonl` is replaced by each batch. With more than `BATCH_SIZE` genomes, the report will only describe the last batch, and `get_sra.py` will reject candidates from earlier batches as "not in tech map". For large taxa, regenerate a full report, e.g. `datasets summary genome taxon "<taxon>" --as-json-lines > assemblies/assembly_data_report.jsonl` (use the same filters as the download).
- **Non-FASTA files in `assemblies/` (Mash backend).** `get_mash_dist.sh` processes every file in the directory, including `assembly_data_report.jsonl` and `dataset_catalog.json`. These fail to sketch silently and produce no rows, but they count toward the progress total.
- **Mash parameters must match.** The reference DB and assembly sketches both use Mash defaults (k=21, s=1000). If you change one, change both.
- **Mash resolution at small distances.** With default sketch size, distances near 0.005 are coarse. Treat Mash as a fast prefilter for relatedness, not a substitute for SNP-level comparison of the contextual sequences you select.
- **Entrez rate limits.** `get_sra.py` has no API key or throttling. Without a key NCBI allows about 3 requests/second, so large runs may be slowed or rejected. Lookups are cached per BioSample.
- **Runtime of the ranking step.** `get_mash_dist.sh` scans the full distance table once per reference, so it scales with (references × assemblies).
- **Rerunning into the same outdir is destructive.** `get_mash_dist.sh` deletes and recreates its output CSVs, and `make_ref_db.sh` deletes the old DB, in whichever directory `OUTDIR` points at. Use a fresh `--outdir` to keep earlier results.
- **`xargs` warning.** `get_mash_dist.sh` prints `options --max-args and --replace/-I/-i are mutually exclusive`. It's harmless (the script uses one assembly per job either way) and comes from the original `xargs -n1 -I{}` call.
- `get-taxon-genomes.sh`'s header comment says the default safety margin is 1.3; the code default is 1.2.

## Contributors
Big thank you to the following for their guidance and feedback:
- [Shawn Hawken](DOH-SEH2303)
- [Marcela Torres](DOH-LMT2303)
- [Dahlia Walters](DOH-DEW3303)
