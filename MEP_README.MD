# Quickstart for Molecular Epi
**0. Set up once**

```bash
conda env create -f environment.yml
conda activate episcope
chmod +x *.sh
```

Add your `config.csv` file to the head directory of this repository

**1. Download public assemblies for your taxon** (`get-taxon-genomes.sh`)

This will download all Morganella_morganii released on NCBI on on or after 2016-11-02 and will exclude all a-typical genomes that have USA as a substring as geo-location. Atypical genomes should exlude genomes that are misattributed to your taxa, but you may with she remove this flag if you don't see closely related genomes to your reference (production) seqs. 582 is the taxon ID for Morganella_morganii as can be see on [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/)
```bash
bash get-taxon-genomes.sh 582  --include genome  --mag exclude  --assembly-version latest  --released-after 2016-11-02  --exclude-atypical --geo-location "USA"
```

This will check to see how many genomes pass this filter and see if I have enough local memory to complete the download. The following Shigella query will not begin downloading seqs because it will be too large for my local memory to handle - it will recommend you aler your query OR repeate the analysis in a location with larger memory.
```bash
get-taxon-genomes.sh "Shigella"  --include genome  --mag exclude  --assembly-version latest  --released-after 2016-11-02  --exclude-atypical --geo-location "USA"
```

**2. Collect your cluster reference genomes** (`populate_bb_references.sh`)
If your references live in BigBacter v1 or v2 databases (local or on S3), `populate_bb_references.sh` copies them into `bb_cluster_references/` for you.

```bash
populate_bb_references.sh --config config.csv --bbv2 --species Morganella_morganii
```

**3. Find the closest assembly with SRA data for each reference** (`episcope.sh`)

```bash
episcope.sh --outdir results --tool sourmash --threads 8 --email you@example.org --db Morganella_morganii --illumina --matrix
```

Also useful for if your production seqs have been uploaded you can provide a tsv with one entry per line of any Genbank or SRA ID to exclude from the analysis of finding the closest genomes to your reference dataset with the `--exclude exclude.txt` arg

**4. Examine Output files**
- `results/bigbacter_queue.<bigbacter_db>.closest_to_SRA_YYYYMMDD.csv`: BigBacter manifest
- `results/bigbacter.Morganella_morganii.info.csv`: Contains info on which which reference genomes had the smallest smash distance to which dataset on NCBI that contained both genomes and SRA seqs. Some reference genomes may align to the same genome on NCBI, so there may be fewer lines on the BigBacter manifest than there are reference genomes in order to avoid redundancy
- `distance_matrix.csv`: Long format distance matrix to review mash distance