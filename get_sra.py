import csv
import json
import os
import re
import sys
from argparse import ArgumentParser, BooleanOptionalAction
from collections import defaultdict
from datetime import date

from Bio import Entrez

parser = ArgumentParser(
    description=(
        "For each reference genome, walk its ranked list of non-reference "
        "assemblies (closest first, from ref_nonref_ranked.csv) and take the "
        "first one that has an available SRA run. References whose "
        "candidates *all* lack SRA data are skipped, not just the closest one."
    )
)
parser.add_argument("--email", dest="email", help="your.email@something.com as required by Entrez")
parser.add_argument("--db", dest="organism", help="In the format <Genus_species>")
parser.add_argument("--illumina", dest="illumina", action=BooleanOptionalAction)
parser.add_argument("--ont", action=BooleanOptionalAction)
parser.add_argument("--pacbio", action=BooleanOptionalAction)
parser.add_argument(
    "--infile",
    default="ref_nonref_ranked.csv",
    dest="infile",
    help="input csv with 'reference,rank,assembly,mash_distance' in header",
)
parser.add_argument(
    "--max-candidates",
    dest="max_candidates",
    type=int,
    default=None,
    help="optional cap on how many ranked candidates to try per reference before giving up",
)
parser.add_argument(
    "--exclude",
    dest="exclude",
    action="append",
    default=[],
    metavar="IDS_OR_FILE",
    help=(
        "IDs to exclude from the selection. Either a comma/space-separated list or "
        "the path to a text file (one or more IDs per line, '#' starts a comment). "
        "Accepted IDs: GenBank/RefSeq assembly accessions (GCA_/GCF_, with or without "
        "the .version), SRA runs/experiments/samples/studies (SRR/ERR/DRR, SRX/ERX/DRX, "
        "SRS/ERS/DRS, SRP/ERP/DRP) and BioSamples (SAMN/SAMEA/SAMD). Repeat the flag to "
        "give several lists. An excluded candidate is skipped and the next-closest "
        "assembly is considered instead."
    ),
)
parser.add_argument(
    "--validate-exclude-only",
    dest="validate_exclude_only",
    action="store_true",
    help="check the --exclude values for unrecognised IDs, print a summary and exit",
)
parser.add_argument(
    "--outdir",
    dest="outdir",
    default=".",
    help="directory for output csvs (created if needed); used for the default --outfile/--summary-outfile names",
)
parser.add_argument(
    "--outfile",
    dest="outfile",
    default=None,
    help="optional explicit path for the full info csv; defaults to bigbacter.<db>.info.csv",
)
parser.add_argument(
    "--summary-outfile",
    dest="summary_outfile",
    default=None,
    help=(
        "optional explicit path for the sample/taxa/genbank/sra-only csv; "
        "defaults to bigbacter_queue.<db>.closest_to_SRA_<date>.csv"
    ),
)
parser.add_argument(
    "--report-file",
    dest="report_file",
    default="assemblies/assembly_data_report.jsonl",
    help="path to the NCBI assembly_data_report.jsonl used for sequencing-tech lookups",
)

args = parser.parse_args()
Entrez.email = args.email
taxa = args.organism
illumina = args.illumina
ont = args.ont
pacbio = args.pacbio
infile = args.infile
max_candidates = args.max_candidates

tech_list = []
if illumina:
    tech_list.append("illumina")
if pacbio:
    tech_list.append("pacbio")
if ont:
    tech_list.append("ont")

today = str(date.today()).replace("-", "")
assembly_report_file = args.report_file

# Cache SRA lookups per BioSample so we never hit Entrez twice for the same
# assembly, even if it's the nearest candidate for multiple references.
_srr_cache = {}

# GCA_ <-> GCF_ pairs from the assembly report, so excluding one accession of a
# paired assembly also excludes the other.
assembly_pairs = {}


def check_for_any_substring(target_string, substring_list):
    """
    Checks if any substring from a list is present in a target string.
    """
    return any(sub in target_string for sub in substring_list)


def get_seq_technology(report_file=assembly_report_file, tech_list=tech_list):
    """
    Builds a map of genbank accession -> [tech_clean, biosample_accn] by
    parsing the assembly_data_report.jsonl file and keeping only records
    whose sequencing technology matches tech_list.
    """
    assembly_tech_map = {}

    try:
        with open(report_file, "r") as f:
            for line in f:
                record = json.loads(line)
                _acc, _paired = record.get("accession"), record.get("pairedAccession")
                if _acc and _paired:
                    assembly_pairs[_acc] = _paired
                    assembly_pairs[_paired] = _acc
                if "sequencingTech" in record.keys() or "sequencingTech" in record.get("assemblyInfo", {}).keys():

                    if "sequencingTech" in record.keys():
                        tech = record["sequencingTech"].split("; ")
                    else:
                        tech = record["assemblyInfo"]["sequencingTech"].split("; ")
                    tech = [w.lower() for w in tech]
                    tech_clean = []
                    for s in tech:
                        keep = [a for a in tech_list if a in s.lower()]
                        tech_clean.append(" ".join(keep))

                    if set(tech_list) & set(tech_clean):
                        genbank = record["accession"]
                        biosample_accn = record.get("assemblyInfo", {}).get("biosample", {}).get("accession")
                        if genbank not in assembly_tech_map:
                            assembly_tech_map[genbank] = [tech_clean, biosample_accn]

                        paired_accession = record.get("pairedAccession")
                        if paired_accession and paired_accession not in assembly_tech_map:
                            assembly_tech_map[paired_accession] = [tech_clean, biosample_accn]

    except FileNotFoundError:
        raise SystemExit(
            f"Error: Assembly report file not found at '{report_file}'.\n"
            f"This file is required to look up sequencing technology and BioSample "
            f"accessions — without it, every candidate will be rejected. Pass the "
            f"correct location with --report-file, or generate it first "
            f"(e.g. via `datasets summary genome ... --as-json-lines`)."
        )
    except json.JSONDecodeError as e:
        raise SystemExit(f"Error: Could not parse JSON in '{report_file}'. {e}")

    if not assembly_tech_map:
        raise SystemExit(
            f"Error: '{report_file}' was read but produced no usable records for "
            f"tech_list={tech_list!r}. Check that --illumina/--ont/--pacbio matches "
            f"how sequencingTech is recorded in the report, and that the file isn't empty."
        )

    return assembly_tech_map


def get_srr_from_biosample(biosample_accn):
    """
    Return [(run_id, {identifiers}), ...] for the SRR/DRR/ERR runs linked to a
    BioSample accession, filtered by tech_list. The identifier set holds every
    field on the run's runinfo line (run, experiment, study, sample, BioSample,
    ...) in upper case, so the run can be matched against --exclude IDs of any
    of those kinds.
    """
    if biosample_accn in _srr_cache:
        return _srr_cache[biosample_accn]

    srrs = []
    handle = Entrez.esearch(db="sra", term=biosample_accn)
    record = Entrez.read(handle)
    handle.close()
    if not record["IdList"]:
        _srr_cache[biosample_accn] = srrs
        return srrs

    ids = ",".join(record["IdList"])
    runinfo = Entrez.efetch(db="sra", id=ids, rettype="runinfo", retmode="text")
    lines = runinfo.read().splitlines()
    runinfo.close()

    runs = {}
    for line in lines:
        line = line.decode("utf-8")
        if check_for_any_substring(str(line).lower(), tech_list):
            if any(tag in line for tag in ("SRR", "DRR", "ERR")):
                fields = next(csv.reader([line]))
                runs[fields[0]] = frozenset(f.strip().upper() for f in fields if f.strip())

    srrs = sorted(runs.items())
    _srr_cache[biosample_accn] = srrs
    return srrs


def gcf_id_from_filename(filename):
    """e.g. GCF_002120665.1_ASM212066v1_genomic.fna -> GCF_002120665.1"""
    return "_".join(filename.split("_")[:2])


ASSEMBLY_RE = re.compile(r"^GC[AF]_\d+(\.\d+)?$")
SRA_RE = re.compile(r"^[SED]R[RXSP]\d+$")  # run / experiment / sample / study
BIOSAMPLE_RE = re.compile(r"^SAM(N|EA|D)\d+$")


def parse_exclusions(values):
    """
    Turn the --exclude arguments into lookup sets. Each value is either a path to
    a text file or a comma/space-separated list of IDs. Exits with a message if
    any ID is not recognised, so a typo can't silently exclude nothing.
    """
    tokens = []  # (token, where it came from)
    for value in values:
        if os.path.isfile(value):
            with open(value) as fh:
                for lineno, line in enumerate(fh, 1):
                    for tok in re.split(r"[,\s]+", line.split("#", 1)[0].strip()):
                        if tok:
                            tokens.append((tok, f"{value}:{lineno}"))
        else:
            for tok in re.split(r"[,\s]+", value.strip()):
                if tok:
                    tokens.append((tok, "--exclude"))

    excl = {"assembly_exact": set(), "assembly_base": set(), "sra": set(), "biosample": set()}
    bad = []
    for tok, where in tokens:
        t = tok.upper()
        if ASSEMBLY_RE.match(t):
            # with a .version: that exact version only; without: every version
            excl["assembly_exact" if "." in t else "assembly_base"].add(t)
        elif SRA_RE.match(t):
            excl["sra"].add(t)
        elif BIOSAMPLE_RE.match(t):
            excl["biosample"].add(t)
        else:
            bad.append((tok, where))

    if bad:
        lines = [f"  {tok!r} ({where})" for tok, where in bad[:20]]
        more = f"\n  ... and {len(bad) - 20} more" if len(bad) > 20 else ""
        raise SystemExit(
            "Error: these --exclude values are not recognised IDs:\n"
            + "\n".join(lines) + more + "\n"
            "Expected GenBank/RefSeq accessions (GCA_/GCF_), SRA IDs (SRR/ERR/DRR, SRX/ERX/DRX, "
            "SRS/ERS/DRS, SRP/ERP/DRP) or BioSamples (SAMN/SAMEA/SAMD). If one of these is meant "
            "to be a file name, the file was not found."
        )
    excl["all"] = (
        excl["assembly_exact"] | excl["assembly_base"] | excl["sra"] | excl["biosample"]
    )
    return excl


def matched_assembly_exclusion(excl, accession):
    """The exclusion token that matches this assembly (or its GCA_/GCF_ pair), else None."""
    for acc in (accession, assembly_pairs.get(accession)):
        if not acc:
            continue
        a = acc.upper()
        if a in excl["assembly_exact"]:
            return a
        if a.split(".")[0] in excl["assembly_base"]:
            return a.split(".")[0]
    return None


def load_ranked_candidates(path):
    """
    Returns {reference: [(rank, assembly, mash_distance), ...]} sorted by rank,
    skipping rows where no non-reference assembly was found at all.
    """
    candidates = defaultdict(list)
    with open(path) as f:
        reader = csv.reader(f)
        header = next(reader, None)
        for row in reader:
            if not row:
                continue
            reference, rank, assembly, mash_distance = row
            if not assembly:
                continue
            candidates[reference].append((int(rank), assembly, mash_distance))

    for ref in candidates:
        candidates[ref].sort(key=lambda x: x[0])

    return candidates


if __name__ == "__main__":
    exclude = parse_exclusions(args.exclude)  # exits on unrecognised IDs
    if args.validate_exclude_only:
        print(
            f"Exclusions OK: {len(exclude['assembly_exact']) + len(exclude['assembly_base'])} assembly "
            f"accession(s), {len(exclude['sra'])} SRA ID(s), {len(exclude['biosample'])} BioSample(s)"
        )
        sys.exit(0)
    if exclude["all"]:
        print(f"Excluding {len(exclude['all'])} ID(s) from the selection")

    os.makedirs(args.outdir, exist_ok=True)
    outfile = args.outfile or os.path.join(args.outdir, f"bigbacter.{taxa}.info.csv")
    summary_outfile = args.summary_outfile or os.path.join(
        args.outdir, f"bigbacter_queue.{taxa}.closest_to_SRA_{today}.csv"
    )

    candidates_by_ref = load_ranked_candidates(infile)
    data_assembly_json = get_seq_technology()

    n_refs_written = 0
    n_refs_skipped = 0
    n_summary_dupes_skipped = 0
    n_excluded_candidates = 0
    n_excluded_runs = 0
    exclusion_hits = set()
    summary_seen = set()

    with open(outfile, "w", newline="") as out_f, open(summary_outfile, "w", newline="") as summary_f:
        writer = csv.writer(out_f)
        writer.writerow(["reference", "sample", "taxa", "genbank", "sra", "mash_distance", "candidate_rank"])

        summary_writer = csv.writer(summary_f)
        summary_writer.writerow(["sample", "taxa", "genbank", "sra"])

        for reference, candidate_list in candidates_by_ref.items():
            found = False
            reasons = {"not_in_tech_map": 0, "no_biosample": 0, "no_srr": 0, "excluded": 0}
            tried = 0

            for rank, assembly, mash_distance in candidate_list:
                if max_candidates is not None and rank > max_candidates:
                    break
                tried += 1

                gcf_id = gcf_id_from_filename(assembly)

                hit = matched_assembly_exclusion(exclude, gcf_id)
                if hit:
                    exclusion_hits.add(hit)
                    reasons["excluded"] += 1
                    n_excluded_candidates += 1
                    continue

                if gcf_id not in data_assembly_json:
                    reasons["not_in_tech_map"] += 1
                    continue

                biosample = data_assembly_json[gcf_id][1]
                if not biosample:
                    reasons["no_biosample"] += 1
                    continue

                if biosample.upper() in exclude["biosample"]:
                    exclusion_hits.add(biosample.upper())
                    reasons["excluded"] += 1
                    n_excluded_candidates += 1
                    continue

                runs = get_srr_from_biosample(biosample)
                if not runs:
                    reasons["no_srr"] += 1
                    continue

                srrs = []
                for run, ids in runs:
                    matched = ids & exclude["sra"]
                    if matched:
                        exclusion_hits |= matched
                        n_excluded_runs += 1
                    else:
                        srrs.append(run)
                if not srrs:  # every run of this assembly was excluded
                    reasons["excluded"] += 1
                    n_excluded_candidates += 1
                    continue

                for srr in srrs:
                    print(reference, gcf_id, biosample, srr, mash_distance, f"rank={rank}")
                    writer.writerow([reference, biosample, taxa, gcf_id, srr, mash_distance, rank])

                    summary_row = (biosample, taxa, gcf_id, srr)
                    if summary_row not in summary_seen:
                        summary_seen.add(summary_row)
                        summary_writer.writerow(summary_row)
                    else:
                        n_summary_dupes_skipped += 1

                found = True
                break  # stop walking candidates for this reference once one succeeds

            if found:
                n_refs_written += 1
            else:
                n_refs_skipped += 1
                print(
                    f"{reference}: no usable SRA among {tried} candidate(s) tried "
                    f"(of {len(candidate_list)} total) — "
                    f"not_in_tech_map={reasons['not_in_tech_map']}, "
                    f"no_biosample={reasons['no_biosample']}, "
                    f"no_srr_runs={reasons['no_srr']}, "
                    f"excluded={reasons['excluded']}"
                )

    print(f"Wrote full info to {outfile}")
    print(f"Wrote sample/taxa/genbank/sra summary to {summary_outfile}")
    print(f"References with an SRA-backed match: {n_refs_written}")
    print(f"References skipped (no candidate had SRA): {n_refs_skipped}")
    print(f"Redundant rows omitted from summary csv: {n_summary_dupes_skipped}")

    if exclude["all"]:
        unmatched = sorted(exclude["all"] - exclusion_hits)
        print(
            f"Exclusions: {n_excluded_candidates} candidate assembly(ies) skipped, "
            f"{n_excluded_runs} run(s) dropped; {len(exclusion_hits)} of {len(exclude['all'])} "
            f"excluded ID(s) matched a candidate that was examined"
        )
        if unmatched:
            shown = ", ".join(unmatched[:20]) + (" ..." if len(unmatched) > 20 else "")
            print(
                f"  Not encountered ({len(unmatched)}): {shown}\n"
                f"  (An ID is only matched if its genome was reached while walking a reference's "
                f"ranked candidates; genomes that rank behind a selected one, or aren't in this "
                f"dataset, are never examined.)"
            )