"""
16S amplicon sequence processing for the meta-analysis (58 datasets
across 31 studies). Applies an identical processing framework to this
study's own data (see process_study_sequences.py), adapted to handle
per-study differences in primer set and single- vs. paired-end reads.

Each study directory (under base_dir) is processed independently:
create manifest -> import -> trim primers -> DADA2 denoise -> prevalence
filter (>=2 samples) -> taxonomy (vsearch consensus, full-length SILVA)
-> export -> clean taxonomy. Demonstrated below for one representative
study; re-run with study_dir changed for each of the other studies.
"""
import os
import re
import csv
import subprocess
from pathlib import Path
from typing import Tuple

import pandas as pd

# ---- Setup -------------------------------------------------------------
base_dir = Path("")
excluded_dirs = {"scripts"}
ref_seqs = Path("/silva-138-99-seqs.qza")  # full-length SILVA
ref_tax = Path("/silva-138-99-tax.qza")
n_threads = 20

study_dirs = sorted([
    d for d in base_dir.iterdir()
    if d.is_dir() and not d.name.startswith(".") and d.name not in excluded_dirs
])
print(f"Found {len(study_dirs)} study folders for processing:")
for d in study_dirs:
    print(" -", d.name)

# ---- Per-study primer sets ------------------------------------------------
# Read1 carries the reverse complement of the reverse primer, and read2 the
# reverse complement of the forward primer, at the 3' end in some datasets.
primer_settings = {
    "Zheng_2024" : {"fwd": "AAACTYAAAKGAATTGACGG", "rev": "TTACCGCGGCKGCTGGCAC"}, #515F/907R
    "Zhao_2019": {"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "TTACCGCGGCKGCTGGCAC"}, #515F/806R
    "Zhang_2023": {"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "CTGCWGCCNCCCGTAGG"}, #341F/806R
    "Zecchin_2023": {"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "CTGCWGCCNCCCGTAGG"}, #341F/806R
    "Yuan_2022": {"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "TTACCGCGGCKGCTGGCAC"}, #515F/806R
    "Yu_2020": {"fwd": "AAACTYAAAKGAATTGACGG", "rev": "TTACCGCGGCKGCTGGCAC"}, #515F/907R
    "Yang_2022": {"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "TTACCGCGGCKGCTGGCAC"}, #515F/806R
    "Yamamoto_2018": {"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "CTGCWGCCNCCCGTAGG"}, #341F/806R
    "Xu_2021": {"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "CTGCTGCCTCCCGTAGGAGT"}, #338F/806R
    "Xin_2021": {"fwd": "ACTYAAAKGAATTGRCGGGG", "rev" : "TTACCGCGGCKGCTGGCAC"}, #515F/909R
    "Xiao_2023": {"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "TTACCGCGGCKGCTGGCAC"}, #515F/806R
    "Wood_2018": {"fwd": "GGATTAGATACCCBDGTAGTC", "rev": "CTGCWGCCNCCCGTAGG"}, #341F/785R
    "Wolters_2018": {"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "CTGCWGCCNCCCGTAGG"}, #341F/806R
    "Wei_2023": {"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "CTGCTGCCTCCCGTAGGAGT"}, #338F/806R
    "Wang_2024": {"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "TTACCGCGGCKGCTGGCAC"}, #515F/806R
    "Wang_2023b": {"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "CTGCTGCCTCCCGTAGGAGT"}, #338F/806R
    "Wang_2022": {"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "CTGCTGCCTCCCGTAGGAGT"}, #338F/806R
    "Voulgari_2024": {"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "CTGCWGCCNCCCGTAGG"}, #341F/806R
    "Visioli_2018": {"fwd": "AGCAGCCGCGGTAAT", "rev" : "CTGCTGCSYCCCGTAGG"}, #Probio_Uni / Probio_rev
    "Szoboszlay_2019": {"fwd": "GGATTAGATACCCBDGTAGTC", "rev" : "TTACCGCGGCKGCTG"}, #519F/785R
    "Sun_2018": {"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "CTGCTGCCTCCCGTAGGAGT"}, #338F/806R
    "Su_2021" : {"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "CTGCWGCCNCCCGTAGG"}, #341F/806R
    "Song_2020": {"fwd": "AGGATTAGATACCCT", "rev" : "CTGCTGCCTYCCGTA"}, #343F/798R
    "Rilling_2025": {"fwd": "GGATTAGATACCCBDGTAGTC", "rev": "CTGCWGCCNCCCGTAGG"}, #341F/805R
    "Ren_2021": {"fwd": "AAACTYAAAKGAATTGACGG", "rev": "TTACCGCGGCKGCTGGCAC"}, #515F/907R
    "Rathore_2017": {"fwd": "GGATTAGATACCCBDGTAGTC", "rev": "CTGCWGCCNCCCGTAGG"}, #341F/785R
    "Pang_2022": {"fwd": "GGATTAGATACCCBDGTAGTC", "rev": "CTGCWGCCNCCCGTAGG"}, #341F/805R
    "Mapelli_2020": {"fwd": "GGATTAGATACCCBDGTAGTC", "rev": "CTGCWGCCNCCCGTAGG"}, #341F/785R
    "Lu_2017": {"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "TTACCGCGGCKGCTGGCAC"}, #515F/806R
    "Lopes_2023":{"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "TTACCGCGGCKGCTGGCAC"}, #515F/806R
    "Lopes_2016": {"fwd": "GGATTAGATACCCBDGTAGTC", "rev": "CTGCWGCCNCCCGTAGG"}, #341F/785R
    "Liu_2020b": {"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "TTACCGCGGCKGCTGGCAC"}, #515F/806R
    "Liu_2020": {"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "TTACCGCGGCKGCTGGCAC"}, #515F/806R
    "Liu_2019": {"fwd": "AAACTYAAAKGAATTGACGG", "rev": "TTACCGCGGCKGCTGGCAC"}, #515F/907R
    "Liu_2018b": {"fwd": "AAACTYAAAKGAATTGACGG", "rev": "TTACCGCGGCKGCTGGCAC"}, #515F/907R
    "Lee_2021": {"fwd": "TCGTCGGCAGCGTC", "rev" : "AGTCTCGTGGGCTCGG"}, #502F/701R
    "Kong_2020": {"fwd": "CAACGAGCGCAACCCT", "rev" : "CMGGGTATCTAATCCKGTT"},#799F/1115R
    "Jia_2024": {"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "TTACCGCGGCKGCTGGCAC"}, #515F/806R
    "Iniesta-Pallares_2023": {"fwd": "AAACTYAAAKGAATTGACGG", "rev": "TTACCGCGGCKGCTGGCAC"}, #515F/907R
    "Hinsu_2021": {"fwd": "GGATTAGATACCCBDGTAGTC", "rev": "CTGCWGCCNCCCGTAGG"}, #341F/785R
    "Guo_2025": {"fwd": "AAACTYAAAKGAATTGACGG", "rev": "TTACCGCGGCKGCTGGCAC"}, #515F/907R
    "Gagnon_2020": {"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "TTACCGCGGCKGCTGGCAC"}, #515F/806R
    "Fu_2023": {"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "CTGCWGCCNCCCGTAGG"}, #341F/806R
    "Fu_2020": {"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "TTACCGCGGCKGCTGGCAC"}, #515F/806R - incorrect SRA IDS
    "Fu_2020_new": {"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "CTGCWGCCNCCCGTAGG"}, #341F/806R
    "Fernandez_2016": {"fwd": "AGGTGNTGCATGGYYGTCG", "rev" : "GGGTATCTAATCCY"}, #784R/1046R
    "Fan_2017": {"fwd": "AAACTYAAAKGAATTGACGG", "rev": "TTACCGCGGCKGCTGGCAC"}, #515F/907R
    "Efthymiou_2023": {"fwd": "GGATTAGATACCCBDGTAGTC", "rev": "CTGCWGCCNCCCGTAGG"}, #341F/805R
    "Deng_2024": {"fwd": "CAACGAGCGCAACCCT", "rev" : "CMGGGTATCTAATCCKGTT"}, #799F / 1115R
    "Deng_2022_run2": {"fwd": "AAACTYAAAKGAATTGACGG", "rev": "TTACCGCGGCKGCTGGCAC"}, #515F/907R
    "Deng_2022_run1": {"fwd": "AAACTYAAAKGAATTGACGG", "rev": "TTACCGCGGCKGCTGGCAC"}, #515F/907R
    "Das_2025": {"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "TTACCGCGGCKGCTGGCAC"}, #515F/806R
    "Chiniquy_2021": {"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "TTACCGCGGCKGCTGGCAC"}, #515F/806R
    "Cheng_2019": {"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "CTGCWGCCNCCCGTAGG"}, #341F/806R
    "Bledsoe_2020": {"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "TTACCGCGGCKGCTGGCAC"}, #515F/806R
    "Blau_2019": {"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "CTGCWGCCNCCCGTAGG"}, #341F/806R
    "Ayiti_2022": {"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "TTACCGCGGCKGCTGGCAC"}, #341F/785R > 515F/806R
    "Aponte_2022": {"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "TTACCGCGGCKGCTGGCAC"}, #515F/806R
    "Allard_2016": {"fwd": "ATTAGAWACCCBDGTAGTCC", "rev" : "TTACCGCGGCKGCTGGCAC"} #515F/806R
}

# ---- Step 1: manifest creation (auto-detects single- vs. paired-end) -----
def create_manifest(study_dir: Path) -> Path:
    """Build a QIIME2 manifest from fastq.gz files, auto-detecting read
    type from filename patterns (_R1/_R2 or _1/_2)."""
    fastqs = sorted(study_dir.glob("*.fastq.gz"))
    r1_files = [f for f in fastqs if "_R1" in f.name or "_1" in f.name]
    r2_files = [f for f in fastqs if "_R2" in f.name or "_2" in f.name]
    is_paired = len(r1_files) > 0 and len(r2_files) > 0

    if is_paired:
        data = []
        for fwd in r1_files:
            sample_id = fwd.stem.split("_")[0]
            rev = study_dir / fwd.name.replace("_R1", "_R2").replace("_1", "_2")
            if rev.exists():
                data.append([sample_id, str(fwd.resolve()), "forward"])
                data.append([sample_id, str(rev.resolve()), "reverse"])
            else:
                print(f"WARNING: reverse read not found for {fwd.name}")
        manifest_df = pd.DataFrame(data, columns=["sample-id", "absolute-filepath", "direction"])
    else:
        data = [[f.stem.replace(".fastq", "").replace(".gz", ""), str(f.resolve())] for f in fastqs]
        manifest_df = pd.DataFrame(data, columns=["sample-id", "absolute-filepath"])

    manifest_path = study_dir / "study_manifest.txt"
    manifest_df.to_csv(manifest_path, index=False, sep="\t", quoting=csv.QUOTE_NONE)
    print(f"Manifest written to: {manifest_path} (paired={is_paired})")
    return manifest_path


def detect_read_type_from_manifest(manifest_path: Path) -> str:
    """Reads the manifest header to determine 'paired' or 'single'."""
    with open(manifest_path) as f:
        return "paired" if "direction" in f.readline() else "single"


# ---- Step 2: import reads -------------------------------------------------
def import_reads_from_manifest(study_dir: Path, manifest_path: Path) -> Path:
    read_type = detect_read_type_from_manifest(manifest_path)
    output_dir = study_dir / "qiime_output"
    output_dir.mkdir(exist_ok=True)
    output_path = output_dir / "reads.qza"

    import_type = "SampleData[PairedEndSequencesWithQuality]" if read_type == "paired" else "SampleData[SequencesWithQuality]"
    fmt = "PairedEndFastqManifestPhred33" if read_type == "paired" else "SingleEndFastqManifestPhred33V2"

    subprocess.run([
        "qiime", "tools", "import", "--type", import_type,
        "--input-path", str(manifest_path), "--output-path", str(output_path),
        "--input-format", fmt,
    ], check=True)
    print(f"Imported {read_type}-end reads for {study_dir.name}")
    return output_path


# ---- Step 3: primer trimming (cutadapt) ------------------------------------
def trim_primers(study_dir: Path, qza_reads: Path, read_type: str, primers: dict) -> Path:
    output_dir = study_dir / "qiime_output"
    trimmed = output_dir / "trimmed.qza"
    fwd, rev = primers["fwd"], primers.get("rev")

    if read_type == "paired":
        cmd = ["qiime", "cutadapt", "trim-paired",
               "--i-demultiplexed-sequences", str(qza_reads),
               "--p-adapter-f", fwd, "--p-adapter-r", rev,
               "--p-minimum-length", "100", "--o-trimmed-sequences", str(trimmed)]
    else:
        cmd = ["qiime", "cutadapt", "trim-single",
               "--i-demultiplexed-sequences", str(qza_reads),
               "--p-adapter", fwd, "--p-anywhere", rev,
               "--p-minimum-length", "100", "--o-trimmed-sequences", str(trimmed)]

    subprocess.run(cmd, check=True)
    subprocess.run(["qiime", "demux", "summarize", "--i-data", str(trimmed),
                    "--o-visualization", str(output_dir / "trimmed.qzv")], check=True)
    print(f"Primer trimming complete for {study_dir.name} ({read_type}-end)")
    return trimmed


# ---- Step 4: DADA2 denoising -----------------------------------------------
# Trim/truncation parameters were set per study by inspecting each study's
# demux quality plot (qiime_output/trimmed.qzv) -- left blank here; set
# before calling denoise() for each study. max_ee/trunc_q are common
# starting points, not fixed requirements.
def denoise(study_dir: Path, trimmed_reads: Path, read_type: str, n_threads: int,
            trim_left_f=None, trunc_len_f=None, trim_left_r=None, trunc_len_r=None,
            trim_left=None, trunc_len=None, max_ee: int = 3, trunc_q: int = 20) -> Tuple[Path, Path]:
    output_dir = study_dir / "qiime_output"
    table_qza, repseqs_qza, stats_qza = (output_dir / f"{n}.qza" for n in ["table", "rep_seqs", "stats"])

    if read_type == "paired":
        subprocess.run([
            "qiime", "dada2", "denoise-paired",
            "--i-demultiplexed-seqs", str(trimmed_reads),
            "--p-trim-left-f", str(trim_left_f), "--p-trim-left-r", str(trim_left_r),
            "--p-trunc-len-f", str(trunc_len_f), "--p-trunc-len-r", str(trunc_len_r),  # needs >=12bp overlap
            "--p-max-ee-f", str(max_ee), "--p-max-ee-r", str(max_ee), "--p-trunc-q", str(trunc_q),
            "--p-n-threads", str(n_threads),
            "--o-table", str(table_qza), "--o-representative-sequences", str(repseqs_qza),
            "--o-denoising-stats", str(stats_qza),
        ], check=True)
    else:
        subprocess.run([
            "qiime", "dada2", "denoise-single",
            "--i-demultiplexed-seqs", str(trimmed_reads),
            "--p-trim-left", str(trim_left), "--p-trunc-len", str(trunc_len),
            "--p-trunc-q", str(trunc_q), "--p-max-ee", str(max_ee),
            "--p-n-threads", str(n_threads),
            "--o-table", str(table_qza), "--o-representative-sequences", str(repseqs_qza),
            "--o-denoising-stats", str(stats_qza),
        ], check=True)

    subprocess.run(["qiime", "feature-table", "summarize", "--i-table", str(table_qza),
                    "--o-visualization", str(output_dir / "table.qzv")], check=True)
    subprocess.run(["qiime", "feature-table", "tabulate-seqs", "--i-data", str(repseqs_qza),
                    "--o-visualization", str(output_dir / "rep_seqs.qzv")], check=True)
    subprocess.run(["qiime", "metadata", "tabulate", "--m-input-file", str(stats_qza),
                    "--o-visualization", str(output_dir / "stats.qzv")], check=True)

    print(f"Denoising complete for {study_dir.name}. Outputs in: {output_dir}")
    return table_qza, repseqs_qza


# ---- Step 5: prevalence filter (ASVs present in >=2 samples) ---------------
def filter_low_prevalence(study_dir: Path) -> Tuple[Path, Path]:
    output_dir = study_dir / "qiime_output"
    table, rep_seqs = output_dir / "table.qza", output_dir / "rep_seqs.qza"
    table_final, rep_seqs_final = output_dir / "table-final.qza", output_dir / "rep-seqs-final.qza"

    subprocess.run(["qiime", "feature-table", "filter-features", "--i-table", str(table),
                    "--p-min-samples", "2", "--o-filtered-table", str(table_final)], check=True)
    subprocess.run(["qiime", "feature-table", "filter-seqs", "--i-data", str(rep_seqs),
                    "--i-table", str(table_final), "--o-filtered-data", str(rep_seqs_final)], check=True)
    print(f"Prevalence filtering complete for {study_dir.name}")
    return table_final, rep_seqs_final


# ---- Step 6: taxonomy (vsearch consensus classifier, full-length SILVA) ----
def assign_taxonomy(study_dir: Path, ref_seqs: Path, ref_tax: Path) -> Path:
    output_dir = study_dir / "qiime_output"
    taxonomy_qza, hits_qza = output_dir / "taxonomy.qza", output_dir / "vsearch_hits.qza"

    subprocess.run([
        "qiime", "feature-classifier", "classify-consensus-vsearch",
        "--i-query", str(output_dir / "rep-seqs-final.qza"),
        "--i-reference-reads", str(ref_seqs), "--i-reference-taxonomy", str(ref_tax),
        "--p-perc-identity", "0.80", "--p-threads", "12",
        "--p-top-hits-only", "TRUE", "--p-maxaccepts", "10",
        "--o-classification", str(taxonomy_qza), "--o-search-results", str(hits_qza),
    ], check=True)
    print(f"Taxonomy classification complete for {study_dir.name}")
    return taxonomy_qza


# ---- Step 7: export ---------------------------------------------------------
def export_results(study_dir: Path):
    output_dir = study_dir / "qiime_output"
    final_dir = output_dir / "qiime_final"
    final_dir.mkdir(exist_ok=True)

    subprocess.run(["qiime", "tools", "export", "--input-path", str(output_dir / "taxonomy.qza"),
                    "--output-path", str(final_dir)], check=True)
    subprocess.run(["qiime", "tools", "export", "--input-path", str(output_dir / "table-final.qza"),
                    "--output-path", str(final_dir)], check=True)
    subprocess.run(["biom", "convert", "-i", str(final_dir / "feature-table.biom"),
                    "-o", str(final_dir / "feature-table.tsv"), "--to-tsv"], check=True)
    print(f"Export complete for {study_dir.name}")
    return final_dir


# ---- Step 8: clean taxonomy (split ranks, strip prefixes, blank -> NA) -----
def clean_taxonomy(study_dir: Path):
    final_dir = study_dir / "qiime_output" / "qiime_final"
    taxonomy_tsv, out_tsv = final_dir / "taxonomy.tsv", final_dir / "taxonomy_clean.tsv"

    ranks = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    prefix_re = re.compile(r"^(d|p|c|o|f|g|s)__\s*", flags=re.IGNORECASE)

    df = pd.read_csv(taxonomy_tsv, sep="\t", dtype=str)
    df["Consensus"] = pd.to_numeric(df["Consensus"], errors="coerce")

    split = df["Taxon"].fillna("").str.split(";", expand=True).iloc[:, :len(ranks)]
    while split.shape[1] < len(ranks):
        split[split.shape[1]] = None
    split.columns = ranks

    for c in ranks:
        split[c] = (
            split[c].astype(str).str.strip()
            .str.replace(prefix_re, "", regex=True).str.strip()
            .replace({"": pd.NA, "unassigned": pd.NA, "Unassigned": pd.NA,
                      "unclassified": pd.NA, "unknown": pd.NA, "NA": pd.NA, "None": pd.NA})
        )

    out = pd.concat([df[["Feature ID", "Consensus"]], split], axis=1)
    out.to_csv(out_tsv, sep="\t", index=False, na_rep="NA")
    print(f"Taxonomy cleaned for {study_dir.name}: {out.shape[0]} ASVs -> {out_tsv}")


# ---- Full pipeline for one study (repeat with study_dir changed for each) --
study_dir = Path("/home/marie/predators/Sun_2018")  # <- change per study

manifest_path = create_manifest(study_dir)
read_type = detect_read_type_from_manifest(manifest_path)
primers = primer_settings.get(study_dir.name)
if primers is None:
    raise ValueError(f"No primer sequences found for study: {study_dir.name}")

reads_qza = import_reads_from_manifest(study_dir, manifest_path)
trimmed_qza = trim_primers(study_dir, reads_qza, read_type, primers)

# Inspect qiime_output/trimmed.qzv (https://view.qiime2.org) and set these
# based on this study's read length/quality before denoising -- left blank,
# adjusted per study rather than fixed.
trim_left_f = trim_left_r = trunc_len_f = trunc_len_r = None  # paired-end
trim_left = trunc_len = None                                   # single-end

table_qza, repseqs_qza = denoise(
    study_dir, trimmed_qza, read_type, n_threads,
    trim_left_f=trim_left_f, trunc_len_f=trunc_len_f, trim_left_r=trim_left_r, trunc_len_r=trunc_len_r,
    trim_left=trim_left, trunc_len=trunc_len,
)
table_final, repseqs_final = filter_low_prevalence(study_dir)
taxonomy_qza = assign_taxonomy(study_dir, ref_seqs, ref_tax)
export_results(study_dir)
clean_taxonomy(study_dir)
