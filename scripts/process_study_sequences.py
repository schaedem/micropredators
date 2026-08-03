"""
16S amplicon sequence processing (this study, isotope pulse-chase greenhouse experiment).

Two sequencing runs were processed separately through DADA2 (organic:
"nifa"; conventional: "nifa_conv"), then merged and re-classified as one
combined dataset. Run inside the qiime2-amplicon-2023.9 conda
environment (see nifa_dada2.sh for environment activation).

"""
import os, re

project = "/path"  

processors = 20
db = "Silva"
treename = "NIFA"
min_support = 0.9

# One entry per sequencing run: [name, working dir, raw data dir, manifest path]
libraries = [
    ["nifa",      f"{project}/raw_data/processed",      f"{project}/raw_data",      f"{project}/raw_data/nifa.manifest.csv"],
    ["nifa_conv", f"{project}/raw_data_conv/processed", f"{project}/raw_data_conv", f"{project}/raw_data_conv/nifa.conv.manifest.csv"],
]

# Per-library DADA2 trimming parameters: [trim_left_f, trunc_len_f, trim_left_r, trunc_len_r]
trim_dict = {
    "nifa": [10, 200, 10, 200],
    "nifa_conv": [10, 200, 10, 200],
}

classifier_db = "~/databases/silva/silva-138-99-515-806-nb-classifier.qza"


def format_taxonomy(tax_dirty, min_support):
    """Reformat QIIME2 taxonomy.tsv into a phyloseq-compatible tab-separated
    table with one column per rank, filling unclassified ranks and flagging
    low-confidence assignments as 'putative'."""
    output_path = re.sub(".tsv", "-fixed.tsv", tax_dirty)
    full_rank_length = 7
    with open(output_path, "w") as output, open(tax_dirty) as f:
        output.write("\t".join(["ASV", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]) + "\n")
        next(f)
        for line in f:
            read_id, tax_string, confidence = line.strip().split("\t")[:3]
            if db == "Silva":
                tax_string = re.sub("[a-z]__", "", tax_string)
                tax_string = re.sub("_", " ", tax_string)
            full_rank = tax_string.split(";")

            last_classified = full_rank[-1]
            count = 1
            while last_classified == " ":
                last_classified = full_rank[len(full_rank) - count]
                count += 1

            if float(confidence) < float(min_support):
                full_rank[full_rank.index(last_classified)] = "putative " + last_classified
                last_classified = "putative " + last_classified

            for n in range(full_rank.index(last_classified) + 1, full_rank_length):
                if n < len(full_rank):
                    full_rank[n] = "unclassified " + last_classified
                else:
                    full_rank.append("unclassified " + last_classified)

            output.write(read_id + "\t" + "\t".join(full_rank) + "\n")


# ---- Step 1: import demultiplexed paired-end reads -------------------------
for name, proc, raw, manifest in libraries:
    os.system(" ".join([
        "qiime tools import",
        "--type 'SampleData[PairedEndSequencesWithQuality]'",
        f"--input-path {manifest}",
        f"--output-path {proc}/{name}.qza",
        "--input-format PairedEndFastqManifestPhred33",
    ]))
    os.system(" ".join([
        "qiime demux summarize",
        f"--i-data {proc}/{name}.qza",
        f"--o-visualization {proc}/{name}.demux.qzv",
    ]))

# ---- Step 2: remove primers (cutadapt) --------------------------------------
# forward primer (515F): ATTAGAWACCCBDGTAGTCC | reverse primer (806R, revcomp): TTACCGCGGCKGCTGGCAC
for name, proc, raw, manifest in libraries:
    os.system(" ".join([
        "qiime cutadapt trim-paired",
        f"--i-demultiplexed-sequences {proc}/{name}.qza",
        "--p-adapter-f ATTAGAWACCCBDGTAGTCC",
        "--p-adapter-r TTACCGCGGCKGCTGGCAC",
        "--p-minimum-length 100",
        "--p-discard-untrimmed",
        f"--o-trimmed-sequences {proc}/trimmed",
    ]))
    os.system(" ".join([
        "qiime demux summarize",
        f"--i-data {proc}/trimmed.qza",
        f"--o-visualization {proc}/trimmed.qzv",
    ]))

# ---- Step 3: denoise with DADA2 ---------------------------------------------
for name, proc, raw, manifest in libraries:
    tl_f, trunc_f, tl_r, trunc_r = trim_dict[name]
    os.system(" ".join([
        "qiime dada2 denoise-paired",
        f"--i-demultiplexed-seqs {proc}/trimmed.qza",
        f"--o-representative-sequences {proc}/{name}.rep-seqs.qza",
        f"--o-table {proc}/{name}.table.qza",
        f"--p-trim-left-f {tl_f}", f"--p-trim-left-r {tl_r}",
        f"--p-trunc-len-f {trunc_f}", f"--p-trunc-len-r {trunc_r}",
        "--p-max-ee-f 3", "--p-max-ee-r 3",
        "--p-n-threads", str(processors),
    ]))

# ---- Step 4: summarize ASVs (view outputs at https://view.qiime2.org) ------
for name, proc, raw, manifest in libraries:
    os.system(" ".join(["qiime feature-table summarize",
                        f"--i-table {proc}/{name}.table.qza",
                        f"--o-visualization {proc}/{name}.table.qzv"]))
    os.system(" ".join(["qiime feature-table tabulate-seqs",
                        f"--i-data {proc}/{name}.rep-seqs.qza",
                        f"--o-visualization {proc}/{name}.rep-seqs.qzv"]))

# ---- Step 5: classify taxonomy (SILVA v138, pre-trained classifier) --------
for name, proc, raw, manifest in libraries:
    os.system(" ".join([
        "qiime feature-classifier classify-sklearn",
        f"--i-classifier {classifier_db}",
        f"--i-reads {proc}/{name}.rep-seqs.qza",
        f"--o-classification {proc}/{name}.taxonomy.qza",
        "--p-n-jobs", str(processors),
    ]))
    os.system(" ".join(["qiime metadata tabulate",
                        f"--m-input-file {proc}/{name}.taxonomy.qza",
                        f"--o-visualization {proc}/{name}.taxonomy-summary.qzv"]))
    os.system(" ".join(["qiime taxa barplot",
                        f"--i-table {proc}/{name}.table.qza",
                        f"--i-taxonomy {proc}/{name}.taxonomy.qza",
                        f"--o-visualization {proc}/{name}.taxa-bar-plots.qzv"]))

# ---- Step 6: merge the two sequencing runs ----------------------------------
os.makedirs(f"{project}/joined_runs", exist_ok=True)
out_table      = f"{project}/joined_runs/merged_table.qza"
out_rep_seqs   = f"{project}/joined_runs/merged_rep-seqs.qza"
out_taxonomy   = f"{project}/joined_runs/merged_taxonomy.qza"

run1_table, run2_table       = f"{project}/raw_data/processed/nifa.table.qza", f"{project}/raw_data_conv/processed/nifa_conv.table.qza"
run1_rep_seqs, run2_rep_seqs = f"{project}/raw_data/processed/nifa.rep-seqs.qza", f"{project}/raw_data_conv/processed/nifa_conv.rep-seqs.qza"

os.system(f"qiime feature-table merge --i-tables {run1_table} --i-tables {run2_table} --o-merged-table {out_table}")
os.system(f"qiime feature-table merge-seqs --i-data {run1_rep_seqs} --i-data {run2_rep_seqs} --o-merged-data {out_rep_seqs}")

# Re-classify taxonomy on the merged rep-seqs (per-run classifications
# aren't simply concatenable since ASV sets can differ between runs)
os.system(" ".join([
    "qiime feature-classifier classify-sklearn",
    f"--i-classifier {classifier_db}",
    f"--i-reads {out_rep_seqs}",
    f"--o-classification {out_taxonomy}",
    "--p-n-jobs 2",
]))

# ---- Step 7: phylogenetic tree (built from the merged rep-seqs) ------------
tree_dir = f"{project}/tree"
os.makedirs(tree_dir, exist_ok=True)
os.system(f"qiime alignment mafft --i-sequences {out_rep_seqs} --o-alignment {tree_dir}/{treename}.rep-seqs-merged-aligned.qza --p-n-threads {processors}")
os.system(f"qiime alignment mask --i-alignment {tree_dir}/{treename}.rep-seqs-merged-aligned.qza --o-masked-alignment {tree_dir}/{treename}.rep-seqs-merged-aligned-masked.qza")
os.system(f"qiime phylogeny fasttree --i-alignment {tree_dir}/{treename}.rep-seqs-merged-aligned-masked.qza --o-tree {tree_dir}/{treename}.tree-unrooted.qza --p-n-threads {processors}")
os.system(f"qiime phylogeny midpoint-root --i-tree {tree_dir}/{treename}.tree-unrooted.qza --o-rooted-tree {tree_dir}/{treename}.tree-rooted.qza")
os.system(f"qiime tools export --input-path {tree_dir}/{treename}.tree-rooted.qza --output-path {tree_dir}/")
os.rename(f"{tree_dir}/tree.nwk", f"{tree_dir}/{treename}.tree-final.nwk")

# ---- Step 8: export merged ASV table, taxonomy, and rep-seqs ---------------
final_dir = f"{project}/joined_runs/final"
os.makedirs(final_dir, exist_ok=True)

os.system(f"qiime tools export --input-path {out_table} --output-path {final_dir}")
os.system(f"qiime tools export --input-path {out_taxonomy} --output-path {final_dir}")
os.system(f"qiime tools export --input-path {out_rep_seqs} --output-path {final_dir}")

format_taxonomy(f"{final_dir}/taxonomy.tsv", min_support)

# Convert feature table from .biom to .tsv, with ASV-labeled header
count_biom = f"{final_dir}/feature-table.biom"
count_tsv  = f"{final_dir}/feature-table.tsv"
tmp_tsv, tmp2_tsv = f"{final_dir}/tmp.tsv", f"{final_dir}/tmp2.tsv"

os.system(f"biom convert -i {count_biom} -o {tmp_tsv} --to-tsv")
os.system(f"tail -n +2 {tmp_tsv} > {tmp2_tsv}")          # drop the comment header line
os.system(f"sed 's/#OTU ID/ASV/g' {tmp2_tsv} > {count_tsv}")
os.remove(tmp_tsv)
os.remove(tmp2_tsv)

print(f"Done. Final feature table, taxonomy, and rep-seqs written to {final_dir}")
