# Figure 5 (meta-analysis) -- Step 1: per-study differential abundance
#
# Run once per study in the meta-analysis (31 studies total). Reads one
# study's ASV feature table + taxonomy + sample metadata, aggregates ASVs
# to the order level (trophic-group classification below), and runs
# DESeq2 (rhizosphere vs. bulk soil) with ashr LFC shrinkage. Output feeds
# Step 2 (fig5_meta_analysis.R), which pools these per-study results
# across studies via a hierarchical meta-regression.
#
# Shown here for one representative study (Ayiti_2022); the same script
# was re-run per study with the working directory, feature table, and
# metadata file names adjusted accordingly.

library(tidyverse)
library(DESeq2)

# ---- Shared helper functions (identical across all study runs) ------------

# Classify ASVs into an order-level analysis group and predator/non-predator
# type, based on taxonomy previously associated with predatory lifestyles.
clean_taxonomy <- function(tax_df) {
  bad_terms <- "uncultured|unidentified|putative|Incertae_Sedis|NA|unknown|\\bAmbiguous\\b"

  tax_df %>%
    dplyr::rename(ASV = `Feature ID`) %>%
    mutate(
      # fall back to the next-highest rank when Order is missing/ambiguous
      CleanOrder = case_when(
        !str_detect(Order, bad_terms) & !is.na(Order)   ~ Order,
        !str_detect(Class, bad_terms) & !is.na(Class)   ~ paste0(Class, "_uncl"),
        !str_detect(Phylum, bad_terms) & !is.na(Phylum) ~ paste0(Phylum, "_uncl"),
        TRUE ~ paste0(Domain, "_uncl")
      ),
      AnalysisGroup = ifelse(Genus == "Lysobacter", "g__Lysobacter", paste0("o__", CleanOrder)),
      PlotLabel = ifelse(Genus == "Lysobacter", paste0(Phylum, " | Lysobacter"), paste0(Phylum, " | ", CleanOrder)),
      pred_type = case_when(
        Phylum == "Myxococcota" ~ "Facultative",
        Order %in% c("Streptomycetales", "Cytophagales") | Genus == "Lysobacter" ~ "Facultative",
        Order == "Vampirovibrionales" | Phylum == "Bdellovibrionota" ~ "Obligate",
        TRUE ~ "None"
      )
    )
}

# Aggregate ASV counts to AnalysisGroup, run DESeq2 (rhizosphere vs. bulk),
# and shrink LFC estimates with ashr.
run_deseq_meta <- function(count_subset, metadata_subset, tax_map, study_name, dataset_name) {

  otu_long <- count_subset %>%
    rownames_to_column("ASV") %>%
    inner_join(tax_map, by = "ASV")
  if (nrow(otu_long) == 0) stop(paste("No ASVs from", study_name, "matched the taxonomy map!"))

  group_counts <- otu_long %>%
    group_by(AnalysisGroup) %>%
    summarise(across(where(is.numeric), sum)) %>%
    filter(!is.na(AnalysisGroup) & AnalysisGroup != "") %>%
    column_to_rownames("AnalysisGroup")

  meta_df <- as.data.frame(metadata_subset) %>%
    mutate(SampleID = trimws(as.character(SampleID))) %>%
    filter(!is.na(SampleID) & SampleID != "")

  common_samples <- intersect(colnames(group_counts), meta_df$SampleID)
  if (length(common_samples) < 2) stop(paste("Insufficient overlapping samples (<2) for", study_name, dataset_name))

  group_counts <- group_counts[, common_samples, drop = FALSE]
  meta_df <- meta_df[meta_df$SampleID %in% common_samples, ]
  rownames(meta_df) <- meta_df$SampleID
  meta_df <- meta_df[colnames(group_counts), , drop = FALSE]

  dds <- DESeqDataSetFromMatrix(
    countData = as.matrix(round(group_counts)),
    colData = meta_df,
    design = ~ isolation_source
  )
  if (!"bulk" %in% dds$isolation_source) stop(paste("Reference level 'bulk' not found in", study_name))
  dds$isolation_source <- factor(dds$isolation_source, levels = c("bulk", "rhizosphere"))

  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq(dds, fitType = "local")
  res_shrunk <- lfcShrink(dds, contrast = c("isolation_source", "rhizosphere", "bulk"), type = "ashr")

  as.data.frame(res_shrunk) %>%
    rownames_to_column("AnalysisGroup") %>%
    left_join(tax_map %>% select(AnalysisGroup, PlotLabel, Phylum, pred_type) %>% distinct(), by = "AnalysisGroup") %>%
    mutate(Study = study_name, Dataset = dataset_name, vi = lfcSE^2)
}

# ---- Per-study run ----------------------------------------------------

# Feature (ASV count) table
features <- read_tsv("feature-table.tsv") %>%
  dplyr::rename(OTU = `OTU ID`) %>%
  column_to_rownames("OTU")

read_depth <- data.frame(SampleID = colnames(features), TotalReads = colSums(features))

# Sample metadata: classify each sample as bulk or rhizosphere
metadata <- read_csv("Ayiti_2022.csv") %>%
  dplyr::select(Run, BioProject, BioSample, Experiment, source, lat_lon, `Sample Name`) %>%
  dplyr::rename(SampleID = Run) %>%
  mutate(isolation_source = ifelse(grepl("bulk", source), "bulk", "rhizosphere")) %>%
  merge(read_depth)

write_csv(metadata, "Ayiti_2022_meta.csv")

counts <- features %>% select(metadata$SampleID)

# Taxonomy, with rank fallback for unclassified ASVs
tax_df <- read_tsv("taxonomy_clean.tsv") %>%
  mutate(Phylum = ifelse(is.na(Phylum), Domain, Phylum),
         Class = ifelse(is.na(Class), Phylum, Class),
         Order = ifelse(is.na(Order), Class, Order),
         Family = ifelse(is.na(Family), Order, Family),
         Genus = ifelse(is.na(Genus), Family, Genus),
         Species = ifelse(is.na(Species), Genus, Species))
tax_processed <- clean_taxonomy(tax_df)

# Run DESeq2 for this study and save the per-study result
res_data <- run_deseq_meta(counts, metadata, tax_processed, "Ayiti_2022", "data") %>%
  filter(baseMean > 0)

write_csv(res_data, "master_output_Ayiti_2022_ashr.csv")
