# Figure 4 -- Step 1: per-timepoint differential abundance (this study)
#
# DESeq2 (rhizosphere vs. bulk) run separately for each of the 8
# soil-background x plant x timepoint combinations, ASVs aggregated to
# the Order level (Lysobacter kept at genus level), LFC shrunk via ashr.
# Output (4 CSVs, one per plant x treatment, each with all 4 timepoints
# stacked) feeds Step 2 (fig4_rhizosphere_effect.R).

library(tidyverse)
library(phyloseq)
library(DESeq2)

physeq <- readRDS("nifa_joint.phyloseq.RDS")

# Order-level grouping (genus-level for Lysobacter), with rank fallback
# for unclassified ASVs, and predator-strategy classification.
bad_terms <- "uncultured|unidentified|putative|Incertae_Sedis|NA|unknown|\\bAmbiguous\\b"

tax_df <- as.data.frame(tax_table(physeq)) %>%
  rownames_to_column("ASV") %>%
  mutate(
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

# Aggregate to AnalysisGroup, run DESeq2 (rhizosphere vs. bulk) for one
# plant x treatment x timepoint subset, shrink LFC with ashr.
run_deseq_from_phyloseq <- function(obj, plant_val, treatment_val, timepoint_val,
                                    group_var = "sample_type",
                                    base_level = "bulk",
                                    contrast_level = "rhizo") {

  ps_sub <- eval(parse(text = paste0(
    "subset_samples(obj, plant == '", plant_val,
    "' & treatment == '", treatment_val,
    "' & timepoint == '", timepoint_val, "')"
  )))
  sample_data(ps_sub)$sample_type <- factor(sample_data(ps_sub)$sample_type)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)

  count_subset <- as.data.frame(as(otu_table(ps_sub), "matrix"))
  if (!taxa_are_rows(ps_sub)) count_subset <- as.data.frame(t(count_subset))

  otu_long <- count_subset %>%
    rownames_to_column("ASV") %>%
    inner_join(tax_df, by = "ASV")
  if (nrow(otu_long) == 0) stop("No ASVs matched the taxonomy map!")

  group_counts <- otu_long %>%
    group_by(AnalysisGroup) %>%
    summarise(across(where(is.numeric), sum)) %>%
    filter(!is.na(AnalysisGroup) & AnalysisGroup != "") %>%
    column_to_rownames("AnalysisGroup")

  meta_df <- as(sample_data(ps_sub), "data.frame")
  meta_df$sample_type <- factor(meta_df$sample_type, levels = c(base_level, contrast_level))

  common_samples <- intersect(colnames(group_counts), rownames(meta_df))
  if (length(common_samples) < 2) stop("Insufficient overlapping samples (<2)")
  group_counts <- group_counts[, common_samples, drop = FALSE]
  meta_df <- meta_df[common_samples, , drop = FALSE]

  dds <- DESeqDataSetFromMatrix(countData = as.matrix(round(group_counts)),
                                colData = meta_df, design = ~ sample_type)
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq(dds, fitType = "local")
  res_shrunk <- lfcShrink(dds, contrast = c(group_var, contrast_level, base_level), type = "ashr")

  as.data.frame(res_shrunk) %>%
    rownames_to_column("AnalysisGroup") %>%
    left_join(tax_df %>% select(AnalysisGroup, PlotLabel, Phylum, pred_type) %>% distinct(), by = "AnalysisGroup") %>%
    mutate(vi = lfcSE^2)
}

# Run all 4 timepoints for one plant x treatment combination, stack, save.
run_plant_treatment <- function(plant_val, treatment_val, out_file) {
  results <- purrr::map_dfr(c("A", "B", "C", "D"), function(tp) {
    run_deseq_from_phyloseq(physeq, plant_val, treatment_val, tp) %>%
      mutate(plant = plant_val, treatment = treatment_val, timepoint = tp)
  })
  write.csv(results, out_file)
  results
}

org_vetch_results  <- run_plant_treatment("vetch", "org",  "org_vetch_results_order_level.csv")
org_wheat_results  <- run_plant_treatment("wheat", "org",  "org_wheat_results_order_level.csv")
conv_vetch_results <- run_plant_treatment("vetch", "conv", "conv_vetch_results_order_level.csv")
conv_wheat_results <- run_plant_treatment("wheat", "conv", "conv_wheat_results_order_level.csv")

all_results <- rbind(conv_wheat_results, conv_vetch_results, org_vetch_results, org_wheat_results)
write_csv(all_results, "rhizo_enriched_order_level.csv")
