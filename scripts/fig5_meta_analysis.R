# Figure 5 -- Step 2: across-study meta-regression
#
# Pools the per-study DESeq2 results (fig5_per_study_deseq.R output,
# *_ashr.csv, one file per of the 31 studies) into a single dataset,
# aggregates order-level log2FoldChange estimates into taxonomic trophic
# groups via inverse-variance weighting, then fits a hierarchical
# mixed-effects meta-regression (metafor::rma.mv) to estimate global
# rhizosphere enrichment for each trophic group. Produces Figure 5
# (combined_plots_by_order_ashr.png).

library(tidyverse)
library(metafor)
library(cowplot)
library(ggtext)

# ---- 1. Compile per-study results --------------------------------------

study_dirs <- list.dirs()

combined <- study_dirs %>%
  map_df(~ {
    files <- list.files(.x, pattern = "*_ashr\\.csv$", recursive = TRUE, full.names = TRUE)
    if (length(files) > 0) map_df(files, read_csv, show_col_types = FALSE) else NULL
  })

# Deng_2022 was sequenced across two runs; treat as one study, two datasets
combined <- combined %>%
  mutate(
    Dataset = ifelse(Study %in% c("run1", "run2"), Study, Dataset),
    Study   = ifelse(Study %in% c("run1", "run2"), "Deng_2022", Study),
    predator = ifelse(pred_type %in% c("Obligate", "Facultative"), "Predator", "Non-predator")
  )

# Join primer set / sample size metadata (needed for random effects below)
sample_metadata <- read_csv("meta_analysis/metadata/study_sample_numbers.csv", show_col_types = FALSE)

master_df <- combined %>%
  left_join(sample_metadata, by = c("Study", "Dataset")) %>%
  filter(!is.na(log2FoldChange) & !is.na(vi))

write_csv(master_df, "final_predator_meta_analysis_data_ashr.csv")

# ---- 2. Aggregate order-level LFC into trophic groups ------------------
# Inverse-variance weighted mean within each Study x Dataset x trophic
# group (Supplemental Methods Eq. 1-2).

master_aggregated <- master_df %>%
  mutate(Unit = case_when(
    grepl("Myxococcota", PlotLabel)        ~ "Myxococcota",
    grepl("Streptomycetales", AnalysisGroup) ~ "Streptomycetales",
    grepl("Bdellovibrionota", PlotLabel)   ~ "Bdellovibrionota",
    grepl("Lysobacter", AnalysisGroup)     ~ "Lysobacter",
    grepl("Agromyces", AnalysisGroup)     ~ "Agromyces",
    grepl("Cytophagales", AnalysisGroup)   ~ "Cytophagales",
    grepl("Vampirovibrionales", AnalysisGroup) ~ "Vampirovibrionales",
    TRUE ~ "Non-Predator"
  )) %>%
  group_by(Study, Dataset, primers16S, pred_type, Unit, plant_functional_group) %>%
  summarise(
    yi = sum(log2FoldChange / vi) / sum(1 / vi),
    V  = 1 / sum(1 / vi),
    .groups = "drop"
  )

# Cap precision weights at 500 to prevent high-depth/high-abundance datasets
# from dominating the model (Supplemental Methods). A sensitivity check
# across thresholds 100-10,000 confirmed estimates are stable at >=500.
compare_caps <- function(cap_val, data) {
  dat <- data %>% mutate(V_temp = 1 / pmin(1 / V, cap_val))
  mod <- rma.mv(yi ~ Unit - 1, V = V_temp,
                random = list(~ 1 | primers16S, ~ 1 | Study / Dataset), data = dat)
  as.vector(mod$beta)
}
cap_sensitivity <- sapply(c(100, 300, 500, 1000, 10000), compare_caps, data = master_aggregated)
# max absolute difference from the chosen cap (500) across all thresholds:
max(abs(cap_sensitivity - cap_sensitivity[, 3]))

master_weighted <- master_aggregated %>%
  filter(!is.na(yi) & !is.infinite(1 / V)) %>%
  mutate(V_stable = 1 / pmin(1 / V, 500),
         Unit = relevel(factor(Unit), ref = "Non-Predator"))

# ---- 3. Hierarchical meta-regression (metafor::rma.mv) -----------------
# Random effects: primer set (amplification bias) and dataset nested
# within study (within- vs. between-study variance).

# Non-predator as reference level -- for the "relative to non-predators" panel
res_meta_final <- rma.mv(
  yi = yi, V = V_stable, mods = ~ Unit,
  random = list(~ 1 | primers16S, ~ 1 | Study / Dataset),
  data = master_weighted, method = "ML"
)

# Does adding plant functional group improve the model? (reported in Results)
res_meta_plant <- rma.mv(
  yi = yi, V = V_stable, mods = ~ Unit + plant_functional_group,
  random = list(~ 1 | primers16S, ~ 1 | Study / Dataset),
  data = master_weighted, method = "ML"
)
anova(res_meta_final, res_meta_plant)

# No-intercept specification -- for the "relative to bulk soil" panel
res_global <- rma.mv(
  yi = yi, V = V_stable, mods = ~ Unit - 1,
  random = list(~ 1 | primers16S, ~ 1 | Study / Dataset),
  data = master_weighted, method = "REML"
)

# ---- 4. Figure 5: forest plots ------------------------------------------

strategy_map <- data.frame(
  Unit = c("Cytophagales", "Streptomycetales", "Lysobacter","Agromyces",
           "Bdellovibrionota", "Vampirovibrionales", "Myxococcota", "Non-Predator"),
  pred_type = c("Facultative", "Facultative", "Facultative",
                "Obligate", "Obligate", "Facultative", "Non-Predator")
)

study_counts <- master_weighted %>%
  group_by(Unit) %>%
  summarise(n_studies = n_distinct(paste0(Study, Dataset)))

sig_label <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")

# Panel B: enrichment relative to non-predators
plot_data_final <- data.frame(
  Unit = gsub("Unit", "", rownames(res_meta_final$beta)),
  LFC = as.vector(res_meta_final$beta),
  ci.lb = res_meta_final$ci.lb, ci.ub = res_meta_final$ci.ub, pval = res_meta_final$pval
) %>%
  filter(Unit != "intrcpt") %>%
  left_join(study_counts, by = "Unit") %>%
  inner_join(strategy_map, by = "Unit") %>%
  mutate(pred_type = factor(pred_type, levels = c("Obligate", "Facultative")),
         sig_stars = sig_label(pval),
         Unit = paste0("<i>", Unit),
         Unit_label = reorder(Unit, LFC),
         star_pos = ifelse(LFC >= 0, ci.ub + 0.01, ci.lb - 0.01),
         star_hjust = ifelse(LFC >= 0, 0, 1))

relative_plot <- ggplot(plot_data_final, aes(x = LFC, y = Unit_label, fill = pred_type)) +
  geom_col() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 2) +
  geom_errorbarh(aes(xmin = ci.lb, xmax = ci.ub), width = 0.1) +
  geom_point(shape = 18, size = 5, color = "black", show.legend = FALSE) +
  geom_text(aes(x = star_pos, label = sig_stars, hjust = star_hjust), vjust = 0.5, size = 8, color = "black") +
  scale_fill_manual(values = c("Obligate" = "#2c7bb6", "Facultative" = "#d7191c")) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.15)), limits = c(-0.3, 0.55), breaks = c(-0.25, 0, 0.25, 0.5)) +
  labs(x = "Rhizosphere effect relative to non-predators", y = NULL) +
  theme_linedraw() +
  theme(axis.title = element_text(size = 15), axis.text.y = element_markdown(size = 15),
        axis.text.x = element_text(size = 15), panel.grid.minor = element_blank(),
        legend.title = element_blank(), legend.text = element_text(size = 15),
        legend.position = c(0.99, 0.01), legend.justification = c(1, 0),
        legend.background = element_rect(color = "black", fill = "white"))

# Panel A: enrichment relative to bulk soil (no-intercept model)
plot_global <- data.frame(
  Unit = gsub("Unit", "", rownames(res_global$beta)),
  LFC = as.vector(res_global$beta),
  ci.lb = res_global$ci.lb, ci.ub = res_global$ci.ub, pval = res_global$pval
) %>%
  left_join(study_counts, by = "Unit") %>%
  left_join(strategy_map %>% rename(pred_type_all = pred_type), by = "Unit") %>%
  mutate(sig_stars = sig_label(pval),
         Unit = paste0("<i>", Unit),
         Unit_label = reorder(Unit, LFC),
         pred_type_all = factor(pred_type_all, levels = c("Obligate", "Facultative", "Non-Predator")),
         star_pos = ifelse(LFC >= 0, ci.ub + 0.01, ci.lb - 0.01),
         star_hjust = ifelse(LFC >= 0, 0, 1))

global_plot <- ggplot(plot_global, aes(x = LFC, y = Unit_label, fill = pred_type_all)) +
  geom_col() +
  geom_vline(xintercept = 0, color = "black") +
  geom_errorbarh(aes(xmin = ci.lb, xmax = ci.ub), width = 0.1, color = "black", show.legend = FALSE) +
  geom_point(shape = 18, size = 5, color = "black", show.legend = FALSE) +
  geom_text(aes(x = star_pos, label = sig_stars, hjust = star_hjust), vjust = 0.5, size = 8, color = "black") +
  scale_fill_manual(values = c("#1F78B4", "#E31A1C", "darkgrey"),
                     labels = c("Obligate", "Facultative", "Non-Predator")) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.15)), limits = c(-0.3, 0.55), breaks = c(-0.25, 0, 0.25, 0.5)) +
  labs(x = "Rhizosphere effect", y = NULL) +
  theme_linedraw() +
  theme(axis.title = element_text(size = 15), axis.text.y = element_markdown(size = 15),
        axis.text.x = element_text(size = 15), panel.grid.minor = element_blank(),
        legend.title = element_blank(), legend.text = element_text(size = 15),
        legend.position = c(0.99, 0.01), legend.justification = c(1, 0),
        legend.background = element_rect(color = "black", fill = "white"))

# Figure 5: (A) no-intercept model, (B) relative to non-predators
combined_figure <- plot_grid(global_plot, relative_plot, labels = "AUTO", align = "hv", axis = "tb")
combined_figure
ggsave("Fig5_meta_analysis.png", combined_figure, height = 5, width = 14)
