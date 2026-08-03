setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Marie/Postdoc_Cornell/SOM_Greenhouse/final_data")

library(tidyverse)
library(cowplot)
library(ggpubr)   # stat_cor()
library(emmeans)
library(lme4)
library(car)
library(ggtext)   # geom_richtext()/element_markdown()

# Figure 3, Supplemental Table 4, Supplemental Figure 3.
#
# Predator:prey ratio = ratio of labeled predator abundance to labeled
# non-predator abundance in rarefied rhizosphere communities (abundance
# adjusted for rrn copy number and 16S qPCR, see Supplemental Methods).

# ---- 1. Data preparation ------------------------------------------------
incorps_soil <- readRDS("incorp_soil_adj_abun.RDS")

tax <- read_tsv("taxonomy-fixed.tsv")

conv_meta <- read_csv("conv_sample_metadata.csv")
org_meta <- read_csv("org_sample_metadata_new.csv")

irms <- read_csv("irms_results_all_final.csv") %>%
  mutate(timepoint = case_when(timepoint == 1 ~ "A", timepoint == 3 ~ "B",
                               timepoint == 7 ~ "C", timepoint == 14 ~ "D")) %>%
  rename(sample = `Sample ID`) %>%
  mutate(C_N = perc_C / perc_N)

meta_full <- rbind(conv_meta, org_meta) %>%
  dplyr::rename(sample = sample_names) %>%
  left_join(irms) %>%
  filter(sample %in% incorps_soil$sample) %>%
  column_to_rownames("sample")

# NOTE: not explicitly filtered to soil == "rhizo" here; assumes
# incorps_soil/meta_full are already rhizosphere-only upstream.
incorp_preds <- incorps_soil %>%
  merge(meta_full) %>%
  filter(ASV %in% tax$ASV) %>%
  merge(tax) %>%
  mutate(plant = ifelse(plant == "wheat", "Wheat", "Vetch")) %>%
  mutate(predator = ifelse(grepl("Bdello", Phylum) |
                             grepl("Cytophag", Order) |
                             grepl("Lyso", Genus) |
                             grepl("Myxo", Phylum) |
                             grepl("Streptomycetales", Order) |
                             grepl("Vampiro", Order), "Predator", "Non-predator")) %>%
  mutate(timepoint = str_trim(timepoint)) %>%
  mutate(timepoint = ifelse(timepoint == "A", 1,
                            ifelse(timepoint == "B", 3,
                                   ifelse(timepoint == "C", 7, 14)))) %>%
  mutate(timepoint = as.double(timepoint),
         plant = factor(plant),
         treatment = factor(treatment)) %>%
  filter(log_adj_abund > 0) %>%
  distinct()

## ---- Canonical, complete rhizosphere sample roster ------------------------
# incorp_preds only has rows for ASV-sample combos actually detected as
# labeled, so a sample with zero labeled ASVs in a functional group would
# otherwise disappear rather than showing a true zero. Roster built from
# sample metadata directly; labeled abundances are left-joined on and
# zero-filled below.
rhizo_sample_roster <- rownames_to_column(meta_full, "sample") %>%
  mutate(plant = ifelse(plant == "wheat", "Wheat", "Vetch")) %>%
  dplyr::select(sample, treatment, plant, timepoint, rep, label) %>%
  distinct() %>%
  rename(timepoint_factor = timepoint) %>%
  mutate(timepoint = case_when(timepoint_factor == "A" ~ 1, timepoint_factor == "B" ~ 3,
                               timepoint_factor == "C" ~ 7, timepoint_factor == "D" ~ 14))

# --- Helper: outlier-robust "top of the data" for annotation placement --
# Tukey whisker convention (Q3 + 1.5*IQR, capped at the highest point
# still within that fence).
robust_upper <- function(x) {
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  fence <- q3 + 1.5 * IQR(x, na.rm = TRUE)
  max(x[x <= fence], na.rm = TRUE)
}

# --- Helper: explicit p-value threshold text (p<0.0001 ... p<0.1, ns)
p_label <- function(p) {
  dplyr::case_when(
    p < 0.0001 ~ "p<0.0001", p < 0.001 ~ "p<0.001", p < 0.01 ~ "p<0.01",
    p < 0.05 ~ "p<0.05", p < 0.1 ~ "p<0.1", TRUE ~ "ns"
  )
}
# --- Helper: markdown-bold only if p < 0.05 (rendered via geom_richtext)
bold_if_sig <- function(label, p) ifelse(!is.na(p) & p < 0.05, paste0("**", label, "**"), label)

# --- Helper: a significance label for a single 2-level factor's marginal
# effect from a fitted lmer model (emmeans), rather than a raw two-sample
# test on unmodeled data.
get_marginal_signif <- function(mod, term) {
  emm <- emmeans::emmeans(mod, as.formula(paste("~", term)))
  p <- as.data.frame(pairs(emm))$p.value[1]
  bold_if_sig(p_label(p), p)
}

# --- Helper: strip any previously-added annotation layer before adding a
# new one (idempotent re-runs).
strip_annotations <- function(p) {
  is_annotation_layer <- vapply(p$layers, function(l) {
    inherits(l$geom, "GeomText") || inherits(l$geom, "GeomRichText")
  }, logical(1))
  p$layers <- p$layers[!is_annotation_layer]
  p
}

#### Predator:Prey ratio (abundance-based) ###################################
# Summed per functional group per sample using the untransformed adj_abund
# (not log_adj_abund -- summing already-logged values isn't a meaningful
# total), log2 applied once to the final ratio.
incorp_preds3_abund <- rhizo_sample_roster %>%
  left_join(incorp_preds, by = c("sample", "plant", "treatment", "timepoint", "label", "rep")) %>%
  group_by(sample, plant, treatment, timepoint, label, rep, predator) %>%
  summarise(sum_adj_abund = sum(adj_abund, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(id_cols = c(sample, plant, treatment, timepoint, rep, label),
              names_from = predator, values_from = sum_adj_abund) %>%
  mutate(Predator = tidyr::replace_na(Predator, 0),
         `Non-predator` = tidyr::replace_na(`Non-predator`, 0),
         predator_prey_ratio_abund = Predator / `Non-predator`)

# ---- Figure 3A: soil background (treatment) --------------------------------
soil_background_abund <- incorp_preds3_abund %>%
  mutate(treatment = ifelse(treatment == "conv", "Conventional", "Organic")) %>%
  ggplot(aes(y = predator_prey_ratio_abund, x = treatment, color = treatment, fill = treatment)) +
  theme_linedraw() +
  geom_boxplot(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.6,
               outlier.shape = NA, color = "black") +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), size = 3) +
  labs(y = "Labeled Predator : Labeled Prey abundance") +
  xlab("") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15),
        legend.position = "none") +
  scale_color_manual(values = c("#FED439FF", "#709AE1FF")) +
  scale_fill_manual(values = c("#FED439FF", "#709AE1FF"))

# ---- Figure 3B: plant host --------------------------------------------------
plant_background_abund <- ggplot(incorp_preds3_abund, aes(y = predator_prey_ratio_abund, x = plant, color = plant, fill = plant)) +
  theme_linedraw() +
  geom_boxplot(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.6,
               outlier.shape = NA, color = "black") +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), size = 3) +
  labs(y = "Labeled Predator : Labeled Prey abundance") +
  xlab("") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15),
        legend.position = "none") +
  scale_color_manual(values = c("#7570B3FF", "#1B9E77FF")) +
  scale_fill_manual(values = c("#7570B3FF", "#1B9E77FF"))

# Time-course base (Figure 3C), brackets/trend line added below
ratios_time_abund_base <- ggplot(incorp_preds3_abund,
                                 aes(x = factor(timepoint), y = predator_prey_ratio_abund, fill = treatment)) +
  theme_linedraw() +
  geom_boxplot(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = treatment), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), size = 3) +
  labs(y = "Labeled Predator : Labeled Prey abundance", x = "Days post-label") +
  scale_fill_manual(values = c("#FED439FF", "#709AE1FF")) +
  scale_color_manual(values = c("#FED439FF", "#709AE1FF")) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15),
        legend.position = "none") +
  facet_wrap(~plant)

# ---- Linear mixed-effects model (Supplemental Table 4) ---------------------
mod_abund <- lme4::lmer(
  predator_prey_ratio_abund ~ treatment * plant * as.factor(timepoint) + (1 | rep),
  data = incorp_preds3_abund
)
car::Anova(mod_abund) # Supplemental Table 4

treatment_signif_abund <- get_marginal_signif(mod_abund, "treatment")
plant_signif_abund <- get_marginal_signif(mod_abund, "plant")
y_annot_abund <- robust_upper(incorp_preds3_abund$predator_prey_ratio_abund) * 1.08

soil_background_abund <- strip_annotations(soil_background_abund) +
  annotate("richtext", x = 1.5, y = y_annot_abund, label = treatment_signif_abund,
           size = 6, fill = NA, label.color = NA)
plant_background_abund <- strip_annotations(plant_background_abund) +
  annotate("richtext", x = 1.5, y = y_annot_abund, label = plant_signif_abund,
           size = 6, fill = NA, label.color = NA)

# ---- Directional change over time: emtrends slope test ---------------------
mod_abund_trend <- lme4::lmer(
  predator_prey_ratio_abund ~ treatment * plant * timepoint + (1 | rep),
  data = incorp_preds3_abund
)

# Is the linear timepoint trend significant within each plant x treatment
# stratum? Bonferroni-corrected (not BH, unlike the pairwise brackets
# below) since a fitted trend line is a stronger claim and there are only
# 4 simultaneous tests.
get_sig_trends <- function(mod, group_vars = c("treatment", "plant"), method = "bonferroni") {
  trends <- emmeans::emtrends(mod, specs = group_vars, var = "timepoint")
  trend_df <- as.data.frame(summary(trends, infer = c(TRUE, TRUE))) %>%
    mutate(p.adj = p.adjust(p.value, method = method))
  message("trend tests -- raw p-values: ", paste(round(trend_df$p.value, 4), collapse = ", "),
          " | ", method, "-adjusted: ", paste(round(trend_df$p.adj, 4), collapse = ", "))
  trend_df %>% dplyr::filter(p.adj < 0.05)
}

# --- Helper: significance-bracket coordinates for all significant
# pairwise comparisons within each plant facet -- (a) across timepoints
# within the same treatment, and (b) across treatments within the same
# timepoint. Both types are pooled and BH-corrected together per plant
# facet, since both appear in the same panel.
build_sig_brackets <- function(mod, data, value_col, dodge_offset = 0.175) {
  trt_levels <- levels(factor(data$treatment))
  tp_levels  <- sort(unique(data$timepoint))

  plant_max_y <- data %>%
    group_by(plant) %>%
    summarize(base_y = robust_upper(.data[[value_col]]), .groups = "drop")

  # (a) within-treatment, across timepoints -- explicit, self-named
  # contrast weight vectors rather than parsing emmeans' auto-generated
  # contrast label text, which isn't reliable across emmeans versions
  emm_time <- emmeans::emmeans(mod, ~ timepoint | plant * treatment)
  tp_pairs <- utils::combn(tp_levels, 2, simplify = FALSE)
  contrast_list <- lapply(tp_pairs, function(pr) as.numeric(tp_levels == pr[1]) - as.numeric(tp_levels == pr[2]))
  names(contrast_list) <- vapply(tp_pairs, function(pr) paste(pr[1], "vs", pr[2]), character(1))

  brackets_time <- emmeans::contrast(emm_time, method = contrast_list, by = c("plant", "treatment"), adjust = "none") %>%
    as.data.frame() %>%
    mutate(
      t1 = as.numeric(sub(" vs .*", "", contrast)),
      t2 = as.numeric(sub(".* vs ", "", contrast)),
      offset = ifelse(treatment == trt_levels[1], -dodge_offset, dodge_offset),
      xmin = match(pmin(t1, t2), tp_levels) + offset,
      xmax = match(pmax(t1, t2), tp_levels) + offset,
      type = "time"
    ) %>%
    dplyr::select(plant, xmin, xmax, p.value, type)

  # (b) within-timepoint, across treatments
  emm_trt <- emmeans::emmeans(mod, ~ treatment | plant * timepoint)
  brackets_trt <- pairs(emm_trt, by = c("plant", "timepoint"), adjust = "none") %>%
    as.data.frame() %>%
    mutate(xpos = match(timepoint, tp_levels), xmin = xpos - dodge_offset, xmax = xpos + dodge_offset, type = "trt") %>%
    dplyr::select(plant, xmin, xmax, p.value, type)

  # pool raw p-values from both comparison types per plant facet, then
  # FDR-correct once across the combined family
  sig_brackets <- dplyr::bind_rows(brackets_time, brackets_trt) %>%
    group_by(plant) %>%
    mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
    ungroup() %>%
    dplyr::filter(p.adj < 0.05)

  sig_brackets %>%
    mutate(annotation = case_when(p.adj < 0.001 ~ "***", p.adj < 0.01 ~ "**", TRUE ~ "*")) %>%
    left_join(plant_max_y, by = "plant") %>%
    group_by(plant) %>%
    arrange(xmax - xmin, .by_group = TRUE) %>%
    # stack brackets within a facet so they don't overlap
    mutate(y_position = base_y * (1.03 + 0.09 * (dplyr::row_number() - 1)),
           tick_height = base_y * 0.015) %>%
    ungroup()
}

sig_trends_abund <- get_sig_trends(mod_abund_trend)
sig_brackets_abund <- build_sig_brackets(mod_abund, incorp_preds3_abund, "predator_prey_ratio_abund")

# only strata with a significant linear trend get a fitted regression line
trend_data_abund <- incorp_preds3_abund %>%
  inner_join(sig_trends_abund %>% dplyr::select(treatment, plant), by = c("treatment", "plant"))

dodge_offset <- 0.175 # 2 treatment levels dodged at width=0.7 -> half-offset = 0.175

ratios_time_abund <- ratios_time_abund_base +
  geom_smooth(data = trend_data_abund,
              aes(x = as.numeric(factor(timepoint)), y = predator_prey_ratio_abund, group = treatment, color = treatment),
              method = "lm", se = TRUE, linewidth = 1) +
  # brackets drawn with geom_segment()/geom_text() (not geom_signif()),
  # which don't reliably handle per-facet data subsetting
  geom_segment(data = sig_brackets_abund, aes(x = xmin, xend = xmax, y = y_position, yend = y_position),
               inherit.aes = FALSE, linewidth = 0.4, color = "black") +
  geom_segment(data = sig_brackets_abund, aes(x = xmin, xend = xmin, y = y_position, yend = y_position - tick_height),
               inherit.aes = FALSE, linewidth = 0.4, color = "black") +
  geom_segment(data = sig_brackets_abund, aes(x = xmax, xend = xmax, y = y_position, yend = y_position - tick_height),
               inherit.aes = FALSE, linewidth = 0.4, color = "black") +
  geom_text(data = sig_brackets_abund, aes(x = (xmin + xmax) / 2, y = y_position, label = annotation),
            inherit.aes = FALSE, vjust = -0.3, size = 6)

# Figure 3: combined A/B/C
fig3_ratios <- cowplot::plot_grid(soil_background_abund, plant_background_abund, ratios_time_abund,
                                   labels = "AUTO", nrow = 2)
fig3_ratios
ggsave("Fig3_predator_prey_ratios.png", fig3_ratios, width = 14, height = 11, units = "in", dpi = 600)

#### Relationships with soil carbon and 16S abundance (simple correlations) ###
# Merge covariates needed below
incorp_preds4_abund <- incorp_preds3_abund %>%
  left_join(
    rownames_to_column(meta_full, "sample") %>%
      dplyr::select(sample, perc_C, perc_N, C_N, delta_13C, copy_number_16S_g_soil_gwc) %>%
      distinct(),
    by = "sample"
  ) %>%
  mutate(log_16S = log10(copy_number_16S_g_soil_gwc))

cor.test(incorp_preds4_abund$perc_C, incorp_preds4_abund$predator_prey_ratio_abund, method = "spearman")
cor.test(incorp_preds4_abund$copy_number_16S_g_soil_gwc, incorp_preds4_abund$predator_prey_ratio_abund, method = "spearman")

#### Does each covariate itself differ by treatment/plant? (Supplemental Figure 3) ####
# Covariate as the OUTCOME: covariate ~ treatment*plant*timepoint + (1|rep)
test_covariate_as_outcome <- function(covariate, data = incorp_preds4_abund) {
  f <- as.formula(paste(covariate, "~ treatment * plant * as.factor(timepoint) + (1 | rep)"))
  mod <- lme4::lmer(f, data = data)
  cat("\n====", covariate, "(as outcome) ====\n")
  print(car::Anova(mod))
  invisible(mod)
}

mod_perc_C_outcome <- test_covariate_as_outcome("perc_C")
mod_perc_N_outcome <- test_covariate_as_outcome("perc_N")
mod_CN_outcome     <- test_covariate_as_outcome("C_N")
mod_16S_outcome    <- test_covariate_as_outcome("log_16S")
mod_delta_outcome  <- test_covariate_as_outcome("delta_13C")

extract_term <- function(mod, covariate, term_name) {
  tbl <- as.data.frame(car::Anova(mod))
  tbl$term <- rownames(tbl)
  tbl %>%
    dplyr::filter(term == term_name) %>%
    dplyr::rename(p.value = `Pr(>Chisq)`) %>%
    dplyr::mutate(covariate = covariate) %>%
    dplyr::select(covariate, term, Chisq, Df, p.value)
}

models <- list(perc_C = mod_perc_C_outcome, perc_N = mod_perc_N_outcome, C_N = mod_CN_outcome,
               log_16S = mod_16S_outcome, delta_13C = mod_delta_outcome)

# treatment, plant, and treatment:plant terms, each BH-corrected across
# the 5 covariates (5 genuinely distinct hypotheses per term)
build_effect_table <- function(term_name) {
  purrr::imap_dfr(models, ~ extract_term(.x, .y, term_name)) %>%
    mutate(p.adj = p.adjust(p.value, method = "BH"), signif = p_label(p.adj)) %>%
    arrange(p.adj)
}
treatment_effect_pvals <- build_effect_table("treatment")
plant_effect_pvals <- build_effect_table("plant")
interaction_effect_pvals <- build_effect_table("treatment:plant")

print(treatment_effect_pvals)
print(plant_effect_pvals)
print(interaction_effect_pvals)

# Boxplots (treatment on x, plant as fill), each panel annotated with all
# three terms
plot_covariate_boxplot_combined <- function(covariate, ylab, data = incorp_preds4_abund) {
  p_trt   <- treatment_effect_pvals$p.adj[treatment_effect_pvals$covariate == covariate][1]
  p_plant <- plant_effect_pvals$p.adj[plant_effect_pvals$covariate == covariate][1]
  p_int   <- interaction_effect_pvals$p.adj[interaction_effect_pvals$covariate == covariate][1]

  label_trt   <- paste0("Soil: ", bold_if_sig(p_label(p_trt), p_trt))
  label_plant <- paste0("Plant: ", bold_if_sig(p_label(p_plant), p_plant))
  label_int   <- paste0("Soil", "\u00d7", "Plant: ", bold_if_sig(p_label(p_int), p_int))

  data <- data %>% mutate(treatment = ifelse(treatment == "conv", "Conventional", "Organic"))

  # three separate single-line annotations (not one multi-line block) so
  # all three share the same x/hjust and left-align consistently
  ggplot(data, aes(x = treatment, y = .data[[covariate]], color = plant, fill = plant)) +
    theme_linedraw() +
    geom_boxplot(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.6, outlier.shape = NA, color = "black") +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7), size = 2) +
    annotate("richtext", x = -Inf, y = Inf, label = label_trt, hjust = -0.05, vjust = 1.3, size = 4, fill = NA, label.color = NA) +
    annotate("richtext", x = -Inf, y = Inf, label = label_plant, hjust = -0.05, vjust = 2.6, size = 4, fill = NA, label.color = NA) +
    annotate("richtext", x = -Inf, y = Inf, label = label_int, hjust = -0.05, vjust = 3.9, size = 4, fill = NA, label.color = NA) +
    labs(x = NULL, y = ylab, fill = "Plant", color = "Plant") +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15), legend.position = "none") +
    scale_color_manual(values = c("#7570B3FF", "#1B9E77FF")) +
    scale_fill_manual(values = c("#7570B3FF", "#1B9E77FF"))
}

perc_C_box   <- plot_covariate_boxplot_combined("perc_C", "% Carbon")
perc_N_box   <- plot_covariate_boxplot_combined("perc_N", "% Nitrogen")
CN_box       <- plot_covariate_boxplot_combined("C_N", "C:N ratio")
copy16S_box  <- plot_covariate_boxplot_combined("log_16S", expression("log"[10]~"16S copies g"^{-1}*" soil"))
delta_box    <- plot_covariate_boxplot_combined("delta_13C", expression(delta^13*"C"))

# unified legend, extracted from a throwaway plot with the same plant mapping
legend_source_plot <- ggplot(incorp_preds4_abund, aes(x = plant, y = perc_C, color = plant, fill = plant)) +
  geom_boxplot() +
  labs(fill = "Plant", color = "Plant") +
  scale_color_manual(values = c("#7570B3FF", "#1B9E77FF")) +
  scale_fill_manual(values = c("#7570B3FF", "#1B9E77FF")) +
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 13))
covariate_legend <- cowplot::get_legend(legend_source_plot)

fig_covariates_combined <- cowplot::plot_grid(
  perc_C_box, perc_N_box, CN_box, copy16S_box, delta_box, covariate_legend,
  labels = c("A", "B", "C", "D", "E", ""), ncol = 3, nrow = 2
)
fig_covariates_combined
ggsave("SuppFig3_covariates_treatment_by_plant.png", fig_covariates_combined,
       width = 16, height = 9, units = "in", dpi = 600)

#### Does the ratio relate to these covariates, controlling for design? ######
# ratio ~ covariate*plant*treatment*timepoint + (1|rep); the covariate's
# own main-effect term tests whether it relates to the ratio after
# accounting for plant, treatment, and time since labeling.
test_covariate_interaction <- function(covariate, data = incorp_preds4_abund) {
  f <- as.formula(paste("predator_prey_ratio_abund ~", covariate, "* plant * treatment * as.factor(timepoint) + (1 | rep)"))
  mod <- lme4::lmer(f, data = data)
  cat("\n====", covariate, "====\n")
  print(car::Anova(mod))
  invisible(mod)
}

mod_perc_C_interact <- test_covariate_interaction("perc_C")
mod_perc_N_interact <- test_covariate_interaction("perc_N")
mod_CN_interact     <- test_covariate_interaction("C_N")
mod_16S_interact    <- test_covariate_interaction("log_16S")
mod_delta_interact  <- test_covariate_interaction("delta_13C")
