# Figure 4 -- Step 2: aggregation + meta-regression + plot
#
# Pools the per-timepoint DESeq2 results (fig4_per_timepoint_deseq.R
# output, 4 CSVs) across the 4 timepoints, aggregates Order-level LFC
# into taxonomic trophic groups via inverse-variance weighting, then fits
# a hierarchical mixed-effects meta-regression (metafor::rma.mv) to
# estimate the rhizosphere effect for each trophic group, separately for
# each plant. 

library(tidyverse)
library(metafor)
library(cowplot)
library(ggtext)

# ---- 1. Compile and aggregate to trophic-group level -----------------------
org_wheat  <- read_csv("org_wheat_results_order_level.csv")
org_vetch  <- read_csv("org_vetch_results_order_level.csv")
conv_wheat <- read_csv("conv_wheat_results_order_level.csv")
conv_vetch <- read_csv("conv_vetch_results_order_level.csv")
all_results <- rbind(org_wheat, org_vetch, conv_wheat, conv_vetch)

# Order-level LFC -> trophic group, inverse-variance weighted mean
# (Supplemental Methods Eq. 1-2, same approach as the meta-analysis).
master_aggregated <- all_results %>%
  mutate(Unit = case_when(
    grepl("Myxococcota", PlotLabel) ~ "Myxococcota",
    grepl("Streptomycetales", AnalysisGroup) ~ "Streptomycetales",
    grepl("Bdellovibrionota", PlotLabel) ~ "Bdellovibrionota",
    grepl("Lysobacter", AnalysisGroup) ~ "Lysobacter",
    grepl("Cytophagales", AnalysisGroup) ~ "Cytophagales",
    grepl("Vampirovibrionales", AnalysisGroup) ~ "Vampirovibrionales",
    TRUE ~ "Non-Predator"
  )) %>%
  group_by(pred_type, Unit, timepoint, plant, treatment) %>%
  summarise(yi = sum(log2FoldChange / vi) / sum(1 / vi), V = 1 / sum(1 / vi), .groups = "drop")

# Cap precision weights at 500 to prevent high-depth datasets from
# dominating the model; sensitivity check across thresholds confirms
# stability (Supplemental Methods).
compare_caps <- function(cap_val, data) {
  dat <- data %>% mutate(V_temp = 1 / pmin(1 / V, cap_val))
  mod <- rma.mv(yi ~ Unit - 1, V = V_temp, random = list(~ 1 | timepoint), data = dat)
  as.vector(mod$beta)
}
cap_sensitivity <- sapply(c(100, 300, 500, 1000, 10000), compare_caps,
                          data = master_aggregated %>% mutate(Unit = relevel(factor(Unit), ref = "Non-Predator")))
max(abs(cap_sensitivity - cap_sensitivity[, 3])) # max deviation from the chosen cap (500)

master_weighted <- master_aggregated %>%
  mutate(W = 1 / V) %>%
  filter(!is.na(yi) & !is.infinite(W)) %>%
  mutate(V_stable = 1 / pmin(W, 500),
         Unit = relevel(factor(Unit), ref = "Non-Predator"))

# ---- 2. Hierarchical meta-regression (metafor::rma.mv) ---------------------
# Non-predator as reference level -- for the "relative to non-predators" panels
res_meta_final <- rma.mv(yi = yi, V = V_stable, mods = ~ Unit,
                         random = list(~ 1 | timepoint), data = master_weighted, method = "ML")

# Does adding plant/soil background improve the model? (reported in Results)
res_meta_plant <- rma.mv(yi = yi, V = V_stable, mods = ~ Unit + plant,
                         random = list(~ 1 | timepoint), data = master_weighted, method = "ML")
res_meta_plant2 <- rma.mv(yi = yi, V = V_stable, mods = ~ plant,
                          random = list(~ 1 | timepoint), data = master_weighted, method = "ML")
anova(res_meta_plant, res_meta_plant2)

res_meta_soil <- rma.mv(yi = yi, V = V_stable, mods = ~ Unit + treatment,
                        random = list(~ 1 | timepoint), data = master_weighted, method = "ML")
anova(res_meta_final, res_meta_soil)

# ---- 3. Figure 4: forest plots, per plant -----------------------------------
strategy_map <- data.frame(
  Unit = c("Cytophagales", "Streptomycetales", "Lysobacter", "Bdellovibrionota", "Vampirovibrionales", "Myxococcota", "Non-Predator"),
  pred_type = c("Facultative", "Facultative", "Facultative", "Obligate", "Obligate", "Facultative", "Non-Predator")
)
display_order <- c("Streptomycetales", "Myxococcota", "Bdellovibrionota", "Vampirovibrionales", "Cytophagales", "Lysobacter")

sig_label <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", p < 0.1 ~ ".", TRUE ~ "")

# Panel C/D: enrichment relative to non-predators, one model per plant
make_plant_forest <- function(plant_name, data = master_weighted) {
  plant_data <- data %>%
    filter(plant == plant_name) %>%
    mutate(Unit = relevel(droplevels(factor(Unit)), ref = "Non-Predator"))

  res_plant <- rma.mv(yi, V = V_stable, mods = ~ Unit, random = ~ 1 | timepoint, data = plant_data)

  level_order <- paste0("<i>", rev(display_order))

  plot_data <- data.frame(
    Unit = gsub("Unit", "", rownames(res_plant$beta)), LFC = as.vector(res_plant$beta),
    ci.lb = res_plant$ci.lb, ci.ub = res_plant$ci.ub, pval = res_plant$pval
  ) %>%
    filter(Unit != "intrcpt") %>%
    inner_join(strategy_map, by = "Unit") %>%
    mutate(pred_type = factor(pred_type, levels = c("Obligate", "Facultative")),
           Unit = paste0("<i>", Unit), Unit_label = factor(Unit, levels = level_order),
           sig_stars = sig_label(pval),
           star_pos = ifelse(LFC >= 0, ci.ub + 0.01, ci.lb - 0.01),
           star_hjust = ifelse(LFC >= 0, 0, 1))

  ggplot(plot_data, aes(x = LFC, y = Unit_label, fill = pred_type)) +
    geom_col() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 2) +
    geom_errorbarh(aes(xmin = ci.lb, xmax = ci.ub), width = 0.1) +
    geom_point(shape = 18, size = 5, color = "black", show.legend = FALSE) +
    geom_text(aes(x = star_pos, label = sig_stars, hjust = star_hjust), vjust = 0.5, size = 8, color = "black") +
    scale_fill_manual(values = c("Obligate" = "#2c7bb6", "Facultative" = "#d7191c")) +
    theme_linedraw() +
    theme(axis.title = element_text(size = 15), axis.text.y = element_markdown(size = 15),
          axis.text.x = element_text(size = 15), panel.grid.minor = element_blank(),
          legend.title = element_blank(), legend.text = element_text(size = 15),
          legend.position = c(0.99, 0.01), legend.justification = c(1, 0),
          legend.background = element_rect(color = "black", fill = "white")) +
    labs(x = "Rhizosphere effect relative to non-predators", y = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.15)), limits = c(-0.5, 0.55),
                       breaks = c(-0.5, -0.25, 0, 0.25, 0.5)) +
    scale_y_discrete(drop = FALSE)
}

# Panel A/B: enrichment relative to bulk soil (no-intercept model), per plant
make_global_forest <- function(plant_name, data = master_weighted) {
  plot_dat <- data %>% filter(plant == plant_name) %>% mutate(Unit = droplevels(factor(Unit)))
  res_global <- rma.mv(yi, V = V_stable, mods = ~ Unit - 1, random = ~ 1 | timepoint, data = plot_dat)

  level_order <- paste0("<i>", rev(c(display_order, "Non-Predator")))
  strategy_map_all <- strategy_map %>% mutate(pred_type = ifelse(pred_type == "Non-Predator", "None", pred_type))

  plot_data <- data.frame(
    Unit = gsub("Unit", "", rownames(res_global$beta)), LFC = as.vector(res_global$beta),
    ci.lb = res_global$ci.lb, ci.ub = res_global$ci.ub, pval = res_global$pval
  ) %>%
    left_join(strategy_map_all, by = "Unit") %>%
    mutate(pred_type = factor(pred_type, levels = c("Obligate", "Facultative", "None")),
           Unit = paste0("<i>", Unit), Unit_label = factor(Unit, levels = level_order),
           sig_stars = sig_label(pval),
           star_pos = ifelse(LFC >= 0, ci.ub + 0.01, ci.lb - 0.01),
           star_hjust = ifelse(LFC >= 0, 0, 1))

  ggplot(plot_data, aes(x = LFC, y = Unit_label, fill = pred_type)) +
    geom_col() +
    geom_vline(xintercept = 0, color = "black") +
    geom_errorbarh(aes(xmin = ci.lb, xmax = ci.ub), width = 0.1, color = "black", show.legend = FALSE) +
    geom_point(shape = 18, size = 5, color = "black", show.legend = FALSE) +
    geom_text(aes(x = star_pos, label = sig_stars, hjust = star_hjust), vjust = 0.5, size = 8, color = "black") +
    scale_fill_manual(values = c("Obligate" = "#1F78B4", "Facultative" = "#E31A1C", "None" = "darkgrey"),
                      labels = c("Obligate", "Facultative", "Non-Predator")) +
    theme_linedraw() +
    theme(axis.title = element_text(size = 15), axis.text.y = element_markdown(size = 15),
          axis.text.x = element_text(size = 15), panel.grid.minor = element_blank(),
          legend.title = element_blank(), legend.text = element_text(size = 15),
          legend.position = c(0.99, 0.01), legend.justification = c(1, 0),
          legend.background = element_rect(color = "black", fill = "white")) +
    labs(x = "Rhizosphere effect", y = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.15)), limits = c(-0.5, 0.55),
                       breaks = c(-0.5, -0.25, 0, 0.25, 0.5)) +
    scale_y_discrete(drop = FALSE)
}

vetch_plot_overall <- make_global_forest("vetch")
wheat_plot_overall <- make_global_forest("wheat")
vetch_plot <- make_plant_forest("vetch")
wheat_plot <- make_plant_forest("wheat")

# Figure 4: (A) vetch no-intercept, (B) wheat no-intercept,
# (C) vetch relative-to-non-predator, (D) wheat relative-to-non-predator
fig4 <- cowplot::plot_grid(vetch_plot_overall, wheat_plot_overall, vetch_plot, wheat_plot, labels = "AUTO")
fig4
ggsave("Figure4_rhizosphere_effect.png", fig4, dpi = 600, height = 6, width = 12, units = "in")
