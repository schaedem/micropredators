# Figure 1, Supplemental Figure 1, Supplemental Tables 2-3.
#
# Models the day of first 13C-label detection (an ordinal outcome: 1, 3,
# 7, or 14) as a function of treatment, plant host, and predator
# taxonomic group, via a cumulative link (proportional-odds ordinal
# logistic regression) model. Two versions of the trophic grouping are
# fit: a 7-level taxonomic grouping (CustomGroup; Supplemental Table 3)
# and a simpler 3-level Obligate/Facultative/Non-Predator grouping
# (pred_group; Supplemental Table 2). Predicted probabilities and
# bootstrapped SEs from the CustomGroup interaction model are plotted
# two ways: marginalized over treatment x plant (Figure 1) and shown
# separately for each treatment x plant combination (Supplemental
# Figure 1).

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Marie/Postdoc_Cornell/SOM_Greenhouse/final_data")

library(tidyverse)
library(boot)
library(ggtext)
library(ordinal)
library(emmeans)
library(multcomp)

# ---- 1. Data preparation ------------------------------------------------
incorps_soil <- readRDS("incorp_soil_adj_abun.RDS")

l2fc <- readRDS("joint_incorps_all.RDS") %>%
  mutate(padj_adj = p.adjust(padj, method = "BH")) %>%
  mutate(timepoint = ifelse(timepoint == "A", 1,
                            ifelse(timepoint == "B", 3,
                                   ifelse(timepoint == "C", 7, 14)))) %>%
  mutate(timepoint = as.numeric(timepoint)) %>%
  mutate(predator = ifelse(
    grepl("Bdello", Phylum) |
      grepl("Cytophag", Order) |
      grepl("Lyso", Genus) |
      grepl("Myxo", Phylum) |
      grepl("Streptomycetales", Order) |
      grepl("Vampiro", Order), "Predator", "Non-predator")) %>%
  mutate(pred_group = ifelse(
    grepl("Bdello", Phylum) | grepl("Vampiro", Order), "Obligate",
    ifelse(predator == "Predator" &
             !grepl("Bdello", Phylum) &
             !grepl("Vampiro", Order), "Facultative", "Non-Predator"))) %>%
  filter(ASV %in% incorps_soil$ASV) %>%
  group_by(ASV, treatment, plant, pred_group) %>%
  mutate(first_label = min(timepoint)) %>%
  ungroup() %>%
  mutate(timepoint   = as.factor(timepoint),
         first_label = as.factor(first_label)) %>%
  mutate(pred_group = factor(pred_group,
                             levels = c("Non-Predator", "Facultative", "Obligate"))) %>%
  mutate(CustomGroup = case_when(
    Phylum == "Bdellovibrionota"  ~ "*Bdellovibrionota*",
    Phylum == "Myxococcota"       ~ "*Myxococcota*",
    Genus  == "Lysobacter"        ~ "*Lysobacter*",
    Order  == "Streptomycetales"  ~ "*Streptomycetales*",
    Order  == "Cytophagales"      ~ "*Cytophagales*",
    Order  == "Vampirovibrionales"~ "*Vampirovibrionales*",
    TRUE                          ~ "Non-predator"
  )) %>%
  dplyr::select(ASV, CustomGroup, first_label, pred_group, treatment, plant) %>%
  unique()

l2fc$CustomGroup <- as.factor(l2fc$CustomGroup)
l2fc$treatment   <- as.factor(l2fc$treatment)
l2fc$plant       <- as.factor(l2fc$plant)
l2fc$first_label <- factor(l2fc$first_label, levels = c("1", "3", "7", "14"))

# ---- 2. Model fitting (Supplemental Tables 2-3) -------------------------

# CustomGroup (7-level taxonomic grouping) -- Supplemental Table 3
clm_interaction <- clm(first_label ~ treatment * plant * CustomGroup, data = l2fc)
clm_final       <- clm(first_label ~ treatment + plant + CustomGroup, data = l2fc)
anova(clm_final, clm_interaction) # model selection: additive vs. full interaction
anova(clm_interaction)            # Supplemental Table 3

# pred_group (3-level Obligate/Facultative/Non-Predator grouping) -- Supplemental Table 2
clm_interaction2 <- clm(first_label ~ treatment * plant * pred_group, data = l2fc)
clm_pred_group   <- clm(first_label ~ treatment + plant + pred_group, data = l2fc)
anova(clm_interaction2, clm_pred_group)
anova(clm_interaction2) # Supplemental Table 2

# ---- 3. Post-hoc contrasts (CustomGroup model) --------------------------

# Global contrasts, marginalizing over treatment x plant
emm_global <- emmeans(clm_interaction, ~ CustomGroup)
pairs(emm_global, adjust = "BH")

cld_int <- cld(emm_global, adjust = "BH", Letters = letters)
print(cld_int)

# compact-letter-display labels for the legend, from the contrasts above
cld_labels_int <- c(
  "*Bdellovibrionota*"   = "*Bdellovibrionota* (ab)",
  "*Cytophagales*"       = "*Cytophagales* (b)",
  "*Lysobacter*"         = "*Lysobacter* (ab)",
  "*Myxococcota*"        = "*Myxococcota* (a)",
  "*Streptomycetales*"   = "*Streptomycetales* (a)",
  "*Vampirovibrionales*" = "*Vampirovibrionales* (ab)",
  "Non-predator"         = "Non-predator (a)"
)

# Per-facet contrasts (differ from the global ones only in the interaction model)
emm_facet   <- emmeans(clm_interaction, ~ CustomGroup | treatment + plant)
pairs_facet <- pairs(emm_facet, adjust = "fdr")
summary(pairs_facet)

# Contrasts vs. non-predator specifically
contrast_vs_ref <- contrast(emm_facet, method = "trt.vs.ctrl", ref = "Non-predator", adjust = "fdr")
summary(contrast_vs_ref)

# ============================================================
# SUPPLEMENTAL FIGURE 1: predicted probabilities per treatment x plant
# ============================================================

# ---- 4. Bootstrap SEs for the per-facet (interaction) predictions -------
newdat <- expand.grid(
  CustomGroup = levels(l2fc$CustomGroup),
  treatment   = levels(l2fc$treatment),
  plant       = levels(l2fc$plant)
)
n_rows       <- nrow(newdat)
n_timepoints <- nlevels(l2fc$first_label)

boot_int_fn <- function(data, indices) {
  d <- data[indices, ]
  if (nlevels(droplevels(d$CustomGroup)) < nlevels(l2fc$CustomGroup)) return(rep(NA, n_rows * n_timepoints))
  if (nlevels(droplevels(d$treatment))   < nlevels(l2fc$treatment))   return(rep(NA, n_rows * n_timepoints))
  if (nlevels(droplevels(d$plant))       < nlevels(l2fc$plant))       return(rep(NA, n_rows * n_timepoints))
  if (nlevels(droplevels(d$first_label)) < nlevels(l2fc$first_label)) return(rep(NA, n_rows * n_timepoints))

  fit <- tryCatch(
    clm(first_label ~ treatment + plant + CustomGroup + treatment:plant:CustomGroup,
        nominal = ~ treatment + plant, data = d),
    error = function(e) NULL
  )
  if (is.null(fit)) return(rep(NA, n_rows * n_timepoints))

  probs <- tryCatch(
    suppressWarnings(predict(fit, newdata = newdat, type = "prob")$fit),
    error = function(e) NULL
  )
  if (is.null(probs)) return(rep(NA, n_rows * n_timepoints))
  as.vector(t(probs))
}

set.seed(123)
boot_int <- boot(data = l2fc, statistic = boot_int_fn, R = 500)
cat("Interaction model bootstrap failure rate:",
    round(mean(is.na(boot_int$t[, 1])) * 100, 1), "%\n")

# ---- 5. Predicted probabilities + SEs, long format ----------------------
pred_probs <- predict(clm_interaction, newdata = newdat, type = "prob")$fit
colnames(pred_probs) <- levels(l2fc$first_label)
se_main <- apply(boot_int$t, 2, sd, na.rm = TRUE)

pred_probs_long <- data.frame(
  newdat[rep(1:n_rows, each = n_timepoints), ],
  first_label = rep(as.numeric(levels(l2fc$first_label)), times = n_rows),
  prob        = as.vector(t(pred_probs)),
  SE          = se_main
) %>%
  mutate(
    # Bdellovibrionota in organic wheat has only n=3 observations -- SE
    # dropped for that cell rather than shown with misleadingly low
    # uncertainty
    low_confidence = treatment == "org" & plant == "wheat" & CustomGroup == "*Bdellovibrionota*",
    SE = ifelse(low_confidence, NA, SE),
    pred_group = case_when(
      CustomGroup == "*Bdellovibrionota*"   ~ "Obligate predator",
      CustomGroup == "*Vampirovibrionales*" ~ "Obligate predator",
      CustomGroup %in% c("*Cytophagales*", "*Lysobacter*", "*Myxococcota*", "*Streptomycetales*") ~ "Facultative predator",
      TRUE ~ "Non-predator"
    ),
    CustomGroup_label = recode(CustomGroup, !!!cld_labels_int),
    treatment_label = recode(treatment, "conv" = "Conventional", "org" = "Organic"),
    plant_label = recode(plant, "vetch" = "Vetch", "wheat" = "Wheat")
  )

# ---- 6. Plot -------------------------------------------------------------
group_colors <- c(
  "Obligate predator"    = "#1F78B4",
  "Facultative predator" = "#E31A1C",
  "Non-predator"         = "darkgrey"
)

interaction_plot <-
  ggplot(pred_probs_long,
       aes(x = first_label, y = prob * 100, ymin = (prob - SE) * 100, ymax = (prob + SE) * 100,
           group = CustomGroup, color = pred_group, shape = CustomGroup, alpha = low_confidence)) +
  geom_line(linewidth = 0.5, position = position_dodge(width = 1)) +
  geom_point(size = 4, position = position_dodge(width = 1)) +
  geom_errorbar(width = 0.3, position = position_dodge(width = 1)) +
  facet_grid(plant_label ~ treatment_label) +
  scale_x_continuous(breaks = as.numeric(levels(l2fc$first_label))) +
  scale_color_manual(values = group_colors) +
  scale_shape_manual(values = c(
    "*Bdellovibrionota*"   = 16, "*Cytophagales*" = 17, "*Lysobacter*" = 15,
    "*Myxococcota*" = 18, "*Streptomycetales*" = 8, "Non-predator" = 1
  )) +
  scale_alpha_manual(values = c("FALSE" = 1, "TRUE" = 0.6), guide = "none") +
  labs(x = "Day", y = expression("Probability (%) of first " ^13 * "C label"),
       color = "Predator group", shape = "Trophic group") +
  theme_linedraw() +
  theme(axis.title = element_text(size = 15), strip.text = element_text(size = 15),
        axis.text = element_text(size = 12), legend.text = ggtext::element_markdown(size = 12),
        legend.title = element_blank())
interaction_plot
ggsave("SuppFig1_label_probability_interaction.png", width = 10, height = 7, units = "in", dpi = 600)

# ============================================================
# FIGURE 1: predicted probabilities marginalized over treatment x plant
# ============================================================

# ---- 7. Bootstrap SEs, averaging predictions over treatment x plant ----
newdat_marginal <- data.frame(CustomGroup = levels(l2fc$CustomGroup))
n_rows_marginal       <- nrow(newdat_marginal)
n_timepoints_marginal <- nlevels(l2fc$first_label)

boot_marginal_fn <- function(data, indices) {
  d <- data[indices, ]
  if (nlevels(droplevels(d$CustomGroup)) < nlevels(l2fc$CustomGroup)) return(rep(NA, n_rows_marginal * n_timepoints_marginal))
  if (nlevels(droplevels(d$treatment))   < nlevels(l2fc$treatment))   return(rep(NA, n_rows_marginal * n_timepoints_marginal))
  if (nlevels(droplevels(d$plant))       < nlevels(l2fc$plant))       return(rep(NA, n_rows_marginal * n_timepoints_marginal))
  if (nlevels(droplevels(d$first_label)) < nlevels(l2fc$first_label)) return(rep(NA, n_rows_marginal * n_timepoints_marginal))

  fit <- tryCatch(
    clm(first_label ~ treatment + plant + CustomGroup + treatment:plant:CustomGroup,
        nominal = ~ treatment + plant, data = d),
    error = function(e) NULL
  )
  if (is.null(fit)) return(rep(NA, n_rows_marginal * n_timepoints_marginal))

  probs_full <- tryCatch(
    suppressWarnings(predict(fit, newdata = newdat, type = "prob")$fit),
    error = function(e) NULL
  )
  if (is.null(probs_full)) return(rep(NA, n_rows_marginal * n_timepoints_marginal))

  probs_df <- as.data.frame(probs_full)
  probs_df$CustomGroup <- newdat$CustomGroup
  probs_marginal <- probs_df %>%
    group_by(CustomGroup) %>%
    summarise(across(everything(), mean), .groups = "drop") %>%
    arrange(match(CustomGroup, levels(l2fc$CustomGroup))) %>%
    dplyr::select(-CustomGroup)

  as.vector(t(as.matrix(probs_marginal)))
}

set.seed(123)
boot_marginal <- boot(data = l2fc, statistic = boot_marginal_fn, R = 500)
cat("Marginal bootstrap failure rate:",
    round(mean(is.na(boot_marginal$t[, 1])) * 100, 1), "%\n")

# ---- 8. Predicted probabilities + SEs, long format -----------------------
probs_full_orig <- predict(clm_interaction, newdata = newdat, type = "prob")$fit
colnames(probs_full_orig) <- levels(l2fc$first_label)

probs_full_df <- as.data.frame(probs_full_orig)
probs_full_df$CustomGroup <- newdat$CustomGroup

probs_marginal_orig <- probs_full_df %>%
  group_by(CustomGroup) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  arrange(match(CustomGroup, levels(l2fc$CustomGroup)))

se_marginal <- apply(boot_marginal$t, 2, sd, na.rm = TRUE)

pred_probs_marginal <- data.frame(
  CustomGroup = rep(levels(l2fc$CustomGroup), each = n_timepoints_marginal),
  first_label = rep(as.numeric(levels(l2fc$first_label)), times = n_rows_marginal),
  prob        = as.vector(t(as.matrix(probs_marginal_orig %>% dplyr::select(-CustomGroup)))),
  SE          = se_marginal
) %>%
  mutate(
    pred_group = case_when(
      CustomGroup == "*Bdellovibrionota*"   ~ "Obligate predator",
      CustomGroup == "*Vampirovibrionales*" ~ "Obligate predator",
      CustomGroup %in% c("*Cytophagales*", "*Lysobacter*", "*Myxococcota*", "*Streptomycetales*") ~ "Facultative predator",
      TRUE ~ "Non-predator"
    ),
    CustomGroup_label = recode(CustomGroup, !!!cld_labels_int)
  )

# ---- 9. Plot --------------------------------------------------------------
ggplot(pred_probs_marginal,
       aes(x = first_label, y = prob, ymin = prob - SE, ymax = prob + SE,
           group = CustomGroup_label, color = pred_group, shape = CustomGroup_label)) +
  geom_line(linewidth = 0.5, position = position_dodge(width = 0.5)) +
  geom_point(size = 5, position = position_dodge(width = 0.5)) +
  geom_errorbar(width = 0.3, position = position_dodge(width = 0.5)) +
  scale_x_continuous(breaks = as.numeric(levels(l2fc$first_label))) +
  scale_color_manual(values = group_colors) +
  scale_shape_manual(values = c(
    "*Bdellovibrionota* (ab)" = 16, "*Cytophagales* (b)" = 17, "*Lysobacter* (ab)" = 15,
    "*Myxococcota* (a)" = 18, "*Streptomycetales* (a)" = 8, "Non-predator (a)" = 1
  )) +
  labs(x = "Day", y = expression("Probability of first " ^13 * "C label"),
       color = "Trophic group", shape = "Predator group") +
  guides(color = guide_legend(order = 1, ncol = 1), shape = guide_legend(order = 2, ncol = 1)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12),
        legend.text = ggtext::element_markdown(size = 12), legend.title = element_text(size = 12, face = "bold"),
        legend.position = "inside", legend.position.inside = c(0.65, 0.011),
        legend.justification = c(1, 0), legend.direction = "vertical",
        legend.box = "horizontal", legend.box.just = "top",
        legend.box.background = element_rect(color = "black", linewidth = 0.3),
        legend.box.margin = margin(5, 5, 5, 5))
ggsave("Fig1_label_probability_marginal.png", width = 10, height = 7, units = "in", dpi = 600)
