# Figure 2. Micropredator composition, abundance, richness, and targeted
# genus detection in the rhizosphere.
#
# (A) Bubble plot of the most abundant labeled taxa (top 10 predator ASVs,
#     top 5 non-predator ASVs per treatment x plant), ordered by
#     hierarchical clustering (binary distance on presence/absence across
#     treatment x plant x timepoint) within each of the two panels.
# (B) Labeled predator abundance (log-adjusted) in the rhizosphere.
# (C) Labeled predator richness in the rhizosphere, with ALL significant
#     pairwise (treatment x plant) comparisons shown as brackets.
# (D) Binary detection grid for chosen genera, built from ALL labeling
#     events in `all_incorps` (not just the top-ASV display subset),
#     restricted to ASVs present in soil_incorp_rare/final_data.
# Panels combined via patchwork; saved as one figure.

library(tidyverse)
library(phyloseq)
library(viridis)
library(lme4)
library(lmerTest)
library(ggsignif)
library(ggpubr)
library(emmeans)
library(Polychrome)
library(car)
library(rstatix)
library(patchwork)
library(ggdendro)
library(ggtext)
library(colorspace)
library(cowplot)
library(ggnewscale)

# 1. LOAD DATA --------------------------------------------------------------
org_incorp <- readRDS("source_incorp_by_day.RDS")
conv_incorp <- readRDS("conv_source_incorp_by_day.RDS")
all_incorps <- rbind(org_incorp, conv_incorp)

#soil_incorp_rare <- readRDS("nifa_joint.phyloseq.RDS")
soil_incorp_rare <- readRDS("joint_incorps_soil.RDS")

org_pap <- readRDS("nifa.paprica.RDS")
conv_pap <- readRDS("nifa_conv.paprica.RDS")
pap_ref <- bind_rows(org_pap, conv_pap) %>%
  distinct(ASV, .keep_all = TRUE) %>%
  dplyr::select(ASV, n16S)

# 2. FILTERING KEYS -----------------------------------------------------------
incorp_keys <- all_incorps %>%
  mutate(key = paste(ASV, treatment, timepoint, plant, sep = "_")) %>%
  pull(key)

# 3. DATA PREP & NORMALIZATION ------------------------------------------------
rare_depth <- 23993

final_data <- psmelt(soil_incorp_rare) %>%
  dplyr::rename(ASV = OTU) %>%
  mutate(match_key = paste(ASV, treatment, timepoint, plant, sep = "_")) %>%
  filter(match_key %in% incorp_keys) %>%
  inner_join(pap_ref, by = "ASV") %>%
  mutate(
    rel_abund = Abundance / rare_depth,
    adj_abund = (rel_abund * (copy_number_16S_g_soil_gwc / 2)) / n16S,
    log_adj_abund = ifelse(adj_abund <= 0, 0, log10(adj_abund))
  ) %>%
  mutate(
    predator = ifelse(
      grepl("Bdellovib|Myxococc", Phylum) |
        grepl("Streptomycetales|Vampiro|Cytophag", Order) |
        grepl("Lysobacter", Genus) |
      grepl("Agromyces", Genus),
      "Predator", "Non-predator"),
    Taxa_Label = case_when(
      !is.na(Genus) & Genus != "unclassified" ~ Genus,
      !is.na(Family) & Family != "unclassified" ~ paste("f_", Family),
      !is.na(Order) & Order != "unclassified" ~ paste("o_", Order),
      TRUE ~ "Unclassified")
  )

final_data_rhizo <- final_data %>%
  filter(sample_type == "rhizo")

# 4. SUMMARIZE & SELECT TOP ASVs ----------------------------------------------
mean_plot_data <- final_data_rhizo %>%
  group_by(ASV, Taxa_Label, predator, treatment, plant, timepoint, Phylum, Order, Family, Genus) %>%
  summarize(mean_log_adj = mean(log_adj_abund, na.rm = TRUE), .groups = "drop")

top_asvs <- mean_plot_data %>%
  group_by(ASV, predator, treatment, plant) %>%
  summarize(max_mean = max(mean_log_adj), .groups = "drop") %>%
  group_split(predator) %>%
  map_dfr(~ {
    n_to_slice <- if (unique(.x$predator) == "Predator") 10 else 5
    .x %>%
      group_by(treatment, plant) %>%
      slice_max(max_mean, n = n_to_slice, with_ties = FALSE)
  }) %>%
  pull(ASV)

# 5. HIERARCHICAL CLUSTERING ORDER --------------------------------------------
binary_input <- all_incorps %>%
  mutate(sample_id = paste(treatment, plant, timepoint, sep = "_"),
         labeled = 1) %>%
  distinct(ASV, sample_id, labeled) %>%
  filter(ASV %in% top_asvs)

asv_predator <- final_data %>%
  distinct(ASV, predator)

wide_binary <- binary_input %>%
  pivot_wider(id_cols = ASV, names_from = sample_id,
              values_from = labeled, values_fill = 0) %>%
  left_join(asv_predator, by = "ASV")

cluster_results <- wide_binary %>%
  group_split(predator) %>%
  purrr::map(function(df) {
    mat <- as.matrix(dplyr::select(df, -ASV, -predator))
    rownames(mat) <- df$ASV
    grp <- unique(df$predator)
    if (nrow(mat) < 3) {
      return(list(predator = grp, hc = NULL, order = df$ASV))
    }
    d  <- dist(mat, method = "binary")
    hc <- hclust(d, method = "complete")
    list(predator = grp, hc = hc, order = df$ASV[hc$order])
  }) %>%
  purrr::set_names(purrr::map_chr(., "predator"))

ordered_asvs <- purrr::map(cluster_results, "order") %>% unlist(use.names = FALSE)

plot_subset <- mean_plot_data %>%
  filter(ASV %in% top_asvs) %>%
  mutate(
    ASV = factor(ASV, levels = ordered_asvs),
    timepoint = fct_recode(factor(as.character(timepoint)), "1"="A", "3"="B", "7"="C", "14"="D"),
    timepoint = fct_relevel(timepoint, "1", "3", "7", "14"),
    treatment = case_when(treatment == "org" ~ "Organic",
                          treatment == "conv" ~ "Conventional",
                          TRUE ~ treatment),
    plant = case_when(plant == "vetch" ~ "Vetch",
                      plant == "wheat" ~ "Wheat",
                      TRUE ~ plant)) %>%
  mutate(mean_log_adj = ifelse(mean_log_adj == 0, NA, mean_log_adj)) %>%
  mutate(Taxa_Label = ifelse(Taxa_Label == "uncultured", paste(Taxa_Label, Order, sep = " "), Taxa_Label)) %>%
  mutate(Order_key = paste(predator, Order, sep = "__"))

# 6. COLOR MAPPING & LABELS ----------------------------------------------------
make_distinct_hcl <- function(n, h_range, l_bands = c(80, 60, 38), c_bands = c(95, 65)) {
  hues   <- seq(h_range[1], h_range[2], length.out = n)
  l_vals <- rep(l_bands, length.out = n)
  c_vals <- rep(c_bands, length.out = n)
  colorspace::hex(colorspace::polarLUV(L = l_vals, C = c_vals, H = hues), fixup = TRUE)
}

strip_order_key <- function(x) sub("^(Predator|Non-predator)__", "", x)

#Set colors for preds and non-preds
nonpred_orders <- sort(unique(plot_subset$Order_key[plot_subset$predator == "Non-predator"]))
pred_orders    <- sort(unique(plot_subset$Order_key[plot_subset$predator == "Predator"]))

warm_colors <- make_distinct_hcl(length(nonpred_orders), h_range = c(10, 90))
cool_colors <- make_distinct_hcl(length(pred_orders),    h_range = c(185, 295))

order_color_map <- c(setNames(warm_colors, nonpred_orders),
                     setNames(cool_colors, pred_orders))

recode_order_label <- function(x) {
  # operates on the plain (already-stripped) Order name
  x <- gsub("^unclassified\\s+KD4-96$", "Chloroflexota unclassified KD4-96", x)
  x <- gsub("^unclassified\\s+bacteriap25$", "Myxococcota unclassified bacteriap25", x)
  x <- gsub("Xanthomonadales$", "Xanthomonadales, Lysobacter", x)
  x <- gsub("Micrococcales$", "Micrococcales, Agromyces", x)
}

show_taxa_labels <- FALSE

# color lookup now keyed on Order_key (unambiguous); visible text still
# uses the plain Order-derived Taxa_Label
asv_label_lookup <- plot_subset %>%
  distinct(ASV, Taxa_Label, Order, Order_key) %>%
  mutate(label_html = paste0(
    "<span style='color:", order_color_map[Order_key], "; font-weight:bold; font-style:italic'>",
    Taxa_Label, "</span>"
  ))

asv_labels <- setNames(asv_label_lookup$label_html, asv_label_lookup$ASV)

size_limits <- range(plot_subset$mean_log_adj, na.rm = TRUE)
size_breaks <- pretty(size_limits, n = 4)
size_breaks <- size_breaks[size_breaks >= size_limits[1] & size_breaks <= size_limits[2]]

# 7. MEDIAN FIRST-LABEL DAY (annotation, currently off) -----------------------
show_median_first_label <- FALSE

asv_predator_map <- distinct(final_data, ASV, predator)
dup_asvs <- asv_predator_map %>% count(ASV) %>% filter(n > 1)
if (nrow(dup_asvs) > 0) {
  warning(nrow(dup_asvs), " ASV(s) have inconsistent predator classification across samples: ",
          paste(dup_asvs$ASV, collapse = ", "))
}

unmatched_asvs <- setdiff(unique(all_incorps$ASV), unique(asv_predator_map$ASV))
if (length(unmatched_asvs) > 0) {
  warning(length(unmatched_asvs), " ASV(s) in all_incorps have no match in final_data.")
}

median_label_day <- all_incorps %>%
  select(-any_of(c("Order", "predator"))) %>%
  mutate(
    timepoint = case_when(
      timepoint == "A" ~ 1, timepoint == "B" ~ 3,
      timepoint == "C" ~ 7, timepoint == "D" ~ 14,
      TRUE ~ suppressWarnings(as.numeric(timepoint))
    ),
    treatment = case_when(treatment == "org" ~ "Organic",
                          treatment == "conv" ~ "Conventional",
                          TRUE ~ treatment)
  ) %>%
  left_join(asv_predator_map, by = "ASV") %>%
  group_by(ASV, plant, treatment, predator) %>%
  summarize(first_label = min(timepoint), .groups = "drop") %>%
  group_by(predator, treatment) %>%
  summarize(median_day = median(first_label, na.rm = TRUE), n_events = n(), .groups = "drop") %>%
  mutate(label_text = paste0("Median first label: day ", median_day))

# 8. VISUALIZATION: BUBBLE PLOT + DENDROGRAM -----------------------------------
make_dendro_plot <- function(hc) {
  if (is.null(hc)) return(plot_spacer())
  dd <- ggdendro::dendro_data(hc, type = "rectangle")
  ggplot() +
    geom_segment(data = dd$segments,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 linewidth = 0.4, color = "grey35") +
    scale_x_continuous(expand = expansion(add = 0.6)) +
    scale_y_reverse(expand = c(0, 0)) +
    coord_flip() +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

dendro_nonpred <- make_dendro_plot(cluster_results[["Non-predator"]]$hc)
dendro_pred    <- make_dendro_plot(cluster_results[["Predator"]]$hc)

build_group_bubble <- function(pred_grp, show_column_strip, show_x_axis,
                               fill_values, fill_breaks, legend_title) {
  df <- plot_subset %>% filter(predator == pred_grp)
  p <- ggplot(df, aes(x = timepoint, y = ASV, fill = Order_key)) +
    geom_point(aes(size = mean_log_adj, alpha = plant), shape = 21, position = position_dodge(0.5)) +
    facet_grid(predator ~ treatment) +
    scale_fill_manual(
      values = fill_values, breaks = fill_breaks,
      labels = recode_order_label(strip_order_key(fill_breaks)),
      name = legend_title,
      guide = guide_legend(
        override.aes = list(size = 5), ncol = 1,
        label.theme = element_text(face = "italic", size = 12)
      )
    ) +
    scale_size_continuous(
      range = c(1, 10), limits = size_limits, breaks = size_breaks,
      name = "Mean Log\nAdjusted Abundance"
    ) +
    scale_y_discrete(labels = if (show_taxa_labels) asv_labels else NULL, drop = TRUE) +
    scale_alpha_manual(
      values = c("Wheat" = 1, "Vetch" = 0.5), name = "Plant Host",
      guide = guide_legend(
        title.theme = element_text(size = 12, angle = 0),
        override.aes = list(color = "#3B3B3B", size = 5, shape = 16)
      )
    ) +
    theme_linedraw() +
    labs(x = if (show_x_axis) "Days since pulse label" else NULL, y = NULL) +
    theme(
      strip.background = element_rect(fill = "black"),
      strip.text = element_text(color = "white", size = 15),
      strip.text.x = if (show_column_strip) element_text(color = "white", size = 15) else element_blank(),
      strip.background.x = if (show_column_strip) element_rect(fill = "black") else element_blank(),
      axis.text.y = if (show_taxa_labels) ggtext::element_markdown(size = 12) else element_blank(),
      axis.text.x = if (show_x_axis) element_text(size = 12) else element_blank(),
      axis.ticks.x = if (show_x_axis) element_line() else element_blank(),
      axis.title.x = if (show_x_axis) element_text(size = 15) else element_blank(),
      axis.ticks.y = element_blank(),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      plot.margin = margin(t = 2, r = 5, b = 2, l = 0)
    )
  
  if (show_median_first_label) {
    p <- p + geom_text(
      data = median_label_day %>% filter(predator == pred_grp),
      aes(x = Inf, y = Inf, label = label_text),
      inherit.aes = FALSE, hjust = 1.05, vjust = 1.5, size = 3.2, fontface = "italic"
    )
  }
  
  p
}

bubble_nonpred <- build_group_bubble(
  "Non-predator", show_column_strip = TRUE, show_x_axis = FALSE,
  fill_values = order_color_map[nonpred_orders], fill_breaks = nonpred_orders,
  legend_title = "Non-Predators"
)
bubble_pred <- build_group_bubble(
  "Predator", show_column_strip = FALSE, show_x_axis = TRUE,
  fill_values = order_color_map[pred_orders], fill_breaks = pred_orders,
  legend_title = "Predators"
)

extract_legend <- function(p) {
  g <- ggplot2::ggplotGrob(p)
  idx <- grep("^guide-box", g$layout$name)
  if (length(idx) == 0) return(grid::nullGrob())
  g$grobs[[idx[1]]]
}

legend_plot <- ggplot(plot_subset, aes(x = timepoint, y = ASV)) +
  geom_point(
    data = filter(plot_subset, predator == "Non-predator"),
    aes(fill = Order_key, size = mean_log_adj, alpha = plant), shape = 21
  ) +
  scale_fill_manual(
    values = order_color_map[nonpred_orders], breaks = nonpred_orders,
    labels = recode_order_label(strip_order_key(nonpred_orders)), name = "Non-Predators",
    guide = guide_legend(order = 1, override.aes = list(size = 4), ncol = 1,
                         label.theme = element_text(face = "italic", size = 10.5))
  ) +
  ggnewscale::new_scale_fill() +
  geom_point(
    data = filter(plot_subset, predator == "Predator"),
    aes(fill = Order_key, size = mean_log_adj, alpha = plant), shape = 21
  ) +
  scale_fill_manual(
    values = order_color_map[pred_orders], breaks = pred_orders,
    labels = recode_order_label(strip_order_key(pred_orders)), name = "Predators",
    guide = guide_legend(order = 2, override.aes = list(size = 4), ncol = 1,
                         label.theme = element_text(face = "italic", size = 10.5))
  ) +
  scale_size_continuous(
    range = c(1, 10), limits = size_limits, breaks = size_breaks,
    name = "Mean Log\nAdjusted Abundance", guide = guide_legend(order = 3)
  ) +
  scale_alpha_manual(
    values = c("Wheat" = 1, "Vetch" = 0.5), name = "Plant Host",
    guide = guide_legend(
      order = 4, title.theme = element_text(size = 11, angle = 0),
      override.aes = list(color = "#3B3B3B", size = 4, shape = 16)
    )
  ) +
  theme_linedraw() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10.5),
    legend.key.size = unit(0.4, "cm")
  )

legend_grob <- extract_legend(legend_plot)

design <- "
ABE
CDE
"

n_pred    <- length(cluster_results[["Predator"]]$order)
n_nonpred <- length(cluster_results[["Non-predator"]]$order)

bubble_plot <- dendro_nonpred +
  (bubble_nonpred + theme(legend.position = "none")) +
  dendro_pred +
  (bubble_pred + theme(legend.position = "none")) +
  wrap_elements(full = legend_grob) +
  plot_layout(design = design,
              heights = c(n_nonpred, n_pred),
              widths = c(0.12, 1, 0.5))

bubble_plot
ggsave("bubble_abund_combined_hierarchical.png", height = 10, width = 12, units = "in", dpi = 600)

# 9. PREDATOR ABUNDANCE (Panel B) ----------------------------------------------
plot_df <- final_data %>%
  filter(predator == "Predator") %>%
  filter(sample_type == "rhizo") %>%
  filter(log_adj_abund > 0) %>%
  mutate(treatment = case_when(treatment == "conv" ~ "Conventional",
                               treatment == "org" ~ "Organic",
                               TRUE ~ treatment),
         plant = case_when(plant == "vetch" ~ "Vetch",
                           plant == "wheat" ~ "Wheat",
                           TRUE ~ plant))

abund_plot <-
  ggplot(plot_df, aes(x = plant, y = log_adj_abund, fill = treatment)) +
  geom_point(aes(alpha = plant, color = treatment), size = 5, position = position_jitterdodge()) +
  theme_linedraw() +
  scale_alpha_manual(values = c("Wheat" = 1, "Vetch" = 0.5), name = "Plant") +
  scale_color_manual(values = c("#FED439FF", "#709AE1FF")) +
  labs(x = "", y = "Log Adjusted\nLabeled Predator Abundance") +
  geom_boxplot(outliers = FALSE, alpha = 0.6, position = position_dodge()) +
  scale_fill_manual(values = c("#FED439FF", "#709AE1FF")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.16))) +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15),
        plot.margin = margin(t = 10, r = 10, b = 10, l = -10),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  ) +
  # NOTE: brackets left as manually-placed literals (unchanged) -- only
  # panel C was asked to switch to auto-generated significant brackets.
  geom_signif(
    y_position = c(9, 9),
    xmin = c(0.8, 1.8),
    xmax = c(1.2, 2.2),
    annotation = c("***", "***"),
    tip_length = 0.01
  ) +
  geom_signif(
    y_position = c(9.4),
    xmin = c(1.8),
    xmax = c(1.2),
    annotation = c("***"),
    tip_length = 0.01
  ) +
  geom_signif(
    y_position = c(9.8),
    xmin = c(0.8),
    xmax = c(2.2),
    annotation = c("***"),
    tip_length = 0.01
  ) +
  geom_signif(
    y_position = c(10.2),
    xmin = c(0.8),
    xmax = c(1.8),
    annotation = c("ns"),
    tip_length = 0.01
  ) +
  geom_signif(
    y_position = c(10.6),
    xmin = c(1.2),
    xmax = c(2.2),
    annotation = c("***"),
    tip_length = 0.01
  )

abund_plot
ggsave("pred_abund_plot.png", height = 3.5, width = 6, units = "in", dpi = 600)

# 10. PREDATOR RICHNESS (Panel C) ----------------------------------------------
plot_df2 <- final_data %>%
  filter(predator == "Predator") %>%
  filter(sample_type == "rhizo") %>%
  filter(log_adj_abund > 0) %>%
  mutate(treatment = case_when(treatment == "conv" ~ "Conventional",
                               treatment == "org" ~ "Organic",
                               TRUE ~ treatment),
         plant = case_when(plant == "vetch" ~ "Vetch",
                           plant == "wheat" ~ "Wheat",
                           TRUE ~ plant)) %>%
  group_by(plant, treatment, timepoint, rep, label) %>%
  summarize(total_abund = log10(sum(adj_abund)),
            total_richness = sum(length(unique(ASV)))) %>%
  ungroup()

test_df2 <- plot_df2 %>%
  mutate(rep = as.factor(timepoint), treatment = as.factor(treatment),
         plant = as.factor(plant), timepoint = as.factor(timepoint))
mod2 <- lmer(total_richness ~ treatment * plant + timepoint + (1 | rep), data = test_df2)
car::Anova(mod2)
summary(mod2)

richness_contrasts <- as.data.frame(emmeans(mod2, pairwise ~ treatment * plant)$contrasts) %>%
  mutate(padj = adjust_pvalue(p.value, method = "BH"),
         signif = cut(padj, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                      labels = c("***", "**", "*", "ns"), right = TRUE))
richness_contrasts

richness_plot <-
  ggplot(plot_df2, aes(x = plant, y = total_richness, fill = treatment)) +
  geom_jitter(aes(alpha = plant, color = treatment), size = 5, position = position_jitterdodge()) +
  theme_linedraw() +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05 + 0.1 * max(1, nrow(bracket_df))))) +
  scale_alpha_manual(values = c("Wheat" = 1, "Vetch" = 0.5), name = "Plant") +
  scale_color_manual(values = c("#FED439FF", "#709AE1FF")) +
  labs(x = "", y = "Labeled\nPredator Richness") +
  geom_boxplot(outliers = FALSE, alpha = 0.6, position = position_dodge()) +
  scale_fill_manual(values = c("#FED439FF", "#709AE1FF")) +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15),
        plot.margin = margin(t = 10, r = 10, b = 10, l = -10),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )
  geom_signif(
    y_position = c(15, 14.5),
    xmin = c(1.2, 1.8),
    xmax = c(1.8, 2.2),
    annotation = c("**", "**"),
    tip_length = 0.01
  ) 



richness_plot
ggsave("pred_richness_plot.png", height = 3.5, width = 6, units = "in", dpi = 600)

# PANEL D -----------------------------------------------------------------------

chosen_genus <- c("Haliangium", "Agromyces", "Streptomyces", "Lysobacter")


top_tax <- unique(plot_subset$Genus[!is.na(plot_subset$Genus) & plot_subset$Genus != "unclassified"])
chosen_genus <- intersect(chosen_genus, top_tax)

genus_detection <- plot_subset %>%
  mutate(Genus = ifelse(is.na(Genus) | Genus == "unclassified",
                        paste0(Order, " (order)"), Genus),
         Genus = ifelse(Genus == "uncultured", paste0(Order, " (uncultured) "), Genus)) %>%
  mutate(label = ifelse(is.na(mean_log_adj), 0, 1)) %>%
  filter(Genus %in% chosen_genus) %>%
  filter(label == 1) %>%                                   # only actually-labeled detections
  mutate(timepoint = as.numeric(as.character(timepoint))) %>%  # factor -> numeric for x axis
  distinct(plant, treatment, timepoint, Genus) %>%
  mutate(treatment_label = dplyr::recode(treatment, "conv" = "Conventional", "org" = "Organic"),
         plant = factor(plant, levels = c("Vetch", "Wheat")),
         Genus_label = paste0("*", Genus, "*"))

detection_grid_plot <- ggplot(genus_detection,
                              aes(x = timepoint, y = plant, fill = treatment_label)) +
  geom_point(shape = 21, color = "black", size = 7, position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("Conventional" = "#FED439FF", "Organic" = "#709AE1FF"),
                    name = "Soil background") +
  scale_x_continuous(breaks = c(1, 3, 7, 14), limits = c(-0.5, 15)) +
  labs(x = "Days post label", y = NULL) +
  facet_wrap(~Genus_label, nrow = 1) +
  theme_linedraw() +
  theme(axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 15),
        strip.text = element_markdown(size = 15),
        panel.grid.minor = element_blank(),
        legend.position = "none")

detection_grid_plot

plot_grid(detection_grid_plot, labels = c("D"), label_size = 20)
ggsave("detection_grid.png", height = 4.5, width = 8.5, units = "in", dpi = 600)

# 11. STATISTICAL TESTING -------------------------------------------------------
test_df <- final_data %>%
  filter(predator == "Predator") %>%
  filter(sample_type == "rhizo") %>%
  filter(log_adj_abund > 0) %>%
  mutate(rep = as.factor(timepoint),
         treatment = as.factor(treatment),
         plant = as.factor(plant),
         timepoint = as.factor(timepoint))

mod <- lmer(log_adj_abund ~ treatment * plant + timepoint + (1 | rep), data = test_df)
car::Anova(mod)
summary(mod)

# treatment       147.8328  1  < 2.2e-16 ***
# plant             3.8371  1    0.05013 .
# timepoint         3.3009  3    0.34751
# treatment:plant  24.3926  1  7.857e-07 ***

specs <- emmeans(mod, pairwise ~ treatment * plant)
results_df <- as.data.frame(specs$contrasts) %>%
  mutate(padj = adjust_pvalue(p.value, method = "BH")) %>%
  mutate(signif = ifelse(padj<0.05, "*", ""))
results_df

# Panel C's model (mod2) was already fit in section 10 so its contrasts
# could drive the plotted significance brackets -- reprinted here for the
# record rather than refit.
car::Anova(mod2)
summary(mod2)

# treatment       14.2182  1  0.0001628 ***
# plant            0.9983  1  0.3177270
# timepoint        1.2320  3  0.7453341
# treatment:plant  2.1447  1  0.1430617

richness_contrasts

# 12. COMBINE PLOTS -------------------------------------------------------------
tag_shift <- theme(plot.tag.position = c(-0.06, 0.95))

legend_bc_plot <- ggplot(plot_df, aes(x = plant, y = log_adj_abund, fill = treatment, color = treatment)) +
  geom_point(size = 5, shape = 21) +
  scale_fill_manual(values = c("Conventional" = "#FED439FF", "Organic" = "#709AE1FF"), name = "Treatment",
                    guide = guide_legend(title.position = "top", title.hjust = 0.5)) +
  scale_color_manual(values = c("Conventional" = "#FED439FF", "Organic" = "#709AE1FF"), name = "Treatment") +
  theme_linedraw() +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10.5)
  ) 

#Extract legend
extract_legend_bc <- function(p) {
  g <- ggplot2::ggplotGrob(p)
  idx <- grep("^guide-box", g$layout$name)
  if (length(idx) == 0) return(grid::nullGrob())
  non_empty <- idx[!vapply(g$grobs[idx], function(x) inherits(x, "zeroGrob"), logical(1))]
  if (length(non_empty) == 0) return(grid::nullGrob())
  g$grobs[[non_empty[1]]]
}

legend_bc_grob <- extract_legend_bc(legend_bc_plot)

# B and C explicitly forced to equal widths (1:1) in one row.
bc_row <- (abund_plot + labs(tag = "B") | richness_plot + labs(tag = "C")) +
  plot_layout(widths = c(1, 1))

right_col <- (
  wrap_elements(full = legend_bc_grob) /
    bc_row /
    (detection_grid_plot + labs(tag = "D")) /
    plot_spacer()
) +
  plot_layout(heights = c(1, 1.5, 1, 1.5))

combined_figure <- (wrap_elements(full = bubble_plot) + labs(tag = "A") | plot_spacer() | right_col) +
  plot_layout(widths = c(3, 0.02, 1.7)) &
  theme(
    plot.tag = element_text(face = "bold", size = 15),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )

combined_figure
ggsave("bubble_abund_combined2.png", combined_figure, height = 12, width = 17, units = "in", dpi = 600)
