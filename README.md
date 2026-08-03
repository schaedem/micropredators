# micropredators
 Data files and scripts for an experiment involving soils collected from long-term management regimes (72-year corn monoculture or 17-year diverse organic rotation). These soils were collected and used in a greenhouse experiment involving ¹³CO₂ pulse-labeling to track the incorporation of plant carbon within microbial communities in the rhizosphere of two different plant species (hairy vetch or winter wheat). 

 This repository contains supporting information for Schaedel et al. in prep: *Plant species and land management structure bacterial trophic dynamics in the rhizosphere*. DOI: unassigned

## Research objectives
- Apply a pulse-chase methodology to determine whether the trophic dynamics of putative bacterial micropredators in the rhizosphere are structured by plant host identity and soil management legacy.
- Test whether micropredator enrichment in the rhizosphere is a general feature of plant-microbe interactions, by comparing this study's findings against a meta-analysis of 16S rRNA amplicon data from 31 previously published studies.

## Methodology
- Greenhouse experiment involved two carbon isotopic labeling treatments, two soil legacies, and two plant species:
  - ¹²CO₂ / ¹³CO₂
  - Conventional (conv) / organic (org)
  - Hairy vetch / winter wheat
- Isopycnic centrifugation and density fractionation was performed on genomic DNA extracted from rhizosphere soil samples collected at 4 timepoints (1/3/7/14 days) after labeling
- Amplicon sequencing of bacterial 16S rRNA was conducted, followed by multiple-window-high-resolution-stable-isotope-probing (MW-HR-SIP) to identify isotopically labeled taxa

## Availability of sequence data
16S rRNA sequence data are available in the NCBI Short Read Archives (BioProject PRJNA1496949).

## Data files
- **Feature tables & taxonomy** (this study and each of the 31 meta-analysis studies): ASV feature table (`feature-table.tsv`/`.biom`) and cleaned taxonomy (`taxonomy-fixed.tsv` / `taxonomy_clean.tsv`), output of `process_study_sequences.py` / `process_meta_analysis_sequences.py`.
- **Sample metadata**: `conv_sample_metadata.csv`, `org_sample_metadata_new.csv` (treatment, plant, timepoint, replicate, isotope-labeling status), and per-study metadata for the meta-analysis (soil type, primer set, sample counts).
- **IRMS data**: `irms_results_all_final.csv` (percent C, percent N, δ¹³C enrichment per sample).
- **Phyloseq objects**: `nifa_joint.phyloseq.RDS` (combined ASV table, taxonomy, and sample metadata for this study).
- **Incorporator ASVs**:
  - `incorp_soil_adj_abun.RDS`: rhizosphere, sample-by-ASV table of isotopically labeled ASVs with rrn/qPCR-adjusted abundance (`adj_abund`, `log_adj_abund`) and sample metadata (treatment, plant, timepoint, replicate). Used in Figures 1 and 3.
  - `joint_incorps_all.RDS`: labeled-ASV table across all samples with taxonomy (Phylum/Order/Genus) and results from MW-HR-SIP. Used in Figure 1.
  - `source_incorp_by_day.RDS`,`conv_source_incorp_by_day.RDS` : per-ASV record of which day(s) each ASV was detected as labeled, split by soil background (organic / conventional); used to build the presence/absence matrix for hierarchical clustering and to restrict the ASV table to labeled ASV-sample combinations. Used in Figure 2.

## Scripts
- **Figure 1** — `fig1_polr_first_label.R`: Cumulative link (proportional-odds ordinal logistic regression) model of the day of first ¹³C-label detection as a function of treatment, plant host, and predator taxonomic group. Produces Figure 1, Supplemental Figure 1, and Supplemental Tables 2–3.

- **Figure 2** — `fig2_composition_abundance_richness.R`: Bubble plot of the most abundant labeled taxa (hierarchically clustered within predator/non-predator groups) alongside labeled predator abundance and richness in the rhizosphere. Produces Figure 2.

- **Figure 3** — `fig3_predator_prey_ratios.R`: Linear mixed-effects models of the labeled predator:prey abundance ratio by soil management, plant host, and time, including trend tests and relationships with soil carbon/nitrogen/16S abundance. Produces Figure 3, Supplemental Table 4, and Supplemental Figure 3.

- **Figure 4, step 1** — `fig4_per_timepoint_deseq.R`: DESeq2 differential abundance (rhizosphere vs. bulk soil) for this study, run separately for each of the 8 soil-management x plant x timepoint combinations, with ASVs aggregated to the Order level.

- **Figure 4, step 2** — `fig4_rhizosphere_effect.R`: Aggregates the per-timepoint DESeq2 results into taxonomic trophic groups via inverse-variance weighting and fits a hierarchical mixed-effects meta-regression (`metafor::rma.mv`) to estimate the rhizosphere effect for each group, separately by plant. Produces Figure 4.

- **Figure 5, step 1** — `fig5_per_study_deseq.R`: DESeq2 differential abundance template applied once per meta-analysis study (rhizosphere vs. bulk soil, Order-level aggregation, ashr LFC shrinkage); shown for one representative study.

- **Figure 5, step 2** — `fig5_meta_analysis.R`: Compiles per-study DESeq2 results across all 31 studies, aggregates to trophic groups, and fits the global hierarchical meta-regression. Produces Figure 5.

- **Sequence processing (this study)** — `process_study_sequences.py`: QIIME2/DADA2 pipeline importing, primer-trimming, denoising, taxonomically classifying, and building a phylogenetic tree from this study's raw paired-end 16S amplicon reads (two sequencing runs, merged).

- **Sequence processing (meta-analysis)** — `process_meta_analysis_sequences.py`: QIIME2/DADA2 pipeline applied per study across the 31 meta-analysis studies, handling per-study primer sets and single- vs. paired-end reads, with taxonomy assigned via a vsearch consensus classifier against full-length SILVA.
