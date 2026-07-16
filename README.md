# micropredators
 Data files and scripts for an experiment involving soils collected from long-term management regimes (72-year corn monoculture or 17-year diverse organic rotation). These soils were collected and used in a greenhouse experiment involving 13CO2 pulse-labeling to track the incorporation of plant carbon within microbial communities in the rhizosphere of two different plant species (hairy vetch or winter wheat). 

 This repository contains supporting information for Schaedel et al. in prep: *Plant species and land management structure bacterial trophic dynamics in the rhizosphere*. DOI: unassigned

## Research objectives
- A
- B

## Methodology
- Microcosm experiments involved two carbon isotopic labeling treatments, two soil legacies, and two plant species:
  - 12CO2 / 13CO2
  - Conventional (conv) / organic (org)
  - Hairy vetch / winter wheat
- Isopycnic centrifugation and density fractionation was performed on genomic DNA extracted from rhizosphere soil samples collected at 4 timepoints (1/3/7/14 days) after labeling
- Amplicon sequencing of bacterial 16S rRNA was conducted, followed by multiple-window-high-resolution-stable-isotope-probing (MW-HR-SIP) to identify isotopically labeled taxa

## Availability of sequence data
16S rRNA sequence data are available in the NCBI Short Read Archives (BioProject Unassigned).

## Data files
- c_min_df.csv
  - Carbon mineralization data from GC-MS analayis, including respiration of 12C/13C xylose and cellulose 
- dna_yield.csv
  - DNA yield data for microcosm cummunities, measured with Picogreen 
- seq_l2fc.csv
  - HTSSIP results identifying taxa with significant shifts in buoyant density due to incorporation of isotope label
- incorp_bulk_rare_abund.RDS
  - Incorporator taxa present in rarefied microcosm communities
  - Includes associated metadata (taxonomy, rrn, maximum log2-foldchange, latency) for each incorporator
  - R data object

## Scripts
- Figures2_3.R
  - Carbon mineralization rates and cumulative mineralization
- Figure4.R
  - Relative abundance of bacterial communities in microcosms
  - Beta-diversity shifts in response to carbon addition
- Figure5.R
  - Incorporator venn diagram (detected vs labeled)
- Figure6.R
  - Carbon assimilation dynamics of top incorporator taxa 
- Figure7.R
  - Incorporator growth responses by tillage regime  
