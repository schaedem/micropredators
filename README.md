# micropredators
 Data files and scripts for an experiment involving soils collected from long-term management regimes (72-year corn monoculture or 17-year diverse organic rotation). These soils were collected and used in a greenhouse experiment involving 13CO2 pulse-labeling to track the incorporation of plant carbon within microbial communities in the rhizosphere of two different plant species (hairy vetch or winter wheat). 

 This repository contains supporting information for Schaedel et al. 2025: *Long term tillage regime structures bacterial carbon assmilation*. DOI: 10.1101/2025.01.07.631211

## Research objectives
- Determine the impact of a 42 year legacy of no-till or moldboard plow tillage on bacterial carbon assimilation
- Investigate the role of divergent life histories resulting from management legacy in contributing to carbon dynamics

## Methodology
- Microcosm experiments involved two carbon isotopic labeling treatments and two soil legacies:
  - 12C-xylose + 13C-cellulose / 13C-xylose + 12C-cellulose
  - No-till biomass harvested (NTH) / plow-till biomass harvested (PTH)
- Isopycnic centrifugation and density fractionation was performed on genomic DNA extracted from microcosm soil samples
- Amplicon sequencing of bacterial 16S rRNA was conducted, followed by multiple-window-high-resolution-stable-isotope-probing (MW-HR-SIP) to identify isotopically labeled taxa

## Availability of sequence data
16S rRNA sequence data are available in the NCBI Short Read Archives (BioProject PRJNA1170979).

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
