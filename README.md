# ACID Microbiome
### Gut enterotype- and body mass index (BMI)-dependent effects of anthocyanin supplementation on gut microbiota composition in individuals at risk for cognitive decline: a randomized placebo-controlled trial
## Scripts
* ***DADA2:*** Processing of raw sequences
* ***pre_processing.R:*** Performs data cleaning, preparation and rarefaction
* ***top10_taxa.R:*** Generates plots of the top 10 taxa at genus and phylum levels
* ***alpha_diversity.R:*** Calculates and statistically tests various alpha diversity metrics
* ***beta_diversity.R:*** Computes beta diversity measures
* ***enterotype_dmm.R:*** Classifies microbiome samples into distinct enterotypes using Dirichlet Multinomial Mixture (DMM) models
* ***BMI_age_differential_abundance.R:*** Performs a differential abundance analysis stratified by body mass index (BMI) and age (quartile)
* ***effectofintervention_microbiome.R:*** Assess the overall effect of intervention on microbiome composition, as well as stratified effects by enterotype, baseline BMI, and age quartiles
* ***effectofintervention_cognition.R:*** Assess the overall and stratified effects of intervention on cognition, incorporating microbiome and mediation analysis
