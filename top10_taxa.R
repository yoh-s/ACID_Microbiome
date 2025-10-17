# ================================================
# Baseline Microbiome Composition Visualization
# ================================================
# Libraries
library(magrittr)
library(stringr)
library(phyloseq)
library(microViz)
library(microbiome)
library(tidyverse)

# ================================================
# Prepare Final Dataset
# ================================================
# Use ps_study8_10 for differential abundance analysis
ps_study8_10

# Remove low-depth samples (<12,358 reads)
ps_final <- subset_samples(
  ps_study8_10,
  !(Screening_code %in% c("P-076_24", "B-147_12", "P-330_24", "P-073_0", "P-136_0", "P-173_0", "P-073_24"))
)

# Remove zero-sum taxa
ps_final <- prune_taxa(taxa_sums(ps_final) > 0, ps_final)
summarize_phyloseq(ps_final)

# Replace sample metadata with updated one (e.g., rarefied metadata)
sample_data(ps_final) <- df_ps_study810_rare_f3

# Check taxonomic composition
length(unique(tax_table(ps_final)[, 2]))  # 12 phyla (2 Archaea + 10 Bacteria)

# Collapse at Genus level
tax_glom(ps_final, taxrank = "Genus")  # ~295 genera

# ================================================
# Filter Low-Prevalence Taxa and Transform
# ================================================
ps_study5_filter <- microViz::tax_filter(ps_final, min_prevalence = 0.05)  # 5% prevalence
psstudy5_filtercomp <- microbiome::transform(ps_study5_filter, transform = "compositional")

# Subset baseline samples
base_filtercomp <- subset_samples(psstudy5_filtercomp, Time == "Baseline")
base_filtercomp <- prune_taxa(taxa_sums(base_filtercomp) > 0, base_filtercomp)  # 1240 ASVs

# ================================================
# PHYLUM-LEVEL PLOTS
# ================================================
# Aggregate to Phylum
base_comp_phy <- tax_glom(base_filtercomp, taxrank = "Phylum")
tax_table(base_comp_phy) <- tax_table(base_comp_phy)[, 1:2]
base_comp_phy <- prune_taxa(taxa_sums(base_comp_phy) > 0, base_comp_phy)

# Subset top 6 phyla
base_top6phyla <- names(sort(taxa_sums(base_comp_phy), decreasing = TRUE)[1:6])
base_phyla6_subset <- prune_taxa(base_top6phyla, base_comp_phy)

# Prepare sorting info
phyla_to_sort_b <- data.frame(
  id = 1:6,
  phyla = as.character(tax_table(base_phyla6_subset)[, "Phylum"]),
  otu = as.character(taxa_names(base_phyla6_subset))
)
rownames(phyla_to_sort_b) <- phyla_to_sort_b$otu
base_phylum_ranks <- base_phyla6_subset %>% otu_table() %>% rowSums() %>% sort(TRUE) %>% names()
phyla_to_sort_b <- phyla_to_sort_b[base_phylum_ranks, ]

# Convert to proportions
base_prop_phylum6 <- transform_sample_counts(base_phyla6_subset, function(x) x / sum(x))

# Melt for plotting
base_bardat <- psmelt(base_prop_phylum6) %>%
  mutate(
    Sample = as.numeric(factor(Sample)),
    OTU = factor(OTU, levels = phyla_to_sort_b$otu, labels = phyla_to_sort_b$phyla)
  )

# Sort by Firmicutes abundance
base_firmicutes_order <- base_bardat %>%
  filter(OTU == "Firmicutes") %>%
  arrange(Abundance) %>%
  select(Sample)

base_bardat %<>%
  mutate(Sample = as.numeric(factor(Sample, levels = factor(base_firmicutes_order$Sample))))

# Function to remove % symbol in Y-axis labels
remove_percent_sign <- function(x) as.character(as.numeric(x) * 100)

# Define colors for major phyla
area_phylum <- c(
  "Actinobacteriota" = "#ee3333",
  "Bacteroidota" = "#ee7722",
  "Euryarchaeota" = "#ffee33",
  "Firmicutes" = "deepskyblue2",
  "Proteobacteria" = "#3366aa",
  "Verrucomicrobiota" = "#992288"
)

# Ensure treatment factor order
base_bardat$Treatment_Type <- factor(base_bardat$Treatment_Type, levels = c("Intervention", "Placebo"))

# ================================================
# GENUS-LEVEL PLOTS
# ================================================
base_comp_gen <- tax_glom(base_filtercomp, taxrank = "Genus")
tax_table(base_comp_gen) <- tax_table(base_comp_gen)[, 1:6]
base_comp_gen <- prune_taxa(taxa_sums(base_comp_gen) > 0, base_comp_gen)

# Get top 10 genera
top10genus_base <- names(sort(taxa_sums(base_comp_gen), decreasing = TRUE)[1:10])
base_genus10 <- prune_taxa(top10genus_base, base_comp_gen)

# Sort genera by abundance
genus_to_sort <- data.frame(
  id = 1:10,
  phyla = as.character(tax_table(base_genus10)[, "Genus"]),
  otu = as.character(taxa_names(base_genus10))
)
rownames(genus_to_sort) <- genus_to_sort$otu
genus_ranks <- base_genus10 %>% otu_table() %>% rowSums() %>% sort(TRUE) %>% names()
genus_to_sort <- genus_to_sort[genus_ranks, ]

# Proportion transform and melt
base_genusprop <- transform_sample_counts(base_genus10, function(x) x / sum(x))
base_genus_10 <- psmelt(base_genusprop) %>%
  mutate(
    Sample = as.numeric(factor(Sample)),
    OTU = factor(OTU, levels = genus_to_sort$otu, labels = genus_to_sort$phyla)
  )

# Sort by Bacteroides abundance
bacteroides_order <- base_genus_10 %>%
  filter(OTU == "Bacteroides") %>%
  arrange(Abundance) %>%
  select(Sample)

base_genus_10 %<>%
  mutate(Sample = as.numeric(factor(Sample, levels = factor(bacteroides_order$Sample))),
         Treatment_Type = factor(Treatment_Type, levels = c("Intervention", "Placebo")))

# Genus color palette
area_genus <- c(
  "Agathobacter" = "#bcbd22", "Akkermansia" = "darkorchid",
  "Clostridia_UCG-014" = "deeppink1", "Bacteroides" = "#66aa55",
  "Christensenellaceae_R-7_group" = "darkslategray", "Coprococcus" = "gray",
  "Faecalibacterium" = "#d62728", "Ruminococcus" = "#17becf",
  "Subdoligranulum" = "#0000ff", "UCG-002" = "saddlebrown"
)

# Clean legend labels (with line breaks)
genus_labels <- c(
  "Agathobacter", "Akkermansia", "Bacteroides",
  "Christensenellaceae\nR-7_group", "Clostridia_UCG-014",
  "Coprococcus", "Faecalibacterium", "Ruminococcus",
  "Subdoligranulum", "UCG-002"
)
