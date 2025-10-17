# ===============================
# Load Required Libraries
# ===============================
library(decontam)
library(DT)
library(qiime2R)
library(phyloseq)
library(microbiome)
library(tidyverse)
library(ape)
library(microbiomeutilities)
library(Biostrings)
library(microViz)
library(readr)
library(dada2)
library(biomformat)
library(dplyr)

# ===============================
# Create Phyloseq Object
# ===============================
pscontam_july <- qza_to_phyloseq(
  features = "unfiltered_qiime/maxeeon_pseudotable_contam.qza",
  taxonomy = "unfiltered_qiime/silva341_785_finaltaxonomy.qza",
  metadata = "unfiltered_qiime/check_blank_metadata.txt"
)

# Save initial object
dir.create("decontam", showWarnings = FALSE)
saveRDS(pscontam_july, "decontam/pscontam_july.rds")

# ===============================
# Preprocessing
# ===============================
# Remove mock samples
pscontam_july <- subset_samples(pscontam_july, Description != "MockCocktail")

# Remove zero-sum taxa and samples
pscontam_julya <- prune_taxa(taxa_sums(pscontam_july) > 1, pscontam_july)
pscontam_julya <- prune_samples(sample_sums(pscontam_julya) > 0, pscontam_julya)

# Summarize phyloseq object
summarize_phyloseq(pscontam_julya)

# ===============================
# Library Size Check
# ===============================
df <- as.data.frame(sample_data(pscontam_julya))
df$LibrarySize <- sample_sums(pscontam_julya)
df <- df[order(df$LibrarySize), ]
df$Index <- seq(nrow(df))

ggplot(df, aes(x = Index, y = LibrarySize, color = Description)) +
  geom_point() +
  labs(x = "Sample index", y = "Sequencing depth") +
  theme_bw() +
  theme(legend.title = element_blank())

ggsave("decontam/library_size_pscontam_julya.jpeg", dpi = 600)

# ===============================
# Identify Contaminants
# ===============================
trial <- isContaminant(pscontam_julya, method = "prevalence", neg = "is.neg")
hist(trial$p, 100, xlab = "Decontam score", ylab = "Number of ASVs", main = "Prevalence method")

# Prevalence threshold = 0.5
contamdf.prev.05 <- isContaminant(pscontam_julya, method = "prevalence", neg = "is.neg", threshold = 0.5)
table(contamdf.prev.05$contaminant)

row_indices05 <- which(contamdf.prev.05$contaminant)

# Extract taxonomy of contaminants
otu_taxonomy <- phyloseq::tax_table(pscontam_julya)
classification <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
taxonomy_table05 <- tibble()

for (i in row_indices05) {
  tax_value <- otu_taxonomy[i, ]
  taxonomy_table05 <- rbind(taxonomy_table05, tax_value)
}

names(taxonomy_table05) <- classification
datatable(taxonomy_table05)

write.csv(taxonomy_table05, "decontam/taxonomy_table05july.csv")
write.csv(data.frame(phyloseq::tax_table(pscontam_julya)), "decontam/tax_pscontamf_july.csv")

# ===============================
# Remove Contaminants
# ===============================
contaminant_taxa <- rownames(taxonomy_table05)
pscontam_julya_names <- phyloseq::taxa_names(pscontam_julya)
pscontam_julya_names_1 <- pscontam_julya_names[!(pscontam_julya_names %in% contaminant_taxa)]
pscontam_julya1 <- prune_taxa(pscontam_julya_names_1, pscontam_julya)

# Remove Eukaryota and Cyanobacteria
pscontam_julya2 <- subset_taxa(pscontam_julya1, !(Kingdom %in% c("d__Eukaryota")) &
                                 !(Phylum == "Cyanobacteria"))

# Subset study samples only
ps_study <- subset_samples(pscontam_julya2, Description == "Study")
ps_study <- prune_taxa(taxa_sums(ps_study) > 0, ps_study)

# ===============================
# Update Metadata
# ===============================
trialfinal_metadata <- read_delim("unfiltered_qiime/trialfinal_metadata070323.tsv", delim = "\t", trim_ws = TRUE)
study_replicate <- sample_data(trialfinal_metadata)
sample_names(study_replicate) <- study_replicate$`sample-id`

ps_study1 <- phyloseq(sample_data(study_replicate), otu_table(ps_study), tax_table(ps_study))

# Subset high sequencing depth samples
ps_study2 <- subset_samples(ps_study1, Highsequencingdepth == "1")
ps_study2 <- prune_taxa(taxa_sums(ps_study2) > 0, ps_study2)

# ===============================
# Taxonomic Cleaning
# ===============================
ps_study3 <- subset_taxa(ps_study2, !Phylum %in% c("Elusimicrobiota", "Spirochaetota"))
ps_study3 <- prune_taxa(taxa_sums(ps_study3) > 2, ps_study3)
ps_study3 <- subset_taxa(ps_study3, !is.na(Class) & Class != "uncultured")

# Save cleaned taxonomy
write.csv(data.frame(tax_table(ps_study3)), "output_results/tax_psstudy3.csv")

# ===============================
# Merge Eukaryotic Phyloseq
# ===============================
psbermeta_euk1 <- readRDS("output_results/psbermeta_euk1.rds")
ps_study4 <- merge_phyloseq(ps_study3, psbermeta_euk1)
saveRDS(ps_study4, "output_results/ps_study4.rds")

# ===============================
# Update Sample Names and Final Metadata
# ===============================
sample_names(ps_study4) <- sample_data(ps_study4)$Screening_code
final_metadata24 <- read_csv("final_metadata2024may.csv")
final_metadata24 <- sample_data(final_metadata24)
sample_names(final_metadata24) <- final_metadata24$...1

ps_study5 <- phyloseq(sample_data(final_metadata24), otu_table(ps_study4), tax_table(ps_study4))
ps_study5 <- prune_taxa(taxa_sums(ps_study5) > 2, ps_study5)
ps_study5 <- subset_taxa(ps_study5, Kingdom != "Eukaryota")

# ===============================
# Add Reference Sequences
# ===============================
repseqs <- readDNAStringSet("unfiltered_qiime/all_dnasequences.fasta", format = "fasta")
ps_study6 <- phyloseq(otu_table(ps_study5), tax_table(ps_study5), sample_data(ps_study5), refseq(repseqs))
saveRDS(ps_study6, "output_results/ps_study6.rds")

# ===============================
# Build Phylogenetic Tree
# ===============================
study6_tree <- read.tree("output_results/tree2024/mafft_fasttree/tree_rooted/tree.nwk")
ps_study7 <- phyloseq(otu_table(ps_study6), sample_data(ps_study6), tax_table(ps_study6),
                      refseq(ps_study6), phy_tree(study6_tree))
ps_study7 <- prune_taxa(taxa_sums(ps_study7) > 2, ps_study7)
ps_study8 <- subset_taxa(ps_study7, Kingdom != "Eukaryota")
saveRDS(ps_study8, "output_results/ps_study8.rds")

# ===============================
# Rarefaction and Genus-level Aggregation
# ===============================
ps_study8_rare <- rarefy_even_depth(ps_study8, sample.size = 12358, rngseed = 1988, replace = FALSE)
ps_study8_10 <- prune_taxa(taxa_sums(ps_study8) > 10, ps_study8)

# Recode treatment groups
sample_data(ps_study8_10) <- sample_data(ps_study8_10) %>% as.data.frame() %>%
  mutate(Treatment_group = fct_recode(Treatment_Type,
                                      "Treatment A" = "Intervention",
                                      "Treatment B" = "Placebo"))

# Merge at genus level and export tree
ps_study8genus <- tax_glom(ps_study8, taxrank = "Genus")
treegenus <- phy_tree(ps_study8genus)
ape::write.tree(treegenus, "treegenus.tree")
