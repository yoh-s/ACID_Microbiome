#load library
library(decontam)
library(DT)
library(qiime2R)
library(phyloseq)
library(microbiome)
library(tidyverse)
library(ape)
library(microbiomeutilities)
library(Biostrings)
#library(phangorn)
#library(DECIPHER)
library(microViz)
library(readr)
library(dplyr)
library(dada2)
library(biomformat)
#library(speedyseq) ###orient taxa
library(readr)
library(phyloseq)
library(microViz)

##Decontamination

#create a phyloseq object
pscontam_july <- qza_to_phyloseq(
  features="unfiltered_qiime/maxeeon_pseudotable_contam.qza",
  taxonomy="unfiltered_qiime/silva341_785_finaltaxonomy.qza",
  metadata = "unfiltered_qiime/check_blank_metadata.txt"
)

pscontam_july   ###implying a phyloseq object that contains negative control
dir.create("decontam")
saveRDS(pscontam_july, "decontam/pscontam_july.rds")

##one zero sum taxa
pscontam_july <- subset_samples(pscontam_july, Description != "MockCocktail")

pscontam_julya <- prune_taxa(taxa_sums(pscontam_july) > 1, pscontam_july)
pscontam_julya <- prune_samples(sample_sums(pscontam_julya) > 0, pscontam_julya)
pscontam_julya  #contains replicate one and two #mockcocktail

summarize_phyloseq(pscontam_julya)

#check the library size of samples and negative controls
df <- as.data.frame(sample_data(pscontam_julya))
df$LibrarySize <- sample_sums(pscontam_julya)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Description)) + geom_point() +
  labs(x = "Sample index", y = "Sequencing depth") +
  theme_bw() +
  theme(legend.title = element_blank()) 
ggsave("decontam/library_size_pscontam_julya.jpeg", dpi = 600)

trial <- isContaminant(pscontam_julya, method = "prevalence", neg = "is.neg")
hist(trial$p, 100, xlab = "decontam score", ylab = "Number of ASVs", main = "Prevalence method")

#prevalence threshold of 0.1
contamdf.prev.05 <- isContaminant(pscontam_julya, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev.05$contaminant)
##FALSE  TRUE 
##6060    34

head(which(contamdf.prev.05$contaminant))
##[1]  143  519 1182 1448 1584 1703

row_indices05 <- which(contamdf.prev.05$contaminant) #grab the row indices that correspond with identified contaminants to locate taxonomic information in the corresponding OTU file
otu_taxonomy <- phyloseq::tax_table(pscontam_julya)
classification <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

taxonomy_table05 <- tibble()
for (i in row_indices05){
  loc <-  contamdf.prev.05[i, 0]
  tax_key <- row.names(loc)
  tax_value <- otu_taxonomy[tax_key, ]
  taxonomy_table05 <- rbind(taxonomy_table05, tax_value)
}
names(taxonomy_table05) <- classification
datatable(taxonomy_table05)   ####9 bacteria data  table

write.csv(taxonomy_table05, "decontam/taxonomy_table05july.csv")

tax_pscontamf_july <- data.frame(phyloseq::tax_table(pscontam_julya))
write.csv(tax_pscontamf_july, "decontam/tax_pscontamf_july.csv")

#### p = 0.5 prevalence threshold will be used
datatable(taxonomy_table05)   ####9 list of contaminant taxa

contaminant_taxa = rownames(taxonomy_table05)

##use decontam phyloseq object
pscontam_julya

##Extract taxa names from the phyloseq object
pscontam_julya_names <- phyloseq::taxa_names(pscontam_julya)

##Extract taxa names not included in the contaminant taxa
pscontam_julya_names_1 <- pscontam_julya_names[!(pscontam_julya_names %in% contaminant_taxa)]

##create a new phyloseq object without contaminant
pscontam_julya1 = prune_taxa(pscontam_julya_names_1, pscontam_julya)

##Export taxa_table as a dataframe
table_pscontam_julya1 <- data.frame(phyloseq::tax_table(pscontam_julya1))
write.csv(table_pscontam_julya1, "decontam//table_pscontam_julya1.csv")

##Remove Eukaryota, Archaea, and Cyanobacteria
pscontam_julya2 = subset_taxa(pscontam_julya1, !(Kingdom %in% c("d__Eukaryota")) &
                             !(Phylum == "Cyanobacteria"))

ps_study = subset_samples(pscontam_julya2, Description == "Study")

sum(taxa_sums(ps_study) == 0)  ##23 zero otu

#Remove any OTUs that are 0
ps_study <- prune_taxa(taxa_sums(ps_study) > 0, ps_study)

##Read new metadata from final blueberry project
trialfinal_metadata070323 <- read_delim("unfiltered_qiime/trialfinal_metadata070323.tsv", 
                                        delim = "\t", escape_double = FALSE, 
                                        trim_ws = TRUE)
View(trialfinal_metadata070323)

##add the new metadata to the phyloseq object
study_replicate <- sample_data(trialfinal_metadata070323)

#new phyloseq object
ps_study1 <- phyloseq(sample_data(study_replicate), otu_table(ps_study), tax_table(ps_study))
ps_samplingdepth <- sample_sums(ps_study1)
study_depth <- ps_samplingdepth %>% as.data.frame()
write.csv(study_depth, "decontam/study_depth.csv")
##Check sequence depth

##Select samples with high sequencing depth out of the duplicates
##Read new metadata from final blueberry project
trialfinal_metadata070323 <- read_delim("unfiltered_qiime/trialfinal_metadata070323.tsv", 
                                        delim = "\t", escape_double = FALSE, 
                                        trim_ws = TRUE)
View(trialfinal_metadata070323)

##add the new metadata to the phyloseq object
study_replicate <- sample_data(trialfinal_metadata070323)
sample_names(study_replicate) <- study_replicate$`sample-id`

#new phyloseq object
ps_study1 <- phyloseq(sample_data(study_replicate), otu_table(ps_study), tax_table(ps_study))
sample_sums(ps_study1)
View(sample_data(ps_study1))
sample_variables(ps_study1)
############################
##Single samples without duplicate
ps_study2 <- subset_samples(ps_study1, Highsequencingdepth == "1") ##6014

#Remove any OTUs that are 0
ps_study2 <- prune_taxa(taxa_sums(ps_study2) > 0, ps_study2)  ##5528
sample_sums(ps_study2) %>% range() ##7046 157418
##7046, 7396, 9804, below 10,000 reads

summarize_phyloseq(ps_study2)  ##Sparsity = 0.9413

variability <- plot_taxa_cv(ps_study2, plot.type = "scatter")
variability + scale_x_log10()

# Create table, number of features for each phyla
table(tax_table(ps_study2)[, "Phylum"], exclude = NULL)

# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps_study2),
               MARGIN = ifelse(taxa_are_rows(ps_study2), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps_study2),
                    tax_table(ps_study2))

phyla_distribution <- plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
dir.create("output_results")
write.csv(phyla_distribution, "output_results/phyla_distribution.csv")

# Filter entries with unidentified Phylum.
ps_study3 = subset_taxa(ps_study2, !Phylum %in% c("Elusimicrobiota", "Spirochaetota"))
sum(taxa_sums(ps_study3) == 0) ##486

#Remove any doubleton
ps_study3 <- prune_taxa(taxa_sums(ps_study3) > 2, ps_study3)

#Remove unclassified and uncultured taxa at the Class level
ps_study3 <- subset_taxa(ps_study3, !is.na(Class) & Class != "uncultured")

tax_psstudy3 <- data.frame(tax_table(ps_study3))
write.csv(tax_psstudy3, "output_results/tax_psstudy3.csv")
##ASV13_5360, ASV13_3337, and ASV13_4871 had uncultured at the class level

##Clean taxonomy file
edit_tax_psstudy3 <- read.csv("output_results/tax_psstudy3.csv")

#change rownames to column
edit_tax_psstudy3 <- edit_tax_psstudy3 %>% column_to_rownames("X")
tax_table(ps_study3) <- as.matrix(edit_tax_psstudy3)
View(tax_table(ps_study3))

##read eukaryotic phyloseq
psbermeta_euk1 <- readRDS("output_results/psbermeta_euk1.rds")

##merge two phyloseq objects
ps_study4 <- merge_phyloseq(ps_study3, psbermeta_euk1)
View(tax_table(ps_study4))

##Save
saveRDS(ps_study4, "output_results/ps_study4.rds")

sampling_depth <- sample_sums(ps_study4) %>% as.data.frame() %>% dplyr::rename(depth = 1) %>% 
  rownames_to_column("sample_id")

##change sample names of the phyloseq object
ps_study4a <- ps_study4

##Assign the sample names to screening code
sample_names(ps_study4a) <- sample_data(ps_study4a)$Screening_code

##Add the final metadata file to the phyloseq object
#Read the metadata
final_metadata24 <- read_csv("final_metadata2024may.csv")
final_metadata24 <- sample_data(final_metadata24)
sample_names(final_metadata24) <- final_metadata24$...1

#new phyloseq object
ps_study5 <- phyloseq(sample_data(final_metadata24), otu_table(ps_study4a), tax_table(ps_study4a))

#Remove any doubleton
ps_study5 <- prune_taxa(taxa_sums(ps_study5) > 2, ps_study5)

#remove eukaryotes
ps_study5 <- subset_taxa(ps_study5, Kingdom != "Eukaryota")

sum((taxa_sums(ps_study5) < 10)) ##680

##read the rep.seqs
repseqs <- Biostrings::readDNAStringSet("unfiltered_qiime/all_dnasequences.fasta", format = "fasta")
saveRDS(repseqs, "output_results/repseqs.rds")

##all included
ps_study6 <- phyloseq(otu_table(ps_study5), tax_table(ps_study5), sample_data(ps_study5),
                      refseq(repseqs))
saveRDS(ps_study6, "output_results/ps_study6.rds")

sample_sums(ps_study6) %>% range()
##build a phylogenetic tree

study6_tree <- read.tree("output_results/tree2024/mafft_fasttree/tree_rooted/tree.nwk")

ps_study7 <- phyloseq(otu_table(ps_study6),sample_data(ps_study6),tax_table(ps_study6),
                      refseq(ps_study6),phy_tree(study6_tree))
View(tax_table(ps_study7))

#Remove any doubleton
ps_study7 <- prune_taxa(taxa_sums(ps_study7) > 2, ps_study7)

plot_tree(ps_study7)

##remove eukaryotes
ps_study8 <- subset_taxa(ps_study7, Kingdom != "Eukaryota")
saveRDS(ps_study8, "output_results/ps_study8.rds")

##microViz::tax_filter(ps_study8, min_prevalence = 0.001)
##the final phyloseq object passed 0.01 and 0.1% filter

##12358 rarefy even depth
ps_study8_rare <- rarefy_even_depth(ps_study8, sample.size = 12358, rngseed = 1988, replace = F)
##5237
sum(taxa_sums(ps_study8_rare) > 10)  ##3934

##number of ASV below 10 reads in total samples are about 1303

##Taxa sum above 10 is 4570, lower than 10 reads were 749
ps_study8_10 <- prune_taxa(taxa_sums(ps_study8) > 10, ps_study8)
5319-4570 ##749

sum(taxa_sums(ps_study8_10) > 10)  ##4570

##Recode the treatment groups to Anthocyanins and Placebo groups
sample_data(ps_study8_10) <- sample_data(ps_study8_10) %>% as.data.frame() %>% 
  mutate(Treatment_group = fct_recode(Treatment_Type,
                             "Treatment A" = "Intervention",
                             "Treatment B" = "Placebo")) 

df_ps_study8_10 <- data.frame(sample_data(ps_study8_10))

df_ps_study8_10mod <- df_ps_study8_10 %>% mutate(Treatment_Type = fct_recode(Treatment_Type, "Intervention" = "Treatment A",
                                     "Placebo" = "Treatment B"))

sample_data(ps_study8_10) <- df_ps_study8_10mod

ps_study810_rare <- rarefy_even_depth(ps_study8_10, sample.size = 12358, rngseed = 1988, replace = F) ##4564
sum(taxa_sums(ps_study810_rare) > 2)  ##4535

##number of ASV below 10 reads in total samples are about 620

##merge at the genus level
ps_study8genus<- tax_glom(ps_study8, taxrank = "Genus")

####QIIME2 EXPORT
####phylogenetic tree
treegenus = phy_tree(ps_study8genus)
ape::write.tree(treegenus, "treegenus.tree")

##feature table export to qiime2
otu_ps_studygenus <- as(otu_table(ps_study8genus),"matrix")
otu_genusbiom<-make_biom(data=otu_ps_studygenus)
write_biom(otu_genusbiom,"otu_genusbiom.biom")


##Export for picrust2
#use ps_study8

dir.create("picrust2_october24")
#feature table export to qiime2
otu_ps_study8 <- as(otu_table(ps_study8),"matrix")
otu_psstudy8_biom<-make_biom(data=otu_ps_study8)
write_biom(otu_psstudy8_biom,"picrust2_october24/otu_psstudy8.biom")

##read the rep.seqs
##Export to generate phylogenetic tree using qiime2
ps_study8 %>%
  refseq() %>%
  Biostrings::writeXStringSet("picrust2_october24/psstudy8_seq.fna", append=FALSE,
                              compress=FALSE, compression_level=NA, format="fasta")

##Subset baseline
base_psstudy8 <- subset_samples(ps_study8, Time == "Baseline")
df_base_psstudy8 <- data.frame(sample_data(base_psstudy8))

##set levels for age group
df_base_psstudy8$Age_quartile <- factor(df_base_psstudy8$Age_quartile, levels = c("60-64_years", "65-68_years", "68-73_years", "73-80_years"))

##set levels for age group
df_base_psstudy8$Treatment_Type <- factor(df_base_psstudy8$Treatment_Type, levels = c("Treatment A", "Treatment B"))

chisq_test(df_base_psstudy8$Age_quartile, df_base_psstudy8$Treatment_Type)

##create contigency table
table(df_base_psstudy8$Age_quartile, df_base_psstudy8$Treatment_Type)

#plot age group
df_base_psstudy8 %>% 
  ggplot(aes(x = Treatment_Type, fill = Age_quartile)) +
  geom_bar(position = "fill")

# load packages
library(ggstatsplot)
library(ggplot2)
# plot
age_treatment <- ggbarstats(
  data = df_base_psstudy8,
  x = Treatment_Type,
  y = Age_quartile,
  sample.size.label.args = list(size =3)) +
  labs(caption = NULL, x = "") +
  scale_x_discrete(labels = c("60 - 64 years", "65 - 68 years", "68 - 73 years", "73 - 80 years")) +
  scale_fill_discrete(labels = c("Placebo", "Intervention")) +
  guides(fill = guide_legend(title = "")) +
  theme(axis.text.x = element_text(size = 11, colour = "black"),
        legend.text = element_text(size = 11))

ggsave(plot = age_treatment,"age_treatment.tiff",dpi = 600,width = 6,height = 5)

age_treatment <- ggbarstats(
  data = df_base_psstudy8,
  x = Age_quartile,
  y = Treatment_Type,
  sample.size.label.args = list(size =3)) +
  labs(caption = NULL, x = "") +
  scale_x_discrete(labels = c("Intervention", "Placebo")) +
  guides(fill = guide_legend(title = "")) +
  theme(axis.text.x = element_text(size = 11, colour = "black"),
        legend.text = element_text(size = 11)) 
dir.create("age_category")  
ggsave(plot = age_treatment,"age_category/age_treatment_final.tiff",dpi = 600,width = 6,height = 5)

library(ggstatsplot)
library(ggplot2)

ggbarstats(
  data = df_base_psstudy8,
  x = Age_quartile,
  y = Treatment_Type,
  sample.size.label.args = list(size = 3)
) +
  labs(caption = NULL, x = "") +
  scale_x_discrete(labels = c("60 - 64 years", "65 - 68 years", "68 - 73 years", "73 - 80 years")) +
  scale_fill_discrete(labels = c("Placebo", "Intervention")) +
  guides(fill = guide_legend(title = "")) +
  theme(axis.text.x = element_text(size = 11, colour = "black"),
        legend.text = element_text(size = 11)) +
  theme_void()


##set levels for age group
df_base_psstudy8$BMI_3_class <- factor(df_base_psstudy8$BMI_3_class, levels = c("Normal", "Overweight", "Obesity"))


##set levels for age group
df_base_psstudy8$Screening_Diagnosis <- factor(df_base_psstudy8$Screening_Diagnosis, levels = c("CMD", "MCI"))

##CMD and BMI
bmi_screening <- ggbarstats(
  data = df_base_psstudy8,
  x = Screening_Diagnosis,
  y = BMI_3_class,
  sample.size.label.args = list(size =3)
) +
  labs(caption = NULL, x = "") +
  scale_fill_discrete(labels = c("Mild cognitive impairment", "Cardiometabolic diseases")) +
  guides(fill = guide_legend(title = "")) +
  theme(axis.text.x = element_text(size = 11, colour = "black"),
        legend.text = element_text(size = 11))

ggsave(plot = bmi_screening,"bmi_screening.tiff",dpi = 600,width = 6,height = 5)

