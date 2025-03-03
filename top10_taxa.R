library(magrittr)
library(stringr)
ps_study8_10 ###use for differential abundance analysis

##Remove sample below the 12,358 reads
ps_final <- subset_samples(ps_study8_10, !(Screening_code %in% c("P-076_24","B-147_12","P-330_24","P-073_0","P-136_0","P-173_0","P-073_24")))                            

##Remove zero sum
ps_final <- prune_taxa(taxa_sums(ps_final) > 0, ps_final) ##4569 ASVs
summarize_phyloseq(ps_final)

##this should be the one used for picruts2
dir.create("picrust2_october24")
#feature table export to qiime2
otu_ps_final <- as(otu_table(ps_final),"matrix")
otu_psfinal_biom<-make_biom(data=otu_ps_final)
write_biom(otu_psfinal_biom,"picrust2_october24/otu_psfinal.biom")

##read the rep.seqs
##Export to generate phylogenetic tree using qiime2
ps_final %>%
  refseq() %>%
  Biostrings::writeXStringSet("picrust2_october24/psfinal_seq.fna", append=FALSE,
                              compress=FALSE, compression_level=NA, format="fasta")





##change the metadata to the new one that calculated the metadata
sample_data(ps_final) <- df_ps_study810_rare_f3

length(unique(tax_table(ps_final)[,2])) ##12 phylum
##Two Archaea and ten bacteria

tax_glom(ps_final, taxrank = "Genus")  ##295

##filter for 5% prevalence using microViz
ps_study5_filter <- microViz::tax_filter(ps_final, min_prevalence = 0.05) ##1240

##Calculate relative abundance
psstudy5_filtercomp <- microbiome::transform(ps_study5_filter, transform = "compositional")

##Subset baseline 
base_filtercomp <- subset_samples(psstudy5_filtercomp, Time == "Baseline")

##Remove zero sum
base_filtercomp <- prune_taxa(taxa_sums(base_filtercomp) > 0, base_filtercomp) ##1240 ASVs

##Plot at the phylum level
base_comp_phy <- tax_glom(base_filtercomp, taxrank = "Phylum")

##remove NA taxonomic levels
tax_table(base_comp_phy) <- tax_table(base_comp_phy)[,1:2]

##Remove zero sum
base_comp_phy <- prune_taxa(taxa_sums(base_comp_phy) > 0, base_comp_phy)
##taxa are rows is TRUE

##Baseline phylum level 
#subset most abundant 8 phyla
base_top6phyla = names(sort(taxa_sums(base_comp_phy), TRUE)[1:6])
base_phyla6_subset <- prune_taxa(base_top6phyla, base_comp_phy)

#sort by phylum abundance
phyla_to_sort_b <- data.frame(id=1:6, phyla = as.character(tax_table(base_phyla6_subset)[,"Phylum"]), 
                              otu = as.character(taxa_names(base_phyla6_subset)))
rownames(phyla_to_sort_b) <- phyla_to_sort_b$otu
base_phylum_ranks <- base_phyla6_subset %>% otu_table %>% rowSums %>% sort(T) %>% names
phyla_to_sort_b <- phyla_to_sort_b[base_phylum_ranks, ]

#calculate abundance as a proportion of top 8
base_prop_phylum6 <- transform_sample_counts(base_phyla6_subset, function(i) i / sum(i))

#melt for use in plotting
base_bardat <- psmelt(base_prop_phylum6) %>% 
  mutate(Sample = as.numeric(factor(Sample)),    # geom_area only seems to work with a numeric x
         OTU = factor(OTU,                         # sort by phyla
                      levels=phyla_to_sort_b$otu,   
                      labels=phyla_to_sort_b$phyla)) 

base_firmicutes_order <- base_bardat %>%     # get order by firmicutes abundance
  filter(OTU=="Firmicutes") %>% 
  arrange(Abundance) %>% 
  select(Sample)

base_bardat %<>% mutate(Sample = as.numeric(factor(Sample, levels=factor(base_firmicutes_order$Sample)))) # apply firmicutes order to data.frame

# Check for duplicated Sample IDs
duplicated_samples <- base_bardat[duplicated(base_bardat$Sample), "Sample"]
duplicated_samples

##remove percent sign function
# Custom label formatter function to remove "%" sign
remove_percent_sign <- function(x) {
  as.character(as.numeric(x) * 100)
}

##sort(unique(base_bardat$Phylum))
##"Actinobacteriota"  "Bacteroidota"      "Euryarchaeota"     "Firmicutes"        "Proteobacteria"   
##"Verrucomicrobiota"
 

##Define phylum colors
area_phylum <- c("Actinobacteriota" = "#ee3333", "Bacteroidota" = "#ee7722", "Euryarchaeota" = "#ffee33", 
                 "Firmicutes" = "deepskyblue2","Proteobacteria" = "#3366aa", "Verrucomicrobiota" = "#992288")
base_bardat$Treatment_Type <- factor(base_bardat$Treatment_Type, levels = c("Intervention",
                                                                            "Placebo"))

# Barplot by treatment group
set.seed(14)
dir.create("area_plot")
ggplot(arrange(base_bardat, desc(OTU), Abundance), aes(x = Sample, y = Abundance, fill = OTU)) +
  geom_area() +
  scale_fill_manual(values = area_phylum, breaks = c("Actinobacteriota", "Bacteroidota", 
                                                     "Euryarchaeota", "Firmicutes", 
                                                     "Proteobacteria", "Verrucomicrobiota")) +
  labs(fill = "Phylum") +
  facet_wrap( ~ Treatment_Type, scales = "free_x") + 
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Sample", y = "Relative abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(), legend.position = "right", axis.line = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 10), legend.title.align = 0.5, legend.title = element_blank(),
        legend.text = element_text(size = 10), strip.background = element_blank(), strip.text = element_text(size = 10),
        legend.key.size = unit(0.4, "cm"), legend.margin = margin(t = 0, b = 0, r = 0, l = -3)) +
  scale_y_continuous(expand = c(0, 0), labels = remove_percent_sign)
ggsave("area_plot/finaltop10/top5_phylum_baseline_treatment.tiff", dpi = 600, bg = "white", width = 5.5,height = 3)

base_bardat %>% select(Phylum, Abundance) %>% group_by(Phylum, Abundance) %>% 
  summarise(average_rv = mean(Abundance)) %>% View()

##Overall mean prevalence
base_bardat %>% 
  group_by(Phylum) %>%
  summarize(mean_abundance = mean(Abundance, na.rm = TRUE))
  
##average by treatment groups
base_bardat %>% 
  group_by(Phylum, Treatment_Type) %>%
  summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>% View()

##Plot all samples together
ggplot(arrange(base_bardat, desc(OTU), Abundance), aes(x = Sample, y = Abundance, fill = OTU)) +
  geom_area() +
  scale_fill_manual(values = area_phylum, breaks = c("Actinobacteriota", "Bacteroidota", 
                                                     "Euryarchaeota", "Firmicutes", 
                                                     "Proteobacteria", "Verrucomicrobiota")) +
  labs(fill = "Phylum") +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Sample", y = "Relative abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(), legend.position = "right", axis.line = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 10), legend.title.align = 0.5, legend.title = element_blank(),
        legend.text = element_text(size = 10), strip.background = element_blank(), strip.text = element_text(face = "bold", size = 10),
        legend.key.size = unit(0.4, "cm"), legend.margin = margin(t = 0, b = 0, r = 0, l = -3)) +
  scale_y_continuous(expand = c(0, 0), labels = remove_percent_sign)
ggsave("area_plot/top5_phylum_baseline_all.tiff", dpi = 600, bg = "white")


##Plot at the genus level
base_comp_gen <- tax_glom(base_filtercomp, taxrank = "Genus")

##remove NA taxonomic levels
tax_table(base_comp_gen) <- tax_table(base_comp_gen)[,1:6]

##Remove zero sum
base_comp_gen <- prune_taxa(taxa_sums(base_comp_gen) > 0, base_comp_gen)

base_genusall <- transform_sample_counts(base_comp_gen, function(i) i / sum(i))

df_base_genusall <- psmelt(base_genusall)

##Overall mean prevalence
df_base_genusall %>% 
  group_by(Genus) %>%
  summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>% View()

##average by treatment groups
base_genus_10 %>% 
  group_by(Genus, Treatment_Type) %>%
  summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>% View()





#subset most abundant 8 phyla
top10genus_base = names(sort(taxa_sums(base_comp_gen), TRUE)[1:10])
base_genus10 <- prune_taxa(top10genus_base, base_comp_gen)

#sort by phylum abundance
genus_to_sort <- data.frame(id=1:10, phyla = as.character(tax_table(base_genus10)[,"Genus"]), 
                            otu = as.character(taxa_names(base_genus10)))
rownames(genus_to_sort) <- genus_to_sort$otu
genus_ranks <- base_genus10 %>% otu_table %>% rowSums %>% sort(T) %>% names
genus_to_sort <- genus_to_sort[genus_ranks, ]

#calculate abundance as a proportion of top 8
base_genusprop <- transform_sample_counts(base_genus10, function(i) i / sum(i))

#melt for use in plotting
base_genus_10 <- psmelt(base_genusprop) %>%  
  mutate(Sample = as.numeric(factor(Sample)),    # geom_area only seems to work with a numeric x
         OTU = factor(OTU,                         # sort by phyla
                      levels=genus_to_sort$otu,   
                      labels=genus_to_sort$phyla)) 


bacteroides_order <- base_genus_10 %>%     # get order by firmicutes abundance
  filter(OTU=="Bacteroides") %>% 
  arrange(Abundance) %>% 
  select(Sample)

base_genus_10$Treatment_Type <- factor(base_genus_10$Treatment_Type,
                                       levels = c("Intervention", "Placebo"))

base_genus_10 %<>% mutate(Sample = as.numeric(factor(Sample, levels=factor(bacteroides_order$Sample)))) # apply Bacteroides order to data.frame

area_genus <- c("Agathobacter" = "#bcbd22", "Akkermansia" = "darkorchid", "Clostridia_UCG-014" = "deeppink1", "Bacteroides" = "#66aa55", 
                "Christensenellaceae_R-7_group" = "darkslategray", "Coprococcus" = "gray", "Faecalibacterium" = "#d62728", 
                "Ruminococcus" = "#17becf", "Subdoligranulum" = "#0000ff","UCG-002" = "saddlebrown")
sort(unique(base_genus_10$Genus))

# Modify the labels using str_wrap to create line breaks in the legend
genus_labels <- c("Agathobacter", "Akkermansia", "Bacteroides", 
                  "Christensenellaceae\nR-7_group",  # Insert a line break with \n
                  "Clostridia_UCG-014", "Coprococcus", "Faecalibacterium", 
                  "Ruminococcus", "Subdoligranulum", "UCG-002")


##all samples together
ggplot(arrange(base_genus_10, desc(OTU), Abundance), aes(x = Sample, y = Abundance, fill = OTU)) +
  geom_area() +
  scale_fill_manual(
    values = area_genus, 
    labels = genus_labels,
    breaks = c("Agathobacter", "Akkermansia", "Bacteroides", "Christensenellaceae_R-7_group",
               "Clostridia_UCG-014", "Coprococcus", "Faecalibacterium", "Ruminococcus", 
               "Subdoligranulum", "UCG-002")) +  
  labs(fill = "Genus") +
  facet_wrap(~ Treatment_Type, scales = "free_x") + 
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Sample", y = "Relative abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(), legend.position = "right", axis.line = element_blank(), 
        axis.ticks.x = element_blank(), axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10), axis.text.y = element_text(size = 10, face = "bold"), 
        legend.title = element_blank(), legend.text = element_text(face = "italic"),
        legend.key.height = unit(0.4, "cm"), 
        legend.key.width = unit(0.4, "cm"),
        strip.background = element_blank(), 
        legend.margin = margin(t = 0, b = 0, r = 0, l = -3)) +
  scale_y_continuous(expand = c(0, 0), labels = remove_percent_sign)

ggsave("area_plot/finaltop10/base_top10_genustreatment.tiff", dpi = 600, bg = "white", width = 5.5, height = 3)

##divided by treatment type
ggplot(arrange(base_genus_10, desc(OTU), Abundance), aes(x = Sample, y = Abundance, fill = OTU)) +
  geom_area() +
  scale_fill_manual(values = area_genus, breaks = c("Agathobacter", "Akkermansia", "Bacteroides", "Christensenellaceae_R-7_group",
                                                    "Clostridia_UCG-014", "Coprococcus", "Faecalibacterium", "Ruminococcus", "Subdoligranulum",
                                                    "UCG-002")) +  # Corrected line
  labs(fill = "Genus") +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Sample", y = "Relative abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(), legend.position = "right", axis.line = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10), legend.text = element_text(face = "italic"),
        axis.text.y = element_text(size = 10), legend.title.align = 0.5, legend.title = element_blank(),
        legend.key.size = unit(0.4, "cm"), legend.margin = margin(t = 0, b = 0, r = 0, l = -3)) +
  scale_y_continuous(expand = c(0, 0), labels = remove_percent_sign)
ggsave("area_plot/base_top10_genus_all.tiff", dpi = 600, bg = "white")

##Overall mean prevalence
base_genus_10 %>% 
  group_by(Genus) %>%
  summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>% View()

##average by treatment groups
base_genus_10 %>% 
  group_by(Genus, Treatment_Type) %>%
  summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>% View()




##Subset midline
wk12_filtercomp <- subset_samples(psstudy5_filtercomp, Time == "12-wks")
##1240
##Remove zero sum
wk12_filtercomp <- prune_taxa(taxa_sums(wk12_filtercomp) > 0, wk12_filtercomp) ##1240 ASVs

##Plot at the phylum level
wk12_comp_phy <- tax_glom(wk12_filtercomp, taxrank = "Phylum")

##remove NA taxonomic levels
tax_table(wk12_comp_phy) <- tax_table(wk12_comp_phy)[,1:2]

##Week 12 phylum level 
#subset most abundant 8 phyla
wk12_top6phyla = names(sort(taxa_sums(wk12_comp_phy), TRUE)[1:6])
wk12_phyla6_subset <- prune_taxa(wk12_top6phyla, wk12_comp_phy)

#sort by phylum abundance
phyla_to_sort_wk12 <- data.frame(id=1:6, phyla = as.character(tax_table(wk12_phyla6_subset)[,"Phylum"]), 
                              otu = as.character(taxa_names(wk12_phyla6_subset)))
rownames(phyla_to_sort_wk12) <- phyla_to_sort_wk12$otu
wk12_phylum_ranks <- wk12_phyla6_subset %>% otu_table %>% rowSums %>% sort(T) %>% names
phyla_to_sort_wk12 <- phyla_to_sort_wk12[wk12_phylum_ranks, ]

#calculate abundance as a proportion of top 8
wk12_prop_phylum6 <- transform_sample_counts(wk12_phyla6_subset, function(i) i / sum(i))

#melt for use in plotting
wk12_bardat <- psmelt(wk12_prop_phylum6) %>% 
  mutate(Sample = as.numeric(factor(Sample)),    # geom_area only seems to work with a numeric x
         OTU = factor(OTU,                         # sort by phyla
                      levels=phyla_to_sort_wk12$otu,   
                      labels=phyla_to_sort_wk12$phyla)) 

wk12_firmicutes_order <- wk12_bardat %>%     # get order by firmicutes abundance
  filter(OTU=="Firmicutes") %>% 
  arrange(Abundance) %>% 
  select(Sample)

wk12_bardat %<>% mutate(Sample = as.numeric(factor(Sample, levels=factor(wk12_firmicutes_order$Sample)))) # apply firmicutes order to data.frame

# Check for duplicated Sample IDs
duplicated_samples <- wk12_bardat[duplicated(wk12_bardat$Sample), "Sample"]
duplicated_samples

##remove percent sign function
# Custom label formatter function to remove "%" sign
remove_percent_sign <- function(x) {
  as.character(as.numeric(x) * 100)
}

##Define phylum colors
area_phylum <- c("Actinobacteriota" = "#ee3333", "Bacteroidota" = "#ee7722", "Euryarchaeota" = "#ffee33", 
                 "Firmicutes" = "deepskyblue2","Proteobacteria" = "#3366aa", "Verrucomicrobiota" = "#992288")
wk12_bardat$Treatment_Type <- factor(wk12_bardat$Treatment_Type, levels = c("Placebo",

sort(unique(wk12_bardat$Phylum))
                                                                                                                                                 "Intervention"))
#Barplot by treatment group
set.seed(14)
dir.create("area_plot")
ggplot(arrange(wk12_bardat, desc(OTU), Abundance), aes(x = Sample, y = Abundance, fill = OTU)) +
  geom_area() +
  scale_fill_manual(values = area_phylum, breaks = c("Actinobacteriota", "Bacteroidota", 
                                                     "Euryarchaeota", "Firmicutes", 
                                                     "Proteobacteria", "Verrucomicrobiota")) +
  labs(fill = "Phylum") +
  facet_wrap( ~ Treatment_Type, scales = "free_x") + 
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Sample", y = "Relative abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(), legend.position = "right", axis.line = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 10), legend.title.align = 0.5, legend.title = element_blank(),
        legend.text = element_text(size = 10), strip.background = element_blank(), strip.text = element_text(size = 10),
        legend.key.size = unit(0.4, "cm"), legend.margin = margin(t = 0, b = 0, r = 0, l = -3)) +
  scale_y_continuous(expand = c(0, 0), labels = remove_percent_sign)
ggsave("area_plot/top5_phylum_week12_treatment.tiff", dpi = 600, bg = "white")

##Plot all samples together
ggplot(arrange(wk12_bardat, desc(OTU), Abundance), aes(x = Sample, y = Abundance, fill = OTU)) +
  geom_area() +
  scale_fill_manual(values = area_phylum, breaks = c("Actinobacteriota", "Bacteroidota", 
                                                     "Euryarchaeota", "Firmicutes", 
                                                     "Proteobacteria", "Verrucomicrobiota")) +
  labs(fill = "Phylum") +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Sample", y = "Relative abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(), legend.position = "right", axis.line = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 10), legend.title.align = 0.5, legend.title = element_blank(),
        legend.text = element_text(size = 10), strip.background = element_blank(), strip.text = element_text(face = "bold", size = 10),
        legend.key.size = unit(0.4, "cm"), legend.margin = margin(t = 0, b = 0, r = 0, l = -3)) +
  scale_y_continuous(expand = c(0, 0), labels = remove_percent_sign)
ggsave("area_plot/top5_phylum_week12_all.tiff", dpi = 600, bg = "white")

##Plot at the genus level
wk12_comp_gen <- tax_glom(wk12_filtercomp, taxrank = "Genus")

##remove NA taxonomic levels
tax_table(wk12_comp_gen) <- tax_table(wk12_comp_gen)[,1:6]

##Remove zero sum
wk12_comp_gen <- prune_taxa(taxa_sums(wk12_comp_gen) > 0, wk12_comp_gen)

wk12_genall <- transform_sample_counts(wk12_comp_gen, function(i) i / sum(i))

df_wk12_genall <- psmelt(wk12_genall)

##Overall mean prevalence
df_wk12_genall %>% 
  group_by(Genus) %>%
  summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>% View()

##average by treatment groups
df_wk12_genall %>% 
  group_by(Genus, Treatment_Type) %>%
  summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>% View()




#subset most abundant 8 phyla
top10genus_wk12 = names(sort(taxa_sums(wk12_comp_gen), TRUE)[1:10])
wk12_genus10 <- prune_taxa(top10genus_wk12, wk12_comp_gen)

#sort by phylum abundance
genus_to_sort <- data.frame(id=1:10, phyla = as.character(tax_table(wk12_genus10)[,"Genus"]), 
                            otu = as.character(taxa_names(wk12_genus10)))
rownames(genus_to_sort) <- genus_to_sort$otu
wk12_genus_ranks <- wk12_genus10 %>% otu_table %>% rowSums %>% sort(T) %>% names
genus_to_sort <- genus_to_sort[wk12_genus_ranks, ]

#calculate abundance as a proportion of top 8
wk12_genus10prop <- transform_sample_counts(wk12_genus10, function(i) i / sum(i))

#melt for use in plotting
wk12_genus_10 <- psmelt(wk12_genus10prop) %>%  
  mutate(Sample = as.numeric(factor(Sample)),    # geom_area only seems to work with a numeric x
         OTU = factor(OTU,                         # sort by phyla
                      levels=genus_to_sort$otu,   
                      labels=genus_to_sort$phyla)) 

bacteroides_order <- wk12_genus_10 %>%     # get order by firmicutes abundance
  filter(OTU=="Bacteroides") %>% 
  arrange(Abundance) %>% 
  select(Sample)

sort(unique(wk12_genus_10$Genus))

wk12_genus_10$Treatment_Type <- factor(wk12_genus_10$Treatment_Type,
                                       levels = c("Placebo", "Intervention"))

wk12_genus_10 %<>% mutate(Sample = as.numeric(factor(Sample, levels=factor(bacteroides_order$Sample)))) # apply Bacteroides order to data.frame

area_genus <- c("Agathobacter" = "#bcbd22", "Akkermansia" = "darkorchid", "Clostridia_UCG-014" = "deeppink1", "Bacteroides" = "#66aa55", 
                "Christensenellaceae_R-7_group" = "darkslategray", "Coprococcus" = "gray", "Faecalibacterium" = "#d62728", 
                "Ruminococcus" = "#17becf", "Subdoligranulum" = "#0000ff","UCG-002" = "saddlebrown")
sort(unique(wk12_genus_10$Genus))

##all samples together
ggplot(arrange(wk12_genus_10, desc(OTU), Abundance), aes(x = Sample, y = Abundance, fill = OTU)) +
  geom_area() +
  scale_fill_manual(values = area_genus, breaks = c("Agathobacter", "Akkermansia", "Bacteroides", "Christensenellaceae_R-7_group",
                                                    "Clostridia_UCG-014", "Coprococcus", "Faecalibacterium", "Ruminococcus", "Subdoligranulum",
                                                    "UCG-002")) +  # Corrected line
  labs(fill = "Genus") +
  facet_wrap( ~ Treatment_Type, scales = "free_x") + 
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Sample", y = "Relative abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(), legend.position = "right", axis.line = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 10, face = "bold"), legend.title = element_blank(), legend.text = element_text(face = "italic"),
        legend.key.size = unit(0.4, "cm"),strip.background = element_blank(),legend.margin = margin(t = 0, b = 0, r = 0, l = -3)) +
  scale_y_continuous(expand = c(0, 0), labels = remove_percent_sign)
ggsave("area_plot/wk12_top10_genustreatment.tiff", dpi = 600, bg = "white")

##all samples together
ggplot(arrange(wk12_genus_10, desc(OTU), Abundance), aes(x = Sample, y = Abundance, fill = OTU)) +
  geom_area() +
  scale_fill_manual(values = area_genus, breaks = c("Agathobacter", "Akkermansia", "Bacteroides", "Christensenellaceae_R-7_group",
                                                    "Clostridia_UCG-014", "Coprococcus", "Faecalibacterium", "Ruminococcus", "Subdoligranulum",
                                                    "UCG-002")) +  # Corrected line
  labs(fill = "Genus") +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Sample", y = "Relative abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(), legend.position = "right", axis.line = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 10, face = "bold"), legend.title = element_blank(), legend.text = element_text(face = "italic"),
        legend.key.size = unit(0.4, "cm"),strip.background = element_blank(),legend.margin = margin(t = 0, b = 0, r = 0, l = -3)) +
  scale_y_continuous(expand = c(0, 0), labels = remove_percent_sign)
ggsave("area_plot/wk12_top10_genus_all.tiff", dpi = 600, bg = "white")

##Subset endline
wk24_filtercomp <- subset_samples(psstudy5_filtercomp, Time == "24-wks")

##Plot at the phylum level
wk24_comp_phy <- tax_glom(wk24_filtercomp, taxrank = "Phylum")

##remove NA taxonomic levels
tax_table(wk24_comp_phy) <- tax_table(wk24_comp_phy)[,1:2]

##Week 12 phylum level 
#subset most abundant 8 phyla
wk24_top6phyla = names(sort(taxa_sums(wk24_comp_phy), TRUE)[1:6])
wk24_phyla6_subset <- prune_taxa(wk24_top6phyla, wk24_comp_phy)

#sort by phylum abundance
phyla_to_sort_wk24 <- data.frame(id=1:6, phyla = as.character(tax_table(wk24_phyla6_subset)[,"Phylum"]), 
                                 otu = as.character(taxa_names(wk24_phyla6_subset)))
rownames(phyla_to_sort_wk24) <- phyla_to_sort_wk24$otu
wk24_phylum_ranks <- wk24_phyla6_subset %>% otu_table %>% rowSums %>% sort(T) %>% names
phyla_to_sort_wk24 <- phyla_to_sort_wk24[wk24_phylum_ranks, ]

#calculate abundance as a proportion of top 8
wk24_prop_phylum6 <- transform_sample_counts(wk24_phyla6_subset, function(i) i / sum(i))

#melt for use in plotting
wk24_bardat <- psmelt(wk24_prop_phylum6) %>% 
  mutate(Sample = as.numeric(factor(Sample)),    # geom_area only seems to work with a numeric x
         OTU = factor(OTU,                         # sort by phyla
                      levels=phyla_to_sort_wk24$otu,   
                      labels=phyla_to_sort_wk24$phyla)) 
wk24_firmicutes_order <- wk24_bardat %>%     # get order by firmicutes abundance
  filter(OTU=="Firmicutes") %>% 
  arrange(Abundance) %>% 
  dplyr::select(Sample)

wk24_bardat %<>% mutate(Sample = as.numeric(factor(Sample, levels=factor(wk24_firmicutes_order$Sample)))) # apply firmicutes order to data.frame

# Check for duplicated Sample IDs
duplicated_samples <- wk24_bardat[duplicated(wk24_bardat$Sample), "Sample"]
duplicated_samples

##remove percent sign function
# Custom label formatter function to remove "%" sign
remove_percent_sign <- function(x) {
  as.character(as.numeric(x) * 100)
}

##Define phylum colors
area_phylum <- c("Actinobacteriota" = "#ee3333", "Bacteroidota" = "#ee7722", "Euryarchaeota" = "#ffee33", 
                 "Firmicutes" = "deepskyblue2","Proteobacteria" = "#3366aa", "Verrucomicrobiota" = "#992288")
wk24_bardat$Treatment_Type <- factor(wk24_bardat$Treatment_Type, levels = c("Intervention",sort(unique(wk12_bardat$Phylum)),
                                                                            "Placebo"))


#Barplot by treatment group
set.seed(14)
dir.create("area_plot")
ggplot(arrange(wk24_bardat, desc(OTU), Abundance), aes(x = Sample, y = Abundance, fill = OTU)) +
  geom_area() +
  scale_fill_manual(values = area_phylum, breaks = c("Actinobacteriota", "Bacteroidota", 
                                                     "Euryarchaeota", "Firmicutes", 
                                                     "Proteobacteria", "Verrucomicrobiota")) +
  labs(fill = "Phylum") +
  facet_wrap( ~ Treatment_Type, scales = "free_x") + 
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Sample", y = "Relative abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(), legend.position = "right", axis.line = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 10), legend.title.align = 0.5, legend.title = element_blank(),
        legend.text = element_text(size = 10), strip.background = element_blank(), strip.text = element_text(size = 10),
        legend.key.size = unit(0.4, "cm"), legend.margin = margin(t = 0, b = 0, r = 0, l = -3)) +
  scale_y_continuous(expand = c(0, 0), labels = remove_percent_sign)
ggsave("area_plot/finaltop10/top5_phylum_week24_treatment.tiff", dpi = 600, bg = "white", width = 5.5, height = 3)

##Phylum level for all treatment type
ggplot(arrange(wk24_bardat, desc(OTU), Abundance), aes(x = Sample, y = Abundance, fill = OTU)) +
  geom_area() +
  scale_fill_manual(values = area_phylum, breaks = c("Actinobacteriota", "Bacteroidota", 
                                                     "Euryarchaeota", "Firmicutes", 
                                                     "Proteobacteria", "Verrucomicrobiota")) +
  labs(fill = "Phylum") +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Sample", y = "Relative abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(), legend.position = "right", axis.line = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 10), legend.title.align = 0.5, legend.title = element_blank(),
        legend.text = element_text(size = 10), strip.background = element_blank(), strip.text = element_text(size = 10),
        legend.key.size = unit(0.4, "cm"), legend.margin = margin(t = 0, b = 0, r = 0, l = -3)) +
  scale_y_continuous(expand = c(0, 0), labels = remove_percent_sign)
ggsave("area_plot/top5_phylum_week24_all.tiff", dpi = 600, bg = "white")


##Plot at the genus level
wk24_comp_gen <- tax_glom(wk24_filtercomp, taxrank = "Genus")

##Remove zero sum
wk24_comp_gen <- prune_taxa(taxa_sums(wk24_comp_gen) > 0, wk24_comp_gen)

wk24_genall <- transform_sample_counts(wk24_comp_gen, function(i) i / sum(i))

df_wk24_genall <- psmelt(wk24_genall)

##Overall mean prevalence
df_wk24_genall %>% 
  group_by(Genus) %>%
  summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>% View()

##average by treatment groups
df_wk12_genall %>% 
  group_by(Genus, Treatment_Type) %>%
  summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>% View()


#subset most abundant 8 phyla
top10genus_wk24 = names(sort(taxa_sums(wk24_comp_gen), TRUE)[1:10])
wk24_genus10 <- prune_taxa(top10genus_wk24, wk24_comp_gen)

#sort by phylum abundance
wk24genus_to_sort <- data.frame(id=1:10, phyla = as.character(tax_table(wk24_genus10)[,"Genus"]), 
                            otu = as.character(taxa_names(wk24_genus10)))
rownames(wk24genus_to_sort) <- wk24genus_to_sort$otu
wk24genus_ranks <- wk24_genus10 %>% otu_table %>% rowSums %>% sort(T) %>% names
wk24genus_to_sort <- wk24genus_to_sort[wk24genus_ranks, ]

#calculate abundance as a proportion of top 8
wk24_genusprop <- transform_sample_counts(wk24_genus10, function(i) i / sum(i))

#melt for use in plotting
wk24_bardat <- psmelt(wk24_genusprop) %>%  
  mutate(Sample = as.numeric(factor(Sample)),    # geom_area only seems to work with a numeric x
         OTU = factor(OTU,                         # sort by phyla
                      levels=wk24genus_to_sort$otu,   
                      labels=wk24genus_to_sort$phyla)) 


bacteroides_order <- wk24_bardat %>%     # get order by firmicutes abundance
  filter(OTU=="Bacteroides") %>% 
  arrange(Abundance) %>% 
  dplyr::select(Sample)

wk24_bardat$Treatment_Type <- factor(wk24_bardat$Treatment_Type,
                                       levels = c("Intervention", "Placebo"))


wk24_bardat %<>% mutate(Sample = as.numeric(factor(Sample, levels=factor(bacteroides_order$Sample)))) # apply Bacteroides order to data.frame

area_genus <- c("Agathobacter" = "#bcbd22", "Akkermansia" = "darkorchid", "Clostridia_UCG-014" = "deeppink1", "Bacteroides" = "#66aa55", 
                "Christensenellaceae_R-7_group" = "darkslategray", "Coprococcus" = "gray", "Faecalibacterium" = "#d62728", 
                "Ruminococcus" = "#17becf", "Subdoligranulum" = "#0000ff","UCG-002" = "saddlebrown")

# Modify the labels using str_wrap to create line breaks in the legend
genus_labels <- c("Agathobacter", "Akkermansia", "Bacteroides", 
                  "Christensenellaceae\nR-7_group",  # Insert a line break with \n
                  "Clostridia_UCG-014", "Coprococcus", "Faecalibacterium", 
                  "Ruminococcus", "Subdoligranulum", "UCG-002")

##all samples together
ggplot(arrange(wk24_bardat, desc(OTU), Abundance), aes(x = Sample, y = Abundance, fill = OTU)) +
  geom_area() +
  scale_fill_manual(values = area_genus, 
                    labels = genus_labels,
                    breaks = c("Agathobacter", "Akkermansia", "Bacteroides", "Christensenellaceae_R-7_group",
                                                    "Clostridia_UCG-014", "Coprococcus", "Faecalibacterium", "Ruminococcus", "Subdoligranulum",
                                                    "UCG-002")) +  # Corrected line
  labs(fill = "Genus") +
  facet_wrap( ~ Treatment_Type, scales = "free_x") + 
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Sample", y = "Relative abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(), legend.position = "right", axis.line = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 10, face = "bold"), legend.title = element_blank(), legend.text = element_text(face = "italic"),
        legend.key.size = unit(0.4, "cm"),strip.background = element_blank(),legend.margin = margin(t = 0, b = 0, r = 0, l = -3)) +
  scale_y_continuous(expand = c(0, 0), labels = remove_percent_sign)
ggsave("area_plot/finaltop10/wk24_top10_genustreatment.tiff", dpi = 600, bg = "white", width = 5.5, height = 3)

##all samples together
ggplot(arrange(wk24_bardat, desc(OTU), Abundance), aes(x = Sample, y = Abundance, fill = OTU)) +
  geom_area() +
  scale_fill_manual(values = area_genus, breaks = c("Agathobacter", "Akkermansia", "Bacteroides", "Christensenellaceae_R-7_group",
                                                    "Clostridia_UCG-014", "Coprococcus", "Faecalibacterium", "Ruminococcus", "Subdoligranulum",
                                                    "UCG-002")) +  # Corrected line
  labs(fill = "Genus") +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Sample", y = "Relative abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(), legend.position = "right", axis.line = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 10, face = "bold"), legend.title = element_blank(), legend.text = element_text(face = "italic"),
        legend.key.size = unit(0.4, "cm"),strip.background = element_blank(),legend.margin = margin(t = 0, b = 0, r = 0, l = -3)) +
  scale_y_continuous(expand = c(0, 0), labels = remove_percent_sign)
ggsave("area_plot/wk24_top10_genus_all.tiff", dpi = 600, bg = "white")
