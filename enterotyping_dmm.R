library(DirichletMultinomial)
library(microViz)
library(ggplot2)
library(scales)
library(rstatix)

ps_study810_rare  ##rarefied phyloseq object

ps_study810_rare_f2 <- ps_join(ps_study810_rare_f2, metadatapaired, by = "sampleid")


sample_data(ps_study810_rare_f2) <- data.frame(sample_data(ps_study810_rare_f2)) %>% dplyr::select(-c("Enterotype_BIC...3.x", "Enterotype_BIC...4.x","Baseline_enterotype_BIC.y",
                                                               "Baseline_enterotype_BIC.x","Enterotype_BIC...3.y","Enterotype_BIC...4.y","Enterotype_BIC...4.x"))

##Glom at the genus level
ps810rare_genus <- tax_glom(ps_study810_rare, taxrank = "Genus")

###change to data frame
ps810raregenus_count <-  psmelt(ps810rare_genus)  ##wide data.frame

##Select sample, abundance, and genus
long_study810rare_count <- ps810raregenus_count %>% dplyr::select(Sample, Genus, Abundance)   ##data.frame
View(long_study810rare_count)

rownames(long_study810rare_count) <- NULL
long_study810rare_count <- as_tibble(long_study810rare_count)

# ##Change to wide and transpose ------------------------------------------
# Adding a unique identifier for each Genus occurrence within each Sample
long_study810rare_count_1 <- long_study810rare_count %>%
  group_by(Sample, Genus) %>%
  mutate(id = row_number()) %>%
  ungroup()

# Pivoting the data to wide format
wide_study810rare_count <- long_study810rare_count_1 %>%
  pivot_wider(
    names_from = c(Genus, id),
    values_from = Abundance,
    names_sep = "_",
    id_cols = Sample
  )
# Removing '_id' and numeric suffixes from column names if they were generated during the pivot
names(wide_study810rare_count) <- gsub("_\\d+$", "", names(wide_study810rare_count)) # Removes any '_id' and following digits

colnames(wide_study810rare_count)[1] <- "" ###rename sample on the first column names with ""

##change to matrix
wide_study810rare_count

##Convert to dataframe
wide_study810rare_count <- data.frame(wide_study810rare_count)

#Set the first column as row names
rownames(wide_study810rare_count) <- wide_study810rare_count[,1] ##change id column to rownames

#Drop the duplicated sampleid column
wide_study810rare_count <- wide_study810rare_count[-1]

#Change to matrix
wide_study810rare_count <- as.matrix(wide_study810rare_count) 

# Add pseudocount to the OTU table
count_study810rare_final <- wide_study810rare_count + 0.000000001

##Distribution of reads from each taxon, on a log scale
cnts <- log10(colSums(count_study810rare_final))
densityplot(cnts, xlim=range(cnts), xlab="Taxon representation (log 10 count)")

# Compute dirichlet multinomial for different numbers of components
all_dmns <- 6 # how many DMMs should be checked?
dmn_list <- numeric(all_dmns)
for (i in 1:all_dmns) {
  print(i)
  assign(paste0("dmn_", i), dmn(as.matrix(count_study810rare_final), i, verbose = F))
}
dmn_list <- list(dmn_1, dmn_2, dmn_3, dmn_4, dmn_5, dmn_6)

##the effect of alpha diversity/

# Find the DMM with the minimum number of Dirichlet components based on laplace
lplc <- sapply(dmn_list, laplace)
BIC <- sapply(dmn_list, BIC)
AIC <- sapply(dmn_list, AIC)
dmn_min_lplc <- dmn_list[[which.min(lplc)]]
dmn_min_BIC <- dmn_list[[which.min(BIC)]]
dmn_min_AIC <- dmn_list[[which.min(AIC)]]
# Print the description
print(dmn_min_lplc)
#class: DMN 
#k: 3 
#samples x taxa: 283 x 295 
#Laplace: 159364.4 BIC: 161669.7 AIC: 160053
print(dmn_min_BIC)
#class: DMN 
#k: 2 
#samples x taxa: 283 x 295 
#Laplace: 159891.3 BIC: 161328 AIC: 160250.8 
print(dmn_min_AIC)
#class: DMN 
#k: 3 
#samples x taxa: 283 x 295 
#Laplace: 159364.4 BIC: 161669.7 AIC: 160053

# Plot model fit for different numbers of Dirichlet components
pdf("blueberry_dirichlet_components.pdf", onefile = TRUE)
plot(lplc, type = "b", xlab = "Number of Dirichlet Components", ylab = "Model Fit, Laplace") # 3
plot(BIC, type = "b", xlab = "Number of Dirichlet Components", ylab = "Model Fit, BIC") # 2
plot(AIC, type = "b", xlab = "Number of Dirichlet Components", ylab = "Model Fit, AIC") # 3
dev.off()

# Extract the DMM results for different numbers of components
dmn2 <- dmn_list[[2]]
Dirichlet_multinomial_2 <- mixture(dmn2, assign = TRUE)

dmn3 <- dmn_list[[3]]
Dirichlet_multinomial_3 <- mixture(dmn3, assign = TRUE)

# Combine the DMM results into one dataframe
Dirichlet_multinomial_all <- cbind(Dirichlet_multinomial_2, Dirichlet_multinomial_3)
colnames(Dirichlet_multinomial_all) <- c("k2", "k3")

##Save the dirichlet multinomial all
write.csv(Dirichlet_multinomial_all, "output_results/Dirichlet_multinomial_all.csv")

fit <- lapply(1:3, dmn, count = count_study810rare_final, verbose=TRUE)
lplc <- base::sapply(fit, DirichletMultinomial::laplace) # AIC / BIC / Laplace
best <- fit[[which.min(unlist(lplc))]]
ass <- apply(mixture(best), 1, which.max)

library(reshape2)
library(ggplot2)

# Create an empty list to store the plots
d <- reshape2::melt(fitted(best))
colnames(d) <- c("OTU", "cluster", "value")

cluster1 <- subset(d, cluster == "1") %>% 
  arrange(value) %>%
  mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
  # Only show the most important drivers
  filter(abs(value) > quantile(abs(value), 0.8))  
  
cluster1 <- cluster1 %>% arrange(desc(value)) %>% 
  filter(value > 1.663406)
ggplot(cluster1, aes(x = OTU, y = value)) +
    geom_bar(stat = "identity", fill = "#1B9E77") +
    coord_flip() +
    labs(title = paste("Top drivers: Enterotype", "one", x = "")) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0))
ggsave("topdriversenterotypeone.pdf",dpi = 600,width = 6,height = 4)


cluster2 <- subset(d, cluster == "2") %>% 
  arrange(value) %>%
  mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
  # Only show the most important drivers
  filter(abs(value) > quantile(abs(value), 0.8)) 

cluster2 <- cluster2 %>% arrange(desc(value)) %>% 
  filter(value > 1.53)

ggplot(cluster2, aes(x = OTU, y = value)) +
  geom_bar(stat = "identity", fill = "#D95F02") +
  coord_flip() +
  labs(title = paste("Top drivers: Enterotype", "one", x = "")) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0))
ggsave("topdriversenterotypetwo.pdf",dpi = 600,width = 6,height = 4)

cluster3 <- subset(d, cluster == "3") %>% 
  arrange(value) %>%
  mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
  # Only show the most important drivers
  filter(abs(value) > quantile(abs(value), 0.8)) 

cluster3 <- cluster3 %>% arrange(desc(value)) %>% 
  filter(value > 0.4741727)

ggplot(cluster3, aes(x = OTU, y = value)) +
  geom_bar(stat = "identity", fill = "#006BA3") +
  coord_flip() +
  labs(title = paste("Top drivers: Enterotype", "three", x = "")) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0))
ggsave("topdriversenterotypethree.pdf",dpi = 600,width = 6,height = 4)


##Add Enterotype to the phyloseq object
##Read enterotype assignment file
enterotype_2024_final <- read_csv("enterotype_2024_final2.csv")
View(enterotype_2024_final)

#change sample names to screening code
sample_names(enterotype_2024_final) <- enterotype_2024_final$sampleid

##Add enterotype to rarefie phyloseq object
ps_study810_rare_f <- ps_join(ps_study810_rare, enterotype_2024_final, match_sample_names = "sampleid")

##Remove variables from phyloseq object (Screening code and ...1)
sample_data(ps_rarefied_nte) <- sample_data(ps_rarefied_nte)[, -c(1,2)]

##Final phyloseq object
ps_study810_rare_f2
df_study810raref2 <- data.frame(sample_data(ps_study810_rare_f2))

##Enterotype levels
df_study810raref2$Time <- factor(df_study810raref2$Time, levels = c("Baseline", "12-wks",
                                                                    "24-wks"))
df_study810raref2$Enterotype <- factor(df_study810raref2$Enterotype, levels = c("Enterotype one",
                                                                                "Enterotype two",
                                                                                "Enterotype three"))

##plot the proprotion of enterotype at each timepoint
# Plot with integrated proportion calculation and text labels
df_study810raref2 %>%
  group_by(Time, Enterotype) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(Time) %>%
  mutate(Proportion = Count / sum(Count)) %>% 
  ggplot(aes(x = Time, fill = Enterotype)) +  
  geom_bar(position = "fill", aes(y = after_stat(prop)), stat = "count") +  # Now using after_stat()
  labs(y = "", x = "") +
  theme_bw() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.title = element_blank(),plot.margin = margin(t=0,r=0,b=0,l=0)) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#006BA3")) +
  geom_text(aes(label = scales::percent(Proportion, accuracy = 1), y = Proportion),
            position = position_stack(vjust = 0.5), stat = "identity", show.legend = FALSE)
ggsave("output_results/enterotype_time.tiff", dpi = 600)

# First, calculate the proportions manually
df_enterotype <- df_study810raref2 %>%
  group_by(Time, Enterotype) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(Time) %>%
  mutate(Proportion = Count / sum(Count),
         CumulativeProp = cumsum(Proportion) - (Proportion / 2)) %>%  # Calculate the midpoint of each stack
  ungroup()

# Then, plot the data
ggplot(df_enterotype, aes(x = Time, y = Proportion, fill = Enterotype)) + 
  geom_bar(position = "fill", stat = "identity") +  # Use manually calculated proportions
  labs(y = "") +
  theme_bw() +
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#006BA3")) +
  geom_text(aes(label = scales::percent(Proportion, accuracy = 1), y = CumulativeProp),
            position = position_fill(vjust = 0.5), show.legend = FALSE, color = "white")  # You can change color as needed


###Kruskal wallis test at the Baseline between Enterotypes
##Subset baseline 

psstudy5_filtercomp_1 <- ps_join(psstudy5_filtercomp, metadatapaired, by = "sampleid")

##
base_filtercomp <- subset_samples(psstudy5_filtercomp, Time == "Baseline")

##Glom at the genus level
base_filtercompGenus <- tax_glom(base_filtercomp, taxrank = "Genus")
sum(taxa_sums(base_filtercomp) == 0) ##no zero sum

##Kurskal wallis test for the difference between Enterotypes
enterokruskal_base <- data.frame(t(otu_table(base_filtercompGenus)))
enterokruskal_base$Enterotype <- phyloseq::sample_data(base_filtercompGenus)$Enterotype
enterokruskal_base$Enterotype <- factor(enterokruskal_base$Enterotype, c("Enterotype one", "Enterotype two", "Enterotype three"))

##Adding taxonomic labels
base_taxa <- data.frame(tax_table(base_filtercompGenus))
base_taxa <- base_taxa %>% select(-Species) %>% rownames_to_column(var = "ASV")

##Relative abundance per ASV##RelatSpeciesive abundance per ASV
rel_kruskalbase <- enterokruskal_base %>% 
  gather(key = ASV, value = relabun, -Enterotype) %>% 
  group_by(ASV)

mean_kruskalbase <- rel_kruskalbase %>% 
  dplyr::group_by(Enterotype, ASV) %>% 
  summarise(mean_relabun = mean(relabun))
write.csv(mean_kruskalbase, "output_results/mean_kruskalbase.csv")

##add taxa_info to the mean_kruskalbase
basefinal_meankruskal <- full_join(mean_kruskalbase, base_taxa)
write.csv(basefinal_meankruskal, "output_results/basefinal_meankruskal.csv")

write.csv(rel_kruskalbase, "output_results/rel_kruskalbase.csv")

##Kruskal wallis test of Enterotype
kruskalentero_basebh <- rel_kruskalbase %>% 
  group_by(ASV) %>% 
  kruskal_test(relabun~Enterotype) %>% 
  adjust_pvalue(method = "BH") %>% 
  add_significance()

##Kruskal wallis without p value adjustment
kruskalentero_base <- rel_kruskalbase %>% 
  group_by(ASV) %>% 
  kruskal_test(relabun~Enterotype) %>% 
  add_significance()

##Test only prevotella
prevotella_enterotype <- read.csv("prevotella_enterotype.csv")
prevotella_enterotype$relabun <- as.numeric(prevotella_enterotype$relabun)
prevotella_enterotype$Enterotype <- factor(prevotella_enterotype$Enterotype,
                                           c("Enterotype one", "Enterotype two",
                                             "Enterotype three"))
prevoentero_base <- prevotella_enterotype %>% 
  group_by(ASV) %>% 
  kruskal_test(relabun~Enterotype) %>% 
  add_significance()

##filter the significant results and use dunn posthoc
sig_kruskalenterobase <- kruskalentero_base %>% filter(p.signif != "ns")

##Subset for only significant from the original
# Filter rel_kruskalbase to only include ASVs in sig_kruskalenterobase
sig_relkruskalbase <- rel_kruskalbase %>%
  semi_join(sig_kruskalenterobase, by = "ASV")

##Posthoc dunn test
posthoc_dunnentero <- sig_relkruskalbase %>% 
  group_by(ASV) %>% 
  dunn_test(relabun ~ Enterotype, p.adjust.method = "BH")%>% 
  add_significance()
write.csv(posthoc_dunnentero, "output_results/posthoc_dunnentero.csv")

##filter only signficant results
sig_posthoc_dunnentero<- posthoc_dunnentero %>% filter(p.adj.signif != "ns")

sig_posthoc_dunnentero$ASV %in% basefinal_meankruskal$ASV
baseline_dunnentero <- inner_join(sig_posthoc_dunnentero, base_taxa, by = "ASV", 
          relationship = "many-to-many")
write.csv(baseline_dunnentero, "output_results/baseline_dunnentero.csv")

##Show only unique asvs
uniqueasv_enterobase <- unique(baseline_dunnentero$ASV)
sigasv_relativeabund <- mean_kruskalbase %>% filter(ASV %in% uniqueasv_enterobase)

##Add taxonomic information
tax_sigasv_relativeabund <- inner_join(sigasv_relativeabund, base_taxa, by = "ASV")

write.csv(tax_sigasv_relativeabund, "output_results/tax_sigasv_relativeabund.csv")

dunn_entero_base <- ps_kruskal_base %>% 
  gather(key = ASV, value = relabun, -Enterotype) %>% 
  group_by(ASV) %>%
  dunn_test(relabun ~ Enterotype, p.adjust.method = "BH")%>% 
  adjust_pvalue(method = "BH") %>%
  add_significance()

df_base_filtercompGenus <- psmelt(base_filtercompGenus)
sele_df_base_filtercomp <- df_base_filtercompGenus %>% select(OTU, Sample,Genus,Treatment_Type, Enterotype, Abundance)

##Factor
sele_df_base_filtercomp$Enterotype <- factor(sele_df_base_filtercomp$Enterotype,
                                             levels = c("Enterotype one", "Enterotype two",
                                                        "Enterotype three"))


sele_df_base_filtercomp %>% 
  select(OTU, Genus, Abundance, Enterotype) %>% 
  dplyr::group_by(Enterotype, OTU, Genus) %>% 
  summarise(mean_relabun = mean(Abundance), .groups = 'drop')


##Unnest the results, grab the ASV names and p-values, add the taxonomic labels, and calculate the FDR adjusted
##p-values
#Unnesting
kruskal_results_1 <- kruskal_results %>%
  dplyr::select(ASV, p_value) %>%
  tidyr::unnest(cols = c(p_value))

##check the data
head(kruskal_results_1)

##plot enterotypes
##first merge 
tax_table(psstudy5_filtercomp)
##Glom at the genus level
ps_filtercompGenus <- tax_glom(psstudy5_filtercomp, taxrank = "Genus")
View(tax_table(ps_filtercompGenus))
##Remove species from the tax table
tax_table(ps_filtercompGenus) <- tax_table(ps_filtercompGenus)[,1:6]
ps_prevotella <- ps_filtercompGenus

tax_ps_filtercompGenus <- data.frame(tax_table(ps_filtercompGenus))
tax_ps_filtercompGenus <- tax_ps_filtercompGenus %>% rownames_to_column(var = "ASV") %>% select(-Species)
View(tax_ps_filtercompGenus)
##Extract OTU table
otu_prevotella <- as.data.frame(otu_table(ps_prevotella))

# Specify the ASVs to rename
asvs_to_rename <- c("ASV13_126", "ASV13_65", "ASV13_139", "ASV13_1333")

dim(otu_prevotella)
#198 283

# Rename the ASVs to a common name
otu_prevotella_1 <- otu_prevotella %>%
  rownames_to_column("ASV") %>%
  mutate(ASV = ifelse(ASV %in% asvs_to_rename, "ASV13_126", ASV)) %>%
  group_by(ASV) %>%
  summarise_all(sum) %>%
  column_to_rownames("ASV")

#write.csv(otu_prevotella_1, "otu_prevotella_1.csv")
#write.csv(otu_prevotella, "otu_prevotella.csv")

prevotella_taxa <- data.frame(tax_table(ps_prevotella))
dim(prevotella_taxa)

##Rename the ASV and Genus in the tax_table
prevotella_taxa_1 <- prevotella_taxa %>%
  rownames_to_column("ASV") %>%
  mutate(ASV = ifelse(ASV %in% asvs_to_rename, "ASV13_126", ASV)) %>%
  mutate(Genus = ifelse(Genus %in% c("Paraprevotella",
                                     "Prevotella_9", "Prevotella_7"), "Prevotella", Genus)) %>% 
  distinct(ASV, .keep_all = TRUE) %>% 
  column_to_rownames("ASV")
  
##change the otu table and tax_table
otu_prevotella_1matrix <- as.matrix(otu_prevotella_1)
prevotella_taxa_1matrix <- as.matrix(prevotella_taxa_1)

otu_table(ps_prevotella) <- otu_table(otu_prevotella_1matrix, taxa_are_rows = TRUE)
tax_table(ps_prevotella) <- tax_table(prevotella_taxa_1matrix)

##ps_prevotella is ready for plotting
##subset for baseline
base_prevotella <- subset_samples(ps_prevotella, Time == "Baseline")

##Remove zero sum
base_prevotella <- prune_taxa(taxa_sums(base_prevotella) > 0, base_prevotella) ##195 ASVs

df_baseprevotella <- psmelt(base_prevotella)

##Target bacteria for plotting
target_bacteria <- c("Faecalibacterium", "Bacteroides", "Akkermansia", "Alistipes", "Prevotella", "Methanobrevibacter", "Ruminococcus",
                     "Subdoligranulum")

##Rename everything to other
df_baseprevotella$Genus <- ifelse(df_baseprevotella$Genus %in% target_bacteria, df_baseprevotella$Genus, "Other")

other_baseprevotella <- df_baseprevotella %>% filter(Genus == "Other")

other_baseprevotella1 <- other_baseprevotella %>% 
  group_by(Sample, Genus) %>%
  summarise(Abundance = sum(Abundance))

##add sample metadata
df_base_prevotella <- data.frame(sample_data(base_prevotella))

##Rename sampleid to Sample
df_base_prevotella <- df_base_prevotella %>% rename(Sample = sampleid)

##Join sample data to the Other dataframe
other_baselinecomplete <- right_join(other_baseprevotella1, df_base_prevotella, by = "Sample")

other_baselinecomplete1 <- other_baselinecomplete %>% 
  select(c("Sample","Abundance","Genus","Enterotype"))
notother_baseprevotella <- df_baseprevotella %>% filter(Genus != "Other")

notother_baseprevotella1 <- notother_baseprevotella %>% 
  select(c("Sample","Abundance","Genus","Enterotype"))

entero_baselinerel <- bind_rows(notother_baseprevotella1,other_baselinecomplete1)

##Factorize enterotype
entero_baselinerel$Enterotype <- factor(entero_baselinerel$Enterotype,
                                        levels = c("Enterotype one", "Enterotype two",
                                                   "Enterotype three"))
entero_baselinerel <- entero_baselinerel %>% 
  mutate(Genus = factor(Genus, levels = c("Other", "Methanobrevibacter", "Alistipes",
                                          "Ruminococcus", "Faecalibacterium", 
                                          "Akkermansia","Subdoligranulum", "Prevotella",
                                          "Bacteroides")))



entgenus_color <- c("Other" = "#a5cee2", "Methanobrevibacter" = "#6a3d9a", "Alistipes" = "#1f78b4",
                    "Ruminococcus" = "#b2d98d", "Faecalibacterium" = "#de1e1e", "Subdoligranulum" = "#ffee33",
                    "Akkermansia" = "#33a42b","Prevotella" = "#fc9a99",
                    "Bacteroides" = "#fdbe6f")
##plot stacked bar plot
ggplot(entero_baselinerel, aes(fill = Genus, y=Abundance, x=Enterotype)) +
  geom_bar(position = "fill",stat = "identity") +
  scale_fill_manual(values = entgenus_color, breaks = c("Akkermansia", 
                                                        "Alistipes", 
                                                        "Bacteroides", "Faecalibacterium", "Ruminococcus", "Methanobrevibacter", "Prevotella",
                                                        "Subdoligranulum", "Other")) +
  labs(y = "Mean Relative Abundance (%)", x = "") +
  theme_classic() +
  theme(legend.title = element_blank(), axis.line = element_blank(),
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(face = "italic"), legend.key.size = unit(0.4, "cm"),
        legend.margin = margin(t = 0, b = 0, r = 0, l = -3)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0))
ggsave("area_plot/baseline_enterotype.pdf", dpi = 600, bg = "white",
       width = 6, height = 4)
ggsave(plot = enterotype_all,"enterotype_all.pdf",dpi = 600,width = 6,height = 4)

##Beta diversity of Enterotypes

##Phyloseq object
ps_study810_rare_f2

##gloom taxa at the genus level
psstudy810rare_f2genus <- tax_glom(ps_study810_rare_f2, taxrank = "Genus")

##Extract the baseline
baseline_f2genus <- subset_samples(psstudy810rare_f2genus, Time == "Baseline") ##295

baseline_f2genus <- prune_taxa(taxa_sums(baseline_f2genus) > 0, baseline_f2genus)  ##283

##Calculate relative abundance
baselinef2genus_comp <- microbiome::transform(baseline_f2genus, transform = "compositional")

##Baseline
df_baselinef2genus_comp <- data.frame(sample_data(baselinef2genus_comp))

#Bray-Curtis at genus level
basef2genus_bray <- phyloseq::distance(baselinef2genus_comp, method = "bray")

set.seed(1988)
braygenus_enterotype <- adonis2(formula = basef2genus_bray ~ Enterotype, data = df_baselinef2genus_comp, permutations = 999,
                               na.action = na.exclude)

##Pairwise comparison
bray_enterotype_pair <- pairwise.adonis2(basef2genus_bray ~ Enterotype, data=df_baselinef2genus_comp)

#Dispersion test and plot
braydispr_enterotype <- vegan::betadisper(basef2genus_bray,df_baselinef2genus_comp$Enterotype)
anova(braydispr_enterotype)
###Weighted UniFrac
baseline_f2genus

#Bray-Curtis at genus level
basef2genus_wunifrac <- phyloseq::distance(baseline_f2genus, method = "wunifrac")

set.seed(1988)
wunifrac_enterotype <- adonis2(formula = basef2genus_wunifrac ~ Enterotype, data = df_baselinef2genus_comp, permutations = 999,
                                na.action = na.exclude)

##Pairwise comparison
wunifrac_enterotype_pair <- pairwise.adonis2(basef2genus_wunifrac ~ Enterotype, data=df_baselinef2genus_comp)

##plot everything together
##gloom taxa at the genus level
psstudy810rare_f2genus <- tax_glom(ps_study810_rare_f2, taxrank = "Genus")

##Calculate relative abundance
psstudyf2genus_comp <- microbiome::transform(psstudy810rare_f2genus, transform = "compositional")

##Baseline
df_psstudyf2genus_comp <- data.frame(sample_data(psstudyf2genus_comp))

#Bray-Curtis at genus level
f2genus_bray <- phyloseq::distance(psstudyf2genus_comp, method = "bray")

##
set.seed(1988)
f2braygenus_enterotype <- adonis2(formula = f2genus_bray ~ Enterotype, data = df_psstudyf2genus_comp, permutations = 999,
                                na.action = na.exclude, strata = df_psstudyf2genus_comp$Time)

##Pairwise comparison
f2bray_enterotype_pair <- pairwise.adonis2(f2genus_bray ~ Enterotype, data=df_psstudyf2genus_comp)

#Bray-Curtis at genus level
f2asv_wunifrac <- phyloseq::distance(psstudy810rare_f2genus, method = "wunifrac")

##
set.seed(1988)
f2asvwunifrac_enterotype <- adonis2(formula = f2asv_wunifrac ~ Enterotype, data = df_psstudyf2genus_comp, permutations = 999,
                                  na.action = na.exclude, strata = df_psstudyf2genus_comp$Time)

##Pairwise comparison
f2asvwunifrac_enterotype_pair <- pairwise.adonis2(f2asv_wunifrac ~ Enterotype, data=df_psstudyf2genus_comp)

##plot bray curtis
f2genus_bray

pcoa <- pcoa(f2genus_bray,correction = "none", rn = NULL) 

groups <- df_psstudyf2genus_comp ##sample data
PC1 <- pcoa$vectors[,1]
PC2 <- pcoa$vectors[,2]

groups<-groups[match(rownames(pcoa$vectors),rownames(groups)),]
pcoadata <- data.frame(rownames(pcoa$vectors),
                       PC1,PC2,groups$Enterotype, groups$Time)
colnames(pcoadata) <- c("sample","PC1","PC2","group", "group_1")

pcoadata$group <- factor(pcoadata$group,levels = c("Enterotype one","Enterotype two","Enterotype three"))

pcoadata$group_1 <- factor(pcoadata$group_1,levels = c("Baseline","12-wks","24-wks"))

#Adonis test
set.seed(1988)
otu.adonis <- adonis2(base_rare_wunifrac~group,data = pcoadata,permutations = 999)
write.table(otu.adonis,'beta_bray-curtis_adonis.tsv',sep = '\t',quote = F)

#Boxplot
df <- pcoadata
df1 <- df %>% group_by(group) %>% summarise(Max = max(PC1))
df2 <- df %>% group_by(group) %>% summarise(Max = max(PC2))
df1$Max <- df1$Max + max(df1$Max)*0.1
df2$Max <- df2$Max + max(df2$Max)*0.1

tuk1 <- aov(PC1~group,data = pcoadata) %>% 
  glht(alternative = 'two.sided',linfct = mcp(group = "Tukey")) %>% cld(alpah = 0.05)
#mod<-aov(PC1~group,data = pcoadata)
#glt<-glht(mod, alternative = 'two.sided', linfct = mcp(group = 'Tukey'))
#summary(glt)
tuk2 <- aov(PC2~group,data = pcoadata) %>% 
  glht(alternative = 'two.sided',linfct = mcp(group = "Tukey")) %>% cld(alpah = 0.05)

res <- data.frame(PC1 = tuk1$mcletters$Letters,PC2 = tuk2$mcletters$Letters,
                  df1 = df1$Max,df2 = df2$Max,group = df1$group)

enterotype_all <- ggplot(pcoadata, aes(PC1, PC2,colour = group, shape = group_1)) +
  geom_point(aes(colour = group),size = 2)+
  labs(x = (floor(pcoa$values$Relative_eig[1]*100)) %>% 
         paste0("PCoA1 ( ", ., "%", " )"),
       y = (floor(pcoa$values$Relative_eig[2]*100)) %>% 
         paste0("PCoA2 ( ", ., "%", " )"), colour = "Enterotype", shape = "Time") +
  scale_colour_manual(values = c("#1B9E77", "#D95F02", "#006BA3")) +
  scale_shape_manual(values = c(3, 18, 19)) +
  theme(text = element_text(size = 15))+
  geom_vline(aes(xintercept = 0),linetype = "dotted")+
  geom_hline(aes(yintercept = 0),linetype = "dotted")+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.title.x = element_text(colour = 'black', size = 10),
        axis.title.y = element_text(colour = 'black', size = 10),
        axis.text = element_text(colour = 'black',size = 10),
        legend.title = element_text(size = 9, hjust = 0.5),
        legend.margin = margin(t = 0, b = 0, r = 0, l = -7),
        legend.text = element_text(size = 8),
        legend.key.height = unit(0.4,"cm"),
        legend.key = element_blank(),
        panel.grid = element_blank())

p <- p2 + plot_spacer() + p1 + p3 + 
  plot_layout(heights = c(1, 4), widths = c(4, 1), ncol = 2, nrow = 2)
ggsave(plot = enterotype_all,"enterotype_all.pdf",dpi = 600,width = 6,height = 4)


##
base_filtercomp <- subset_samples(psstudy5_filtercomp, Time == "Baseline")

##Glom at the genus level
psstudy5filtercompG <- tax_glom(psstudy5_filtercomp, taxrank = "Genus")
sum(taxa_sums(base_filtercomp) == 0) ##no zero sum

##Kurskal wallis test for the difference between Enterotypes
enterokruskal_all <- data.frame(t(otu_table(psstudy5filtercompG)))
enterokruskal_all$Enterotype <- phyloseq::sample_data(psstudy5filtercompG)$Enterotype
enterokruskal_all$Enterotype <- factor(enterokruskal_all$Enterotype, c("Enterotype one", "Enterotype two", "Enterotype three"))

##Adding taxonomic labels
all_taxa <- data.frame(tax_table(psstudy5filtercompG))
all_taxa <- all_taxa %>% select(-Species) %>% rownames_to_column(var = "ASV")

##Relative abundance per ASV##RelatSpeciesive abundance per ASV
rel_kruskalall <- enterokruskal_all %>% 
  gather(key = ASV, value = relabun, -Enterotype) %>% 
  group_by(ASV)

mean_kruskalall <- rel_kruskalall %>% 
  dplyr::group_by(Enterotype, ASV) %>% 
  summarise(mean_relabun = mean(relabun))
write.csv(mean_kruskalall, "output_results/mean_kruskalall.csv")

##add taxa_info to the mean_kruskalbase
all_meankruskal <- full_join(mean_kruskalall, all_taxa)
write.csv(all_meankruskal, "output_results/all_meankruskal.csv")

write.csv(rel_kruskalbase, "output_results/rel_kruskalbase.csv")

##Kruskal wallis test of Enterotype
kruskalentero_basebh <- rel_kruskalall %>% 
  group_by(ASV) %>% 
  kruskal_test(relabun~Enterotype) %>% 
  adjust_pvalue(method = "BH") %>% 
  add_significance()

##Kruskal wallis without p value adjustment
kruskalentero_all <- rel_kruskalall %>% 
  group_by(ASV) %>% 
  kruskal_test(relabun~Enterotype) %>% 
  add_significance()

##filter the significant results and use dunn posthoc
sig_rel_kruskalall <- kruskalentero_all %>% filter(p.signif != "ns")

##Subset for only significant from the original
# Filter rel_kruskalbase to only include ASVs in sig_kruskalenterobase
sig_relkruskalall <- rel_kruskalall %>%
  semi_join(sig_rel_kruskalall, by = "ASV")

##Posthoc dunn test
posthoc_enteroall <- sig_relkruskalall %>% 
  group_by(ASV) %>% 
  dunn_test(relabun ~ Enterotype, p.adjust.method = "BH")%>% 
  add_significance()
write.csv(posthoc_enteroall, "output_results/posthoc_enteroall.csv")

##filter the significant results and use dunn posthoc
sig_posthoc_enteroall <- posthoc_enteroall %>% filter(p.adj.signif != "ns")
write.csv(sig_posthoc_enteroall, "output_results/sig_posthoc_enteroall.csv")
sig_posthoc_enteroall
all_taxa

sigposthocenterotaxa <- right_join(sig_posthoc_enteroall, all_taxa, by = "ASV", relationship = "many-to-many")
write.csv(sigposthocenterotaxa, "output_results/sigposthocenterotaxa.csv")