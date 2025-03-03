require(knitr)
require(tidyverse)
require(RColorBrewer)
require(cowplot)
require(reshape2)
require(ggdendro)
require(vegan)
require(ape)
library(compositions)
library(PCAtools)
library(tidyverse)
library(car)
library(lme4)
library(nlme)

metadatapaired <- read_csv("metadata_pairedsample.csv")
metadatapaired <- metadatapaired[,-1]

ps_study810_rare_f2 <- ps_join(ps_study810_rare_f2, metadatapaired, by = "sampleid")

View(sample_data(ps_study810_rare_f2))
##4564

##procrsustes anlaysis
##filter pair baseline and endline paired sample
psrare_pair <- subset_samples(ps_study810_rare_f2, Paired_final == "1")
pspair <- psrare_pair

##Remove zero sum
pspair <- prune_taxa(taxa_sums(pspair) > 0, pspair) ##195 ASVs

##Bray-Curtis dissimilarity
##glom at the genus level
pspairGenus <- tax_glom(pspair, taxrank = "Genus") ##293

##Calculate relative abundance
psfilcomp_pairGenus_comp <- microbiome::transform(pspairGenus, transform = "compositional")

df_psfilcomp_pairGenus_comp <- data.frame(sample_data(psfilcomp_pairGenus_comp))

df_psfilcomp_pairGenus_comp$Baseline_enterotype_BIC <- factor(df_psfilcomp_pairGenus_comp$Baseline_enterotype_BIC,
                                                              levels = c("Enterotype one", "Enterotype two"))

df_psfilcomp_pairGenus_comp$Time <- factor(df_psfilcomp_pairGenus_comp$Time, levels = c("Baseline","24-wks"))

enterotypebray <- phyloseq::distance(psfilcomp_pairGenus_comp, method = "bray")

# PERMANOVA using the adonis function on the Bray-Curtis distance matrix

permanova_bray_terms <- adonis2(enterotypebray ~ Treatment_Type + Baseline_enterotype_BIC + Time + 
                                  Treatment_Type:Baseline_enterotype_BIC + 
                                  Treatment_Type:Time + 
                                  Baseline_enterotype_BIC:Time + 
                                  Treatment_Type:Baseline_enterotype_BIC:Time, 
                                data = df_psfilcomp_pairGenus_comp, 
                                permutations = 999, 
                                strata = df_psfilcomp_pairGenus_comp$subject,
                                by = "terms")

# Print the results
permanova_bray_terms

# Perform PCoA
pcoa_result <- cmdscale(vegdist(enterotypebray, method = "bray"), k = 2)

# Create a data frame for plotting
plot_data <- data.frame(pcoa_result, Treatment_Type = df_psfilcomp_pairGenus_comp$Treatment_Type, 
                        Baseline_enterotype_BIC = df_psfilcomp_pairGenus_comp$Baseline_enterotype_BIC, 
                        Time = df_psfilcomp_pairGenus_comp$Time)
# Plot the PCoA results
library(ggplot2)
ggplot(plot_data, aes(x = X1, y = X2, color = Treatment_Type, shape = Baseline_enterotype_BIC)) +
  geom_point(size = 3) +
  facet_wrap(~ Time) +
  labs(title = "PCoA of Enterotype Bray-Curtis Distances", x = "PCoA1", y = "PCoA2") +
  theme_minimal()

##Baseline_enterotype_BIC
##subset based on baseline enterotype
enterotype_one <- subset_samples(psfilcomp_pairGenus_comp, sample_data(psfilcomp_pairGenus_comp)$Baseline_enterotype_BIC == "Enterotype one")

##Remove zero sum
enterotype_one <- prune_taxa(taxa_sums(enterotype_one) > 0, enterotype_one) ##195 ASVs

##subset treatment A
enterotype_oneBase <- subset_samples(enterotype_one, Time == "Baseline")

##Remove zero sum
enterotype_oneBase <- prune_taxa(taxa_sums(enterotype_oneBase) > 0, enterotype_oneBase) ##264 Genus

##subset treatment A
enterotype_oneBaseINT <- subset_samples(enterotype_oneBase, Treatment_Type == "Intervention")

##Remove zero sum
enterotype_oneBaseINT <- prune_taxa(taxa_sums(enterotype_oneBaseINT) > 0, enterotype_oneBaseINT) ##238 ASVs

##Endline in intervention group
enterotype_oneEnd <- subset_samples(enterotype_one, Time == "24-wks")

##Remove zero sum
enterotype_oneEnd <- prune_taxa(taxa_sums(enterotype_oneEnd) > 0, enterotype_oneEnd) ##195 ASVs

##subset treatment A
enterotype_oneEndINT <- subset_samples(enterotype_oneEnd, Treatment_Type == "Intervention")

##Remove zero sum
enterotype_oneEndINT <- prune_taxa(taxa_sums(enterotype_oneEndINT) > 0, enterotype_oneEndINT) ##195 ASVs

##Compare baseline and endline intervention in Enterotype one group
enterotype_oneBaseINT

##########ENTEROTYPE ONE INTERVENTION
#Bray-Curtis at genus level
enteroOneBaseINT_bray <- phyloseq::distance(enterotype_oneBaseINT, method = "bray")

pcoa_enteroOneBaseINT <- cmdscale(enteroOneBaseINT_bray, k = 2)
head(pcoa_enteroOneBaseINT)

rownames(pcoa_enteroOneBaseINT) <- sub("_0$", "", rownames(pcoa_enteroOneBaseINT))

enterotype_oneEndINT
enterotypeOneEndINT_bray <- phyloseq::distance(enterotype_oneEndINT, method = "bray")

pcoa_enterotypeOneEndINT <- cmdscale(enterotypeOneEndINT_bray, k = 2)
head(pcoa_enterotypeOneEndINT)

rownames(pcoa_enterotypeOneEndINT) <- sub("_24$", "", rownames(pcoa_enterotypeOneEndINT))

align_pcoa <- function(pcoa1, pcoa2) {
  common_samples <- intersect(rownames(pcoa1), rownames(pcoa2))
  list(pcoa1[common_samples, ], pcoa2[common_samples, ])
}

aligned_intervention <- align_pcoa(pcoa_enteroOneBaseINT, pcoa_enterotypeOneEndINT)
# Repeat for Enterotype two and three

# procrustes
pro_intervention <- procrustes(aligned_intervention[[1]], aligned_intervention[[2]], symmetric = TRUE, scale = TRUE)

protest_ENTERO1intervention <- protest(aligned_intervention[[1]], aligned_intervention[[2]])

eigen_interone <- sqrt(pro_intervention$svd$d)
percent_var_interone <- signif(eigen_interone/sum(eigen_interone), 4)*100

OneBaseINT <- data.frame(pro_intervention$X)
OneBaseINT$UserName <- rownames(OneBaseINT)
OneBaseINT$type <- "Time(Bray_curtis)"

##sample data
df_enterotype_oneBaseINT <- data.frame(sample_data(enterotype_oneBaseINT))

Time <- df_enterotype_oneBaseINT[,4]

OneBaseINT_pro=cbind(OneBaseINT,Time)

OneEndINT <- data.frame(pro_intervention$Yrot)
OneEndINT$UserName <- rownames(OneEndINT)
OneEndINT$type <- "Time(Bray_curtis)"

##sample data
df_enterotype_oneEndINT <- data.frame(sample_data(enterotype_oneEndINT))
Time <- df_enterotype_oneEndINT[,4]

OneEndINT_pro=cbind(OneEndINT,Time)

# Open a TIFF device with 600 dpi resolution
tiff("area_plot/finaltop10/pro_intervention_plot.tiff", width = 6, height = 6, units = "in", res = 600)
# Plot the Procrustes analysis without points and arrows
plot(pro_intervention, kind = 1, type = "none", 
     main = "", 
     xlab = "PC1", 
     ylab = "PC2")

# Add blue dots for the points
points(pro_intervention, display = "target", pch = 16, col = "#136F53") ##baseline
points(pro_intervention, display = "rotated", pch = 16, col = "#C13C3E") ##endline

# Add arrows with closed triangle heads
arrows(pro_intervention$X[,1], pro_intervention$X[,2], 
       pro_intervention$Yrot[,1], pro_intervention$Yrot[,2], 
       col = "black", length = 0.1, lwd = 1)

# Add a legend without a box
legend("topright", legend = c("Baseline", "Week-24"), 
       col = c("#136F53", "#C13C3E"), pch = 16, bty = "n")

# Close the PDF device
dev.off()

#############################ENTEROTYPE ONE PLACEBO
enterotype_oneBasePLAC <- subset_samples(enterotype_oneBase, Treatment_Type == "Placebo")

##Remove zero sum
enterotype_oneBasePLAC <- prune_taxa(taxa_sums(enterotype_oneBasePLAC) > 0, enterotype_oneBasePLAC) ##195 ASVs

#Bray-Curtis at genus level
enterotype_oneBasePLAC_bray <- phyloseq::distance(enterotype_oneBasePLAC, method = "bray")

pcoa_enteroOneBasePLAC <- cmdscale(enterotype_oneBasePLAC_bray, k = 2)
head(pcoa_enteroOneBasePLAC)

rownames(pcoa_enteroOneBasePLAC) <- sub("_0$", "", rownames(pcoa_enteroOneBasePLAC))

enterotype_oneEnd
enterotype_oneEndPLAC <- subset_samples(enterotype_oneEnd, Treatment_Type == "Placebo")

##Remove zero sum
enterotype_oneEndPLAC <- prune_taxa(taxa_sums(enterotype_oneEndPLAC) > 0, enterotype_oneEndPLAC) ##195 ASVs

#Bray-Curtis at genus level
enterotype_oneEndPLAC_bray <- phyloseq::distance(enterotype_oneEndPLAC, method = "bray")

pcoa_enteroOneEndPLAC <- cmdscale(enterotype_oneEndPLAC_bray, k = 2)
head(pcoa_enteroOneEndPLAC)

rownames(pcoa_enteroOneEndPLAC) <- sub("_24$", "", rownames(pcoa_enteroOneEndPLAC))

aligned_placebo <- align_pcoa(pcoa_enteroOneBasePLAC, pcoa_enteroOneEndPLAC)


pro_placebo <- procrustes(aligned_placebo[[1]], aligned_placebo[[2]], symmetric = TRUE, scale = TRUE)
# Repeat for Enterotype two and three
print(pro_placebo)

protest_placebo <- protest(aligned_placebo[[1]], aligned_placebo[[2]])
print(protest_placebo)


# Open a TIFF device with 600 dpi resolution
tiff("area_plot/finaltop10/pro_placebo_plot.tiff", width = 6, height = 6, units = "in", res = 600)

# Plot the Procrustes analysis without points and arrows
plot(pro_placebo, kind = 1, type = "none", 
     main = "", 
     xlab = "PC1", 
     ylab = "PC2")

# Add blue dots for the points
points(pro_placebo, display = "target", pch = 16, col = "#136F53") ##baseline
points(pro_placebo, display = "rotated", pch = 16, col = "#C13C3E") ##endline

# Add arrows with closed triangle heads
arrows(pro_placebo$X[,1], pro_placebo$X[,2], 
       pro_placebo$Yrot[,1], pro_placebo$Yrot[,2], 
       col = "black", length = 0.1, lwd = 1.5, code = 2, angle = 30)

# Add a legend
# Add a legend without a box
legend("topright", legend = c("Baseline", "Week-24"), 
       col = c("#136F53", "#C13C3E"), pch = 16, bty = "n")


# Close the PDF device
dev.off()




#####################ENTEROTYPE TWO######################################
enterotype_two <- subset_samples(psfilcomp_pairGenus_comp, sample_data(psfilcomp_pairGenus_comp)$Baseline_enterotype_BIC == "Enterotype two")

##Remove zero sum
enterotype_two <- prune_taxa(taxa_sums(enterotype_two) > 0, enterotype_two) ##195 ASVs

##subset treatment A
enterotype_twoBase <- subset_samples(enterotype_two, Time == "Baseline")

##Remove zero sum
enterotype_twoBase <- prune_taxa(taxa_sums(enterotype_twoBase) > 0, enterotype_twoBase) ##195 ASVs

##subset treatment A
enterotype_twoBaseINT <- subset_samples(enterotype_twoBase, Treatment_Type == "Intervention")

##Remove zero sum
enterotype_twoBaseINT <- prune_taxa(taxa_sums(enterotype_twoBaseINT) > 0, enterotype_twoBaseINT) ##195 ASVs

##Endline one entertype treatment A
enterotype_twoEnd <- subset_samples(enterotype_two, Time == "24-wks")

##Remove zero sum
enterotype_twoEnd <- prune_taxa(taxa_sums(enterotype_twoEnd) > 0, enterotype_twoEnd) ##195 ASVs

##subset treatment A
enterotype_twoEndINT <- subset_samples(enterotype_twoEnd, Treatment_Type == "Intervention")

##Remove zero sum
enterotype_twoEndINT <- prune_taxa(taxa_sums(enterotype_twoEndINT) > 0, enterotype_twoEndINT) ##195 ASVs

##Compare baseline and endline intervention in Enterotype one group
enterotype_twoBaseINT

##########ENTEROTYPE ONE INTERVENTION
#Bray-Curtis at genus level
enterotwoBaseINT_bray <- phyloseq::distance(enterotype_twoBaseINT, method = "bray")

pcoa_enteroTwoBaseINT <- cmdscale(enterotwoBaseINT_bray, k = 2)
head(pcoa_enteroTwoBaseINT)

rownames(pcoa_enteroTwoBaseINT) <- sub("_0$", "", rownames(pcoa_enteroTwoBaseINT))

enterotype_twoEndINT
enterotypeTwoEndINT_bray <- phyloseq::distance(enterotype_twoEndINT, method = "bray")

pcoa_enterotypeTwoEndINT <- cmdscale(enterotypeTwoEndINT_bray, k = 2)
head(pcoa_enterotypeTwoEndINT)

rownames(pcoa_enterotypeTwoEndINT) <- sub("_24$", "", rownames(pcoa_enterotypeTwoEndINT))

align_pcoa <- function(pcoa1, pcoa2) {
  common_samples <- intersect(rownames(pcoa1), rownames(pcoa2))
  list(pcoa1[common_samples, ], pcoa2[common_samples, ])
}

aligned_intervention_two <- align_pcoa(pcoa_enteroTwoBaseINT, pcoa_enterotypeTwoEndINT)
# Repeat for Enterotype two and three

# procrustes
pro_intervention_two <- procrustes(aligned_intervention_two[[1]], aligned_intervention_two[[2]])

protest_ENTERO2intervention <- protest(aligned_intervention_two[[1]], aligned_intervention_two[[2]])

# Open a TIFF device with 600 dpi resolution
tiff("area_plot/finaltop10/pro_Intervention_two_plot.tiff", width = 6, height = 6, units = "in", res = 600)

# Plot the Procrustes analysis without points and arrows
plot(pro_intervention_two, kind = 1, type = "none", 
     main = "", 
     xlab = "PC1", 
     ylab = "PC2")

# Add blue dots for the points
points(pro_intervention_two, display = "target", pch = 16, col = "#136F53") ##baseline
points(pro_intervention_two, display = "rotated", pch = 16, col = "#C13C3E") ##endline

# Add arrows with closed triangle heads
arrows(pro_intervention_two$X[,1], pro_intervention_two$X[,2], 
       pro_intervention_two$Yrot[,1], pro_intervention_two$Yrot[,2], 
       col = "black", length = 0.1, lwd = 1.5)

# Add a legend without a box
legend("topright", legend = c("Baseline", "Week-24"), 
       col = c("#136F53", "#C13C3E"), pch = 16, bty = "n")

# Close the PDF device
dev.off()


#############################ENTEROTYPE two PLACEBO
enterotype_twoBasePLAC <- subset_samples(enterotype_twoBase, Treatment_Type == "Placebo")

##Remove zero sum
enterotype_twoBasePLAC <- prune_taxa(taxa_sums(enterotype_twoBasePLAC) > 0, enterotype_twoBasePLAC) ##195 ASVs

#Bray-Curtis at genus level
enterotype_twoBasePLAC_bray <- phyloseq::distance(enterotype_twoBasePLAC, method = "bray")

pcoa_enteroTwoBasePLAC <- cmdscale(enterotype_twoBasePLAC_bray, k = 2)
head(pcoa_enteroTwoBasePLAC)

rownames(pcoa_enteroTwoBasePLAC) <- sub("_0$", "", rownames(pcoa_enteroTwoBasePLAC))

enterotype_twoEnd
enterotype_twoEndPLAC <- subset_samples(enterotype_twoEnd, Treatment_Type == "Placebo")

##Remove zero sum
enterotype_twoEndPLAC <- prune_taxa(taxa_sums(enterotype_twoEndPLAC) > 0, enterotype_twoEndPLAC) ##195 ASVs

#Bray-Curtis at genus level
enterotype_twoEndPLAC_bray <- phyloseq::distance(enterotype_twoEndPLAC, method = "bray")

pcoa_enteroTwoEndPLAC <- cmdscale(enterotype_twoEndPLAC_bray, k = 2)
head(pcoa_enteroTwoEndPLAC)

rownames(pcoa_enteroTwoEndPLAC) <- sub("_24$", "", rownames(pcoa_enteroTwoEndPLAC))

aligned_placebo_two <- align_pcoa(pcoa_enteroTwoBasePLAC, pcoa_enteroTwoEndPLAC)

pro_placebo_two <- procrustes(aligned_placebo_two[[1]], aligned_placebo_two[[2]])

# Repeat for Enterotype two and three
print(pro_placebo_two)

protest_placebo_two <- protest(aligned_placebo_two[[1]], aligned_placebo_two[[2]])
print(protest_placebo_two)

# Open a TIFF device with 600 dpi resolution
tiff("area_plot/finaltop10/pro_placebo_two_plot.tiff", width = 6, height = 6, units = "in", res = 600)
# Plot the Procrustes analysis without points and arrows
plot(pro_placebo_two, kind = 1, type = "none", 
     main = "", 
     xlab = "PC1", 
     ylab = "PC2")

# Add blue dots for the points
points(pro_placebo_two, display = "target", pch = 16, col = "#136F53") ##baseline
points(pro_placebo_two, display = "rotated", pch = 16, col = "#C13C3E") ##endline

# Add arrows with closed triangle heads
arrows(pro_placebo_two$X[,1], pro_placebo_two$X[,2], 
       pro_placebo_two$Yrot[,1], pro_placebo_two$Yrot[,2], 
       col = "black", length = 0.1, lwd = 1)

# Add a legend without a box
legend("topright", legend = c("Baseline", "Week-24"), 
       col = c("#136F53", "#C13C3E"), pch = 16, bty = "n")

# Close the PDF device
dev.off()



###############################WEIGHTED UNIFRAC###################
##Remove zero sum
pspair <- prune_taxa(taxa_sums(pspair) > 0, pspair) ##195 ASVs


##Baseline_enterotype_BIC
##subset based on baseline enterotype
enterotype_one_asv <- subset_samples(pspair, sample_data(pspair)$Baseline_enterotype_BIC == "Enterotype one")

##Remove zero sum
enterotype_one_asv <- prune_taxa(taxa_sums(enterotype_one_asv) > 0, enterotype_one_asv) ##195 ASVs

##subset treatment A
enterotype_oneASVBase <- subset_samples(enterotype_one_asv, Time == "Baseline")

##Remove zero sum
enterotype_oneASVBase <- prune_taxa(taxa_sums(enterotype_oneASVBase) > 0, enterotype_oneASVBase) ##264 Genus

##subset treatment A
enterotype_oneASVBaseINT <- subset_samples(enterotype_oneASVBase, Treatment_Type == "Intervention")

##Remove zero sum
enterotype_oneASVBaseINT <- prune_taxa(taxa_sums(enterotype_oneASVBaseINT) > 0, enterotype_oneASVBaseINT) ##238 ASVs

##Endline in intervention group
enterotype_oneASVEnd <- subset_samples(enterotype_one_asv, Time == "24-wks")

##Remove zero sum
enterotype_oneASVEnd <- prune_taxa(taxa_sums(enterotype_oneASVEnd) > 0, enterotype_oneASVEnd) ##195 ASVs

##subset treatment A
enterotype_oneASVEndINT <- subset_samples(enterotype_oneASVEnd, Treatment_Type == "Intervention")

##Remove zero sum
enterotype_oneASVEndINT <- prune_taxa(taxa_sums(enterotype_oneASVEndINT) > 0, enterotype_oneASVEndINT) ##195 ASVs

##Compare baseline and endline intervention in Enterotype one group
enterotype_oneASVBaseINT

##########ENTEROTYPE ONE INTERVENTION
#Bray-Curtis at genus level
enteroOneBaseINT_wunifrac <- phyloseq::distance(enterotype_oneASVBaseINT, method = "wunifrac")

pcoa_OneBaseINT_wunif <- cmdscale(enteroOneBaseINT_wunifrac, k = 2)
head(pcoa_OneBaseINT_wunif)

rownames(pcoa_OneBaseINT_wunif) <- sub("_0$", "", rownames(pcoa_OneBaseINT_wunif))

enterotype_oneEndINT
enterotypeOneEndINT_wunifrac <- phyloseq::distance(enterotype_oneEndINT, method = "wunifrac")

pcoa_OneEndINT_wunif <- cmdscale(enterotypeOneEndINT_wunifrac, k = 2)
head(pcoa_OneEndINT_wunif)

rownames(pcoa_OneEndINT_wunif) <- sub("_24$", "", rownames(pcoa_OneEndINT_wunif))

align_pcoa <- function(pcoa1, pcoa2) {
  common_samples <- intersect(rownames(pcoa1), rownames(pcoa2))
  list(pcoa1[common_samples, ], pcoa2[common_samples, ])
}

aligned_intervention_wunif <- align_pcoa(pcoa_OneBaseINT_wunif, pcoa_OneEndINT_wunif)
# Repeat for Enterotype two and three

# procrustes
pro_intervention_wunifrac <- procrustes(aligned_intervention_wunif[[1]], aligned_intervention_wunif[[2]])

protest_ENTERO1interv_wunifrac <- protest(aligned_intervention_wunif[[1]], aligned_intervention_wunif[[2]])

####Plot
plot(pro_intervention, kind = 1, main = "Procrustes Analysis",
     xlab = "PC1", ylab = "PC2", ar.col = "blue", pch = 16)

# Add blue dots for the points
points(pro_intervention, display = "target", pch = 16, col = "darkgreen")
points(pro_intervention, display = "rotated", pch = 16, col = "#C13C3E")

# Open a PDF device
pdf("enterotype_final/ProIntervention_One.pdf", width = 6, height = 5)  # Specify width and height as needed
# Plot the Procrustes analysis without points and arrows
plot(pro_intervention, kind = 1, type = "none", 
     main = "", 
     xlab = "PC1", 
     ylab = "PC2")

# Add blue dots for the points
points(pro_intervention, display = "target", pch = 16, col = "#136F53") ##baseline
points(pro_intervention, display = "rotated", pch = 16, col = "#C13C3E") ##endline

# Add arrows with closed triangle heads
arrows(pro_intervention$X[,1], pro_intervention$X[,2], 
       pro_intervention$Yrot[,1], pro_intervention$Yrot[,2], 
       col = "black", length = 0.1, lwd = 1)

# Close the PDF device
dev.off()

#############################ENTEROTYPE ONE PLACEBO
enterotype_oneBasePLAC <- subset_samples(enterotype_oneBase, Treatment_Type == "Placebo")

##Remove zero sum
enterotype_oneBasePLAC <- prune_taxa(taxa_sums(enterotype_oneBasePLAC) > 0, enterotype_oneBasePLAC) ##195 ASVs

#Bray-Curtis at genus level
enterotype_oneBasePLAC_wunifrac <- phyloseq::distance(enterotype_oneBasePLAC, method = "wunifrac")

pcoa_enteroOneBasePLAC_wunifrac <- cmdscale(enterotype_oneBasePLAC_wunifrac, k = 2)
head(pcoa_enteroOneBasePLAC_wunifrac)

rownames(pcoa_enteroOneBasePLAC_wunifrac) <- sub("_0$", "", rownames(pcoa_enteroOneBasePLAC_wunifrac))

enterotype_oneEnd
enterotype_oneEndPLAC <- subset_samples(enterotype_oneEnd, Treatment_Type == "Placebo")

##Remove zero sum
enterotype_oneEndPLAC <- prune_taxa(taxa_sums(enterotype_oneEndPLAC) > 0, enterotype_oneEndPLAC) ##195 ASVs

#Bray-Curtis at genus level
enterotype_oneEndPLAC_wunifrac <- phyloseq::distance(enterotype_oneEndPLAC, method = "wunifrac")

pcoa_enteroOneEndPLAC_wunifrac <- cmdscale(enterotype_oneEndPLAC_wunifrac, k = 2)
head(pcoa_enteroOneEndPLAC_wunifrac)

rownames(pcoa_enteroOneEndPLAC_wunifrac) <- sub("_24$", "", rownames(pcoa_enteroOneEndPLAC_wunifrac))

aligned_placebo_wunifrac <- align_pcoa(pcoa_enteroOneBasePLAC_wunifrac, pcoa_enteroOneEndPLAC_wunifrac)


pro_placebo_wunifrac <- procrustes(aligned_placebo_wunifrac[[1]], aligned_placebo_wunifrac[[2]])
# Repeat for Enterotype two and three
print(pro_placebo_wunifrac)

plot(pro_placebo, kind = 1)


protest_placebo_wunifrac <- protest(aligned_placebo_wunifrac[[1]], aligned_placebo_wunifrac[[2]])
print(protest_placebo_wunifrac)

#####################ENTEROTYPE TWO######################################
enterotype_two <- subset_samples(psfilcomp_pairGenus_comp, sample_data(psfilcomp_pairGenus_comp)$Baseline_enterotype_BIC == "Enterotype two")

##Remove zero sum
enterotype_two <- prune_taxa(taxa_sums(enterotype_two) > 0, enterotype_two) ##195 ASVs

##subset treatment A
enterotype_twoBase <- subset_samples(enterotype_two, Time == "Baseline")

##Remove zero sum
enterotype_twoBase <- prune_taxa(taxa_sums(enterotype_twoBase) > 0, enterotype_twoBase) ##195 ASVs

##subset treatment A
enterotype_twoBaseINT <- subset_samples(enterotype_twoBase, Treatment_Type == "Intervention")

##Remove zero sum
enterotype_twoBaseINT <- prune_taxa(taxa_sums(enterotype_twoBaseINT) > 0, enterotype_twoBaseINT) ##195 ASVs

##Endline one entertype treatment A
enterotype_twoEnd <- subset_samples(enterotype_two, Time == "24-wks")

##Remove zero sum
enterotype_twoEnd <- prune_taxa(taxa_sums(enterotype_twoEnd) > 0, enterotype_twoEnd) ##195 ASVs

##subset treatment A
enterotype_twoEndINT <- subset_samples(enterotype_twoEnd, Treatment_Type == "Intervention")

##Remove zero sum
enterotype_twoEndINT <- prune_taxa(taxa_sums(enterotype_twoEndINT) > 0, enterotype_twoEndINT) ##195 ASVs

##Compare baseline and endline intervention in Enterotype one group
enterotype_twoBaseINT

##########ENTEROTYPE ONE INTERVENTION
#Bray-Curtis at genus level
enterotwoBaseINT_wunifrac <- phyloseq::distance(enterotype_twoBaseINT, method = "wunifrac")

pcoa_enteroTwoBaseINT_wunifrac <- cmdscale(enterotwoBaseINT_wunifrac, k = 2)
head(pcoa_enteroTwoBaseINT_wunifrac)

rownames(pcoa_enteroTwoBaseINT_wunifrac) <- sub("_0$", "", rownames(pcoa_enteroTwoBaseINT_wunifrac))

enterotype_twoEndINT
enterotypeTwoEndINT_wunifrac <- phyloseq::distance(enterotype_twoEndINT, method = "wunifrac")

pcoa_enterotypeTwoEndINT_wunifrac <- cmdscale(enterotypeTwoEndINT_wunifrac, k = 2)
head(pcoa_enterotypeTwoEndINT_wunifrac)

rownames(pcoa_enterotypeTwoEndINT_wunifrac) <- sub("_24$", "", rownames(pcoa_enterotypeTwoEndINT_wunifrac))

align_pcoa <- function(pcoa1, pcoa2) {
  common_samples <- intersect(rownames(pcoa1), rownames(pcoa2))
  list(pcoa1[common_samples, ], pcoa2[common_samples, ])
}

aligned_interventiontwo_wunifrac <- align_pcoa(pcoa_enteroTwoBaseINT_wunifrac, pcoa_enterotypeTwoEndINT_wunifrac)
# Repeat for Enterotype two and three

# procrustes
pro_intervention_two_wunifrac <- procrustes(aligned_interventiontwo_wunifrac[[1]], aligned_interventiontwo_wunifrac[[2]])

protest_ENTERO2intervention_wunifrac <- protest(aligned_interventiontwo_wunifrac[[1]], aligned_interventiontwo_wunifrac[[2]])

# Open a PDF device
pdf("enterotype_final/ProIntervention_Two.pdf", width = 6, height = 5)  # Specify width and height as needed
# Plot the Procrustes analysis without points and arrows
plot(pro_intervention_two, kind = 1, type = "none", 
     main = "", 
     xlab = "PC1", 
     ylab = "PC2")

# Add blue dots for the points
points(pro_intervention_two, display = "target", pch = 16, col = "#136F53") ##baseline
points(pro_intervention_two, display = "rotated", pch = 16, col = "#C13C3E") ##endline

# Add arrows with closed triangle heads
arrows(pro_intervention_two$X[,1], pro_intervention_two$X[,2], 
       pro_intervention_two$Yrot[,1], pro_intervention_two$Yrot[,2], 
       col = "black", length = 0.1, lwd = 1.5)

# Close the PDF device
dev.off()


#############################ENTEROTYPE two PLACEBO
enterotype_twoBasePLAC <- subset_samples(enterotype_twoBase, Treatment_Type == "Placebo")

##Remove zero sum
enterotype_twoBasePLAC <- prune_taxa(taxa_sums(enterotype_twoBasePLAC) > 0, enterotype_twoBasePLAC) ##195 ASVs

#Bray-Curtis at genus level
enterotype_twoBasePLAC_wunifrac <- phyloseq::distance(enterotype_twoBasePLAC, method = "wunifrac")

pcoa_enteroTwoBasePLAC_wunifrac <- cmdscale(enterotype_twoBasePLAC_wunifrac, k = 2)
head(pcoa_enteroTwoBasePLAC_wunifrac)

rownames(pcoa_enteroTwoBasePLAC_wunifrac) <- sub("_0$", "", rownames(pcoa_enteroTwoBasePLAC_wunifrac))

enterotype_twoEnd
enterotype_twoEndPLAC <- subset_samples(enterotype_twoEnd, Treatment_Type == "Placebo")

##Remove zero sum
enterotype_twoEndPLAC <- prune_taxa(taxa_sums(enterotype_twoEndPLAC) > 0, enterotype_twoEndPLAC) ##195 ASVs

#Bray-Curtis at genus level
enterotype_twoEndPLAC_wunifrac <- phyloseq::distance(enterotype_twoEndPLAC, method = "wunifrac")

pcoa_enteroTwoEndPLAC_wunifrac <- cmdscale(enterotype_twoEndPLAC_wunifrac, k = 2)
head(pcoa_enteroTwoEndPLAC_wunifrac)

rownames(pcoa_enteroTwoEndPLAC_wunifrac) <- sub("_24$", "", rownames(pcoa_enteroTwoEndPLAC_wunifrac))

aligned_placebo_two_wunifrac <- align_pcoa(pcoa_enteroTwoBasePLAC_wunifrac, pcoa_enteroTwoEndPLAC_wunifrac)

pro_placebo_two_wunifrac <- procrustes(aligned_placebo_two_wunifrac[[1]], aligned_placebo_two_wunifrac[[2]])

# Repeat for Enterotype two and three
print(pro_placebo_two)

protest_placebo_two_wunifrac <- protest(aligned_placebo_two_wunifrac[[1]], aligned_placebo_two_wunifrac[[2]])
print(protest_placebo_two_wunifrac)


#############################JACCARD######################################################
##Remove zero sum
pspair <- prune_taxa(taxa_sums(pspair) > 0, pspair) ##195 ASVs


##Baseline_enterotype_BIC
##subset based on baseline enterotype
enterotype_one_asv <- subset_samples(pspair, sample_data(pspair)$Baseline_enterotype_BIC == "Enterotype one")

##Remove zero sum
enterotype_one_asv <- prune_taxa(taxa_sums(enterotype_one_asv) > 0, enterotype_one_asv) ##195 ASVs

##subset treatment A
enterotype_oneASVBase <- subset_samples(enterotype_one_asv, Time == "Baseline")

##Remove zero sum
enterotype_oneASVBase <- prune_taxa(taxa_sums(enterotype_oneASVBase) > 0, enterotype_oneASVBase) ##264 Genus

##subset treatment A
enterotype_oneASVBaseINT <- subset_samples(enterotype_oneASVBase, Treatment_Type == "Intervention")

##Remove zero sum
enterotype_oneASVBaseINT <- prune_taxa(taxa_sums(enterotype_oneASVBaseINT) > 0, enterotype_oneASVBaseINT) ##238 ASVs

##Endline in intervention group
enterotype_oneASVEnd <- subset_samples(enterotype_one_asv, Time == "24-wks")

##Remove zero sum
enterotype_oneASVEnd <- prune_taxa(taxa_sums(enterotype_oneASVEnd) > 0, enterotype_oneASVEnd) ##195 ASVs

##subset treatment A
enterotype_oneASVEndINT <- subset_samples(enterotype_oneASVEnd, Treatment_Type == "Intervention")

##Remove zero sum
enterotype_oneASVEndINT <- prune_taxa(taxa_sums(enterotype_oneASVEndINT) > 0, enterotype_oneASVEndINT) ##195 ASVs

##Compare baseline and endline intervention in Enterotype one group
enterotype_oneASVBaseINT

##########ENTEROTYPE ONE INTERVENTION
#Bray-Curtis at genus level
enteroOneBaseINT_jaccard <- phyloseq::distance(enterotype_oneASVBaseINT, method = "jaccard", binary = TRUE)

pcoa_OneBaseINT_jaccard <- cmdscale(enteroOneBaseINT_jaccard, k = 2)
head(pcoa_OneBaseINT_jaccard)

rownames(pcoa_OneBaseINT_jaccard) <- sub("_0$", "", rownames(pcoa_OneBaseINT_jaccard))

enterotype_oneEndINT
enterotypeOneEndINT_jaccard <- phyloseq::distance(enterotype_oneEndINT, method = "jaccard", binary = TRUE)

pcoa_OneEndINT_jaccard <- cmdscale(enterotypeOneEndINT_jaccard, k = 2)
head(pcoa_OneEndINT_jaccard)

rownames(pcoa_OneEndINT_jaccard) <- sub("_24$", "", rownames(pcoa_OneEndINT_jaccard))

align_pcoa <- function(pcoa1, pcoa2) {
  common_samples <- intersect(rownames(pcoa1), rownames(pcoa2))
  list(pcoa1[common_samples, ], pcoa2[common_samples, ])
}

aligned_intervention_jaccard <- align_pcoa(pcoa_OneBaseINT_jaccard, pcoa_OneEndINT_jaccard)
# Repeat for Enterotype two and three

# procrustes
pro_intervention_jaccard <- procrustes(aligned_intervention_jaccard[[1]], aligned_intervention_jaccard[[2]])

protest_ENTERO1interv_jaccard <- protest(aligned_intervention_jaccard[[1]], aligned_intervention_jaccard[[2]])

# Open a PDF device
pdf("enterotype_final/ProInterve_JACCARD_One.pdf", width = 6, height = 5)  # Specify width and height as needed
# Plot the Procrustes analysis without points and arrows
plot(pro_intervention_jaccard, kind = 1, type = "none", 
     main = "", 
     xlab = "PC1", 
     ylab = "PC2")

# Add blue dots for the points
points(pro_intervention_jaccard, display = "target", pch = 16, col = "#136F53") ##baseline
points(pro_intervention_jaccard, display = "rotated", pch = 16, col = "#C13C3E") ##endline

# Add arrows with closed triangle heads
arrows(pro_intervention_jaccard$X[,1], pro_intervention_jaccard$X[,2], 
       pro_intervention_jaccard$Yrot[,1], pro_intervention_jaccard$Yrot[,2], 
       col = "black", length = 0.1, lwd = 1)

# Close the PDF device
dev.off()

#############################ENTEROTYPE ONE PLACEBO
enterotype_oneBasePLAC <- subset_samples(enterotype_oneBase, Treatment_Type == "Placebo")

##Remove zero sum
enterotype_oneBasePLAC <- prune_taxa(taxa_sums(enterotype_oneBasePLAC) > 0, enterotype_oneBasePLAC) ##195 ASVs

#Bray-Curtis at genus level
enterotype_oneBasePLAC_jaccard <- phyloseq::distance(enterotype_oneBasePLAC, method = "jaccard",
                                                      binary = TRUE)

pcoa_enteroOneBasePLAC_jaccard <- cmdscale(enterotype_oneBasePLAC_jaccard, k = 2)
head(pcoa_enteroOneBasePLAC_jaccard)

rownames(pcoa_enteroOneBasePLAC_jaccard) <- sub("_0$", "", rownames(pcoa_enteroOneBasePLAC_jaccard))

enterotype_oneEnd
enterotype_oneEndPLAC <- subset_samples(enterotype_oneEnd, Treatment_Type == "Placebo")

##Remove zero sum
enterotype_oneEndPLAC <- prune_taxa(taxa_sums(enterotype_oneEndPLAC) > 0, enterotype_oneEndPLAC) ##195 ASVs

#Bray-Curtis at genus level
enterotype_oneEndPLAC_jaccard <- phyloseq::distance(enterotype_oneEndPLAC, method = "jaccard",
                                                     binary = TRUE)

pcoa_enteroOneEndPLAC_jaccard <- cmdscale(enterotype_oneEndPLAC_jaccard, k = 2)
head(pcoa_enteroOneEndPLAC_jaccard)

rownames(pcoa_enteroOneEndPLAC_jaccard) <- sub("_24$", "", rownames(pcoa_enteroOneEndPLAC_jaccard))

aligned_placebo_jaccard <- align_pcoa(pcoa_enteroOneBasePLAC_jaccard, pcoa_enteroOneEndPLAC_jaccard)


pro_placebo_jaccard <- procrustes(aligned_placebo_jaccard[[1]], aligned_placebo_jaccard[[2]])
# Repeat for Enterotype two and three
print(pro_placebo_wunifrac)

plot(pro_placebo, kind = 1)

protest_placebo_jaccard <- protest(aligned_placebo_jaccard[[1]], aligned_placebo_jaccard[[2]])
print(protest_placebo_jaccard)

# Open a PDF device
pdf("enterotype_final/ProPlacebo_JACCARD_One.pdf", width = 6, height = 5)  # Specify width and height as needed
# Plot the Procrustes analysis without points and arrows
plot(pro_placebo_jaccard, kind = 1, type = "none", 
     main = "", 
     xlab = "PC1", 
     ylab = "PC2")

# Add blue dots for the points
points(pro_placebo_jaccard, display = "target", pch = 16, col = "#136F53") ##baseline
points(pro_placebo_jaccard, display = "rotated", pch = 16, col = "#C13C3E") ##endline

# Add arrows with closed triangle heads
arrows(pro_placebo_jaccard$X[,1], pro_placebo_jaccard$X[,2], 
       pro_placebo_jaccard$Yrot[,1], pro_placebo_jaccard$Yrot[,2], 
       col = "black", length = 0.1, lwd = 1)

# Close the PDF device
dev.off()



#####################ENTEROTYPE TWO######################################
enterotype_two <- subset_samples(psfilcomp_pairGenus_comp, sample_data(psfilcomp_pairGenus_comp)$Baseline_enterotype_BIC == "Enterotype two")

##Remove zero sum
enterotype_two <- prune_taxa(taxa_sums(enterotype_two) > 0, enterotype_two) ##195 ASVs

##subset treatment A
enterotype_twoBase <- subset_samples(enterotype_two, Time == "Baseline")

##Remove zero sum
enterotype_twoBase <- prune_taxa(taxa_sums(enterotype_twoBase) > 0, enterotype_twoBase) ##195 ASVs

##subset treatment A
enterotype_twoBaseINT <- subset_samples(enterotype_twoBase, Treatment_Type == "Intervention")

##Remove zero sum
enterotype_twoBaseINT <- prune_taxa(taxa_sums(enterotype_twoBaseINT) > 0, enterotype_twoBaseINT) ##195 ASVs

##Endline one entertype treatment A
enterotype_twoEnd <- subset_samples(enterotype_two, Time == "24-wks")

##Remove zero sum
enterotype_twoEnd <- prune_taxa(taxa_sums(enterotype_twoEnd) > 0, enterotype_twoEnd) ##195 ASVs

##subset treatment A
enterotype_twoEndINT <- subset_samples(enterotype_twoEnd, Treatment_Type == "Intervention")

##Remove zero sum
enterotype_twoEndINT <- prune_taxa(taxa_sums(enterotype_twoEndINT) > 0, enterotype_twoEndINT) ##195 ASVs

##Compare baseline and endline intervention in Enterotype one group
enterotype_twoBaseINT

# Calculate effect sizes
effect_size_intervention_1 <- procrustes(aligned_intervention[[1]], aligned_intervention[[2]])$ss
effect_size_placebo_1 <- procrustes(aligned_placebo[[1]], aligned_placebo[[2]])$ss
effect_size_intervention_2 <- procrustes(aligned_intervention_two[[1]], aligned_intervention_two[[2]])$ss
effect_size_placebo_2 <- procrustes(aligned_placebo_two[[1]], aligned_placebo_two[[2]])$ss

# Combine data into a dataframe
data_effectsize_braycurtis <- data.frame(
  Group = rep(c("Intervention", "Placebo"), each = 2),
  Enterotype = rep(c("One", "Two"), times = 2),
  EffectSize = c(effect_size_intervention_1, effect_size_placebo_1, effect_size_intervention_2, effect_size_placebo_2)
)

# Perform ANOVA
anova_result <- aov(EffectSize ~ Group * Enterotype, data = data_effectsize_braycurtis)
summary(anova_result)

# Perform ANOVA using car package
anova_result <- aov(EffectSize ~ Group * Enterotype, data = data_effectsize_braycurtis)
anova_summary <- Anova(anova_result, type = "III")  # Type III sums of squares
anova_summary





##########ENTEROTYPE ONE INTERVENTION
#Bray-Curtis at genus level
enterotwoBaseINT_jaccard <- phyloseq::distance(enterotype_twoBaseINT, method = "jaccard",
                                                binary = TRUE)

pcoa_enteroTwoBaseINT_jaccard <- cmdscale(enterotwoBaseINT_jaccard, k = 2)
head(pcoa_enteroTwoBaseINT_jaccard)

rownames(pcoa_enteroTwoBaseINT_jaccard) <- sub("_0$", "", rownames(pcoa_enteroTwoBaseINT_jaccard))

enterotype_twoEndINT
enterotypeTwoEndINT_jaccard <- phyloseq::distance(enterotype_twoEndINT, method = "jaccard",
                                                   binary = TRUE)

pcoa_enterotypeTwoEndINT_jaccard <- cmdscale(enterotypeTwoEndINT_jaccard, k = 2)
head(pcoa_enterotypeTwoEndINT_jaccard)

rownames(pcoa_enterotypeTwoEndINT_jaccard) <- sub("_24$", "", rownames(pcoa_enterotypeTwoEndINT_jaccard))

align_pcoa <- function(pcoa1, pcoa2) {
  common_samples <- intersect(rownames(pcoa1), rownames(pcoa2))
  list(pcoa1[common_samples, ], pcoa2[common_samples, ])
}

aligned_interventiontwo_jaccard <- align_pcoa(pcoa_enteroTwoBaseINT_jaccard, pcoa_enterotypeTwoEndINT_jaccard)
# Repeat for Enterotype two and three

# procrustes
pro_intervention_two_jaccard <- procrustes(aligned_interventiontwo_jaccard[[1]], aligned_interventiontwo_jaccard[[2]])

protest_ENTERO2intervention_jaccard <- protest(aligned_interventiontwo_jaccard[[1]], aligned_interventiontwo_jaccard[[2]])

# Open a PDF device
pdf("enterotype_final/ProInterv_JACCARD_Two.pdf", width = 6, height = 5)  # Specify width and height as needed
# Plot the Procrustes analysis without points and arrows
plot(pro_intervention_two_jaccard, kind = 1, type = "none", 
     main = "", 
     xlab = "PC1", 
     ylab = "PC2",)

# Add blue dots for the points
points(pro_intervention_two_jaccard, display = "target", pch = 16, col = "#136F53") ##baseline
points(pro_intervention_two_jaccard, display = "rotated", pch = 16, col = "#C13C3E") ##endline

# Add arrows with closed triangle heads
arrows(pro_intervention_two_jaccard$X[,1], pro_intervention_two_jaccard$X[,2], 
       pro_intervention_two_jaccard$Yrot[,1], pro_intervention_two_jaccard$Yrot[,2], 
       col = "black", length = 0.1, lwd = 1.5)

# Close the PDF device
dev.off()


#############################ENTEROTYPE two PLACEBO
enterotype_twoBasePLAC <- subset_samples(enterotype_twoBase, Treatment_Type == "Placebo")

##Remove zero sum
enterotype_twoBasePLAC <- prune_taxa(taxa_sums(enterotype_twoBasePLAC) > 0, enterotype_twoBasePLAC) ##195 ASVs

#Bray-Curtis at genus level
enterotype_twoBasePLAC_wunifrac <- phyloseq::distance(enterotype_twoBasePLAC, method = "wunifrac")

pcoa_enteroTwoBasePLAC_wunifrac <- cmdscale(enterotype_twoBasePLAC_wunifrac, k = 2)
head(pcoa_enteroTwoBasePLAC_wunifrac)

rownames(pcoa_enteroTwoBasePLAC_wunifrac) <- sub("_0$", "", rownames(pcoa_enteroTwoBasePLAC_wunifrac))

enterotype_twoEnd
enterotype_twoEndPLAC <- subset_samples(enterotype_twoEnd, Treatment_Type == "Placebo")

##Remove zero sum
enterotype_twoEndPLAC <- prune_taxa(taxa_sums(enterotype_twoEndPLAC) > 0, enterotype_twoEndPLAC) ##195 ASVs

#Bray-Curtis at genus level
enterotype_twoEndPLAC_jaccard <- phyloseq::distance(enterotype_twoEndPLAC, method = "jaccard",
                                                     binary = TRUE)

pcoa_enteroTwoEndPLAC_jaccard <- cmdscale(enterotype_twoEndPLAC_jaccard, k = 2)
head(pcoa_enteroTwoEndPLAC_jaccard)

rownames(pcoa_enteroTwoEndPLAC_jaccard) <- sub("_24$", "", rownames(pcoa_enteroTwoEndPLAC_jaccard))

aligned_placebo_two_jaccard <- align_pcoa(pcoa_enteroTwoBasePLAC_wunifrac, pcoa_enteroTwoEndPLAC_jaccard)

pro_placebo_two_jaccard <- procrustes(aligned_placebo_two_jaccard[[1]], aligned_placebo_two_jaccard[[2]])

# Repeat for Enterotype two and three
print(pro_placebo_two)

protest_placebo_two_jaccard <- protest(aligned_placebo_two_jaccard[[1]], aligned_placebo_two_jaccard[[2]])
print(protest_placebo_two_jaccard)

# Open a PDF device
pdf("enterotype_final/PropLACEBO_JACCARD_Two.pdf", width = 6, height = 5)  # Specify width and height as needed
# Plot the Procrustes analysis without points and arrows
plot(pro_placebo_two_jaccard, kind = 1, type = "none", 
     main = "", 
     xlab = "PC1", 
     ylab = "PC2")

# Add blue dots for the points
points(pro_placebo_two_jaccard, display = "target", pch = 16, col = "#136F53") ##baseline
points(pro_placebo_two_jaccard, display = "rotated", pch = 16, col = "#C13C3E") ##endline

# Add arrows with closed triangle heads
arrows(pro_placebo_two_jaccard$X[,1], pro_placebo_two_jaccard$X[,2], 
       pro_placebo_two_jaccard$Yrot[,1], pro_placebo_two_jaccard$Yrot[,2], 
       col = "black", length = 0.1, lwd = 1.5)

# Close the PDF device
dev.off()

# Load required libraries
library(phyloseq)
library(vegan)

##Calculate relative abundance
psfilcomp_pairGenus_comp <- microbiome::transform(pspairGenus, transform = "compositional")

df_psfilcomp_pairGenus_comp <- data.frame(sample_data(psfilcomp_pairGenus_comp))

df_psfilcomp_pairGenus_comp$Baseline_enterotype_BIC <- factor(df_psfilcomp_pairGenus_comp$Baseline_enterotype_BIC,
                                                              levels = c("Enterotype one", "Enterotype two"))

df_psfilcomp_pairGenus_comp$Time <- factor(df_psfilcomp_pairGenus_comp$Time, levels = c("Baseline","24-wks"))

genus_composition_bray <- phyloseq::distance(psfilcomp_pairGenus_comp, method = "bray")

# Run db-RDA with Group, Time, and Enterotype, including interaction terms
# Adjust the formula as needed for interaction terms
db_rda_baseentero <- capscale(genus_composition_bray ~ Treatment_Type*Time*Baseline_enterotype_BIC, data = df_psfilcomp_pairGenus_comp)

# Check the results
summary(db_rda_baseentero)

# Test for significance using permutation test
anova(db_rda_baseentero, permutations = 999)

# Perform post-hoc test for specific contrasts of db-RDA axes
permute_rda <- vegan::anova.cca(db_rda_baseentero, permutations = 999)
print(permute_rda)

# Visualize ordination (CAP or CCA)
plot(db_rda_baseentero, display = c("sites", "species"), main = "CAP/CCA Plot")

# If you want to color by Time (Baseline vs. 24-wks)
ordiplot(db_rda_baseentero, type = "t", display = "sites")
points(db_rda_baseentero, col = df_psfilcomp_pairGenus_comp$Time, pch = 16)

# Add legends or customizations as needed
legend("topright", legend = levels(df_psfilcomp_pairGenus_comp$Time), col = 1:2, pch = 16)

otu_pairGenuscomp <- as.matrix(t(otu_table(psfilcomp_pairGenus_comp)))

# Run distance-based RDA
dbrda_entero <- dbrda(otu_pairGenuscomp ~ Baseline_enterotype_BIC * Treatment_Type * Time + Condition(subject), data = df_psfilcomp_pairGenus_comp, distance = "bray")

# Summary of the dbRDA result
summary(dbrda_entero)

# Load necessary libraries
library(vegan)
library(car)

# Check for multicollinearity
vif_model <- lm(otu_pairGenuscomp ~ Baseline_enterotype_BIC * Treatment_Type + Time, data = df_psfilcomp_pairGenus_comp)

# Check for multicollinearity
vif_values <- vif(vif_model, type = "predictor")
print(vif_values)

# Simplify the model by removing collinear terms
dbrda_entero <- dbrda(otu_pairGenuscomp ~ Baseline_enterotype_BIC * Treatment_Type * Time + Condition(subject), 
                      data = df_psfilcomp_pairGenus_comp, distance = "bray")
# Summary of the dbRDA result
summary(dbrda_entero)

# Plot the dbRDA results
plot(dbrda_entero, display = "sites", type = "n")
points(dbrda_entero, display = "sites", col = as.factor(df_psfilcomp_pairGenus_comp$Baseline_enterotype_BIC))
text(dbrda_entero, display = "bp", col = "red")

# Add legend
legend("topright", legend = levels(as.factor(df_psfilcomp_pairGenus_comp$Baseline_enterotype_BIC)), col = 1:length(levels(as.factor(df_psfilcomp_pairGenus_comp$Baseline_enterotype_BIC))), pch = 1)


# Partition variance
varpart_result <- varpart(otu_pairGenuscomp, ~ Baseline_enterotype_BIC, ~ Treatment_Type, ~ subject, data = df_psfilcomp_pairGenus_comp)

# Plot the results
plot(varpart_result)

adonis2(genus_composition_bray ~ Treatment_Type * Time * Baseline_enterotype_BIC, strata = df_psfilcomp_pairGenus_comp$subject, data = df_psfilcomp_pairGenus_comp, permutations = 999)

# Alternatively, you can use the pairwise.adonis function to perform pairwise comparisons directly

# For example, for all Time comparisons within Enterotype one / Intervention
time_entero_treat <- pairwise.adonis2(genus_composition_bray ~  * Baseline_enterotype_BIC * Treatment_Type, 
                                   data = df_psfilcomp_pairGenus_comp,strata = 'subject', perm = 999)

# Print results
print(time_entero_treat)









library(lme4)
model <- lmer(genus_composition_bray ~ Treatment_Type * Time * Baseline_enterotype_BIC + (1|subject), data = distances_data)
summary(model)

