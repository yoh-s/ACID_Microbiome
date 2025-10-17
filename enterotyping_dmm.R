###############################################################################
# ENTEROTYPING PIPELINE — Two Enterotypes (BIC-based, Wilcoxon Test)
# Author: [Your Name]
# Date: 2025-10-15
###############################################################################

# ─────────────────────────────
# 1. Load libraries
# ─────────────────────────────
library(phyloseq)
library(DirichletMultinomial)
library(microViz)
library(microbiome)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rstatix)
library(scales)
library(reshape2)
library(multcomp)
library(patchwork)
library(readr)
library(vegan)
library(pairwiseAdonis)
library(ape)

#Prepare phyloseq data
ps <- ps_study810_rare
ps <- ps_join(ps, metadatapaired, by = "sampleid")

sample_data(ps) <- data.frame(sample_data(ps)) %>%
  select(-c(
    "Enterotype_BIC...3.x", "Enterotype_BIC...4.x",
    "Baseline_enterotype_BIC.y", "Baseline_enterotype_BIC.x",
    "Enterotype_BIC...3.y", "Enterotype_BIC...4.y"
  ))

#Dirichlet Multinomial Model (Genus level)
ps_genus <- tax_glom(ps, taxrank = "Genus")

count_mat <- as(otu_table(ps_genus), "matrix")
count_mat <- count_mat + 1e-9

n_comp <- 6
dmn_list <- lapply(1:n_comp, function(k) dmn(count_mat, k, verbose = FALSE))
names(dmn_list) <- paste0("dmn_", 1:n_comp)

lplc <- sapply(dmn_list, laplace)
bic  <- sapply(dmn_list, BIC)
aic  <- sapply(dmn_list, AIC)

best_bic <- dmn_list[[which.min(bic)]]
best_k <- which.min(bic)
cat("Best number of enterotypes based on BIC:", best_k, "\n")

pdf("output_results/dirichlet_model_fit.pdf", onefile = TRUE)
plot(lplc, type = "b", main = "Laplace", xlab = "No. of Components")
plot(bic,  type = "b", main = "BIC", xlab = "No. of Components")
plot(aic,  type = "b", main = "AIC", xlab = "No. of Components")
dev.off()

#Assign Enterotypes (Two groups)
best_fit <- best_bic
mix <- mixture(best_fit, assign = TRUE)
enterotype_df <- data.frame(
  SampleID = rownames(mix),
  Enterotype = apply(mix, 1, which.max)
)

enterotype_df <- enterotype_df %>%
  filter(Enterotype %in% c(1, 2)) %>%
  mutate(Enterotype = factor(Enterotype, labels = c("Enterotype one", "Enterotype two")))

write.csv(enterotype_df, "output_results/enterotype_assignments.csv", row.names = FALSE)

#Add Enterotype labels to phyloseq object
ps <- ps_join(ps, enterotype_df, by = "SampleID")

#Enterotype distribution over time
meta <- data.frame(sample_data(ps))
meta$Time <- factor(meta$Time, levels = c("Baseline", "12-wks", "24-wks"))

enterotype_prop <- meta %>%
  group_by(Time, Enterotype) %>%
  summarise(Count = n(), .groups = "drop_last") %>%
  mutate(Proportion = Count / sum(Count))

ggplot(enterotype_prop, aes(x = Time, y = Proportion, fill = Enterotype)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = percent(Proportion, accuracy = 1)),
            position = position_fill(vjust = 0.5)) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02")) +
  theme_bw() +
  labs(x = "", y = "Proportion") +
  theme(legend.title = element_blank())

ggsave("output_results/enterotype_time_distribution.tiff", dpi = 600)

#Wilcoxon Test for Differential Abundance (Baseline)
ps_baseline <- subset_samples(ps, Time == "Baseline")
ps_baseline_genus <- tax_glom(ps_baseline, "Genus")

otu_rel <- microbiome::transform(ps_baseline_genus, "compositional")
df_rel <- psmelt(otu_rel)

# Wilcoxon rank-sum test per OTU (two groups only)
wilcox_results <- df_rel %>%
  group_by(OTU) %>%
  wilcox_test(Abundance ~ Enterotype) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

write.csv(wilcox_results, "output_results/wilcoxon_enterotype_baseline.csv", row.names = FALSE)

#Beta diversity (Bray & UniFrac) and PERMANOVA
ps_genus_all <- tax_glom(ps, "Genus")
ps_genus_comp <- microbiome::transform(ps_genus_all, "compositional")
meta_all <- data.frame(sample_data(ps_genus_comp))

dist_bray <- phyloseq::distance(ps_genus_comp, method = "bray")
set.seed(1988)
perm_bray <- adonis2(dist_bray ~ Enterotype, data = meta_all, permutations = 999, strata = meta_all$Time)

dist_wuni <- phyloseq::distance(ps_genus_all, method = "wunifrac")
set.seed(1988)
perm_wuni <- adonis2(dist_wuni ~ Enterotype, data = meta_all, permutations = 999, strata = meta_all$Time)

#PCoA plot (Bray–Curtis)
pcoa_res <- ape::pcoa(dist_bray)
scores <- data.frame(pcoa_res$vectors[, 1:2])
scores <- cbind(scores, meta_all)

ggplot(scores, aes(Axis.1, Axis.2, color = Enterotype, shape = Time)) +
  geom_point(size = 2) +
  theme_classic() +
  scale_color_manual(values = c("#1B9E77", "#D95F02")) +
  labs(
    x = sprintf("PCoA1 (%.1f%%)", pcoa_res$values$Relative_eig[1] * 100),
    y = sprintf("PCoA2 (%.1f%%)", pcoa_res$values$Relative_eig[2] * 100)
  )

ggsave("output_results/pcoa_enterotypes.pdf", dpi = 600, width = 6, height = 4)