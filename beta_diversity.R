# Load libraries
library(phyloseq)
library(vegan)
library(microbiome)
library(dplyr)
library(ggplot2)
library(multcomp)
library(multcompView)
library(patchwork)
library(ape)
library(pairwiseAdonis)
library(readr)

set.seed(562)

# Prepare data (Genus level)
ps_genus_rare <- tax_glom(ps_study810_rare, taxrank = "Genus")
df_ps_genus_rare <- data.frame(sample_data(ps_genus_rare))

factor_vars <- c(
  "APOE_status","BMI_3_class","Screening_Diagnosis","BMI_SD","Gender",
  "Age_quartile","Family_dementia","Smoker","Coronary_heart_disease",
  "Hypertension_treatment","Diabetes_mellitus_all","Antithrombotic_agents",
  "Calcium_channel_blockers","Thyroid_therapy","Drug_acid_related_disorders",
  "Hypercholesterolemia_drug","Residence"
)

df_ps_genus_rare <- df_ps_genus_rare %>%
  mutate(across(all_of(factor_vars), as.factor))

# Subset baseline samples
base_genus_rare <- subset_samples(ps_genus_rare, Time == "Baseline")
base_genus_rare <- prune_taxa(taxa_sums(base_genus_rare) > 0, base_genus_rare)
base_ps_genus_rare <- df_ps_genus_rare %>% filter(Time == "Baseline")

# ASV-level subset for UniFrac and Bray
base_asv <- subset_samples(ps_study8_rare, Time == "Baseline")
base_asv <- prune_taxa(taxa_sums(base_asv) > 0, base_asv)

get_distance <- function(ps_obj, method) phyloseq::distance(ps_obj, method = method)

run_permanova <- function(dist_mat, metadata, vars, seed = 1988) {
  set.seed(seed)
  results <- lapply(vars, function(var) {
    model <- as.formula(paste("dist_mat ~", var))
    res <- adonis2(model, data = metadata, permutations = 999, na.action = na.exclude)
    data.frame(Variable = var, R2 = res$R2[1], p_value = res$`Pr(>F)`[1], stringsAsFactors = FALSE)
  })
  do.call(rbind, results)
}

# Run PERMANOVA
meta_vars <- c(
  "Screening_Diagnosis","Residence","Age_quartile","Gender","BMI_3_class",
  "Coronary_heart_disease","Hypercholesterolemia_drug","Thyroid_therapy",
  "Family_dementia","Hypertension_treatment","Diabetes_mellitus_all",
  "Antithrombotic_agents","Drug_acid_related_disorders","BMI_SD",
  "Calcium_channel_blockers"
)

dist_methods <- c("bray", "jaccard", "wunifrac", "uunifrac")
results_list <- list()

for (method in dist_methods) {
  message("Running ", method, " ...")
  dist_asv <- get_distance(base_asv, method)
  res_asv <- run_permanova(dist_asv, base_ps_genus_rare, meta_vars)
  res_asv$Level <- "ASV"
  res_asv$Metric <- method
  results_list[[paste0("asv_", method)]] <- res_asv

  if (method %in% c("bray", "jaccard")) {
    dist_genus <- get_distance(microbiome::transform(base_genus_rare, "compositional"), method)
    res_genus <- run_permanova(dist_genus, base_ps_genus_rare, meta_vars)
    res_genus$Level <- "Genus"
    res_genus$Metric <- method
    results_list[[paste0("genus_", method)]] <- res_genus
  }
}

permanova_results <- do.call(rbind, results_list)
write.csv(permanova_results, "permanova_summary.csv", row.names = FALSE)

#Bray–Curtis + Weighted UniFrac (ASV) plot
plot_pcoa_permanova <- function(distance_matrix, metadata, group_var, metric, output_path) {
  library(multcomp)
  library(multcompView)
  library(ggplot2)
  library(patchwork)
  library(ape)
  
  pcoa_res <- ape::pcoa(distance_matrix)
  groups <- metadata
  groups <- groups[match(rownames(pcoa_res$vectors), rownames(groups)), ]
  
  df_pcoa <- data.frame(
    sample = rownames(pcoa_res$vectors),
    PC1 = pcoa_res$vectors[, 1],
    PC2 = pcoa_res$vectors[, 2],
    group = groups[[group_var]]
  )
  
  df_pcoa$group <- factor(df_pcoa$group, 
                          levels = c("60-64_years", "65-68_years", "68-73_years", "73-80_years"))
  
  # PERMANOVA
  set.seed(1988)
  model <- as.formula(paste("distance_matrix ~ group"))
  perm_res <- adonis2(model, data = df_pcoa, permutations = 9999)
  write.table(perm_res, paste0(output_path, "_adonis.tsv"), sep = "\t", quote = FALSE)
  
  # Tukey tests
  tuk1 <- aov(PC1 ~ group, data = df_pcoa) %>% 
    glht(linfct = mcp(group = "Tukey")) %>% 
    cld(alpha = 0.05)
  tuk2 <- aov(PC2 ~ group, data = df_pcoa) %>% 
    glht(linfct = mcp(group = "Tukey")) %>% 
    cld(alpha = 0.05)
  
  df1 <- df_pcoa %>% group_by(group) %>% summarise(Max = max(PC1)) %>% 
    mutate(Max = Max + max(Max) * 0.1)
  df2 <- df_pcoa %>% group_by(group) %>% summarise(Max = max(PC2)) %>% 
    mutate(Max = Max + max(Max) * 0.1)
  
  res <- data.frame(
    PC1 = tuk1$mcletters$Letters,
    PC2 = tuk2$mcletters$Letters,
    df1 = df1$Max,
    df2 = df2$Max,
    group = df1$group
  )
  
  # Plot assembly
  df_pcoa$group <- factor(gsub("_", " ", df_pcoa$group))
  
  p1 <- ggplot(df_pcoa, aes(PC1, PC2, color = group)) +
    geom_point(size = 2) +
    labs(
      x = paste0("PCoA1 (", floor(pcoa_res$values$Relative_eig[1] * 100), "%)"),
      y = paste0("PCoA2 (", floor(pcoa_res$values$Relative_eig[2] * 100), "%)")
    ) +
    scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD")) +
    theme_minimal(base_size = 14) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    theme(legend.position = "none")
  
  p2 <- ggplot(df_pcoa, aes(group, PC1, fill = group)) +
    geom_boxplot(outlier.colour = NA) +
    coord_flip() +
    scale_fill_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD")) +
    theme_void()
  
  p3 <- ggplot(df_pcoa, aes(group, PC2, fill = group)) +
    geom_boxplot(outlier.colour = NA) +
    scale_fill_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD")) +
    theme_void()
  
  p4 <- ggplot() +
    geom_text(
      aes(x = 0, y = 0,
          label = paste0(
            "PERMANOVA:\nR² = ", round(perm_res$R2[1], 4),
            "\np = ", signif(perm_res$`Pr(>F)`[1], 3)
          )),
      size = 4
    ) + theme_void()
  
  final_plot <- p2 + plot_spacer() + p1 + p3 + 
    plot_layout(heights = c(1, 4), widths = c(4, 1), ncol = 2, nrow = 2)
  
  ggsave(plot = final_plot, filename = paste0(output_path, ".tiff"), dpi = 600, width = 6, height = 5)
}

#Run visualizations
base_asv_bray <- get_distance(base_asv, "bray")
base_asv_wunifrac <- get_distance(base_asv, "wunifrac")

dir.create("area_plot/finaltop10", recursive = TRUE, showWarnings = FALSE)

plot_pcoa_permanova(base_asv_bray, base_ps_genus_rare, "Age_quartile", 
                    "BrayCurtis", "area_plot/finaltop10/bray_curtis_agequartile")

plot_pcoa_permanova(base_asv_wunifrac, base_ps_genus_rare, "Age_quartile", 
                    "WeightedUniFrac", "area_plot/finaltop10/weightedunifrac_agequartile")
