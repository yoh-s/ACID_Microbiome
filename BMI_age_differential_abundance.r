#Load library

library(phyloseq)
library(microbiome)
library(dplyr)
library(tidyr)
library(rstatix)
library(pheatmap)

# Agglomerate to Genus level and remove zero-sum taxa
dat_ancom_prepared <- tax_glom(ps_study5_filter_2, taxrank = "Genus")
dat_ancom_prepared <- prune_taxa(taxa_sums(dat_ancom_prepared) > 0, dat_ancom_prepared)

# Subset Baseline samples
dat_baseline <- subset_samples(dat_ancom_prepared, Time == "Baseline")
dat_baseline <- prune_taxa(taxa_sums(dat_baseline) > 0, dat_baseline)

# Transform to compositional data
dat_baseline_comp <- microbiome::transform(dat_baseline, "compositional")

# Extract taxonomy table without Species column
taxa_info <- tax_table(dat_baseline_comp) %>%
  as.data.frame() %>%
  dplyr::select(-Species) %>%
  tibble::rownames_to_column(var = "ASV")

# =================== BMI Analysis ===================
BMI_df <- data.frame(t(otu_table(dat_baseline_comp)))
BMI_df$BMI_3_class <- factor(sample_data(dat_baseline_comp)$BMI_3_class,
                             levels = c("Normal", "Overweight", "Obesity"))

BMI_long <- BMI_df %>%
  pivot_longer(cols = -BMI_3_class, names_to = "ASV", values_to = "relabun") %>%
  mutate(percent = relabun * 100)

# Mean relative abundance
mean_BMI <- BMI_long %>%
  group_by(BMI_3_class, ASV) %>%
  summarise(mean_relabun = mean(relabun) * 100, .groups = "drop") %>%
  pivot_wider(names_from = BMI_3_class, values_from = mean_relabun)

mean_BMIwide <- left_join(taxa_info, mean_BMI, by = "ASV") %>%
  dplyr::select(-Class, -Order, -Family)

dir.create("BMIcategory_final", showWarnings = FALSE)
write.csv(mean_BMIwide, "BMIcategory_final/mean_BMIwide_all.csv", row.names = FALSE)

# Wilcoxon test per ASV across BMI groups
wilcox_BMI <- BMI_long %>%
  group_by(ASV) %>%
  wilcox_test(percent ~ BMI_3_class, p.adjust.method = "BH") %>%
  add_significance() %>%
  ungroup()

wilcox_BMI_sig <- wilcox_BMI %>% filter(p.adj.signif != "ns")

# Merge results with taxonomy and mean abundance
final_bmi_wilcoxon <- wilcox_BMI_sig %>%
  left_join(taxa_info, by = "ASV") %>%
  left_join(mean_BMIwide %>% dplyr::select(-Kingdom, -Phylum, -Genus), by = "ASV") %>%
  dplyr::select(-Kingdom, -Class, -Order, -Family)

write.csv(final_bmi_wilcoxon, "BMIcategory_final/final_bmi_wilcoxon.csv", row.names = FALSE)




# =================== Age Quartile Analysis ===================
age_df <- data.frame(t(otu_table(dat_baseline_comp)))
age_df$Age_quartile <- factor(sample_data(dat_baseline_comp)$Age_quartile,
                              levels = c("60-64_years", "65-68_years", "68-73_years", "73-80_years"))

age_long <- age_df %>%
  pivot_longer(cols = -Age_quartile, names_to = "ASV", values_to = "relabun") %>%
  mutate(percent = relabun * 100)

# Mean relative abundance per ASV and age quartile
mean_age <- age_long %>%
  group_by(Age_quartile, ASV) %>%
  summarise(mean_relabun = mean(relabun) * 100, .groups = "drop") %>%
  pivot_wider(names_from = Age_quartile, values_from = mean_relabun)

mean_agewide <- left_join(taxa_info, mean_age, by = "ASV") %>%
  dplyr::select(-Class, -Order, -Family)

dir.create("Agequartile_final", showWarnings = FALSE)
write.csv(mean_agewide, "Agequartile_final/mean_Agequartilewide_all.csv", row.names = FALSE)

# Wilcoxon test per ASV across age quartiles
wilcox_age <- age_long %>%
  filter(!is.na(percent) & !is.na(Age_quartile)) %>%
  group_by(ASV) %>%
  wilcox_test(percent ~ Age_quartile, p.adjust.method = "BH") %>%
  add_significance() %>%
  ungroup()

wilcox_age_sig <- wilcox_age %>% filter(p.adj.signif != "ns")

final_age_wilcoxon <- wilcox_age %>%
  left_join(taxa_info, by = "ASV") %>%
  left_join(mean_agewide %>% dplyr::select(-Kingdom, -Phylum, -Genus), by = "ASV") %>%
  dplyr::select(-Kingdom, -Class, -Order, -Family)

write.csv(final_age_wilcoxon, "Agequartile_final/final_agequartile_wilcoxon.csv", row.names = FALSE)

#Heatmap for Significant Age ASVs
sig_asvs <- wilcox_age_sig$ASV %>% unique()

if(length(sig_asvs) > 0){
  heat_data <- mean_agewide %>%
    filter(ASV %in% sig_asvs) %>%
    tibble::column_to_rownames("Genus") %>%
    dplyr::select(`60-64_years`, `65-68_years`, `68-73_years`, `73-80_years`)

  pheatmap(heat_data, cluster_rows = FALSE, cluster_cols = FALSE,
           scale = "row",
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
}