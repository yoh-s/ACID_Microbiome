# --- Load Libraries ---
library(phyloseq)
library(phyloseq.extended) # for Faith PD
library(microbiome)
library(rstatix)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(magrittr)

# --- Compute Alpha Diversity Metrics ---
alpha_all <- microbiome::alpha(ps_study810_rare_f, index = "all") %>%
  rownames_to_column(var = "sampleid") %>%
  select(sampleid, observed, diversity_shannon, evenness_pielou)

faith_phylo <- phylodiv(ps_study810_rare_f) # Faith’s PD

alpha_diversity <- right_join(faith_phylo, alpha_all, by = "sampleid") %>%
  select(-contains("..."), -contains("Time.y")) %>%
  rename(pd = pd) %>%
  relocate(sampleid, pd, observed, diversity_shannon, evenness_pielou)

write.csv(alpha_diversity, "alpha_diversity_f.csv", row.names = FALSE)

# --- Merge with Phyloseq Object ---
ps_study810_rare_f1 <- ps_join(ps_study810_rare_f, alpha_diversity, by = "sampleid")

# --- Merge Classification Metadata ---
alphadiversity_class <- read_csv("alphadiversity_classification.csv")
ps_study810_rare_f2 <- ps_join(ps_study810_rare_f1, alphadiversity_class, by = "sampleid")

# --- Clean & Rename Metadata ---
sample_data(ps_study810_rare_f2) <- data.frame(sample_data(ps_study810_rare_f2)) %>%
  select(-contains("..."), -contains("Time.y"), -Screening_code) %>%
  rename(Time = Time.x)

# --- Prepare Data for Statistical Tests ---
alphadiversity_category <- alpha_diversity %>%
  right_join(sample_data(ps_study810_rare_f2) %>% as.data.frame(), by = "sampleid") %>%
  mutate(across(where(is.character), factor))

# --- Baseline Subset ---
baseline_alpha <- alphadiversity_category %>% filter(Time == "Baseline")

# Variables to test
krus_vars <- c("Screening_Diagnosis", "Residence", "Age_quartile", "Gender",
               "BMI_3_class", "Coronary_heart_disease", "Hypercholesterolemia_drug",
               "Thyroid_therapy", "Family_dementia", "Hypertension_treatment",
               "Calcium_channel_blockers", "Diabetes_mellitus_all",
               "Antithrombotic_agents", "Drug_acid_related_disorders",
               "BMI_SD", "Enterotype", "Baseline_enterotype")

# --- Helper: Kruskal–Wallis Loop ---
kruskal_batch <- function(df, variable_list, metric) {
  map_dfr(variable_list, function(var) {
    formula <- as.formula(paste(metric, "~", var))
    df %>%
      kruskal_test(formula) %>%
      adjust_pvalue(method = "BH") %>%
      add_significance("p.adj") %>%
      mutate(variable = var)
  })
}

# Run Kruskal-Wallis for all alpha metrics
kruskal_results <- list(
  observed = kruskal_batch(baseline_alpha, krus_vars, "observed"),
  pd = kruskal_batch(baseline_alpha, krus_vars, "pd"),
  shannon = kruskal_batch(baseline_alpha, krus_vars, "diversity_shannon"),
  evenness = kruskal_batch(baseline_alpha, krus_vars, "evenness_pielou")
)

# Save results
write.csv(kruskal_results$observed, "output_results/kruskal_baseobserved_df.csv", row.names = FALSE)
write.csv(kruskal_results$pd, "output_results/kruskal_basepd_df.csv", row.names = FALSE)
write.csv(kruskal_results$shannon, "output_results/kruskal_baseshannon_df.csv", row.names = FALSE)
write.csv(kruskal_results$evenness, "output_results/kruskal_baseevenness_df.csv", row.names = FALSE)

# --- Helper: Dunn’s Test for Significant Variables ---
dunn_compare <- function(df, metric, variable) {
  df %>%
    dunn_test(as.formula(paste(metric, "~", variable)), p.adjust.method = "BH")
}

# Example:
# dunn_compare(baseline_alpha, "observed", "BMI_3_class")

# --- Plot Function for Alpha Diversity ---
plot_alpha_box <- function(df, xvar, yvar, title, fill_colors, comparisons) {
  ggplot(df, aes(x = !!sym(xvar), y = !!sym(yvar), fill = !!sym(xvar))) +
    geom_boxplot(outlier.shape = NA) +
    labs(y = "", title = paste0(title, "\n")) +
    theme_classic(base_size = 9) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 8.5),
      axis.title.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ) +
    scale_fill_manual(values = fill_colors) +
    stat_compare_means(
      comparisons = comparisons,
      method = "wilcox.test",
      label = "p.sgnif",
      symnum.args = list(
        cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        symbols = c("****", "***", "**", "*", "ns")
      )
    )
}

# --- BMI Plots (Baseline) ---
baseline_alpha$BMI_3_class <- factor(
  baseline_alpha$BMI_3_class,
  levels = c("Normal", "Overweight", "Obesity"),
  labels = c("Healthy weight", "Overweight", "Obesity")
)

bmi_comparisons <- list(c("Healthy weight", "Obesity"), c("Overweight", "Obesity"))
bmi_colors <- c("#4daf4a", "#377eb8", "#e41a1c")

BMI_observed <- plot_alpha_box(baseline_alpha, "BMI_3_class", "observed", "Observed richness", bmi_colors, bmi_comparisons)
BMI_pd <- plot_alpha_box(baseline_alpha, "BMI_3_class", "pd", "Faith PD", bmi_colors, list(c("Overweight", "Obesity")))
BMI_shannon <- plot_alpha_box(baseline_alpha, "BMI_3_class", "diversity_shannon", "Shannon diversity", bmi_colors, bmi_comparisons)

(BMI_observed | BMI_shannon | BMI_pd)
ggsave("output_results/BMI_alpha_diversity_baseline.tiff", dpi = 600, width = 8.5, height = 4, bg = "white")

# --- Enterotype Plots (Baseline, Week12, Week24) ---
enterotype_levels <- c("Enterotype one", "Enterotype two", "Enterotype three")
enterotype_labels <- c("Enterotype\none", "Enterotype\ntwo", "Enterotype\nthree")
enterotype_colors <- c("#1B9E77", "#D95F02", "#006BA3")
enterotype_comparisons <- list(
  c("Enterotype\none", "Enterotype\ntwo"),
  c("Enterotype\none", "Enterotype\nthree"),
  c("Enterotype\ntwo", "Enterotype\nthree")
)

plot_enterotype_set <- function(df, timepoint, filename) {
  df$Enterotype <- factor(df$Enterotype, levels = enterotype_levels, labels = enterotype_labels)

  p1 <- plot_alpha_box(df, "Enterotype", "observed", "Observed richness", enterotype_colors, enterotype_comparisons)
  p2 <- plot_alpha_box(df, "Enterotype", "diversity_shannon", "Shannon diversity", enterotype_colors, enterotype_comparisons)
  p3 <- plot_alpha_box(df, "Enterotype", "pd", "Faith PD", enterotype_colors, enterotype_comparisons)

  combined_plot <- p1 | p2 | p3
  ggsave(paste0("output_results/ENTEROTYPE_alpha_diversity_", timepoint, ".tiff"),
         plot = combined_plot, dpi = 600, width = 6.8, height = 3.7, bg = "white")
}

plot_enterotype_set(baseline_alpha, "baseline")
plot_enterotype_set(week12_alphadiversity, "week12")
plot_enterotype_set(week24_alphadiversity, "week24")
