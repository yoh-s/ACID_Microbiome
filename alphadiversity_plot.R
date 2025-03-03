library(ggplot2)
library(ggpubr)
library(magrittr)

#################################BASELINE#################################
baseline_alphadiversity$BMI_3_class <- factor(baseline_alphadiversity$BMI_3_class, levels = c("Normal", "Overweight", "Obesity"))

levels(baseline_alphadiversity$BMI_3_class)[levels(baseline_alphadiversity$BMI_3_class) == "Normal"] <- "Healthy weight"


bmi_comparison <- list(c("Healthy weight", "Overweight"), c("Healthy weight", "Obesity"), c("Overweight", "Obesity"))
sig_bmi_bmi <- list(c("Healthy weight", "Obesity"), c("Overweight", "Obesity"))
p_value_cut = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

#observed richness (BMI)
BMI_observed_all <- baseline_alphadiversity %>%
  ggplot(aes(x=BMI_3_class, y=observed, fill=BMI_3_class)) +
  geom_boxplot(outlier.shape = NA) +
  labs(y = "", title = "Observed richness\n") +
  theme_classic() +
  theme(legend.position="none", axis.text.x = element_text(size = 8.5), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8.5), plot.title = element_text(hjust = 0.5, size = 9), plot.margin = margin(t=0,r=0,b=0,l=-30)) +
  scale_fill_manual(values = c("#4daf4a", "#377eb8", "#e41a1c")) +
  stat_compare_means(comparisons = sig_bmi_bmi, method = "wilcox.test",
                     label = "p.sgnif", symnum.args = p_value_cut)

sig_faithpd_bmi <- list(c("Overweight", "Obesity"))
#faith pd
BMI_faithpd_all <- baseline_alphadiversity %>%
  ggplot(aes(x=BMI_3_class, y=pd, fill=BMI_3_class)) +
  geom_boxplot(outlier.shape = NA) +
  labs(y = "", title = "Faith PD\n") +
  theme_classic() +
  theme(legend.position="none", axis.text.x = element_text(size = 8.5), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8.5), plot.title = element_text(hjust = 0.5, size = 9), plot.margin = margin(t=0,r=0,b=0,l=-30)) +
  scale_fill_manual(values = c("#4daf4a", "#377eb8", "#e41a1c")) +
  stat_compare_means(comparisons = sig_faithpd_bmi, method = "wilcox.test",
                     label = "p.sgnif", symnum.args = p_value_cut)

#Shannon diversity
BMI_shannon_all <- baseline_alphadiversity %>%
  ggplot(aes(x=BMI_3_class, y=diversity_shannon, fill=BMI_3_class)) +
  geom_boxplot(outlier.shape = NA) +
  labs(y = "", title = "Shannon diversity\n") +
  theme_classic() +
  theme(legend.position="none", axis.text.x = element_text(size = 8.5), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8.5), plot.title = element_text(hjust = 0.5, size = 9), plot.margin = margin(t=0,r=0,b=0,l=-30)) +
  scale_fill_manual(values = c("#4daf4a", "#377eb8", "#e41a1c")) +
  stat_compare_means(comparisons = sig_bmi_bmi, method = "wilcox.test",
                     label = "p.sgnif", symnum.args = p_value_cut)
##plot only three diversity metrics
BMI_observed_all | BMI_shannon_all | BMI_faithpd_all

ggsave("output_results/BMI_alpha_diversity_baseline.tiff", dpi = 600, height = 4, width = 8.5)

# Define the comparisons between enterotypes
entero_comparison <- list(c("Enterotype\none", "Enterotype\ntwo"), 
                          c("Enterotype\none", "Enterotype\nthree"),
                          c("Enterotype\ntwo", "Enterotype\nthree"))


baseline_alphadiversity$Enterotype <- factor(baseline_alphadiversity$Enterotype, 
                                             levels = c("Enterotype one", "Enterotype two", "Enterotype three"))

# Modify Enterotype levels
baseline_alphadiversity$Enterotype <- factor(baseline_alphadiversity$Enterotype,
                                             levels = c("Enterotype one", "Enterotype two", "Enterotype three"),
                                             labels = c("Enterotype\none", "Enterotype\ntwo", "Enterotype\nthree"))

#observed richness (Enterotype)
ENTEROTYPE_observed_all <- baseline_alphadiversity %>%
  ggplot(aes(x=Enterotype, y=observed, fill=Enterotype)) +
  geom_boxplot(outlier.shape = NA) +
  labs(y = "", title = "Observed richness\n") +
  theme_classic() +
  theme(legend.position="none", axis.text.x = element_text(size = 8.5), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8.5), plot.title = element_text(hjust = 0.5, size = 9), plot.margin = margin(t=0,r=0,b=0,l=-30)) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#006BA3")) +
  stat_compare_means(comparisons = entero_comparison, method = "wilcox.test",
                     label = "p.sgnif", symnum.args = p_value_cut)

  
sig_faithpd_enterotype <- list(c("Enterotype one", "Enterotype two"), c("Enterotype one", "Enterotype three"))

# Define the comparisons between enterotypes
sig_faithpd_enterotype <- list(c("Enterotype\none", "Enterotype\ntwo"),
                               c("Enterotype\none", "Enterotype\nthree"))


#faith pd
ENTEROTYPE_faithpd_all <- baseline_alphadiversity %>%
  ggplot(aes(x=Enterotype, y=pd, fill=Enterotype)) +
  geom_boxplot(outlier.shape = NA) +
  labs(y = "", title = "Faith PD\n") +
  theme_classic() +
  theme(legend.position="none", axis.text.x = element_text(size = 8.5), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8.5), plot.title = element_text(hjust = 0.5, size = 9), plot.margin = margin(t=0,r=0,b=0,l=-30)) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#006BA3")) +
  stat_compare_means(comparisons = sig_faithpd_enterotype, method = "wilcox.test",
                     label = "p.sgnif", symnum.args = p_value_cut)

#Shannon diversity
ENTEROTYPE_shannon_all <- baseline_alphadiversity %>%
  ggplot(aes(x=Enterotype, y=diversity_shannon, fill=Enterotype)) +
  geom_boxplot(outlier.shape = NA) +
  labs(y = "", title = "Shannon diversity\n") +
  theme_classic() +
  theme(legend.position="none", axis.text.x = element_text(size = 8.5), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8.5), plot.title = element_text(hjust = 0.5, size = 9), plot.margin = margin(t=0,r=0,b=0,l=-30)) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#006BA3")) +
  stat_compare_means(comparisons = entero_comparison, method = "wilcox.test",
                     label = "p.sgnif", symnum.args = p_value_cut)
##plot only three diversity metrics
ENTEROTYPE_observed_all | ENTEROTYPE_shannon_all | ENTEROTYPE_faithpd_all

##Save plot
ggsave("output_results/ENTEROTYPE_alpha_diversity_baseline.tiff",dpi = 600)

##Week 12
week12_alphadiversity

week12_alphadiversity$Enterotype <- factor(week12_alphadiversity$Enterotype,
                                             levels = c("Enterotype one", "Enterotype two", "Enterotype three"),
                                             labels = c("Enterotype\none", "Enterotype\ntwo", "Enterotype\nthree"))


#observed richness (Enterotype)
ENTEROTYPE_observed_wk12 <- week12_alphadiversity %>%
  ggplot(aes(x=Enterotype, y=observed, fill=Enterotype)) +
  geom_boxplot(outlier.shape = NA) +
  labs(y = "", title = "Observed richness\n") +
  theme_classic() +
  theme(legend.position="none", axis.text.x = element_text(size = 8.5), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8.5), plot.title = element_text(hjust = 0.5, size = 9), plot.margin = margin(t=0,r=0,b=0,l=-30)) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#006BA3")) +
  stat_compare_means(comparisons = entero_comparison, method = "wilcox.test",
                     label = "p.sgnif", symnum.args = p_value_cut)


sig_faithpd_enterotype <- list(c("Enterotype one", "Enterotype two"), c("Enterotype one", "Enterotype three"))

# Define the comparisons between enterotypes
sig_faithpd_enterotype <- list(c("Enterotype\none", "Enterotype\ntwo"),
                               c("Enterotype\none", "Enterotype\nthree"))


#faith pd
ENTEROTYPE_faithpd_wk12 <- week12_alphadiversity %>%
  ggplot(aes(x=Enterotype, y=pd, fill=Enterotype)) +
  geom_boxplot(outlier.shape = NA) +
  labs(y = "", title = "Faith PD\n") +
  theme_classic() +
  theme(legend.position="none", axis.text.x = element_text(size = 8.5), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8.5), plot.title = element_text(hjust = 0.5, size = 9), plot.margin = margin(t=0,r=0,b=0,l=-30)) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#006BA3")) +
  stat_compare_means(comparisons = entero_comparison, method = "wilcox.test",
                     label = "p.sgnif", symnum.args = p_value_cut)

#Shannon diversity
ENTEROTYPE_shannon_wk12 <- week12_alphadiversity %>%
  ggplot(aes(x=Enterotype, y=diversity_shannon, fill=Enterotype)) +
  geom_boxplot(outlier.shape = NA) +
  labs(y = "", title = "Shannon diversity\n") +
  theme_classic() +
  theme(legend.position="none", axis.text.x = element_text(size = 8.5), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8.5), plot.title = element_text(hjust = 0.5, size = 9), plot.margin = margin(t=0,r=0,b=0,l=-30)) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#006BA3")) +
  stat_compare_means(comparisons = entero_comparison, method = "wilcox.test",
                     label = "p.sgnif", symnum.args = p_value_cut)
##plot only three diversity metrics
ENTEROTYPE_observed_wk12 | ENTEROTYPE_shannon_wk12 | ENTEROTYPE_faithpd_wk12

##Save plot
ggsave("output_results/ENTEROTYPE_alpha_diversity_week12.tiff",dpi = 600)
##6.81 x 3.71
##Week 24
week24_alphadiversity
week24_alphadiversity$Enterotype <- factor(week24_alphadiversity$Enterotype,
                                           levels = c("Enterotype one", "Enterotype two", "Enterotype three"),
                                           labels = c("Enterotype\none", "Enterotype\ntwo", "Enterotype\nthree"))


#observed richness (Enterotype)
ENTEROTYPE_observed_wk24 <- week24_alphadiversity %>%
  ggplot(aes(x=Enterotype, y=observed, fill=Enterotype)) +
  geom_boxplot(outlier.shape = NA) +
  labs(y = "", title = "Observed richness\n") +
  theme_classic() +
  theme(legend.position="none", axis.text.x = element_text(size = 8.5), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8.5), plot.title = element_text(hjust = 0.5, size = 9), plot.margin = margin(t=0,r=0,b=0,l=-30)) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#006BA3")) +
  stat_compare_means(comparisons = entero_comparison, method = "wilcox.test",
                     label = "p.sgnif", symnum.args = p_value_cut)


sig_faithpd_enterotype <- list(c("Enterotype one", "Enterotype two"), c("Enterotype one", "Enterotype three"))

# Define the comparisons between enterotypes
sig_faithpd_enterotype <- list(c("Enterotype\none", "Enterotype\ntwo"),
                               c("Enterotype\none", "Enterotype\nthree"))


#faith pd
ENTEROTYPE_faithpd_wk24 <- week24_alphadiversity %>%
  ggplot(aes(x=Enterotype, y=pd, fill=Enterotype)) +
  geom_boxplot(outlier.shape = NA) +
  labs(y = "", title = "Faith PD\n") +
  theme_classic() +
  theme(legend.position="none", axis.text.x = element_text(size = 8.5), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8.5), plot.title = element_text(hjust = 0.5, size = 9), plot.margin = margin(t=0,r=0,b=0,l=-30)) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#006BA3")) +
  stat_compare_means(comparisons = entero_comparison, method = "wilcox.test",
                     label = "p.sgnif", symnum.args = p_value_cut)

#Shannon diversity
ENTEROTYPE_shannon_wk24 <- week24_alphadiversity %>%
  ggplot(aes(x=Enterotype, y=diversity_shannon, fill=Enterotype)) +
  geom_boxplot(outlier.shape = NA) +
  labs(y = "", title = "Shannon diversity\n") +
  theme_classic() +
  theme(legend.position="none", axis.text.x = element_text(size = 8.5), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8.5), plot.title = element_text(hjust = 0.5, size = 9), plot.margin = margin(t=0,r=0,b=0,l=-30)) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#006BA3")) +
  stat_compare_means(comparisons = entero_comparison, method = "wilcox.test",
                     label = "p.sgnif", symnum.args = p_value_cut)
##plot only three diversity metrics
ENTEROTYPE_observed_wk24 | ENTEROTYPE_shannon_wk24 | ENTEROTYPE_faithpd_wk24

##Save plot
ggsave("output_results/ENTEROTYPE_alpha_diversity_week24.tiff",dpi = 600)

ps_study810_rare_f2

###select a paired
