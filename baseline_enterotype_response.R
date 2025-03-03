library(glmmTMB)
library(emmeans)
library(multcomp)
library(nortest)
library(lme4)
library(tidyr)
library(dplyr)
library(broom)
library(vegan)

##read the metadata file
enterotype_BIC <- read_csv("enterotype_BIC.csv")

ps_study5_filter

enterotype_BIC <- enterotype_BIC %>% dplyr::select(-Time)
ps_study5_filter_1 <- ps_join(ps_study5_filter, enterotype_BIC, by = "sampleid")

metadatapaired1 <- metadatapaired %>% dplyr::select(-Baseline_enterotypefinal)

ps_study5_filter_2 <- ps_join(ps_study5_filter_1, metadatapaired1, by = "sampleid")

##Glom at the genus level
ps_study5filter2genus <- tax_glom(ps_study5_filter_2, taxrank = "Genus")
sum(taxa_sums(ps_study5filter2genus) == 0) ##no zero sum

##Enterotype baseline
sample_data(ps_study5filter2genus)$Baseline_enterotype_BIC <- factor(sample_data(ps_study5filter2genus)$Baseline_enterotype_BIC,
                                                levels = c("Enterotype one","Enterotype two"))

##Intervention groups
sample_data(ps_study5filter2genus)$Treatment_Type <- factor(sample_data(ps_study5filter2genus)$Treatment_Type,
                                                       levels = c("Intervention","Placebo"))

##Intervention groups
sample_data(ps_study5filter2genus)$Time <- factor(sample_data(ps_study5filter2genus)$Time,
                                                             levels = c("Baseline","24-wks"))

pspairgenus <- subset_samples(ps_study5filter2genus, Paired_final == "1")

df_pspairgenus <- psmelt(pspairgenus)

##transform clr
pspairgenus_compositional <- microbiome::transform(pspairgenus, transform = 'compositional')

# Prepare unrarefied count data
pspairgenus_counts <- pspairgenus_clr %>% 
  otu_table() %>% t() %>% 
  as.data.frame() %>% 
  rownames_to_column('id')

mdata_pairgenus_counts <- sample_data(pspairgenus_clr) %>%
  data.frame() %>%
  rownames_to_column('id') 

fulldata <- left_join(mdata_pairgenus_counts, pspairgenus_counts)

# Centered log-ratio transformation
fulldata_clr <- fulldata %>%
  pivot_longer(cols = starts_with('ASV'), names_to = 'ASV', values_to = 'value') %>%
  group_by(id) %>% 
  ungroup()

fulldata_clrselect <- fulldata_clr %>% dplyr::select(c("id", "subject", "Treatment_Type","Time",
                                                       "Baseline_enterotype_BIC","value",matches("ASV")))

# We want to pivot the data to long format and include treatment, time, and enterotype
lmm_data <- fulldata_clrselect %>%
  group_by(ASV) %>%
  nest() %>%
  mutate(model = map(data, ~ glmmTMB(value ~ Treatment_Type * Time * Baseline_enterotype_BIC + (1|subject), 
                                     family = gaussian, data = .x)))

# Step 4: Post-hoc Comparisons using `emmeans`
# For each model, we want to get specific pairwise comparisons
# Enterotype 1 - Baseline vs Endline for Intervention and Placebo
# Enterotype 2 - Baseline vs Endline for Intervention and Placebo

# Step 4: Post-hoc Comparisons using `emmeans`
comparison_results <- lmm_data %>%
  mutate(emmeans_intervention = map(model, ~ emmeans(.x, ~ Time | Baseline_enterotype_BIC | Treatment_Type, at = list(Time = c("Baseline", "24-wks")))),
         
         # Enterotype one: Intervention - Baseline vs Endline (Estimate and P-value)
         E1_Intervention_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>% 
                                          summary() %>% 
                                          filter(Baseline_enterotype_BIC == 'Enterotype one', Treatment_Type == 'Intervention')),
         E1_Intervention_estimate = map_dbl(E1_Intervention_contrast, ~ pull(., estimate)),
         E1_Intervention_pvalue = map_dbl(E1_Intervention_contrast, ~ pull(., p.value)),
         
         # Enterotype one: Placebo - Baseline vs Endline (Estimate and P-value)
         E1_Placebo_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>% 
                                     summary() %>% 
                                     filter(Baseline_enterotype_BIC == 'Enterotype one', Treatment_Type == 'Placebo')),
         E1_Placebo_estimate = map_dbl(E1_Placebo_contrast, ~ pull(., estimate)),
         E1_Placebo_pvalue = map_dbl(E1_Placebo_contrast, ~ pull(., p.value)),
         
         # Enterotype two: Intervention - Baseline vs Endline (Estimate and P-value)
         E2_Intervention_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>% 
                                          summary() %>% 
                                          filter(Baseline_enterotype_BIC == 'Enterotype two', Treatment_Type == 'Intervention')),
         E2_Intervention_estimate = map_dbl(E2_Intervention_contrast, ~ pull(., estimate)),
         E2_Intervention_pvalue = map_dbl(E2_Intervention_contrast, ~ pull(., p.value)),
         
         # Enterotype two: Placebo - Baseline vs Endline (Estimate and P-value)
         E2_Placebo_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>% 
                                     summary() %>% 
                                     filter(Baseline_enterotype_BIC == 'Enterotype two', Treatment_Type == 'Placebo')),
         E2_Placebo_estimate = map_dbl(E2_Placebo_contrast, ~ pull(., estimate)),
         E2_Placebo_pvalue = map_dbl(E2_Placebo_contrast, ~ pull(., p.value))
  )

# Step 5: Adjust p-values for multiple comparisons (BH method)
comparison_results_1 <- comparison_results %>%
  mutate(E1_Intervention_pvalue_adj = p.adjust(E1_Intervention_pvalue, method = 'BH'),
         E1_Placebo_pvalue_adj = p.adjust(E1_Placebo_pvalue, method = 'BH'),
         E2_Intervention_pvalue_adj = p.adjust(E2_Intervention_pvalue, method = 'BH'),
         E2_Placebo_pvalue_adj = p.adjust(E2_Placebo_pvalue, method = 'BH'))
# Step 6: View results
enterotype_individual_taxa <- comparison_results_1 %>%
  dplyr::select(ASV, E1_Intervention_estimate, E1_Intervention_pvalue, E1_Placebo_estimate,
                E1_Placebo_pvalue, E2_Intervention_estimate, E2_Intervention_pvalue,
                E2_Placebo_estimate, E2_Placebo_pvalue, E1_Intervention_pvalue_adj, E1_Placebo_pvalue_adj,
                E2_Intervention_pvalue_adj, E2_Placebo_pvalue_adj)
##export the tax table 

tax_filtercomp1Gen <- data.frame(tax_table(pspairgenus_clr))

tax_filtercomp1Gen1 <- tax_filtercomp1Gen %>% dplyr::select(-Species) %>% rownames_to_column("ASV")

##Semi join to filter only those useful
tax_filtercomp1Gen125 <- tax_filtercomp1Gen1 %>% semi_join(enterotype_individual_taxa, by = "ASV")

effectenterotypeIndividual <- full_join(tax_filtercomp1Gen125, enterotype_individual_taxa, by = "ASV")

write.csv(effectenterotypeIndividual, "enterotype_final/effectenterotypeIndividual.csv")


###week-24 placebo vs intervention at enterotype one and enterotype two

# Step 4: Post-hoc Comparisons for Week-24 Placebo vs Intervention in Each Enterotype
comparison_endline <- lmm_data %>%
  mutate(emmeans_week24 = map(model, ~ emmeans(.x, ~ Treatment_Type | Baseline_enterotype_BIC | Time, 
                                               at = list(Time = "24-wks"))),
         # Enterotype one: Week-24 Placebo vs Intervention (Estimate and p-value)
         E1_Week24_Placebo_vs_Intervention = map(emmeans_week24, ~ contrast(., "pairwise") %>% 
                                                   summary() %>% 
                                                   filter(Baseline_enterotype_BIC == 'Enterotype one')),
         E1_Week24_estimate = map_dbl(E1_Week24_Placebo_vs_Intervention, ~ pull(., estimate)),
         E1_Week24_pvalue = map_dbl(E1_Week24_Placebo_vs_Intervention, ~ pull(., p.value)),
         
         # Enterotype two: Week-24 Placebo vs Intervention (Estimate and p-value)
         E2_Week24_Placebo_vs_Intervention = map(emmeans_week24, ~ contrast(., "pairwise") %>% 
                                                   summary() %>% 
                                                   filter(Baseline_enterotype_BIC == 'Enterotype two')),
         E2_Week24_estimate = map_dbl(E2_Week24_Placebo_vs_Intervention, ~ pull(., estimate)),
         E2_Week24_pvalue = map_dbl(E2_Week24_Placebo_vs_Intervention, ~ pull(., p.value))
  )

# Step 5: Adjust p-values for multiple comparisons (BH method)
comparison_endline_final <- comparison_endline %>%
  mutate(E1_Week24_pvalue_adj = p.adjust(E1_Week24_pvalue, method = 'BH'),
         E2_Week24_pvalue_adj = p.adjust(E2_Week24_pvalue, method = 'BH'))

# Step 6: View results
comparison_endline_final12 <- comparison_endline_final %>%
  dplyr::select("ASV", "E1_Week24_estimate", "E1_Week24_pvalue", "E2_Week24_estimate", "E2_Week24_pvalue",
                "E1_Week24_pvalue_adj", "E2_Week24_pvalue_adj")


##Semi join to filter only those useful
tax_filtercomp1Gen125 <- tax_filtercomp1Gen1 %>% semi_join(comparison_endline_final12, by = "ASV")

effectenterotypeweek24 <- full_join(tax_filtercomp1Gen125, comparison_endline_final12, by = "ASV")

write.csv(effectenterotypeweek24, "enterotype_final/effectenterotypeweek24.csv")

##Effect on alpha diversity 
fulldata_alphadiversity<- fulldata_clr %>% dplyr::select(c("id", "subject", "Treatment_Type","Time",
                                                       "Baseline_enterotype_BIC", "diversity_shannon"))

# We want to pivot the data to long format and include treatment, time, and enterotype
##Extract sample data from the phyloseq object
df_pspairgenus_clr <- data.frame(sample_data(pspairgenus_clr))

df_alphadiversity <- df_pspairgenus_clr %>% dplyr::select(c("sampleid", "subject", "Treatment_Type","Time", 
                                       "Baseline_enterotype_BIC", "diversity_shannon", "observed", "pd",
                                       "evenness_pielou"))

##Linear model on alpha diversity
lmm_alpha <- df_alphadiversity %>%
  nest() %>%
  mutate(model = map(data, ~ glmmTMB(evenness_pielou ~ Treatment_Type * Time * Baseline_enterotype_BIC + (1|subject), 
                                     family = gaussian, data = .x)))

# Step 4: Post-hoc Comparisons using `emmeans`
# For each model, we want to get specific pairwise comparisons
# Enterotype 1 - Baseline vs Endline for Intervention and Placebo
# Enterotype 2 - Baseline vs Endline for Intervention and Placebo

# Step 4: Post-hoc Comparisons using `emmeans`
comparison_endline_alpha <- lmm_alpha %>%
  mutate(emmeans_week24 = map(model, ~ emmeans(.x, ~ Treatment_Type | Baseline_enterotype_BIC | Time, 
                                               at = list(Time = "24-wks"))),
         # Enterotype one: Week-24 Placebo vs Intervention (Estimate and p-value)
         E1_Week24_Placebo_vs_Intervention = map(emmeans_week24, ~ contrast(., "pairwise") %>% 
                                                   summary() %>% 
                                                   filter(Baseline_enterotype_BIC == 'Enterotype one')),
         E1_Week24_estimate = map_dbl(E1_Week24_Placebo_vs_Intervention, ~ pull(., estimate)),
         E1_Week24_pvalue = map_dbl(E1_Week24_Placebo_vs_Intervention, ~ pull(., p.value)),
         
         # Enterotype two: Week-24 Placebo vs Intervention (Estimate and p-value)
         E2_Week24_Placebo_vs_Intervention = map(emmeans_week24, ~ contrast(., "pairwise") %>% 
                                                   summary() %>% 
                                                   filter(Baseline_enterotype_BIC == 'Enterotype two')),
         E2_Week24_estimate = map_dbl(E2_Week24_Placebo_vs_Intervention, ~ pull(., estimate)),
         E2_Week24_pvalue = map_dbl(E2_Week24_Placebo_vs_Intervention, ~ pull(., p.value))
  )

# Step 5: Adjust p-values for multiple comparisons (BH method)
comparison_endline_alpha1 <- comparison_endline_alpha %>%
  mutate(E1_Week24_pvalue_adj = p.adjust(E1_Week24_pvalue, method = 'BH'),
         E2_Week24_pvalue_adj = p.adjust(E2_Week24_pvalue, method = 'BH'))

# Step 6: View results
comparison_endline_final12 <- comparison_endline_alpha1 %>%
  dplyr::select("E1_Week24_estimate", "E1_Week24_pvalue", "E2_Week24_estimate", "E2_Week24_pvalue",
                "E1_Week24_pvalue_adj", "E2_Week24_pvalue_adj")
View(comparison_endline_final12)

###Linear mixed models for alpha diversity
##Linear model on alpha diversity
lmm_alpha <- df_alphadiversity %>%
  nest() %>%
  mutate(model = map(data, ~ glmmTMB(diversity_shannon ~ Treatment_Type * Time * Baseline_enterotype_BIC + (1|subject), 
                                     family = gaussian, data = .x)))

# Step 4: Post-hoc Comparisons using `emmeans`
# For each model, we want to get specific pairwise comparisons
# Enterotype 1 - Baseline vs Endline for Intervention and Placebo
# Enterotype 2 - Baseline vs Endline for Intervention and Placebo

# Step 4: Post-hoc Comparisons using `emmeans`
comparison_endline_alpha <- lmm_alpha %>%
  mutate(emmeans_week24 = map(model, ~ emmeans(.x, ~ Treatment_Type | Baseline_enterotype_BIC | Time, 
                                               at = list(Time = "24-wks"))),
         # Enterotype one: Week-24 Placebo vs Intervention (Estimate and p-value)
         E1_Week24_Placebo_vs_Intervention = map(emmeans_week24, ~ contrast(., "pairwise") %>% 
                                                   summary() %>% 
                                                   filter(Baseline_enterotype_BIC == 'Enterotype one')),
         E1_Week24_estimate = map_dbl(E1_Week24_Placebo_vs_Intervention, ~ pull(., estimate)),
         E1_Week24_pvalue = map_dbl(E1_Week24_Placebo_vs_Intervention, ~ pull(., p.value)),
         
         # Enterotype two: Week-24 Placebo vs Intervention (Estimate and p-value)
         E2_Week24_Placebo_vs_Intervention = map(emmeans_week24, ~ contrast(., "pairwise") %>% 
                                                   summary() %>% 
                                                   filter(Baseline_enterotype_BIC == 'Enterotype two')),
         E2_Week24_estimate = map_dbl(E2_Week24_Placebo_vs_Intervention, ~ pull(., estimate)),
         E2_Week24_pvalue = map_dbl(E2_Week24_Placebo_vs_Intervention, ~ pull(., p.value))
  )

# Step 5: Adjust p-values for multiple comparisons (BH method)
comparison_endline_alpha1 <- comparison_endline_alpha %>%
  mutate(E1_Week24_pvalue_adj = p.adjust(E1_Week24_pvalue, method = 'BH'),
         E2_Week24_pvalue_adj = p.adjust(E2_Week24_pvalue, method = 'BH'))

# Step 6: View results
comparison_endline_final12 <- comparison_endline_alpha1 %>%
  dplyr::select("E1_Week24_estimate", "E1_Week24_pvalue", "E2_Week24_estimate", "E2_Week24_pvalue",
                "E1_Week24_pvalue_adj", "E2_Week24_pvalue_adj")
View(comparison_endline_final12)


