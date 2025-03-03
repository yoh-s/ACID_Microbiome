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
pspairgenus_clr <- microbiome::transform(pspairgenus, transform = 'clr')

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

fulldata_clrselect <- fulldata_clr %>% dplyr::select(c("id", "subject", "Treatment_Type","Time", "Age_quartile",
                                                       "Baseline_enterotype_BIC","value",matches("ASV")))

# We want to pivot the data to long format and include treatment, time, and enterotype
effect_intervention <- fulldata_clrselect %>%
  group_by(ASV) %>%
  nest() %>%
  mutate(model = map(data, ~ glmmTMB(value ~ Treatment_Type * Time + (1|subject), 
                                     family = gaussian, data = .x)))

# Create contrast matrix
contrast_results <- effect_intervention %>%
  mutate(emmeans_intervention = map(model, ~ emmeans(.x, ~ Time | Treatment_Type)),
         
         # Comparison: Intervention - Baseline vs Endline
         Intervention_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>% 
                                       summary() %>% 
                                       filter(Treatment_Type == 'Intervention')),
         Intervention_estimate = map_dbl(Intervention_contrast, ~ ifelse(nrow(.) > 0, pull(., estimate), NA)),
         Intervention_pvalue = map_dbl(Intervention_contrast, ~ ifelse(nrow(.) > 0, pull(., p.value), NA)),
         
         # Comparison: Placebo - Baseline vs Endline
         Placebo_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>% 
                                  summary() %>% 
                                  filter(Treatment_Type == 'Placebo')),
         Placebo_estimate = map_dbl(Placebo_contrast, ~ ifelse(nrow(.) > 0, pull(., estimate), NA)),
         Placebo_pvalue = map_dbl(Placebo_contrast, ~ ifelse(nrow(.) > 0, pull(., p.value), NA)),
         
         # Comparison: Endline - Intervention vs Placebo
         Endline_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Treatment_Type") %>% 
                                  summary() %>% 
                                  filter(Time == '24-wks')),
         Endline_estimate = map_dbl(Endline_contrast, ~ ifelse(nrow(.) > 0, pull(., estimate), NA)),
         Endline_pvalue = map_dbl(Endline_contrast, ~ ifelse(nrow(.) > 0, pull(., p.value), NA))
  )
# View results
contrast_results

# Adjust p-values using BH method
contrast_results1 <- contrast_results %>%
  mutate(
    Intervention_adjusted_pvalue = p.adjust(Intervention_pvalue, method = "BH"),
    Placebo_adjusted_pvalue = p.adjust(Placebo_pvalue, method = "BH"),
    Endline_adjusted_pvalue = p.adjust(Endline_pvalue, method = "BH")
  )

# View adjusted results
View(contrast_results1)

contrast_results2 <- contrast_results1 %>% dplyr::select(-c("model", "data", "emmeans_intervention", "Placebo_contrast", "Intervention_contrast", 
                                      "Endline_contrast"))
##Semi join to filter only those useful
contrast_BMI_ALL_2_taxa <- tax_filtercomp1Gen1 %>% semi_join(contrast_results2, by = "ASV")

intervention_difference_genus <- full_join(contrast_BMI_ALL_2_taxa, contrast_results2, by = "ASV")

write.csv(intervention_difference_genus, "BMI_category/intervention_difference_genus.csv")

#############################################################################################################
fulldata_clrselect <- fulldata_clr %>% dplyr::select(c("id", "subject", "Treatment_Type","Time", "Age_quartile",
                                                       "Baseline_enterotype_BIC","value",matches("ASV")))

# We want to pivot the data to long format and include treatment, time, and enterotype
effect_intervention_agequartile <- fulldata_clrselect %>%
  group_by(ASV) %>%
  nest() %>%
  mutate(model = map(data, ~ glmmTMB(value ~ Treatment_Type * Time * Age_quartile + (1|subject), 
                                     family = gaussian, data = .x)))

# Create contrast matrix
contrast_results_age <- effect_intervention_agequartile %>%
  mutate(emmeans_intervention = map(model, ~ emmeans(.x, ~ Time | Treatment_Type + Age_quartile)),
         
         # Comparison: Intervention - Baseline vs Endline
         Intervention_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>% 
                                       summary() %>% 
                                       filter(Treatment_Type == 'Intervention')),
         Intervention_estimate = map_dbl(Intervention_contrast, ~ ifelse(nrow(.) > 0, pull(., estimate), NA)),
         Intervention_pvalue = map_dbl(Intervention_contrast, ~ ifelse(nrow(.) > 0, pull(., p.value), NA)),
         
         # Comparison: Placebo - Baseline vs Endline
         Placebo_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>% 
                                  summary() %>% 
                                  filter(Treatment_Type == 'Placebo')),
         Placebo_estimate = map_dbl(Placebo_contrast, ~ ifelse(nrow(.) > 0, pull(., estimate), NA)),
         Placebo_pvalue = map_dbl(Placebo_contrast, ~ ifelse(nrow(.) > 0, pull(., p.value), NA)),
         
         # Comparison: Endline - Intervention vs Placebo
         Endline_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Treatment_Type") %>% 
                                  summary() %>% 
                                  filter(Time == '24-wks')),
         Endline_estimate = map_dbl(Endline_contrast, ~ ifelse(nrow(.) > 0, pull(., estimate), NA)),
         Endline_pvalue = map_dbl(Endline_contrast, ~ ifelse(nrow(.) > 0, pull(., p.value), NA))
  )

# View results
contrast_results_age

# Adjust p-values using BH method
contrast_results_age_1 <- contrast_results_age %>%
  mutate(
    Intervention_adjusted_pvalue = p.adjust(Intervention_pvalue, method = "BH"),
    Placebo_adjusted_pvalue = p.adjust(Placebo_pvalue, method = "BH"),
    Endline_adjusted_pvalue = p.adjust(Endline_pvalue, method = "BH")
  )

# View adjusted results
View(contrast_results_age_1)

contrast_results_age_2 <- contrast_results_age_1 %>% dplyr::select(-c("model", "data", "emmeans_intervention", "Placebo_contrast", "Intervention_contrast", 
                                                            "Endline_contrast"))
##Semi join to filter only those useful
contrast_BMI_ALL_2_taxa <- tax_filtercomp1Gen1 %>% semi_join(contrast_results_age_2, by = "ASV")

intervention_agequartile_genus <- full_join(contrast_BMI_ALL_2_taxa, contrast_results_age_2, by = "ASV")

write.csv(intervention_agequartile_genus, "BMI_category/intervention_agequartile_genus.csv")



###########################################################################################################
# We want to pivot the data to long format and include treatment, time, and enterotype
effect_intervention <- fulldata_clrselect %>%
  group_by(ASV) %>%
  nest() %>%
  mutate(model = map(data, ~ glmmTMB(value ~ Treatment_Type * Time + (1|subject), 
                                     family = gaussian, data = .x)))

# Create contrast matrix
contrast_results <- effect_intervention %>%
  mutate(emmeans_intervention = map(model, ~ emmeans(.x, ~ Time | Treatment_Type)),
         
         # Comparison: Intervention - Baseline vs Endline
         Intervention_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>% 
                                       summary() %>% 
                                       filter(Treatment_Type == 'Intervention')),
         Intervention_estimate = map_dbl(Intervention_contrast, ~ ifelse(nrow(.) > 0, pull(., estimate), NA)),
         Intervention_pvalue = map_dbl(Intervention_contrast, ~ ifelse(nrow(.) > 0, pull(., p.value), NA)),
         Intervention_confint = map(emmeans_intervention, ~ confint(contrast(., interaction = "pairwise", simple = "Time")) %>% 
                                      filter(Treatment_Type == 'Intervention')),
         Intervention_lower = map_dbl(Intervention_confint, ~ ifelse(nrow(.) > 0, pull(., lower.CL), NA)),
         Intervention_upper = map_dbl(Intervention_confint, ~ ifelse(nrow(.) > 0, pull(., upper.CL), NA)),
         
         # Comparison: Placebo - Baseline vs Endline
         Placebo_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>% 
                                  summary() %>% 
                                  filter(Treatment_Type == 'Placebo')),
         Placebo_estimate = map_dbl(Placebo_contrast, ~ ifelse(nrow(.) > 0, pull(., estimate), NA)),
         Placebo_pvalue = map_dbl(Placebo_contrast, ~ ifelse(nrow(.) > 0, pull(., p.value), NA)),
         Placebo_confint = map(emmeans_intervention, ~ confint(contrast(., interaction = "pairwise", simple = "Time")) %>% 
                                 filter(Treatment_Type == 'Placebo')),
         Placebo_lower = map_dbl(Placebo_confint, ~ ifelse(nrow(.) > 0, pull(., lower.CL), NA)),
         Placebo_upper = map_dbl(Placebo_confint, ~ ifelse(nrow(.) > 0, pull(., upper.CL), NA)),
         
         # Comparison: Endline - Intervention vs Placebo
         Endline_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Treatment_Type") %>% 
                                  summary() %>% 
                                  filter(Time == '24-wks')),
         Endline_estimate = map_dbl(Endline_contrast, ~ ifelse(nrow(.) > 0, pull(., estimate), NA)),
         Endline_pvalue = map_dbl(Endline_contrast, ~ ifelse(nrow(.) > 0, pull(., p.value), NA)),
         Endline_confint = map(emmeans_intervention, ~ confint(contrast(., interaction = "pairwise", simple = "Treatment_Type")) %>% 
                                 filter(Time == '24-wks')),
         Endline_lower = map_dbl(Endline_confint, ~ ifelse(nrow(.) > 0, pull(., lower.CL), NA)),
         Endline_upper = map_dbl(Endline_confint, ~ ifelse(nrow(.) > 0, pull(., upper.CL), NA))
  )

# View results
contrast_results

# Adjust p-values using BH method
contrast_results1 <- contrast_results %>%
  mutate(
    Intervention_adjusted_pvalue = p.adjust(Intervention_pvalue, method = "BH"),
    Placebo_adjusted_pvalue = p.adjust(Placebo_pvalue, method = "BH"),
    Endline_adjusted_pvalue = p.adjust(Endline_pvalue, method = "BH")
  )

# View adjusted results
View(contrast_results1)

contrast_results2 <- contrast_results1 %>% dplyr::select(-c("model", "data", "emmeans_intervention", "Placebo_contrast", "Intervention_contrast", 
                                                            "Endline_contrast"))

# View final results
View(contrast_results2)
####################################################################################################################
effect_intervention <- fulldata_clrselect %>%
  group_by(ASV) %>%
  nest() %>%
  mutate(model = map(data, ~ glmmTMB(value ~ Treatment_Type * Time + (1|subject), 
                                     family = gaussian, data = .x)))

emmeans_resultsINTERVENTION <- effect_intervention %>%
  mutate(emmeans_intervention = map(model, ~ emmeans(.x, ~ Time | Treatment_Type)))

# Extract EMMs for each combination
emmeans_dataINTERVENTION_final <- emmeans_resultsINTERVENTION %>%
  mutate(emmeans_df = map(emmeans_intervention, as.data.frame)) %>%
  dplyr::select(ASV, emmeans_df) %>%
  unnest(emmeans_df)

# View the extracted EMMs
print(emmeans_dataINTERVENTION_final)

write.csv(emmeans_dataINTERVENTION_final, "output_results/emmeans_dataINTERVENTION_final.csv")

emmeans_time_INTERVENTION <- effect_intervention %>%
  mutate(emmeans_intervention = map(model, ~ emmeans(.x, ~ Treatment_Type | Time)))

# Extract EMMs for each combination
emmeans_time_INTERVENTION_final <- emmeans_time_INTERVENTION %>%
  mutate(emmeans_df = map(emmeans_intervention, as.data.frame)) %>%
  dplyr::select(ASV, emmeans_df) %>%
  unnest(emmeans_df)


# Assuming contrast_results is your data frame
contrastresults_intervention2 <- contrastresults_intervention1 %>%
  select_if(~ !is.list(.))


##Semi join to filter only those useful
contrastresults_intervention2_taxa <- tax_filtercomp1Gen1 %>% semi_join(contrastresults_intervention2, by = "ASV")

intervention_estimate_genus <- full_join(contrastresults_intervention2_taxa, contrastresults_intervention2, by = "ASV")

write.csv(intervention_estimate_genus, "output_results/intervention_estimate_genus.csv")





#######################################################################################################################

df_metadatapspair <- data.frame(sample_data(pspair))
df_metadatapspair$BMI_3_class <- factor(df_metadatapspair$BMI_3_class, levels = 
                                          c("Normal", "Overweight", "Obesity"))
df_metadatapspair$Time <- factor(df_metadatapspair$Time, levels = 
                                   c("Baseline", "24-wks"))
##Effect of BMI observed
observed_bmi_treatment <- lmer(observed ~ Time*Treatment_Type + (1 | subject),
                     data = df_pspairGenusComp_2)

summary(observed_bmi_treatment)

##Effect of BMI shannon diversity
diversity_shannon_treatment <- lmer(diversity_shannon ~ Time*Treatment_Type + (1 | subject),
                              data = df_metadatapspair)

summary(diversity_shannon_treatment)

##Effect of BMI faith pd
pd_treatment <- lmer(pd ~ Time*Treatment_Type + (1 | subject),
               data = df_metadatapspair)

summary(pd_treatment)

##Effect of BMI faith pd
evenness_pielou_treatment <- lmer(evenness_pielou ~ Time*Treatment_Type + (1 | subject),
                            data = df_metadatapspair)

summary(evenness_pielou_treatment)

##Age quartile effect on alpha diversity after treatment
df_metadatapspair <- data.frame(sample_data(pspair))
df_metadatapspair$BMI_3_class <- factor(df_metadatapspair$BMI_3_class, levels = 
                                          c("Normal", "Overweight", "Obesity"))
df_metadatapspair$Time <- factor(df_metadatapspair$Time, levels = 
                                   c("Baseline", "24-wks"))

df_metadatapspair$Age_quartile <- factor(df_metadatapspair$Age_quartile, levels = c("60-64_years", "65-68_years", "68-73_years","73-80_years"))


##Effect of BMI observed
observed_age_treatment <- lmer(observed ~ Age_quartile*Time*Treatment_Type + (1 | subject),
                               data = df_metadatapspair)

summary(observed_age_treatment)

##Effect of BMI shannon diversity
shannon_age_treatment <- lmer(diversity_shannon ~ Age_quartile*Time*Treatment_Type + (1 | subject),
                                    data = df_metadatapspair)

summary(shannon_age_treatment)

##Effect of BMI faith pd
pd_age_treatment <- lmer(pd ~ Age_quartile*Time*Treatment_Type + (1 | subject),
                     data = df_pspairGenusComp_2)

summary(pd_age_treatment)

##Effect of BMI faith pd
age_pielou_treatment <- lmer(evenness_pielou ~ Age_quartile*Time*Treatment_Type + (1 | subject),
                                  data = df_metadatapspair)

summary(age_pielou_treatment)

##Beta diversity 
pspair

##Gloom at the genus level
pspairGenus <- tax_glom(pspair, taxrank = "Genus")

##relative abundance
pspairGenusComp <- microbiome::transform(pspairGenus, transform = "compositional")

##Bray-Curtis dissimilarity
pspairGenusCompbray <- phyloseq::distance(pspairGenusComp, method = "bray")

bray_genus <- adonis2(pspairGenusCompbray~Treatment_Type*Time, data = df_pspairGenusComp_2, strata = df_pspairGenusComp$subject, by= "terms")

##pspair
##Bray-Curtis dissimilarity (treatment only)
pspair_asv_bray <- phyloseq::distance(pspair, method = "bray")

bray_asv <- adonis2(pspair_asv_bray~Treatment_Type*Time, data = df_pspairGenusComp_2, strata = df_pspairGenusComp$subject, by= "terms")

##Jaccard distance
pspair_asv_jaccard <- phyloseq::distance(pspair, method = "jaccard", binary = TRUE)

jaccard_asv <- adonis2(pspair_asv_jaccard~Treatment_Type*Time, data = df_pspairGenusComp_2, strata = df_pspairGenusComp$subject, by= "terms")

##Weighted UniFrac
pspair_asv_wunifrac <- phyloseq::distance(pspair, method = "wunifrac")

wunifrac_asv <- adonis2(pspair_asv_wunifrac~Treatment_Type*Time, data = df_pspairGenusComp_2, strata = df_pspairGenusComp$subject, by= "terms")

##Unweighted UniFrac
pspair_asv_unifrac <- phyloseq::distance(pspair, method = "unifrac")

unifrac_asv <- adonis2(pspair_asv_unifrac~Treatment_Type*Time, data = df_pspairGenusComp_2, strata = df_pspairGenusComp$subject, by= "terms")

##pspair
##Bray-Curtis dissimilarity (treatment only and age)
pspairGenusCompbray <- phyloseq::distance(pspairGenusComp, method = "bray")

bray_genus_age <- adonis2(pspairGenusCompbray~Age_quartile*Treatment_Type*Time, data = df_pspairGenusComp_2, strata = df_pspairGenusComp$subject, by= "terms")

pspair_asv_bray <- phyloseq::distance(pspair, method = "bray")

bray_asvage <- adonis2(pspair_asv_bray~Age_quartile*Treatment_Type*Time, data = df_pspairGenusComp_2, strata = df_pspairGenusComp$subject, by= "terms")

##Jaccard distance
pspair_asv_jaccard <- phyloseq::distance(pspair, method = "jaccard", binary = TRUE)

jaccard_asv_age <- adonis2(pspair_asv_jaccard~Age_quartile*Treatment_Type*Time, data = df_pspairGenusComp_2, strata = df_pspairGenusComp$subject, by= "terms")

##Weighted UniFrac
pspair_asv_wunifrac <- phyloseq::distance(pspair, method = "wunifrac")

wunifrac_asv <- adonis2(pspair_asv_wunifrac~Age_quartile*Treatment_Type*Time, data = df_pspairGenusComp_2, strata = df_pspairGenusComp$subject, by= "terms")

##Unweighted UniFrac
pspair_asv_unifrac <- phyloseq::distance(pspair, method = "unifrac")

unifrac_asv <- adonis2(pspair_asv_unifrac~Age_quartile*Treatment_Type*Time, data = df_pspairGenusComp_2, strata = df_pspairGenusComp$subject, by= "terms")
