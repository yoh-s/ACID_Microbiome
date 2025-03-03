library(lme4)
library(lmerTest)
library(pairwiseAdonis)

ps_study810_rare_f2 <- ps_join(ps_study810_rare_f2, metadatapaired, by = "sampleid")
##4564

##procrsustes anlaysis
##filter pair baseline and endline paired sample
psrare_pair <- subset_samples(ps_study810_rare_f2, Paired_final == "1")
pspair <- psrare_pair

##Remove zero sum
pspair <- prune_taxa(taxa_sums(pspair) > 0, pspair) ##195 ASVs

df_metadatapspair <- data.frame(sample_data(pspair))
df_metadatapspair$BMI_3_class <- factor(df_metadatapspair$BMI_3_class, levels = 
                                          c("Normal", "Overweight", "Obesity"))
df_metadatapspair$Time <- factor(df_metadatapspair$Time, levels = 
                                   c("Baseline", "24-wks"))
##Effect of BMI observed
observed_bmi <- lmer(observed ~ Time*Treatment_Type*BMI_3_class + (1 | subject),
     data = df_pspairGenusComp_2)

summary(observed_bmi)

##Effect of BMI shannon diversity
diversity_shannon_bmi <- lmer(diversity_shannon ~ Time*Treatment_Type*BMI_3_class + (1 | subject),
                     data = df_pspairGenusComp_2)

summary(diversity_shannon_bmi)

##Effect of BMI faith pd
pd_bmi <- lmer(pd ~ Time*Treatment_Type*BMI_3_class + (1 | subject),
                              data = df_pspairGenusComp_2)

summary(pd_bmi)

##Effect of BMI faith pd
evenness_pielou_bmi <- lmer(evenness_pielou ~ Time*Treatment_Type*BMI_3_class + (1 | subject),
               data = df_pspairGenusComp_2)

summary(evenness_pielou_bmi)

################################################################################################
bmi_alphadiversity <- df_metadatapspair %>% dplyr::select(c("sampleid", "subject", "Treatment_Type","Time", 
                                                            "Baseline_BMI3", "diversity_shannon", "observed", "pd", "evenness_pielou"))

bmi_alphadiversity$Baseline_BMI3 <- factor(bmi_alphadiversity$Baseline_BMI3, levels = c("Normal", "Overweight", "Obesity"))  

##Linear model on alpha diversity
lmm_bmialpha <- bmi_alphadiversity %>%
  nest() %>%
  mutate(model = map(data, ~ glmmTMB(diversity_shannon ~ Treatment_Type * Time * Baseline_BMI3 + (1|subject), 
                                     family = gaussian, data = .x)))

# Step 4: Post-hoc Comparisons using `emmeans`
bmi_shannon_base24 <- lmm_bmialpha %>%
  mutate(emmeans_intervention = map(model, ~ emmeans(.x, ~ Time | Baseline_BMI3 | Treatment_Type, at = list(Time = c("Baseline", "24-wks")))),
         
         # BMI Normal: Intervention - Baseline vs Endline (Estimate and P-value)
         Normal_Intervention_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>% 
                                              summary() %>% 
                                              filter(Baseline_BMI3 == 'Normal', Treatment_Type == 'Intervention')),
         Normal_Intervention_estimate = map_dbl(Normal_Intervention_contrast, ~ pull(., estimate)),
         Normal_Intervention_pvalue = map_dbl(Normal_Intervention_contrast, ~ pull(., p.value)),
         
         # BMI Normal: Placebo - Baseline vs Endline (Estimate and P-value)
         Normal_Placebo_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>% 
                                         summary() %>% 
                                         filter(Baseline_BMI3 == 'Normal', Treatment_Type == 'Placebo')),
         Normal_Placebo_estimate = map_dbl(Normal_Placebo_contrast, ~ pull(., estimate)),
         Normal_Placebo_pvalue = map_dbl(Normal_Placebo_contrast, ~ pull(., p.value)),
         
         # BMI Overweight: Intervention - Baseline vs Endline (Estimate and P-value)
         Overweight_Intervention_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>% 
                                                  summary() %>% 
                                                  filter(Baseline_BMI3 == 'Overweight', Treatment_Type == 'Intervention')),
         Overweight_Intervention_estimate = map_dbl(Overweight_Intervention_contrast, ~ pull(., estimate)),
         Overweight_Intervention_pvalue = map_dbl(Overweight_Intervention_contrast, ~ pull(., p.value)),
         
         # BMI Overweight: Placebo - Baseline vs Endline (Estimate and P-value)
         Overweight_Placebo_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>% 
                                             summary() %>% 
                                             filter(Baseline_BMI3 == 'Overweight', Treatment_Type == 'Placebo')),
         Overweight_Placebo_estimate = map_dbl(Overweight_Placebo_contrast, ~ pull(., estimate)),
         Overweight_Placebo_pvalue = map_dbl(Overweight_Placebo_contrast, ~ pull(., p.value)),
         
         # BMI Obesity: Intervention - Baseline vs Endline (Estimate and P-value)
         Obesity_Intervention_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>% 
                                               summary() %>% 
                                               filter(Baseline_BMI3 == 'Obesity', Treatment_Type == 'Intervention')),
         Obesity_Intervention_estimate = map_dbl(Obesity_Intervention_contrast, ~ pull(., estimate)),
         Obesity_Intervention_pvalue = map_dbl(Obesity_Intervention_contrast, ~ pull(., p.value)),
         
         # BMI Obesity: Placebo - Baseline vs Endline (Estimate and P-value)
         Obesity_Placebo_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>% 
                                          summary() %>% 
                                          filter(Baseline_BMI3 == 'Obesity', Treatment_Type == 'Placebo')),
         Obesity_Placebo_estimate = map_dbl(Obesity_Placebo_contrast, ~ pull(., estimate)),
         Obesity_Placebo_pvalue = map_dbl(Obesity_Placebo_contrast, ~ pull(., p.value))
  ) %>%
  # Adjust p-values using the BH method
  mutate(across(ends_with("pvalue"), ~ p.adjust(., method = "BH"), .names = "adj_{col}"))

##no signficant effect due to baseline bmi group
View(bmi_shannon_base24)
#################Week 24 (Intervention vs placebo)
# Step 4: Post-hoc Comparisons for Week-24 Placebo vs Intervention in Each BMI Category
shannon_endline_bmi <- lmm_bmialpha %>%
  mutate(emmeans_week24 = map(model, ~ emmeans(.x, ~ Treatment_Type | Baseline_BMI3 | Time, 
                                               at = list(Time = "24-wks"))),
         
         # BMI Normal: Week-24 Placebo vs Intervention (Estimate and p-value)
         Normal_Week24_Placebo_vs_Intervention = map(emmeans_week24, ~ contrast(., "pairwise") %>% 
                                                       summary() %>% 
                                                       filter(Baseline_BMI3 == 'Normal')),
         Normal_Week24_estimate = map_dbl(Normal_Week24_Placebo_vs_Intervention, ~ pull(., estimate)),
         Normal_Week24_pvalue = map_dbl(Normal_Week24_Placebo_vs_Intervention, ~ pull(., p.value)),
         
         # BMI Overweight: Week-24 Placebo vs Intervention (Estimate and p-value)
         Overweight_Week24_Placebo_vs_Intervention = map(emmeans_week24, ~ contrast(., "pairwise") %>% 
                                                           summary() %>% 
                                                           filter(Baseline_BMI3 == 'Overweight')),
         Overweight_Week24_estimate = map_dbl(Overweight_Week24_Placebo_vs_Intervention, ~ pull(., estimate)),
         Overweight_Week24_pvalue = map_dbl(Overweight_Week24_Placebo_vs_Intervention, ~ pull(., p.value)),
         
         # BMI Obesity: Week-24 Placebo vs Intervention (Estimate and p-value)
         Obesity_Week24_Placebo_vs_Intervention = map(emmeans_week24, ~ contrast(., "pairwise") %>% 
                                                        summary() %>% 
                                                        filter(Baseline_BMI3 == 'Obesity')),
         Obesity_Week24_estimate = map_dbl(Obesity_Week24_Placebo_vs_Intervention, ~ pull(., estimate)),
         Obesity_Week24_pvalue = map_dbl(Obesity_Week24_Placebo_vs_Intervention, ~ pull(., p.value))
  )

# Step 5: Adjust p-values for multiple comparisons (BH method)
shannon_endline_bmi_final <- shannon_endline_bmi %>%
  mutate(Normal_Week24_pvalue_adj = p.adjust(Normal_Week24_pvalue, method = 'BH'),
         Overweight_Week24_pvalue_adj = p.adjust(Overweight_Week24_pvalue, method = 'BH'),
         Obesity_Week24_pvalue_adj = p.adjust(Obesity_Week24_pvalue, method = 'BH')) %>%
  mutate(
    Normal_Week24_significant = Normal_Week24_pvalue_adj < 0.05,
    Overweight_Week24_significant = Overweight_Week24_pvalue_adj < 0.05,
    Obesity_Week24_significant = Obesity_Week24_pvalue_adj < 0.05
  )


###################################################################################################
bmi_alphadiversity <- df_metadatapspair %>% dplyr::select(c("sampleid", "subject", "Treatment_Type","Time", 
                                                            "Baseline_BMI3", "diversity_shannon", "observed", "pd", "evenness_pielou"))

bmi_alphadiversity$Baseline_BMI3 <- factor(bmi_alphadiversity$Baseline_BMI3, levels = c("Normal", "Overweight", "Obesity"))  

##Linear model on alpha diversity
lmm_observedalpha <- bmi_alphadiversity %>%
  nest() %>%
  mutate(model = map(data, ~ glmmTMB(observed ~ Treatment_Type * Time * Baseline_BMI3 + (1|subject), 
                                     family = gaussian, data = .x)))

# Step 4: Post-hoc Comparisons using `emmeans`
bmi_observed_base24 <- lmm_observedalpha %>%
  mutate(emmeans_intervention = map(model, ~ emmeans(.x, ~ Time | Baseline_BMI3 | Treatment_Type, at = list(Time = c("Baseline", "24-wks")))),
         
         # BMI Normal: Intervention - Baseline vs Endline (Estimate and P-value)
         Normal_Intervention_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>% 
                                              summary() %>% 
                                              filter(Baseline_BMI3 == 'Normal', Treatment_Type == 'Intervention')),
         Normal_Intervention_estimate = map_dbl(Normal_Intervention_contrast, ~ pull(., estimate)),
         Normal_Intervention_pvalue = map_dbl(Normal_Intervention_contrast, ~ pull(., p.value)),
         
         # BMI Normal: Placebo - Baseline vs Endline (Estimate and P-value)
         Normal_Placebo_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>% 
                                         summary() %>% 
                                         filter(Baseline_BMI3 == 'Normal', Treatment_Type == 'Placebo')),
         Normal_Placebo_estimate = map_dbl(Normal_Placebo_contrast, ~ pull(., estimate)),
         Normal_Placebo_pvalue = map_dbl(Normal_Placebo_contrast, ~ pull(., p.value)),
         
         # BMI Overweight: Intervention - Baseline vs Endline (Estimate and P-value)
         Overweight_Intervention_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>% 
                                                  summary() %>% 
                                                  filter(Baseline_BMI3 == 'Overweight', Treatment_Type == 'Intervention')),
         Overweight_Intervention_estimate = map_dbl(Overweight_Intervention_contrast, ~ pull(., estimate)),
         Overweight_Intervention_pvalue = map_dbl(Overweight_Intervention_contrast, ~ pull(., p.value)),
         
         # BMI Overweight: Placebo - Baseline vs Endline (Estimate and P-value)
         Overweight_Placebo_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>% 
                                             summary() %>% 
                                             filter(Baseline_BMI3 == 'Overweight', Treatment_Type == 'Placebo')),
         Overweight_Placebo_estimate = map_dbl(Overweight_Placebo_contrast, ~ pull(., estimate)),
         Overweight_Placebo_pvalue = map_dbl(Overweight_Placebo_contrast, ~ pull(., p.value)),
         
         # BMI Obesity: Intervention - Baseline vs Endline (Estimate and P-value)
         Obesity_Intervention_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>% 
                                               summary() %>% 
                                               filter(Baseline_BMI3 == 'Obesity', Treatment_Type == 'Intervention')),
         Obesity_Intervention_estimate = map_dbl(Obesity_Intervention_contrast, ~ pull(., estimate)),
         Obesity_Intervention_pvalue = map_dbl(Obesity_Intervention_contrast, ~ pull(., p.value)),
         
         # BMI Obesity: Placebo - Baseline vs Endline (Estimate and P-value)
         Obesity_Placebo_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>% 
                                          summary() %>% 
                                          filter(Baseline_BMI3 == 'Obesity', Treatment_Type == 'Placebo')),
         Obesity_Placebo_estimate = map_dbl(Obesity_Placebo_contrast, ~ pull(., estimate)),
         Obesity_Placebo_pvalue = map_dbl(Obesity_Placebo_contrast, ~ pull(., p.value))
  ) %>%
  # Adjust p-values using the BH method
  mutate(across(ends_with("pvalue"), ~ p.adjust(., method = "BH"), .names = "adj_{col}"))

##no signficant effect due to baseline bmi group
View(bmi_observed_base24)

##Step 4: Post-hoc Comparisons for Week-24 Placebo vs Intervention in Each BMI Category
observed_endline_bmi <- lmm_observedalpha %>%
  mutate(emmeans_week24 = map(model, ~ emmeans(.x, ~ Treatment_Type | Baseline_BMI3 | Time, 
                                               at = list(Time = "24-wks"))),
         
         # BMI Normal: Week-24 Placebo vs Intervention (Estimate and p-value)
         Normal_Week24_Placebo_vs_Intervention = map(emmeans_week24, ~ contrast(., "pairwise") %>% 
                                                       summary() %>% 
                                                       filter(Baseline_BMI3 == 'Normal')),
         Normal_Week24_estimate = map_dbl(Normal_Week24_Placebo_vs_Intervention, ~ pull(., estimate)),
         Normal_Week24_pvalue = map_dbl(Normal_Week24_Placebo_vs_Intervention, ~ pull(., p.value)),
         
         # BMI Overweight: Week-24 Placebo vs Intervention (Estimate and p-value)
         Overweight_Week24_Placebo_vs_Intervention = map(emmeans_week24, ~ contrast(., "pairwise") %>% 
                                                           summary() %>% 
                                                           filter(Baseline_BMI3 == 'Overweight')),
         Overweight_Week24_estimate = map_dbl(Overweight_Week24_Placebo_vs_Intervention, ~ pull(., estimate)),
         Overweight_Week24_pvalue = map_dbl(Overweight_Week24_Placebo_vs_Intervention, ~ pull(., p.value)),
         
         # BMI Obesity: Week-24 Placebo vs Intervention (Estimate and p-value)
         Obesity_Week24_Placebo_vs_Intervention = map(emmeans_week24, ~ contrast(., "pairwise") %>% 
                                                        summary() %>% 
                                                        filter(Baseline_BMI3 == 'Obesity')),
         Obesity_Week24_estimate = map_dbl(Obesity_Week24_Placebo_vs_Intervention, ~ pull(., estimate)),
         Obesity_Week24_pvalue = map_dbl(Obesity_Week24_Placebo_vs_Intervention, ~ pull(., p.value))
  )

# Step 5: Adjust p-values for multiple comparisons (BH method)
observed_endline_bmi_final <- observed_endline_bmi %>%
  mutate(Normal_Week24_pvalue_adj = p.adjust(Normal_Week24_pvalue, method = 'BH'),
         Overweight_Week24_pvalue_adj = p.adjust(Overweight_Week24_pvalue, method = 'BH'),
         Obesity_Week24_pvalue_adj = p.adjust(Obesity_Week24_pvalue, method = 'BH')) %>%
  mutate(
    Normal_Week24_significant = Normal_Week24_pvalue_adj < 0.05,
    Overweight_Week24_significant = Overweight_Week24_pvalue_adj < 0.05,
    Obesity_Week24_significant = Obesity_Week24_pvalue_adj < 0.05
  )


#############################################################################################################
bmi_alphadiversity <- df_metadatapspair %>% dplyr::select(c("sampleid", "subject", "Treatment_Type","Time", 
                                                            "Baseline_BMI3", "diversity_shannon", "observed", "pd", "evenness_pielou"))

bmi_alphadiversity$Baseline_BMI3 <- factor(bmi_alphadiversity$Baseline_BMI3, levels = c("Normal", "Overweight", "Obesity"))  

##Linear model on alpha diversity
lmm_pd_alpha <- bmi_alphadiversity %>%
  nest() %>%
  mutate(model = map(data, ~ glmmTMB(pd ~ Treatment_Type * Time * Baseline_BMI3 + (1|subject), 
                                     family = gaussian, data = .x)))

# Step 4: Post-hoc Comparisons using `emmeans`
bmi_pd_base24 <- lmm_pd_alpha %>%
  mutate(emmeans_intervention = map(model, ~ emmeans(.x, ~ Time | Baseline_BMI3 | Treatment_Type, at = list(Time = c("Baseline", "24-wks")))),
         
         # BMI Normal: Intervention - Baseline vs Endline (Estimate and P-value)
         Normal_Intervention_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>% 
                                              summary() %>% 
                                              filter(Baseline_BMI3 == 'Normal', Treatment_Type == 'Intervention')),
         Normal_Intervention_estimate = map_dbl(Normal_Intervention_contrast, ~ pull(., estimate)),
         Normal_Intervention_pvalue = map_dbl(Normal_Intervention_contrast, ~ pull(., p.value)),
         
         # BMI Normal: Placebo - Baseline vs Endline (Estimate and P-value)
         Normal_Placebo_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>% 
                                         summary() %>% 
                                         filter(Baseline_BMI3 == 'Normal', Treatment_Type == 'Placebo')),
         Normal_Placebo_estimate = map_dbl(Normal_Placebo_contrast, ~ pull(., estimate)),
         Normal_Placebo_pvalue = map_dbl(Normal_Placebo_contrast, ~ pull(., p.value)),
         
         # BMI Overweight: Intervention - Baseline vs Endline (Estimate and P-value)
         Overweight_Intervention_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>% 
                                                  summary() %>% 
                                                  filter(Baseline_BMI3 == 'Overweight', Treatment_Type == 'Intervention')),
         Overweight_Intervention_estimate = map_dbl(Overweight_Intervention_contrast, ~ pull(., estimate)),
         Overweight_Intervention_pvalue = map_dbl(Overweight_Intervention_contrast, ~ pull(., p.value)),
         
         # BMI Overweight: Placebo - Baseline vs Endline (Estimate and P-value)
         Overweight_Placebo_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>% 
                                             summary() %>% 
                                             filter(Baseline_BMI3 == 'Overweight', Treatment_Type == 'Placebo')),
         Overweight_Placebo_estimate = map_dbl(Overweight_Placebo_contrast, ~ pull(., estimate)),
         Overweight_Placebo_pvalue = map_dbl(Overweight_Placebo_contrast, ~ pull(., p.value)),
         
         # BMI Obesity: Intervention - Baseline vs Endline (Estimate and P-value)
         Obesity_Intervention_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>% 
                                               summary() %>% 
                                               filter(Baseline_BMI3 == 'Obesity', Treatment_Type == 'Intervention')),
         Obesity_Intervention_estimate = map_dbl(Obesity_Intervention_contrast, ~ pull(., estimate)),
         Obesity_Intervention_pvalue = map_dbl(Obesity_Intervention_contrast, ~ pull(., p.value)),
         
         # BMI Obesity: Placebo - Baseline vs Endline (Estimate and P-value)
         Obesity_Placebo_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>% 
                                          summary() %>% 
                                          filter(Baseline_BMI3 == 'Obesity', Treatment_Type == 'Placebo')),
         Obesity_Placebo_estimate = map_dbl(Obesity_Placebo_contrast, ~ pull(., estimate)),
         Obesity_Placebo_pvalue = map_dbl(Obesity_Placebo_contrast, ~ pull(., p.value))
  ) %>%
  # Adjust p-values using the BH method
  mutate(across(ends_with("pvalue"), ~ p.adjust(., method = "BH"), .names = "adj_{col}"))

##no signficant effect due to baseline bmi group
View(bmi_pd_base24)

# Step 5: Adjust p-values for multiple comparisons (BH method)
observed_endline_bmi_final <- observed_endline_bmi %>%
  mutate(Normal_Week24_pvalue_adj = p.adjust(Normal_Week24_pvalue, method = 'BH'),
         Overweight_Week24_pvalue_adj = p.adjust(Overweight_Week24_pvalue, method = 'BH'),
         Obesity_Week24_pvalue_adj = p.adjust(Obesity_Week24_pvalue, method = 'BH')) %>%
  mutate(
    Normal_Week24_significant = Normal_Week24_pvalue_adj < 0.05,
    Overweight_Week24_significant = Overweight_Week24_pvalue_adj < 0.05,
    Obesity_Week24_significant = Obesity_Week24_pvalue_adj < 0.05
  )



###############################WEEK24 (Placebo vs Intervention)
##Step 4: Post-hoc Comparisons for Week-24 Placebo vs Intervention in Each BMI Category
observed_endline_bmi <- lmm_observedalpha %>%
  mutate(emmeans_week24 = map(model, ~ emmeans(.x, ~ Treatment_Type | Baseline_BMI3 | Time, 
                                               at = list(Time = "24-wks"))),
         
         # BMI Normal: Week-24 Placebo vs Intervention (Estimate and p-value)
         Normal_Week24_Placebo_vs_Intervention = map(emmeans_week24, ~ contrast(., "pairwise") %>% 
                                                       summary() %>% 
                                                       filter(Baseline_BMI3 == 'Normal')),
         Normal_Week24_estimate = map_dbl(Normal_Week24_Placebo_vs_Intervention, ~ pull(., estimate)),
         Normal_Week24_pvalue = map_dbl(Normal_Week24_Placebo_vs_Intervention, ~ pull(., p.value)),
         
         # BMI Overweight: Week-24 Placebo vs Intervention (Estimate and p-value)
         Overweight_Week24_Placebo_vs_Intervention = map(emmeans_week24, ~ contrast(., "pairwise") %>% 
                                                           summary() %>% 
                                                           filter(Baseline_BMI3 == 'Overweight')),
         Overweight_Week24_estimate = map_dbl(Overweight_Week24_Placebo_vs_Intervention, ~ pull(., estimate)),
         Overweight_Week24_pvalue = map_dbl(Overweight_Week24_Placebo_vs_Intervention, ~ pull(., p.value)),
         
         # BMI Obesity: Week-24 Placebo vs Intervention (Estimate and p-value)
         Obesity_Week24_Placebo_vs_Intervention = map(emmeans_week24, ~ contrast(., "pairwise") %>% 
                                                        summary() %>% 
                                                        filter(Baseline_BMI3 == 'Obesity')),
         Obesity_Week24_estimate = map_dbl(Obesity_Week24_Placebo_vs_Intervention, ~ pull(., estimate)),
         Obesity_Week24_pvalue = map_dbl(Obesity_Week24_Placebo_vs_Intervention, ~ pull(., p.value))
  )

# Step 5: Adjust p-values for multiple comparisons (BH method)
observed_endline_bmi_final <- observed_endline_bmi %>%
  mutate(Normal_Week24_pvalue_adj = p.adjust(Normal_Week24_pvalue, method = 'BH'),
         Overweight_Week24_pvalue_adj = p.adjust(Overweight_Week24_pvalue, method = 'BH'),
         Obesity_Week24_pvalue_adj = p.adjust(Obesity_Week24_pvalue, method = 'BH')) %>%
  mutate(
    Normal_Week24_significant = Normal_Week24_pvalue_adj < 0.05,
    Overweight_Week24_significant = Overweight_Week24_pvalue_adj < 0.05,
    Obesity_Week24_significant = Obesity_Week24_pvalue_adj < 0.05
  )

###Beta diversity
pspair

##Gloom at the genus level
pspairGenus <- tax_glom(pspair, taxrank = "Genus")

##relative abundance
pspairGenusComp <- microbiome::transform(pspairGenus, transform = "compositional")
df_pspairGenusComp <- data.frame(sample_data(pspairGenusComp))
write.csv(df_pspairGenusComp, "df_pspairGenusComp.csv")

df_pspairGenusComp_1 <- read_csv("df_pspairGenusComp.csv")

df_pspairGenusComp_2 <- df_pspairGenusComp_1 %>% dplyr::select(-...1)

##Factorize the variables
df_pspairGenusComp_2$Time <- factor(df_pspairGenusComp_2$Time, levels = c("Baseline", "24-wks"))
df_pspairGenusComp_2$BMI_3_class <- factor(df_pspairGenusComp_2$BMI_3_class, levels = c("Normal", "Overweight", "Obesity"))
df_pspairGenusComp_2$Treatment_Type <- factor(df_pspairGenusComp_2, levels = c("Intervention", "Placebo"))

##Bray-Curtis dissimilarity
pspairGenusCompbray <- phyloseq::distance(pspairGenusComp, method = "bray")

bray_genus <- adonis2(pspairGenusCompbray~Treatment_Type*Time*BMI_3_class, data = df_pspairGenusComp_2, strata = df_pspairGenusComp$subject, by= "terms")

##pspair
##Bray-Curtis dissimilarity
pspair_asv_bray <- phyloseq::distance(pspair, method = "bray")

bray_asv <- adonis2(pspair_asv_bray~Treatment_Type*Time*BMI_3_class, data = df_pspairGenusComp_2, strata = df_pspairGenusComp$subject, by= "terms")

# Perform pairwise comparisons
pairwise_time <- pairwise.adonis(pspair_asv_bray ~ Time, data = df_pspairGenusComp_2, strata = df_pspairGenusComp$subject)
print(pairwise_results)


##Jaccard distance
pspair_asv_jaccard <- phyloseq::distance(pspair, method = "jaccard", binary = TRUE)

jaccard_asv <- adonis2(pspair_asv_jaccard~BMI_3_class*Treatment_Type*Time, data = df_pspairGenusComp_2, strata = df_pspairGenusComp$subject, by= "terms")

##Weighted UniFrac
pspair_asv_wunifrac <- phyloseq::distance(pspair, method = "wunifrac")

wunifracasv_bmi <- adonis2(pspair_asv_wunifrac~BMI_3_class*Treatment_Type*Time, data = df_pspairGenusComp_2, strata = df_pspairGenusComp$subject, by= "terms")

##Unweighted UniFrac
pspair_asv_unifrac <- phyloseq::distance(pspair, method = "unifrac")

unifrac_asvbmi <- adonis2(pspair_asv_unifrac~BMI_3_class*Treatment_Type*Time, data = df_pspairGenusComp_2, strata = df_pspairGenusComp$subject, by= "terms")
df_pspairGenusComp_2$subject <- factor(df_pspairGenusComp_2$subject)
df_pspairGenusComp_2 <- data.frame(df_pspairGenusComp_2)

# Perform pairwise comparisons
pairwise_bmitreat <- pairwise.adonis2(pspair_asv_unifrac ~ BMI_3_class*Treatment_Type*Time, 
                                    data = df_pspairGenusComp_2, 
                                    strata = 'subject')


#######Effect of intervention/BMI on response
# Centered log-ratio transformation
fulldata_clr <- fulldata %>%
  pivot_longer(cols = starts_with('ASV'), names_to = 'ASV', values_to = 'value') %>%
  group_by(id) %>% 
  ungroup()

fulldata_clrselect <- fulldata_clr %>% dplyr::select(c("id", "subject", "Treatment_Type","Time",
                                                       "BMI_3_class","value",matches("ASV")))

fulldata_clrselect$BMI_3_class <- factor(fulldata_clrselect$BMI_3_class, levels = c("Normal", "Overweight", "Obesity"))

fulldata_clrselect$Treatment_Type <- factor(fulldata_clrselect$Treatment_Type, levels = c("Intervention", "Placebo"))

# We want to pivot the data to long format and include treatment, time, and enterotype
effect_intervention_bmi <- fulldata_clrselect %>%
  group_by(ASV) %>%
  nest() %>%
  mutate(model = map(data, ~ glmmTMB(value ~ BMI_3_class*Treatment_Type * Time + (1|subject), 
                                     family = gaussian, data = .x)))

# Create contrast matrix with BMI as covariate
contrast_BMI_ALL <- effect_intervention_bmi %>%
  mutate(
    # Generate emmeans for the interaction of Time, Treatment_Type, and BMI_3_class
    emmeans_intervention = map(model, ~ emmeans(.x, ~ Time | Treatment_Type * BMI_3_class)),
    
    ### Contrasts for each BMI group ###
    
    # Normal BMI group - Intervention: Baseline vs Endline
    NormalBMI_Intervention_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>%
                                            summary() %>%
                                            filter(Treatment_Type == 'Intervention', BMI_3_class == 'Normal')),
    
    NormalBMI_Intervention_estimate = map_dbl(NormalBMI_Intervention_contrast, ~ ifelse(nrow(.) > 0, pull(., estimate), NA)),
    NormalBMI_Intervention_pvalue = map_dbl(NormalBMI_Intervention_contrast, ~ ifelse(nrow(.) > 0, pull(., p.value), NA)),
    
    # Normal BMI group - Placebo: Baseline vs Endline
    NormalBMI_Placebo_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>%
                                       summary() %>%
                                       filter(Treatment_Type == 'Placebo', BMI_3_class == 'Normal')),
    
    NormalBMI_Placebo_estimate = map_dbl(NormalBMI_Placebo_contrast, ~ ifelse(nrow(.) > 0, pull(., estimate), NA)),
    NormalBMI_Placebo_pvalue = map_dbl(NormalBMI_Placebo_contrast, ~ ifelse(nrow(.) > 0, pull(., p.value), NA)),
    
    # Normal BMI group - Treatment: Intervention vs Placebo at Endline
    NormalBMI_Endline_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Treatment_Type") %>%
                                       summary() %>%
                                       filter(Time == '24-wks', BMI_3_class == 'Normal')),
    
    NormalBMI_Endline_estimate = map_dbl(NormalBMI_Endline_contrast, ~ ifelse(nrow(.) > 0, pull(., estimate), NA)),
    NormalBMI_Endline_pvalue = map_dbl(NormalBMI_Endline_contrast, ~ ifelse(nrow(.) > 0, pull(., p.value), NA)),
    
    # Overweight BMI group - Intervention: Baseline vs Endline
    OverweightBMI_Intervention_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>%
                                                summary() %>%
                                                filter(Treatment_Type == 'Intervention', BMI_3_class == 'Overweight')),
    
    OverweightBMI_Intervention_estimate = map_dbl(OverweightBMI_Intervention_contrast, ~ ifelse(nrow(.) > 0, pull(., estimate), NA)),
    OverweightBMI_Intervention_pvalue = map_dbl(OverweightBMI_Intervention_contrast, ~ ifelse(nrow(.) > 0, pull(., p.value), NA)),
    
    # Overweight BMI group - Placebo: Baseline vs Endline
    OverweightBMI_Placebo_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>%
                                           summary() %>%
                                           filter(Treatment_Type == 'Placebo', BMI_3_class == 'Overweight')),
    
    OverweightBMI_Placebo_estimate = map_dbl(OverweightBMI_Placebo_contrast, ~ ifelse(nrow(.) > 0, pull(., estimate), NA)),
    OverweightBMI_Placebo_pvalue = map_dbl(OverweightBMI_Placebo_contrast, ~ ifelse(nrow(.) > 0, pull(., p.value), NA)),
    
    # Overweight BMI group - Treatment: Intervention vs Placebo at Endline
    OverweightBMI_Endline_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Treatment_Type") %>%
                                           summary() %>%
                                           filter(Time == '24-wks', BMI_3_class == 'Overweight')),
    
    OverweightBMI_Endline_estimate = map_dbl(OverweightBMI_Endline_contrast, ~ ifelse(nrow(.) > 0, pull(., estimate), NA)),
    OverweightBMI_Endline_pvalue = map_dbl(OverweightBMI_Endline_contrast, ~ ifelse(nrow(.) > 0, pull(., p.value), NA)),
    
    # Obese BMI group - Intervention: Baseline vs Endline
    ObeseBMI_Intervention_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>%
                                           summary() %>%
                                           filter(Treatment_Type == 'Intervention', BMI_3_class == 'Obesity')),
    
    ObeseBMI_Intervention_estimate = map_dbl(ObeseBMI_Intervention_contrast, ~ ifelse(nrow(.) > 0, pull(., estimate), NA)),
    ObeseBMI_Intervention_pvalue = map_dbl(ObeseBMI_Intervention_contrast, ~ ifelse(nrow(.) > 0, pull(., p.value), NA)),
    
    # Obese BMI group - Placebo: Baseline vs Endline
    ObeseBMI_Placebo_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Time") %>%
                                      summary() %>%
                                      filter(Treatment_Type == 'Placebo', BMI_3_class == 'Obesity')),
    
    ObeseBMI_Placebo_estimate = map_dbl(ObeseBMI_Placebo_contrast, ~ ifelse(nrow(.) > 0, pull(., estimate), NA)),
    ObeseBMI_Placebo_pvalue = map_dbl(ObeseBMI_Placebo_contrast, ~ ifelse(nrow(.) > 0, pull(., p.value), NA)),
    
    # Obese BMI group - Treatment: Intervention vs Placebo at Endline
    ObeseBMI_Endline_contrast = map(emmeans_intervention, ~ contrast(., interaction = "pairwise", simple = "Treatment_Type") %>%
                                      summary() %>%
                                      filter(Time == '24-wks', BMI_3_class == 'Obesity')),
    
    ObeseBMI_Endline_estimate = map_dbl(ObeseBMI_Endline_contrast, ~ ifelse(nrow(.) > 0, pull(., estimate), NA)),
    ObeseBMI_Endline_pvalue = map_dbl(ObeseBMI_Endline_contrast, ~ ifelse(nrow(.) > 0, pull(., p.value), NA)))



effect_intervention_bmi <- fulldata_long %>%
  group_by(ASV) %>%
  nest() %>%
  mutate(model = map(data, ~ glmmTMB(value ~ BMI_3_class * Treatment_Type * Time + (1|subject), 
                                     family = gaussian, data = .x)))


emmeans_resultsBMI <- effect_intervention_bmi %>%
  mutate(emmeans_intervention = map(model, ~ emmeans(.x, ~ Time | Treatment_Type * BMI_3_class)))

# Extract EMMs for each combination
emmeans_dataBMI <- emmeans_resultsBMI %>%
  mutate(emmeans_df = map(emmeans_intervention, as.data.frame)) %>%
  dplyr::select(ASV, emmeans_df) %>%
  unnest(emmeans_df)

# View the extracted EMMs
print(emmeans_dataBMI)

emmeans_differences <- emmeans_resultsBMI %>%
  mutate(contrast_baseline_endline = map(emmeans_intervention, ~ contrast(.x, interaction = "pairwise", simple = "Time") %>%
                                           summary() %>%
                                           filter(contrast == "Baseline - 24-wks")))

# Extract the differences
emmeans_diff_data <- emmeans_differences %>%
  mutate(contrast_df = map(contrast_baseline_endline, as.data.frame)) %>%
  select(ASV, contrast_df) %>%
  unnest(contrast_df)

# View the differences
print(emmeans_diff_data)



##Semi join to filter only those useful
emmeans_dataBMI_2_taxa <- tax_filtercomp1Gen1 %>% semi_join(emmeans_dataBMI, by = "ASV")

emmeans_dataBMI_genus <- full_join(emmeans_dataBMI_2_taxa, emmeans_dataBMI, by = "ASV")
write.csv(emmeans_dataBMI_genus, "BMIcategory_final/emmeans_dataBMI_genus.csv")


# Adjust p-values using BH method
contrast_BMI_ALL_1 <- contrast_BMI_ALL %>%
  mutate(
    ObeseBMI_Intervention_adj = p.adjust(ObeseBMI_Intervention_pvalue, method = "BH"),
    ObeseBMI_Placebo_adj = p.adjust(ObeseBMI_Placebo_pvalue, method = "BH"),
    ObeseBMI_Endline_adj = p.adjust(ObeseBMI_Endline_pvalue, method = "BH"),
    OverweightBMI_Endline_adj = p.adjust(OverweightBMI_Endline_pvalue, method = "BH"),
    OverweightBMI_Placebo_adj = p.adjust(OverweightBMI_Placebo_pvalue, method = "BH"),
    OverweightBMI_Intervention_adj = p.adjust(OverweightBMI_Intervention_pvalue, method = "BH"),
    NormalBMI_Endline_adj = p.adjust(NormalBMI_Endline_pvalue, method = "BH"),
    NormalBMI_Placebo_adj = p.adjust(NormalBMI_Placebo_pvalue, method = "BH"),
    NormalBMI_Intervention_adj = p.adjust(NormalBMI_Intervention_pvalue, method = "BH"))

contrast_BMI_ALL_2 <- contrast_BMI_ALL_1 %>% dplyr::select(-c("data", "OverweightBMI_Placebo_contrast", "model", "NormalBMI_Intervention_contrast", "NormalBMI_Placebo_contrast",
                                                                          "NormalBMI_Endline_contrast", "OverweightBMI_Intervention_contrast", "OverweightBMI_Placebo_contrast",
                                                                  "OverweightBMI_Endline_contrast", "ObeseBMI_Intervention_contrast", "ObeseBMI_Placebo_contrast",
                                                                  "ObeseBMI_Endline_contrast", "emmeans_intervention"))

View(contrast_BMI_ALL_2)

tax_filtercomp1Gen <- data.frame(tax_table(pspairgenus_clr))

tax_filtercomp1Gen1 <- tax_filtercomp1Gen %>% dplyr::select(-Species) %>% rownames_to_column("ASV")

##Semi join to filter only those useful
contrast_BMI_ALL_2_taxa <- tax_filtercomp1Gen1 %>% semi_join(contrast_BMI_ALL_2, by = "ASV")

contrast_BMI_ALL_genus <- full_join(contrast_BMI_ALL_2_taxa, contrast_BMI_ALL_2, by = "ASV")

write.csv(contrast_BMI_ALL_genus, "BMI_category/contrast_BMI_ALL_genus.csv")

#############################################################################################################
fulldata_clrnewselect <- fulldata_clr %>% dplyr::select(c("id", "subject", "Treatment_Type","Time", "baseline_shannon", "Baseline_enterotype_BIC",
                                                       "BMI_3_class","value",matches("ASV")))

fulldata_clrnewselect$BMI_3_class <- factor(fulldata_clrnewselect$BMI_3_class, levels = c("Normal", "Overweight", "Obesity"))

fulldata_clrnewselect$Treatment_Type <- factor(fulldata_clrnewselect$Treatment_Type, levels = c("Intervention", "Placebo"))

fulldata_clrnewselect$Baseline_enterotype_BIC ##Levels: Enterotype one Enterotype two

fulldata_clrnewselect$baseline_shannon <- factor(fulldata_clrnewselect$baseline_shannon, levels = c("Low", "Medium", "High"))
# We want to pivot the data to long format and include treatment, time, and enterotype
effect_intervention_bmi <- fulldata_clrnewselect %>%
  group_by(ASV) %>%
  nest() %>%
  mutate(model = map(data, ~ glmmTMB(value ~ BMI_3_class*Treatment_Type * Time + (1|subject), 
                                     family = gaussian, data = .x)))

# Load necessary libraries
library(lme4)
library(interactions)
library(ggplot2)
model_all <- lmer(value ~ Treatment_Type * Time * Baseline_enterotype_BIC * BMI_3_class + (1|subject), data = fulldata_clrnewselect)
summary(model)

# Fit a linear model to check VIF for main effects and lower-order terms
lm_model <- lm(value ~ Treatment_Type + Time + Baseline_enterotype_BIC + BMI_3_class + baseline_shannon, data = fulldata_clrnewselect)
vif_values <- vif(lm_model, type = "predictor")
print(vif_values)


# Load necessary libraries
library(lme4)
library(interactions)
library(ggplot2)

# Convert categorical variables to factors
fulldata_clrnewselect$Treatment_Type <- as.factor(fulldata_clrnewselect$Treatment_Type)
fulldata_clrnewselect$Time <- as.factor(fulldata_clrnewselect$Time)
fulldata_clrnewselect$baseline_shannon <- as.factor(fulldata_clrnewselect$baseline_shannon)
fulldata_clrnewselect$Baseline_enterotype_BIC <- as.factor(fulldata_clrnewselect$Baseline_enterotype_BIC)
fulldata_clrnewselect$BMI_3_class <- as.factor(fulldata_clrnewselect$BMI_3_class)
fulldata_clrnewselect$age_group <- as.factor(fulldata_clrnewselect$age_group)

# Function to fit the simplified model for each ASV
fit_simplified_model_for_asv <- function(asv) {
  # Subset the data for the given ASV
  asv_data <- subset(fulldata_clrnewselect, ASV == asv)
  
  # Fit the simplified linear mixed-effects model
  model <- lmer(value ~ Treatment_Type * Time * Baseline_enterotype_BIC * BMI_3_class + (1|subject), data = asv_data)
  
  # Return the model summary
  return(summary(model))
}

# Get a list of unique ASVs
unique_asvs <- unique(fulldata_clrnewselect$ASV)

# Fit the simplified model for each ASV and store the results
simplified_model_results <- lapply(unique_asvs, fit_simplified_model_for_asv)

# Example of how to access the results for the first ASV
print(simplified_model_results[])

# Apply ANOVA to each model in the list and store the results
anova_results_list <- lapply(simplified_model_results, anova)

# Example of how to access the ANOVA results for the first ASV
print(anova_results_list[])

# If you want to save all ANOVA results to a file
save(anova_results_list, file = "anova_results_list.RData")




# Visualization: Interaction plot for the first ASV
interact_plot(simplified_model_results[], pred = "BMI_3_class", modx = "Baseline_enterotype_BIC")

# Save the simplified model results to a file
save(simplified_model_results, file = "simplified_model_results.RData")

# Example of creating a heatmap for the first ASV
heatmap_data <- asv_data[, c("Treatment_Type", "Time", "Baseline_enterotype_BIC", "BMI_3_class", "age_group", "baseline_shannon", "value")]
heatmap_data <- reshape2::melt(heatmap_data, id.vars = c("Treatment_Type", "Time", "Baseline_enterotype_BIC", "BMI_3_class", "age_group", "baseline_shannon"))
ggplot(heatmap_data, aes(x = variable, y = value, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(title = "Heatmap of Genus Abundance", x = "Factors", y = "Abundance")

# Model selection and validation
# Example using AIC for model comparison
aic_values <- sapply(simplified_model_results, AIC)
print(aic_values)

# Cross-validation (example using caret package)
library(caret)
train_control <- trainControl(method = "cv", number = 10)
cv_model <- train(value ~ Treatment_Type * Time * Baseline_enterotype_BIC * BMI_3_class, data = asv_data, method = "lm", trControl = train_control)
print(cv_model)
