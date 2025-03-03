library(geepack)
#ps_study810_rare_f2

baseline_study810raref2 <- subset_samples(ps_study810_rare_f2, Time == "Baseline") %>% 
  sample_data() %>% data.frame()
class(baseline_study810raref2)

##Factorize
baseline_study810raref2 <- baseline_study810raref2 %>% mutate(APOE_status = factor(APOE_status), BMI_3_class = factor(BMI_3_class),
         Screening_Diagnosis = factor(Screening_Diagnosis), BMI_SD = factor(BMI_SD), Gender = factor(Gender), Age_quartile = factor(Age_quartile), 
         Family_dementia = factor(Family_dementia), Smoker = factor(Smoker), Coronary_heart_disease = factor(Coronary_heart_disease), 
         Hypertension_treatment = factor(Hypertension_treatment), Diabetes_mellitus_all = factor(Diabetes_mellitus_all),
         Antithrombotic_agents = factor(Antithrombotic_agents), Calcium_channel_blockers = factor(Calcium_channel_blockers), 
         Thyroid_therapy = factor(Thyroid_therapy), Drug_acid_related_disorders = factor(Drug_acid_related_disorders), 
         Hypercholesterolemia_drug = factor(Hypercholesterolemia_drug), Residence = factor(Residence), Enterotype_BIC = factor(Enterotype_BIC), 
         faithpd_tertile = factor(faithpd_tertile, levels = c("Low", "Medium", "High")), observed_tertile = factor(observed_tertile, levels = c("Low", "Medium", "High")), 
         shannon_tertile = factor(shannon_tertile, levels = c("Low", "Medium", "High")), evenness_tertile = factor(evenness_tertile, levels = c("Low", "Medium", "High")))

# List of variables to test
variables <- c("APOE_status", "BMI_3_class", "Screening_Diagnosis", "BMI_SD", "Gender", "Age_quartile", 
               "Family_dementia", "Smoker", "Coronary_heart_disease", "Hypertension_treatment", 
               "Diabetes_mellitus_all", "Antithrombotic_agents", "Calcium_channel_blockers", 
               "Thyroid_therapy", "Drug_acid_related_disorders", "Hypercholesterolemia_drug", 
               "Residence", "faithpd_tertile", "observed_tertile", "shannon_tertile", "evenness_tertile")

# Function to determine significance level
get_significance <- function(p_value) {
  if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}

# Loop through each variable and perform Fisher's exact test
results_fisher <- lapply(variables, function(var) {
  test <- fisher.test(table(baseline_study810raref2$Enterotype_BIC, baseline_study810raref2[[var]]))
  significance <- get_significance(test$p.value)
  return(list(variable = var, p.value = test$p.value, significance = significance))
})

# Convert results to a data frame for easy viewing
results_fisher_df <- do.call(rbind, lapply(results_fisher, as.data.frame))
print(results_fisher_df) ##only alpha diversity groups are the only ones that are associated with enterotypes

###compare BMI and Age with Enterotype
##BMI
bmi_enterotype <- kruskal.test(BMI ~ Enterotype_BIC, data = baseline_study810raref2) 

##Age
age_enterotype <- kruskal.test(Age ~ Enterotype_BIC, data = baseline_study810raref2) 
##Faith PD and Enterotype
table_enterofaithpd <- table(baseline_study810raref2$Enterotype_BIC, baseline_study810raref2$faithpd_tertile)

# Perform Fisher's exact test
fisher_result_faithpdentero <- fisher.test(table_enterofaithpd)

pairwise_fisher_test(table_enterofaithpd, p.adjust.method = "bonferroni")

# group1 group2     n          p      p.adj p.adj.signif
# Low    Medium    72 0.00000122 0.00000366 ****        
# Low    High      48 0.0000065  0.0000195  ****        
# Medium High      72 0.739      1          ns 

##Faith PD and Enterotype
table_entero_observed <- table(baseline_study810raref2$Enterotype_BIC, baseline_study810raref2$observed_tertile)

# Perform Fisher's exact test
fisher_result_observedentero <- fisher.test(table_entero_observed)

pairwise_fisher_test(table_entero_observed, p.adjust.method = "bonferroni")

# group1 group2     n             p        p.adj p.adj.signif
# Low    Medium    73 0.00000000493 0.0000000148 ****        
# Low    High      47 0.0000000365  0.00000011   ****        
# Medium High      72 0.485         1            ns 


##Shannon diversity and Enterotype
table_entero_shannon <- table(baseline_study810raref2$Enterotype_BIC, baseline_study810raref2$shannon_tertile)

# Perform Fisher's exact test
fisher_result_shannonentero <- fisher.test(table_entero_shannon)

pairwise_fisher_test(table_entero_shannon, p.adjust.method = "bonferroni")

# group1 group2     n        p    p.adj p.adj.signif
# Low    Medium    72 1.08e-11 3.24e-11 ****        
# Low    High      48 3.58e-11 1.07e-10 ****        
# Medium High      72 2.55e- 1 7.65e- 1 ns


##Pielou's evenness and Enterotype
table_entero_evenness <- table(baseline_study810raref2$Enterotype_BIC, baseline_study810raref2$evenness_tertile)

# Perform Fisher's exact test
fisher_result_evennessentero <- fisher.test(table_entero_evenness)

pairwise_fisher_test(table_entero_evenness, p.adjust.method = "bonferroni")

# group1 group2     n             p        p.adj p.adj.signif
# Low    Medium    72 0.0000000213  0.0000000639 ****        
# Low    High      48 0.00000000357 0.0000000107 ****        
# Medium High      72 0.149         0.447        ns

##Effect of treatment and time on enterotype distribution
df_paired_metadata <- subset_samples(pspair) %>%
  sample_data() %>% data.frame()

###Time as a factor
df_paired_metadata$Time <- factor(df_paired_metadata$Time, levels = c("Baseline", "24-wks"))
df_paired_metadata$Treatment_Type <- factor(df_paired_metadata$Treatment_Type)
df_paired_metadata$subject <- as.factor(df_paired_metadata$subject)
# Convert Enterotype_BIC to numeric
df_paired_metadata$Enterotype_BIC <- as.numeric(df_paired_metadata$Enterotype_BIC) - 1

# Fit the GEE model
enterotimetreatment_model <- geeglm(Enterotype_BIC ~ Treatment_Type * Time, family = binomial(link = "logit"), 
                                    data = df_paired_metadata, id = subject, corstr = "exchangeable")
summary(enterotimetreatment_model)

df_study810_rare_f2 <- ps_study810_rare_f2 %>% sample_data() %>% data.frame()
class(df_study810_rare_f2) ##data.frame
df_study810_rare_f2$Enterotype_BIC <- factor(df_study810_rare_f2$Enterotype_BIC, levels = c("Enterotype one",
                                                                                            "Enterotype two"))

df_study810_rare_f2$Time <- factor(df_study810_rare_f2$Time, levels = c("Baseline", "12-wks", "24-wks"))


##plot the proprotion of enterotype at each timepoint
# Plot with integrated proportion calculation and text labels
enterotype_prop <- df_study810_rare_f2 %>%
  group_by(Time, Enterotype_BIC) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(Time) %>%
  mutate(Proportion = Count / sum(Count)) %>% 
  ggplot(aes(x = Time, fill = Enterotype_BIC)) +  
  geom_bar(position = "fill", aes(y = after_stat(prop)), stat = "count") +  # Now using after_stat()
  labs(y = "", x = "") +
  theme_bw() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.title = element_blank(),
        legend.margin = margin(t=0,r=0,b=0,l=-6)) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02")) +
  geom_text(aes(label = scales::percent(Proportion, accuracy = 1), y = Proportion),
            position = position_stack(vjust = 0.5), stat = "identity", show.legend = FALSE)
ggsave(plot = enterotype_prop,"enterotype_final/enterotypeprop.pdf",dpi = 600,width = 5,height = 4)

