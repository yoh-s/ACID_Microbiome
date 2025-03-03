library(phyloseq)
library(phyloseq.extended) ##FaithPD
library(microbiome)
library(rstatix)

##Calculate alpha diversity
ps_study810_rare_f

alpha_all <- microbiome::alpha(ps_study810_rare_f, index = "all")
sample_variables(alpha_all)
alpha_all1 <- alpha_all %>% select("observed", "diversity_shannon", "evenness_pielou")

##Change the rownames to sampleid
alpha_all1 <- rownames_to_column(alpha_all1, var = "sampleid")
##faith phylogenetic diversity
faith_phylo <- phylodiv(ps_study810_rare_f) #theta = 0 is faith phylogenetic diversity

##merge two alpha diversity metrics
alpha_diversity <- right_join(faith_phylo, alpha_all1, by = "sampleid")
sample_variables(alpha_diversity)
alpha_diversity <- alpha_diversity %>% select(-c("...1","Time.y"))
alpha_diversity_f <- alpha_diversity %>% select("Screening_code","pd","observed","diversity_shannon",
                                                "evenness_pielou")
alpha_diversity_f <- alpha_diversity_f %>% rename(sampleid = Screening_code)
write.csv(alpha_diversity_f, "alpha_diversity_f.csv")

#Merge to phyloseq object
ps_study810_rare_f1 <- ps_join(ps_study810_rare_f, alpha_diversity_f, by = "sampleid")

sample_variables(ps_study810_rare_f1) 
##   "...1"    "Time.x"  "Time.y"

alphadiversity_class <- read_csv("alphadiversity_classification.csv")

ps_study810_rare_f2 <- ps_join(ps_study810_rare_f1,alphadiversity_class,by = "sampleid")

sample_variables(ps_study810_rare_f2)
df_ps_study810_rare_f2 <- data.frame(sample_data(ps_study810_rare_f2))
df_ps_study810_rare_f3 <- df_ps_study810_rare_f2 %>% select(-c("...1", "Time.y", "Screening_code"))

##Rename Time.x
df_ps_study810_rare_f3 <- df_ps_study810_rare_f3 %>% rename(Time = Time.x) 

sample_data(ps_study810_rare_f2) <- df_ps_study810_rare_f3
sample_data(ps_study810_rare_f2) %>% as.data.frame() %>% View()

##Final phyloseq object
ps_study810_rare_f2

##use the data frame used to calculate alpha diversity
alpha_diversity
sample_variables(alpha_diversity)
alphadiversity_categoryk <- alpha_diversity %>% 
  select(c("sampleid", "Treatment_Type", "Time.x", "Paired_all", "subject", "BMI", "Age", "APOE_status","BMI_3_class","Screening_Diagnosis","BMI_SD","Gender",
           "Age_quartile","Family_dementia","Smoker","Coronary_heart_disease",
           "Hypertension_treatment", "Diabetes_mellitus_all", "Antithrombotic_agents",
           "Calcium_channel_blockers", "Thyroid_therapy","Drug_acid_related_disorders",
           "Hypercholesterolemia_drug","Residence", "Enterotype", "Baseline_enterotype",
           "pd","observed","diversity_shannon","evenness_pielou"))

##Rename Time.x
alphadiversity_categoryk <- alphadiversity_categoryk %>% rename(Time = Time.x) 

##Factorize the comparison variables
alphadiversity_categoryk <- alphadiversity_categoryk %>% mutate(APOE_status = factor(APOE_status), Time = factor(Time), Treatment_Type = factor(Treatment_Type),
                                                                BMI_3_class = factor(BMI_3_class),Screening_Diagnosis = factor(Screening_Diagnosis),
                                                                BMI_SD = factor(BMI_SD),Gender = factor(Gender),Age_quartile = factor(Age_quartile), Family_dementia = factor(Family_dementia),Smoker = factor(Smoker),
                                                                Coronary_heart_disease = factor(Coronary_heart_disease), Hypertension_treatment = factor(Hypertension_treatment),
                                                                Diabetes_mellitus_all = factor(Diabetes_mellitus_all), Antithrombotic_agents = factor(Antithrombotic_agents),
                                                                Calcium_channel_blockers = factor(Calcium_channel_blockers), Thyroid_therapy = factor(Thyroid_therapy),
                                                                Drug_acid_related_disorders = factor(Drug_acid_related_disorders), Hypercholesterolemia_drug = factor(Hypercholesterolemia_drug),
                                                                Residence = factor(Residence), BMI_SD = factor(BMI_SD),
                                                                Enterotype = factor(Enterotype), Baseline_enterotype = factor(Baseline_enterotype))

##kruskal wallis test
baseline_alphadiversity <- alphadiversity_categoryk %>% filter(Time == "Baseline")

##select variable for kruskal wallis
##kruskal variables (used the same variable than beta diversity)
krus_variables <- c("Screening_Diagnosis", "Residence", "Age_quartile", "Gender", 
                    "BMI_3_class", "Coronary_heart_disease", "Hypercholesterolemia_drug",
                    "Thyroid_therapy", "Family_dementia", "Hypertension_treatment", "Calcium_channel_blockers",
                    "Diabetes_mellitus_all", "Antithrombotic_agents", 
                    "Drug_acid_related_disorders", "BMI_SD", "Enterotype", "Baseline_enterotype")

# Initialize an empty list to store results
base_observed_kruskal <- list()

# Loop through each variable and perform the Kruskal-Wallis test on observed richness
for (var in krus_variables) {
  # Create a formula
  formula <- as.formula(paste("observed ~", var))
  
  # Perform Kruskal-Wallis test
  result <- baseline_alphadiversity %>% 
    kruskal_test(formula) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj")
  
  # Add the variable name to the result
  result <- result %>% mutate(variable = var)
  
  # Store the result in the list
  base_observed_kruskal[[var]] <- result
}

# Combine all results into one data frame
kruskal_baseobserved_df <- bind_rows(base_observed_kruskal)
write.csv(kruskal_baseobserved_df, "output_results/kruskal_baseobserved_df.csv")

# Print the combined results
View(kruskal_baseobserved_df)
##BMI_3_class
##pairwise comparison for BMI_3_class
bmi3observed_pair <- baseline_alphadiversity %>% 
  dunn_test(observed ~ BMI_3_class, p.adjust.method = "BH") 
View(bmi3observed_pair)

##BMI_SD
BMISD_observed_pair <- baseline_alphadiversity %>% 
  dunn_test(observed ~ BMI_SD, p.adjust.method = "BH") 
View(BMISD_observed_pair)

##Enterotype
Enterotype_observed_pair <- baseline_alphadiversity %>% 
  dunn_test(observed ~ Enterotype, p.adjust.method = "BH") 
View(Enterotype_observed_pair)


#############################
# Initialize an empty list to store results
base_pd_kruskal <- list()

# Loop through each variable and perform the Kruskal-Wallis test on observed richness
for (var in krus_variables) {
  # Create a formula
  formula <- as.formula(paste("pd ~", var))
  
  # Perform Kruskal-Wallis test
  result <- baseline_alphadiversity %>% 
    kruskal_test(formula) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj")
  
  # Add the variable name to the result
  result <- result %>% mutate(variable = var)
  
  # Store the result in the list
  base_pd_kruskal[[var]] <- result
}

# Combine all results into one data frame
kruskal_basepd_df <- bind_rows(base_pd_kruskal)

# Print the combined results
View(kruskal_basepd_df)
##BMI_3_class
##BMI_SD
##Enterotype

##BMI_3_class
##pairwise comparison for BMI_3_class
bmi3pd_pair <- baseline_alphadiversity %>% 
  dunn_test(pd ~ BMI_3_class, p.adjust.method = "BH") 
View(bmi3pd_pair)

##BMI_SD
BMISD_pd_pair <- baseline_alphadiversity %>% 
  dunn_test(pd ~ BMI_SD, p.adjust.method = "BH") 
View(BMISD_pd_pair)

##Enterotype
Enterotype_pd_pair <- baseline_alphadiversity %>% 
  dunn_test(pd ~ Enterotype, p.adjust.method = "BH") 
View(Enterotype_pd_pair)

##################################
# Initialize an empty list to store results
base_shannon_kruskal <- list()

# Loop through each variable and perform the Kruskal-Wallis test on shannon diversity
for (var in krus_variables) {
  # Create a formula
  formula <- as.formula(paste("diversity_shannon ~", var))
  
  # Perform Kruskal-Wallis test
  result <- baseline_alphadiversity %>% 
    kruskal_test(formula) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj")
  
  # Add the variable name to the result
  result <- result %>% mutate(variable = var)
  
  # Store the result in the list
  base_shannon_kruskal[[var]] <- result
}

# Combine all results into one data frame
kruskal_baseshannon_df <- bind_rows(base_shannon_kruskal)

# Print the combined results
View(kruskal_baseshannon_df)
##BMI_3_class
##BMI_SD
##Enterotype

##BMI_3_class
##pairwise comparison for BMI_3_class
bmi3shannon_pair <- baseline_alphadiversity %>% 
  dunn_test(diversity_shannon ~ BMI_3_class, p.adjust.method = "BH") 
View(bmi3shannon_pair)

##BMI_SD
BMISD_shannon_pair <- baseline_alphadiversity %>% 
  dunn_test(diversity_shannon ~ BMI_SD, p.adjust.method = "BH") 
BMISD_shannon_pair

##Enterotype
Enterotype_shannon_pair <- baseline_alphadiversity %>% 
  dunn_test(diversity_shannon ~ Enterotype, p.adjust.method = "BH") 
View(Enterotype_shannon_pair)


###################################################
# Initialize an empty list to store results
base_evenness_kruskal <- list()

# Loop through each variable and perform the Kruskal-Wallis test on shannon diversity
for (var in krus_variables) {
  # Create a formula
  formula <- as.formula(paste("evenness_pielou ~", var))
  
  # Perform Kruskal-Wallis test
  result <- baseline_alphadiversity %>% 
    kruskal_test(formula) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj")
  
  # Add the variable name to the result
  result <- result %>% mutate(variable = var)
  
  # Store the result in the list
  base_evenness_kruskal[[var]] <- result
}

# Combine all results into one data frame
kruskal_basevenness_df <- bind_rows(base_evenness_kruskal)

# Print the combined results
View(kruskal_basevenness_df)

##Enterotype
##Baseline_enterotype

##BMI_3_class
##pairwise comparison for BMI_3_class
bmi3evenness_pair <- baseline_alphadiversity %>% 
  dunn_test(evenness_pielou ~ BMI_3_class, p.adjust.method = "BH") 
View(bmi3evenness_pair)

##BMI_SD
BMISD_evenness_pair <- baseline_alphadiversity %>% 
  dunn_test(evenness_pielou ~ BMI_SD, p.adjust.method = "BH") 
View(BMISD_evenness_pair)

##Enterotype
Enterotype_evenness_pair <- baseline_alphadiversity %>% 
  dunn_test(evenness_pielou ~ Enterotype, p.adjust.method = "BH") 
Enterotype_evenness_pair

####
##kruskal wallis test
week12_alphadiversity <- alphadiversity_categoryk %>% filter(Time == "12-wks")

################################################################################
##Week 12 variables
kruswe12_variables <- c("Enterotype", "BMI_3_class")

# Initialize an empty list to store results
week12_observed_kruskal <- list()

# Loop through each variable and perform the Kruskal-Wallis test on shannon diversity
for (var in kruswe12_variables) {
  # Create a formula
  formula <- as.formula(paste("observed ~", var))
  
  # Perform Kruskal-Wallis test
  result <- week12_alphadiversity %>% 
    kruskal_test(formula) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj")
  
  # Add the variable name to the result
  result <- result %>% mutate(variable = var)
  
  # Store the result in the list
  week12_observed_kruskal[[var]] <- result
}

# Combine all results into one data frame
kruskal_week12evenness_df <- bind_rows(week12_observed_kruskal)

# Print the combined results
View(kruskal_week12evenness_df)
##Enterotype
##Baseline_enterotype

##Enterotype
Enterotype_evenness_pair <- week12_alphadiversity %>% 
  dunn_test(evenness_pielou ~ Enterotype, p.adjust.method = "BH") 
View(Enterotype_evenness_pair)

##Shannon diversity
# Initialize an empty list to store results
week12_shannon_kruskal <- list()

# Loop through each variable and perform the Kruskal-Wallis test on shannon diversity
for (var in kruswe12_variables) {
  # Create a formula
  formula <- as.formula(paste("diversity_shannon ~", var))
  
  # Perform Kruskal-Wallis test
  result <- week12_alphadiversity %>% 
    kruskal_test(formula) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj")
  
  # Add the variable name to the result
  result <- result %>% mutate(variable = var)
  
  # Store the result in the list
  week12_shannon_kruskal[[var]] <- result
}

# Combine all results into one data frame
kruskal_week12shannon_df <- bind_rows(week12_shannon_kruskal)

# Print the combined results
View(kruskal_week12shannon_df)

##Enterotype
##pairwise comparison for BMI_3_class
enterotypeshannon_pair <- week12_alphadiversity %>% 
  dunn_test(diversity_shannon ~ Enterotype, p.adjust.method = "BH") 
View(enterotypeshannon_pair)

##Faith pd
# Initialize an empty list to store results
week12_faithpd_kruskal <- list()

# Loop through each variable and perform the Kruskal-Wallis test on shannon diversity
for (var in kruswe12_variables) {
  # Create a formula
  formula <- as.formula(paste("pd ~", var))
  
  # Perform Kruskal-Wallis test
  result <- week12_alphadiversity %>% 
    kruskal_test(formula) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj")
  
  # Add the variable name to the result
  result <- result %>% mutate(variable = var)
  
  # Store the result in the list
  week12_faithpd_kruskal[[var]] <- result
}

# Combine all results into one data frame
kruskal_week12faithpd_df <- bind_rows(week12_faithpd_kruskal)

# Print the combined results
View(kruskal_week12faithpd_df)

##Enterotype
##pairwise comparison for BMI_3_class
enterotypefaithpd_pair <- week12_alphadiversity %>% 
  dunn_test(pd ~ Enterotype, p.adjust.method = "BH") 
View(enterotypefaithpd_pair)

##BMI_3_class
##pairwise comparison for BMI_3_class
bmi3faithpd_pair <- week12_alphadiversity %>% 
  dunn_test(pd ~ BMI_3_class, p.adjust.method = "BH") 
View(bmi3faithpd_pair)

##kruskal wallis test
week24_alphadiversity <- alphadiversity_categoryk %>% filter(Time == "24-wks")

##Faith pd
# Initialize an empty list to store results
week24_faithpd_kruskal <- list()

# Loop through each variable and perform the Kruskal-Wallis test on shannon diversity
for (var in kruswe12_variables) {
  # Create a formula
  formula <- as.formula(paste("pd ~", var))
  
  # Perform Kruskal-Wallis test
  result <- week24_alphadiversity %>% 
    kruskal_test(formula) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj")
  
  # Add the variable name to the result
  result <- result %>% mutate(variable = var)
  
  # Store the result in the list
  week24_faithpd_kruskal[[var]] <- result
}

# Combine all results into one data frame
kruskal_week24faithpd_df <- bind_rows(week24_faithpd_kruskal)

# Print the combined results
View(kruskal_week24faithpd_df)

##Enterotype
##pairwise comparison for BMI_3_class
wk24enterotypefaithpd_pair <- week24_alphadiversity %>% 
  dunn_test(pd ~ Enterotype, p.adjust.method = "BH") 
View(wk24enterotypefaithpd_pair)

#################################################
##kruskal wallis test
week24_alphadiversity <- alphadiversity_categoryk %>% filter(Time == "24-wks")

##Faith pd
# Initialize an empty list to store results
week24_shannon_kruskal <- list()

# Loop through each variable and perform the Kruskal-Wallis test on shannon diversity
for (var in kruswe12_variables) {
  # Create a formula
  formula <- as.formula(paste("diversity_shannon ~", var))
  
  # Perform Kruskal-Wallis test
  result <- week24_alphadiversity %>% 
    kruskal_test(formula) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj")
  
  # Add the variable name to the result
  result <- result %>% mutate(variable = var)
  
  # Store the result in the list
  week24_shannon_kruskal[[var]] <- result
}

# Combine all results into one data frame
kruskal_week24shannon_df <- bind_rows(week24_shannon_kruskal)

# Print the combined results
View(kruskal_week24shannon_df)

##Enterotype
##pairwise comparison for Enterotype
wk24enterotypeshannon_pair <- week24_alphadiversity %>% 
  dunn_test(diversity_shannon ~ Enterotype, p.adjust.method = "BH") 
View(wk24enterotypeshannon_pair)

##Observed
# Initialize an empty list to store results
week24_observed_kruskal <- list()

# Loop through each variable and perform the Kruskal-Wallis test on shannon diversity
for (var in kruswe12_variables) {
  # Create a formula
  formula <- as.formula(paste("observed ~", var))
  
  # Perform Kruskal-Wallis test
  result <- week24_alphadiversity %>% 
    kruskal_test(formula) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj")
  
  # Add the variable name to the result
  result <- result %>% mutate(variable = var)
  
  # Store the result in the list
  week24_observed_kruskal[[var]] <- result
}

# Combine all results into one data frame
kruskal_week24observed_df <- bind_rows(week24_observed_kruskal)

# Print the combined results
View(kruskal_week24observed_df)

##Enterotype
##pairwise comparison for Enterotype
wk24enterotypeobserved_pair <- week24_alphadiversity %>% 
  dunn_test(observed ~ Enterotype, p.adjust.method = "BH") 
View(wk24enterotypeobserved_pair)
