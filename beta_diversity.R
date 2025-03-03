set.seed(562)
library(phyloseq)
library(vegan)
library(data.table)
library(ggplot2)
library(readr)
library(pairwiseAdonis)

c("#1B9E77", "#D95F02", "#006BA3")

##merge at the genus level
ps_genus_rare<- tax_glom(ps_study810_rare, taxrank = "Genus")

df_ps_genus_rare <- data.frame(sample_data(ps_genus_rare))

df_ps_genus_rare %>% count(Family_dementia)

##Factorize selected variables
# Convert the specified variables to factors
df_ps_genus_rare <- df_ps_genus_rare %>%
  mutate(APOE_status = factor(APOE_status), BMI_3_class = factor(BMI_3_class),Screening_Diagnosis = factor(Screening_Diagnosis),
         BMI_SD = factor(BMI_SD),Gender = factor(Gender),Age_quartile = factor(Age_quartile),Family_dementia = factor(Family_dementia),Smoker = factor(Smoker),
         Coronary_heart_disease = factor(Coronary_heart_disease),Hypertension_treatment = factor(Hypertension_treatment),Diabetes_mellitus_all = factor(Diabetes_mellitus_all),
         Antithrombotic_agents = factor(Antithrombotic_agents),Calcium_channel_blockers = factor(Calcium_channel_blockers),Thyroid_therapy = factor(Thyroid_therapy),
         Drug_acid_related_disorders = factor(Drug_acid_related_disorders), Hypercholesterolemia_drug = factor(Hypercholesterolemia_drug),
         Residence = factor(Residence), BMI_SD = factor(BMI_SD))

##metadata from rarefied phyloseq object
ps_genus_rare

base_genus_rare <- subset_samples(ps_genus_rare, Time == "Baseline") ##295

base_genus_rare <- prune_taxa(taxa_sums(base_genus_rare) > 0, base_genus_rare)  ##283

prune_taxa(taxa_sums(base_genus_rare) > 1, base_genus_rare)  ##279 (number of singletons is 4)

prune_taxa(taxa_sums(base_genus_rare) > 2, base_genus_rare)  ##273 (number of doubletons is 10)

##Calculate relative abundance
basegenusrare_comp <- microbiome::transform(base_genus_rare, transform = "compositional")

##Baseline
base_ps_genus_rare <- df_ps_genus_rare %>% filter(Time == "Baseline")

#Bray-Curtis at genus level
base_genus_bray <- phyloseq::distance(basegenusrare_comp, method = "bray")

basegenus_ordi <- ordinate(basegenusrare_comp, distance = "bray")

braygenusfit <- envfit(basegenus_ordi, envfit_baselinemetadata, perm = 999, na.rm = TRUE)

set.seed(1988)
bray_genus_baseline_terms <- adonis2(formula = base_genus_bray ~ Screening_Diagnosis + Residence + 
                                       Age_quartile + Gender + BMI_3_class + Coronary_heart_disease + 
                                       Hypercholesterolemia_drug +Thyroid_therapy + Family_dementia + Hypertension_treatment + Metabolic_disease + Diabetes_mellitus_all + Antithrombotic_agents + Drug_acid_related_disorders, 
                                     data = base_ps_genus_rare, permutations = 999, by = "terms", na.action = na.exclude)

##Bray curtis at the genus level relative abundance
set.seed(1988)
braygenus_screening <- adonis2(formula = base_genus_bray ~ Screening_Diagnosis, data = base_ps_genus_rare, permutations = 999,
                               na.action = na.exclude)
braygenus_residence <- adonis2(formula = base_genus_bray ~ Residence, data = base_ps_genus_rare, permutations = 999,
                               na.action = na.exclude)
braygenus_agequartile <- adonis2(formula = base_genus_bray ~ Age_quartile, data = base_ps_genus_rare, permutations = 999,
                                 na.action = na.exclude)
braygenus_gender <- adonis2(formula = base_genus_bray ~ Gender, data = base_ps_genus_rare, permutations = 999,
                            na.action = na.exclude) ##Only age partitioned by quartile is significantly different interms of microbiome structure
braygenus_bmi3 <- adonis2(formula = base_genus_bray ~ BMI_3_class, data = base_ps_genus_rare, permutations = 999,
                          na.action = na.exclude) 
braygenus_coronary <- adonis2(formula = base_genus_bray ~ Coronary_heart_disease, data = base_ps_genus_rare, permutations = 999,
                              na.action = na.exclude) 
braygenus_hypercholest <- adonis2(formula = base_genus_bray ~ Hypercholesterolemia_drug, data = base_ps_genus_rare, permutations = 999,
                                  na.action = na.exclude) 
braygenus_thyroid <- adonis2(formula = base_genus_bray ~ Thyroid_therapy, data = base_ps_genus_rare, permutations = 999,
                             na.action = na.exclude) 
braygenus_dementia <- adonis2(formula = base_genus_bray ~ Family_dementia, data = base_ps_genus_rare, permutations = 999,
                              na.action = na.exclude) 
braygenus_hypertension <- adonis2(formula = base_genus_bray ~ Hypertension_treatment, data = base_ps_genus_rare, permutations = 999,
                                  na.action = na.exclude) 
braygenus_diabetes <- adonis2(formula = base_genus_bray ~ Diabetes_mellitus_all, data = base_ps_genus_rare, permutations = 999,na.action = na.exclude) 

braygenus_antithromboitic <- adonis2(formula = base_genus_bray ~ Antithrombotic_agents, data = base_ps_genus_rare, permutations = 999,
                                     na.action = na.exclude) 
braygenus_drugacid <- adonis2(formula = base_genus_bray ~ Drug_acid_related_disorders, data = base_ps_genus_rare, permutations = 999,
                              na.action = na.exclude) 

braygenus_bmisd <- adonis2(formula = base_genus_bray ~ BMI_SD, data = base_ps_genus_rare, permutations = 999,
                              na.action = na.exclude) 

braygenus_calcium <- adonis2(formula = base_genus_bray ~ Calcium_channel_blockers, data = base_ps_genus_rare, permutations = 999,
                           na.action = na.exclude)

base_genus_bray
#Bray-Curtis at genus level
base_genus_bray <- phyloseq::distance(basegenusrare_comp, method = "bray")

##ordination
basegenus_ordi <- ordinate(basegenusrare_comp, distance = "bray")

##plot
# Calculate centroids by Gender
plot_ordination(basegenusrare_comp,basegenus_ordi, color = 'Age_quartile') + 
  geom_point(size = 2) +  # Change the size of points to 3 (you can adjust the value as needed)
  stat_ellipse(linewidth = 1, linetype = 'solid') + 
  theme_classic() +
  labs(x = "", y = "") +
  theme(legend.position = 'right', axis.text = element_text(size = 10), axis.title = element_text(face = "bold", size = 10),
        legend.title = element_text(hjust = 0.5, face = "bold"),
        legend.margin = margin(0,0,0,-25)) +
  labs(x = "PCoA 1 (21.1%)", y = "PCoA 2 (15.3%)", color = "Biopsy approach")


##Pairwise comparison
brayagequartile_pair <- pairwise.adonis2(base_genus_bray ~ Age_quartile, data=base_ps_genus_rare)

###Bray-Curtis at the ASV level
ps_study8_rare

##subset Baseline level
base_study8_rare <- subset_samples(ps_study8_rare, Time == "Baseline") ##321

base_study8_rare <- prune_taxa(taxa_sums(base_study8_rare) > 0, base_study8_rare)  ##296

#Bray-Curtis at the ASV level
base_rare_bray <- phyloseq::distance(base_study8_rare, method = "bray")
baseasv_ordi <- ordinate(base_study8_rare, distance = "bray", "PCoA")


envfit_baselinemetadata <- base_ps_genus_rare %>% select(Treatment_Type, APOE_status, BMI_3_class, Screening_Diagnosis,
                              BMI_SD,Residence, Gender, Age_quartile, Age_median, Family_dementia,
                              Coronary_heart_disease, Hypertension_treatment, Hypercholesterolemia_drug,
                              Diabetes_mellitus_all,Antithrombotic_agents, Calcium_channel_blockers,
                              Antiinflammatory_antirheumatic_drugs, Thyroid_therapy, Drug_acid_related_disorders)

brayasvfit <- envfit(baseasv_ordi, envfit_baselinemetadata, perm = 999, na.rm = TRUE)

##Bray curtis at the genus level relative abundance
set.seed(1989)
brayasv_screening <- adonis2(formula = base_rare_bray ~ Screening_Diagnosis, data = base_ps_genus_rare, permutations = 999,
                               na.action = na.exclude)
brayasv_residence <- adonis2(formula = base_rare_bray ~ Residence, data = base_ps_genus_rare, permutations = 999,
                               na.action = na.exclude)
brayasv_agequartile <- adonis2(formula = base_rare_bray ~ Age_quartile, data = base_ps_genus_rare, permutations = 999,
                                 na.action = na.exclude)
brayasv_gender <- adonis2(formula = base_rare_bray ~ Gender, data = base_ps_genus_rare, permutations = 999,
                            na.action = na.exclude) ##Only age partitioned by quartile is significantly different interms of microbiome structure
brayasv_bmi3 <- adonis2(formula = base_rare_bray ~ BMI_3_class, data = base_ps_genus_rare, permutations = 999,
                          na.action = na.exclude) 
brayasv_coronary <- adonis2(formula = base_rare_bray ~ Coronary_heart_disease, data = base_ps_genus_rare, permutations = 999,
                              na.action = na.exclude) 
brayasv_hypercholest <- adonis2(formula = base_rare_bray ~ Hypercholesterolemia_drug, data = base_ps_genus_rare, permutations = 999,
                                  na.action = na.exclude) 
brayasv_thyroid <- adonis2(formula = base_rare_bray ~ Thyroid_therapy, data = base_ps_genus_rare, permutations = 999,
                             na.action = na.exclude) 
brayasv_dementia <- adonis2(formula = base_rare_bray ~ Family_dementia, data = base_ps_genus_rare, permutations = 999,
                              na.action = na.exclude) 
brayasv_hypertension <- adonis2(formula = base_rare_bray ~ Hypertension_treatment, data = base_ps_genus_rare, permutations = 999,
                                  na.action = na.exclude) 
brayasv_diabetes <- adonis2(formula = base_rare_bray ~ Diabetes_mellitus_all, data = base_ps_genus_rare, permutations = 999,na.action = na.exclude) 

brayasv_antithromboitic <- adonis2(formula = base_rare_bray ~ Antithrombotic_agents, data = base_ps_genus_rare, permutations = 999,
                                     na.action = na.exclude) 
brayasv_drugacid <- adonis2(formula = base_rare_bray ~ Drug_acid_related_disorders, data = base_ps_genus_rare, permutations = 999,
                              na.action = na.exclude) 
brayasv_bmisd <- adonis2(formula = base_rare_bray ~ BMI_SD, data = base_ps_genus_rare, permutations = 999,
                           na.action = na.exclude) 

##Jaccard at the ASV level
#Bray-Curtis at the ASV level
base_rare_jaccard <- phyloseq::distance(base_study8_rare, method = "jaccard")

##Bray curtis at the genus level relative abundance
set.seed(1990)
jaccard_screening <- adonis2(formula = base_rare_jaccard ~ Screening_Diagnosis, data = base_ps_genus_rare, permutations = 999,
                             na.action = na.exclude)
jaccard_residence <- adonis2(formula = base_rare_jaccard ~ Residence, data = base_ps_genus_rare, permutations = 999,
                             na.action = na.exclude)
jaccard_agequartile <- adonis2(formula = base_rare_jaccard ~ Age_quartile, data = base_ps_genus_rare, permutations = 999,
                               na.action = na.exclude)
jaccard_gender <- adonis2(formula = base_rare_jaccard ~ Gender, data = base_ps_genus_rare, permutations = 999,
                          na.action = na.exclude) ##Only age partitioned by quartile is significantly different interms of microbiome structure
jaccard_bmi3 <- adonis2(formula = base_rare_jaccard ~ BMI_3_class, data = base_ps_genus_rare, permutations = 999,
                        na.action = na.exclude) 
jaccard_coronary <- adonis2(formula = base_rare_jaccard ~ Coronary_heart_disease, data = base_ps_genus_rare, permutations = 999,
                            na.action = na.exclude) 
jaccard_hypercholest <- adonis2(formula = base_rare_jaccard ~ Hypercholesterolemia_drug, data = base_ps_genus_rare, permutations = 999,
                                na.action = na.exclude) 
jaccard_thyroid <- adonis2(formula = base_rare_jaccard ~ Thyroid_therapy, data = base_ps_genus_rare, permutations = 999,
                           na.action = na.exclude) 
jaccard_dementia <- adonis2(formula = base_rare_jaccard ~ Family_dementia, data = base_ps_genus_rare, permutations = 999,
                            na.action = na.exclude) 
jaccard_hypertension <- adonis2(formula = base_rare_jaccard ~ Hypertension_treatment, data = base_ps_genus_rare, permutations = 999,
                                na.action = na.exclude) 
jaccard_diabetes <- adonis2(formula = base_rare_jaccard ~ Diabetes_mellitus_all, data = base_ps_genus_rare, permutations = 999,na.action = na.exclude) 
jaccard_antithromboitic <- adonis2(formula = base_rare_jaccard ~ Antithrombotic_agents, data = base_ps_genus_rare, permutations = 999,
                                   na.action = na.exclude) 
jaccard_drugacid <- adonis2(formula = base_rare_jaccard ~ Drug_acid_related_disorders, data = base_ps_genus_rare, permutations = 999,
                            na.action = na.exclude) 
jaccard_bmisd <- adonis2(formula = base_rare_jaccard ~ BMI_SD, data = base_ps_genus_rare, permutations = 999,
                         na.action = na.exclude) 

##Jaccard at the genus level 
##Calculate relative abundance
basegenusrare_comp <- microbiome::transform(base_genus_rare, transform = "compositional")

##Baseline
base_ps_genus_rare <- df_ps_genus_rare %>% filter(Time == "Baseline")
trial <- envfit_baselinemetadata %>% select(-Age_median)

#Bray-Curtis at genus level
base_genus_jaccard <- phyloseq::distance(basegenusrare_comp, method = "jaccard")
basegenus_jacc <- ordinate(base_study8_rare, distance = "jaccard")

braygenujaccfit <- envfit(basegenus_jacc, trial, perm = 999, na.rm = TRUE)

##Bray curtis at the genus level relative abundance
set.seed(1998)
jaccgenus_screening <- adonis2(formula = base_genus_jaccard ~ Screening_Diagnosis, data = base_ps_genus_rare, permutations = 999,
                             na.action = na.exclude)
jaccgenus_residence <- adonis2(formula = base_genus_jaccard ~ Residence, data = base_ps_genus_rare, permutations = 999,
                             na.action = na.exclude)
jaccgenus_agequartile <- adonis2(formula = base_genus_jaccard ~ Age_quartile, data = base_ps_genus_rare, permutations = 999,
                               na.action = na.exclude)
jaccgenus_gender <- adonis2(formula = base_genus_jaccard ~ Gender, data = base_ps_genus_rare, permutations = 999,
                          na.action = na.exclude) ##Only age partitioned by quartile is significantly different interms of microbiome structure
jaccgenus_bmi3 <- adonis2(formula = base_genus_jaccard ~ BMI_3_class, data = base_ps_genus_rare, permutations = 999,
                        na.action = na.exclude) 
jaccgenus_coronary <- adonis2(formula = base_genus_jaccard ~ Coronary_heart_disease, data = base_ps_genus_rare, permutations = 999,
                            na.action = na.exclude) 
jaccgenus_hypercholest <- adonis2(formula = base_genus_jaccard ~ Hypercholesterolemia_drug, data = base_ps_genus_rare, permutations = 999,
                                na.action = na.exclude) 
jaccgenus_thyroid <- adonis2(formula = base_genus_jaccard ~ Thyroid_therapy, data = base_ps_genus_rare, permutations = 999,
                           na.action = na.exclude) 
jaccgenus_dementia <- adonis2(formula = base_genus_jaccard ~ Family_dementia, data = base_ps_genus_rare, permutations = 999,
                            na.action = na.exclude) 
jaccgenus_hypertension <- adonis2(formula = base_genus_jaccard ~ Hypertension_treatment, data = base_ps_genus_rare, permutations = 999,
                                na.action = na.exclude) 
jaccgenus_diabetes <- adonis2(formula = base_genus_jaccard ~ Diabetes_mellitus_all, data = base_ps_genus_rare, permutations = 999,na.action = na.exclude) 
jaccgenus_antithromboitic <- adonis2(formula = base_genus_jaccard ~ Antithrombotic_agents, data = base_ps_genus_rare, permutations = 999,
                                   na.action = na.exclude) 
jaccgenus_drugacid <- adonis2(formula = base_genus_jaccard ~ Drug_acid_related_disorders, data = base_ps_genus_rare, permutations = 999,
                            na.action = na.exclude) 
jaccgenus_bmisd <- adonis2(formula = base_genus_jaccard ~ BMI_SD, data = base_ps_genus_rare, permutations = 999,
                         na.action = na.exclude) 

##export the refseq to qiime2 and make tree

new_tre <- ape::multi2di(phy_tree(base_study8_rare))
phy_tree(base_study8_rare) <- new_tre

##Weighted UniFrac
#Bray-Curtis at the ASV level
base_rare_wunifrac <- phyloseq::distance(base_study8_rare, method = "wunifrac")


##Weighted UniFrac at the ASV level
set.seed(1990)
wunifrac_screening <- adonis2(formula = base_rare_wunifrac ~ Screening_Diagnosis, data = base_ps_genus_rare, permutations = 999,
                             na.action = na.exclude)
wunifrac_residence <- adonis2(formula = base_rare_wunifrac ~ Residence, data = base_ps_genus_rare, permutations = 999,
                             na.action = na.exclude)
wunifrac_agequartile <- adonis2(formula = base_rare_wunifrac ~ Age_quartile, data = base_ps_genus_rare, permutations = 999,
                               na.action = na.exclude)
wunifrac_gender <- adonis2(formula = base_rare_wunifrac ~ Gender, data = base_ps_genus_rare, permutations = 999,
                          na.action = na.exclude) ##Only age partitioned by quartile is significantly different interms of microbiome structure
wunifrac_bmi3 <- adonis2(formula = base_rare_wunifrac ~ BMI_3_class, data = base_ps_genus_rare, permutations = 999,
                        na.action = na.exclude) 
wunifrac_coronary <- adonis2(formula = base_rare_wunifrac ~ Coronary_heart_disease, data = base_ps_genus_rare, permutations = 999,
                            na.action = na.exclude) 
wunifrac_hypercholest <- adonis2(formula = base_rare_wunifrac ~ Hypercholesterolemia_drug, data = base_ps_genus_rare, permutations = 999,
                                na.action = na.exclude) 
wunifrac_thyroid <- adonis2(formula = base_rare_wunifrac ~ Thyroid_therapy, data = base_ps_genus_rare, permutations = 999,
                           na.action = na.exclude) 
wunifrac_dementia <- adonis2(formula = base_rare_wunifrac ~ Family_dementia, data = base_ps_genus_rare, permutations = 999,
                            na.action = na.exclude) 
wunifrac_hypertension <- adonis2(formula = base_rare_wunifrac ~ Hypertension_treatment, data = base_ps_genus_rare, permutations = 999,
                                na.action = na.exclude) 
wunifrac_diabetes <- adonis2(formula = base_rare_wunifrac ~ Diabetes_mellitus_all, data = base_ps_genus_rare, permutations = 999,na.action = na.exclude) 
wunifrac_antithromboitic <- adonis2(formula = base_rare_wunifrac ~ Antithrombotic_agents, data = base_ps_genus_rare, permutations = 999,
                                   na.action = na.exclude) 
wunifrac_drugacid <- adonis2(formula = base_rare_wunifrac ~ Drug_acid_related_disorders, data = base_ps_genus_rare, permutations = 999,
                            na.action = na.exclude) 
wunifrac_bmisd <- adonis2(formula = base_rare_wunifrac ~ BMI_SD, data = base_ps_genus_rare, permutations = 999,
                         na.action = na.exclude) 

wunifasvage_quartile_pair <- pairwise.adonis2(base_rare_wunifrac ~ Age_quartile, data=base_ps_genus_rare)

##Jaccard at the genus level 
##Calculate relative abundance
basegenusrare_comp <- microbiome::transform(base_genus_rare, transform = "compositional")

##Baseline
base_ps_genus_rare <- df_ps_genus_rare %>% filter(Time == "Baseline")

#Bray-Curtis at genus level
base_genus_wunifrac <- phyloseq::distance(basegenusrare_comp, method = "wunifrac")

##Weighted UniFrac at the ASV level
set.seed(1999)
wunifracgenus_screening <- adonis2(formula = base_genus_wunifrac ~ Screening_Diagnosis, data = base_ps_genus_rare, permutations = 999,
                              na.action = na.exclude)
wunifracgenus_residence <- adonis2(formula = base_genus_wunifrac ~ Residence, data = base_ps_genus_rare, permutations = 999,
                              na.action = na.exclude)
wunifracgenus_agequartile <- adonis2(formula = base_genus_wunifrac ~ Age_quartile, data = base_ps_genus_rare, permutations = 999,
                                na.action = na.exclude)
wunifracgenus_gender <- adonis2(formula = base_genus_wunifrac ~ Gender, data = base_ps_genus_rare, permutations = 999,
                           na.action = na.exclude) ##Only age partitioned by quartile is significantly different interms of microbiome structure
wunifracgenus_bmi3 <- adonis2(formula = base_genus_wunifrac ~ BMI_3_class, data = base_ps_genus_rare, permutations = 999,
                         na.action = na.exclude) 
wunifracgenus_coronary <- adonis2(formula = base_genus_wunifrac ~ Coronary_heart_disease, data = base_ps_genus_rare, permutations = 999,
                             na.action = na.exclude) 
wunifracgenus_hypercholest <- adonis2(formula = base_genus_wunifrac ~ Hypercholesterolemia_drug, data = base_ps_genus_rare, permutations = 999,
                                 na.action = na.exclude) 
wunifracgenus_thyroid <- adonis2(formula = base_genus_wunifrac ~ Thyroid_therapy, data = base_ps_genus_rare, permutations = 999,
                            na.action = na.exclude) 
wunifracgenus_dementia <- adonis2(formula = base_genus_wunifrac ~ Family_dementia, data = base_ps_genus_rare, permutations = 999,
                             na.action = na.exclude) 
wunifracgenus_hypertension <- adonis2(formula = base_genus_wunifrac ~ Hypertension_treatment, data = base_ps_genus_rare, permutations = 999,
                                 na.action = na.exclude) 
wunifracgenus_diabetes <- adonis2(formula = base_genus_wunifrac ~ Diabetes_mellitus_all, data = base_ps_genus_rare, permutations = 999,na.action = na.exclude) 
wunifracgenus_antithromboitic <- adonis2(formula = base_genus_wunifrac ~ Antithrombotic_agents, data = base_ps_genus_rare, permutations = 999,
                                    na.action = na.exclude) 
wunifracgenus_drugacid <- adonis2(formula = base_genus_wunifrac ~ Drug_acid_related_disorders, data = base_ps_genus_rare, permutations = 999,
                             na.action = na.exclude) 
wunifracgenus_bmisd <- adonis2(formula = base_genus_wunifrac ~ BMI_SD, data = base_ps_genus_rare, permutations = 999,
                          na.action = na.exclude) 

wunifgenusage_quartile_pair <- pairwise.adonis2(base_genus_wunifrac ~ Age_quartile, data=base_ps_genus_rare)

##Unweighted UniFrac
uunifrac
base_rare_uunifrac <- phyloseq::distance(base_study8_rare, method = "uunifrac")

##Weighted UniFrac at the ASV level
set.seed(1991)
uunifrac_screening <- adonis2(formula = base_rare_uunifrac ~ Screening_Diagnosis, data = base_ps_genus_rare, permutations = 999,
                              na.action = na.exclude)
uunifrac_residence <- adonis2(formula = base_rare_uunifrac ~ Residence, data = base_ps_genus_rare, permutations = 999,
                              na.action = na.exclude)
uunifrac_agequartile <- adonis2(formula = base_rare_uunifrac ~ Age_quartile, data = base_ps_genus_rare, permutations = 999,
                                na.action = na.exclude)
uunifrac_gender <- adonis2(formula = base_rare_uunifrac ~ Gender, data = base_ps_genus_rare, permutations = 999,
                           na.action = na.exclude) ##Only age partitioned by quartile is significantly different interms of microbiome structure
uunifrac_bmi3 <- adonis2(formula = base_rare_uunifrac ~ BMI_3_class, data = base_ps_genus_rare, permutations = 999,
                         na.action = na.exclude) 
uunifrac_coronary <- adonis2(formula = base_rare_uunifrac ~ Coronary_heart_disease, data = base_ps_genus_rare, permutations = 999,
                             na.action = na.exclude) 
uunifrac_hypercholest <- adonis2(formula = base_rare_uunifrac ~ Hypercholesterolemia_drug, data = base_ps_genus_rare, permutations = 999,
                                 na.action = na.exclude) 
uunifrac_thyroid <- adonis2(formula = base_rare_uunifrac ~ Thyroid_therapy, data = base_ps_genus_rare, permutations = 999,
                            na.action = na.exclude) 
uunifrac_dementia <- adonis2(formula = base_rare_uunifrac ~ Family_dementia, data = base_ps_genus_rare, permutations = 999,
                             na.action = na.exclude) 
uunifrac_hypertension <- adonis2(formula = base_rare_uunifrac ~ Hypertension_treatment, data = base_ps_genus_rare, permutations = 999,
                                 na.action = na.exclude) 
uunifrac_diabetes <- adonis2(formula = base_rare_uunifrac ~ Diabetes_mellitus_all, data = base_ps_genus_rare, permutations = 999,na.action = na.exclude) 
uunifrac_antithromboitic <- adonis2(formula = base_rare_uunifrac ~ Antithrombotic_agents, data = base_ps_genus_rare, permutations = 999,
                                    na.action = na.exclude) 
uunifrac_drugacid <- adonis2(formula = base_rare_uunifrac ~ Drug_acid_related_disorders, data = base_ps_genus_rare, permutations = 999,
                             na.action = na.exclude) 
uunifrac_bmisd <- adonis2(formula = base_rare_uunifrac ~ BMI_SD, data = base_ps_genus_rare, permutations = 999,
                          na.action = na.exclude) 
##Calculate relative abundance
basegenusrare_comp <- microbiome::transform(base_genus_rare, transform = "compositional")

##Baseline
base_ps_genus_rare <- df_ps_genus_rare %>% filter(Time == "Baseline")

#Bray-Curtis at genus level
base_genus_uunifrac <- phyloseq::distance(basegenusrare_comp, method = "uunifrac")

##Weighted UniFrac at the ASV level
set.seed(2010)
uunifracgenus_screening <- adonis2(formula = base_genus_uunifrac ~ Screening_Diagnosis, data = base_ps_genus_rare, permutations = 999,
                                   na.action = na.exclude)
uunifracgenus_residence <- adonis2(formula = base_genus_uunifrac ~ Residence, data = base_ps_genus_rare, permutations = 999,
                                   na.action = na.exclude)
uunifracgenus_agequartile <- adonis2(formula = base_genus_uunifrac ~ Age_quartile, data = base_ps_genus_rare, permutations = 999,
                                     na.action = na.exclude)
uunifracgenus_gender <- adonis2(formula = base_genus_uunifrac ~ Gender, data = base_ps_genus_rare, permutations = 999,
                                na.action = na.exclude) ##Only age partitioned by quartile is significantly different interms of microbiome structure
uunifracgenus_bmi3 <- adonis2(formula = base_genus_uunifrac ~ BMI_3_class, data = base_ps_genus_rare, permutations = 999,
                              na.action = na.exclude) 
uunifracgenus_coronary <- adonis2(formula = base_genus_uunifrac ~ Coronary_heart_disease, data = base_ps_genus_rare, permutations = 999,
                                  na.action = na.exclude) 
uunifracgenus_hypercholest <- adonis2(formula = base_genus_uunifrac ~ Hypercholesterolemia_drug, data = base_ps_genus_rare, permutations = 999,
                                      na.action = na.exclude) 
uunifracgenus_thyroid <- adonis2(formula = base_genus_uunifrac ~ Thyroid_therapy, data = base_ps_genus_rare, permutations = 999,
                                 na.action = na.exclude) 
uunifracgenus_dementia <- adonis2(formula = base_genus_uunifrac ~ Family_dementia, data = base_ps_genus_rare, permutations = 999,
                                  na.action = na.exclude) 
uunifracgenus_hypertension <- adonis2(formula = base_genus_uunifrac ~ Hypertension_treatment, data = base_ps_genus_rare, permutations = 999,
                                      na.action = na.exclude) 
uunifracgenus_diabetes <- adonis2(formula = base_genus_uunifrac ~ Diabetes_mellitus_all, data = base_ps_genus_rare, permutations = 999,na.action = na.exclude) 
uunifracgenus_antithromboitic <- adonis2(formula = base_genus_uunifrac ~ Antithrombotic_agents, data = base_ps_genus_rare, permutations = 999,
                                         na.action = na.exclude) 
uunifracgenus_drugacid <- adonis2(formula = base_genus_uunifrac ~ Drug_acid_related_disorders, data = base_ps_genus_rare, permutations = 999,
                                  na.action = na.exclude) 
uunifracgenus_bmisd <- adonis2(formula = base_genus_uunifrac ~ BMI_SD, data = base_ps_genus_rare, permutations = 999,
                               na.action = na.exclude) 


wunifgenubmisd_quartile_pair <- pairwise.adonis2(base_genus_uunifrac ~ BMI_SD, data=base_ps_genus_rare, na.action = na.omit)



##Read the distance matrix calculate using qiime2 boots
##Bray curtis at the ASV level

braycurtis_matrix <- read_delim("q2boot_diversity/braycurtis-matrix.tsv", 
                                delim = "\t", escape_double = FALSE, 
                                trim_ws = TRUE)
braycurtis_matrix <- as.matrix(braycurtis_matrix)

rownames(braycurtis_matrix) <- braycurtis_matrix[, 1]  # Assign the first column as row names
braycurtis_matrix <- braycurtis_matrix[, -1]           # Remove the first column


##Baseline
base_ps_genus_rare <- df_ps_genus_rare %>% filter(Time == "Baseline")

#Filter the distance matrix
base_sample <- base_ps_genus_rare$...1 

# Filter the matrix rows and columns by these subject names
base_braycurtis_matrix <- braycurtis_matrix[base_sample, base_sample]

set.seed(1988)
bray_screening <- adonis2(formula = base_braycurtis_matrix ~ Screening_Diagnosis, data = base_ps_genus_rare, permutations = 999,
                          na.action = na.exclude)
bray_agequartile <- adonis2(formula = base_braycurtis_matrix ~ Age_quartile, data = base_ps_genus_rare, permutations = 999,
                            na.action = na.exclude)
bray_gender <- adonis2(formula = base_braycurtis_matrix ~ Gender, data = base_ps_genus_rare, permutations = 999,
                       na.action = na.exclude)
bray_bmi3 <- adonis2(formula = base_braycurtis_matrix ~ BMI_3_class, data = base_ps_genus_rare, permutations = 999,
                       na.action = na.exclude) 
bray_coronary <- adonis2(formula = base_braycurtis_matrix ~ Coronary_heart_disease, data = base_ps_genus_rare, permutations = 999,
                     na.action = na.exclude) 
bray_hypercholest <- adonis2(formula = base_braycurtis_matrix ~ Hypercholesterolemia_drug, data = base_ps_genus_rare, permutations = 999,
                         na.action = na.exclude) 
bray_thyroid <- adonis2(formula = base_braycurtis_matrix ~ Thyroid_therapy, data = base_ps_genus_rare, permutations = 999,
                             na.action = na.exclude) 
bray_dementia <- adonis2(formula = base_braycurtis_matrix ~ Family_dementia, data = base_ps_genus_rare, permutations = 999,
                        na.action = na.exclude) 
bray_hypertension <- adonis2(formula = base_braycurtis_matrix ~ Hypertension_treatment, data = base_ps_genus_rare, permutations = 999,
                         na.action = na.exclude) 
bray_diabetes <- adonis2(formula = base_braycurtis_matrix ~ Diabetes_mellitus_all, data = base_ps_genus_rare, permutations = 999,
                             na.action = na.exclude) 
bray_antithromboitic <- adonis2(formula = base_braycurtis_matrix ~ Antithrombotic_agents, data = base_ps_genus_rare, permutations = 999,
                         na.action = na.exclude) 
bray_drugacid <- adonis2(formula = base_braycurtis_matrix ~ Drug_acid_related_disorders, data = base_ps_genus_rare, permutations = 999,
                                na.action = na.exclude) 

###Bray curtis at the genus level
braygenus_matrix <- read_delim("q2boot_diversity/genus_diveristy/genus_braymatrix.tsv", 
                                delim = "\t", escape_double = FALSE, 
                                trim_ws = TRUE)
braygenus_matrix <- as.matrix(braygenus_matrix)

rownames(braygenus_matrix) <- braygenus_matrix[, 1]  # Assign the first column as row names
braygenus_matrix <- braygenus_matrix[, -1]           # Remove the first column

##Baseline
base_ps_genus_rare

#Filter the distance matrix
base_sample <- base_ps_genus_rare$...1 

# Filter the matrix rows and columns by these subject names
base_braygenus_matrix <- braygenus_matrix[base_sample, base_sample]

set.seed(1988)
braygenus_screening <- adonis2(formula = base_braygenus_matrix ~ Screening_Diagnosis, data = base_ps_genus_rare, permutations = 999,
                          na.action = na.exclude)
braygenus_agequartile <- adonis2(formula = base_braygenus_matrix ~ Age_quartile, data = base_ps_genus_rare, permutations = 999,
                            na.action = na.exclude)
braygenus_gender <- adonis2(formula = base_braygenus_matrix ~ Gender, data = base_ps_genus_rare, permutations = 999,
                       na.action = na.exclude)
braygenus_bmi3 <- adonis2(formula = base_braygenus_matrix ~ BMI_3_class, data = base_ps_genus_rare, permutations = 999,
                     na.action = na.exclude) 
braygenus_coronary <- adonis2(formula = base_braygenus_matrix ~ Coronary_heart_disease, data = base_ps_genus_rare, permutations = 999,
                         na.action = na.exclude) 
braygenus_hypercholest <- adonis2(formula = base_braygenus_matrix ~ Hypercholesterolemia_drug, data = base_ps_genus_rare, permutations = 999,
                             na.action = na.exclude) 
braygenus_thyroid <- adonis2(formula = base_braygenus_matrix ~ Thyroid_therapy, data = base_ps_genus_rare, permutations = 999,
                        na.action = na.exclude) 
braygenus_dementia <- adonis2(formula = base_braygenus_matrix ~ Family_dementia, data = base_ps_genus_rare, permutations = 999,
                         na.action = na.exclude) 
braygenus_hypertension <- adonis2(formula = base_braygenus_matrix ~ Hypertension_treatment, data = base_ps_genus_rare, permutations = 999,
                             na.action = na.exclude) 
braygenus_diabetes <- adonis2(formula = base_braygenus_matrix ~ Diabetes_mellitus_all, data = base_ps_genus_rare, permutations = 999,na.action = na.exclude) 

braygenus_antithromboitic <- adonis2(formula = base_braygenus_matrix ~ Antithrombotic_agents, data = base_ps_genus_rare, permutations = 999,
                                na.action = na.exclude) 
braygenus_drugacid <- adonis2(formula = base_braygenus_matrix ~ Drug_acid_related_disorders, data = base_ps_genus_rare, permutations = 999,
                         na.action = na.exclude) 