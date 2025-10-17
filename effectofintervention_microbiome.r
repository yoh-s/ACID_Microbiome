set.seed(562)

#Libraries 
library(phyloseq)
library(microbiome)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(vegan)
library(broom)
library(broom.mixed)
library(glmmTMB)
library(emmeans)
library(lme4)
library(lmerTest)
library(ggplot2)

#create output directories 
root_out <- "output_results"
dirs <- c("overall","BMI","enterotype","agequartile")
for (d in dirs) dir.create(file.path(root_out, d), recursive = TRUE, showWarnings = FALSE)

#safe adonis2 wrapper
# dist_mat: a dist object
# meta: dataframe with sample metadata, rownames = sample names
# formula_rhs: character, RHS of formula (e.g. "Treatment_Type * Time")
run_permanova <- function(dist_mat, meta, formula_rhs, permutations = 999, strata = NULL) {
  # prepare formula
  f <- as.formula(paste0("dist_mat ~ ", formula_rhs))
  # run adonis2: adonis2 accepts a distance matrix and a data frame
  ad <- tryCatch({
    adonis2(dist_mat ~ ., data = meta %>% dplyr::select(all_of(all.vars(as.formula(formula_rhs)))), permutations = permutations, by = "margin", strata = if (!is.null(strata)) meta[[strata]] else NULL)
  }, error = function(e) {
    warning("PERMANOVA failed: ", conditionMessage(e))
    return(NULL)
  })
  ad
}

#ANCOVA per feature (baseline-adjusted linear model)
# df_wide: columns for Baseline and Week24 named exactly, plus Treatment_Type and optionally stratifier
# returns tibble of results per ASV
ancova_per_feature <- function(df_wide, predictor = "Treatment_Type", baseline_col = "Baseline", outcome_col = "Week24") {
  # df_wide: must contain columns: ASV, Baseline, Week24, Treatment_Type (and possibly stratifier)
  res <- tryCatch({
    fit <- lm(as.formula(paste0(outcome_col, " ~ ", predictor, " + ", baseline_col)), data = df_wide)
    tidy_fit <- broom::tidy(fit)
    # return the tidy row for Treatment_Type if exists
    term_name <- paste0(predictor, "Intervention")
    row <- tidy_fit %>% filter(grepl(term_name, term))
    if (nrow(row) == 0) {
      # maybe predictor has different naming, return NA row
      return(tibble(term = term_name, estimate = NA_real_, std.error = NA_real_, statistic = NA_real_, p.value = NA_real_))
    }
    row
  }, error = function(e) {
    tibble(term = paste0(predictor,"_error"), estimate = NA_real_, std.error = NA_real_, statistic = NA_real_, p.value = NA_real_)
  })
  res
}
# NOTE: adapt these variables to your loaded objects; the script expects:
#- pspairgenus (phyloseq genus-level phyloseq used for paired analyses)
#- pspair (phyloseq at ASV/OTU level)
#- pspairgenus_clr (CLR-transformed genus-level phyloseq) OR script will create it
#- pspairgenus and pspair should have sample_data with columns:
#  sampleid (or rownames), subject, Time (Baseline/24-wks), Treatment_Type

# Ensure sample_data rownames are sample names
ensure_sample_names <- function(ps) {
  sdf <- data.frame(sample_data(ps))
  rownames(sdf) <- sample_names(ps)
  sample_data(ps) <- sdf
  ps
}

#ensure objects exist (if not, the script will stop and tell you)
if (!exists("pspair") | !exists("pspairgenus")) stop("Please load 'pspair' and 'pspairgenus' phyloseq objects into the environment before running this script.")

pspair <- ensure_sample_names(pspair)
pspairgenus <- ensure_sample_names(pspairgenus)

#Create CLR-transformed genus phyloseq if not present
if (!exists("pspairgenus_clr")) {
  pspairgenus_clr <- microbiome::transform(pspairgenus, "clr")
  # strip Species column if present in tax_table
  if (ncol(tax_table(pspairgenus_clr)) >= 7) tax_table(pspairgenus_clr) <- tax_table(pspairgenus_clr)[,1:6]
}

#Common metadata frame
meta_all <- data.frame(sample_data(pspair)) %>% rownames_to_column(var = "id")

#Overall analysis
#- baseline-adjusted ANCOVA per genus (wide -> Baseline & 24-wks)
#- mixed models across time (Treatment_Type * Time with random subject) optional
#- PERMANOVA for several distances (Bray, Jaccard, Weighted & Unweighted UniFrac)
{
  outdir <- file.path(root_out, "overall")

  #ANCOVA (genus-level): Baseline -> Week24
  #Prepare wide data: take CLR genus table for modeling (pspairgenus_clr)
  genus_counts <- as.data.frame(t(otu_table(pspairgenus_clr)))
  genus_counts <- genus_counts %>% rownames_to_column(var = "id")
  meta <- data.frame(sample_data(pspairgenus_clr)) %>% rownames_to_column(var = "id")
  fulldata <- left_join(meta, genus_counts, by = "id")

  # Keep samples that have both Baseline and 24-wks per subject per ASV: pivot-wide by ASV is performed below
  # Pivot long for per-ASV ancova building
  fulldata_long <- fulldata %>%
    pivot_longer(cols = starts_with("ASV"), names_to = "ASV", values_to = "value") %>%
    mutate(Time = factor(Time, levels = c("Baseline","24-wks")),
           Treatment_Type = factor(Treatment_Type, levels = c("Placebo", "Intervention")),
           subject = factor(subject))

  genus_wide <- fulldata_long %>%
    filter(Time %in% c("Baseline","24-wks")) %>%
    dplyr::select(subject, Treatment_Type, Time, ASV, value) %>%
    pivot_wider(names_from = Time, values_from = value) %>%
    dplyr::rename(Baseline = `Baseline`, Week24 = `24-wks`) %>%
    drop_na(Baseline, Week24) %>%
    mutate(Treatment_Type = factor(Treatment_Type, levels = c("Placebo","Intervention")))

  ancova_genus_res <- genus_wide %>%
    group_by(ASV) %>%
    nest() %>%
    mutate(ancova = map(data, ~ ancova_per_feature(.x, predictor = "Treatment_Type", baseline_col = "Baseline", outcome_col = "Week24"))) %>%
    unnest(ancova)

  # Attach taxonomy and BH adjust across ASVs
  tax_g <- data.frame(tax_table(pspairgenus_clr)) %>% rownames_to_column(var = "ASV")
  ancova_genus_res2 <- ancova_genus_res %>%
    left_join(tax_g, by = "ASV") %>%
    mutate(p.adj = p.adjust(p.value, method = "BH"))

  write_csv(ancova_genus_res2, file.path(outdir, "ancova_genus_overall.csv"))

  #PERMANOVA: define distance matrices and run adonis2 (with subject strata for longitudinal design)
  meta_for_dist <- data.frame(sample_data(pspair)) %>% rownames_to_column(var = "id") %>%
    mutate(Treatment_Type = factor(Treatment_Type, levels = c("Placebo","Intervention")),
           Time = factor(Time, levels = c("Baseline","24-wks")),
           subject = factor(subject))

  # Distances to compute
  # Bray (ASV)
  dist_bray_asv <- phyloseq::distance(pspair, method = "bray")
  ad_bray_asv <- adonis2(dist_bray_asv ~ Treatment_Type * Time, data = meta_for_dist, permutations = 999, strata = meta_for_dist$subject)
  write_csv(as.data.frame(ad_bray_asv), file.path(outdir, "permanova_bray_asv_overall.csv"))

  # Bray (genus)
  dist_bray_genus <- phyloseq::distance(pspairgenus, method = "bray")
  ad_bray_genus <- adonis2(dist_bray_genus ~ Treatment_Type * Time, data = meta_for_dist, permutations = 999, strata = meta_for_dist$subject)
  write_csv(as.data.frame(ad_bray_genus), file.path(outdir, "permanova_bray_genus_overall.csv"))

  # Jaccard (ASV) - presence/absence
  dist_jaccard <- phyloseq::distance(pspair, method = "jaccard", binary = TRUE)
  ad_jaccard <- adonis2(dist_jaccard ~ Treatment_Type * Time, data = meta_for_dist, permutations = 999, strata = meta_for_dist$subject)
  write_csv(as.data.frame(ad_jaccard), file.path(outdir, "permanova_jaccard_asv_overall.csv"))

  # Weighted UniFrac
  dist_wunifrac <- phyloseq::distance(pspair, method = "wunifrac")
  ad_wunifrac <- adonis2(dist_wunifrac ~ Treatment_Type * Time, data = meta_for_dist, permutations = 999, strata = meta_for_dist$subject)
  write_csv(as.data.frame(ad_wunifrac), file.path(outdir, "permanova_wunifrac_overall.csv"))

  # Unweighted UniFrac
  dist_unifrac <- phyloseq::distance(pspair, method = "unifrac")
  ad_unifrac <- adonis2(dist_unifrac ~ Treatment_Type * Time, data = meta_for_dist, permutations = 999, strata = meta_for_dist$subject)
  write_csv(as.data.frame(ad_unifrac), file.path(outdir, "permanova_unifrac_overall.csv"))

  #Alpha diversity between-group (mixed model)
  #compute alpha metrics into meta_for_dist if not present (observed, shannon, pd, pielou)
  alpha_metrics <- c("observed","diversity_shannon","pd","evenness_pielou")
  alpha_df <- meta_for_dist
  # if metrics missing, compute from pspair (ASV) or pspairgenus
  if (!all(alpha_metrics %in% colnames(alpha_df))) {
    # compute using estimate_richness for observed, shannon; pd from picante/phyloseq? here assume available
    ar <- estimate_richness(pspair, measures = c("Observed","Shannon"))
    ar <- ar %>% dplyr::rename(observed = Observed, diversity_shannon = Shannon)
    alpha_df <- alpha_df %>% left_join(ar %>% rownames_to_column("id"), by = "id")
    # add placeholders for pd and pielou if missing
    if (!"pd" %in% colnames(alpha_df)) alpha_df$pd <- NA_real_
    if (!"evenness_pielou" %in% colnames(alpha_df)) alpha_df$evenness_pielou <- NA_real_
  }
  # mixed model per metric
  alpha_results <- map_df(alpha_metrics, function(metric) {
    formula <- as.formula(paste0(metric, " ~ Treatment_Type * Time + (1|subject)"))
    # keep only columns needed
    df <- alpha_df %>% dplyr::select(all_of(c("id","subject","Time","Treatment_Type",metric)))
    colnames(df)[colnames(df) == metric] <- "metric_value"
    mod <- tryCatch(lmer(metric_value ~ Treatment_Type * Time + (1|subject), data = df), error = function(e) NULL)
    if (is.null(mod)) return(tibble(metric = metric, term = NA_character_, estimate = NA_real_, p.value = NA_real_, p.adj = NA_real_))
    # emmeans pairwise contrasts
    emm <- emmeans(mod, ~ Treatment_Type | Time)
    contr <- contrast(emm, "pairwise") %>% summary()
    contr$metric <- metric
    contr$p.adj <- p.adjust(contr$p.value, method = "BH")
    contr
  })

  write_csv(alpha_results, file.path(outdir, "alpha_mixedmodels_overall.csv"))

} #end overall

#By Baseline BMI
{
  outdir <- file.path(root_out, "BMI")

  # Ensure Baseline_BMI3 exists in metadata
  meta_df <- data.frame(sample_data(pspair)) %>% rownames_to_column("id")
  if (!"Baseline_BMI3" %in% colnames(meta_df)) stop("Baseline_BMI3 not found in sample data.")

  meta_df <- meta_df %>% mutate(Baseline_BMI3 = factor(Baseline_BMI3, levels = c("Normal","Overweight","Obesity")))

  # PERMANOVA stratified by BMI groups (Bray genus & ASV, Jaccard, W/Unifrac)
  bmi_groups <- levels(meta_df$Baseline_BMI3)
  perma_list <- list()

  for (g in bmi_groups) {
    meta_g <- meta_df %>% filter(Baseline_BMI3 == g)
    ids <- meta_g$id
    if (length(ids) < 3) {
      perma_list[[g]] <- NA
      next
    }
    # prune phyloseq to these samples
    ps_g_asv <- prune_samples(ids, pspair)
    ps_g_genus <- prune_samples(ids, pspairgenus)
    # distances
    dist_bray_a <- phyloseq::distance(ps_g_asv, method = "bray")
    dist_bray_g <- phyloseq::distance(ps_g_genus, method = "bray")
    dist_jacc <- phyloseq::distance(ps_g_asv, method = "jaccard", binary = TRUE)
    dist_wuni <- phyloseq::distance(ps_g_asv, method = "wunifrac")
    dist_uuni <- phyloseq::distance(ps_g_asv, method = "unifrac")

    meta_local <- data.frame(sample_data(ps_g_asv)) %>% rownames_to_column("id")
    # run adonis
    perma_list[[g]] <- list(
      bray_asv = adonis2(dist_bray_a ~ Treatment_Type * Time, data = meta_local, permutations = 999, strata = meta_local$subject),
      bray_genus = adonis2(dist_bray_g ~ Treatment_Type * Time, data = data.frame(sample_data(ps_g_genus)) %>% rownames_to_column("id"), permutations = 999, strata = meta_local$subject),
      jaccard = adonis2(dist_jacc ~ Treatment_Type * Time, data = meta_local, permutations = 999, strata = meta_local$subject),
      wunifrac = adonis2(dist_wuni ~ Treatment_Type * Time, data = meta_local, permutations = 999, strata = meta_local$subject),
      unifrac = adonis2(dist_uuni ~ Treatment_Type * Time, data = meta_local, permutations = 999, strata = meta_local$subject)
    )
  }

  # Save results to csv files per group
  for (g in names(perma_list)) {
    if (is.list(perma_list[[g]])) {
      for (nm in names(perma_list[[g]])) {
        write_csv(as.data.frame(perma_list[[g]][[nm]]), file.path(outdir, paste0("permanova_", g, "_", nm, ".csv")))
      }
    } else {
      # write placeholder
      write_lines("insufficient_samples", file.path(outdir, paste0("permanova_", g, "_NA.txt")))
    }
  }

  # ANCOVA per genus stratified by BMI is heavy — run ANCOVA_genus_bmi analogous to overall but filtered
  # Reuse fulldata_long from overall (it has Time, subject, Treatment_Type, ASV, value)
  if (exists("fulldata_long")) {
    genus_wide_bmi <- fulldata_long %>%
      filter(Time %in% c("Baseline","24-wks")) %>%
      dplyr::select(subject, Treatment_Type, Time, ASV, value, Baseline_BMI3) %>%
      pivot_wider(names_from = Time, values_from = value) %>%
      dplyr::rename(Baseline = `Baseline`, Week24 = `24-wks`) %>%
      drop_na(Baseline, Week24) %>%
      mutate(Treatment_Type = factor(Treatment_Type, levels = c("Placebo","Intervention")),
             Baseline_BMI3 = factor(Baseline_BMI3, levels = c("Normal","Overweight","Obesity")))

    ancova_genus_bmi <- genus_wide_bmi %>%
      group_by(ASV) %>%
      nest() %>%
      mutate(model = map(data, ~ tryCatch(lm(Week24 ~ Treatment_Type * Baseline_BMI3 + Baseline, data = .x), error = function(e) NULL))) %>%
      mutate(
        # Extract contrasts per BMI level
        res_normal = map(model, ~ if (!is.null(.x)) { emmeans(.x, ~ Treatment_Type | Baseline_BMI3, at = list(Baseline_BMI3 = "Normal")) %>% contrast("pairwise") %>% summary() } else NULL),
        res_over = map(model, ~ if (!is.null(.x)) { emmeans(.x, ~ Treatment_Type | Baseline_BMI3, at = list(Baseline_BMI3 = "Overweight")) %>% contrast("pairwise") %>% summary() } else NULL),
        res_obese = map(model, ~ if (!is.null(.x)) { emmeans(.x, ~ Treatment_Type | Baseline_BMI3, at = list(Baseline_BMI3 = "Obesity")) %>% contrast("pairwise") %>% summary() } else NULL)
      )

    # extract into data frame
    extract_res <- function(lst, colname) {
      map_df(lst, ~ if (is.null(.x)) tibble(estimate = NA_real_, p.value = NA_real_) else as_tibble(.x)[1, c("estimate","p.value")]) %>%
        rename_with(~ paste0(colname,"_",.), everything())
    }
    # assembly
    ancova_out <- ancova_genus_bmi %>%
      transmute(ASV,
                normal_est = map_dbl(res_normal, ~ if (is.null(.x)) NA_real_ else .x$estimate[1]),
                normal_p = map_dbl(res_normal, ~ if (is.null(.x)) NA_real_ else .x$p.value[1]),
                over_est = map_dbl(res_over, ~ if (is.null(.x)) NA_real_ else .x$estimate[1]),
                over_p = map_dbl(res_over, ~ if (is.null(.x)) NA_real_ else .x$p.value[1]),
                obese_est = map_dbl(res_obese, ~ if (is.null(.x)) NA_real_ else .x$estimate[1]),
                obese_p = map_dbl(res_obese, ~ if (is.null(.x)) NA_real_ else .x$p.value[1])
                ) %>%
      mutate(normal_padj = p.adjust(normal_p, method = "BH"),
             over_padj = p.adjust(over_p, method = "BH"),
             obese_padj = p.adjust(obese_p, method = "BH")) %>%
      left_join(data.frame(tax_table(pspairgenus_clr)) %>% rownames_to_column("ASV"), by = "ASV")

    write_csv(ancova_out, file.path(outdir, "ancova_genus_by_BMI.csv"))
  }

} # end BMI

#Enterotype (two groups ONLY)
{
  outdir <- file.path(root_out, "enterotype")

  meta <- data.frame(sample_data(pspairgenus)) %>% rownames_to_column("id")
  # Check column name: Baseline_enterotype_BIC or Enterotype
  et_col <- if ("Baseline_enterotype_BIC" %in% colnames(meta)) "Baseline_enterotype_BIC" else if ("Enterotype" %in% colnames(meta)) "Enterotype" else stop("Enterotype column not found in sample_data.")

  # keep only two enterotypes: recode if needed; drop samples not in the two chosen groups
  meta <- meta %>% mutate(!!et_col := as.character(.data[[et_col]]))
  # If there is "Enterotype three", drop it. Keep only one and two.
  keep_levels <- c("Enterotype one","Enterotype two")
  meta <- meta %>% filter(.data[[et_col]] %in% keep_levels)
  # update phyloseq objects to same sample set
  common_ids <- meta$id
  pspair_et <- prune_samples(common_ids, pspair)
  pspairgenus_et <- prune_samples(common_ids, pspairgenus)
  # run PERMANOVA similar to overall but only samples kept
  meta_local <- meta %>% mutate(Baseline_enterotype_BIC = factor(.data[[et_col]], levels = keep_levels),
                                Treatment_Type = factor(Treatment_Type, levels = c("Placebo","Intervention")),
                                Time = factor(Time, levels = c("Baseline","24-wks")),
                                subject = factor(subject))
  # distances
  dist_bray_genus <- phyloseq::distance(pspairgenus_et, method = "bray")
  ad_bray_genus <- adonis2(dist_bray_genus ~ Baseline_enterotype_BIC * Treatment_Type * Time, data = meta_local, permutations = 999, strata = meta_local$subject)
  write_csv(as.data.frame(ad_bray_genus), file.path(outdir, "permanova_bray_genus_enterotype2.csv"))

  # ASV Bray
  dist_bray_asv <- phyloseq::distance(pspair_et, method = "bray")
  ad_bray_asv <- adonis2(dist_bray_asv ~ Baseline_enterotype_BIC * Treatment_Type * Time, data = meta_local, permutations = 999, strata = meta_local$subject)
  write_csv(as.data.frame(ad_bray_asv), file.path(outdir, "permanova_bray_asv_enterotype2.csv"))

  # other distances
  dist_jacc <- phyloseq::distance(pspair_et, method = "jaccard", binary = TRUE)
  write_csv(as.data.frame(adonis2(dist_jacc ~ Baseline_enterotype_BIC * Treatment_Type * Time, data = meta_local, permutations = 999, strata = meta_local$subject)), file.path(outdir, "permanova_jaccard_enterotype2.csv"))
  dist_wunifrac <- phyloseq::distance(pspair_et, method = "wunifrac")
  write_csv(as.data.frame(adonis2(dist_wunifrac ~ Baseline_enterotype_BIC * Treatment_Type * Time, data = meta_local, permutations = 999, strata = meta_local$subject)), file.path(outdir, "permanova_wunifrac_enterotype2.csv"))
  dist_unifrac <- phyloseq::distance(pspair_et, method = "unifrac")
  write_csv(as.data.frame(adonis2(dist_unifrac ~ Baseline_enterotype_BIC * Treatment_Type * Time, data = meta_local, permutations = 999, strata = meta_local$subject)), file.path(outdir, "permanova_unifrac_enterotype2.csv"))

  # ANCOVA per genus with enterotype stratification (baseline-adjusted)
  # Reuse fulldata_long (ensure it includes enterotype) — prepare fulldata_long_et
  meta_counts <- data.frame(sample_data(pspairgenus_clr)) %>% rownames_to_column("id")
  genus_counts <- as.data.frame(t(otu_table(pspairgenus_clr))) %>% rownames_to_column("id")
  fulldata_all <- left_join(meta_counts, genus_counts, by = "id")
  fulldata_long_et <- fulldata_all %>%
    pivot_longer(cols = starts_with("ASV"), names_to = "ASV", values_to = "value") %>%
    mutate(Time = factor(Time, levels = c("Baseline","24-wks")),
           Baseline_enterotype_BIC = factor(.data[[et_col]], levels = keep_levels),
           Treatment_Type = factor(Treatment_Type, levels = c("Placebo","Intervention")),
           subject = factor(subject))

  genus_wide_et <- fulldata_long_et %>%
    filter(Time %in% c("Baseline","24-wks")) %>%
    dplyr::select(subject, Treatment_Type, Baseline_enterotype_BIC, Time, ASV, value) %>%
    pivot_wider(names_from = Time, values_from = value) %>%
    rename(Baseline = `Baseline`, Week24 = `24-wks`) %>%
    drop_na(Baseline, Week24) %>%
    mutate(Treatment_Type = factor(Treatment_Type, levels = c("Placebo","Intervention")),
           Baseline_enterotype_BIC = factor(Baseline_enterotype_BIC, levels = keep_levels))

  ancova_genus_et <- genus_wide_et %>%
    group_by(ASV) %>%
    nest() %>%
    mutate(model = map(data, ~ tryCatch(lm(Week24 ~ Treatment_Type * Baseline_enterotype_BIC + Baseline, data = .x), error = function(e) NULL))) %>%
    mutate(
      res_e1 = map(model, ~ if (!is.null(.x)) { emmeans(.x, ~ Treatment_Type | Baseline_enterotype_BIC, at = list(Baseline_enterotype_BIC = "Enterotype one")) %>% contrast("pairwise") %>% summary() } else NULL),
      res_e2 = map(model, ~ if (!is.null(.x)) { emmeans(.x, ~ Treatment_Type | Baseline_enterotype_BIC, at = list(Baseline_enterotype_BIC = "Enterotype two")) %>% contrast("pairwise") %>% summary() } else NULL)
    )

  ancova_out_et <- ancova_genus_et %>% transmute(
    ASV,
    e1_est = map_dbl(res_e1, ~ if (is.null(.x)) NA_real_ else .x$estimate[1]),
    e1_p = map_dbl(res_e1, ~ if (is.null(.x)) NA_real_ else .x$p.value[1]),
    e2_est = map_dbl(res_e2, ~ if (is.null(.x)) NA_real_ else .x$estimate[1]),
    e2_p = map_dbl(res_e2, ~ if (is.null(.x)) NA_real_ else .x$p.value[1])
  ) %>%
    mutate(e1_padj = p.adjust(e1_p, method = "BH"), e2_padj = p.adjust(e2_p, method = "BH")) %>%
    left_join(data.frame(tax_table(pspairgenus_clr)) %>% rownames_to_column("ASV"), by = "ASV")

  write_csv(ancova_out_et, file.path(outdir, "ancova_genus_enterotype2.csv"))

} # end enterotype

#Age quartile 
# Use Age_quartile levels "60-64_years", "65-68_years", "68-73_years", "73-80_years"
{
  outdir <- file.path(root_out, "agequartile")

  meta <- data.frame(sample_data(pspair)) %>% rownames_to_column("id")
  if (!"Age_quartile" %in% colnames(meta)) stop("Age_quartile not found in sample data.")
  meta <- meta %>% mutate(Age_quartile = factor(Age_quartile, levels = c("60-64_years","65-68_years","68-73_years","73-80_years")))

  age_groups <- levels(meta$Age_quartile)

  # PERMANOVA stratified by age
  for (age in age_groups) {
    meta_age <- meta %>% filter(Age_quartile == age)
    ids <- meta_age$id
    if (length(ids) < 3) {
      write_lines("insufficient_samples", file.path(outdir, paste0("permanova_age_", age, ".txt")))
      next
    }
    ps_a_asv <- prune_samples(ids, pspair)
    ps_a_genus <- prune_samples(ids, pspairgenus)
    meta_local <- data.frame(sample_data(ps_a_asv)) %>% rownames_to_column("id")
    # distances and adonis2
    write_csv(as.data.frame(adonis2(phyloseq::distance(ps_a_genus, method = "bray") ~ Treatment_Type * Time, data = meta_local, permutations = 999, strata = meta_local$subject)), file.path(outdir, paste0("permanova_bray_genus_age_", age, ".csv")))
    write_csv(as.data.frame(adonis2(phyloseq::distance(ps_a_asv, method = "bray") ~ Treatment_Type * Time, data = meta_local, permutations = 999, strata = meta_local$subject)), file.path(outdir, paste0("permanova_bray_asv_age_", age, ".csv")))
    write_csv(as.data.frame(adonis2(phyloseq::distance(ps_a_asv, method = "jaccard", binary = TRUE) ~ Treatment_Type * Time, data = meta_local, permutations = 999, strata = meta_local$subject)), file.path(outdir, paste0("permanova_jaccard_age_", age, ".csv")))
    write_csv(as.data.frame(adonis2(phyloseq::distance(ps_a_asv, method = "wunifrac") ~ Treatment_Type * Time, data = meta_local, permutations = 999, strata = meta_local$subject)), file.path(outdir, paste0("permanova_wunifrac_age_", age, ".csv")))
    write_csv(as.data.frame(adonis2(phyloseq::distance(ps_a_asv, method = "unifrac") ~ Treatment_Type * Time, data = meta_local, permutations = 999, strata = meta_local$subject)), file.path(outdir, paste0("permanova_unifrac_age_", age, ".csv")))
  }

  # ANCOVA per-genus by age quartile (baseline-adjusted) similar to age section in original code
  if (exists("fulldata_long")) {
    genus_wide_age <- fulldata_long %>%
      filter(Time %in% c("Baseline","24-wks")) %>%
      dplyr::select(subject, Treatment_Type, Time, ASV, value, Age_quartile) %>%
      pivot_wider(names_from = Time, values_from = value) %>%
      dplyr::rename(Baseline = `Baseline`, Week24 = `24-wks`) %>%
      drop_na(Baseline, Week24) %>%
      mutate(Treatment_Type = factor(Treatment_Type, levels = c("Placebo", "Intervention")),
             Age_quartile = factor(Age_quartile, levels = c("60-64_years","65-68_years","68-73_years","73-80_years")))

    ancova_genus_age <- genus_wide_age %>%
      group_by(ASV) %>%
      nest() %>%
      mutate(model = map(data, ~ tryCatch(lm(Week24 ~ Treatment_Type * Age_quartile + Baseline, data = .x), error = function(e) NULL))) %>%
      mutate(
        a1 = map(model, ~ if (!is.null(.x)) emmeans(.x, ~ Treatment_Type | Age_quartile, at = list(Age_quartile = "60-64_years")) %>% contrast("pairwise") %>% summary() else NULL),
        a2 = map(model, ~ if (!is.null(.x)) emmeans(.x, ~ Treatment_Type | Age_quartile, at = list(Age_quartile = "65-68_years")) %>% contrast("pairwise") %>% summary() else NULL),
        a3 = map(model, ~ if (!is.null(.x)) emmeans(.x, ~ Treatment_Type | Age_quartile, at = list(Age_quartile = "68-73_years")) %>% contrast("pairwise") %>% summary() else NULL),
        a4 = map(model, ~ if (!is.null(.x)) emmeans(.x, ~ Treatment_Type | Age_quartile, at = list(Age_quartile = "73-80_years")) %>% contrast("pairwise") %>% summary() else NULL)
      )

    ancova_out_age <- ancova_genus_age %>% transmute(
      ASV,
      a1_est = map_dbl(a1, ~ if (is.null(.x)) NA_real_ else .x$estimate[1]),
      a1_p = map_dbl(a1, ~ if (is.null(.x)) NA_real_ else .x$p.value[1]),
      a2_est = map_dbl(a2, ~ if (is.null(.x)) NA_real_ else .x$estimate[1]),
      a2_p = map_dbl(a2, ~ if (is.null(.x)) NA_real_ else .x$p.value[1]),
      a3_est = map_dbl(a3, ~ if (is.null(.x)) NA_real_ else .x$estimate[1]),
      a3_p = map_dbl(a3, ~ if (is.null(.x)) NA_real_ else .x$p.value[1]),
      a4_est = map_dbl(a4, ~ if (is.null(.x)) NA_real_ else .x$estimate[1]),
      a4_p = map_dbl(a4, ~ if (is.null(.x)) NA_real_ else .x$p.value[1])
    ) %>%
      mutate(a1_padj = p.adjust(a1_p, method = "BH"),
             a2_padj = p.adjust(a2_p, method = "BH"),
             a3_padj = p.adjust(a3_p, method = "BH"),
             a4_padj = p.adjust(a4_p, method = "BH")) %>%
      left_join(data.frame(tax_table(pspairgenus_clr)) %>% rownames_to_column("ASV"), by = "ASV")

    write_csv(ancova_out_age, file.path(outdir, "ancova_genus_by_agequartile.csv"))
  }

} # end age quartile

# Plotting and Visualization
########################################

# Overall fold-change plot
ancova_genus_fc <- ancova_genus %>%
  filter(p.adj < 0.05) %>% 
  mutate(
    fold_change = exp(estimate),
    fc_low = exp(estimate - 1.96 * std.error),
    fc_high = exp(estimate + 1.96 * std.error),
    log2FC = estimate / log(2),
    sig = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      p.adj < 0.10  ~ "•",
      TRUE          ~ ""
    ),
    genus_label = coalesce(as.character(Genus), ASV)
  )

readr::write_csv(ancova_genus_fc, "output_results/ancova_genus_with_foldchange.csv")

ggplot(heat_df, aes(x = "Intervention vs Placebo", y = genus_label, fill = estimate)) +
  geom_tile(color = "grey90", size = 0.1) +
  geom_text(aes(label = sig), size = 2.5, color = "black", vjust = 0.5) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick3", midpoint = 0,
                       limits = max(abs(heat_df$estimate)) * c(-1,1),
                       breaks = c(-1,0,1), labels = c("-1","0","1")) +
  facet_grid(Phylum ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = "") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    strip.background = element_blank(),
    axis.text.y = element_text(face = "italic", size = 8),
    strip.text.y = element_text(angle = 0, size = 8),
    legend.position = "bottom",
    legend.title = element_text(size = 7.5),
    legend.text = element_text(size = 7.5),
    legend.key.size = unit(0.3, "cm")
  ) +
  scale_x_discrete(expand = expansion(mult = 0)) +
  scale_y_discrete(expand = expansion(mult = 0))

ggsave(filename = "output_results/effectofintervention.tiff", height = 4, width = 4, dpi = 600, bg = "white")

# Enterotype heatmap
enterotypesignficant <-  ancova_genus_enterotype %>%
  # keep both enterotype estimates & p-values
  dplyr::select(ASV, Genus, Phylum, E1_estimate, E1_adj_p, E2_estimate, E2_adj_p) %>%
  pivot_longer(
    cols = starts_with("E"),
    names_to = c("Enterotype", ".value"),
    names_pattern = "(E[12])_(.*)"
  ) %>%
  mutate(
    Enterotype = recode(Enterotype,
                        E1 = "Enterotype one",
                        E2 = "Enterotype two"),
    # Keep estimate only if significant, else NA (so it plots white)
    estimate = ifelse(adj_p < 0.05, estimate, NA),
    sig = case_when(
      adj_p < 0.001 ~ "***",
      adj_p < 0.01  ~ "**",
      adj_p < 0.05  ~ "*",
      TRUE          ~ ""
    ),
    genus_label = Genus
  ) %>%
  # keep only rows where genus was significant in at least one enterotype
  group_by(ASV) %>%
  filter(any(!is.na(estimate))) %>%
  ungroup()

ggplot(enterotypesignficant, aes(x = Enterotype, y = genus_label, fill = estimate)) +
  geom_tile(size = 0.1) +   # no `color=...` → no boxes
  geom_text(aes(label = sig), size = 2.5, color = "black") +
  scale_fill_gradient2(
    name = "Estimate",
    low = "navy", mid = "white", high = "firebrick3",
    midpoint = 0,
    limits = max(abs(enterotypesignficant$estimate), na.rm = TRUE) * c(-1, 1),
    breaks = c(-1, 0, 1),
    labels = c("-1", "0", "1"),
    na.value = "white"
  ) +
  facet_grid(Phylum ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = "") +
  theme_classic(base_size = 11) +
  theme(
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    strip.background = element_blank(),
    axis.text.y = element_text(face = "italic", size = 8),
    strip.text.y = element_text(angle = 0, size = 8),
    legend.position = "bottom",
    legend.title = element_text(size = 7.5),
    legend.text = element_text(size = 7.5),
    legend.key.size = unit(0.3, "cm")
  ) + scale_x_discrete(expand = expansion(mult = 0)) +
  scale_y_discrete(expand = expansion(mult = 0))


ggsave(filename = "output_results/effectofintervention_enterotype.tiff", height = 4.5, width = 5, dpi = 600, bg = "white")

# Filter significant genera across age quartiles
# Reshape age ANCOVA results to long format
age_quartile_significant <- ancova_genus_age %>%
  pivot_longer(
    cols = matches("A[1-4]_(estimate|adj_p)"),
    names_to = c("Quartile", ".value"),
    names_pattern = "(A[1-4])_(.*)"
  ) %>%
  mutate(
    Quartile = recode(Quartile,
                      A1 = "60-64_years",
                      A2 = "65-68_years",
                      A3 = "68-73_years",
                      A4 = "73-80_years"),
    # Keep estimate only if significant
    estimate = ifelse(adj_p < 0.05, estimate, NA),
    sig = case_when(
      adj_p < 0.001 ~ "***",
      adj_p < 0.01  ~ "**",
      adj_p < 0.05  ~ "*",
      TRUE          ~ ""
    ),
    genus_label = Genus
  ) %>%
  # Keep only ASVs significant in at least one quartile
  group_by(ASV) %>%
  filter(any(!is.na(estimate))) %>%
  ungroup()

ggplot(age_quartile_significant, aes(x = Quartile, y = genus_label, fill = estimate)) +
  geom_tile(size = 0.1) +   # no `color=...` → no boxes
  geom_text(aes(label = sig), size = 2.5, color = "black") +
  scale_fill_gradient2(
    name = "Estimate",
    low = "navy", mid = "white", high = "firebrick3",
    midpoint = 0,
    limits = max(abs(enterotypesignficant$estimate), na.rm = TRUE) * c(-1, 1),
    breaks = c(-1, 0, 1),
    labels = c("-1", "0", "1"),
    na.value = "white"
  ) +
  facet_grid(Phylum ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = "") +
  theme_classic(base_size = 11) +
  theme(
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    strip.background = element_blank(),
    axis.text.y = element_text(face = "italic", size = 8),
    strip.text.y = element_text(angle = 0, size = 8),
    legend.position = "bottom",
    legend.title = element_text(size = 7.5),
    legend.text = element_text(size = 7.5),
    axis.text.x = element_text(size = 7.5),
    legend.key.size = unit(0.3, "cm")
  ) + scale_x_discrete(expand = expansion(mult = 0)) +
  scale_y_discrete(expand = expansion(mult = 0))

ggsave(filename = "output_results/effectofintervention_age.tiff", height = 4.8, width = 5.5, dpi = 600, bg = "white")

# BMI heatmap
# Pivot longer for plotting
bmi_significant <- ancova_genus_bmi %>%
  pivot_longer(
    cols = matches("(Normal|Overweight|Obesity)_(estimate|adj_p)"),
    names_to = c("BMI_group", ".value"),
    names_pattern = "(Normal|Overweight|Obesity)_(.*)"
  ) %>%
  mutate(
    # Keep estimate only if significant
    estimate = ifelse(adj_p < 0.05, estimate, NA),
    sig = case_when(
      adj_p < 0.001 ~ "***",
      adj_p < 0.01  ~ "**",
      adj_p < 0.05  ~ "*",
      TRUE          ~ ""
    ),
    genus_label = Genus
  ) %>%
  # Keep only ASVs significant in at least one BMI group
  group_by(ASV) %>%
  filter(any(!is.na(estimate))) %>%
  ungroup()

ggplot(bmi_significant, aes(x = BMI_group, y = genus_label, fill = estimate)) +
  geom_tile(size = 0.1) +
  geom_text(aes(label = sig), size = 2.5, color = "black") +
  scale_fill_gradient2(
    name = "Estimate",
    low = "navy", mid = "white", high = "firebrick3",
    midpoint = 0,
    limits = max(abs(bmi_significant$estimate), na.rm = TRUE) * c(-1, 1),
    breaks = c(-1, 0, 1),
    labels = c("-1", "0", "1"),
    na.value = "white"
  ) +
  facet_grid(Phylum ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = "") +
  theme_classic(base_size = 11) +
  theme(
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    strip.background = element_blank(),
    axis.text.y = element_text(face = "italic", size = 8),
    strip.text.y = element_text(angle = 0, size = 8),
    legend.position = "bottom",
    legend.title = element_text(size = 7.5),
    legend.text = element_text(size = 7.5),
    legend.key.size = unit(0.3, "cm")
  ) + scale_x_discrete(expand = expansion(mult = 0)) +
  scale_y_discrete(expand = expansion(mult = 0))
# Save plot
ggsave(filename = "output_results/effectofintervention_BMI3.tiff",
       height = 4.5, width = 5, dpi = 600, bg = "white")