# Assess effect availability, recover paired precision, and describe scientific
# compatibility and uncertainty under assumed within-person correlations.
# These calculations do not fit pooled models.

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
data_dir <- file.path(root, "meta_analysis", "data")
output_dir <- file.path(root, "meta_analysis", "outputs", "evidence_assessment")
table_dir <- file.path(output_dir, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

read_table <- function(name) {
  read.csv(
    file.path(data_dir, name), stringsAsFactors = FALSE, check.names = FALSE,
    na.strings = c(""), encoding = "UTF-8"
  )
}

studies <- read_table("studies.csv")
conditions <- read_table("conditions.csv")
outcomes <- read_table("outcomes.csv")

numeric_columns <- c("window_min", "timepoint_min", "n", "value", "dispersion_value", "ci_low", "ci_high")
for (column in intersect(numeric_columns, names(outcomes))) {
  outcomes[[column]] <- suppressWarnings(as.numeric(outcomes[[column]]))
}

stopifnot(
  length(unique(studies$study_id)) == 18L,
  nrow(studies) == 20L,
  nrow(conditions) == 60L,
  nrow(outcomes) == 263L,
  all(unique(outcomes$study_id) %in% unique(studies$study_id)),
  all(unique(conditions$population_id) %in% unique(studies$population_id)),
  all(unique(outcomes$population_id) %in% unique(studies$population_id)),
  all(unique(outcomes$condition_id[!is.na(outcomes$condition_id)]) %in% conditions$condition_id)
)

se_from_ci <- function(lower, upper, level = 0.95) {
  z <- qnorm(1 - (1 - level) / 2)
  (upper - lower) / (2 * z)
}

paired_from_exact_p <- function(mean_difference, n, two_sided_p) {
  t_value <- qt(1 - two_sided_p / 2, df = n - 1)
  se_difference <- abs(mean_difference) / t_value
  sd_difference <- se_difference * sqrt(n)
  c(t_value = t_value, sd_difference = sd_difference, se_difference = se_difference)
}

paired_se_upper_from_p_bound <- function(mean_difference, n, two_sided_p_upper) {
  # For p < p_upper, |t| exceeds this boundary and the paired SE is smaller
  # than the conservative upper bound returned here.
  t_boundary <- qt(1 - two_sided_p_upper / 2, df = n - 1)
  se_upper <- abs(mean_difference) / t_boundary
  c(t_boundary = t_boundary, se_upper = se_upper)
}

# Directly reported or exactly recoverable later-minus-earlier effects.
bandin <- subset(outcomes, study_id == "Bandin_2015" & metric == "mean difference")
bo <- subset(outcomes, study_id == "Bo_2015" & metric == "adjusted mean difference")
gu <- subset(outcomes, study_id == "Gu_2020" & metric == "adjusted mean difference")
stutz <- subset(outcomes, study_id == "Stutz_2024" & metric == "2-h iAUC model contrast")
lopez <- subset(outcomes, study_id == "Lopez_Minguez_2018" & contrast_id == "Lopez_all" & metric == "mean difference")
enomoto <- subset(outcomes, study_id == "Enomoto_2026" & contrast_id == "Enomoto_day2" &
                   result_form == "within_person_contrast" & data_status == "author_provided")
stopifnot(nrow(bandin) == 1L, nrow(bo) == 1L, nrow(gu) == 1L, nrow(stutz) == 2L, nrow(lopez) == 1L,
          nrow(enomoto) == 1L)

bandin_df <- data.frame(
  study_effect = "Bandin 2015", study_id = bandin$study_id, population_id = bandin$population_id,
  construct = "same meal occasion: lunch", outcome = "glucose AUC above baseline, 150 min",
  later_minus_earlier = bandin$value, se = bandin$dispersion_value,
  ci_low = bandin$value - qt(0.975, df = bandin$n - 1) * bandin$dispersion_value,
  ci_high = bandin$value + qt(0.975, df = bandin$n - 1) * bandin$dispersion_value,
  unit = bandin$units_original, precision_source = "derived from exact paired p=0.002",
  stringsAsFactors = FALSE
)

bo_df <- data.frame(
  study_effect = "Bo 2015", study_id = bo$study_id, population_id = bo$population_id,
  construct = "daypart: morning versus evening", outcome = "glucose AUC, 180 min",
  later_minus_earlier = -bo$value, se = se_from_ci(bo$ci_low, bo$ci_high),
  ci_low = -bo$ci_high, ci_high = -bo$ci_low, unit = bo$units_original,
  precision_source = "reported adjusted contrast and 95% CI; sign reversed",
  stringsAsFactors = FALSE
)

gu_df <- data.frame(
  study_effect = "Gu 2020", study_id = gu$study_id, population_id = gu$population_id,
  construct = "same meal occasion: dinner, extended window", outcome = "glucose total AUC, 240 min",
  later_minus_earlier = gu$value, se = se_from_ci(gu$ci_low, gu$ci_high),
  ci_low = gu$ci_low, ci_high = gu$ci_high, unit = gu$units_original,
  precision_source = "reported adjusted contrast and 95% CI", stringsAsFactors = FALSE
)

stutz_df <- data.frame(
  study_effect = ifelse(stutz$population_id == "Stutz_2024_early", "Stutz 2024: early chronotype", "Stutz 2024: late chronotype"),
  study_id = stutz$study_id, population_id = stutz$population_id,
  construct = "daypart: morning versus evening", outcome = "CGM glucose iAUC, 120 min",
  later_minus_earlier = stutz$value, se = se_from_ci(stutz$ci_low, stutz$ci_high),
  ci_low = stutz$ci_low, ci_high = stutz$ci_high, unit = stutz$units_original,
  precision_source = "reported multilevel-model contrast and 95% CI", stringsAsFactors = FALSE
)

lopez_precision <- paired_from_exact_p(lopez$value, lopez$n, 0.004)
lopez_df <- data.frame(
  study_effect = "Lopez-Minguez 2018: total sample", study_id = lopez$study_id,
  population_id = lopez$population_id,
  construct = "relative to sleep: early versus late mixed dinner", outcome = "glucose total AUC, 120 min",
  later_minus_earlier = lopez$value, se = unname(lopez_precision["se_difference"]),
  ci_low = lopez$value - qt(0.975, df = lopez$n - 1) * unname(lopez_precision["se_difference"]),
  ci_high = lopez$value + qt(0.975, df = lopez$n - 1) * unname(lopez_precision["se_difference"]),
  unit = lopez$units_original, precision_source = "derived from exact paired p=0.004",
  stringsAsFactors = FALSE
)

enomoto_df <- data.frame(
  study_effect = "Enomoto 2026: Day 2", study_id = enomoto$study_id,
  population_id = enomoto$population_id,
  construct = "relative to sleep: dinner 1 h versus 5 h before bedtime",
  outcome = "CGM dinner glucose iAUC, 180 min",
  later_minus_earlier = enomoto$value,
  se = enomoto$dispersion_value / sqrt(enomoto$n),
  ci_low = enomoto$ci_low, ci_high = enomoto$ci_high,
  unit = enomoto$units_original,
  precision_source = "author-provided paired SD, SE, and 95% CI received 2026-08-31",
  stringsAsFactors = FALSE
)

direct_effects <- rbind(bandin_df, bo_df, gu_df, stutz_df, lopez_df, enomoto_df)
direct_effects$se <- round(direct_effects$se, 4)
direct_effects$ci_low <- round(direct_effects$ci_low, 4)
direct_effects$ci_high <- round(direct_effects$ci_high, 4)
write.csv(direct_effects, file.path(table_dir, "direct_or_recoverable_effects.csv"), row.names = FALSE, na = "")

# Additional exact or conservatively bounded paired precision. Bounds remain
# labeled as bounds and are not treated as exact SEs in any pooled model.
precision_specs <- data.frame(
  study_id = c("Gibbs_2014", "Gibbs_2014", "Haldar_2020", "Haldar_2020", "Nakamura_2021", "Sulaimani_2025", "Sulaimani_2025"),
  comparison = c("low GI: evening minus morning", "high GI: evening minus morning", "high GI: dinner minus breakfast", "low GI: dinner minus breakfast", "dinner day 2: late minus early", "control rice: evening minus morning", "GTE rice: evening minus morning"),
  n = c(10, 10, 34, 34, 12, 14, 14),
  mean_difference = c(167, 133, 284.40, 207.10, 2100, 165, 200),
  p_relation = c("<", "=", "<", "<", "<", "<", "<"),
  p_value_or_bound = c(0.001, 0.003, 0.0001, 0.0001, 0.001, 0.01, 0.01),
  data_status = c("reported table", "reported table", "reported model contrast", "reported model contrast", "graph-digitized mean difference", "graph-digitized mean difference", "graph-digitized mean difference"),
  source = c("Gibbs 2014 Table 2", "Gibbs 2014 Table 2", "Haldar 2020 Table 2", "Haldar 2020 Table 2", "Nakamura 2021 Figure 3B", "Sulaimani 2025 Figure 3B", "Sulaimani 2025 Figure 3B"),
  stringsAsFactors = FALSE
)
precision_rows <- lapply(seq_len(nrow(precision_specs)), function(i) {
  x <- precision_specs[i, ]
  if (x$p_relation == "=") {
    exact <- paired_from_exact_p(x$mean_difference, x$n, x$p_value_or_bound)
    se_value <- unname(exact["se_difference"])
    precision_type <- "exact paired SE derived from reported two-sided p"
  } else {
    bound <- paired_se_upper_from_p_bound(x$mean_difference, x$n, x$p_value_or_bound)
    se_value <- unname(bound["se_upper"])
    precision_type <- "conservative upper bound on paired SE from reported p threshold"
  }
  critical <- qt(0.975, df = x$n - 1)
  data.frame(x, precision_type = precision_type, paired_se_or_upper_bound = se_value,
             conservative_ci_low = x$mean_difference - critical * se_value,
             conservative_ci_high = x$mean_difference + critical * se_value,
             stringsAsFactors = FALSE)
})
paired_precision_recovery <- do.call(rbind, precision_rows)
enomoto_author <- outcomes[
  outcomes$study_id == "Enomoto_2026" & outcomes$result_form == "within_person_contrast" &
    outcomes$data_status == "author_provided",
]
stopifnot(nrow(enomoto_author) == 4L, all(enomoto_author$dispersion_type == "SD (paired difference)"))
enomoto_author_rows <- data.frame(
  study_id = enomoto_author$study_id,
  comparison = paste0(enomoto_author$modifier, ": 1 h minus 5 h before bedtime"),
  n = enomoto_author$n,
  mean_difference = enomoto_author$value,
  p_relation = "=",
  p_value_or_bound = as.numeric(enomoto_author$p_value),
  data_status = "author-provided paired t-test",
  source = enomoto_author$source_provenance,
  precision_type = "author-provided paired SD, SE, and 95% CI",
  paired_se_or_upper_bound = enomoto_author$dispersion_value / sqrt(enomoto_author$n),
  conservative_ci_low = enomoto_author$ci_low,
  conservative_ci_high = enomoto_author$ci_high,
  stringsAsFactors = FALSE
)
paired_precision_recovery <- rbind(paired_precision_recovery, enomoto_author_rows)
paired_precision_recovery[c("paired_se_or_upper_bound", "conservative_ci_low", "conservative_ci_high")] <-
  lapply(paired_precision_recovery[c("paired_se_or_upper_bound", "conservative_ci_low", "conservative_ci_high")], round, 4)
write.csv(paired_precision_recovery, file.path(table_dir, "paired_precision_recovery.csv"), row.names = FALSE, na = "")

# Study-level precision_availability classes. These are evidence-structure judgements, not
# inclusion/exclusion decisions and not risk-of-bias ratings.
precision_availability <- data.frame(
  study_id = c(
    "Bandin_2015", "Bo_2015", "Bravo_Garcia_2024", "Enomoto_2026", "Garaulet_2022",
    "Gibbs_2014", "Gu_2020", "Haldar_2020", "Jarrett_1972", "Lopez_Minguez_2018",
    "Nakamura_2021", "Pizinger_2018", "Saad_2012", "Sato_2011", "Service_1983",
    "Stutz_2024", "Sulaimani_2025", "Yadav_2023"
  ),
  precision_availability_class = c(
    "direct/recoverable paired precision", "direct/recoverable paired precision",
    "reported contrast without usable precision", "direct/recoverable paired precision",
    "figure-derived summaries; covariance missing", "reported p-value precision recovered or bounded",
    "direct/recoverable paired precision", "reported p-value precision recovered or bounded",
    "timed means; AUC contrast precision unavailable", "direct/recoverable paired precision",
    "reported p-value precision recovered or bounded", "condition summaries; covariance/model precision missing",
    "condition summaries; covariance/model precision missing", "condition summaries; covariance/model precision missing",
    "condition summaries; covariance/model precision missing", "direct/recoverable paired precision",
    "reported p-value precision recovered or bounded", "condition summaries; covariance/model precision missing"
  ),
  remaining_limitation = c(
    "unit QA only", "condition means absent but model contrast usable", "SE/CI for breakfast pairwise contrasts",
    "none for prespecified Day-2 contrast; Days 1, 3, and 4 retained as sensitivities", "paired precision and exact figure values",
    "GI-specific paired precision recovered; overall timing contrast and unit QA remain", "none for 240-min model contrast",
    "Overall and GI-specific timing contrasts recovered with p-threshold precision bounds; unit QA remains", "participant-level AUC or covariance across timed measures",
    "unit QA only for total sample", "p-threshold precision bound recovered; exact figure values remain",
    "factorial timing contrast precision and unit QA", "pairwise LS-mean contrast precision", "paired covariance",
    "meal-size estimand, paired covariance, and unit QA", "none for chronotype-specific 120-min contrasts",
    "treatment-specific p-threshold bounds recovered; exact figure values/main timing contrast remain", "pairwise timing precision and unit QA"
  ),
  stringsAsFactors = FALSE
)
stopifnot(setequal(precision_availability$study_id, unique(studies$study_id)))
write.csv(precision_availability, file.path(table_dir, "study_precision_availability.csv"), row.names = FALSE)

# Scientifically coherent hierarchical synthesis map. Participant counts are
# descriptive unique participants, not effective crossover sample sizes. A
# principal family is an organizational/inferential unit, not an automatic
# pooled model.
comparison_groups <- data.frame(
  synthesis_id = c("F1", "F1_D120", "F1_D180_IAUC", "F1_D180_TOTAL", "F2", "F2_180_CGM", "F2_DIRECT", "SAME_EXT", "D300_360_ND", "T2D", "OGTT"),
  analysis_level = c("principal family", "nested stratum", "nested stratum", "nested stratum", "principal family", "nested stratum", "nested stratum", "supportive block", "supportive block", "separate population block", "separate challenge block"),
  parent_family = c(NA, "F1", "F1", "F1", NA, "F2", "F2", NA, NA, NA, NA),
  comparison_question = c(
    "Adults without diabetes: daypart mixed-meal glucose AUC over 120–180 min",
    "Healthy adults: morning versus evening mixed meal, 120-min glucose iAUC",
    "Healthy adults: breakfast versus dinner mixed meal, 180-min glucose iAUC",
    "Healthy adults: earlier-day versus later-day mixed meal, 180-min total/reported AUC",
    "Adults without diabetes: same-meal-occasion/relative-to-sleep mixed-meal glucose AUC over 120–180 min",
    "Earlier versus later dinner, 180-min CGM glucose iAUC",
    "Earlier versus later same meal occasion with directly recoverable precision",
    "Earlier versus later dinner, extended 240–300-min glucose AUC",
    "Adults without diabetes: daypart mixed meals with 300–360-min glucose AUC/iAUC",
    "Adults with type 2 diabetes",
    "Isolated oral glucose challenge evidence"
  ),
  studies = c(
    "Bo_2015; Gibbs_2014; Haldar_2020; Pizinger_2018; Stutz_2024; Sulaimani_2025",
    "Gibbs_2014; Stutz_2024", "Haldar_2020; Sulaimani_2025", "Bo_2015; Pizinger_2018",
    "Bandin_2015; Enomoto_2026; Lopez_Minguez_2018; Nakamura_2021",
    "Enomoto_2026; Nakamura_2021", "Bandin_2015; Lopez_Minguez_2018", "Gu_2020; Sato_2011",
    "Saad_2012; Service_1983; Yadav_2023_ND", "Bravo_Garcia_2024; Yadav_2023_T2D",
    "Garaulet_2022; Jarrett_1972"
  ),
  k_studies = c(6, 2, 2, 2, 4, 2, 2, 2, 3, 2, 2),
  approximate_unique_participants = c(129, 55, 48, 26, 75, 25, 50, 30, 45, 30, 869),
  compatibility = c(
    "Common daypart construct and 2–3-h window; AUC definition, measurement method, modifier structure, and units differ",
    "Same window/iAUC label, but CGM versus plasma and chronotype/GI multi-condition structure",
    "Same window/iAUC label, but both studies have modifier-dependent repeated contrasts",
    "Same window, but AUC definition/unit and factorial structure differ",
    "Common delayed-meal construct and 2–3-h window; lunch/dinner, clock/relative-to-sleep timing, AUC definition, and measurement differ",
    "Same window, CGM, incremental AUC, and dinner occasion; Enomoto has exact author-provided paired precision, while Nakamura remains graph-digitized with a conservative p-threshold SE bound",
    "Both are mixed meals with direct precision, but 150 versus 120 min and AUC-above-baseline versus total AUC differ",
    "Different windows and measurement methods; Gu includes reciprocal snack timing",
    "Long-window supportive evidence; 300 versus 360-min estimands and multi-condition structures differ",
    "Different timing questions: breakfast delay versus breakfast/lunch/dinner daypart",
    "Mixed relative-to-sleep versus clock-time constructs; Jarrett lacks reported AUC precision"
  ),
  synthesis_role = c(
    "Reported individually; total and incremental AUC, measurement methods, and modifiers differ",
    "Reported individually; measurement methods and modifier structure differ",
    "Reported individually; modifier structure and bounded precision limit comparability",
    "Reported individually; AUC definitions and units differ",
    "No family-wide pooled estimate; AUC definitions and windows differ",
    "Primary delayed-dinner synthesis; retain digitization and p-threshold precision limitations",
    "Reported individually; AUC definitions and windows differ",
    "Reported individually; windows and co-interventions differ",
    "Reported individually; estimands and paired precision differ",
    "Reported individually; timing questions differ",
    "Reported individually; challenge timing and precision differ"
  ),
  stringsAsFactors = FALSE
)
write.csv(comparison_groups, file.path(table_dir, "comparison_groups.csv"), row.names = FALSE)

# Audit whether a single family-wide effect representation is both computable
# and scientifically comparable. These criteria explain the synthesis choices.
effect_representation_audit <- data.frame(
  family = c(rep("F1", 6), rep("F2", 4)),
  study_id = c("Bo_2015", "Gibbs_2014", "Haldar_2020", "Pizinger_2018", "Stutz_2024", "Sulaimani_2025", "Bandin_2015", "Enomoto_2026", "Lopez_Minguez_2018", "Nakamura_2021"),
  primary_auc_construct = c("reported glucose AUC, 180 min", "incremental AUC, 120 min", "incremental AUC, 180 min", "total AUC, 180 min", "incremental AUC, 120 min", "incremental AUC, 180 min", "AUC above baseline, 150 min", "incremental AUC, 180 min", "total AUC, 120 min", "incremental AUC, 180 min"),
  raw_mean_difference = c("available as adjusted contrast", "available for GI-specific comparisons", "available for GI-specific comparisons", "available from condition means", "available as chronotype-specific model contrasts", "available from graph-derived condition means", "available with exact paired precision", "author-provided for Days 1–4", "available with exact paired precision", "available from graph-derived condition means"),
  paired_precision = c("reported CI", "exact/bounded from reported p by GI", "bounded from p<0.0001 by GI", "missing", "reported CI by chronotype", "bounded from p<0.01 by treatment", "exact from p=0.002", "author-provided paired SD/SE/CI", "exact from p=0.004", "bounded from p<0.001 for dinner"),
  log_ratio_of_means = c("not computable: condition means absent", "computable", "computable", "computable", "computable", "computable from graph-derived means", "computable", "computable", "computable", "computable from graph-derived means"),
  scientific_compatibility_issue = c("AUC definition and printed units uncertain", "incremental plasma AUC and unequal fasting durations", "incremental plasma AUC with GI factorial structure", "total AUC with meal/sleep factorial structure", "CGM incremental AUC with chronotype strata", "incremental plasma AUC with GTE factorial structure", "above-baseline lunch AUC", "CGM incremental dinner AUC across repeated days", "total dinner AUC", "CGM incremental dinner AUC from figure"),
  family_wide_role = "structured display; nested-stratum estimate only when compatible",
  stringsAsFactors = FALSE
)
write.csv(effect_representation_audit, file.path(table_dir, "effect_representation_audit.csv"), row.names = FALSE, na = "")

pooling_compatibility <- data.frame(
  family = c("F1", "F2"),
  k = c(6, 4),
  common_effect_computable_for_all = c("No: Bo lacks condition means for a log ratio; paired standardized effects remain covariance/design dependent", "Technically for log ratios, but not for a scientifically common AUC construct"),
  construct_compatibility = c("Fail: total/reported AUC and incremental AUC, plasma and CGM, and factorial modifiers are mixed", "Fail: total AUC, above-baseline AUC, and incremental AUC are not interchangeable proportional constructs"),
  precision_compatibility = c("Fail for an omnibus effect: several overall paired/model precisions remain unavailable", "Fail for an omnibus effect: family-wide estimates still mix direct, exact, and bounded precision; the Enomoto–Nakamura nested stratum is now estimable"),
  decision = c("No daypart family-wide pooled estimate; report studies individually", "No delayed-meal family-wide pooled estimate; pool only the compatible Enomoto–Nakamura comparison"),
  stringsAsFactors = FALSE
)
write.csv(pooling_compatibility, file.path(table_dir, "pooling_compatibility.csv"), row.names = FALSE, na = "")

# Explicit one-destination map for every extracted study/population row. This
# guards against smaller homogeneous strata making the remaining eligible
# evidence appear to have vanished.
evidence_destinations <- data.frame(
  population_id = c(
    "Bandin_2015_all", "Bo_2015_all", "Bravo_Garcia_2024_T2D", "Enomoto_2026_all",
    "Garaulet_2022_all", "Gibbs_2014_all", "Gu_2020_all", "Haldar_2020_all",
    "Jarrett_1972_all", "Lopez_Minguez_2018_all", "Nakamura_2021_all", "Pizinger_2018_all",
    "Saad_2012_all", "Sato_2011_all", "Service_1983_all", "Stutz_2024_early",
    "Stutz_2024_late", "Sulaimani_2025_all", "Yadav_2023_ND", "Yadav_2023_T2D"
  ),
  primary_destination = c(
    "F2", "F1", "T2D", "F2", "OGTT", "F1", "SAME_EXT", "F1", "OGTT", "F2",
    "F2", "F1", "D300_360_ND", "SAME_EXT", "D300_360_ND", "F1", "F1", "F1",
    "D300_360_ND", "T2D"
  ),
  destination_label = c(
    "Same-meal occasion, 120–180 min", "Daypart mixed meal, 120–180 min",
    "Type 2 diabetes, separate timing questions", "Same-meal occasion, 120–180 min",
    "Isolated oral glucose challenge", "Daypart mixed meal, 120–180 min",
    "Extended-window same dinner, 240–300 min", "Daypart mixed meal, 120–180 min",
    "Isolated oral glucose challenge", "Same-meal occasion, 120–180 min",
    "Same-meal occasion, 120–180 min", "Daypart mixed meal, 120–180 min",
    "Long-window daypart, 300–360 min, no diabetes", "Extended-window same dinner, 240–300 min",
    "Long-window daypart, 300–360 min, no diabetes", "Daypart mixed meal, 120–180 min",
    "Daypart mixed meal, 120–180 min", "Daypart mixed meal, 120–180 min",
    "Long-window daypart, 300–360 min, no diabetes", "Type 2 diabetes, separate timing questions"
  ),
  stringsAsFactors = FALSE
)
evidence_destinations <- merge(
  studies[c("study_id", "population_id", "metabolic_group", "timing_construct", "challenge_type", "n_analyzed")],
  evidence_destinations,
  by = "population_id", all.x = TRUE, sort = FALSE
)
stopifnot(
  nrow(evidence_destinations) == nrow(studies),
  !anyDuplicated(evidence_destinations$population_id),
  !any(is.na(evidence_destinations$primary_destination)),
  setequal(evidence_destinations$population_id, studies$population_id)
)
write.csv(evidence_destinations, file.path(table_dir, "evidence_destination_map.csv"), row.names = FALSE, na = "")

# Correlation-sensitivity calculation for studies with two condition summaries.
# This quantifies recoverability only; it is not a primary variance imputation.
comparison_specs <- data.frame(
  comparison_id = c(
    "Garaulet_OGTT", "Gibbs_low", "Gibbs_high", "Haldar_high", "Haldar_low",
    "Nakamura_dinner_day2", "Pizinger_normal_sleep",
    "Pizinger_late_sleep", "Sato_dinner", "Service_medium", "Service_large",
    "Sulaimani_control", "Sulaimani_GTE", "Yadav_ND_BD", "Yadav_T2D_BD"
  ),
  study_id = c(
    "Garaulet_2022", "Gibbs_2014", "Gibbs_2014", "Haldar_2020", "Haldar_2020",
    "Nakamura_2021", "Pizinger_2018", "Pizinger_2018",
    "Sato_2011", "Service_1983", "Service_1983", "Sulaimani_2025", "Sulaimani_2025",
    "Yadav_2023", "Yadav_2023"
  ),
  population_id = c(
    "Garaulet_2022_all", "Gibbs_2014_all", "Gibbs_2014_all", "Haldar_2020_all",
    "Haldar_2020_all", "Nakamura_2021_all",
    "Pizinger_2018_all", "Pizinger_2018_all", "Sato_2011_all", "Service_1983_all",
    "Service_1983_all", "Sulaimani_2025_all", "Sulaimani_2025_all", "Yadav_2023_ND",
    "Yadav_2023_T2D"
  ),
  metric = c(
    "120-min glucose AUC", "glucose iAUC", "glucose iAUC", "test-meal glucose iAUC",
    "test-meal glucose iAUC", "3-h glucose iAUC",
    "MTT total glucose AUC", "MTT total glucose AUC", "dinner glucose AUC",
    "integrated glucose response", "integrated glucose response", "glucose iAUC", "glucose iAUC",
    "glucose iAUC", "glucose iAUC"
  ),
  modifier = c(
    NA, "low GI", "high GI", "high GI", "low GI", "dinner day 2",
    "normal sleep", "late sleep", NA, "medium meal", "large meal", "control", "GTE",
    "no diabetes", "T2D"
  ),
  early_condition = c(
    "Garaulet_early", "Gibbs_low_AM", "Gibbs_high_AM", "Haldar_high_breakfast",
    "Haldar_low_breakfast", "Nakamura_early", "Piz_NsNm",
    "Piz_LsNm", "Sato_normal", "Service_medium_08", "Service_large_08",
    "Sulaimani_control_AM", "Sulaimani_GTE_AM", "Yadav_2023_ND_B", "Yadav_2023_T2D_B"
  ),
  late_condition = c(
    "Garaulet_late", "Gibbs_low_PM", "Gibbs_high_PM", "Haldar_high_dinner",
    "Haldar_low_dinner", "Nakamura_late", "Piz_NsLm",
    "Piz_LsLm", "Sato_late", "Service_medium_18", "Service_large_18",
    "Sulaimani_control_PM", "Sulaimani_GTE_PM", "Yadav_2023_ND_D", "Yadav_2023_T2D_D"
  ),
  stringsAsFactors = FALSE
)

get_condition_row <- function(spec, condition_id) {
  rows <- outcomes[
    outcomes$study_id == spec$study_id & outcomes$population_id == spec$population_id &
      outcomes$condition_id == condition_id & outcomes$metric == spec$metric,
  ]
  if (!is.na(spec$modifier)) rows <- rows[rows$modifier == spec$modifier, ]
  if (nrow(rows) != 1L) stop("Expected one outcome for ", spec$comparison_id, " / ", condition_id, "; got ", nrow(rows))
  rows
}

dispersion_to_sd <- function(row) {
  if (row$dispersion_type == "SD") return(row$dispersion_value)
  if (row$dispersion_type == "SEM") return(row$dispersion_value * sqrt(row$n))
  stop("Unsupported dispersion type in sensitivity calculation: ", row$dispersion_type)
}

correlations <- c(0.25, 0.50, 0.75)
sensitivity_rows <- list()
counter <- 1L
for (i in seq_len(nrow(comparison_specs))) {
  spec <- comparison_specs[i, ]
  early <- get_condition_row(spec, spec$early_condition)
  late <- get_condition_row(spec, spec$late_condition)
  n_pair <- min(early$n, late$n)
  sd_early <- dispersion_to_sd(early)
  sd_late <- dispersion_to_sd(late)
  md <- late$value - early$value
  for (r in correlations) {
    sd_difference <- sqrt(sd_early^2 + sd_late^2 - 2 * r * sd_early * sd_late)
    se_difference <- sd_difference / sqrt(n_pair)
    critical <- qt(0.975, df = n_pair - 1)
    sensitivity_rows[[counter]] <- data.frame(
      comparison_id = spec$comparison_id, study_id = spec$study_id,
      population_id = spec$population_id, metric = spec$metric, modifier = spec$modifier,
      n = n_pair, later_minus_earlier = md, unit = late$units_original,
      early_sd = sd_early, late_sd = sd_late, assumed_r = r,
      assumed_sd_difference = sd_difference, assumed_se = se_difference,
      ci_low = md - critical * se_difference, ci_high = md + critical * se_difference,
      data_status = paste(unique(c(early$data_status, late$data_status)), collapse = ";"),
      stringsAsFactors = FALSE
    )
    counter <- counter + 1L
  }
}
covariance_sensitivity <- do.call(rbind, sensitivity_rows)
round_cols <- c("later_minus_earlier", "early_sd", "late_sd", "assumed_sd_difference", "assumed_se", "ci_low", "ci_high")
covariance_sensitivity[round_cols] <- lapply(covariance_sensitivity[round_cols], round, 4)
write.csv(covariance_sensitivity, file.path(table_dir, "crossover_correlation_sensitivity.csv"), row.names = FALSE, na = "")

robustness_rows <- lapply(split(covariance_sensitivity, covariance_sensitivity$comparison_id), function(x) {
  excludes_zero <- x$ci_low > 0 | x$ci_high < 0
  classification <- if (all(excludes_zero) && all(x$ci_low > 0)) {
    "positive CI excludes zero across r=0.25–0.75"
  } else if (all(excludes_zero) && all(x$ci_high < 0)) {
    "negative CI excludes zero across r=0.25–0.75"
  } else if (any(excludes_zero)) {
    "CI conclusion changes across assumed correlations"
  } else {
    "CI includes zero across r=0.25–0.75"
  }
  data.frame(
    comparison_id = x$comparison_id[1], study_id = x$study_id[1],
    later_minus_earlier = x$later_minus_earlier[1], unit = x$unit[1],
    robustness_class = classification, data_status = x$data_status[1],
    stringsAsFactors = FALSE
  )
})
correlation_robustness <- do.call(rbind, robustness_rows)
write.csv(correlation_robustness, file.path(table_dir, "correlation_robustness_summary.csv"), row.names = FALSE, na = "")

precision_availability_counts <- as.data.frame(table(precision_availability$precision_availability_class), stringsAsFactors = FALSE)
names(precision_availability_counts) <- c("precision_availability_class", "n_studies")
write.csv(precision_availability_counts, file.path(table_dir, "precision_availability_counts.csv"), row.names = FALSE)

source(file.path(root, "meta_analysis", "record_session.R"))
record_session(file.path(output_dir, "session_info.txt"))

cat("Validated:", length(unique(studies$study_id)), "studies;", nrow(outcomes), "outcome rows\n")
cat("Direct/recoverable study IDs:", length(unique(direct_effects$study_id)), "\n")
cat("Additional exact/bounded paired-precision comparisons:", nrow(paired_precision_recovery), "\n")
cat("Correlation-sensitivity comparisons:", length(unique(covariance_sensitivity$comparison_id)), "\n")
cat("Assumption-stable positive comparisons:", sum(grepl("positive CI", correlation_robustness$robustness_class)), "\n")
cat("No pooled model was fitted.\n")
