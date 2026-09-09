# Generate a complete, source-driven map of the principal glucose result for
# every included study/population. The table is descriptive: rows retain their
# original AUC definitions, windows, and units and are never pooled merely
# because a later-minus-earlier difference can be calculated.

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
data_dir <- file.path(root, "meta_analysis", "data")
output_dir <- file.path(root, "meta_analysis", "outputs", "study_summary")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_table <- function(name) read.csv(
  file.path(data_dir, name), stringsAsFactors = FALSE, check.names = FALSE,
  na.strings = "", encoding = "UTF-8"
)
studies <- read_table("studies.csv")
outcomes <- read_table("outcomes.csv")
for (column in c("window_min", "timepoint_min", "n", "value", "dispersion_value", "ci_low", "ci_high")) {
  outcomes[[column]] <- suppressWarnings(as.numeric(outcomes[[column]]))
}

direct_specs <- data.frame(
  study_id = c("Bandin_2015", "Bo_2015", "Bravo_Garcia_2024", "Enomoto_2026", "Garaulet_2022", "Gu_2020", "Haldar_2020", "Lopez_Minguez_2018", "Stutz_2024", "Stutz_2024"),
  population_id = c("Bandin_2015_all", "Bo_2015_all", "Bravo_Garcia_2024_T2D", "Enomoto_2026_all", "Garaulet_2022_all", "Gu_2020_all", "Haldar_2020_all", "Lopez_Minguez_2018_all", "Stutz_2024_early", "Stutz_2024_late"),
  contrast_id = c("Bandin_lunch", "Bo_daypart", "Bravo_nonexercise", "Enomoto_day2", "Garaulet_OGTT", "Gu_dinner", "Haldar_overall_daypart", "Lopez_all", "Stutz_early_2-h_iAUC", "Stutz_late_2-h_iAUC"),
  metric = c("mean difference", "adjusted mean difference", "CGM glucose iAUC difference", "dinner glucose iAUC paired difference", "relative difference", "adjusted mean difference", "overall model timing contrast", "mean difference", "2-h iAUC model contrast", "2-h iAUC model contrast"),
  multiplier = c(1, -1, 1, 1, 1, 1, 1, 1, 1, 1),
  display_label = c("Late versus early lunch", "Evening versus morning mixed meal", "Later versus early breakfast", "Dinner 1 h versus 5 h before sleep, Day 2", "Late versus early evening OGTT", "Late versus routine dinner", "Dinner versus breakfast, GI conditions combined", "Late versus early dinner", "Evening versus morning, early chronotype", "Evening versus morning, late chronotype"),
  stringsAsFactors = FALSE
)

condition_specs <- data.frame(
  study_id = c(
    "Gibbs_2014", "Gibbs_2014", "Nakamura_2021",
    "Pizinger_2018", "Pizinger_2018", "Saad_2012", "Saad_2012",
    "Sato_2011", "Service_1983", "Service_1983", "Service_1983", "Service_1983",
    "Sulaimani_2025", "Sulaimani_2025",
    "Yadav_2023", "Yadav_2023", "Yadav_2023", "Yadav_2023"
  ),
  population_id = c(
    "Gibbs_2014_all", "Gibbs_2014_all", "Nakamura_2021_all",
    "Pizinger_2018_all", "Pizinger_2018_all", "Saad_2012_all", "Saad_2012_all",
    "Sato_2011_all", rep("Service_1983_all", 4), rep("Sulaimani_2025_all", 2),
    rep("Yadav_2023_ND", 2), rep("Yadav_2023_T2D", 2)
  ),
  early_condition = c(
    "Gibbs_low_AM", "Gibbs_high_AM", "Nakamura_early",
    "Piz_NsNm", "Piz_LsNm", "Saad_breakfast", "Saad_breakfast", "Sato_normal",
    "Service_medium_08", "Service_medium_08", "Service_large_08", "Service_large_08",
    "Sulaimani_control_AM", "Sulaimani_GTE_AM",
    "Yadav_2023_ND_B", "Yadav_2023_ND_B", "Yadav_2023_T2D_B", "Yadav_2023_T2D_B"
  ),
  late_condition = c(
    "Gibbs_low_PM", "Gibbs_high_PM", "Nakamura_late",
    "Piz_NsLm", "Piz_LsLm", "Saad_lunch", "Saad_dinner", "Sato_late",
    "Service_medium_13", "Service_medium_18", "Service_large_13", "Service_large_18",
    "Sulaimani_control_PM", "Sulaimani_GTE_PM",
    "Yadav_2023_ND_L", "Yadav_2023_ND_D", "Yadav_2023_T2D_L", "Yadav_2023_T2D_D"
  ),
  contrast_id = c(
    "Gibbs_low", "Gibbs_high", "Nakamura_dinner_day_2",
    "Pizinger_normal_sleep", "Pizinger_late_sleep", rep("Saad_daypart", 2), "Sato_dinner",
    rep("Service_medium", 2), rep("Service_large", 2), "Sulaimani_control", "Sulaimani_GTE",
    rep("Yadav_2023_ND_daypart", 2), rep("Yadav_2023_T2D_daypart", 2)
  ),
  display_label = c(
    "Evening versus morning, low-GI meal", "Evening versus morning, high-GI meal", "Late versus early dinner, Day 2",
    "Late versus normal meal, normal sleep", "Late versus normal meal, late sleep", "Lunch versus breakfast", "Dinner versus breakfast", "Late versus normal dinner",
    "13:00 versus 08:00, medium meal", "18:00 versus 08:00, medium meal", "13:00 versus 08:00, large meal", "18:00 versus 08:00, large meal",
    "Evening versus morning, control rice", "Evening versus morning, GTE rice",
    "Lunch versus breakfast, no diabetes", "Dinner versus breakfast, no diabetes", "Lunch versus breakfast, T2D", "Dinner versus breakfast, T2D"
  ), stringsAsFactors = FALSE
)

make_direct <- function(spec) {
  x <- outcomes[
    outcomes$study_id == spec$study_id & outcomes$population_id == spec$population_id &
      outcomes$contrast_id == spec$contrast_id & outcomes$metric == spec$metric &
      outcomes$analyte == "glucose" &
      outcomes$result_form != "condition_summary",
  ]
  stopifnot(nrow(x) >= 1L)
  lower <- x$ci_low
  upper <- x$ci_high
  from_se <- x$dispersion_type == "SE" & !is.na(x$dispersion_value)
  lower[from_se] <- x$value[from_se] - qt(0.975, df = x$n[from_se] - 1) * x$dispersion_value[from_se]
  upper[from_se] <- x$value[from_se] + qt(0.975, df = x$n[from_se] - 1) * x$dispersion_value[from_se]
  data.frame(
    study_id = x$study_id, population_id = x$population_id,
    contrast = if (nrow(x) == 1L) spec$display_label else paste(spec$display_label, x$direction),
    timing_construct = NA_character_, challenge_type = NA_character_, metabolic_group = NA_character_,
    metric = x$metric, auc_definition = x$auc_definition, window_min = x$window_min,
    early_value = NA_real_, late_value = NA_real_,
    later_minus_earlier = spec$multiplier * x$value,
    ci_low = ifelse(is.na(lower), NA_real_, ifelse(spec$multiplier == 1, lower, -upper)),
    ci_high = ifelse(is.na(upper), NA_real_, ifelse(spec$multiplier == 1, upper, -lower)),
    p_value = x$p_value, units = x$units_original, result_form = x$result_form,
    data_status = x$data_status, precision_status = ifelse(
      x$data_status == "author_provided", "author-provided paired SD, SE, and 95% CI",
      ifelse(from_se, "95% CI derived from paired SE", ifelse(is.na(lower), "no SE/CI reported for this contrast", "reported 95% CI"))
    ),
    source_provenance = x$source_provenance, stringsAsFactors = FALSE
  )
}

make_condition <- function(spec) {
  x <- outcomes[
    outcomes$study_id == spec$study_id & outcomes$population_id == spec$population_id &
      outcomes$contrast_id == spec$contrast_id & outcomes$condition_id %in% c(spec$early_condition, spec$late_condition) &
      outcomes$analyte == "glucose" & outcomes$result_form == "condition_summary" & !is.na(outcomes$auc_definition),
  ]
  early <- x[x$condition_id == spec$early_condition,]
  late <- x[x$condition_id == spec$late_condition,]
  stopifnot(nrow(early) == 1L, nrow(late) == 1L, early$metric == late$metric, early$units_original == late$units_original)
  data.frame(
    study_id = spec$study_id, population_id = spec$population_id, contrast = spec$display_label,
    timing_construct = NA_character_, challenge_type = NA_character_, metabolic_group = NA_character_,
    metric = early$metric, auc_definition = early$auc_definition, window_min = early$window_min,
    early_value = early$value, late_value = late$value, later_minus_earlier = late$value - early$value,
    ci_low = NA_real_, ci_high = NA_real_, p_value = NA_character_, units = early$units_original,
    result_form = "difference of reported condition means", data_status = paste(unique(c(early$data_status, late$data_status)), collapse = "; "),
    precision_status = "paired/model precision unavailable in the report", source_provenance = paste(unique(c(early$source_provenance, late$source_provenance)), collapse = "; "),
    stringsAsFactors = FALSE
  )
}

direct_specs <- subset(direct_specs, study_id %in% studies$study_id)
condition_specs <- subset(condition_specs, study_id %in% studies$study_id)
direct_rows <- do.call(rbind, lapply(seq_len(nrow(direct_specs)), function(i) make_direct(direct_specs[i,])))
condition_rows <- do.call(rbind, lapply(seq_len(nrow(condition_specs)), function(i) make_condition(condition_specs[i,])))

# Jarrett reports timed glucose means rather than a participant-level AUC. Use
# mean-curve trapezoidal AUC only as a descriptive magnitude; no variance is
# invented for the reconstructed contrast.
jarrett <- outcomes[outcomes$study_id == "Jarrett_1972" & outcomes$analyte == "glucose" & outcomes$domain == "timed_glucose",]
auc <- function(condition) {
  x <- jarrett[jarrett$condition_id == condition,]
  x <- x[order(x$timepoint_min),]
  sum(diff(x$timepoint_min) * (head(x$value, -1) + tail(x$value, -1)) / 2)
}
jarrett_row <- data.frame(
  study_id = "Jarrett_1972", population_id = "Jarrett_1972_all", contrast = "20:00 versus 09:00 OGTT",
  timing_construct = NA_character_, challenge_type = NA_character_, metabolic_group = NA_character_,
  metric = "mean-curve glucose AUC", auc_definition = "trapezoidal from reported timed means", window_min = 120,
  early_value = auc("Jarrett_09"), late_value = auc("Jarrett_20"),
  later_minus_earlier = auc("Jarrett_20") - auc("Jarrett_09"), ci_low = NA_real_, ci_high = NA_real_, p_value = NA_character_,
  units = "mg/dL×min", result_form = "descriptive reconstruction", data_status = "derived",
  precision_status = "participant-level AUC variance unavailable", source_provenance = "Jarrett_1972; PDF pp2-3; Tables II-III", stringsAsFactors = FALSE
)

lu_row <- jarrett_row
lu_row[1,] <- NA
lu_row$study_id <- "Lu_2022"; lu_row$population_id <- "Lu_2022_all"
lu_row$contrast <- "14:00 versus 12:00 lunch, no preload"
lu_row$metric <- "Capillary glucose iAUC"; lu_row$auc_definition <- "incremental; reported window ambiguous"
lu_row$p_value <- "<0.05"; lu_row$units <- "mmol/L×min"
lu_row$result_form <- "narrative_contrast"; lu_row$data_status <- "narrative_only"
lu_row$precision_status <- "No numerical contrast/paired precision; AUC label 0-120 versus sampling -30 to +90 min"
lu_row$source_provenance <- "Lu 2022; pp7-8, Results 3.3 and Figure 3"
results <- rbind(direct_rows, condition_rows, jarrett_row, lu_row)
study_fields <- studies[c("study_id", "population_id", "timing_construct", "challenge_type", "metabolic_group")]
results <- merge(results, study_fields, by = c("study_id", "population_id"), all.x = TRUE, suffixes = c("", ".study"), sort = FALSE)
for (field in c("timing_construct", "challenge_type", "metabolic_group")) results[[field]] <- results[[paste0(field, ".study")]]
results <- results[setdiff(names(results), grep("\\.study$", names(results), value = TRUE))]
results$direction <- ifelse(results$later_minus_earlier > 0, "larger response later", ifelse(results$later_minus_earlier < 0, "smaller response later", "no difference"))
results$direction[results$study_id == "Lu_2022"] <- "larger capillary response later (reported P<0.05); no numeric contrast extracted"
results <- results[order(match(results$population_id, studies$population_id), results$contrast),]

stopifnot(
  setequal(unique(results$population_id), studies$population_id),
  length(unique(results$study_id)) == 17L,
  all(is.finite(results$later_minus_earlier) | results$study_id == "Lu_2022"),
  nrow(results) == 25L
)

write.csv(results, file.path(output_dir, "study_primary_result_map.csv"), row.names = FALSE, na = "")
source(file.path(root, "meta_analysis", "record_session.R"))
record_session(file.path(output_dir, "session_info.txt"))
cat("Published-data study result map complete:", nrow(results), "contrasts across", length(unique(results$study_id)), "studies.\n")
