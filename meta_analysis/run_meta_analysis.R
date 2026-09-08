# Primary meta-analysis for the Meal Timing and Glycemia systematic review.
#
# The only pooled analysis is the delayed-dinner comparison reported by
# Enomoto 2026 and Nakamura 2021. Both studies measured 180-minute CGM glucose
# iAUC after the same dinner was delayed by several hours. Effects are later
# minus earlier; positive values indicate a greater response after later dinner.

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
data_file <- file.path(root, "meta_analysis", "data", "outcomes.csv")
output_dir <- file.path(root, "meta_analysis", "outputs", "primary_meta_analysis")
table_dir <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

outcomes <- read.csv(data_file, stringsAsFactors = FALSE, check.names = FALSE)
mgdl_to_mmoll <- 1 / 18

enomoto <- outcomes[
  outcomes$study_id == "Enomoto_2026" &
    outcomes$result_form == "within_person_contrast" &
    outcomes$data_status == "author_provided",
]
enomoto$day <- as.integer(sub("day ", "", enomoto$modifier, fixed = TRUE))
enomoto$dispersion_value <- as.numeric(enomoto$dispersion_value)
enomoto <- enomoto[order(enomoto$day),]
stopifnot(
  nrow(enomoto) == 4L,
  all(enomoto$dispersion_type == "SD (paired difference)"),
  all(enomoto$units_original == "mg/dL×min")
)

nakamura_conditions <- outcomes[
  outcomes$study_id == "Nakamura_2021" &
    outcomes$contrast_id == "Nakamura_dinner_day_2" &
    outcomes$result_form == "condition_summary",
]
nakamura_conditions <- nakamura_conditions[
  match(c("earlier", "later"), nakamura_conditions$timing_level),
]
stopifnot(
  nrow(nakamura_conditions) == 2L,
  !anyNA(nakamura_conditions$timing_level),
  all(nakamura_conditions$n == 12L),
  all(nakamura_conditions$data_status == "graph_digitized"),
  all(nakamura_conditions$units_original == "mg/3hour/dL (as figure)")
)
nakamura_effect_mg <- diff(nakamura_conditions$value)
nakamura_n <- unique(nakamura_conditions$n)
stopifnot(length(nakamura_n) == 1L, nakamura_effect_mg == 2100)

# Nakamura reports P < 0.001 for this Day 2 dinner iAUC contrast and states
# that early and late dinner conditions were compared with paired t tests.
# Treating a threshold as equality gives the largest SE compatible with it.
nakamura_se_from_p <- function(p_value) {
  abs(nakamura_effect_mg) / qt(1 - p_value / 2, df = nakamura_n - 1)
}

make_study_data <- function(enomoto_day = 2L, nakamura_p = 0.001) {
  e <- enomoto[enomoto$day == enomoto_day,]
  stopifnot(nrow(e) == 1L)
  data.frame(
    study = c("Enomoto 2026", "Nakamura 2021"),
    effect = c(e$value, nakamura_effect_mg) * mgdl_to_mmoll,
    se = c(e$dispersion_value / sqrt(e$n), nakamura_se_from_p(nakamura_p)) * mgdl_to_mmoll,
    ci_low = c(e$ci_low, nakamura_effect_mg - qt(0.975, df = nakamura_n - 1) * nakamura_se_from_p(nakamura_p)) * mgdl_to_mmoll,
    ci_high = c(e$ci_high, nakamura_effect_mg + qt(0.975, df = nakamura_n - 1) * nakamura_se_from_p(nakamura_p)) * mgdl_to_mmoll,
    source = c(
      paste0("Author-provided paired effect and 95% CI; Day ", enomoto_day),
      paste0("Digitized from figure; SE derived by treating P < 0.001 as P = ", format(nakamura_p, scientific = FALSE))
    ),
    unit = "mmol/L × min",
    stringsAsFactors = FALSE
  )
}

fit_meta <- function(dat) {
  yi <- dat$effect
  vi <- dat$se^2
  k <- length(yi)
  stopifnot(k == 2L, all(is.finite(yi)), all(is.finite(vi)), all(vi > 0))

  w_common <- 1 / vi
  common_estimate <- sum(w_common * yi) / sum(w_common)
  common_se <- sqrt(1 / sum(w_common))
  q <- sum(w_common * (yi - common_estimate)^2)
  q_df <- k - 1L
  i2 <- if (q > 0) max(0, (q - q_df) / q) * 100 else 0

  reml_objective <- function(tau2) {
    w <- 1 / (vi + tau2)
    mu <- sum(w * yi) / sum(w)
    0.5 * (sum(log(vi + tau2)) + log(sum(w)) + sum(w * (yi - mu)^2))
  }
  closed_form_tau2 <- max(0, ((yi[1] - yi[2])^2 - sum(vi)) / 2)
  upper <- max(1, closed_form_tau2 * 20, stats::var(yi) * 20, max(vi) * 20)
  tau2 <- optimize(reml_objective, interval = c(0, upper), tol = 1e-12)$minimum
  if (tau2 < 1e-10) tau2 <- 0

  w_random <- 1 / (vi + tau2)
  random_estimate <- sum(w_random * yi) / sum(w_random)
  random_se <- sqrt(1 / sum(w_random))

  q_random <- sum(w_random * (yi - random_estimate)^2)
  hk_scale <- max(1, q_random / q_df)
  hk_se <- sqrt(hk_scale / sum(w_random))
  hk_critical <- qt(0.975, df = q_df)

  results <- data.frame(
    model = c("REML random effects", "Common effect", "REML with HKSJ interval"),
    estimate = c(random_estimate, common_estimate, random_estimate),
    ci_low = c(
      random_estimate - qnorm(0.975) * random_se,
      common_estimate - qnorm(0.975) * common_se,
      random_estimate - hk_critical * hk_se
    ),
    ci_high = c(
      random_estimate + qnorm(0.975) * random_se,
      common_estimate + qnorm(0.975) * common_se,
      random_estimate + hk_critical * hk_se
    ),
    Q = q,
    Q_df = q_df,
    Q_p = pchisq(q, df = q_df, lower.tail = FALSE),
    tau2 = c(tau2, 0, tau2),
    I2_percent = i2,
    stringsAsFactors = FALSE
  )

  list(results = results, closed_form_tau2 = closed_form_tau2)
}

analysis_specs <- data.frame(
  analysis = c(
    "Primary analysis: Enomoto Day 2 + Nakamura",
    "Enomoto Day 1 sensitivity",
    "Enomoto Day 3 sensitivity",
    "Enomoto Day 4 sensitivity",
    "Nakamura P-threshold sensitivity: P = 0.0005",
    "Nakamura P-threshold sensitivity: P = 0.0001"
  ),
  enomoto_day = c(2L, 1L, 3L, 4L, 2L, 2L),
  nakamura_p = c(0.001, 0.001, 0.001, 0.001, 0.0005, 0.0001),
  stringsAsFactors = FALSE
)

all_models <- list()
validation <- list()
primary_data <- NULL
for (i in seq_len(nrow(analysis_specs))) {
  spec <- analysis_specs[i,]
  dat <- make_study_data(spec$enomoto_day, spec$nakamura_p)
  fit <- fit_meta(dat)
  selected <- if (i == 1L) fit$results else fit$results[fit$results$model == "REML random effects",]
  selected$analysis <- spec$analysis
  selected$k <- nrow(dat)
  selected$enomoto_day <- spec$enomoto_day
  selected$nakamura_p_assumed <- spec$nakamura_p
  all_models[[i]] <- selected
  validation[[i]] <- data.frame(
    analysis = spec$analysis,
    optimized_REML_tau2 = fit$results$tau2[fit$results$model == "REML random effects"],
    closed_form_REML_tau2 = fit$closed_form_tau2,
    absolute_difference = abs(fit$results$tau2[fit$results$model == "REML random effects"] - fit$closed_form_tau2),
    stringsAsFactors = FALSE
  )
  if (i == 1L) primary_data <- dat
}

models <- do.call(rbind, all_models)
models <- models[c("analysis", "model", "k", "estimate", "ci_low", "ci_high", "Q", "Q_df", "Q_p", "tau2", "I2_percent", "enomoto_day", "nakamura_p_assumed")]
validation <- do.call(rbind, validation)
validation$status <- ifelse(validation$absolute_difference < 0.05, "PASS", "FAIL")
stopifnot(all(validation$status == "PASS"))

numeric_columns <- c("estimate", "ci_low", "ci_high", "Q", "Q_p", "tau2", "I2_percent")
models[numeric_columns] <- lapply(models[numeric_columns], round, 4)
primary_data[c("effect", "se", "ci_low", "ci_high")] <- lapply(primary_data[c("effect", "se", "ci_low", "ci_high")], round, 4)
validation[c("optimized_REML_tau2", "closed_form_REML_tau2", "absolute_difference")] <-
  lapply(validation[c("optimized_REML_tau2", "closed_form_REML_tau2", "absolute_difference")], round, 4)

write.csv(primary_data, file.path(table_dir, "primary_study_effects.csv"), row.names = FALSE)
write.csv(models, file.path(table_dir, "primary_and_sensitivity_models.csv"), row.names = FALSE)
write.csv(validation, file.path(table_dir, "analysis_validation.csv"), row.names = FALSE)

primary_random <- models[
  models$analysis == "Primary analysis: Enomoto Day 2 + Nakamura" &
    models$model == "REML random effects",
]
stopifnot(nrow(primary_random) == 1L)

draw_forest <- function(device) {
  device()
  par(mar = c(5.2, 10.8, 0.8, 2.0), family = "Times New Roman")
  estimates <- c(primary_data$effect, primary_random$estimate)
  lower <- c(primary_data$ci_low, primary_random$ci_low)
  upper <- c(primary_data$ci_high, primary_random$ci_high)
  labels <- c(primary_data$study, "REML random effects")
  y <- c(3.2, 2.2, 0.8)
  xlim <- range(c(0, lower, upper), finite = TRUE)
  xlim <- xlim + c(-15, 15)
  plot(NA, xlim = xlim, ylim = c(0.2, 3.8), yaxt = "n", ylab = "",
       xlab = "Later minus earlier glucose iAUC (mmol/L × min)", bty = "n")
  abline(v = 0, col = "#6B7280", lty = 2, lwd = 1)
  segments(lower[1:2], y[1:2], upper[1:2], y[1:2], col = "#1B4965", lwd = 2)
  points(estimates[1:2], y[1:2], pch = 15, cex = 1.25, col = "#1B4965")
  segments(lower[3], y[3], upper[3], y[3], col = "#0F766E", lwd = 2.5)
  points(estimates[3], y[3], pch = 18, cex = 1.6, col = "#0F766E")
  axis(2, at = y, labels = labels, las = 1, tick = FALSE, cex.axis = 0.95, line = -2)
  dev.off()
}

png_path <- file.path(figure_dir, "primary_delayed_dinner_forest.png")
pdf_path <- file.path(figure_dir, "primary_delayed_dinner_forest.pdf")
draw_forest(function() png(png_path, width = 1800, height = 760, res = 190, type = "cairo", bg = "white"))
draw_forest(function() cairo_pdf(pdf_path, width = 9.5, height = 4.0, family = "Times New Roman"))


stopifnot(
  all(is.finite(models$estimate)),
  all(is.finite(models$ci_low)),
  all(is.finite(models$ci_high)),
  all(file.exists(c(png_path, pdf_path)))
)

source(file.path(root, "meta_analysis", "record_session.R"))
record_session(file.path(output_dir, "session_info.txt"))
cat("Primary delayed-dinner meta-analysis complete.\n")
cat("Outputs:", output_dir, "\n")
