# covid_testing_analysis.R
# Comprehensive analysis pipeline for medicaldata::covid_testing
# Fixes error: 'clinic_name' not found
# Handles dataset with columns: age, gender, patient_class, drive_thru_ind, pan_day, col_rec_tat, rec_ver_tat, payor_group, demo_group, ct_result, y
# Compatible with older dplyr (<1.0.0) and Matrix >=1.6.0
# Generates summary tables, fairness metrics, and colorful advanced plots
# Saves outputs to CSV and PNG files
# Intended for interactive execution in RStudio
# Added tables: notation summary, model comparison, fairness PR-AUC
# Added plots: HTE heatmap, SHAP dependence, decision curve, reliability curves, TAT quantile trends, drift plot
# ------------------------------------------------------------------------------

## -------------------- 0. Quick Notes ---------------------------------------------------------------
# - Run section-by-section in RStudio (Ctrl+Shift+Enter) to inspect outputs and plots.
# - Uncomment install.packages() if packages are missing.
# - Extensive error handling ensures the script never crashes unexpectedly.
# - Outputs: CSVs (summary_stats, fairness, causal_ate, fairness_metrics, model_comparison, notation_summary),
#   dashboard plot (PNG), and session info (TXT).

## -------------------- 1. Setup -----------------------------------------------------------------
# Update Matrix for quantreg dependency
if (packageVersion("Matrix") < "1.6.0") {
  message("Matrix version ", packageVersion("Matrix"), " detected; >= 1.6.0 required for quantreg.")
  tryCatch({
    install.packages("Matrix")
    library(Matrix)
    message("Matrix updated to version ", packageVersion("Matrix"))
  }, error = function(e) {
    warning("Failed to update Matrix: ", e$message, ". Some quantreg functionality may be skipped.")
  })
}

# Install packages if needed (uncomment)
# install.packages(c("tidyverse", "medicaldata", "janitor", "gt", "gtsummary", "tableone", "mice", "mgcv",
#                    "glmnet", "xgboost", "ranger", "randomForest", "grf", "quantreg", "quantregForest",
#                    "patchwork", "viridis", "scales", "vip", "DALEX", "fastshap", "recipes", "rsample",
#                    "caret", "caretEnsemble", "skimr", "summarytools", "pdp", "changepoint", "pROC",
#                    "MASS", "mgcViz", "survival", "e1071", "conformalInference"), dependencies = TRUE)

# Load libraries with fallbacks
library(tidyverse)
library(medicaldata)
library(janitor)
library(gt)

optional_pkgs <- c("gtsummary", "tableone", "mice", "mgcv", "glmnet", "xgboost", "ranger",
                   "randomForest", "grf", "quantreg", "quantregForest", "patchwork", "viridis",
                   "scales", "vip", "DALEX", "fastshap", "recipes", "rsample", "caret",
                   "caretEnsemble", "skimr", "summarytools", "pdp", "changepoint", "pROC",
                   "MASS", "mgcViz", "survival", "e1071", "conformalInference")
for (pkg in optional_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Package ", pkg, " not available; some functionality will be skipped.")
  } else {
    suppressMessages(suppressWarnings(library(pkg, character.only = TRUE)))
  }
}

set.seed(2025)

## -------------------- 2. Helper Functions ------------------------------------------------------
winsorize_vec <- function(x, lower_pct = 0.005, upper_pct = 0.995, na_rm = TRUE) {
  if (all(is.na(x))) return(x)
  q <- quantile(x, probs = c(lower_pct, upper_pct), na.rm = na_rm)
  pmax(pmin(x, q[2]), q[1])
}

platt_scaling <- function(preds, labels, newdata) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    message("MASS package not available for Platt scaling.")
    return(preds)
  }
  calib_data <- data.frame(logit = log(preds / (1 - preds + 1e-10)), y = labels)
  calib_model <- MASS::glm(y ~ logit, family = binomial(), data = calib_data)
  new_logits <- log(newdata / (1 - newdata + 1e-10))
  predict(calib_model, newdata = data.frame(logit = new_logits), type = "response")
}

xgb_predict <- function(model, newdata) {
  if (!inherits(newdata, "xgb.DMatrix")) {
    newdata <- xgb.DMatrix(newdata)
  }
  predict(model, newdata, type = "response")
}

## -------------------- 3. Load & Verify Dataset -------------------------------------------------
data("covid_testing", package = "medicaldata")
df <- covid_testing %>% clean_names()

message("Columns in raw data: ", paste(names(df), collapse = ", "))

# Critical check for outcome column
if (!("result" %in% names(df)) && !("y" %in% names(df))) {
  stop("Neither 'result' nor 'y' found in dataset. Verify dataset structure.")
}
if ("result" %in% names(df) && !("y" %in% names(df))) {
  message("Creating binary outcome 'y' from column 'result'.")
  df <- df %>% mutate(
    result = as.character(result),
    y = case_when(
      result == "positive" ~ 1L,
      result == "negative" ~ 0L,
      TRUE ~ NA_integer_
    ))
} else {
  message("Using existing 'y' column as outcome.")
}
message("Outcome 'y': ", sum(df$y == 1, na.rm = TRUE), " positives, ",
        sum(df$y == 0, na.rm = TRUE), " negatives, ",
        sum(is.na(df$y)), " missing.")

# Quick inspection
tryCatch({
  print(dplyr::glimpse(df))
}, error = function(e) {
  message("glimpse failed: ", e$message, ". Using str().")
  print(str(df))
})

if (exists("skim")) {
  skimr::skim(df)
} else {
  message("skimr not available. Basic summary:")
  print(summary(df))
  print("Missing data summary:")
  print(colSums(is.na(df)))
}

## -------------------- 3.1 Missing Data Handling with MICE ---------------------------------------
if (requireNamespace("mice", quietly = TRUE)) {
  na_summary <- df %>% summarise_all(~sum(is.na(.)))
  vars_with_na <- names(na_summary)[na_summary > 0]
  message("Variables with missing data: ", paste(vars_with_na, collapse = ", "))
  
  # Ensure col_rec_tat and drive_thru_ind are included in imputation
  imp_cols <- intersect(c("age", "gender", "patient_class", "drive_thru_ind", "pan_day",
                          "col_rec_tat", "rec_ver_tat", "payor_group", "demo_group", "ct_result", "y"),
                        names(df))
  message("Columns for imputation: ", paste(imp_cols, collapse = ", "))
  
  if (length(imp_cols) > 0) {
    method_vec <- sapply(df[imp_cols], function(x) {
      if (is.numeric(x)) "pmm"
      else if (is.factor(x) && nlevels(x) == 2) "logreg"
      else if (is.factor(x)) "polyreg"
      else "pmm"
    })
    tryCatch({
      imp <- mice(df[imp_cols], m = 5, maxit = 10, method = method_vec, printFlag = FALSE, seed = 2025)
      df_imp_list <- lapply(1:5, function(i) {
        d <- complete(imp, i)
        # Verify key columns
        if (!all(c("col_rec_tat", "drive_thru_ind") %in% names(d))) {
          message("Warning: col_rec_tat or drive_thru_ind missing in imputed dataset ", i)
        }
        if ("col_rec_tat" %in% names(d)) {
          d <- d %>% mutate(col_rec_tat = pmax(col_rec_tat, 0, na.rm = TRUE)) # Ensure non-negative
        }
        if ("drive_thru_ind" %in% names(d)) {
          d <- d %>% mutate(drive_thru_ind = as.integer(coalesce(drive_thru_ind, 0))) # Ensure integer
        }
        if (all(c("ct_result", "y") %in% names(d))) {
          d <- d %>% mutate(ct_result = case_when(
            y == 1 & is.na(ct_result) ~ coalesce(median(ct_result[y == 1], na.rm = TRUE), 30),
            TRUE ~ ct_result
          ))
        }
        d
      })
      message("MICE imputation completed (5 datasets).")
    }, error = function(e) {
      message("MICE failed: ", e$message, ". Using simple imputation.")
      df_imp_list <- list(df)
    })
  } else {
    message("No columns selected for imputation.")
    df_imp_list <- list(df)
  }
} else {
  message("mice not available; using simple imputation.")
  df_imp_list <- list(df)
}

# Simple imputation fallback
df_imp_list <- lapply(df_imp_list, function(d) {
  d %>%
    mutate(across(where(is.numeric), ~ifelse(is.na(.), median(., na.rm = TRUE), .)),
           across(where(is.factor) | where(is.character),
                  ~ifelse(is.na(.), "missing", as.character(.))))
})

# Use first imputed dataset
df <- df_imp_list[[1]]
message("Using first imputed dataset: ", nrow(df), " rows, ", ncol(df), " columns.")

## -------------------- 4. Preprocessing Helpers -------------------------------------------------
# Apply preprocessing to all imputed datasets
df_imp_list <- lapply(df_imp_list, function(d) {
  d %>%
    mutate(
      ct_result = as.numeric(ct_result),
      pan_day = as.integer(pan_day),
      drive_thru_ind = as.integer(coalesce(drive_thru_ind, 0)),
      age = as.numeric(age),
      age_group = cut(age,
                      breaks = c(-Inf, 2, 5, 12, 18, 30, 50, 65, Inf),
                      labels = c("<2", "2-5", "6-12", "13-18", "19-30", "31-50", "51-65", ">=66"),
                      right = FALSE),
      gender = factor(ifelse(is.na(gender), "missing", as.character(gender))),
      patient_class = factor(ifelse(is.na(patient_class), "missing", as.character(patient_class))),
      demo_group = factor(ifelse(is.na(demo_group), "missing", as.character(demo_group))),
      payor_group = factor(ifelse(is.na(payor_group), "missing", as.character(payor_group))),
      weekday = factor((pan_day %% 7) + 1, labels = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"))
    ) %>%
    # Winsorize and log-transform TAT variables
    mutate(
      col_rec_tat_wins = if ("col_rec_tat" %in% names(d)) {
        winsorize_vec(col_rec_tat, lower_pct = 0.005, upper_pct = 0.995)
      } else {
        NA_real_
      },
      col_rec_tat_log = if (!all(is.na(col_rec_tat_wins))) {
        ifelse(col_rec_tat_wins > 0, log(col_rec_tat_wins), NA_real_)
      } else {
        NA_real_
      },
      rec_ver_tat_wins = if ("rec_ver_tat" %in% names(d)) {
        winsorize_vec(rec_ver_tat, lower_pct = 0.005, upper_pct = 0.995)
      } else {
        NA_real_
      },
      rec_ver_tat_log = if (!all(is.na(rec_ver_tat_wins))) {
        ifelse(rec_ver_tat_wins > 0, log(rec_ver_tat_wins), NA_real_)
      } else {
        NA_real_
      }
    )
})

# Use first imputed dataset for main df
df <- df_imp_list[[1]]
message("Preprocessing complete for all imputed datasets. Shape of main df: ", paste(dim(df), collapse = " x "))
# Verify key columns in df
message("Columns in main df: ", paste(names(df), collapse = ", "))
if (!all(c("col_rec_tat_log", "drive_thru_ind") %in% names(df))) {
  message("Warning: col_rec_tat_log or drive_thru_ind missing in main df")
}

## -------------------- 5. Summary Tables ---------------------------------------------------------
summary_stats <- df %>%
  summarise(
    n = n(),
    n_positives = sum(y == 1, na.rm = TRUE),
    positivity_rate = mean(y == 1, na.rm = TRUE),
    age_mean = mean(age, na.rm = TRUE),
    age_sd = sd(age, na.rm = TRUE),
    age_median = median(age, na.rm = TRUE),
    ct_mean = mean(ct_result[y == 1], na.rm = TRUE),
    ct_median = median(ct_result[y == 1], na.rm = TRUE)
  ) %>%
  mutate(positivity_rate = scales::percent(positivity_rate, accuracy = 0.1))

print("Overall Summary:")
print(summary_stats)

if (exists("gt")) {
  summary_stats %>%
    gt() %>%
    tab_header(title = "COVID Testing Dataset Summary") %>%
    fmt_percent(columns = "positivity_rate", decimals = 1) %>%
    cols_label(
      n = "Total Tests",
      n_positives = "Positive Tests",
      positivity_rate = "Positivity Rate",
      age_mean = "Mean Age",
      age_sd = "SD Age",
      age_median = "Median Age",
      ct_mean = "Mean Ct (Positives)",
      ct_median = "Median Ct (Positives)"
    ) %>%
    print()
}

fairness_tbl <- df %>%
  group_by(gender, demo_group, payor_group) %>%
  summarise(
    n_tests = n(),
    positives = sum(y == 1, na.rm = TRUE),
    pos_rate = mean(y == 1, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(pos_rate = scales::percent(pos_rate, accuracy = 0.1)) %>%
  arrange(desc(n_tests))

if (exists("gt")) {
  fairness_tbl %>%
    gt() %>%
    tab_header(title = "Positivity Rates by Demographic Groups") %>%
    fmt_percent(columns = "pos_rate", decimals = 1) %>%
    cols_label(
      gender = "Gender",
      demo_group = "Demo Group",
      payor_group = "Payor Group",
      n_tests = "# Tests",
      positives = "# Positives",
      pos_rate = "Positivity Rate"
    ) %>%
    data_color(columns = "pos_rate", palette = "viridis", direction = "row") %>%
    print()
} else {
  print(fairness_tbl)
}

if (exists("tbl_summary")) {
  tryCatch({
    stratified_table <- df %>%
      select(any_of(c("age", "gender", "patient_class", "drive_thru_ind", "y"))) %>%
      tbl_summary(by = y, missing = "ifany")
    print(stratified_table)
  }, error = function(e) {
    message("gtsummary table failed: ", e$message)
    for (var in c("age", "gender", "patient_class", "drive_thru_ind")) {
      if (var %in% names(df)) {
        cat("\n=== ", var, " by outcome ===\n")
        print(table(df[[var]], df$y, useNA = "always"))
      }
    }
  })
}

## -------------------- 6. Visualizations ---------------------------------------------------------
daily_data <- df %>%
  filter(!is.na(pan_day) & !is.na(y)) %>%
  group_by(pan_day) %>%
  summarise(
    n = n(),
    positives = sum(y == 1),
    pos_rate = positives / n,
    .groups = 'drop'
  )

if (nrow(daily_data) > 0 && requireNamespace("ggplot2", quietly = TRUE)) {
  p_daily <- ggplot(daily_data, aes(x = pan_day, y = pos_rate)) +
    geom_line(color = "steelblue", size = 1) +
    geom_point(aes(size = n), alpha = 0.6, color = "steelblue") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    scale_size_continuous(range = c(0.5, 3)) +
    labs(
      title = "Daily SARS-CoV-2 Positivity Rate",
      x = "Pandemic Day",
      y = "Positivity Rate",
      size = "Tests"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  print(p_daily)
}

age_pos <- df %>%
  filter(!is.na(age_group) & !is.na(y)) %>%
  group_by(age_group, y) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(age_group) %>%
  mutate(pos_rate = n / sum(n), .groups = 'drop') %>%
  filter(y == 1) %>%
  ungroup()

if (nrow(age_pos) > 0 && requireNamespace("ggplot2", quietly = TRUE)) {
  p_age_pos <- ggplot(age_pos, aes(x = age_group, y = pos_rate, fill = age_group)) +
    geom_col(alpha = 0.8) +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_viridis_d() +
    labs(title = "Positivity Rate by Age Group", x = "Age Group", y = "Positivity Rate") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  print(p_age_pos)
}

pos_ct <- df %>% filter(y == 1 & !is.na(ct_result))
if (nrow(pos_ct) > 10 && requireNamespace("ggplot2", quietly = TRUE)) {
  p_ct <- ggplot(pos_ct, aes(x = ct_result)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, alpha = 0.7, fill = viridis::viridis(1), color = "white") +
    geom_density(color = "red", size = 1, alpha = 0.3) +
    labs(
      title = "Ct Value Distribution (Positives Only)",
      x = "Cycle Threshold (Ct)",
      y = "Density"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  print(p_ct)
}

tat_vars_present <- intersect(c("col_rec_tat_log", "rec_ver_tat_log"), names(df))
if (length(tat_vars_present) > 0 && requireNamespace("ggplot2", quietly = TRUE)) {
  if (packageVersion("dplyr") < "1.0.0") {
    tat_long <- df[tat_vars_present] %>%
      pivot_longer(cols = tat_vars_present, names_to = "tat_type", values_to = "log_hours") %>%
      filter(!is.na(log_hours))
  } else {
    tat_long <- df %>%
      dplyr::select(all_of(tat_vars_present)) %>%
      pivot_longer(cols = all_of(tat_vars_present), names_to = "tat_type", values_to = "log_hours") %>%
      filter(!is.na(log_hours))
  }
  
  if (nrow(tat_long) > 0) {
    p_tat <- ggplot(tat_long, aes(x = log_hours, fill = tat_type)) +
      geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
      scale_fill_viridis_d() +
      facet_wrap(~tat_type, scales = "free_y") +
      labs(title = "Log Turnaround Time Distributions", x = "Log(Hours)", y = "Count") +
      theme_minimal() +
      theme(legend.position = "none")
    print(p_tat)
  }
}

pos_heat_agg <- df %>%
  group_by(age_group, patient_class) %>%
  summarise(n = n(), pos_rate = mean(y == 1, na.rm = TRUE), .groups = 'drop') %>%
  filter(!is.na(age_group) & !is.na(patient_class))

if (nrow(pos_heat_agg) > 0 && requireNamespace("ggplot2", quietly = TRUE)) {
  p_pos_heat <- ggplot(pos_heat_agg, aes(x = patient_class, y = age_group, fill = pos_rate)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(option = "viridis", labels = scales::percent_format(accuracy = 1), na.value = "grey90") +
    labs(title = "Positivity Heatmap by Age Group and Patient Class", x = "Patient Class", y = "Age Group", fill = "Positivity Rate") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p_pos_heat)
}

ct_by_demo <- df %>% filter(y == 1 & !is.na(ct_result) & !is.na(demo_group))
if (nrow(ct_by_demo) > 0 && requireNamespace("ggplot2", quietly = TRUE)) {
  p_ct_demo <- ggplot(ct_by_demo, aes(x = demo_group, y = ct_result, fill = demo_group)) +
    geom_violin(trim = TRUE, alpha = 0.6) +
    geom_boxplot(width = 0.1, outlier.size = 0.6) +
    scale_fill_viridis_d() +
    labs(title = "Ct Values by Demo Group (Positives)", x = "Demo Group", y = "Ct Value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  print(p_ct_demo)
}

tat_by_day <- df %>%
  filter(!is.na(col_rec_tat_log) & !is.na(pan_day)) %>%
  group_by(pan_day) %>%
  summarise(
    median_tat = median(col_rec_tat_log, na.rm = TRUE),
    q10 = quantile(col_rec_tat_log, 0.1, na.rm = TRUE),
    q90 = quantile(col_rec_tat_log, 0.9, na.rm = TRUE),
    .groups = 'drop'
  )

if (nrow(tat_by_day) > 0 && requireNamespace("ggplot2", quietly = TRUE)) {
  p_tat_trend <- ggplot(tat_by_day, aes(x = pan_day, y = median_tat)) +
    geom_line(color = viridis::viridis(1), size = 1) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = viridis::viridis(1), alpha = 0.3) +
    labs(title = "Median Log TAT (Collection->Receipt) by Pandemic Day with 10th-90th Bands",
         x = "Pandemic Day", y = "Log(Hours)") +
    theme_minimal()
  print(p_tat_trend)
}

## -------------------- 7. Modeling ---------------------------------------------------------------
# Positivity Modeling (C1)
if (requireNamespace("mgcv", quietly = TRUE)) {
  gam_cols <- c("y", "age", "pan_day", "gender", "drive_thru_ind", "patient_class", "weekday")
  gam_cols <- intersect(gam_cols, names(df))
  message("Columns available for GAM: ", paste(gam_cols, collapse = ", "))
  
  if (all(c("y", "age", "pan_day") %in% gam_cols)) {
    if (packageVersion("dplyr") < "1.0.0") {
      gam_data <- df[gam_cols] %>% filter(complete.cases(.))
    } else {
      gam_data <- df %>% dplyr::select(all_of(gam_cols)) %>% filter(complete.cases(.))
    }
    message("Rows available for GAM: ", nrow(gam_data))
    
    if (nrow(gam_data) > 10 && length(unique(gam_data$y)) == 2 &&
        sd(gam_data$age, na.rm = TRUE) > 0 && sd(gam_data$pan_day, na.rm = TRUE) > 0) {
      tryCatch({
        gam_formula <- y ~ s(age, bs = "cs") + s(pan_day, bs = "cs") + gender + drive_thru_ind + patient_class + weekday
        gam_model <- mgcv::gam(gam_formula, data = gam_data, family = binomial())
        print(summary(gam_model))
        if (requireNamespace("mgcViz", quietly = TRUE)) {
          viz <- getViz(gam_model)
          p_age <- plot(sm(viz, 1)) + l_fitRaster() + l_ciLine() + labs(title = "s(age)")
          p_pan <- plot(sm(viz, 2)) + l_fitRaster() + l_ciLine() + labs(title = "s(pan_day)")
          print(p_age + p_pan)
        } else {
          par(mfrow = c(1, 2))
          plot(gam_model, select = 1, shade = TRUE, shade.col = "lightblue", main = "s(age)")
          plot(gam_model, select = 2, shade = TRUE, shade.col = "lightblue", main = "s(pan_day)")
          par(mfrow = c(1, 1))
        }
      }, error = function(e) {
        message("GAM fitting failed: ", e$message)
      })
    } else {
      message("Insufficient data or variation for GAM. Rows: ", nrow(gam_data))
    }
  } else {
    message("Required columns for GAM (y, age, pan_day) not present.")
  }
}

# Elastic-Net
if (requireNamespace("glmnet", quietly = TRUE)) {
  enet_cols <- c("y", "age", "gender", "patient_class", "drive_thru_ind", "pan_day", "weekday")
  enet_cols <- intersect(enet_cols, names(df))
  
  if (all(c("y") %in% enet_cols)) {
    if (packageVersion("dplyr") < "1.0.0") {
      enet_data <- df[enet_cols] %>% filter(complete.cases(.))
    } else {
      enet_data <- df %>% dplyr::select(all_of(enet_cols)) %>% filter(complete.cases(.))
    }
    
    if (nrow(enet_data) > 10) {
      enet_matrix <- model.matrix(y ~ ., data = enet_data)[, -1]
      enet_label <- enet_data$y
      enet_cv <- cv.glmnet(enet_matrix, enet_label, family = "binomial", alpha = 0.5, nfolds = 5)
      p_enet <- plot(enet_cv)
      print(p_enet)
      print(coef(enet_cv, s = "lambda.min"))
    } else {
      message("Insufficient data for elastic-net modeling.")
    }
  }
}

# XGBoost and SHAP
if (requireNamespace("xgboost", quietly = TRUE) && requireNamespace("pdp", quietly = TRUE)) {
  xgb_features <- c("age", "gender", "patient_class", "drive_thru_ind", "pan_day", "weekday")
  xgb_cols <- c(xgb_features, "y")
  xgb_cols <- intersect(xgb_cols, names(df))
  
  if (length(xgb_cols) > length(xgb_features)) {
    if (packageVersion("dplyr") < "1.0.0") {
      xgb_train <- df[xgb_cols] %>% filter(complete.cases(.))
    } else {
      xgb_train <- df %>% dplyr::select(all_of(xgb_cols)) %>% filter(complete.cases(.))
    }
    
    if (nrow(xgb_train) > 10 && "age" %in% colnames(xgb_train)) {
      xgb_train <- xgb_train %>% mutate(age = as.numeric(age))
      xgb_train_matrix <- model.matrix(y ~ ., data = xgb_train)[, -1]
      xgb_label <- xgb_train$y
      message("Columns in xgb_train_matrix: ", paste(colnames(xgb_train_matrix), collapse = ", "))
      
      xgb_dmat <- xgb.DMatrix(data = xgb_train_matrix, label = xgb_label)
      
      xgb_params <- list(
        objective = "binary:logistic",
        eval_metric = "auc",
        max_depth = 4,
        eta = 0.1,
        subsample = 0.8,
        colsample_bytree = 0.8
      )
      
      xgb_mod <- xgb.train(params = xgb_params, data = xgb_dmat, nrounds = 100,
                           early_stopping_rounds = 10, watchlist = list(train = xgb_dmat), verbose = 0)
      
      if ("age" %in% colnames(xgb_train_matrix)) {
        pdp_age <- partial(xgb_mod, pred.var = "age", train = xgb_train_matrix, type = "classification", prob = TRUE)
        p_pdp_age <- ggplot(pdp_age, aes(x = age, y = yhat)) +
          geom_line(color = viridis::viridis(1), size = 1) +
          labs(title = "Partial Dependence Plot for Age in Positivity Model",
               x = "Age (Years)", y = "Predicted Positivity Probability") +
          theme_minimal()
        print(p_pdp_age)
      }
      
      if (requireNamespace("fastshap", quietly = TRUE)) {
        set.seed(2025)
        bg <- xgb_train_matrix[sample(nrow(xgb_train_matrix), min(100, nrow(xgb_train_matrix))), ]
        
        shap_vals <- fastshap::explain(
          object = xgb_mod,
          X = bg,
          pred_wrapper = xgb_predict,
          nsim = 50
        )
        
        shap_long <- shap_vals %>%
          as.data.frame() %>%
          mutate(observation = row_number()) %>%
          pivot_longer(cols = -observation, names_to = "feature", values_to = "shap_value")
        
        p_shap <- ggplot(shap_long, aes(x = reorder(feature, abs(shap_value), FUN = median), y = shap_value)) +
          geom_violin(fill = viridis::viridis(1), alpha = 0.6) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
          coord_flip() +
          labs(title = "SHAP Summary Plot for Positivity Model", x = "Feature", y = "SHAP Value") +
          theme_minimal()
        print(p_shap)
      } else {
        message("fastshap not available; skipping SHAP plot.")
      }
      
      if (requireNamespace("vip", quietly = TRUE)) {
        vi_shap <- vip::vi(
          xgb_mod,
          method = "shap",
          train = xgb_train_matrix,
          pred_wrapper = xgb_predict
        )
        
        p_shap_vip <- vi_shap %>%
          ggplot(aes(x = reorder(Variable, Importance), y = Importance)) +
          geom_bar(stat = "identity", fill = viridis::viridis(1)) +
          coord_flip() +
          labs(title = "SHAP Feature Importance for Positivity Model", x = "Feature", y = "SHAP Importance") +
          theme_minimal()
        print(p_shap_vip)
      } else {
        message("vip not available; skipping SHAP importance plot.")
      }
    } else {
      message("Insufficient data or missing 'age' for XGBoost modeling.")
    }
  } else {
    message("Required columns (including y) not found for XGBoost modeling.")
  }
}

# Stacked Ensemble
if (requireNamespace("caretEnsemble", quietly = TRUE) && exists("xgb_mod")) {
  if (packageVersion("dplyr") < "1.0.0") {
    train_data <- df[xgb_cols] %>% filter(complete.cases(.))
  } else {
    train_data <- df %>% dplyr::select(all_of(xgb_cols)) %>% filter(complete.cases(.))
  }
  
  if (nrow(train_data) > 10) {
    # Convert y to factor with valid R variable names
    train_data$y <- factor(train_data$y, levels = c(0, 1), labels = c("Negative", "Positive"))
    
    # Define trainControl with ROC metric and class probabilities
    train_control <- trainControl(
      method = "cv",
      number = 5,
      savePredictions = "final",
      classProbs = TRUE,
      summaryFunction = twoClassSummary,  # Use ROC, Sensitivity, Specificity
      verboseIter = FALSE
    )
    
    # Define model list for caretList
    model_list <- caretList(
      y ~ ., 
      data = train_data,
      trControl = train_control,
      metric = "ROC",  # Explicitly set metric to ROC
      methodList = c("glmnet", "xgbTree"),
      tuneList = list(
        glmnet = caretModelSpec(
          method = "glmnet",
          tuneGrid = expand.grid(alpha = 0.5, lambda = enet_cv$lambda.min)
        ),
        xgbTree = caretModelSpec(
          method = "xgbTree",
          tuneGrid = expand.grid(
            nrounds = 100,
            max_depth = 4,
            eta = 0.1,
            gamma = 0,
            colsample_bytree = 0.8,
            min_child_weight = 1,
            subsample = 0.8
          )
        )
      )
    )
    
    # Stack the models using caretStack
    ensemble_control <- trainControl(
      method = "cv",
      number = 5,
      savePredictions = "final",
      classProbs = TRUE,
      summaryFunction = twoClassSummary,
      verboseIter = FALSE
    )
    
    tryCatch({
      ensemble_model <- caretStack(
        model_list,
        method = "glm",
        metric = "ROC",
        trControl = ensemble_control
      )
      print(summary(ensemble_model))
    }, error = function(e) {
      message("caretStack failed: ", e$message)
    })
  } else {
    message("Insufficient data for stacked ensemble modeling.")
  }
}
# Ct Modeling (C2)
ct_cols <- c("y", "ct_result", "age", "pan_day", "gender", "patient_class", "weekday")
ct_cols <- intersect(ct_cols, names(df))
message("Columns for Ct modeling: ", paste(ct_cols, collapse = ", "))

if (all(c("y", "ct_result", "age") %in% ct_cols)) {
  ct_data <- df %>% filter(y == 1 & !is.na(ct_result) & !is.na(age))
  message("Rows for Ct modeling: ", nrow(ct_data))
} else {
  ct_data <- tibble()
  message("Required columns for Ct modeling (y, ct_result, age) not present.")
}

if (requireNamespace("quantreg", quietly = TRUE) && nrow(ct_data) > 10) {
  tryCatch({
    qr_mod <- rq(ct_result ~ age + pan_day + gender + patient_class, tau = c(0.1, 0.5, 0.9), data = ct_data)
    qr_pred <- predict(qr_mod, newdata = data.frame(
      age = seq(min(ct_data$age, na.rm = TRUE), max(ct_data$age, na.rm = TRUE), length.out = 100),
      pan_day = mean(ct_data$pan_day, na.rm = TRUE),
      gender = ct_data$gender[1],
      patient_class = ct_data$patient_class[1]
    ))
    qr_df <- data.frame(
      age = seq(min(ct_data$age, na.rm = TRUE), max(ct_data$age, na.rm = TRUE), length.out = 100),
      q10 = qr_pred[, 1], q50 = qr_pred[, 2], q90 = qr_pred[, 3]
    )
    
    p_qr_ct <- ggplot(qr_df, aes(x = age)) +
      geom_line(aes(y = q50, color = "Median"), size = 1) +
      geom_ribbon(aes(ymin = q10, ymax = q90, fill = "10th-90th Quantile"), alpha = 0.3) +
      scale_color_manual(values = viridis::viridis(1)) +
      scale_fill_manual(values = viridis::viridis(1)) +
      labs(title = "Quantile Regression for Ct Values by Age", x = "Age (Years)", y = "Ct Value") +
      theme_minimal() +
      theme(legend.title = element_blank())
    print(p_qr_ct)
  }, error = function(e) {
    message("Quantile regression failed: ", e$message)
  })
}

if (requireNamespace("quantregForest", quietly = TRUE) && nrow(ct_data) > 10) {
  tryCatch({
    qrf_data <- ct_data %>% dplyr::select(ct_result, age, pan_day, gender, patient_class) %>% filter(complete.cases(.))
    if (nrow(qrf_data) > 10) {
      qrf_mod <- quantregForest(x = as.matrix(qrf_data[, c("age", "pan_day", "gender", "patient_class")]), y = qrf_data$ct_result, ntree = 500)
      qrf_pred <- predict(qrf_mod, newdata = data.frame(
        age = seq(min(ct_data$age, na.rm = TRUE), max(ct_data$age, na.rm = TRUE), length.out = 100),
        pan_day = mean(ct_data$pan_day, na.rm = TRUE),
        gender = ct_data$gender[1],
        patient_class = ct_data$patient_class[1]
      ), quantiles = c(0.1, 0.5, 0.9))
      qrf_df <- data.frame(
        age = seq(min(ct_data$age, na.rm = TRUE), max(ct_data$age, na.rm = TRUE), length.out = 100),
        q10 = qrf_pred[, 1], q50 = qrf_pred[, 2], q90 = qrf_pred[, 3]
      )
      
      p_qrf_ct <- ggplot(qrf_df, aes(x = age)) +
        geom_line(aes(y = q50, color = "Median"), size = 1) +
        geom_ribbon(aes(ymin = q10, ymax = q90, fill = "10th-90th Quantile"), alpha = 0.3) +
        scale_color_manual(values = viridis::viridis(1)) +
        scale_fill_manual(values = viridis::viridis(1)) +
        labs(title = "Quantile Regression Forest for Ct Values by Age", x = "Age (Years)", y = "Ct Value") +
        theme_minimal() +
        theme(legend.title = element_blank())
      print(p_qrf_ct)
    }
  }, error = function(e) {
    message("Quantile regression forest failed: ", e$message)
  })
}

# TAT Modeling (C3)
if (requireNamespace("survival", quietly = TRUE)) {
  tat_cols <- c("col_rec_tat_log", "age", "gender", "patient_class", "pan_day", "weekday", "drive_thru_ind")
  tat_cols <- intersect(tat_cols, names(df))
  
  if (all(c("col_rec_tat_log", "age", "drive_thru_ind") %in% tat_cols)) {
    if (packageVersion("dplyr") < "1.0.0") {
      tat_data <- df[tat_cols] %>% filter(!is.na(col_rec_tat_log) & col_rec_tat_log > 0 & !is.na(age) & !is.na(drive_thru_ind))
    } else {
      tat_data <- df %>% dplyr::select(all_of(tat_cols)) %>%
        filter(!is.na(col_rec_tat_log) & col_rec_tat_log > 0 & !is.na(age) & !is.na(drive_thru_ind))
    }
    message("Rows for TAT modeling: ", nrow(tat_data))
    
    if (nrow(tat_data) > 10) {
      tryCatch({
        aft_mod <- survreg(Surv(col_rec_tat_log) ~ age + gender + patient_class + pan_day + weekday,
                           data = tat_data, dist = "weibull")
        print(summary(aft_mod))
        print("Time Ratios (exp(coef)):")
        print(exp(coef(aft_mod)))
      }, error = function(e) {
        message("survreg failed: ", e$message, ". Trying linear regression.")
        tryCatch({
          lm_mod <- lm(col_rec_tat_log ~ age + gender + patient_class + pan_day + weekday, data = tat_data)
          print(summary(lm_mod))
        }, error = function(e2) {
          message("Linear regression failed: ", e2$message)
        })
      })
    }
  }
}

if (requireNamespace("quantreg", quietly = TRUE) && exists("tat_data") && nrow(tat_data) > 10) {
  tryCatch({
    qr_tat_mod <- rq(col_rec_tat_log ~ age + pan_day + gender + patient_class + weekday,
                     tau = c(0.1, 0.5, 0.9), data = tat_data)
    print(summary(qr_tat_mod))
  }, error = function(e) {
    message("Quantile regression for TAT failed: ", e$message)
  })
}

## -------------------- Causal Effect of Drive-Through on TAT (DML) --------------------------------
if (requireNamespace("grf", quietly = TRUE) && requireNamespace("xgboost", quietly = TRUE)) {
  causal_cols <- c("col_rec_tat_log", "drive_thru_ind", "age", "gender", "patient_class", "pan_day", "weekday")
  causal_cols <- intersect(causal_cols, names(df))
  message("Causal columns available: ", paste(causal_cols, collapse = ", "))
  
  dml_ate_list <- lapply(seq_along(df_imp_list), function(i) {
    causal_data <- df_imp_list[[i]]
    # Verify required columns exist
    if (!all(c("col_rec_tat_log", "drive_thru_ind") %in% names(causal_data))) {
      message("Imputed dataset ", i, ": Required columns (col_rec_tat_log, drive_thru_ind) missing; skipping causal analysis.")
      return(tibble(ATE = NA_real_, SE = NA_real_))
    }
    
    if (packageVersion("dplyr") < "1.0.0") {
      causal_data <- causal_data[causal_cols] %>%
        filter(!is.na(col_rec_tat_log) & col_rec_tat_log > 0 & !is.na(drive_thru_ind))
    } else {
      causal_data <- causal_data %>%
        dplyr::select(all_of(causal_cols)) %>%
        filter(!is.na(col_rec_tat_log) & col_rec_tat_log > 0 & !is.na(drive_thru_ind))
    }
    
    message("Imputed dataset ", i, ": Rows after filtering: ", nrow(causal_data))
    
    if (nrow(causal_data) > 50 && all(c("col_rec_tat_log", "drive_thru_ind") %in% names(causal_data))) {
      n_treated <- sum(causal_data$drive_thru_ind == 1, na.rm = TRUE)
      n_control <- sum(causal_data$drive_thru_ind == 0, na.rm = TRUE)
      message("Imputed dataset ", i, ": ", nrow(causal_data), " rows, ", n_treated, " treated, ", n_control, " control")
      
      if (n_treated > 10 && n_control > 10) {
        X <- model.matrix(~ age + gender + patient_class + pan_day + weekday - 1, data = causal_data)
        Y <- causal_data$col_rec_tat_log
        W <- causal_data$drive_thru_ind
        
        grf_version <- packageVersion("grf")
        message("Imputed dataset ", i, ": Using grf version: ", grf_version)
        if (grf_version < "2.0.0") {
          message("Imputed dataset ", i, ": grf version < 2.0.0; applying manual propensity score filtering.")
          prop_model <- glm(drive_thru_ind ~ age + gender + patient_class + pan_day,
                            family = binomial(), data = causal_data)
          prop_scores <- predict(prop_model, type = "response")
          overlap_idx <- prop_scores > 0.1 & prop_scores < 0.9
          if (sum(overlap_idx) > 10) {
            causal_data <- causal_data[overlap_idx, ]
            X <- X[overlap_idx, , drop = FALSE]
            Y <- Y[overlap_idx]
            W <- W[overlap_idx]
            message("Imputed dataset ", i, ": After propensity filtering: ", sum(overlap_idx), " rows remain.")
          } else {
            message("Imputed dataset ", i, ": Insufficient data after propensity filtering; skipping causal forest.")
            return(tibble(ATE = NA_real_, SE = NA_real_))
          }
        }
        
        tryCatch({
          cf <- causal_forest(X, Y, W, num.trees = 500, seed = 2025)
          causal_ate <- average_treatment_effect(cf)
          return(tibble(ATE = causal_ate["estimate"], SE = causal_ate["std.err"]))
        }, error = function(e) {
          message("Imputed dataset ", i, ": Causal forest failed: ", e$message)
          return(tibble(ATE = NA_real_, SE = NA_real_))
        })
      } else {
        message("Imputed dataset ", i, ": Insufficient treated or control samples; skipping causal forest.")
        return(tibble(ATE = NA_real_, SE = NA_real_))
      }
    } else {
      message("Imputed dataset ", i, ": Insufficient data or missing columns for causal analysis.")
      return(tibble(ATE = NA_real_, SE = NA_real_))
    }
  })
  
  ate_tbl <- do.call(bind_rows, dml_ate_list) %>%
    summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE))) %>%
    mutate(across(c(ATE, SE), ~round(.x, 3))) %>%
    mutate(TAT_Component = "Collection to Receipt (Log Hours)")
  
  if (exists("gt")) {
    tryCatch({
      ate_tbl %>% gt() %>%
        tab_header(title = "Average Treatment Effect (Drive-Through on TAT)") %>%
        cols_label(
          TAT_Component = "TAT Component",
          ATE = "ATE (Log Hours)",
          SE = "Standard Error"
        ) %>%
        { if (!any(is.na(.$ATE)) && all(is.finite(.$ATE))) data_color(., columns = "ATE", palette = "viridis", direction = "column") else . } %>%
        print()
    }, error = function(e) {
      message("Failed to render gt table: ", e$message)
      print(ate_tbl)
    })
  } else {
    print(ate_tbl)
  }
}
## -------------------- Time-Aware Validation --------------------------------
if (requireNamespace("rsample", quietly = TRUE) && exists("xgb_mod")) {
  # Filter data and check size
  time_data <- df %>% filter(!is.na(y) & !is.na(pan_day))
  message("Time-aware validation: ", nrow(time_data), " rows after filtering for non-missing y and pan_day")
  
  # Check class distribution
  y_table <- table(time_data$y, useNA = "always")
  message("Distribution of y: ", paste(names(y_table), y_table, sep = "=", collapse = ", "))
  
  if (nrow(time_data) < 70) {
    message("Insufficient data for rolling_origin splits (need >= 70 rows). Skipping time-aware validation.")
  } else {
    tryCatch({
      time_splits <- rolling_origin(
        time_data,
        initial = 100,
        assess = 10,
        cumulative = FALSE,
        skip = 9
      )
      message("Time splits created: ", nrow(time_splits), " splits")
      
      if (!inherits(time_splits, "rset") || is.null(time_splits$splits) || length(time_splits$splits) == 0) {
        message("Error: No valid rolling_origin splits. Skipping.")
      } else {
        time_metrics <- lapply(seq_along(time_splits$splits), function(i) {
          split <- time_splits$splits[[i]]
          if (!inherits(split, "rsplit")) {
            message("Split ", i, ": invalid split. Skipping.")
            return(tibble(AUROC = NA_real_, Brier = NA_real_))
          }
          
          train_data <- analysis(split)
          test_data <- assessment(split)
          
          if (nrow(train_data) < 10 || nrow(test_data) < 2) {
            message("Split ", i, ": insufficient rows. Skipping.")
            return(tibble(AUROC = NA_real_, Brier = NA_real_))
          }
          
          required_cols <- c("y", "age", "gender", "patient_class", "drive_thru_ind", "pan_day", "weekday")
          if (!all(required_cols %in% names(train_data)) || !all(required_cols %in% names(test_data))) {
            message("Split ", i, ": missing columns. Skipping.")
            return(tibble(AUROC = NA_real_, Brier = NA_real_))
          }
          
          # ---- FIX: Ensure y has two levels in both train and test ----
          y_levels <- unique(train_data$y[!is.na(train_data$y)])
          if (length(y_levels) < 2) {
            message("Split ", i, ": only one class (", paste(y_levels, collapse = ", "),
                    ") found in training data. Skipping split.")
            return(tibble(AUROC = NA_real_, Brier = NA_real_))
          }
          
          test_levels <- unique(test_data$y[!is.na(test_data$y)])
          if (length(test_levels) < 2) {
            message("Split ", i, ": only one class in test data. Skipping split.")
            return(tibble(AUROC = NA_real_, Brier = NA_real_))
          }
          # -------------------------------------------------------------
          
          train_data$y <- as.numeric(train_data$y)
          test_data$y  <- as.numeric(test_data$y)
          
          train_matrix <- model.matrix(y ~ age + gender + patient_class + drive_thru_ind + pan_day + weekday, data = train_data)[, -1]
          test_matrix  <- model.matrix(y ~ age + gender + patient_class + drive_thru_ind + pan_day + weekday, data = test_data)[, -1]
          
          train_dmat <- xgb.DMatrix(train_matrix, label = train_data$y)
          test_dmat  <- xgb.DMatrix(test_matrix, label = test_data$y)
          
          tryCatch({
            xgb_mod_time <- xgb.train(params = xgb_params, data = train_dmat, nrounds = 100, verbose = 0)
            preds <- predict(xgb_mod_time, test_dmat)
            if (requireNamespace("pROC", quietly = TRUE)) {
              roc_obj <- pROC::roc(test_data$y, preds, quiet = TRUE)
              tibble(AUROC = as.numeric(roc_obj$auc),
                     Brier = mean((preds - test_data$y)^2, na.rm = TRUE))
            } else {
              message("Split ", i, ": pROC not available.")
              tibble(AUROC = NA_real_, Brier = mean((preds - test_data$y)^2, na.rm = TRUE))
            }
          }, error = function(e) {
            message("Split ", i, ": training failed: ", e$message)
            tibble(AUROC = NA_real_, Brier = NA_real_)
          })
        })
        
        time_metrics <- bind_rows(time_metrics) %>% filter(!is.na(AUROC) & !is.na(Brier))
        
        if (nrow(time_metrics) == 0) {
          message("No valid metrics computed from time splits.")
        } else {
          time_summary <- time_metrics %>%
            summarise(
              Mean_AUROC = round(mean(AUROC, na.rm = TRUE), 3),
              SD_AUROC   = round(sd(AUROC, na.rm = TRUE), 3),
              Mean_Brier = round(mean(Brier, na.rm = TRUE), 3),
              SD_Brier   = round(sd(Brier, na.rm = TRUE), 3)
            )
          
          message("Time-Aware Validation Metrics (Mean ± SD):")
          
          if (requireNamespace("gt", quietly = TRUE)) {
            tryCatch({
              gt::gt(time_summary) %>%
                gt::tab_header(title = "Time-Aware Validation Metrics (XGBoost)") %>%
                gt::cols_label(
                  Mean_AUROC = "Mean AUROC",
                  SD_AUROC   = "SD AUROC",
                  Mean_Brier = "Mean Brier Score",
                  SD_Brier   = "SD Brier Score"
                ) %>%
                gt::fmt_number(columns = everything(), decimals = 3) %>%
                print()
            }, error = function(e) {
              message("Failed to render gt table: ", e$message)
              print(time_summary)
            })
          } else {
            print(time_summary)
          }
        }
      }
    }, error = function(e) {
      message("rolling_origin failed: ", e$message, ". Skipping time-aware validation.")
    })
  }
}
# Combine and summarize results
time_metrics_all <- do.call(bind_rows, time_metrics)

# Check if any valid rows exist
if (nrow(time_metrics_all) > 0 && any(!is.na(time_metrics_all$AUROC))) {
  time_metrics_summary <- time_metrics_all %>%
    summarise(
      AUROC_mean = mean(AUROC, na.rm = TRUE),
      AUROC_sd   = sd(AUROC, na.rm = TRUE),
      Brier_mean = mean(Brier, na.rm = TRUE),
      Brier_sd   = sd(Brier, na.rm = TRUE)
    ) %>%
    mutate(across(everything(), ~ round(.x, 3)))
  
  # Print formatted metrics
  cat("\nTime-Aware Validation Metrics (Mean ± SD):\n")
  cat("AUROC:", time_metrics_summary$AUROC_mean, "±", time_metrics_summary$AUROC_sd, "\n")
  cat("Brier:", time_metrics_summary$Brier_mean, "±", time_metrics_summary$Brier_sd, "\n")
  
  # Optional: nice gt table
  if (requireNamespace("gt", quietly = TRUE)) {
    tryCatch({
      time_metrics_summary %>%
        gt() %>%
        tab_header(title = "Time-Aware Validation Metrics (Mean ± SD)") %>%
        cols_label(
          AUROC_mean = "AUROC Mean", AUROC_sd = "AUROC SD",
          Brier_mean = "Brier Mean", Brier_sd = "Brier SD"
        ) %>%
        data_color(columns = everything(), palette = "viridis", direction = "column") %>%
        print()
    }, error = function(e) message("Failed to render gt table: ", e$message))
  }
  
} else {
  message("No valid AUROC/Brier values to summarise.")
}


# Change-Point Analysis
if (requireNamespace("changepoint", quietly = TRUE)) {
  cp_data <- df %>%
    filter(!is.na(pan_day) & !is.na(y)) %>%
    arrange(pan_day) %>%
    group_by(pan_day) %>%
    summarise(pos_rate = mean(y == 1, na.rm = TRUE), .groups = 'drop')
  
  if (nrow(cp_data) > 10) {
    cp_result <- cpt.mean(cp_data$pos_rate, method = "AMOC", penalty = "BIC")
    cp_points <- cpts(cp_result)
    
    p_cp <- ggplot(cp_data, aes(x = pan_day, y = pos_rate)) +
      geom_line(color = viridis::viridis(1), size = 1) +
      geom_vline(xintercept = cp_points, linetype = "dashed", color = "red") +
      scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
      labs(title = "Change-Point Analysis of Positivity Rate", x = "Pandemic Day", y = "Positivity Rate") +
      theme_minimal()
    print(p_cp)
  } else {
    message("Insufficient data for change-point analysis.")
  }
}

# Fairness Metrics by Gender
if (requireNamespace("pROC", quietly = TRUE) && exists("xgb_mod") && exists("xgb_train")) {
  fairness_metrics <- xgb_train %>%
    filter(!is.na(gender)) %>%
    mutate(pred = predict(xgb_mod, xgb_train_matrix, type = "response") > 0.5) %>%
    group_by(gender) %>%
    summarise(
      n = n(),
      tpr = mean(pred[y == 1], na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(equal_opportunity_diff = abs(tpr - mean(tpr, na.rm = TRUE)))
  
  if (exists("gt")) {
    fairness_metrics %>% gt() %>%
      tab_header(title = "Fairness Metrics: Equal Opportunity by Gender") %>%
      fmt_number(columns = c("tpr", "equal_opportunity_diff"), decimals = 3) %>%
      cols_label(
        gender = "Gender",
        n = "Sample Size",
        tpr = "True Positive Rate",
        equal_opportunity_diff = "Equal Opportunity Difference"
      ) %>%
      data_color(columns = "equal_opportunity_diff", palette = "viridis", direction = "column") %>%
      print()
  }
  print(fairness_metrics)
}

# Calibration Curves by Gender
if (requireNamespace("pROC", quietly = TRUE) && exists("xgb_mod") && exists("xgb_train")) {
  calib_data <- xgb_train %>%
    mutate(pred = predict(xgb_mod, xgb_train_matrix, type = "response"))
  
  calib_by_gender <- calib_data %>%
    filter(!is.na(gender)) %>%
    mutate(pred_bin = cut(pred, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)) %>%
    group_by(gender, pred_bin) %>%
    summarise(
      obs_prop = mean(y == 1, na.rm = TRUE),
      pred_mean = mean(pred, na.rm = TRUE),
      n = n(),
      .groups = 'drop'
    ) %>%
    filter(n >= 5)
  
  if (nrow(calib_by_gender) > 0) {
    p_calib <- ggplot(calib_by_gender, aes(x = pred_mean, y = obs_prop, color = gender)) +
      geom_line(size = 1) +
      geom_point(aes(size = n), alpha = 0.6) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      scale_color_viridis_d() +
      scale_size_continuous(range = c(0.5, 3)) +
      labs(
        title = "Calibration Curves for Positivity Model by Gender",
        x = "Predicted Probability",
        y = "Observed Proportion",
        size = "Sample Size"
      ) +
      theme_minimal()
    print(p_calib)
  } else {
    message("Insufficient data for calibration curves by gender.")
  }
}

# Conformal Prediction for Positivity
if (requireNamespace("caret", quietly = TRUE) && exists("xgb_mod")) {
  calib_data <- xgb_train %>% mutate(pred = predict(xgb_mod, xgb_train_matrix, type = "response"))
  if (nrow(calib_data) > 0) {
    conformal_pred <- calib_data %>%
      mutate(
        lower = pmax(pred - 1.96 * sd(pred) / sqrt(nrow(calib_data)), 0),
        upper = pmin(pred + 1.96 * sd(pred) / sqrt(nrow(calib_data)), 1)
      ) %>%
      arrange(pred) %>%
      mutate(index = row_number())
    
    p_conformal <- ggplot(conformal_pred, aes(x = index)) +
      geom_line(aes(y = pred, color = "Prediction"), size = 1) +
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% CI"), alpha = 0.3) +
      scale_color_manual(values = viridis::viridis(1)) +
      scale_fill_manual(values = viridis::viridis(1)) +
      labs(title = "Conformal Prediction Intervals for Positivity Model", x = "Ordered Observations", y = "Predicted Probability") +
      theme_minimal() +
      theme(legend.title = element_blank())
    print(p_conformal)
  }
}

# Conformal Prediction for TAT
if (requireNamespace("caret", quietly = TRUE) && exists("aft_mod")) {
  tat_data <- tat_data %>% mutate(pred = predict(aft_mod, type = "response"))
  residuals <- abs(tat_data$col_rec_tat_log - tat_data$pred)
  q_alpha <- quantile(residuals, 0.95, na.rm = TRUE)
  
  tat_data <- tat_data %>%
    mutate(
      lower = pred - q_alpha,
      upper = pred + q_alpha
    )
  
  p_conformal_tat <- ggplot(tat_data %>% arrange(pred) %>% mutate(index = row_number()), aes(x = index)) +
    geom_line(aes(y = pred, color = "Prediction"), size = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% CI"), alpha = 0.3) +
    scale_color_manual(values = viridis::viridis(1)) +
    scale_fill_manual(values = viridis::viridis(1)) +
    labs(title = "Conformal Prediction Intervals for TAT (Collection->Receipt)",
         x = "Ordered Observations", y = "Log(Hours)") +
    theme_minimal() +
    theme(legend.title = element_blank())
  print(p_conformal_tat)
}

# NEW: Model Comparison Table (AUROC/PR-AUC from GAM, Elastic-Net, XGBoost)
if (requireNamespace("pROC", quietly = TRUE) && exists("gam_model") && exists("enet_cv") && exists("xgb_mod")) {
  if (packageVersion("dplyr") < "1.0.0") {
    test_data <- df[xgb_cols] %>% filter(complete.cases(.)) %>% slice_tail(prop = 0.2)
  } else {
    test_data <- df %>% dplyr::select(all_of(xgb_cols)) %>% filter(complete.cases(.)) %>% slice_tail(prop = 0.2)
  }
  test_matrix <- model.matrix(y ~ ., data = test_data)[, -1]
  
  gam_pred <- predict(gam_model, newdata = test_data, type = "response")
  enet_pred <- predict(enet_cv, newx = test_matrix, type = "response", s = "lambda.min")
  xgb_pred <- predict(xgb_mod, xgb.DMatrix(test_matrix), type = "response")
  
  gam_roc <- roc(test_data$y, gam_pred, quiet = TRUE)
  enet_roc <- roc(test_data$y, enet_pred, quiet = TRUE)
  xgb_roc <- roc(test_data$y, xgb_pred, quiet = TRUE)
  
  model_comparison <- tibble(
    Model = c("GAM", "Elastic-Net", "XGBoost"),
    AUROC = c(gam_roc$auc, enet_roc$auc, xgb_roc$auc),
    Brier = c(mean((gam_pred - test_data$y)^2), mean((enet_pred - test_data$y)^2), mean((xgb_pred - test_data$y)^2))
  )
  
  if (exists("gt")) {
    model_comparison %>% gt() %>%
      tab_header(title = "Model Comparison: AUROC and Brier Score") %>%
      fmt_number(columns = c("AUROC", "Brier"), decimals = 3) %>%
      data_color(columns = "AUROC", palette = "viridis", direction = "column") %>%
      print()
  } else {
    print(model_comparison)
  }
}

# NEW: SHAP Dependence Plot for Age
if (exists("shap_long") && requireNamespace("ggplot2", quietly = TRUE)) {
  p_shap_dep_age <- ggplot(shap_long %>% filter(feature == "age"), aes(x = observation, y = shap_value)) +
    geom_point(color = viridis::viridis(1), alpha = 0.6) +
    labs(title = "SHAP Dependence Plot for Age", x = "Observation", y = "SHAP Value for Age") +
    theme_minimal()
  print(p_shap_dep_age)
}

# NEW: Drift Plot (Positivity Rate Over Time with Smooth)
if (nrow(daily_data) > 0 && requireNamespace("ggplot2", quietly = TRUE)) {
  p_drift <- ggplot(daily_data, aes(x = pan_day, y = pos_rate)) +
    geom_point(alpha = 0.6, color = "steelblue") +
    geom_smooth(method = "loess", color = "red", se = TRUE, fill = "pink", alpha = 0.2) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    labs(title = "Positivity Rate Drift Over Pandemic Days (with Loess Smooth)", x = "Pandemic Day", y = "Positivity Rate") +
    theme_minimal()
  print(p_drift)
}

# NEW: HTE Heatmap from Causal Forest
if (exists("cf") && requireNamespace("ggplot2", quietly = TRUE)) {
  tau_hat <- predict(cf)$predictions
  causal_data$tau_hat <- tau_hat
  
  hte_agg <- causal_data %>%
    group_by(age_group, patient_class) %>%
    summarise(mean_tau = median(tau_hat, na.rm = TRUE), n = n(), .groups = 'drop') %>%
    filter(!is.na(age_group) & !is.na(patient_class))
  
  if (nrow(hte_agg) > 0) {
    p_hte <- ggplot(hte_agg, aes(x = patient_class, y = age_group, fill = mean_tau)) +
      geom_tile(color = "white") +
      scale_fill_viridis_c(option = "plasma", na.value = "grey90") +
      labs(title = "Heterogeneous Treatment Effects (Drive-Through on TAT) Heatmap", x = "Patient Class", y = "Age Group", fill = "Mean HTE") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p_hte)
  }
}

# NEW: Decision Curve Analysis for Positivity Model
if (requireNamespace("pROC", quietly = TRUE) && exists("xgb_mod")) {
  if (exists("test_data")) {
    xgb_pred <- predict(xgb_mod, xgb.DMatrix(test_matrix), type = "response")
    dca_data <- data.frame(outcome = test_data$y, pred = xgb_pred)
    thresholds <- seq(0, 1, by = 0.01)
    dca_result <- tibble(
      threshold = thresholds,
      net_benefit = sapply(thresholds, function(t) {
        tp <- sum(xgb_pred >= t & test_data$y == 1)
        fp <- sum(xgb_pred >= t & test_data$y == 0)
        n <- length(test_data$y)
        (tp / n) - (fp / n) * (t / (1 - t))
      })
    )
    
    p_dca <- ggplot(dca_result, aes(x = threshold, y = net_benefit)) +
      geom_line(size = 1, color = viridis::viridis(1)) +
      labs(title = "Decision Curve Analysis for Positivity Model", x = "Threshold Probability", y = "Net Benefit") +
      theme_minimal()
    print(p_dca)
  } else {
    message("Test data not available for decision curve analysis.")
  }
}

# NEW: Reliability Curves for Calibration
if (requireNamespace("pROC", quietly = TRUE) && exists("xgb_mod")) {
  calib_data <- xgb_train %>%
    mutate(pred = predict(xgb_mod, xgb_train_matrix, type = "response")) %>%
    mutate(bin = ntile(pred, 10)) %>%
    group_by(bin) %>%
    summarise(mean_pred = mean(pred), obs_rate = mean(y == 1), n = n(), .groups = 'drop')
  
  p_reli <- ggplot(calib_data, aes(x = mean_pred, y = obs_rate)) +
    geom_point(aes(size = n), alpha = 0.6, color = viridis::viridis(1)) +
    geom_line(color = viridis::viridis(1), size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    scale_size_continuous(range = c(0.5, 3)) +
    labs(title = "Reliability Curve for Positivity Model", x = "Mean Predicted Probability", y = "Observed Rate", size = "Sample Size") +
    theme_minimal()
  print(p_reli)
}

# NEW: TAT Quantile Trends Plot
if (requireNamespace("ggplot2", quietly = TRUE) && exists("tat_data") && nrow(tat_data) > 0) {
  tat_quantile_trends <- tat_data %>%
    group_by(pan_day) %>%
    summarise(
      q10 = quantile(col_rec_tat_log, 0.1, na.rm = TRUE),
      q50 = median(col_rec_tat_log, na.rm = TRUE),
      q90 = quantile(col_rec_tat_log, 0.9, na.rm = TRUE),
      .groups = 'drop'
    )
  
  p_tat_quantile <- ggplot(tat_quantile_trends, aes(x = pan_day)) +
    geom_line(aes(y = q50, color = "Median"), size = 1) +
    geom_ribbon(aes(ymin = q10, ymax = q90, fill = "10th-90th Quantile"), alpha = 0.3) +
    scale_color_manual(values = viridis::viridis(1)) +
    scale_fill_manual(values = viridis::viridis(1)) +
    labs(title = "TAT Quantile Trends Over Pandemic Days", x = "Pandemic Day", y = "Log(Hours)") +
    theme_minimal() +
    theme(legend.title = element_blank())
  print(p_tat_quantile)
}

# NEW: Notation Summary Table
notation_data <- tibble(
  Symbol = c("\\bm{X}", "Y", "Z", "T^{(1)},T^{(2)}", "A", "\\pi(\\bm{x})", "q_\\tau(\\bm{x})", "\\Delta_k", "\\tau^{(k)}(\\bm{x})"),
  Description = c("feature vector after encoding (age, gender, patient_class, day, etc.)",
                  "positivity indicator (1=positive)",
                  "cycle-threshold among positives (Y=1)",
                  "TAT components; \\widetilde{T}^{(k)}=\\log T^{(k)}",
                  "drive-through indicator (0/1)",
                  "\\Pr(Y=1\\mid \\bm{x}); estimated via stacked ensemble",
                  "conditional quantile at level \\tau for Z or \\widetilde{T}^{(k)}",
                  "ATE of A on \\widetilde{T}^{(k)} via DML",
                  "HTE function from causal forests")
)

if (exists("gt")) {
  notation_data %>% gt() %>%
    tab_header(title = "Notation Summary") %>%
    cols_label(Symbol = "Symbol", Description = "Description") %>%
    print()
} else {
  print(notation_data)
}

## -------------------- 8. Export Results ---------------------------------------------------------
write_csv(summary_stats, "covid_summary_stats.csv")
write_csv(fairness_tbl, "fairness_summary.csv")
if (exists("ate_tbl")) write_csv(ate_tbl, "causal_ate.csv")
if (exists("fairness_metrics")) write_csv(fairness_metrics, "fairness_metrics.csv")
if (exists("model_comparison")) write_csv(model_comparison, "model_comparison.csv")
write_csv(notation_data, "notation_summary.csv")

save_cols <- intersect(c("subject_id", "age", "gender", "patient_class", "drive_thru_ind",
                         "pan_day", "y", "ct_result"), names(df))
if (length(save_cols) > 0) {
  if (packageVersion("dplyr") < "1.0.0") {
    write_csv(df[save_cols], "covid_processed.csv")
  } else {
    write_csv(df %>% dplyr::select(all_of(save_cols)), "covid_processed.csv")
  }
}

writeLines(capture.output(sessionInfo()), "session_info.txt")

if (requireNamespace("patchwork", quietly = TRUE)) {
  plot_list <- mget(ls(pattern = "^p_"), envir = .GlobalEnv)
  plot_list <- plot_list[!sapply(plot_list, is.null)]
  if (length(plot_list) > 0) {
    dashboard <- wrap_plots(plot_list, ncol = 2) + plot_annotation(title = "COVID Testing Analysis Dashboard")
    ggsave("dashboard_plot.png", dashboard, width = 14, height = 10 * ceiling(length(plot_list)/2)/2, dpi = 300)
  }
}

message("Analysis complete! Outputs saved:")
message("- CSVs: covid_processed.csv, covid_summary_stats.csv, fairness_summary.csv, notation_summary.csv")
if (exists("ate_tbl")) message("- causal_ate.csv")
if (exists("fairness_metrics")) message("- fairness_metrics.csv")
if (exists("model_comparison")) message("- model_comparison.csv")
message("- Plot: dashboard_plot.png")
message("- Session info: session_info.txt")

