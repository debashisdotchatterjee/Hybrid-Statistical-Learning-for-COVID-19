library(tidyverse)
library(medicaldata)
library(janitor)
library(gt)

if (!requireNamespace("gtsummary", quietly = TRUE)) {
  message("gtsummary not available, using alternatives")
} else {
  library(gtsummary)
}

packages_optional <- c("tableone", "mice", "mgcv", "glmnet", "xgboost", "ranger", 
                       "randomForest", "grf", "quantreg", "quantregForest", "patchwork", 
                       "viridis", "scales", "vip", "DALEX", "fastshap", "recipes", 
                       "rsample", "caret", "caretEnsemble", "skimr", "summarytools", 
                       "pdp", "changepoint", "pROC", "MASS", "mgcViz")
for(pkg in packages_optional) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Package", pkg, "not available, some functionality will be skipped"))
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

## -------------------- 3. Load & Quick Inspect -------------------------------------------------
data("covid_testing", package = "medicaldata")
df <- covid_testing %>% clean_names()

tryCatch({
  print(dplyr::glimpse(df))
}, error = function(e) {
  message("glimpse failed, using str(): ", e$message)
  print(str(df))
})

if(exists("skim")) {
  skimr::skim(df)
} else {
  print(paste("Total observations:", nrow(df)))
  print(paste("Variables:", ncol(df)))
  print("Missing data summary:")
  print(colSums(is.na(df)))
}

## -------------------- 3.1 Missing Data Handling with MICE ---------------------------------------
if(requireNamespace("mice", quietly = TRUE)) {
  library(mice)
  na_summary <- df %>% summarise_all(~sum(is.na(.)))
  vars_with_na <- names(na_summary)[na_summary > 0]
  
  if(length(vars_with_na) > 0) {
    message("Variables with missing data: ", paste(vars_with_na, collapse = ", "))
    
    imp_cols <- intersect(c("age", "gender", "patient_class", "drive_thru_ind", "pan_day", "y", "col_rec_tat", "rec_ver_tat"), names(df))
    method_vec <- sapply(df[imp_cols], function(x) {
      if(is.numeric(x)) "pmm"
      else if(is.factor(x) && nlevels(x) == 2) "logreg"
      else if(is.factor(x)) "polyreg"
      else "pmm"
    })
    
    tryCatch({
      imp <- mice(df[imp_cols], m = 1, maxit = 5, method = method_vec, printFlag = FALSE, seed = 2025)
      df_imp <- complete(imp)
      df <- df %>% select(-all_of(imp_cols)) %>% bind_cols(df_imp)
      
      if("ct_result" %in% names(df)) {
        df <- df %>%
          mutate(ct_result = case_when(
            y == 1 & is.na(ct_result) ~ median(ct_result[y == 1], na.rm = TRUE),
            TRUE ~ ct_result
          ))
      }
      
      message("MICE imputation applied successfully.")
    }, error = function(e) {
      message("MICE imputation failed: ", e$message, ". Using simple imputation.")
      df <- df %>%
        mutate(across(where(is.numeric), ~ifelse(is.na(.), median(., na.rm = TRUE), .)),
               across(where(is.factor) | where(is.character), ~ifelse(is.na(.), "missing", as.character(.))))
    })
  } else {
    message("No missing data detected in selected columns.")
  }
} else {
  message("MICE not available; install with install.packages('mice'). Using simple imputation.")
  df <- df %>%
    mutate(across(where(is.numeric), ~ifelse(is.na(.), median(., na.rm = TRUE), .)),
           across(where(is.factor) | where(is.character), ~ifelse(is.na(.), "missing", as.character(.))))
}

## -------------------- 4. Preprocessing Helpers -------------------------------------------------
df <- df %>%
  mutate(
    result = as.character(result),
    y = case_when(
      result == "positive" ~ 1L,
      result == "negative" ~ 0L,
      TRUE ~ NA_integer_
    ),
    ct_result = as.numeric(ct_result),
    pan_day = as.integer(pan_day),
    drive_thru_ind = as.integer(drive_thru_ind),
    age = as.numeric(age),
    age_group = cut(age, 
                    breaks = c(-Inf, 2, 5, 12, 18, 30, 50, 65, Inf),
                    labels = c("<2", "2-5", "6-12", "13-18", "19-30", "31-50", "51-65", ">=66"),
                    right = FALSE),
    gender = factor(ifelse(is.na(gender), "missing", as.character(gender))),
    patient_class = factor(ifelse(is.na(patient_class), "missing", as.character(patient_class))),
    clinic_name = factor(ifelse(is.na(clinic_name), "missing", as.character(clinic_name))),
    demo_group = factor(ifelse(is.na(demo_group), "missing", as.character(demo_group))),
    payor_group = factor(ifelse(is.na(payor_group), "missing", as.character(payor_group))),
    weekday = factor((pan_day %% 7) + 1, labels = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"))
  )

tat_vars <- c("col_rec_tat", "rec_ver_tat")
for(v in tat_vars) {
  if(v %in% names(df)) {
    df[[paste0(v, "_wins")]] <- winsorize_vec(df[[v]], lower_pct = 0.005, upper_pct = 0.995)
    df[[paste0(v, "_log")]] <- ifelse(df[[paste0(v, "_wins")]] > 0, log(df[[paste0(v, "_wins")]]), NA_real_)
  }
}

message("Preprocessing complete. Shape: ", paste(dim(df), collapse = " x "))

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

if(exists("gt")) {
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

clinic_summary <- df %>%
  group_by(clinic_name) %>%
  summarise(
    n_tests = n(),
    positives = sum(y == 1, na.rm = TRUE),
    pos_rate = mean(y == 1, na.rm = TRUE),
    median_ct = median(ct_result[y == 1], na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(desc(n_tests)) %>%
  slice_head(n = 20) %>%
  mutate(pos_rate = scales::percent(pos_rate, accuracy = 0.1))

print("Top 20 Clinics:")
print(clinic_summary)

if(exists("gt")) {
  clinic_summary %>% 
    gt() %>% 
    tab_header(title = "Top 20 Clinics by Test Volume") %>%
    fmt_percent(columns = "pos_rate") %>%
    cols_label(
      clinic_name = "Clinic",
      n_tests = "# Tests",
      positives = "# Positives",
      pos_rate = "Positivity Rate",
      median_ct = "Median Ct"
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

if(exists("gt")) {
  fairness_tbl %>% gt() %>%
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
}

if(exists("tbl_summary")) {
  tryCatch({
    stratified_table <- df %>% 
      select(any_of(c("age", "gender", "patient_class", "drive_thru_ind", "y"))) %>%
      tbl_summary(by = y, missing = "ifany")
    print(stratified_table)
  }, error = function(e) {
    message("gtsummary table failed: ", e$message)
    for(var in c("age", "gender", "patient_class", "drive_thru_ind")) {
      if(var %in% names(df)) {
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

if(nrow(daily_data) > 0) {
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

if(nrow(age_pos) > 0) {
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
if(nrow(pos_ct) > 10) {
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
if(length(tat_vars_present) > 0) {
  if(packageVersion("dplyr") < "1.0.0") {
    tat_long <- df[tat_vars_present] %>%
      pivot_longer(cols = tat_vars_present, names_to = "tat_type", values_to = "log_hours") %>%
      filter(!is.na(log_hours))
  } else {
    tat_long <- df %>%
      dplyr::select(all_of(tat_vars_present)) %>%
      pivot_longer(cols = all_of(tat_vars_present), names_to = "tat_type", values_to = "log_hours") %>%
      filter(!is.na(log_hours))
  }
  
  if(nrow(tat_long) > 0) {
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

top_clinics <- df %>%
  count(clinic_name, sort = TRUE) %>%
  slice_head(n = min(15, nrow(.))) %>%
  pull(clinic_name)

if(length(top_clinics) > 0) {
  clinic_heat <- df %>%
    filter(clinic_name %in% top_clinics) %>%
    group_by(clinic_name, pan_day) %>%
    summarise(pos_rate = mean(y == 1, na.rm = TRUE), n = n(), .groups = 'drop') %>%
    filter(!is.na(pos_rate))
  
  if(nrow(clinic_heat) > 0) {
    p_heat <- ggplot(clinic_heat, aes(x = pan_day, y = reorder(clinic_name, desc(n)), fill = pos_rate)) +
      geom_tile() +
      scale_fill_viridis_c(option = "plasma", labels = scales::percent_format(accuracy = 1), na.value = "grey90") +
      labs(title = "Positivity Heatmap by Clinic and Day", x = "Pandemic Day", y = "Clinic", fill = "Positivity") +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 8))
    print(p_heat)
  }
}

pos_heat_agg <- df %>%
  group_by(age_group, patient_class) %>%
  summarise(n = n(), pos_rate = mean(y == 1, na.rm = TRUE), .groups = 'drop') %>%
  filter(!is.na(age_group) & !is.na(patient_class))

if(nrow(pos_heat_agg) > 0) {
  p_pos_heat <- ggplot(pos_heat_agg, aes(x = patient_class, y = age_group, fill = pos_rate)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(option = "viridis", labels = percent_format(accuracy = 1), na.value = "grey90") +
    labs(title = "Positivity Heatmap by Age Group and Patient Class", x = "Patient Class", y = "Age Group", fill = "Positivity Rate") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p_pos_heat)
}

ct_by_demo <- df %>% filter(y == 1 & !is.na(ct_result) & !is.na(demo_group))
if(nrow(ct_by_demo) > 0) {
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

if(nrow(tat_by_day) > 0) {
  p_tat_trend <- ggplot(tat_by_day, aes(x = pan_day, y = median_tat)) +
    geom_line(color = viridis(1), size = 1) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = viridis(1), alpha = 0.3) +
    labs(title = "Median Log TAT (Collection->Receipt) by Pandemic Day with 10th-90th Bands",
         x = "Pandemic Day", y = "Log(Hours)") +
    theme_minimal()
  print(p_tat_trend)
}

## -------------------- 7. Modeling ---------------------------------------------------------------
if(exists("gam") && requireNamespace("mgcv", quietly = TRUE)) {
  tryCatch({
    gam_formula <- y ~ s(age, bs = "cs") + s(pan_day, bs = "cs") + gender + drive_thru_ind + patient_class
    gam_model <- mgcv::gam(gam_formula, data = df, family = binomial())
    print(summary(gam_model))
    par(mfrow = c(2, 2))
    plot(gam_model, pages = 1, shade = TRUE, shade.col = "lightblue")
    par(mfrow = c(1, 1))
  }, error = function(e) {
    message("GAM fitting failed: ", e$message)
  })
}

if(requireNamespace("xgboost", quietly = TRUE) && requireNamespace("pdp", quietly = TRUE)) {
  library(xgboost)
  library(pdp)
  
  xgb_features <- c("age", "gender", "patient_class", "drive_thru_ind", "pan_day", "clinic_name")
  xgb_features <- intersect(xgb_features, names(df))
  xgb_cols <- c(xgb_features, "y")
  xgb_cols <- intersect(xgb_cols, names(df))
  
  if(length(xgb_cols) > length(xgb_features)) {
    if(packageVersion("dplyr") < "1.0.0") {
      xgb_train <- df[xgb_cols] %>% filter(complete.cases(.))
    } else {
      xgb_train <- df %>% dplyr::select(all_of(xgb_cols)) %>% filter(complete.cases(.))
    }
    
    if(nrow(xgb_train) > 0 && "age" %in% colnames(xgb_train)) {
      xgb_train <- xgb_train %>% mutate(age = as.numeric(age))
      xgb_train_matrix <- model.matrix(y ~ ., data = xgb_train)[, -1]
      xgb_label <- xgb_train$y
      cat("Columns in xgb_train_matrix:", colnames(xgb_train_matrix), "\n")
      
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
      
      if("age" %in% colnames(xgb_train_matrix)) {
        pdp_age <- partial(xgb_mod, pred.var = "age", train = xgb_train_matrix, type = "classification", prob = TRUE)
        p_pdp_age <- ggplot(pdp_age, aes(x = age, y = yhat)) +
          geom_line(color = viridis::viridis(1), size = 1) +
          labs(title = "Partial Dependence Plot for Age in Positivity Model",
               x = "Age (Years)", y = "Predicted Positivity Probability") +
          theme_minimal()
        print(p_pdp_age)
      } else {
        message("Age column not found in xgb_train_matrix for PDP.")
      }
      
      if(requireNamespace("vip", quietly = TRUE) && requireNamespace("xgboost", quietly = TRUE)) {
        library(vip)
        library(xgboost)
        
        xgb_predict <- function(model, newdata) {
          if (!inherits(newdata, "xgb.DMatrix")) {
            newdata <- xgb.DMatrix(newdata)
          }
          predict(model, newdata, type = "response")
        }
        
        vi_shap <- vip::vi(
          xgb_mod,
          method = "shap",
          train = xgb_train_matrix,
          pred_wrapper = xgb_predict
        )
        
        p_shap <- vi_shap %>%
          ggplot(aes(x = reorder(Variable, Importance), y = Importance)) +
          geom_bar(stat = "identity", fill = viridis::viridis(1)) +
          coord_flip() +
          labs(title = "SHAP Feature Importance for Positivity Model", x = "Feature", y = "SHAP Importance") +
          theme_minimal()
        print(p_shap)
      } else {
        message("vip or xgboost not available for SHAP analysis.")
      }
    } else {
      message("Insufficient complete data or age missing for XGBoost modeling.")
    }
  } else {
    message("Required columns (including y) not found for XGBoost modeling.")
  }
}

if(requireNamespace("quantreg", quietly = TRUE)) {
  library(quantreg)
  ct_data <- df %>% filter(y == 1 & !is.na(ct_result) & !is.na(age))
  if(nrow(ct_data) > 10) {
    qr_mod <- rq(ct_result ~ age, tau = c(0.1, 0.5, 0.9), data = ct_data)
    qr_pred <- predict(qr_mod, newdata = data.frame(age = seq(min(ct_data$age, na.rm = TRUE), max(ct_data$age, na.rm = TRUE), length.out = 100)))
    qr_df <- data.frame(age = seq(min(ct_data$age, na.rm = TRUE), max(ct_data$age, na.rm = TRUE), length.out = 100),
                        q10 = qr_pred[, 1], q50 = qr_pred[, 2], q90 = qr_pred[, 3])
    
    p_qr_ct <- ggplot(qr_df, aes(x = age)) +
      geom_line(aes(y = q50, color = "Median"), size = 1) +
      geom_ribbon(aes(ymin = q10, ymax = q90, fill = "10th-90th Quantile"), alpha = 0.3) +
      scale_color_manual(values = viridis::viridis(1)) +
      scale_fill_manual(values = viridis::viridis(1)) +
      labs(title = "Quantile Regression for Ct Values by Age", x = "Age (Years)", y = "Ct Value") +
      theme_minimal() +
      theme(legend.title = element_blank())
    print(p_qr_ct)
  }
}

if(requireNamespace("pROC", quietly = TRUE) && exists("xgb_mod") && exists("xgb_train")) {
  library(pROC)
  
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
  
  if(nrow(calib_by_gender) > 0) {
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

if(requireNamespace("grf", quietly = TRUE)) {
  library(grf)
  causal_cols <- c("col_rec_tat_log", "drive_thru_ind", "age", "gender", "patient_class", "pan_day")
  causal_cols <- intersect(causal_cols, names(df))
  
  if(packageVersion("dplyr") < "1.0.0") {
    causal_data <- df[causal_cols] %>%
      filter(!is.na(col_rec_tat_log) & !is.na(drive_thru_ind))
  } else {
    causal_data <- df %>%
      dplyr::select(all_of(causal_cols)) %>%
      filter(!is.na(col_rec_tat_log) & !is.na(drive_thru_ind))
  }
  
  if(nrow(causal_data) > 0 && all(c("col_rec_tat_log", "drive_thru_ind") %in% names(causal_data))) {
    n_treated <- sum(causal_data$drive_thru_ind == 1, na.rm = TRUE)
    n_control <- sum(causal_data$drive_thru_ind == 0, na.rm = TRUE)
    message("Causal data: ", nrow(causal_data), " rows, ", n_treated, " treated, ", n_control, " control")
    
    if(n_treated > 10 && n_control > 10 && all(colSums(is.na(causal_data)) == 0)) {
      X <- model.matrix(~ age + gender + patient_class + pan_day - 1, data = causal_data)
      Y <- causal_data$col_rec_tat_log
      W <- causal_data$drive_thru_ind
      
      grf_version <- packageVersion("grf")
      message("Using grf version: ", grf_version)
      if(grf_version < "2.0.0") {
        message("grf version < 2.0.0; applying manual propensity score filtering.")
        prop_model <- glm(drive_thru_ind ~ age + gender + patient_class + pan_day, 
                          family = binomial(), data = causal_data)
        prop_scores <- predict(prop_model, type = "response")
        overlap_idx <- prop_scores > 0.1 & prop_scores < 0.9
        if(sum(overlap_idx) > 10) {
          causal_data <- causal_data[overlap_idx, ]
          X <- X[overlap_idx, , drop = FALSE]
          Y <- Y[overlap_idx]
          W <- W[overlap_idx]
          message("After propensity filtering: ", sum(overlap_idx), " rows remain.")
        } else {
          message("Insufficient data after propensity filtering; skipping causal forest.")
          ate_tbl <- tibble(
            TAT_Component = "Collection to Receipt (Log Hours)",
            ATE = NA_real_,
            SE = NA_real_
          ) %>%
            mutate(ATE = round(ATE, 3), SE = round(SE, 3))
          print(ate_tbl)
          next
        }
      }
      
      cf <- causal_forest(X, Y, W, num.trees = 500, seed = 2025)
      causal_ate <- average_treatment_effect(cf)
      
      ate_tbl <- tibble(
        TAT_Component = "Collection to Receipt (Log Hours)",
        ATE = causal_ate["estimate"],
        SE = causal_ate["std.err"]
      ) %>%
        mutate(ATE = round(ATE, 3), SE = round(SE, 3))
      
      if(exists("gt")) {
        gt_table <- ate_tbl %>% gt() %>%
          tab_header(title = "Average Treatment Effect (Drive-Through on TAT)") %>%
          cols_label(
            TAT_Component = "TAT Component",
            ATE = "ATE (Log Hours)",
            SE = "Standard Error"
          )
        
        if(!is.na(ate_tbl$ATE) && is.finite(ate_tbl$ATE)) {
          gt_table <- gt_table %>%
            data_color(columns = "ATE", palette = "viridis", direction = "column")
        } else {
          message("ATE is NA or non-finite; skipping color formatting.")
        }
        print(gt_table)
      }
      print(ate_tbl)
    } else {
      message("Insufficient treated (", n_treated, ") or control (", n_control, ") observations, or missing covariates.")
    }
  } else {
    message("Insufficient data or required columns missing for causal forest analysis.")
  }
}

if(requireNamespace("changepoint", quietly = TRUE)) {
  library(changepoint)
  cp_data <- df %>%
    filter(!is.na(pan_day) & !is.na(y)) %>%
    arrange(pan_day) %>%
    group_by(pan_day) %>%
    summarise(pos_rate = mean(y == 1, na.rm = TRUE), .groups = 'drop')
  
  if(nrow(cp_data) > 10) {
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

if(requireNamespace("caret", quietly = TRUE) && exists("xgb_mod")) {
  library(caret)
  calib_data <- xgb_train %>% mutate(pred = predict(xgb_mod, xgb_train_matrix, type = "response"))
  if(nrow(calib_data) > 0) {
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

if(exists("xgb_mod") && requireNamespace("pROC", quietly = TRUE)) {
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
  
  if(exists("gt")) {
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

## -------------------- 8. Export Results ---------------------------------------------------------
write_csv(summary_stats, "covid_summary_stats.csv")
if(exists("clinic_summary")) write_csv(clinic_summary, "clinic_summary.csv")
write_csv(fairness_tbl, "fairness_summary.csv")
if(exists("ate_tbl")) write_csv(ate_tbl, "causal_ate.csv")
if(exists("fairness_metrics")) write_csv(fairness_metrics, "fairness_metrics.csv")

save_cols <- c("subject_id", "age", "gender", "patient_class", "drive_thru_ind", "pan_day", "y", "ct_result", "clinic_name")
save_cols <- intersect(save_cols, names(df))
if(length(save_cols) > 0) {
  if(packageVersion("dplyr") < "1.0.0") {
    write_csv(df[save_cols], "covid_processed.csv")
  } else {
    write_csv(df %>% dplyr::select(all_of(save_cols)), "covid_processed.csv")
  }
}

writeLines(capture.output(sessionInfo()), "session_info.txt")

message("Analysis complete! Check generated CSV files and plots.")
message("Processed data saved to: covid_processed.csv")
message("Summary statistics: covid_summary_stats.csv")
message("Clinic summary: clinic_summary.csv")
message("Fairness summary: fairness_summary.csv")
if(exists("ate_tbl")) message("Causal ATE: causal_ate.csv")
if(exists("fairness_metrics")) message("Fairness metrics: fairness_metrics.csv")

