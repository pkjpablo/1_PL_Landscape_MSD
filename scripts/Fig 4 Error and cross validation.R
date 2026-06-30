# ============================================================
# SARS-CoV-2 Antigenic Landscape Analysis
# Author: Pablo
# Version: 4.0 (Journal Ready)
# ============================================================

# --- 1. SETUP ---

library(scales)
library(tidyr)
library(dplyr)
library(mgcv)
library(plotly)
library(stringr)
library(ggplot2)
library(forcats)
library(cowplot)

# Set working directory (use relative path for GitHub)
setwd("C:/Users/Pablo/OneDrive/0. MSD Landscape/1_PL_Landscape_MSD")
getwd()

# --- 2. LOAD DATA ---

# Antigenic coordinates (Rössler cartography)
df_coords <- readRDS("dataset/df_coords_processed.rds")

# Longitudinal antibody titer data (Pau da Lima, Salvador)
db_long <- readRDS("dataset/db_long_processed.rds")

# --- 4. CORE FUNCTION: Leave-One-Out Cross-Validation ---

#' Evaluate variant prediction using GAM with leave-one-out cross-validation
#'
#' @param survey_name Character. Survey identifier (e.g., "L45")
#' @param db_long Data frame. Longitudinal antibody titer data
#' @param df_coords Data frame. Antigenic coordinates
#' @return Data frame with prediction errors by variant
evaluate_variant_prediction <- function(survey_name, db_long, df_coords) {
  
  # Subset survey data
  db_sub <- db_long %>% filter(str_starts(sr_name, survey_name))
  merged <- merge(db_sub, df_coords, by = "ag_name", all.x = TRUE)
  
  # Average z (antibody titer) by variant and position
  data_avg <- merged %>%
    mutate(z = as.numeric(titer)) %>%
    group_by(variant = ag_name, x, y) %>%
    summarise(z = mean(z, na.rm = TRUE), .groups = "drop")
  
  # Full GAM with all variants
  gam_full <- gam(z ~ s(x, y, k = 11), data = data_avg)
  
  results <- list()
  variants <- unique(data_avg$variant)
  
  for (v in variants) {
    # Position of the variant
    pos_var <- data_avg %>% filter(variant == v) %>% select(x, y, z)
    
    # Prediction using full model
    pred_full <- predict(gam_full, newdata = pos_var)
    
    # Exclude variant and refit model
    data_train <- data_avg %>% filter(variant != v)
    gam_excl <- gam(z ~ s(x, y, k = 11), data = data_train)
    
    # Prediction at excluded variant position
    pred_excl <- predict(gam_excl, newdata = pos_var)
    
    df_res <- data.frame(
      Survey = survey_name,
      Variant = v,
      Real_z = pos_var$z,
      Predicted_z_full = pred_full,
      Predicted_z_excl = pred_excl
    )
    
    # Calculate errors
    df_res$Abs_Error <- abs(df_res$Predicted_z_excl - df_res$Real_z)
    df_res$dif <- df_res$Predicted_z_excl - df_res$Real_z
    df_res$Rel_Error <- df_res$Abs_Error / df_res$Real_z * 100
    df_res$Rel_Error_signed <- (df_res$dif / df_res$Real_z) * 100
    
    results[[v]] <- df_res
  }
  
  do.call(rbind, results)
}

# --- 5. EXECUTE ANALYSIS ---

# 5.1 Define surveys
surveys <- c("L45", "L46", "L47", "L48", "L49")

# 5.2 Run cross-validation for all surveys
all_results <- lapply(surveys, function(sv) 
  evaluate_variant_prediction(sv, db_long, df_coords))

# 5.3 Combine results
combined_results <- bind_rows(all_results)

# 5.4 Define variant order (for plotting)
variant_order <- c(
  "Ancestral", "Alpha", "Beta", "Gamma", "Delta",
  "BA.1", "BA.2", "BA.5", "BF.7", "BQ.1", "BQ.1.1",
  "XBB.1"
)

combined_results <- combined_results %>%
  mutate(Variant = factor(Variant, levels = variant_order)) %>%
  arrange(Variant)

# --- 6. CLASSIFY VARIANTS BY TEMPORAL CATEGORY ---

variant_classification <- list(
  L45 = list(
    PastCurrent = c("Ancestral", "Alpha", "Beta"),
    Future = c("Gamma", "Delta", "BA.1", "BA.2", "BA.5", "BF.7", "BQ.1", "BQ.1.1", "XBB.1")
  ),
  L46 = list(
    PastCurrent = c("Gamma", "Delta", "Beta"),
    Future = c("BA.1", "BA.2", "BA.5", "BF.7", "BQ.1", "BQ.1.1", "XBB.1")
  ),
  L47 = list(
    PastCurrent = c("BA.1", "BA.2", "Delta"),
    Future = c("BA.5", "BF.7", "BQ.1", "BQ.1.1", "XBB.1")
  ),
  L48 = list(
    PastCurrent = c("BA.5", "BF.7", "BQ.1", "BQ.1.1"),
    Future = c("XBB.1")
  ),
  L49 = list(
    PastCurrent = c("XBB.1", "BQ.1", "BQ.1.1"),
    Future = c()
  )
)

# Assign temporal groups
combined_results$Group <- NA
for (survey_name in names(variant_classification)) {
  past <- variant_classification[[survey_name]]$PastCurrent
  future <- variant_classification[[survey_name]]$Future
  
  combined_results$Group[combined_results$Survey == survey_name &
                           combined_results$Variant %in% past] <- "Past/Current"
  combined_results$Group[combined_results$Survey == survey_name &
                           combined_results$Variant %in% future] <- "Future"
}

combined_results$Group <- factor(combined_results$Group,
                                 levels = c("Past/Current", "Future"))
combined_results$Group[is.na(combined_results$Group)] <- "Past/Current"

# --- 7. RENAME SURVEYS FOR PUBLICATION ---

combined_results$Survey <- factor(
  combined_results$Survey,
  levels = c("L45", "L46", "L47", "L48", "L49"),
  labels = c("Survey 1", "Survey 2", "Survey 3", "Survey 4", "Survey 5")
)

# --- 8. GENERATE FIGURES ---

# 8.1 Figure 1: Prediction error (absolute difference)
plot_error <- ggplot(combined_results,
                     aes(x = fct_rev(Variant), y = dif, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 0.5) +
  coord_flip() +
  facet_wrap(~ Survey, scales = "free_y", nrow = 1) +
  scale_color_manual(values = c("Past/Current" = "#1f77b4",
                                "Future" = "#ff7f0e")) +
  labs(
    title = "",
    y = "Prediction Error (Difference)",
    x = "",
    color = "Variant Category"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    legend.position = "left",
    plot.margin = margin(10, 10, 10, 10)
  )

# 8.2 Figure 2: Relative error (signed percentage)
plot_rel_error <- ggplot(combined_results,
                         aes(x = fct_rev(Variant), y = Rel_Error_signed/100, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 0.5) +
  coord_flip() +
  facet_wrap(~ Survey, scales = "free_y", nrow = 1) +
  scale_color_manual(values = c("Past/Current" = "#1f77b4",
                                "Future" = "#ff7f0e")) +
  scale_y_continuous(
    breaks = seq(
      floor(min(combined_results$Rel_Error_signed/100, na.rm = TRUE)),
      ceiling(max(combined_results$Rel_Error_signed/100, na.rm = TRUE)),
      by = 0.5
    ),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    title = "",
    y = "Relative Prediction Error (%)",
    x = "",
    color = "Variant Category"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    legend.position = "left",
    plot.margin = margin(10, 10, 10, 10)
  )

# 8.3 Figure 3: Observed vs Predicted values
plot_df <- combined_results %>%
  pivot_longer(
    cols = c(Real_z, Predicted_z_excl),
    names_to = "Type",
    values_to = "Z_value"
  )

plot_obs_pred <- ggplot(plot_df,
                        aes(x = fct_rev(Variant), y = Z_value, color = Group, shape = Type)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  coord_flip() +
  facet_wrap(~ Survey, scales = "free_y", nrow = 1) +
  scale_color_manual(values = c("Past/Current" = "#1f77b4",
                                "Future" = "#ff7f0e")) +
  scale_shape_manual(values = c(
    "Real_z" = 16,
    "Predicted_z_excl" = 1
  ), labels = c("Observed", "Predicted (LOOCV)")) +
  labs(
    title = "",
    y = "%",
    x = "",
    color = "Variant Category",
    shape = "Time status"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    legend.position = "left",
    plot.margin = margin(10, 10, 10, 10)
  )

# --- 9. COMBINED FIGURE ---

combined_figure <- cowplot::plot_grid(
  plot_obs_pred,
  plot_error,
  plot_rel_error,
  ncol = 1,
  align = "v",
  labels = c("A", "B", "C"),
  label_size = 14,
  label_fontface = "bold"
)

# --- 10. DISPLAY AND SAVE ---

print(combined_figure)




# ============================================================
# MSE AND 95% CI CALCULATION FOR LOOCV PREDICTIONS
# ============================================================

library(dplyr)
library(boot)

# --- 1. Calculate MSE per Survey ---

calculate_mse <- function(data, survey_name = NULL) {
  
  # Filter by survey if specified
  if (!is.null(survey_name)) {
    data <- data %>% filter(Survey == survey_name)
  }
  
  # Calculate MSE
  mse <- data %>%
    summarise(
      MSE = mean((Real_z - Predicted_z_excl)^2, na.rm = TRUE),
      RMSE = sqrt(MSE),
      MAE = mean(abs(Real_z - Predicted_z_excl), na.rm = TRUE),
      R_squared = 1 - sum((Real_z - Predicted_z_excl)^2, na.rm = TRUE) / 
        sum((Real_z - mean(Real_z, na.rm = TRUE))^2, na.rm = TRUE),
      N = n()
    )
  
  return(mse)
}

# --- 2. Bootstrap function for 95% CI ---

bootstrap_mse <- function(data, n_boot = 1000, survey_name = NULL) {
  
  # Filter by survey if specified
  if (!is.null(survey_name)) {
    data <- data %>% filter(Survey == survey_name)
  }
  
  # Bootstrap function
  boot_fn <- function(data, indices) {
    d <- data[indices, ]
    mean((d$Real_z - d$Predicted_z_excl)^2, na.rm = TRUE)
  }
  
  # Run bootstrap
  boot_results <- boot(data = data, statistic = boot_fn, R = n_boot)
  
  # Calculate 95% CI
  ci <- boot.ci(boot_results, type = "perc")
  
  # Return results
  results <- data.frame(
    MSE = boot_results$t0,
    MSE_lower = ci$percent[4],
    MSE_upper = ci$percent[5],
    Boot_rep = n_boot
  )
  
  return(results)
}

# --- 3. Calculate MSE by Variant Category ---

calculate_mse_by_group <- function(data, survey_name = NULL) {
  
  if (!is.null(survey_name)) {
    data <- data %>% filter(Survey == survey_name)
  }
  
  mse_by_group <- data %>%
    group_by(Group) %>%
    summarise(
      MSE = mean((Real_z - Predicted_z_excl)^2, na.rm = TRUE),
      RMSE = sqrt(MSE),
      MAE = mean(abs(Real_z - Predicted_z_excl), na.rm = TRUE),
      R_squared = 1 - sum((Real_z - Predicted_z_excl)^2, na.rm = TRUE) / 
        sum((Real_z - mean(Real_z, na.rm = TRUE))^2, na.rm = TRUE),
      N = n(),
      .groups = "drop"
    )
  
  return(mse_by_group)
}

# --- 4. Bootstrapped MSE by Group ---

bootstrap_mse_by_group <- function(data, n_boot = 1000) {
  
  groups <- unique(data$Group)
  results <- list()
  
  for (g in groups) {
    group_data <- data %>% filter(Group == g)
    
    boot_fn <- function(d, indices) {
      d_sub <- d[indices, ]
      mean((d_sub$Real_z - d_sub$Predicted_z_excl)^2, na.rm = TRUE)
    }
    
    boot_results <- boot(data = group_data, statistic = boot_fn, R = n_boot)
    ci <- boot.ci(boot_results, type = "perc")
    
    results[[g]] <- data.frame(
      Group = g,
      MSE = boot_results$t0,
      MSE_lower = ci$percent[4],
      MSE_upper = ci$percent[5],
      N = nrow(group_data),
      Boot_rep = n_boot
    )
  }
  
  return(bind_rows(results))
}

# --- 5. Calculate 95% CI for Observed vs Predicted ---

calculate_ci <- function(data, survey_name = NULL) {
  
  if (!is.null(survey_name)) {
    data <- data %>% filter(Survey == survey_name)
  }
  
  # Observed values
  obs_ci <- data %>%
    summarise(
      Observed_Mean = mean(Real_z, na.rm = TRUE),
      Observed_SD = sd(Real_z, na.rm = TRUE),
      Observed_SE = Observed_SD / sqrt(n()),
      Observed_CI_lower = Observed_Mean - 1.96 * Observed_SE,
      Observed_CI_upper = Observed_Mean + 1.96 * Observed_SE,
      N = n()
    )
  
  # Predicted values
  pred_ci <- data %>%
    summarise(
      Predicted_Mean = mean(Predicted_z_excl, na.rm = TRUE),
      Predicted_SD = sd(Predicted_z_excl, na.rm = TRUE),
      Predicted_SE = Predicted_SD / sqrt(n()),
      Predicted_CI_lower = Predicted_Mean - 1.96 * Predicted_SE,
      Predicted_CI_upper = Predicted_Mean + 1.96 * Predicted_SE,
      N = n()
    )
  
  # Combine
  results <- cbind(
    Survey = ifelse(is.null(survey_name), "All", survey_name),
    obs_ci,
    pred_ci
  )
  
  return(results)
}

# --- 6. Calculate MSE by Variant ---

calculate_mse_by_variant <- function(data, survey_name = NULL) {
  
  if (!is.null(survey_name)) {
    data <- data %>% filter(Survey == survey_name)
  }
  
  mse_by_variant <- data %>%
    group_by(Variant) %>%
    summarise(
      MSE = mean((Real_z - Predicted_z_excl)^2, na.rm = TRUE),
      RMSE = sqrt(MSE),
      MAE = mean(abs(Real_z - Predicted_z_excl), na.rm = TRUE),
      R_squared = 1 - sum((Real_z - Predicted_z_excl)^2, na.rm = TRUE) / 
        sum((Real_z - mean(Real_z, na.rm = TRUE))^2, na.rm = TRUE),
      N = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(MSE))
  
  return(mse_by_variant)
}

# ============================================================
# EXECUTE MSE AND CI CALCULATIONS
# ============================================================

# --- 1. Overall MSE ---
mse_overall <- calculate_mse(combined_results)
print("=== OVERALL MSE ===")
print(mse_overall)

# --- 2. MSE by Survey ---
mse_by_survey <- combined_results %>%
  group_by(Survey) %>%
  summarise(
    MSE = mean((Real_z - Predicted_z_excl)^2, na.rm = TRUE),
    RMSE = sqrt(MSE),
    MAE = mean(abs(Real_z - Predicted_z_excl), na.rm = TRUE),
    N = n(),
    .groups = "drop"
  )

print("=== MSE BY SURVEY ===")
print(mse_by_survey)

# --- 3. Bootstrapped MSE by Survey ---
mse_boot_by_survey <- list()
for (s in unique(combined_results$Survey)) {
  mse_boot_by_survey[[s]] <- bootstrap_mse(combined_results, n_boot = 1000, survey_name = s)
}
mse_boot_df <- bind_rows(mse_boot_by_survey, .id = "Survey")

print("=== BOOTSTRAPPED MSE WITH 95% CI BY SURVEY ===")
print(mse_boot_df)

# --- 4. MSE by Temporal Group ---
mse_by_group <- calculate_mse_by_group(combined_results)
print("=== MSE BY TEMPORAL GROUP ===")
print(mse_by_group)

# --- 5. Bootstrapped MSE by Group ---
mse_boot_group <- bootstrap_mse_by_group(combined_results, n_boot = 1000)
print("=== BOOTSTRAPPED MSE WITH 95% CI BY GROUP ===")
print(mse_boot_group)

# --- 6. 95% CI for Observed vs Predicted by Survey ---
ci_by_survey <- list()
for (s in unique(combined_results$Survey)) {
  ci_by_survey[[s]] <- calculate_ci(combined_results, survey_name = s)
}
ci_df <- bind_rows(ci_by_survey)

print("=== 95% CI FOR OBSERVED VS PREDICTED BY SURVEY ===")
print(ci_df)

# --- 7. MSE by Variant ---
mse_by_variant <- calculate_mse_by_variant(combined_results)
print("=== MSE BY VARIANT (TOP 10) ===")
print(head(mse_by_variant, 10))

# --- 8. MSE by Survey and Group ---
mse_survey_group <- combined_results %>%
  group_by(Survey, Group) %>%
  summarise(
    MSE = mean((Real_z - Predicted_z_excl)^2, na.rm = TRUE),
    RMSE = sqrt(MSE),
    MAE = mean(abs(Real_z - Predicted_z_excl), na.rm = TRUE),
    N = n(),
    .groups = "drop"
  )

print("=== MSE BY SURVEY AND GROUP ===")
print(mse_survey_group)