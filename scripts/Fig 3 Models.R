# ============================================================
# MODEL - COMPLETE SCRIPT (CLEANED)
# Author: Pablo
# Version: 4.0
# ============================================================

library(dplyr)
library(gtsummary)
library(here)

# --- 1. LOAD DATA ---
db_v <- readRDS("dataset/db_v_selected.rds") 

# --- 2. CREATE 5 SURVEY-SPECIFIC DATABASES ---
db45 <- db_v %>% filter(survey == "L45")
db46 <- db_v %>% filter(survey == "L46")
db47 <- db_v %>% filter(survey == "L47")
db48 <- db_v %>% filter(survey == "L48")
db49 <- db_v %>% filter(survey == "L49")

# --- 3. LOAD VOLUME FILES ---
dfL45.2 <- readRDS("dataset/dfL45.2.rds")
dfL46.2 <- readRDS("dataset/dfL46.2.rds")
dfL47.2 <- readRDS("dataset/dfL47.2.rds")
dfL48.2 <- readRDS("dataset/dfL48.2.rds")
dfL49.2 <- readRDS("dataset/dfL49.2.rds")

# --- 4. FUNCTION: Merge volume data with survey ---
merge_volume <- function(db, df_vol) {
  df_vol %>%
    select(id, volumen_percent) %>%
    right_join(db, by = c("id" = "idnova"))
}

# --- 5. MERGE ALL SURVEYS WITH VOLUMES ---
db45 <- merge_volume(db45, dfL45.2)
db46 <- merge_volume(db46, dfL46.2)
db47 <- merge_volume(db47, dfL47.2)
db48 <- merge_volume(db48, dfL48.2)
db49 <- merge_volume(db49, dfL49.2)

# --- 6. FUNCTION: Prepare factors ---
prepare_factors <- function(data) {
  data %>%
    mutate(
      categoria_edad = factor(categoria_edad, 
                              levels = c("<18", "18-29", "30-44", "45-59", "≥60")),
      cind_sexo = factor(cind_sexo, levels = c("female", "male"))
    )
}

# --- 7. APPLY FACTORS ---
db45 <- prepare_factors(db45)
db46 <- prepare_factors(db46)
db47 <- prepare_factors(db47)
db48 <- prepare_factors(db48)
db49 <- prepare_factors(db49)

# --- 8. FUNCTION: Create model and table ---
create_model_tbl <- function(data, formula, labels_list, caption_text = "") {
  
  # Fit model
  model <- lm(formula, data = data)
  
  # Create table
  tbl <- tbl_regression(
    model,
    intercept = TRUE,
    label = labels_list,
    estimate_fun = ~style_number(.x, digits = 1),
    pvalue_fun = ~style_pvalue(.x, digits = 2)
  ) %>%
    modify_fmt_fun(
      ci = function(x) paste0("(", gsub(", ", " – ", x), ")")
    ) %>%
    modify_header(
      label = "**Characteristic**",
      estimate = "**Beta**",
      ci = "**95% CI**",
      p.value = "**p-value**"
    ) %>%
    modify_caption(caption_text) %>%
    bold_labels() %>%
    italicize_levels()
  
  return(list(model = model, table = tbl))
}

# ============================================================
# MODEL L45
# ============================================================

result_L45 <- create_model_tbl(
  data = db45,
  formula = volumen_percent ~ cind_sexo + categoria_edad,
  labels_list = list(
    cind_sexo = "Sex",
    categoria_edad = "Age category"
  ),
  caption_text = ""
)

tbl_L45 <- result_L45$table
tbl_L45

# ============================================================
# MODEL L46
# ============================================================

# Convert vaccine doses to factor
db46$n.vacL46 <- factor(db46$n.vacL46)

result_L46 <- create_model_tbl(
  data = db46,
  formula = volumen_percent ~ cind_sexo + categoria_edad + razao_L45 + n.vacL46,
  labels_list = list(
    cind_sexo = "Sex",
    categoria_edad = "Age category",
    razao_L45 = "SARS-CoV-2 IgG OD value (Survey 1)",
    n.vacL46 = "No. of vaccine doses"
  ),
  caption_text = "**Table X. Association between demographic characteristics and SARS-CoV-2 serological markers**"
)

tbl_L46 <- result_L46$table
tbl_L46

# ============================================================
# MODEL L47
# ============================================================

db47$n.vacL47 <- factor(db47$n.vacL47)

result_L47 <- create_model_tbl(
  data = db47,
  formula = volumen_percent ~ cind_sexo + categoria_edad + razao_L45 + razao_L46 + n.vacL47,
  labels_list = list(
    cind_sexo = "Sex",
    categoria_edad = "Age category",
    razao_L45 = "SARS-CoV-2 IgG OD value (Survey 1)",
    razao_L46 = "SARS-CoV-2 IgM OD value (Survey 2)",
    n.vacL47 = "No. of vaccine doses"
  ),
  caption_text = "**Table X. Association between demographic characteristics and SARS-CoV-2 serological markers**"
)

tbl_L47 <- result_L47$table
tbl_L47

# ============================================================
# MODEL L48
# ============================================================

db48$n.vacL48 <- factor(db48$n.vacL48)

result_L48 <- create_model_tbl(
  data = db48,
  formula = volumen_percent ~ cind_sexo + categoria_edad + razao_L45 + razao_L46 + razao_L47 + n.vacL48,
  labels_list = list(
    cind_sexo = "Sex",
    categoria_edad = "Age category",
    razao_L45 = "SARS-CoV-2 IgG OD value (Survey 1)",
    razao_L46 = "SARS-CoV-2 IgM OD value (Survey 2)",
    razao_L47 = "SARS-CoV-2 IgG OD value (Survey 3)",
    n.vacL48 = "No. of vaccine doses"
  ),
  caption_text = "**Table X. Association between demographic characteristics and SARS-CoV-2 serological markers**"
)

tbl_L48 <- result_L48$table
tbl_L48

# ============================================================
# MODEL L49
# ============================================================

db49$n.vacL49 <- factor(db49$n.vacL49)

result_L49 <- create_model_tbl(
  data = db49,
  formula = volumen_percent ~ cind_sexo +categoria_edad + razao_L45 + razao_L46 + razao_L47 + n.vacL49,
  labels_list = list(
    categoria_edad = "Age category",
    razao_L45 = "SARS-CoV-2 IgG OD value (Survey 1)",
    razao_L46 = "SARS-CoV-2 IgM OD value (Survey 2)",
    razao_L47 = "SARS-CoV-2 IgG OD value (Survey 3)",
    n.vacL49 = "No. of vaccine doses"
  ),
  caption_text = "**Table X. Association between demographic characteristics and SARS-CoV-2 serological markers**"
)

tbl_L49 <- result_L49$table
tbl_L49

# ============================================================
# DISPLAY ALL TABLES
# ============================================================

cat("\n", "=", rep("-", 60), "\n")
cat("  TABLES SUMMARY\n")
cat("=", rep("-", 60), "\n\n")

print(tbl_L45)
print(tbl_L46)
print(tbl_L47)
print(tbl_L48)
print(tbl_L49)




# ================================================================
# FOREST PLOT - ORDEN INVERTIDO CON fct_rev()
# Layout: 3 paneles arriba (p1 | p2 | NA) y 3 abajo (p3 | p4 | p5)
# ================================================================

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(patchwork)

# ================================================================
# 1. CREATE DATA FRAME - REFERENCE CATEGORIES CON Beta = 0
# ================================================================

# Survey 1
survey1 <- data.frame(
  Survey = "Survey 1",
  Variable = c("Sex: female (ref)", "Sex: male",
               "Age: <18 (ref)", "Age: 18-29", "Age: 30-44", "Age: 45-59", "Age: ≥60",
               "SARS-CoV-2 IgG OD value (S1)", "SARS-CoV-2 IgG OD value (S2)", "SARS-CoV-2 IgG OD value (S3)",
               "Vaccine doses: 0 (ref)", "Vaccine doses: 1", "Vaccine doses: 2", "Vaccine doses: 3", "Vaccine doses: 4", "Vaccine doses: 5"),
  Beta = c(0, -1.86, 0, -0.24, -0.65, 1.15, -1.50, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  CI_Lower = c(0, -4.05, 0, -3.75, -3.62, -2.47, -4.63, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  CI_Upper = c(0, 0.33, 0, 3.28, 2.31, 4.77, 1.64, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  P_Value = c(NA, "0.10", NA, "0.89", "0.67", "0.53", "0.35", NA, NA, NA, NA, NA, NA, NA, NA, NA)
)

library(ggplot2)
library(dplyr)

# Mantener el orden de las filas
survey1$Variable <- factor(survey1$Variable,
                           levels = rev(survey1$Variable))

# Survey 2
survey2 <- data.frame(
  Survey = "Survey 2",
  Variable = c("Sex: female (ref)", "Sex: male",
               "Age: <18 (ref)", "Age: 18-29", "Age: 30-44", "Age: 45-59", "Age: ≥60",
               "SARS-CoV-2 IgG OD value (S1)", "SARS-CoV-2 IgG OD value (S2)", "SARS-CoV-2 IgG OD value (S3)",
               "Vaccine doses: 0 (ref)", "Vaccine doses: 1", "Vaccine doses: 2", "Vaccine doses: 3", "Vaccine doses: 4", "Vaccine doses: 5"),
  Beta = c(0, -6.70, 0, -2.20, 18.00, 21.60, 17.40, 6.00, NA, NA, 0, 12.60, 0.80, NA, NA, NA),
  CI_Lower = c(0, -13.50, 0, -13.60, 8.90, 7.70, 2.80, 3.80, NA, NA, 0, 1.60, -13.50, NA, NA, NA),
  CI_Upper = c(0, 0.10, 0, 9.20, 27.10, 35.40, 32.00, 8.10, NA, NA, 0, 23.70, 15.10, NA, NA, NA),
  P_Value = c(NA, "0.052", NA, "0.71", "<0.001", "0.002", "0.020", "<0.001", NA, NA, NA, "0.025", "0.91", NA, NA, NA)
)

# Survey 3
survey3 <- data.frame(
  Survey = "Survey 3",
  Variable = c("Sex: female (ref)", "Sex: male",
               "Age: <18 (ref)", "Age: 18-29", "Age: 30-44", "Age: 45-59", "Age: ≥60",
               "SARS-CoV-2 IgG OD value (S1)", "SARS-CoV-2 IgG OD value (S2)", "SARS-CoV-2 IgG OD value (S3)",
               "Vaccine doses: 0 (ref)", "Vaccine doses: 1", "Vaccine doses: 2", "Vaccine doses: 3", "Vaccine doses: 4", "Vaccine doses: 5"),
  Beta = c(0, -6.20, 0, 12.20, 4.10, -6.50, 12.20, 4.20, 2.00, NA, 0, 6.00, 22.10, 22.30, 30.10, NA),
  CI_Lower = c(0, -13.40, 0, 0.90, -6.90, -19.40, -0.60, 1.60, 0.50, NA, 0, -5.60, 12.80, 12.40, 1.50, NA),
  CI_Upper = c(0, 1.00, 0, 23.60, 15.00, 6.50, 25.10, 6.80, 3.50, NA, 0, 17.60, 31.30, 32.10, 58.80, NA),
  P_Value = c(NA, "0.093", NA, "0.035", "0.47", "0.33", "0.062", "0.002", "0.010", NA, NA, "0.31", "<0.001", "<0.001", "0.039", NA)
)

# Survey 4
survey4 <- data.frame(
  Survey = "Survey 4",
  Variable = c("Sex: female (ref)", "Sex: male",
               "Age: <18 (ref)", "Age: 18-29", "Age: 30-44", "Age: 45-59", "Age: ≥60",
               "SARS-CoV-2 IgG OD value (S1)", "SARS-CoV-2 IgG OD value (S2)", "SARS-CoV-2 IgG OD value (S3)",
               "Vaccine doses: 0 (ref)", "Vaccine doses: 1", "Vaccine doses: 2", "Vaccine doses: 3", "Vaccine doses: 4", "Vaccine doses: 5"),
  Beta = c(0, -6.10, 0, -8.60, -14.10, -17.30, -12.50, -0.10, 1.30, 2.70, 0, -4.10, -7.00, 6.10, -0.50, NA),
  CI_Lower = c(0, -12.80, 0, -18.90, -24.30, -28.70, -24.50, -2.50, -0.10, 1.10, 0, -16.20, -17.20, -3.80, -10.60, NA),
  CI_Upper = c(0, 0.60, 0, 1.80, -3.90, -5.90, -0.50, 2.40, 2.70, 4.30, 0, 8.00, 3.20, 16.00, 9.60, NA),
  P_Value = c(NA, "0.074", NA, "0.11", "0.007", "0.003", "0.042", "0.96", "0.062", "<0.001", NA, "0.50", "0.18", "0.23", "0.93", NA)
)

# Survey 5
survey5 <- data.frame(
  Survey = "Survey 5",
  Variable = c("Sex: female (ref)", "Sex: male",
               "Age: <18 (ref)", "Age: 18-29", "Age: 30-44", "Age: 45-59", "Age: ≥60",
               "SARS-CoV-2 IgG OD value (S1)", "SARS-CoV-2 IgG OD value (S2)", "SARS-CoV-2 IgG OD value (S3)",
               "Vaccine doses: 0 (ref)", "Vaccine doses: 1", "Vaccine doses: 2", "Vaccine doses: 3", "Vaccine doses: 4", "Vaccine doses: 5"),
  Beta = c(0, -1.07, 0, 1.34, -11.13, -9.08, -2.22, 0.80, 0.73, 3.27, 0, -6.31, -5.19, 2.29, -3.94, 9.16),
  CI_Lower = c(0, -6.74, 0, -7.64, -19.83, -18.91, -12.41, -1.29, -0.45, 1.92, 0, -16.97, -14.24, -6.38, -13.28, -0.50),
  CI_Upper = c(0, 4.61, 0, 10.31, -2.42, 0.76, 7.96, 2.89, 1.91, 4.61, 0, 4.35, 3.86, 10.95, 5.39, 18.82),
  P_Value = c(NA, "0.71", NA, "0.77", "0.013", "0.070", "0.67", "0.45", "0.23", "<0.001", NA, "0.24", "0.26", "0.60", "0.41", "0.063")
)

# Combine all surveys
forest_data <- bind_rows(survey1, survey2, survey3, survey4, survey5)






library(ggplot2)
library(dplyr)

plot_forest <- function(df, title){
  
  df <- df %>%
    mutate(
      Variable = factor(Variable, levels = rev(Variable)),
      P_num = suppressWarnings(
        as.numeric(gsub("<", "", as.character(P_Value)))
      ),
      Significant = case_when(
        !is.na(P_num) & P_num < 0.05 ~ "Yes",
        TRUE ~ "No"
      )
    ) %>%
    select(-P_num)
  
  ggplot(df, aes(x = Beta, y = Variable)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(
      aes(xmin = CI_Lower, xmax = CI_Upper, colour = Significant),
      height = 0.2,
      na.rm = TRUE
    ) +
    geom_point(
      aes(colour = Significant),
      size = 2.8,
      na.rm = TRUE
    ) +
    scale_color_manual(values = c("Yes" = "red", "No" = "black")) +
    labs(
      title = title,
      x = "Beta (95% CI)",
      y = NULL
    ) +
    coord_cartesian(xlim = c(-30, 50)) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
}


p1 <- plot_forest(survey1, "Survey 1")
p2 <- plot_forest(survey2, "Survey 2")
p3 <- plot_forest(survey3, "Survey 3")
p4 <- plot_forest(survey4, "Survey 4")
p5 <- plot_forest(survey5, "Survey 5")

library(cowplot)

cowplot::plot_grid(
  p1, p2, NULL,
  p3,  p4, p5, 
  ncol = 3
)
