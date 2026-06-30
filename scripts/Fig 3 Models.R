# ============================================================
# MODEL - COMPLETE SCRIPT (CLEANED)
# Author: Pablo
# Version: 4.0
# ============================================================

library(dplyr)
library(gtsummary)
library(here)

# --- 1. LOAD DATA ---
db_v <- readRDS("dataset/db_v_selected.rds")  # Adjust path as needed

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
    estimate_fun = ~style_number(.x, digits = 2),
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

