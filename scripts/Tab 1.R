# ============================================================
# Table 1
# Author: Pablo
# Version: 1.0
# ============================================================

# --- 1. SETUP ---

library(dplyr)     
library(gtsummary) 
library(magrittr)  

# Set working directory (use relative path for GitHub)
setwd("C:/Users/Pablo/OneDrive/0. MSD Landscape/1_PL_Landscape_MSD")
getwd()


db_v <- readRDS(file = "dataset/db_v_selected.rds")
db_v <- db_v %>%
  distinct(idnova, .keep_all = TRUE)
library(dplyr)

db_v <- db_v %>%
  mutate(
    categoria_edad = case_when(
      idade_calc_check <= 5 ~ "≤5",
      idade_calc_check >= 6  & idade_calc_check <= 18 ~ "6-18",
      idade_calc_check >= 19 & idade_calc_check <= 34 ~ "19-34",
      idade_calc_check >= 35 & idade_calc_check <= 64 ~ "35-64",
      idade_calc_check > 64 ~ ">64",
      TRUE ~ NA_character_
    ),
    categoria_edad = factor(
      categoria_edad,
      levels = c("≤5", "6-18", "19-34", "35-64", ">64")
    )
  )


db_v <- db_v %>%
  mutate(
    n.vacL46 = if_else(n.vacL46 >= 1, "Yes", "No"),
    n.vacL47 = if_else(n.vacL47 >= 1, "Yes", "No"),
    n.vacL48 = if_else(n.vacL48 >= 1, "Yes", "No"),
    n.vacL49 = if_else(n.vacL49 >= 1, "Yes", "No")
  ) %>%
  mutate(
    across(
      starts_with("n.vac"),
      ~ factor(.x, levels = c("No", "Yes"))
    )
  )


db_v <- db_v %>%
  mutate(
    across(
      starts_with("razao_"),
      ~ factor(
        if_else(. >= 0.8, "Positive", "Negative"),
        levels = c("Negative", "Positive")
      )
    )
  )

table(db_v$razao_L45)

tbl1 <-
  db_v %>%
  select(
    categoria_edad,
    cind_sexo,
    n.vacL46,
    n.vacL47,
    n.vacL48,
    n.vacL49,
    razao_L45,
    razao_L46,
    razao_L47
  ) %>%
  tbl_summary(
    type = list(
      starts_with("n.vac") ~ "dichotomous",
      starts_with("razao_") ~ "dichotomous"
    ),
    value = list(
      n.vacL46 ~ "Yes",
      n.vacL47 ~ "Yes",
      n.vacL48 ~ "Yes",
      n.vacL49 ~ "Yes",
      razao_L45 ~ "Positive",
      razao_L46 ~ "Positive",
      razao_L47 ~ "Positive"
    ),
    statistic = list(
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "no",
    label = list(
      categoria_edad ~ "Age categories",
      cind_sexo ~ "Sex",
      n.vacL46 ~ "Vaccination in Survey 2",
      n.vacL47 ~ "Vaccination in Survey 3",
      n.vacL48 ~ "Vaccination in Survey 4",
      n.vacL49 ~ "Vaccination in Survey 5",
      razao_L45 ~ "SARS-CoV-2 IgG OD value in Survey 1",
      razao_L46 ~ "SARS-CoV-2 IgG OD value in Survey 2",
      razao_L47 ~ "SARS-CoV-2 IgG OD value in Survey 3"
    )
  ) %>%
  bold_labels()

tbl1
