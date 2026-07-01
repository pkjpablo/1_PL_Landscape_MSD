# ============================================================
# antibody response landscape volume categorized
# Author: Pablo
# Version: 1.0
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

dfL45.2 <- readRDS("dataset/dfL45.2.rds")
dfL46.2 <- readRDS("dataset/dfL46.2.rds")
dfL47.2 <- readRDS("dataset/dfL47.2.rds")
dfL48.2 <- readRDS("dataset/dfL48.2.rds")
dfL49.2 <- readRDS("dataset/dfL49.2.rds")



db_v <- readRDS(file = "dataset/db_v_selected.rds")

db_v <- db_v %>%
  distinct(idnova, .keep_all = TRUE)


dfL45.2 <- dfL45.2 %>%
  mutate(
    volume_group = factor(
      if_else(volumen_percent < 50, "<50%", "≥50%"),
      levels = c("<50%", "≥50%")
    )
  )

db_v <- merge(db_v, dfL45.2, by.x = "idnova",by.y = "id", all.x = TRUE)

names(db_v)

db_v$volume_group <- as.factor(db_v$volume_group)
db_v$cind_sexo <- as.factor(db_v$cind_sexo)
db_v$categoria_edad <- as.factor(db_v$categoria_edad)
db_v$categoria_edad <- factor(
  db_v$categoria_edad,
  levels = c("<18", "18-29", "30-44", "45-59", "≥60")
)

tbl1 <- 
  db_v %>%
  dplyr::select(
    volume_group,
    cind_sexo,
    categoria_edad
      ) %>%
  tbl_summary(
    by = volume_group,
    statistic = list(
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "no",
    digits = list(
      all_categorical() ~ c(0, 0, 1)  # n:0 decimales, N:0, p:1 decimal
    )
  ) %>%
  add_p() %>%                 
  bold_labels()

tbl1


db_v <- readRDS(file = "dataset/db_v_selected.rds")
db_v <- db_v %>%
  distinct(idnova, .keep_all = TRUE)


dfL46.2 <- dfL46.2 %>%
  mutate(
    volume_group = factor(
      if_else(volumen_percent < 50, "<50%", "≥50%"),
      levels = c("<50%", "≥50%")
    )
  )


db_v <- merge(db_v, dfL46.2, by.x = "idnova",by.y = "id", all.x = TRUE)

names(db_v)

db_v$volume_group <- as.factor(db_v$volume_group)
db_v$cind_sexo <- as.factor(db_v$cind_sexo)
db_v$categoria_edad <- as.factor(db_v$categoria_edad)
db_v$categoria_edad <- factor(
  db_v$categoria_edad,
  levels = c("<18", "18-29", "30-44", "45-59", "≥60")
)
db_v$n.vacL46 <- as.factor(db_v$n.vacL46)

tbl2 <- 
  db_v %>%
  dplyr::select(
    volume_group,
    cind_sexo,
    categoria_edad,
    razao_L45,
    n.vacL46
  ) %>%
  tbl_summary(
    by = volume_group,
    statistic = list(
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "no",
    digits = list(
      all_categorical() ~ c(0, 0, 1)  # n:0 decimales, N:0, p:1 decimal
    )
  ) %>%
  add_p() %>%                 
  bold_labels()

tbl2




db_v <- readRDS(file = "dataset/db_v_selected.rds")
db_v <- db_v %>%
  distinct(idnova, .keep_all = TRUE)


dfL47.2 <- dfL47.2 %>%
  mutate(
    volume_group = factor(
      if_else(volumen_percent < 50, "<50%", "≥50%"),
      levels = c("<50%", "≥50%")
    )
  )


db_v <- merge(db_v, dfL47.2, by.x = "idnova",by.y = "id", all.x = TRUE)

names(db_v)

db_v$volume_group <- as.factor(db_v$volume_group)
db_v$cind_sexo <- as.factor(db_v$cind_sexo)
db_v$categoria_edad <- as.factor(db_v$categoria_edad)
db_v$categoria_edad <- factor(
  db_v$categoria_edad,
  levels = c("<18", "18-29", "30-44", "45-59", "≥60")
)
db_v$n.vacL47 <- as.factor(db_v$n.vacL47)

tbl3 <- 
  db_v %>%
  dplyr::select(
    volume_group,
    cind_sexo,
    categoria_edad,
    razao_L45,
    razao_L46,
    n.vacL47
  ) %>%
  tbl_summary(
    by = volume_group,
    statistic = list(
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "no",
    digits = list(
      all_categorical() ~ c(0, 0, 1)  # n:0 decimales, N:0, p:1 decimal
    )
  ) %>%
  add_p() %>%                 
  bold_labels()

tbl3



db_v <- readRDS(file = "dataset/db_v_selected.rds")
db_v <- db_v %>%
  distinct(idnova, .keep_all = TRUE)


dfL48.2 <- dfL48.2 %>%
  mutate(
    volume_group = factor(
      if_else(volumen_percent < 50, "<50%", "≥50%"),
      levels = c("<50%", "≥50%")
    )
  )


db_v <- merge(db_v, dfL48.2, by.x = "idnova",by.y = "id", all.x = TRUE)

names(db_v)

db_v$volume_group <- as.factor(db_v$volume_group)
db_v$cind_sexo <- as.factor(db_v$cind_sexo)
db_v$categoria_edad <- as.factor(db_v$categoria_edad)
db_v$categoria_edad <- factor(
  db_v$categoria_edad,
  levels = c("<18", "18-29", "30-44", "45-59", "≥60")
)
db_v$n.vacL48 <- as.factor(db_v$n.vacL48)

tbl4 <- 
  db_v %>%
  dplyr::select(
    volume_group,
    cind_sexo,
    categoria_edad,
    razao_L45,
    razao_L46,
    razao_L47,
    n.vacL48
  ) %>%
  tbl_summary(
    by = volume_group,
    statistic = list(
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "no",
    digits = list(
      all_categorical() ~ c(0, 0, 1)  # n:0 decimales, N:0, p:1 decimal
    )
  ) %>%
  add_p() %>%                 
  bold_labels()

tbl4



db_v <- readRDS(file = "dataset/db_v_selected.rds")
db_v <- db_v %>%
  distinct(idnova, .keep_all = TRUE)


dfL49.2 <- dfL49.2 %>%
  mutate(
    volume_group = factor(
      if_else(volumen_percent < 50, "<50%", "≥50%"),
      levels = c("<50%", "≥50%")
    )
  )


db_v <- merge(db_v, dfL49.2, by.x = "idnova",by.y = "id", all.x = TRUE)

names(db_v)

db_v$volume_group <- as.factor(db_v$volume_group)
db_v$cind_sexo <- as.factor(db_v$cind_sexo)
db_v$categoria_edad <- as.factor(db_v$categoria_edad)
db_v$categoria_edad <- factor(
  db_v$categoria_edad,
  levels = c("<18", "18-29", "30-44", "45-59", "≥60")
)
db_v$n.vacL49 <- as.factor(db_v$n.vacL49)

tbl5 <- 
  db_v %>%
  dplyr::select(
    volume_group,
    cind_sexo,
    categoria_edad,
    razao_L45,
    razao_L46,
    razao_L47,
    n.vacL49
  ) %>%
  tbl_summary(
    by = volume_group,
    statistic = list(
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "no",
    digits = list(
      all_categorical() ~ c(0, 0, 1)  # n:0 decimales, N:0, p:1 decimal
    )
  ) %>%
  add_p() %>%                 
  bold_labels()

tbl5




### Model 

db_v <- db_v %>%
  mutate(
    # Primero convertir a numérico
    n.vacL49_num = as.numeric(as.character(n.vacL49)),
    
    # Luego recategorizar
    vacunas_categorizadas = case_when(
      n.vacL49_num >= 3 ~ "≥3 doses",
      n.vacL49_num >= 1 & n.vacL49_num <= 2 ~ "1–2 doses",
      n.vacL49_num == 0 ~ "0 doses",
      TRUE ~ NA_character_
    )
  )

table(db_v$vacunas_categorizadas, useNA = "ifany")

contrasts(db_v$volume_group)
db_v$volume_group <- forcats::fct_rev(db_v$volume_group)


m1 <- glm(volume_group ~ cind_sexo + categoria_edad + factor(vacunas_categorizadas), 
          data = db_v, 
          family = binomial())



m2 <- glm(volume_group ~  categoria_edad + factor(vacunas_categorizadas),  
          data = db_v,
          family = binomial())

anova(m1,m2)

sjPlot::tab_model(m2)


library(gtsummary)

tbl_regression(
  m2,
  exponentiate = TRUE,
  estimate_fun = \(x) style_ratio(x, digits = 2),
  conf.level = 0.95
) %>%
  modify_column_merge(
    pattern = "{estimate} ({conf.low}, {conf.high})",
    rows = !is.na(estimate)
  ) %>%
  modify_header(estimate ~ "**OR (95% CI)**") %>%
  bold_labels()


tbl_regression(
  m2,
  exponentiate = TRUE,
  estimate_fun = \(x) style_sigfig(x, digits = 1),  
  pvalue_fun = \(x) style_pvalue(x, digits = 3)      
) %>%
  bold_labels()
