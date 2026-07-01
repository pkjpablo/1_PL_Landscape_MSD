library(readxl)
PL <- readRDS(file = "dataset/main_cohort")


db_v <- readRDS("dataset/db_v_selected.rds")
db <- subset (db_v[,c("idnova","cind_sexo","idade_calc_check","survey")],survey=="L45")



library(dplyr)
library(gtsummary)

PL_sens <- PL %>%
  mutate(
    in_db = if_else(idnova %in% db$idnova, "Included in db", "Not included"),
    age_cat = case_when(
      cens_idade < 18                  ~ "<18",
      cens_idade >= 18 & cens_idade <= 29 ~ "18–29",
      cens_idade >= 30 & cens_idade <= 44 ~ "30–44",
      cens_idade >= 45 & cens_idade <= 59 ~ "45–59",
      cens_idade >= 60                 ~ "≥60"
    )
  )

tbl_sens <- PL_sens %>%
  select(
    in_db,
    cens_idade,
    cind_sexo
  ) %>%
  tbl_summary(
    by = in_db,
    statistic = list(
      cens_idade ~ "{median} ({p25}, {p75})",
      cind_sexo  ~ "{n} ({p}%)"
    ),
    label = list(
      cens_idade ~ "Age (years)",
      cind_sexo  ~ "Sex"
    ),
    missing = "no"
  ) %>%
  add_p(
    test = list(
      cens_idade ~ "wilcox.test",
      cind_sexo  ~ "chisq.test"
    )
  ) %>%
  modify_header(label ~ "**Variable**") %>%
  modify_spanning_header(
    c("stat_1", "stat_2") ~ "**db inclusion status**"
  ) %>%
  bold_labels()

tbl_sens

