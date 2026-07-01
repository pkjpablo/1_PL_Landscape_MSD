# ============================================================
# Longitudinal pseudoneutralization trajectories
# Author: Pablo
# Version: 1
# ============================================================

# ----------------------------
# Load libraries
# ----------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# ----------------------------
# Load and prepare data
# ----------------------------

db_long <- readRDS("dataset/192 long.rds") %>%
  rename(
    pseudoneut = titer,
    variant = ag_name
  ) %>%
  mutate(
    pseudoneut = as.numeric(pseudoneut)
  ) %>%
  separate(
    sr_name,
    into = c("survey", "idjp"),
    sep = "_",
    remove = FALSE
  ) %>%
  mutate(
    survey = recode(
      survey,
      L45 = "S1",
      L46 = "S2",
      L47 = "S3",
      L48 = "S4",
      L49 = "S5"
    ),
    survey = factor(
      survey,
      levels = c("S1", "S2", "S3", "S4", "S5")
    )
  )

# ----------------------------
# Classify change between S4 and S5
# ----------------------------

trend_df <- db_long %>%
  filter(survey %in% c("S4", "S5")) %>%
  select(idjp, variant, survey, pseudoneut) %>%
  pivot_wider(
    names_from = survey,
    values_from = pseudoneut
  ) %>%
  mutate(
    trend = case_when(
      S5 > S4 ~ "Increase",
      S5 < S4 ~ "Decrease",
      TRUE    ~ "Stable"
    )
  )

db_long <- db_long %>%
  left_join(
    trend_df %>%
      select(idjp, variant, trend),
    by = c("idjp", "variant")
  )

# ----------------------------
# Figure
# Grey lines: complete longitudinal trajectory
# Colored segment: change from S4 to S5
# ----------------------------

figure <- ggplot(
  db_long,
  aes(
    x = survey,
    y = pseudoneut,
    group = idjp
  )
) +
  geom_line(
    color = "black",
    linewidth = 0.4,
    alpha = 0.5
  ) +
  geom_point(
    color = "black",
    size = 0.8,
    alpha = 0.5
  ) +
  geom_line(
    data = filter(db_long, survey %in% c("S4", "S5")),
    aes(color = trend),
    linewidth = 0.9
  ) +
  geom_point(
    data = filter(db_long, survey %in% c("S4", "S5")),
    aes(color = trend),
    size = 1.2
  ) +
  facet_wrap(~variant) +
  scale_color_manual(
    values = c(
      Increase = "steelblue",
      Decrease = "firebrick",
      Stable = "black"
    )
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(-0.5, 1)
  ) +
  labs(
    x = "Survey",
    y = "Pseudoneutralization (%)",
    color = "Change from S4 to S5"
  ) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

figure

