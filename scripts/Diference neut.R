# ============================================================
# Longitudinal Changes in Pseudoneutralization Across SARS-CoV-2 Variants
# Figure: Δ % inhibition by variant and survey interval
# Author: Pablo
# Version: 1.0
# ============================================================


library(dplyr)
library(tidyr)
library(ggplot2)

# dataset

db1 <- readRDS("dataset/192 long.rds")
head(db1)

db <- db1 %>%
  select(idjp, idnova, ag_name, survey, titer) %>%
  pivot_wider(
    id_cols = c(idjp, idnova),
    names_from = c(ag_name, survey),
    values_from = titer,
    names_sep = "."
  )

# Diferences by survey
for (v in variantes) {
  
  if(all(c(paste0(v,".L45"), paste0(v,".L46")) %in% names(db)))
    db[[paste0(v,"_diff_L46_L45")]] <- db[[paste0(v,".L46")]] - db[[paste0(v,".L45")]]
  
  if(all(c(paste0(v,".L46"), paste0(v,".L47")) %in% names(db)))
    db[[paste0(v,"_diff_L47_L46")]] <- db[[paste0(v,".L47")]] - db[[paste0(v,".L46")]]
  
  if(all(c(paste0(v,".L47"), paste0(v,".L48")) %in% names(db)))
    db[[paste0(v,"_diff_L48_L47")]] <- db[[paste0(v,".L48")]] - db[[paste0(v,".L47")]]
  
  if(all(c(paste0(v,".L48"), paste0(v,".L49")) %in% names(db)))
    db[[paste0(v,"_diff_L49_L48")]] <- db[[paste0(v,".L49")]] - db[[paste0(v,".L48")]]
}

# long format
db_long <-
  db %>%
  select(matches("_diff_")) %>%
  pivot_longer(
    everything(),
    names_to = "variant_period",
    values_to = "value"
  ) %>%
  mutate(
    period = sub(".*diff_", "", variant_period),
    variant = sub("_diff_.*", "", variant_period)
  )

# order 
db_long$variant <- factor(
  db_long$variant,
  levels = variantes
)

# time
db_long$period <- factor(
  db_long$period,
  levels = c(
    "L46_L45",
    "L47_L46",
    "L48_L47",
    "L49_L48"
  ),
  labels = c(
    "S2-S1",
    "S3-S2",
    "S4-S3",
    "S5-S4"
  )
)

# Paleta
colores_vivos <- c(
  "Ancestral"="#F39C12",
  "Alfa"="#27AE60",
  "Beta"="#117A65",
  "Delta"="#2980B9",
  "Gama"="#F1C40F",
  "Omicron BA.1"="#8E44AD",
  "Omicron BA.2"="#C0392B",
  "Omicron BA.5.2.1"="#8D6E63",
  "Omicron BF.7"="#34495E",
  "Omicron BQ.1"="#D35400",
  "Omicron BQ.1.1"="#2E7D32",
  "Omicron XBB.1"="#922B21"
)

# Figure
ggplot(
  db_long,
  aes(
    interaction(variant, period),
    value,
    colour = variant,
    fill = variant
  )
) +
  geom_boxplot(alpha = .30, outlier.shape = NA) +
  geom_jitter(width = .20, alpha = .20, size = .8) +
  scale_colour_manual(values = colores_vivos) +
  scale_fill_manual(values = colores_vivos) +
  labs(
    x = "",
    y = expression(Delta~"% Inhibition")
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1)
  )

