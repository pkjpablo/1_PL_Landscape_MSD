# ============================================================
# COVID-19 Analysis - Salvador, Brazil
# Author: Juan P. Aguilar Ticona
# Version: 4.0 (English version)
# ============================================================

# --- 1. INITIAL SETUP ---
# If needed: setwd(here::here())

# --- 2. LOAD ALL LIBRARIES ---
library(tidyverse)
library(lubridate)
library(scales)
library(cowplot)
library(readxl)

# --- 3. DEFINE CONSTANTS AND REUSABLE FUNCTIONS ---

# Reference dates for vertical lines
REFERENCE_DATES <- c(
  "2020-11-18", "2021-02-26", "2021-07-14",
  "2021-10-09", "2022-03-25", "2022-08-31", 
  "2022-11-16", "2023-05-09", "2023-10-25", 
  "2024-03-24"
)

# Function to add reference lines to plots
add_reference_lines <- function(plot, dates = REFERENCE_DATES) {
  plot +
    geom_vline(
      xintercept = as.numeric(ymd(dates)),
      linetype = "dashed", 
      color = "#82b0d2", 
      size = 0.8
    )
}

# --- 4. LOAD AND SAVE DATA ---

Brazil <- readRDS("dataset/Brazil.rds")
variants <- readRDS("dataset/variants.rds")
vac <- readRDS("dataset/vac.rds")


# --- 5. FIGURE 1: CASES AND DEATHS ---

create_cases_figure <- function(data) {
  salvador <- data %>%
    filter(municipio == "Salvador") %>%
    mutate(DATE = as.Date(data),
           week = floor_date(DATE, "week")) %>%
    group_by(week) %>%
    summarise(
      cases = sum(casosNovos, na.rm = TRUE),
      deaths = sum(obitosNovos, na.rm = TRUE),
      .groups = "drop"
    )
  
  F1 <- ggplot(salvador, aes(x = week)) +
    geom_col(aes(y = cases), fill = "#e9c46a", alpha = 0.75) +
    geom_col(aes(y = deaths * 8), fill = "#dc0000ff", alpha = 0.4) +
    scale_y_continuous(
      name = "Number of cases",
      sec.axis = sec_axis(~./8, name = "Number of deaths")
    ) +
    scale_x_date(
      limits = as.Date(c('2020-03-01', '2024-04-01')),
      breaks = date_breaks("3 months"),
      labels = date_format("%b %y")
    ) +
    theme_classic(base_size = 12) +
    labs(x = "Date")
  
  # Add reference lines
  F1 <- add_reference_lines(F1)
  
  return(F1)
}

Figure1 <- create_cases_figure(Brazil)

# --- 6. FIGURE 2: VARIANTS ---

colores_pastel <- c(
  "Others" = "#FFB3BA", "Ancestral" = "#FFDFBA", "Gamma" = "#FFFFBA",
  "Alpha" = "#BAFFC9", "Delta" = "#BAE1FF", "BA.1*" = "#D7B3FF",
  "BA.2*" = "#FFB3E6", "BA.4*" = "#B3E6FF", "BA.5*" = "#FFC3E1",
  "BQ.1.1*" = "#C3FF93", "BQ.1" = "#FF9E9D", "BQ.1.22" = "#FFC69D",
  "DL.1" = "#FFD89D", "XBB" = "#FFE39D", "XBB.1" = "#E3FF9D",
  "XBB.1.18*" = "#C5FF9D", "FE.1*" = "#9DFFB3", "XBB.1.15.1" = "#9DFFC5",
  "XBB.1.5*" = "#9DFFE3", "XBB.1.9" = "#9DFFFF", "EG.1 or EG.6.1.1" = "#9DE3FF",
  "GK.1*" = "#9DC5FF", "GK.3" = "#9DAFFF", "JD.1*" = "#B39DFF",
  "FL.1.5.1" = "#D19DFF", "XBB.1.16.6" = "#FF9DFF", "JN.1" = "#FF9DC5"
)

linhas <- as.numeric(ymd(c(
  "2019-09-09","2019-11-09","2020-11-18","2021-02-26","2021-07-14",
  "2021-10-09","2022-03-25","2022-08-31","2022-11-16",
  "2023-05-09","2023-10-25","2024-03-24"
)))

# --- Gráfico de variantes ---
Figure2 <- ggplot(variants, aes(x = month, y = percentage, fill = variant)) +
  geom_area() +
  geom_vline(xintercept = linhas, linetype = "dashed", color = "#82b0d2", size = 0.8) +
  scale_fill_manual(values = colores_pastel) +
  scale_x_date(limits = as.Date(c('2020-03-01','2024-03-30')),
               breaks = date_breaks("3 months"),
               labels = date_format("%b %y")) +
  labs(x = "", y = "Variants (%)") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.ticks.x = element_line(size = 1, linetype = "dashed"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "top"
  )


# --- 7. FIGURE 3: VACCINATION ---


# Define vaccination doses labels
dose_labels <- c(
  "1st Dose",
  "2nd Dose", 
  "3rd Dose",
  "4th Dose"
)

# Convert dates to Date format
vac <- vac %>%
  mutate(
    data1 = as.Date(data1),
    data2 = as.Date(data2),
    data3 = as.Date(data3),
    data4 = as.Date(data4)
  )

# Define total population 
TOTAL_POPULATION <- nrow(vac)  

# FUNCTION TO PROCESS VACCINATION DATA

process_vaccine_data <- function(dates, dose_name, total_pop = TOTAL_POPULATION) {
  data.frame(date = dates) %>%
    count(date) %>%
    arrange(date) %>%
    mutate(
      cumulative = cumsum(n),
      percentage = (cumulative / total_pop) * 100,
      dose = dose_name
    )
}

# PROCESS ALL DOSES 

# Process each dose
data_cum1 <- process_vaccine_data(vac$data1, dose_labels[1])
data_cum2 <- process_vaccine_data(vac$data2, dose_labels[2])
data_cum3 <- process_vaccine_data(vac$data3, dose_labels[3])
data_cum4 <- process_vaccine_data(vac$data4, dose_labels[4])

# Combine all data into one dataframe
all_vaccine_data <- bind_rows(data_cum1, data_cum2, data_cum3, data_cum4)

# --- 6. DEFINE COLORS ---

vaccine_colors <- c(
  "1st Dose" = "skyblue",
  "2nd Dose" = "pink",
  "3rd Dose" = "#107F97",
  "4th Dose" = "#BDCE45"
)

# --- 7. CREATE VACCINATION PLOT ---

create_vaccination_plot <- function(data, colors = vaccine_colors) {
  
  # Define reference dates (only valid dates)
  ref_dates <- c(
    "2020-11-18", "2021-02-26", "2021-07-14",
    "2021-10-09", "2022-03-25", "2022-08-31",
    "2022-11-16", "2023-05-09", "2023-10-25",
    "2024-03-24"
  )
  
  # Create plot
  p <- ggplot(data, aes(x = date, y = percentage, fill = dose, color = dose)) +
    # Area and line for each dose
    geom_area(alpha = 0.7, position = "identity") +
    geom_line(size = 1.2) +
    
    # Colors
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    
    # Scales
    scale_y_continuous(
      limits = c(0, 100),
      labels = percent_format(scale = 1),
      breaks = seq(0, 100, 20)
    ) +
    scale_x_date(
      limits = as.Date(c('2020-03-01', '2024-03-30')),
      breaks = date_breaks("3 months"),
      labels = date_format("%b %y")
    ) +
    
    # Reference lines
    geom_vline(
      xintercept = as.numeric(ymd(ref_dates)),
      linetype = "dashed",
      color = "#82b0d2",
      size = 0.8
    ) +
    
    # Labels and theme
    labs(
      title = "COVID-19 Vaccination Coverage",
      subtitle = "Salvador, Brazil",
      x = "Date",
      y = "Cumulative Vaccination (%)",
      fill = "Dose",
      color = "Dose"
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray50")
    )
  
  return(p)
}

# GENERATE PLOT 

Figure3 <- create_vaccination_plot(all_vaccine_data)


# --- 8. SAVE RESULTS ---

# Create output directory if it doesn't exist
dir.create("output", showWarnings = FALSE)

# Save figures as RDS
saveRDS(Figure1, "output/Figure1_cases_deaths.rds")
saveRDS(Figure2, "output/Figure2_variants.rds")
saveRDS(Figure3, "output/Figure3_vaccination.rds")

# --- 9. DISPLAY FIGURES ---
print(Figure1)
print(Figure2)
print(Figure3)


