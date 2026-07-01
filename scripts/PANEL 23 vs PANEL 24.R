# ================================================================
# CORRELATION ANALYSIS: PANEL 23 vs PANEL 24
# Spike Ancestral and Omicron BA.1 Inhibition
# ================================================================

# Load required libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)

# ================================================================
# 1. DATA IMPORT
# ================================================================
load("dataset/correlation_data.RData")

# ================================================================
# 4. CREATE PLOTS
# ================================================================

# Plot A: Spike Ancestral
plot_spike <- ggplot(df_spike, aes(x = spike_p23, y = spike_p24)) +
  geom_point(size = 3, alpha = 0.7, color = "black") +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", 
              color = "black", fill = "grey80", alpha = 0.3) +
  
  # Regression equation
  stat_regline_equation(
    aes(label = paste(..eq.label.., sep = "~~~")),
    label.x.npc = "left",
    label.y.npc = 0.95,
    size = 4.5
  ) +
  
  # Correlation and p-value
  stat_cor(
    method = "spearman",
    aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
    label.x.npc = "left",
    label.y.npc = 0.80,
    size = 4.5,
    expression = TRUE
  ) +
  
  # Axis labels and theme
  labs(
    x = "Panel 23 – % Inhibition Ancestral Spike",
    y = "Panel 24 – % Inhibition Ancestral Spike",
    title = "A) Ancestral Spike"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

# Plot B: Omicron BA.1
plot_ba1 <- ggplot(df_spike_ba1, aes(x = spikeBA1_p23, y = spikeBA1_p24)) +
  geom_point(size = 3, alpha = 0.7, color = "black") +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", 
              color = "black", fill = "grey80", alpha = 0.3) +
  
  # Regression equation
  stat_regline_equation(
    aes(label = paste(..eq.label.., sep = "~~~")),
    label.x.npc = 0.1,
    label.y.npc = 0.95,
    size = 4.5
  ) +
  
  # Correlation and p-value
  stat_cor(
    method = "spearman",
    aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
    label.x.npc = 0.1,
    label.y.npc = 0.80,
    size = 4.5,
    expression = TRUE
  ) +
  
  # Axis labels and theme
  labs(
    x = "Panel 23 – % Inhibition Omicron BA.1",
    y = "Panel 24 – % Inhibition Omicron BA.1",
    title = "B) Omicron BA.1"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

# ================================================================
# 5. COMBINE PLOTS
# ================================================================

combined_plot <- plot_grid(
  plot_spike, plot_ba1,
  ncol = 1,
  align = "h",
  labels = NULL
)

# Add a single title
title <- ggdraw() + 
  draw_label(
    "",
    fontface = "bold",
    size = 18,
    hjust = 0.5
  )

final_plot <- plot_grid(
  title,
  combined_plot,
  ncol = 1,
  rel_heights = c(0.1, 1)
)


print(final_plot)
