# ============================================================
# MAPA ANTIGÉNICO - SARS-CoV-2 Variants
# Author: Pablo
# Version: 1
# ============================================================

# --- 1. LOAD LIBRARIES ---
library(ggplot2)
library(dplyr)

# --- 2. LOAD DATA ---
df_coords_color <- readRDS("dataset/df_coords_colors")

# --- 3. VERIFICAR DATOS ---
message("\n=== DATOS CARGADOS ===\n")
message("Variantes: ", nrow(df_coords_color))
print(head(df_coords_color, 5))

# --- 4. CREAR MAPA ANTIGÉNICO ---

Fig_AntigenicMap <- ggplot(df_coords_color, aes(x = x, y = y)) +
  
  # --- 4.1 Puntos ---
  geom_point(aes(fill = color), 
             shape = 21, 
             size = 6, 
             color = "black", 
             stroke = 1) +
  
  # --- 4.2 Etiquetas ---
  geom_text(aes(label = ag_name), 
            hjust = -0.2, 
            vjust = 0.5, 
            size = 4.5,
            fontface = "bold") +
  
  # --- 4.3 Escalas ---
  scale_fill_identity() +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  
  # --- 4.4 Etiquetas ---
  labs(
    title = "Antigenic Map of SARS-CoV-2 Variants",
    subtitle = "",
    x = "Antigenic Map X",
    y = "Antigenic Map Y",
    caption = ""
  ) +
  
  # --- 4.5 Tema ---
  theme_minimal(base_size = 14) +
  
  theme(
    # Títulos
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 13, color = "gray50", hjust = 0.5),
    plot.caption = element_text(size = 10, color = "gray50", hjust = 1),
    
    # Grid
    panel.grid = element_line(color = "gray85", size = 0.5),
    panel.grid.minor = element_blank(),
    
    # Ejes
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11),
    
    # Márgenes
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  )

# --- 5. MOSTRAR ---
print(Fig_AntigenicMap)



