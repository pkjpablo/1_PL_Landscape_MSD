# ============================================================
# SARS-CoV-2 Landscape GAM - Production Script
# Author: Pablo
# Version: 1
# ============================================================

# --- 1. LOAD LIBRARIES ---
library(dplyr)
library(mgcv)
library(plotly)
library(stringr)

# --- 2. LOAD PROCESSED DATA ---
db <- readRDS("dataset/192 long.rds")

db <- db %>%
  mutate(titer = if_else(titer < 0, 0, titer))
db <- db[, c("titer", "sr_name", "ag_name")]

if (!file.exists("dataset/db_long_processed.rds")) {
  saveRDS(db, "dataset/db_long_processed.rds")
}

df_coords <- readRDS("dataset/df_coords_processed.rds") ## Rössler cartography 
db_long <- readRDS("dataset/db_long_processed.rds")

head(db_long,n = 10)

message("✓ Datos cargados")
message("  - df_coords: ", nrow(df_coords), " rows")
message("  - db_long: ", nrow(db_long), " rows")

# --- 3. COLORS ---
mapColors <- c(
  "Ancestral" = "#393b79",
  "Alpha"     = "#637939",
  "Beta"      = "#CD9B1D",
  "Gamma"     = "#b471ab",
  "Delta"     = "#d18652",
  "BA.1"      = "#EF3737",
  "BA.2"      = "#5B004C",
  "BA.5"      = "#bfdaa0",
  "BF.7"      = "#647A39",
  "BQ.1"      = "#107f97",
  "BQ.1.1"    = "#107f97",
  "XBB.1"     = "#5B004C"
)

# --- 4. FUNCTION: CREATE FIGURE ---

make_survey_fig <- function(
    survey_name, db_long, df_coords,
    z_range = c(0, 100), text_size = 14,
    expand_xy = 0.5, show_labels = TRUE,
    show_ci = TRUE
) {
  # Filter data for specific survey
  db_sub <- db_long %>% filter(stringr::str_starts(sr_name, survey_name))
  
  # Merge with coordinates
  merged <- merge(db_sub, df_coords, by = "ag_name", all.x = TRUE)
  
  # Prepare data
  data <- data.frame(
    x = merged$x,
    y = merged$y,
    z = as.numeric(merged$titer),
    variant = as.character(merged$ag_name)
  )
  
  # Aggregate data
  data2 <- data %>%
    group_by(x, y, variant) %>%
    summarise(
      z_mean = mean(z, na.rm = TRUE),
      z_sd = sd(z, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(
      se = z_sd / sqrt(n),
      z_low = z_mean - 1.96 * se,
      z_up = z_mean + 1.96 * se,
      hover_label = sprintf(
        "%s: %.1f (%.1f - %.1f)",
        variant, z_mean, z_low, z_up
      )
    )
  
  # GAM fitting
  gam_fit <- mgcv::gam(z ~ s(x, y, k = 12), data = data)
  
  # Create prediction grid
  x_seq <- seq(min(data$x), max(data$x), length.out = 100)
  y_seq <- seq(min(data$y), max(data$y), length.out = 100)
  grid <- expand.grid(x = x_seq, y = y_seq)
  preds <- predict(gam_fit, newdata = grid, se.fit = TRUE)
  
  grid$z_hat <- preds$fit
  grid$z_low <- preds$fit - 1.96 * preds$se.fit
  grid$z_up <- preds$fit + 1.96 * preds$se.fit
  
  # Create matrices for surface plot
  z_matrix <- matrix(grid$z_hat, nrow = length(x_seq), ncol = length(y_seq))
  z_matrix <- t(z_matrix)
  
  z_low_matrix <- matrix(grid$z_low, nrow = length(x_seq), ncol = length(y_seq))
  z_low_matrix <- t(z_low_matrix)
  
  z_up_matrix <- matrix(grid$z_up, nrow = length(x_seq), ncol = length(y_seq))
  z_up_matrix <- t(z_up_matrix)
  
  # Expand axes
  x_min_exp <- min(data$x) - expand_xy
  x_max_exp <- max(data$x) + expand_xy
  y_min_exp <- min(data$y) - expand_xy
  y_max_exp <- max(data$y) + expand_xy
  
  # Create plot
  fig <- plotly::plot_ly()
  
  # Main surface
  fig <- fig %>%
    plotly::add_surface(
      x = x_seq, y = y_seq, z = z_matrix,
      opacity = 0.6, showscale = FALSE, colorscale = "Viridis"
    )
  
  # Confidence intervals (optional)
  if (show_ci) {
    fig <- fig %>%
      plotly::add_surface(
        x = x_seq, y = y_seq, z = z_low_matrix,
        opacity = 0.25, showscale = FALSE, colors = "#999999", name = "CI lower"
      ) %>%
      plotly::add_surface(
        x = x_seq, y = y_seq, z = z_up_matrix,
        opacity = 0.25, showscale = FALSE, colors = "#999999", name = "CI upper"
      )
  }
  
  # Base points (z = 0)
  fig <- fig %>%
    plotly::add_markers(
      data = data2,
      x = ~x, y = ~y, z = 0,
      color = ~variant,
      colors = mapColors,
      marker = list(size = 8, opacity = 0.3, symbol = "circle"),
      hoverinfo = "text",
      text = ~paste("Base -", variant),
      name = "Base (z = 0)"
    )
  
  # 3D points with mean and CI in tooltip
  fig <- fig %>%
    plotly::add_markers(
      data = data2,
      x = ~x, y = ~y, z = ~z_mean,
      color = ~variant,
      colors = mapColors,
      marker = list(size = 6, opacity = 0.9),
      hoverinfo = "text",
      text = ~hover_label,
      name = "3D Points"
    )
  
  # Labels (optional)
  if (show_labels) {
    fig <- fig %>%
      plotly::add_text(
        data = data2,
        x = ~x, y = ~y, z = ~z_mean,
        text = ~variant,
        textposition = "top center",
        textfont = list(size = text_size, color = "black"),
        showlegend = FALSE
      )
  }
  
  # Lines from base plane
  for (i in seq_len(nrow(data2))) {
    fig <- fig %>%
      plotly::add_trace(
        x = c(data2$x[i], data2$x[i]),
        y = c(data2$y[i], data2$y[i]),
        z = c(0, data2$z_mean[i]),
        type = "scatter3d",
        mode = "lines",
        line = list(color = "black", width = 2),
        showlegend = FALSE
      )
  }
  
  # Layout
  fig <- fig %>%
    plotly::layout(
      title = "",
      scene = list(
        xaxis = list(
          title = "", 
          tickfont = list(size = text_size), 
          range = c(x_min_exp, x_max_exp)
        ),
        yaxis = list(
          title = "", 
          tickfont = list(size = text_size), 
          range = c(y_min_exp, y_max_exp)
        ),
        zaxis = list(
          title = "", 
          range = c(0, 1), 
          tickformat = ".0%", 
          tickfont = list(size = text_size)
        ),
        camera = list(eye = list(x = 1.5, y = 1.5, z = 0.7))
      ),
      showlegend = FALSE
    )
  
  return(fig)
}

# --- 5. ADJUST FIGURE VIEW ---

adjust_3d_figure <- function(
    figure, 
    zoom = 1.7, 
    x_eye = 0.1, 
    y_eye = -1.25, 
    z_eye = 0.2, 
    text_size = 14
) {
  figure %>%
    plotly::layout(
      scene = list(
        zaxis = list(
          title = "",
          range = c(0, 1),
          tickmode = "array",
          tickvals = seq(0, 1, 0.2),
          ticktext = paste0(seq(0, 100, 20), "%"),
          tickfont = list(size = text_size)
        ),
        camera = list(
          eye = list(
            x = x_eye * zoom,
            y = y_eye * zoom,
            z = z_eye * zoom
          )
        ),
        xaxis = list(tickfont = list(size = text_size)),
        yaxis = list(tickfont = list(size = text_size))
      )
    )
}

# --- 6. GENERATE ALL FIGURES ---

message("\n=== GENERATING FIGURES ===\n")

# List of surveys
surveys <- c("L45", "L46", "L47", "L48", "L49")

# Generate and adjust figures
figures <- list()
for (s in surveys) {
  message("Creating figure for ", s, "...")
  fig <- make_survey_fig(
    s, db_long, df_coords, 
    show_labels = TRUE, 
    show_ci = FALSE
  )
  figures[[s]] <- adjust_3d_figure(fig)
}

# Assign to variables
fig_L45 <- figures[["L45"]]
fig_L46 <- figures[["L46"]]
fig_L47 <- figures[["L47"]]
fig_L48 <- figures[["L48"]]
fig_L49 <- figures[["L49"]]


# ============================================================
# HIGHEST AND LOWEST POINTS OF EACH GAM SURFACE
# ============================================================

# --- Function to find highest and lowest points of a GAM surface ---
find_surface_extremes <- function(survey_name, db_long, df_coords, grid_size = 100) {
  
  # Extract data for this survey
  db_sub <- db_long %>% filter(stringr::str_starts(sr_name, survey_name))
  merged <- merge(db_sub, df_coords, by = "ag_name", all.x = TRUE)
  
  data <- data.frame(
    x = merged$x,
    y = merged$y,
    z = as.numeric(merged$titer),
    variant = merged$ag_name
  )
  
  # Fit GAM
  gam_fit <- mgcv::gam(z ~ s(x, y, k = 12), data = data)
  
  # Create prediction grid
  x_seq <- seq(min(data$x), max(data$x), length.out = grid_size)
  y_seq <- seq(min(data$y), max(data$y), length.out = grid_size)
  grid <- expand.grid(x = x_seq, y = y_seq)
  grid$z_hat <- predict(gam_fit, newdata = grid)
  
  # Find highest point (peak)
  max_idx <- which.max(grid$z_hat)
  peak <- grid[max_idx, ]
  
  # Find lowest point (valley)
  min_idx <- which.min(grid$z_hat)
  valley <- grid[min_idx, ]
  
  # Find nearest variants
  dists_peak <- sqrt((data$x - peak$x)^2 + (data$y - peak$y)^2)
  nearest_idx_peak <- which.min(dists_peak)
  nearest_variant_peak <- data$variant[nearest_idx_peak]
  
  dists_valley <- sqrt((data$x - valley$x)^2 + (data$y - valley$y)^2)
  nearest_idx_valley <- which.min(dists_valley)
  nearest_variant_valley <- data$variant[nearest_idx_valley]
  
  # Return result
  data.frame(
    Survey = survey_name,
    # Peak (highest)
    Peak_Height = round(peak$z_hat * 100, 2),
    Peak_Height_z = round(peak$z_hat, 4),
    Peak_x = round(peak$x, 3),
    Peak_y = round(peak$y, 3),
    Peak_Variant = nearest_variant_peak,
    # Valley (lowest)
    Valley_Height = round(valley$z_hat * 100, 2),
    Valley_Height_z = round(valley$z_hat, 4),
    Valley_x = round(valley$x, 3),
    Valley_y = round(valley$y, 3),
    Valley_Variant = nearest_variant_valley,
    # Range (difference)
    Range = round((peak$z_hat - valley$z_hat) * 100, 2)
  )
}

# --- Find extremes for all surveys ---
surveys <- c("L45", "L46", "L47", "L48", "L49")
extremes <- do.call(rbind, lapply(surveys, find_surface_extremes, db_long, df_coords))

# --- Display results ---
cat("\n", "=", rep("-", 80), "\n")
cat("  HIGHEST AND LOWEST POINTS OF EACH GAM SURFACE\n")
cat("=", rep("-", 80), "\n")
print(extremes, row.names = FALSE)
cat("=", rep("-", 80), "\n\n")

# --- Summary table (simplified) ---
summary_extremes <- extremes %>%
  select(Survey, Peak_Height, Peak_Variant, Valley_Height, Valley_Variant, Range)

cat("\n", "=", rep("-", 70), "\n")
cat("  SUMMARY - PEAK vs VALLEY\n")
cat("=", rep("-", 70), "\n")
print(summary_extremes, row.names = FALSE)
cat("=", rep("-", 70), "\n\n")

# --- Save results ---
dir.create("output", showWarnings = FALSE)
write.csv(extremes, "output/gam_extremes.csv", row.names = FALSE)
saveRDS(extremes, "output/gam_extremes.rds")

message("\n✓ Extremes saved to 'output/gam_extremes.csv'")


# --- 7. DISPLAY FIGURES ---

message("\n=== DISPLAYING FIGURES ===\n")
print(fig_L45)
print(fig_L46)
print(fig_L47)
print(fig_L48)
print(fig_L49)

# --- 8. SAVE FIGURES ---

message("\n=== SAVING FIGURES ===\n")

# Create output directory
dir.create("output", showWarnings = FALSE)

# Save as HTML (interactive)
htmlwidgets::saveWidget(fig_L45, "output/fig_L45.html")
htmlwidgets::saveWidget(fig_L46, "output/fig_L46.html")
htmlwidgets::saveWidget(fig_L47, "output/fig_L47.html")
htmlwidgets::saveWidget(fig_L48, "output/fig_L48.html")
htmlwidgets::saveWidget(fig_L49, "output/fig_L49.html")



library(htmltools)

# Crear espacio vacío
empty_div <- tags$div(style = "width:100%; height:400px;")

# Función auxiliar para crear un panel con título y figura
panel_with_title <- function(title, figure) {
  tags$div(
    style = "text-align:center;",
    tags$h3(title, style = "margin-bottom:5px; font-family:sans-serif; font-size:18px;"),
    figure
  )
}

# Crear panel HTML con estilo CSS para 2x3 layout
panel_layout <- tags$div(
  style = "
    display: grid;
    grid-template-columns: repeat(3, 33%);
    grid-template-rows: auto auto;
    gap: 10px;
    padding: 10px;
  ",
  
  # Primera fila (2 gráficos + espacio vacío)
  panel_with_title("Survey 1", fig_L45),
  panel_with_title("Survey 2", fig_L46),
  empty_div,
  
  # Segunda fila (3 gráficos)
  panel_with_title("Survey 3", fig_L47),
  panel_with_title("Survey 4", fig_L48),
  panel_with_title("Survey 5", fig_L49)
)

# Mostrar como página HTML
browsable(panel_layout)

save_html(browsable(panel_layout), file = "Figures 3d/Longitudinal_landscapes.html")


# ============================================================
# CÁLCULO DE VOLUMEN - VERSIÓN COMPACTA
# ============================================================

library(dplyr)
library(mgcv)

# --- Extraer datos ---
get_survey_data <- function(survey_name, db_long, df_coords) {
  db_long %>%
    filter(str_starts(sr_name, survey_name)) %>%
    merge(df_coords, by = "ag_name", all.x = TRUE) %>%
    transmute(
      variant = ag_name,
      x = x,
      y = y,
      z = as.numeric(titer)
    )
}

# --- Calcular volumen ---
calc_volume <- function(data, k = 12, grid_size = 100, boot = FALSE, n_boot = 500) {
  # Grilla
  x_seq <- seq(min(data$x), max(data$x), length.out = grid_size)
  y_seq <- seq(min(data$y), max(data$y), length.out = grid_size)
  dA <- diff(range(x_seq)) / (grid_size - 1) * diff(range(y_seq)) / (grid_size - 1)
  total_area <- diff(range(x_seq)) * diff(range(y_seq))
  
  # Función de volumen
  get_volume <- function(gam_fit) {
    grid <- expand.grid(x = x_seq, y = y_seq)
    grid$z_hat <- predict(gam_fit, newdata = grid)
    (sum(grid$z_hat, na.rm = TRUE) * dA / (total_area * 1)) * 100
  }
  
  # Modelo principal
  gam_fit <- gam(z ~ s(x, y, k = k), data = data)
  volume <- get_volume(gam_fit)
  
  if (!boot) return(volume)
  
  # Bootstrap
  boot_vols <- numeric(n_boot)
  for (i in 1:n_boot) {
    boot_data <- data[sample(1:nrow(data), replace = TRUE), ]
    gam_boot <- tryCatch(gam(z ~ s(x, y, k = k), data = boot_data), 
                         error = function(e) NULL)
    if (!is.null(gam_boot)) boot_vols[i] <- get_volume(gam_boot)
  }
  boot_vols <- boot_vols[!is.na(boot_vols)]
  
  list(
    volume = volume,
    mean = mean(boot_vols),
    lower = quantile(boot_vols, 0.025),
    upper = quantile(boot_vols, 0.975)
  )
}

# --- Ejecutar ---
surveys <- c("L45", "L46", "L47", "L48", "L49")
volumes <- list()

for (s in surveys) {
  data <- get_survey_data(s, db_long, df_coords)
  volumes[[s]] <- calc_volume(data, boot = TRUE, n_boot = 500)
  print(volumes[[s]])
}

# --- Guardar ---
saveRDS(volumes, "output/volumes_all.rds")


volumes



# ============================================================
# FIGURA 3 - VOLUMEN (VERSIÓN PROFESIONAL)
# ============================================================

library(ggplot2)
library(dplyr)
library(lubridate)
library(scales)
library(ggthemes)

# --- Datos ---
volume_data <- data.frame(
  survey = c("L45", "L46", "L47", "L48", "L49"),
  date = as.Date(c("2021-01-07", "2021-08-26", "2022-06-12", "2023-02-11", "2024-01-08")),
  volume_mean = c(8.60, 25.86, 69.40, 76.15, 81.7),
  volume_lower = c(8.07, 24.14, 67.55, 74.50, 80.1),
  volume_upper = c(9.20, 27.73, 71.23, 77.82, 83.3)
)

vertical_lines <- as.Date(c(
  "2019-09-09", "2019-11-09",
  "2020-11-18", "2021-02-26",
  "2021-07-14", "2021-10-09",
  "2022-03-25", "2022-08-31",
  "2022-11-16", "2023-05-09",
  "2023-10-25", "2024-03-24"
))

# --- Paleta de colores profesional ---
colors <- c(
  "Volume" = "#006D77",
  "CI" = "#83C5BE"
)

# --- Crear gráfico ---
Fig3_pro <- ggplot(volume_data, aes(x = date)) +
  
  # Área de confianza
  geom_ribbon(aes(ymin = volume_lower, ymax = volume_upper),
              fill = colors["CI"], alpha = 0.3) +
  
  # Línea principal
  geom_line(aes(y = volume_mean),
            color = colors["Volume"], size = 1.5) +
  
  # Puntos
  geom_point(aes(y = volume_mean),
             color = colors["Volume"], size = 4) +
  
  # Barras de error
  geom_errorbar(aes(ymin = volume_lower, ymax = volume_upper),
                color = colors["Volume"], width = 0.2, alpha = 0.5) +
  
  # Etiquetas de valor
  geom_label(aes(y = volume_mean + 5, 
                 label = paste0(round(volume_mean, 1), "%")),
             fill = "white", color = colors["Volume"],
             size = 4, fontface = "bold", alpha = 0.8) +
  
  # Líneas verticales de referencia
  geom_vline(xintercept = as.numeric(vertical_lines),
             linetype = "dashed", color = "#82b0d2", size = 0.4, alpha = 0.6) +
  
  # Escalas
  scale_x_date(
    date_breaks = "3 months",
    date_labels = "%b %Y",
    limits = as.Date(c('2020-03-01', '2024-06-01'))
  ) +
  
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = percent_format(scale = 1),
    expand = expansion(mult = c(0, 0.1))
  ) +
  
  # Título y etiquetas
  labs(
    title = "Evolution of Antigenic Landscape Volume",
    subtitle = "",
    x = "",
    y = "Volume (%)",
    caption = "Error bars represent 95% CI from bootstrap (n=500)"
  ) +
  
  # Tema
  theme_minimal(base_size = 14) +
  
  theme(
    # Títulos
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, 
                              margin = margin(b = 5)),
    plot.subtitle = element_text(size = 13, color = "gray40", hjust = 0.5,
                                 margin = margin(b = 20)),
    plot.caption = element_text(size = 10, color = "gray50", hjust = 0,
                                margin = margin(t = 10)),
    
    # Ejes
    axis.title.y = element_text(size = 13, face = "bold", 
                                margin = margin(r = 10)),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    
    # Grid
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", size = 0.5),
    
    # Márgenes
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  )

# --- Guardar ---
print(Fig3_pro)

ggsave("output/Fig3_volume_professional.png", Fig3_pro, width = 12, height = 7, dpi = 300)
ggsave("output/Fig3_volume_professional.pdf", Fig3_pro, width = 12, height = 7)

