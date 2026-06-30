# ============================================================
# SARS-CoV-2 GAM SURFACE ANALYSIS - COMPLETE SCRIPT
# Author: Pablo
# Version: 3.0
# ============================================================

library(dplyr)
library(tidyr)
library(mgcv)
library(purrr)
library(ggplot2)
library(plotly)
library(htmlwidgets)
library(stringr)

# ============================================================
# PART 1: FUNCIÓN DE VOLUMEN
# ============================================================

calc_gam_volume_percent <- function(data) {
  # Ajustar modelo GAM
  gam_fit <- gam(z ~ s(x, y, k = 12), data = data)
  
  # Crear grilla de predicción
  x_seq <- seq(min(data$x), max(data$x), length.out = 100)
  y_seq <- seq(min(data$y), max(data$y), length.out = 100)
  grid <- expand.grid(x = x_seq, y = y_seq)
  grid$z_hat <- predict(gam_fit, newdata = grid)
  
  # Calcular delta de área por celda
  dx <- diff(range(x_seq)) / (length(x_seq) - 1)
  dy <- diff(range(y_seq)) / (length(y_seq) - 1)
  dA <- dx * dy
  
  # Volumen real: suma de z_hat * área de cada celda
  V_real <- sum(grid$z_hat, na.rm = TRUE) * dA
  
  # Volumen máximo: área total * altura máxima (z = 1)
  V_max <- (diff(range(x_seq)) * diff(range(y_seq))) * 1
  
  # Porcentaje
  perc <- (V_real / V_max) * 100
  
  return(perc)
}


# ============================================================
# PART 2: FUNCIÓN DE ANÁLISIS POR LOTE
# ============================================================

analisis_por_lote <- function(lote) {
  message("Procesando ", lote, "...")
  
  # --- 2.1 Leer datos ---
  df_coords <- readRDS("dataset/df_coords_processed.rds")
  db_long <- readRDS("dataset/db_long_processed.rds")
  
  # --- 2.2 Preprocesamiento ---
  db_long$sr_name <- substr(db_long$sr_name, 1, 13)
  db_long <- db_long %>% separate(sr_name, into = c("survey", "idnova"), sep = "_")
  db_long1 <- db_long %>% filter(str_starts(survey, lote))
  db_long3 <- merge(db_long1, df_coords, by = "ag_name", all.x = TRUE)
  
  db_long3$z <- as.numeric(db_long3$titer)
  data <- db_long3
  
  # --- 2.3 Ajustar GAM por idnova ---
  fits <- data %>%
    group_by(idnova) %>%
    group_split() %>%
    set_names(unique(data$idnova)) %>%
    map(~ gam(z ~ s(x, y, k = 12), data = .x))
  
  # --- 2.4 Grilla común ---
  x_seq <- seq(min(data$x), max(data$x), length.out = 50)
  y_seq <- seq(min(data$y), max(data$y), length.out = 50)
  grid <- expand.grid(x = x_seq, y = y_seq)
  
  # --- 2.5 Predicciones ---
  predictions <- map(fits, function(mod) {
    grid$z_hat <- predict(mod, newdata = grid)
    grid
  })
  
  # --- 2.6 Función auxiliar ---
  make_z_matrix <- function(df) {
    matrix(df$z_hat, nrow = length(x_seq), ncol = length(y_seq))
  }
  
  # --- 2.7 Calcular volumen por GAM ---
  volumen_por_idnova <- data %>%
    group_by(idnova) %>%
    group_split() %>%
    set_names(map_chr(., ~ unique(.x$idnova))) %>%
    map_dbl(calc_gam_volume_percent)
  
  volumen_df <- tibble(
    lote = lote,
    idnova = names(volumen_por_idnova),
    volumen_percent = volumen_por_idnova
  )
  
  # --- 2.8 Superficie plotly ---
  fig <- plot_ly()
  
  for (i in seq_along(predictions)) {
    # Escalar z_hat a porcentaje (0-100%)
    z_scaled <- predictions[[i]]
    z_scaled$z_hat <- z_scaled$z_hat * 100
    
    z_mat <- make_z_matrix(z_scaled)
    
    fig <- fig %>%
      add_surface(
        x = x_seq,
        y = y_seq,
        z = z_mat,
        surfacecolor = matrix(1, nrow = length(x_seq), ncol = length(y_seq)),
        colorscale = list(c(0, 1), c("rgba(80,80,80,0.3)", "rgba(80,80,80,0.3)")),
        showscale = FALSE,
        opacity = 0.3,
        name = paste("idnova", names(predictions)[i]),
        showlegend = FALSE
      )
  }
  
  fig <- fig %>%
    layout(
      title = paste("GAM Surfaces -", lote),
      scene = list(
        xaxis = list(
          title = "",
          showticklabels = FALSE,   # Sin etiquetas en x
          showgrid = FALSE,          # Sin grilla en x
          zeroline = FALSE           # Sin línea cero en x
        ),
        yaxis = list(
          title = "",
          showticklabels = FALSE,   # Sin etiquetas en y
          showgrid = FALSE,          # Sin grilla en y
          zeroline = FALSE           # Sin línea cero en y
        ),
        zaxis = list(
          title = "Z (%)",
          range = c(0, 100),         # Rango 0-100%
          tickmode = "array",
          tickvals = seq(0, 100, 20),
          ticktext = paste0(seq(0, 100, 20), "%"),
          tickfont = list(size = 12))
      )
    )
  
  # --- 2.9 Retornar resultados ---
  list(
    lote = lote,
    fits = fits,
    original_data = data,
    predictions = predictions,
    volumen_df = volumen_df,
    superficie_plotly = fig
  )
}


# ============================================================
# PART 3: EXTRACCIÓN DE PREDICCIONES POR VARIANTE
# ============================================================

extract_predictions <- function(result) {
  
  message("\n=== Extrayendo predicciones para variantes ===")
  
  batch <- result$lote
  fits <- result$fits
  original_data <- result$original_data
  
  # --- 3.1 Identificar columna de variantes ---
  posibles_nombres <- c("ag_name", "variant", "variante", "Variant", "ag", "name", "lineage")
  col_variante <- NULL
  
  for (nombre in posibles_nombres) {
    if (nombre %in% names(original_data)) {
      col_variante <- nombre
      break
    }
  }
  
  if (is.null(col_variante)) {
    exclude <- c("x", "y", "z", "idnova", "survey", "sr_name", "titer", "data")
    candidates <- names(original_data)[!names(original_data) %in% exclude]
    if (length(candidates) > 0) {
      col_variante <- candidates[1]
    } else {
      stop("No se pudo identificar la columna de variantes")
    }
  }
  
  message("  ✓ Columna de variantes: ", col_variante)
  message("  ✓ Variantes: ", paste(unique(original_data[[col_variante]]), collapse = ", "))
  
  # --- 3.2 Coordenadas por variante ---
  variant_coords <- original_data %>%
    group_by(!!sym(col_variante)) %>%
    summarise(
      x = mean(x, na.rm = TRUE),
      y = mean(y, na.rm = TRUE),
      .groups = "drop"
    )
  
  names(variant_coords)[1] <- "variant_name"
  message("  ✓ Coordenadas calculadas para ", nrow(variant_coords), " variantes")
  
  # --- 3.3 Predecir para cada variante en cada idnova ---
  all_predictions <- list()
  
  for (id in names(fits)) {
    message("    Procesando idnova: ", id)
    
    pred_list <- list()
    
    for (i in 1:nrow(variant_coords)) {
      v <- variant_coords$variant_name[i]
      x <- variant_coords$x[i]
      y <- variant_coords$y[i]
      
      newdata <- data.frame(x = x, y = y)
      pred <- predict(fits[[id]], newdata = newdata, se.fit = TRUE)
      
      pred_list[[v]] <- data.frame(
        batch = batch,
        idnova = id,
        variant = v,
        x = x,
        y = y,
        z_pred = pred$fit,
        z_percent = pred$fit * 100,
        se = pred$se.fit,
        lower_ci = (pred$fit - 1.96 * pred$se.fit) * 100,
        upper_ci = (pred$fit + 1.96 * pred$se.fit) * 100
      )
    }
    
    all_predictions[[id]] <- do.call(rbind, pred_list)
  }
  
  # --- 3.4 Combinar ---
  final_df <- do.call(rbind, all_predictions)
  rownames(final_df) <- NULL
  
  message("  ✓ Predicciones completadas! (", nrow(final_df), " filas)")
  
  return(final_df)
}


# ============================================================
# PART 4: VISUALIZACIÓN DE PREDICCIONES
# ============================================================

plot_predictions <- function(all_predictions, variant_filter = NULL) {
  
  # Filtrar si se especifica una variante
  if (!is.null(variant_filter)) {
    plot_data <- all_predictions %>% filter(variant == variant_filter)
    if (nrow(plot_data) == 0) {
      message("⚠️ Variante '", variant_filter, "' no encontrada")
      return(NULL)
    }
    title_suffix <- paste("-", variant_filter)
  } else {
    plot_data <- all_predictions
    title_suffix <- "- Todas las variantes"
  }
  
  # --- 4.1 Gráfico de líneas ---
  p1 <- ggplot(plot_data, aes(x = idnova, y = z_percent, color = variant, group = variant)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = variant),
                alpha = 0.2, color = NA) +
    facet_wrap(~ batch, scales = "free_x") +
    labs(
      title = paste("GAM Predictions by Variant", title_suffix),
      x = "idnova",
      y = "Predicted Z (%)",
      color = "Variant",
      fill = "Variant"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  print(p1)
  
  # --- 4.2 Heatmap ---
  p2 <- ggplot(plot_data, aes(x = idnova, y = variant, fill = z_percent)) +
    geom_tile() +
    geom_text(aes(label = round(z_percent, 1)), size = 2.5) +
    facet_wrap(~ batch, scales = "free_x") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 50) +
    labs(
      title = paste("GAM Predictions Heatmap", title_suffix),
      x = "idnova",
      y = "Variant",
      fill = "Z (%)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  print(p2)
  
  return(list(line_plot = p1, heatmap = p2))
}


# ============================================================
# PART 5: FORMATO ANCHO Y GUARDADO
# ============================================================

save_predictions <- function(all_predictions, output_dir = "output") {
  
  dir.create(output_dir, showWarnings = FALSE)
  
  # --- 5.1 Formato ancho ---
  predictions_wide <- all_predictions %>%
    select(batch, idnova, variant, z_percent) %>%
    pivot_wider(
      names_from = variant,
      values_from = z_percent
    )
  
  # --- 5.2 Formato completo con IC ---
  predictions_full <- all_predictions %>%
    select(batch, idnova, variant, z_percent, lower_ci, upper_ci, x, y)
  
  # --- 5.3 Guardar CSV ---
  write.csv(predictions_wide, file.path(output_dir, "predictions_by_idnova.csv"), row.names = FALSE)
  write.csv(predictions_full, file.path(output_dir, "predictions_all.csv"), row.names = FALSE)
  
  # --- 5.4 Guardar RDS ---
  saveRDS(predictions_wide, file.path(output_dir, "predictions_by_idnova.rds"))
  saveRDS(predictions_full, file.path(output_dir, "predictions_all.rds"))
  
  message("✓ Archivos guardados en '", output_dir, "/'")
  
  return(list(
    wide = predictions_wide,
    full = predictions_full
  ))
}


# ============================================================
# PART 6: VISUALIZAR SUPERFICIES 3D
# ============================================================

display_surfaces <- function(resultados, lotes) {
  
  for (i in 1:length(resultados)) {
    message("\n=== Mostrando superficie 3D - ", lotes[i], " ===")
    
    # Mostrar en el visor de RStudio
    print(resultados[[i]]$superficie_plotly)
    
    # Guardar como HTML
    htmlwidgets::saveWidget(
      resultados[[i]]$superficie_plotly,
      paste0("output/surface_", lotes[i], ".html"),
      selfcontained = TRUE
    )
    
    message("  ✓ Figura guardada: output/surface_", lotes[i], ".html")
  }
}


# ============================================================
# PART 7: VISUALIZAR VOLÚMENES
# ============================================================

plot_volumes <- function(volumenes_totales) {
  
  # --- 7.1 Volumen por idnova ---
  p1 <- ggplot(volumenes_totales, aes(x = idnova, y = volumen_percent, fill = lote)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = round(volumen_percent, 1)),
              position = position_dodge(0.9),
              vjust = -0.5,
              size = 3) +
    facet_wrap(~ lote, scales = "free_x") +
    scale_fill_brewer(palette = "Set1") +
    labs(
      title = "GAM Volume by idnova",
      x = "idnova",
      y = "Volume (%)",
      fill = "Batch"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  print(p1)
  
  # --- 7.2 Volumen promedio por lote ---
  avg_volume <- volumenes_totales %>%
    group_by(lote) %>%
    summarise(
      mean = mean(volumen_percent, na.rm = TRUE),
      sd = sd(volumen_percent, na.rm = TRUE),
      n = n()
    )
  
  p2 <- ggplot(avg_volume, aes(x = lote, y = mean, fill = lote)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                  width = 0.2) +
    geom_text(aes(label = round(mean, 1)),
              vjust = -0.5,
              size = 5,
              fontface = "bold") +
    scale_fill_brewer(palette = "Set1") +
    labs(
      title = "Average GAM Volume by Batch",
      x = "Batch",
      y = "Average Volume (%)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  print(p2)
  
  return(list(by_idnova = p1, by_batch = p2))
}


# ============================================================
# PART 8: EJECUCIÓN COMPLETA
# ============================================================

run_analysis <- function(lotes = c("L45", "L46", "L47", "L48", "L49")) {
  
  message("\n", "=", rep("=", 60), "\n")
  message("  INICIANDO ANÁLISIS GAM")
  message("=", rep("=", 60), "\n")
  
  # --- 8.1 Ejecutar análisis para todos los lotes ---
  message("\n--- PASO 1: Análisis por lote ---\n")
  resultados <- map(lotes, analisis_por_lote)
  names(resultados) <- lotes
  
  # --- 8.2 Unir volúmenes ---
  message("\n--- PASO 2: Procesando volúmenes ---\n")
  volumenes_totales <- bind_rows(map(resultados, "volumen_df"))
  
  message("\nVolúmenes calculados:")
  print(volumenes_totales)
  
  # --- 8.3 Extraer predicciones ---
  message("\n--- PASO 3: Extrayendo predicciones por variante ---\n")
  predicciones_totales <- list()
  
  for (i in 1:length(resultados)) {
    message("\n=== Procesando ", lotes[i], " ===")
    predicciones_totales[[lotes[i]]] <- extract_predictions(resultados[[i]])
  }
  
  # --- 8.4 Combinar todas las predicciones ---
  all_predictions <- bind_rows(predicciones_totales)
  
  message("\n--- PASO 4: Resumen de predicciones ---\n")
  message("  ✓ Total de predicciones: ", nrow(all_predictions))
  message("  ✓ Variantes encontradas: ", 
          paste(unique(all_predictions$variant), collapse = ", "))
  
  # --- 8.5 Guardar predicciones ---
  message("\n--- PASO 5: Guardando resultados ---\n")
  saved <- save_predictions(all_predictions)
  
  # --- 8.6 Mostrar superficies 3D ---
  message("\n--- PASO 6: Mostrando superficies 3D ---\n")
  dir.create("output", showWarnings = FALSE)
  display_surfaces(resultados, lotes)
  
  # --- 8.7 Mostrar volúmenes ---
  message("\n--- PASO 7: Visualizando volúmenes ---\n")
  volume_plots <- plot_volumes(volumenes_totales)
  
  # Guardar gráficos de volumen
  ggsave("output/volumes_by_idnova.png", volume_plots$by_idnova, width = 12, height = 6, dpi = 300)
  ggsave("output/volumes_by_batch.png", volume_plots$by_batch, width = 8, height = 6, dpi = 300)
  
  # --- 8.8 Mostrar predicciones ---
  message("\n--- PASO 8: Visualizando predicciones ---\n")
  plot_pred <- plot_predictions(all_predictions)
  
  if (!is.null(plot_pred)) {
    ggsave("output/predictions_lines.png", plot_pred$line_plot, width = 14, height = 8, dpi = 300)
    ggsave("output/predictions_heatmap.png", plot_pred$heatmap, width = 14, height = 8, dpi = 300)
  }
  
  # --- 8.9 Mostrar Ancestral específico ---
  if ("Ancestral" %in% unique(all_predictions$variant)) {
    message("\n--- PASO 9: Predicciones para Ancestral ---\n")
    ancestral <- all_predictions %>% filter(variant == "Ancestral")
    print(ancestral)
    
    # Guardar Ancestral
    write.csv(ancestral, "output/ancestral_predictions.csv", row.names = FALSE)
    
    # Gráfico de Ancestral
    p_ancestral <- ggplot(ancestral, aes(x = idnova, y = z_percent, fill = batch)) +
      geom_bar(stat = "identity", position = "dodge") +
      geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
                    position = position_dodge(0.9), width = 0.2) +
      geom_text(aes(label = round(z_percent, 1)),
                position = position_dodge(0.9),
                vjust = -0.5, size = 3) +
      facet_wrap(~ batch, scales = "free_x") +
      labs(
        title = "Ancestral Predictions by idnova",
        x = "idnova",
        y = "Predicted Z (%)"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
    
    print(p_ancestral)
    ggsave("output/ancestral_predictions.png", p_ancestral, width = 12, height = 6, dpi = 300)
  }
  
  # --- 8.10 Resumen final ---
  message("\n", "=", rep("=", 60), "\n")
  message("  ANÁLISIS COMPLETADO!")
  message("=", rep("=", 60), "\n")
  message("\nArchivos generados en 'output/':")
  message("  - predictions_by_idnova.csv")
  message("  - predictions_all.csv")
  message("  - surface_*.html (5 archivos)")
  message("  - volumes_by_idnova.png")
  message("  - volumes_by_batch.png")
  message("  - predictions_lines.png")
  message("  - predictions_heatmap.png")
  if ("Ancestral" %in% unique(all_predictions$variant)) {
    message("  - ancestral_predictions.csv")
    message("  - ancestral_predictions.png")
  }
  
  return(list(
    resultados = resultados,
    volumenes = volumenes_totales,
    predicciones = all_predictions,
    predicciones_wide = saved$wide
  ))
}


# ============================================================
# PART 9: EJECUTAR
# ============================================================

# --- Ejecutar el análisis completo ---
resultado_final <- run_analysis()

# --- Ver resumen final ---
message("\nRESUMEN FINAL:")
message("  - Lotes procesados: ", paste(names(resultado_final$resultados), collapse = ", "))
message("  - Total de predicciones: ", nrow(resultado_final$predicciones))
message("  - Variantes: ", paste(unique(resultado_final$predicciones$variant), collapse = ", "))

# ============================================================
# VIOLIN PLOT: ANCESTRAL vs BA.1 - PREDICTED z VALUES
# ============================================================

library(ggplot2)
library(dplyr)

# --- 1. Filtrar datos de Ancestral y BA.1 ---
variantes_interes <- c("Ancestral", "BA.1")

datos_violin <- resultado_final$predicciones %>%
  filter(variant %in% variantes_interes)

# Verificar datos
message("Variantes encontradas:")
print(unique(datos_violin$variant))

message("\nNúmero de observaciones por variante:")
print(table(datos_violin$variant))

# --- 2. Violin Plot básico ---
violin_plot <- ggplot(datos_violin, aes(x = variant, y = z_percent, fill = variant)) +
  
  # Violin
  geom_violin(trim = FALSE, alpha = 0.7, scale = "width") +
  
  # Boxplot interno
  geom_boxplot(width = 0.15, alpha = 0.8, outlier.shape = NA, color = "black") +
  
  # Puntos individuales (jitter)
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.3) +
  
  # Media (punto)
  stat_summary(fun = mean, geom = "point", size = 4, shape = 23, 
               color = "white", fill = "black") +
  
  # Media (label)
  stat_summary(fun = mean, geom = "text", 
               aes(label = round(..y.., 1)),
               vjust = -1.5, size = 4, fontface = "bold") +
  
  # Colores
  scale_fill_manual(values = c("Ancestral" = "#2C7A7B", "BA.1" = "#E74C3C")) +
  
  # Etiquetas
  labs(
    title = "GAM Predicted z Values: Ancestral vs BA.1",
    subtitle = "Distribution across all batches",
    x = "Variant",
    y = "Predicted z (%)",
    fill = "Variant",
    caption = "Diamond = Mean | Box = IQR | Violin = Density"
  ) +
  
  # Tema
  theme_minimal(base_size = 14) +
  
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, color = "gray50", hjust = 0.5),
    plot.caption = element_text(size = 10, color = "gray50", hjust = 1),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

# --- 3. Mostrar ---
print(violin_plot)

# ============================================================
# VIOLIN PLOT POR BATCH - ANCESTRAL vs BA.1
# ============================================================

violin_clean <- ggplot(datos_filtrados, aes(x = variant, y = z_percent, fill = variant)) +
  
  geom_violin(trim = FALSE, alpha = 0.5, scale = "width") +
  geom_boxplot(width = 0.15, alpha = 0.7, outlier.shape = NA, color = "black") +
  
  facet_wrap(~ batch, scales = "free_y", ncol = 5) +
  
  # --- ESCALA Y DE 0 A 100 DE 20 EN 20 ---
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = paste0(seq(0, 100, 20), "%")
  ) +
  
  scale_fill_manual(values = c("Ancestral" = "darkblue", "BA.1" = "skyblue")) +
  
  labs(
    title = "",
    x = "",
    y = "(%)",
    fill = "Variant"
  ) +
  
  theme_minimal(base_size = 12) +
  
  theme(
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "bottom",
    strip.text = element_text(size = 11, face = "bold"),
    panel.grid.major.x = element_blank()
  )

print(violin_clean)




# ============================================================
# SAVE VOLUME DATA BY SURVEY
# ============================================================


library(dplyr)

# --- 1. Load the complete volume data ---
volumenes_df <- as.data.frame(volumenes_totales)

# --- 2. Filter and save each survey ---

# L45
dfL45.2 <- volumenes_df %>%
  filter(lote == "L45") %>%
  select(idnova, volumen_percent) %>%
  rename(id = idnova)


saveRDS(dfL45.2, "dataset/dfL45.2.rds")
message("✓ Saved dfL45.2.rds (", nrow(dfL45.2), " rows)")

# L46
dfL46.2 <- volumenes_df %>%
  filter(lote == "L46") %>%
  select(idnova, volumen_percent) %>%
  rename(id = idnova)

saveRDS(dfL46.2, "dataset/dfL46.2.rds")
message("✓ Saved dfL46.2.rds (", nrow(dfL46.2), " rows)")

# L47
dfL47.2 <- volumenes_df %>%
  filter(lote == "L47") %>%
  select(idnova, volumen_percent) %>%
  rename(id = idnova)

saveRDS(dfL47.2, "dataset/dfL47.2.rds")
message("✓ Saved dfL47.2.rds (", nrow(dfL47.2), " rows)")

# L48
dfL48.2 <- volumenes_df %>%
  filter(lote == "L48") %>%
  select(idnova, volumen_percent) %>%
  rename(id = idnova)

saveRDS(dfL48.2, "dataset/dfL48.2.rds")
message("✓ Saved dfL48.2.rds (", nrow(dfL48.2), " rows)")

# L49
dfL49.2 <- volumenes_df %>%
  filter(lote == "L49") %>%
  select(idnova, volumen_percent) %>%
  rename(id = idnova)

saveRDS(dfL49.2, "dataset/dfL49.2.rds")
message("✓ Saved dfL49.2.rds (", nrow(dfL49.2), " rows)")

