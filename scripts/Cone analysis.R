# ============================================================
# Cone Landscape Analysis
#
# The landscape functions and the reference antigenic map were
# adapted from the work of Roessler et al.
# https://github.com/acorg/roessler_netzl_et_al2022
# ============================================================

# Load required libraries
library(Racmacs)
library(tidyverse)
library(ablandscapes) 
library(r3js)
library(htmlwidgets)

set.seed(100)

# Load custom functions
source("./functions/map_longinfo.R")
source("./functions/sams_landscape_functions.R")

# Load secondary map and extract titer information
full_map_p1_adj <- Racmacs::read.acmap(filename = "maps/Roessler map.ace")
map_jp <- read.acmap("maps/pau da lima p23 and 32 L45-49 all.ace")
map_long <- long_map_info(map_jp)

# Base-2 exponential transformation
map_long <- map_long %>%
  mutate(titer = 10 * 2^(as.numeric(titer) * 10))

# Prepare titer data
titerdata <- map_long %>%
  select(titer, ag_name, sr_name, sr_group) %>%
  mutate(sr_group = substr(sr_name, 1, 3))

# Generate titertables by group
titertables <- titerdata %>%
  group_by(sr_group) %>%
  group_map(get_titertable)

# Fit landscape models
lndscp_fits <- lapply(
  titertables,
  function(titertable) {
    ablandscape.fit(
      titers = titertable,
      bandwidth = 1,
      degree = 1,
      method = "cone",
      error.sd = 1,
      acmap = full_map_p1_adj,
      control = list(
        optimise.cone.slope = TRUE
      )
    )
  }
)

# Get group information - CORRECTED
titertables_groups <- titerdata %>%
  group_by(sr_group) %>%
  group_keys()

# Calculate MTs (mean titers)
mt_data <- titerdata %>%
  mutate(titer = as.numeric(titer)) %>%
  group_by(sr_group, ag_name) %>%
  summarise(mt = mean(titer, na.rm = TRUE), .groups = "drop")

# Modified function for base 3D plot
base_plot_data3js_grs <- function(map, lndscp_fits, highlighted_ags, lims, ag_plot_names, 
                                  add_border = TRUE, add_axis = TRUE){
  
  x_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 1])
  y_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 2])
  z_coords <- rep(0.02, length(highlighted_ags))
  ag_point_size <- c(rep(14, length(highlighted_ags))) / 5
  ag_col <- c(agOutline(map)[agNames(map) %in% highlighted_ags])
  ag_fill <- c(agFill(map)[agNames(map) %in% highlighted_ags])
  labels <- c(ag_plot_names[agNames(map) %in% highlighted_ags])
  border_col <- "grey50"
  
  z_lims <- c(0, 10)
  axis_at <- seq(z_lims[1], z_lims[2], 2)
  
  # Setup 3D plot
  data3js <- ablandscapes:::lndscp3d_setup(
    xlim = lims$xlim,
    ylim = lims$ylim,
    zlim = z_lims,
    aspect.z = 0.5,
    options = list(
      lwd.grid = 0.05,
      sidegrid.lwd = 1,
      sidegrid.col = border_col,
      sidegrid.at = list("z" = axis_at),
      zaxt = "lin"
    ),
    show.axis = FALSE
  )
  
  if(add_axis){
    axis_labels <- seq(0, 1, by = 0.2)
    data3js <- r3js::axis3js(
      data3js,
      side = "z",
      at = axis_at,
      labels = axis_labels,
      cornerside = "f",
      size = 20,
      alignment = "right"
    )
  }
  
  # Add basemap
  data3js <- lndscp3d_map(
    data3js = data3js,
    fit = lndscp_fits[[1]],
    xlim = lims$xlim,
    ylim = lims$ylim,
    zlim = c(0, 10),
    show.map.sera = FALSE,
    options = list(
      opacity.basemap = 0.3
    )
  )
  
  # Add antigen points
  data3js <- r3js::points3js(
    data3js,
    x = x_coords,
    y = y_coords,
    z = z_coords,
    size = ag_point_size,
    col = ag_col,
    fill = ag_fill,
    lwd = 0.5,
    opacity = 1,
    highlight = list(col = "red"),
    label = labels,
    toggle = "Basepoints",
    depthWrite = FALSE,
    shape = "circle filled"
  )
  
  if(add_border){
    data3js <- lines3js(data3js, x = c(lims$xlim[1], lims$xlim[1]), 
                        y = c(lims$ylim[1], lims$ylim[2]), z = c(0, 0),
                        lwd = 1.2, col = border_col)
    data3js <- lines3js(data3js, x = c(lims$xlim[2], lims$xlim[2]), 
                        y = c(lims$ylim[1], lims$ylim[2]), z = c(0, 0),
                        lwd = 1.2, col = border_col)
    
    data3js <- lines3js(data3js, x = c(lims$xlim[1], lims$xlim[2]), 
                        y = c(lims$ylim[1], lims$ylim[1]), z = c(0, 0),
                        lwd = 1.2, col = border_col)
    data3js <- lines3js(data3js, x = c(lims$xlim[1], lims$xlim[2]), 
                        y = c(lims$ylim[2], lims$ylim[2]), z = c(0, 0),
                        lwd = 1.2, col = border_col)
    
    data3js <- r3js::box3js(data3js, col = border_col)
  }
  
  return(data3js)
}

# Create base plot
data3js <- base_plot_data3js_grs(full_map_p1_adj, lndscp_fits, 
                                 agNames(full_map_p1_adj), 
                                 lims = list(xlim = c(-2, 6), ylim = c(-3, 2)), 
                                 agNames(full_map_p1_adj), 
                                 add_axis = TRUE)

# Define landscape colors
land_col <- data.frame(
  Color = c("#FF4040", "#1E90FF", "#FFB90F", "#00688B", "#6A5ACD"), 
  row.names = c("L45", "L46", "L47", "L48", "L49")
)

# Generate final landscape figure - USING titertables_groups CORRECTLY
L45_49 <- plot_landscapes_from_list(data3js,
                                    titertables_groups,  # Now has sr_group column
                                    lndscp_fits,
                                    full_map_p1_adj,
                                    mt_data, 
                                    agNames(full_map_p1_adj),
                                    agNames(full_map_p1_adj),
                                    lndscp_colors = land_col,
                                    show_gmts = FALSE)

# Display figure
L45_49

# Extract fitted values for each survey group
fitted_values_L45 <- lndscp_fits[[1]][["fitted.values"]]
fitted_values_L46 <- lndscp_fits[[2]][["fitted.values"]]
fitted_values_L47 <- lndscp_fits[[3]][["fitted.values"]]
fitted_values_L48 <- lndscp_fits[[4]][["fitted.values"]]
fitted_values_L49 <- lndscp_fits[[5]][["fitted.values"]]

# Display figure
L45_49


# Extract fitted values for each survey group
fitted_values_L45 <- lndscp_fits[[1]][["fitted.values"]]
fitted_values_L46 <- lndscp_fits[[2]][["fitted.values"]]
fitted_values_L47 <- lndscp_fits[[3]][["fitted.values"]]
fitted_values_L48 <- lndscp_fits[[4]][["fitted.values"]]
fitted_values_L49 <- lndscp_fits[[5]][["fitted.values"]]


library(tidyverse)

# Extract fitted values from each antibody landscape
df_fitted_auto <- bind_rows(
  lapply(seq_along(lndscp_fits), function(i) {
    
    # Extract fitted values for landscape i
    fit_values <- lndscp_fits[[i]][["fitted.values"]]
    
    # Create a data frame for this survey
    data.frame(
      Survey = paste0("L", 44 + i),  # L45, L46, L47, L48, L49
      Variant_code = names(fit_values),
      Predicted_z_full = as.numeric(fit_values),
      stringsAsFactors = FALSE
    )
  })
)


# ------------------------------------------------------------------------------
# 1. Mapear nombres de variantes al formato estandarizado
# ------------------------------------------------------------------------------

df_fitted_auto <- df_fitted_auto %>%
  mutate(
    Variant = case_when(
      Variant_code == "D614G"       ~ "Ancestral",
      Variant_code == "B.1.1.7"     ~ "Alpha",
      Variant_code == "B.1.351"     ~ "Beta",
      Variant_code == "B.1.617.2"   ~ "Delta",
      Variant_code == "P.1.1"       ~ "Gamma",
      Variant_code == "BA.1"        ~ "BA.1",
      Variant_code == "BA.2"        ~ "BA.2",
      Variant_code == "BA.5.2.1"    ~ "BA.5",
      Variant_code == "BF.7"        ~ "BF.7",
      Variant_code == "BQ.1.1"      ~ "BQ.1.1",
      Variant_code == "BQ.1.3"      ~ "BQ.1",     
      Variant_code == "XBB.1"       ~ "XBB.1",
      TRUE ~ Variant_code
    )
  )

# ------------------------------------------------------------------------------
# 3. Ordenar y añadir contadores para el formato exacto
# ------------------------------------------------------------------------------

# Definir orden de variantes
orden_variantes <- c("Ancestral", "Alpha", "Beta", "Delta", "Gamma", 
                     "BA.1", "BA.2", "BA.5", "BF.7", "BQ.1.1", "BQ.1", "XBB.1")

df_fitted_auto <- df_fitted_auto %>%
  mutate(Variant = factor(Variant, levels = orden_variantes)) %>%
  arrange(Variant, Survey)

# Añadir contadores (Ancestral...1, Ancestral...2, etc.)
df_final <- df_fitted_auto %>%
  group_by(Variant) %>%
  mutate(
    contador = row_number(),
    variant_label = paste0(Variant, "...", contador)
  ) %>%
  ungroup() %>%
  select(variant_label, Survey, Variant, Predicted_z_full, Variant_code)

# Dividir por 10 como hacías
df_final$Predicted_z_full <- df_final$Predicted_z_full / 10

# Guardar como RDS
saveRDS(df_final, "dataset/cones_pre.rds")


#cones_pre <- readRDS(file = "cones_pre.rds")

cones_pre <- df_final




# ------------------------------------------------------------------------------
# Leave-One-Variant-Out Prediction
# Generate predicted values by excluding each variant from the fit,
# producing a complete table across all variants and antibody landscapes.
# --------------------------------------------------------------------------


# Obtener la lista AUTOMÁTICA de todas las variantes disponibles
lista_variantes <- unique(map_long$ag_name)
print("=== VARIANTES ENCONTRADAS ===")
print(lista_variantes)

# ------------------------------------------------------------------------------
# 2. Función para ajustar modelo con una variante excluida
# ------------------------------------------------------------------------------

ajustar_modelo_sin_variante <- function(variante_na, map_long_data, acmap_object) {
  
  # Preparar datos: convertir la variante seleccionada a NA
  titerdata <- map_long_data %>%
    select(titer, ag_name, sr_name, sr_group) %>%
    mutate(titer = ifelse(ag_name == variante_na, NA, as.numeric(titer))) %>%
    mutate(sr_group = substr(sr_name, 1, 3))
  
  # Crear tablas de títulos por grupo
  titerdata_grouped <- titerdata %>% group_by(sr_group)
  titertables <- titerdata_grouped %>% group_map(get_titertable)
  
  # Ajustar modelos para cada superficie (L45 a L49)
  lndscp_fits <- lapply(titertables, function(titertable) {
    ablandscape.fit(
      titers = titertable,
      bandwidth = 1,
      degree = 1,
      method = "cone",
      error.sd = 1,
      acmap = acmap_object,
      control = list(optimise.cone.slope = TRUE)
    )
  })
  
  # Extraer valores predichos para la variante eliminada
  valores_predichos <- sapply(lndscp_fits, function(fit) {
    valor <- fit[["fitted.values"]][[variante_na]]
    if(is.null(valor)) return(NA) else return(valor)
  })
  
  # Extraer R² de cada modelo
  r_squared <- sapply(lndscp_fits, function(fit) {
    if(!is.null(fit$r.squared)) fit$r.squared else NA
  })
  
  # Crear dataframe con resultados
  resultados <- data.frame(
    variante_excluida = variante_na,
    superficie = c("L45", "L46", "L47", "L48", "L49"),
    valor_predicho = valores_predichos,
    r_squared = r_squared
  )
  
  return(resultados)
}

# ------------------------------------------------------------------------------
# 3. EJECUTAR AUTOMÁTICAMENTE PARA TODAS LAS VARIANTES
# ------------------------------------------------------------------------------

cat("\n=== INICIANDO PROCESAMIENTO DE", length(lista_variantes), "VARIANTES ===\n")

# Lista para almacenar todos los resultados
resultados_totales <- list()

# Barra de progreso opcional
for (i in seq_along(lista_variantes)) {
  variante <- lista_variantes[i]
  cat(sprintf("Procesando variante %d/%d: %s\n", i, length(lista_variantes), variante))
  
  # Ejecutar el modelo
  resultado <- ajustar_modelo_sin_variante(variante, map_long, full_map_p1_adj)
  resultados_totales[[variante]] <- resultado
}

#saveRDS(resultados_totales,file = "dataset/2026 06 25 cones.rds")
# ------------------------------------------------------------------------------
# 4. COMBINAR Y FORMATEAR RESULTADOS
# ------------------------------------------------------------------------------

# Unir todos los resultados en una sola tabla
tabla_final <- bind_rows(resultados_totales)

# Reordenar columnas para mejor visualización
tabla_final <- tabla_final %>%
  select(variante_excluida, superficie, valor_predicho, r_squared) %>%
  arrange(variante_excluida, superficie)

# ------------------------------------------------------------------------------
# 5. CREAR TABLA EN FORMATO PIVOTE (más fácil de leer)
# ------------------------------------------------------------------------------

# Formato ancho: variantes como filas, superficies como columnas
tabla_pivote <- tabla_final %>%
  select(variante_excluida, superficie, valor_predicho) %>%
  pivot_wider(
    names_from = superficie,
    values_from = valor_predicho,
    names_prefix = "predicho_"
  )

# Añadir columna con el promedio de todas las superficies
tabla_pivote$promedio_superficies <- rowMeans(
  tabla_pivote[, c("predicho_L45", "predicho_L46", "predicho_L47", "predicho_L48", "predicho_L49")], 
  na.rm = TRUE
)

# ------------------------------------------------------------------------------
# Reshape to long format
# ------------------------------------------------------------------------------

long_table <- tabla_pivote %>%
  pivot_longer(
    cols = starts_with("predicho_"),
    names_to = "survey",
    values_to = "Predicted_z_excl_cone"
  ) %>%
  mutate(
    # Clean survey names (remove "predicho_" prefix)
    survey = str_remove(survey, "predicho_"),
    
    # Standardize variant names
    Variant = case_when(
      variante_excluida == "B.1.1.7"    ~ "Alpha",
      variante_excluida == "B.1.351"    ~ "Beta",
      variante_excluida == "B.1.617.2"  ~ "Delta",
      variante_excluida == "BA.1"       ~ "BA.1",
      variante_excluida == "BA.2"       ~ "BA.2",
      variante_excluida == "BA.5.2.1"   ~ "BA.5",
      variante_excluida == "BF.7"       ~ "BF.7",
      variante_excluida == "BQ.1.1"     ~ "BQ.1.1",
      variante_excluida == "BQ.1.3"     ~ "BQ.1",
      variante_excluida == "D614G"      ~ "Ancestral",
      variante_excluida == "P.1.1"      ~ "Gamma",
      variante_excluida == "XBB.1"      ~ "XBB.1",
      TRUE ~ variante_excluida
    )
  ) %>%
  arrange(Variant, survey) %>%
  select(variante_excluida, survey, Variant, Predicted_z_excl_cone)

# Convert predictions to antigenic units
long_table$Predicted_z_excl_cone <- long_table$Predicted_z_excl_cone / 10

# Define variant order
variant_order <- c(
  "Ancestral", "Alpha", "Beta", "Delta", "Gamma",
  "BA.1", "BA.2", "BA.5", "BF.7", "BQ.1", "BQ.1.1", "XBB.1"
)

final_table <- long_table %>%
  mutate(Variant = factor(Variant, levels = variant_order)) %>%
  arrange(Variant, survey)

# Save leave-one-out cone predictions
cones_one <- final_table
# saveRDS(cones_one, "cones_one.rds")

# ------------------------------------------------------------------------------
# Merge leave-one-out and full landscape predictions
# ------------------------------------------------------------------------------

cones_one_clean <- cones_one %>%
  rename(Survey = survey)

merged_data <- cones_one_clean %>%
  inner_join(
    cones_pre %>%
      select(Survey, Variant, Predicted_z_full, Variant_code),
    by = c("Survey", "Variant")
  )

cat("=== MERGED DATA ===\n")
print(head(merged_data, 10))

# ------------------------------------------------------------------------------
# Extract mean titers from the antigenic map
# ------------------------------------------------------------------------------

# Read antigenic map
map_jp <- read.acmap("maps/pau da lima p23 and 32 L45-49 all.ace")
map_long <- long_map_info(map_jp)

# Extract titer information
titer_data <- map_long %>%
  select(titer, ag_name, sr_name, sr_group) %>%
  mutate(sr_group = substr(sr_name, 1, 3))

# Calculate mean titers
mt_data <- titer_data %>%
  mutate(titer = as.numeric(titer)) %>%
  group_by(sr_group, ag_name) %>%
  summarise(mt = mean(titer, na.rm = TRUE), .groups = "drop")

# Standardize variant names
mt_data_renamed <- mt_data %>%
  rename(Survey = sr_group) %>%
  mutate(
    Variant = case_when(
      ag_name == "B.1.1.7"   ~ "Alpha",
      ag_name == "B.1.351"   ~ "Beta",
      ag_name == "B.1.617.2" ~ "Delta",
      ag_name == "P.1.1"     ~ "Gamma",
      ag_name == "D614G"     ~ "Ancestral",
      ag_name == "BA.1"      ~ "BA.1",
      ag_name == "BA.2"      ~ "BA.2",
      ag_name == "BA.5.2.1"  ~ "BA.5",
      ag_name == "BF.7"      ~ "BF.7",
      ag_name == "BQ.1.1"    ~ "BQ.1.1",
      ag_name == "BQ.1.3"    ~ "BQ.1",
      ag_name == "XBB.1"     ~ "XBB.1",
      TRUE ~ ag_name
    )
  )

# ------------------------------------------------------------------------------
# Merge predictions with observed mean titers
# ------------------------------------------------------------------------------

merged_complete <- merged_data %>%
  left_join(
    mt_data_renamed %>%
      select(Survey, Variant, mt, ag_name),
    by = c("Survey", "Variant")
  )

# ================================================================
# CALCULATE MSE AND 95% CI FOR CONE MODELS
# Full Model and Exclusion Model (LOOOCV)
# ================================================================

# Assuming merged_complete already has the columns:
# - mt: measured titers
# - Predicted_z_full: predictions from full cone model
# - Predicted_z_excl_cone: predictions from exclusion cone model

# ================================================================
# 1. CALCULATE ERRORS
# ================================================================

merged_complete <- merged_complete %>%
  mutate(
    # Full model errors
    absolute_error = abs(mt - Predicted_z_full),
    relative_error = abs((mt - Predicted_z_full) / mt) * 100,
    squared_error = (mt - Predicted_z_full)^2,
    bias = mt - Predicted_z_full,
    
    # Exclusion model errors (LOOOCV)
    absolute_error_LOOOCV = abs(mt - Predicted_z_excl_cone),
    relative_error_LOOOCV = abs((mt - Predicted_z_excl_cone) / mt) * 100,
    squared_error_LOOOCV = (mt - Predicted_z_excl_cone)^2,
    bias_LOOOCV = mt - Predicted_z_excl_cone
  )

# ================================================================
# 2. CALCULATE MSE AND 95% CI BY SURVEY
# ================================================================

# Function to calculate mean and 95% CI
calc_mse_ci <- function(data) {
  n <- length(data)
  mean_val <- mean(data, na.rm = TRUE)
  se <- sd(data, na.rm = TRUE) / sqrt(n)
  ci_lower <- mean_val - 1.96 * se
  ci_upper <- mean_val + 1.96 * se
  return(c(mean = mean_val, lower = ci_lower, upper = ci_upper))
}

# Calculate for all data
mse_all_full <- calc_mse_ci(merged_complete$squared_error)
mse_all_excl <- calc_mse_ci(merged_complete$squared_error_LOOOCV)

# Calculate by survey group
mse_by_survey <- merged_complete %>%
  group_by(Survey) %>%
  summarise(
    # Full model
    mse_full = mean(squared_error, na.rm = TRUE),
    mse_full_lower = mean(squared_error, na.rm = TRUE) - 1.96 * (sd(squared_error, na.rm = TRUE) / sqrt(n())),
    mse_full_upper = mean(squared_error, na.rm = TRUE) + 1.96 * (sd(squared_error, na.rm = TRUE) / sqrt(n())),
    # Exclusion model
    mse_excl = mean(squared_error_LOOOCV, na.rm = TRUE),
    mse_excl_lower = mean(squared_error_LOOOCV, na.rm = TRUE) - 1.96 * (sd(squared_error_LOOOCV, na.rm = TRUE) / sqrt(n())),
    mse_excl_upper = mean(squared_error_LOOOCV, na.rm = TRUE) + 1.96 * (sd(squared_error_LOOOCV, na.rm = TRUE) / sqrt(n())),
    n = n()
  )

# ================================================================
# 3. CREATE THE TABLE
# ================================================================

# Create results table
results_table <- data.frame(
  Survey = c("all", as.character(sort(unique(merged_complete$Survey[!is.na(merged_complete$Survey)])))),
  Cone_Full_MSE = NA,
  Cone_Full_CI = NA,
  Cone_Excl_MSE = NA,
  Cone_Excl_CI = NA
)

# Add "all" data
results_table[1, "Cone_Full_MSE"] <- round(mse_all_full["mean"], 5)
results_table[1, "Cone_Full_CI"] <- sprintf("(%.5f–%.5f)", 
                                            round(mse_all_full["lower"], 5), 
                                            round(mse_all_full["upper"], 5))
results_table[1, "Cone_Excl_MSE"] <- round(mse_all_excl["mean"], 5)
results_table[1, "Cone_Excl_CI"] <- sprintf("(%.5f–%.5f)", 
                                            round(mse_all_excl["lower"], 5), 
                                            round(mse_all_excl["upper"], 5))

# Add survey data
for(i in 1:nrow(mse_by_survey)) {
  survey_num <- as.character(mse_by_survey$Survey[i])
  row_idx <- which(results_table$Survey == survey_num)
  
  results_table[row_idx, "Cone_Full_MSE"] <- round(mse_by_survey$mse_full[i], 5)
  results_table[row_idx, "Cone_Full_CI"] <- sprintf("(%.5f–%.5f)", 
                                                    round(mse_by_survey$mse_full_lower[i], 5), 
                                                    round(mse_by_survey$mse_full_upper[i], 5))
  results_table[row_idx, "Cone_Excl_MSE"] <- round(mse_by_survey$mse_excl[i], 5)
  results_table[row_idx, "Cone_Excl_CI"] <- sprintf("(%.5f–%.5f)", 
                                                    round(mse_by_survey$mse_excl_lower[i], 5), 
                                                    round(mse_by_survey$mse_excl_upper[i], 5))
}

# ================================================================
# 4. PRINT THE TABLE
# ================================================================

print(results_table)

# Create output directory if it doesn't exist
if(!dir.exists("output")) {
  dir.create("output")
}

# Save as CSV
write.csv(results_table, "output/mse_results_table.csv", row.names = FALSE)
