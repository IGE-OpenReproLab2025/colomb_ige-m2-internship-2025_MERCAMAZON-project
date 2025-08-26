EDGAR_Map_Dep_05_08_2025
================
Martin Colomb
2025-08-26

# Chargement des packages

``` r
library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyverse)
library(ncdf4)
library(stringr)
library(tibble)
library(writexl)
```

``` r
chem_hist_B100<-"C:/Users/colom/Desktop/STAGE/data/clean_mod_data/14_3_1/HIST_B100"
chem_hist<-"C:/Users/colom/Desktop/STAGE/data/clean_mod_data/14_3_1"
chem_EDGAR <- "C:/Users/colom/Desktop/STAGE/data/clean_mod_data/14_3_1/EDGAR"
con_hg_B100<-nc_open(file.path(chem_hist_B100, "GEOSChem.SpeciesConc.2015_m.nc4"))
statemet_B100<-nc_open(file.path(chem_hist_B100, "GEOSChem.StateMet.2015_m.nc4"))
area<-ncvar_get(con_hg_B100, "AREA")

Met_AIRDEN<-ncvar_get(statemet_B100, "Met_AIRDEN")  #Dry air density units: kg m-3
range(Met_AIRDEN)
```

    ## [1] 4.299234e-05 1.493199e+00

``` r
dim(Met_AIRDEN)
```

    ## [1] 144  91  47  12

``` r
Met_AIRDEN <- Met_AIRDEN[, , 1, ]
dim(Met_AIRDEN)
```

    ## [1] 144  91  12

``` r
time_raw <- ncvar_get(statemet_B100, "time")
origin_time <- as.POSIXct("2015-01-01 00:00:00", tz = "UTC")
time <- origin_time + time_raw * 60
days_in_months <- days_in_month(time)
weights <- days_in_months / sum(days_in_months)
weighted_avg <- function(x) sum(x * weights)
mean_Met_AIRDEN <- apply(Met_AIRDEN, MARGIN = c(1, 2), FUN = weighted_avg)
max(mean_Met_AIRDEN) 
```

    ## [1] 1.412798

``` r
nc_mask<-nc_open(file.path("C:/Users/colom/Desktop/STAGE/data/clean_mod_data", "Amazon_basin_mask_2x25.nc"))
Am_mask<-ncvar_get(nc_mask, "MASK")



s_in_yr = 365.2425 * 24 * 3600
unit_conv = s_in_yr
kg_Mg = 1e-3
MW_Hg = 200.59
avo = 6.02e23
g_kg = 1e-3
cm2_m2 = 1e4
unit_conv_dd = MW_Hg / avo * g_kg * cm2_m2 * s_in_yr * area
M_Hg <- 0.20059         # Masse molaire du Hg (kg/mol)
M_air <- 0.028964       # Masse molaire de l'air sec (kg/mol)
ng_per_kg <- 1e12       # Conversion kg → ng
```

``` r
# # Paramètres généraux
# simulations <- c("V03A03", "V20A03")  #       , "V20A20", "V03A20"
# lon <- seq(-180, 180, length.out = 144)
# lat <- seq(-90, 90, length.out = 91)
# df_grid <- expand.grid(x = lon, y = lat)
# 
# # Conversion
# kg_Mg <- 1e-3  # kg → Mg
# unit_conv_dd <- M_Hg / M_air * avo * 1e-9 / s_in_yr / cm2_m2  # à définir avant si pas déjà fait
# 
# # Breaks et couleurs pour dépôts
# breaks_dep <- c(0.01, 0.05, 0.1, 0.3, 0.5, 0.8, 1, 1.5, 2, Inf)
# labels_dep_chr <- c(
#   "0.01–0.05", "0.05–0.1", "0.1–0.3", "0.3–0.5", "0.5–0.8",
#   "0.8–1.0", "1.0–1.5", "1.5–2.0", ">2.0"
# )
# labels_dep_expr <- c(
#   expression("0.01–0.05"), expression("0.05–0.1"), expression("0.1–0.3"),
#   expression("0.3–0.5"), expression("0.5–0.8"), expression("0.8–1.0"),
#   expression("1.0–1.5"), expression("1.5–2.0"), expression("">"~2.0")
# )
# colors_dep <- c(
#   "#005F73", "#0A9396", "#94D2BD", "#E9D8A6",
#   "#EE9B00", "#CA6702", "#BB3E03", "#AE2012", "#9B2226"
# )
# 
# # Fonction de plotting
# plot_drydep_map <- function(matrix_dep, sim_name, species, save_path) {
#   df <- df_grid
#   df$dep_Mg <- as.vector(matrix_dep)
#   df$dep_class <- cut(df$dep_Mg,
#                       breaks = breaks_dep,
#                       labels = labels_dep_chr,
#                       include.lowest = TRUE)
#   df <- na.omit(df)
#   
#   plot <- ggplot(df, aes(x = x, y = y, fill = dep_class)) +
#     geom_raster() +
#     scale_fill_manual(
#       name = paste0("Dry Dep ", species, " [Mg]"),
#       values = colors_dep,
#       labels = labels_dep_expr,
#       na.value = "transparent",
#       guide = guide_legend(reverse = TRUE)
#     ) +
#     coord_equal() +
#     theme_minimal() +
#     labs(
#       title = paste("Annual Dry Deposition of", species, "in 2015"),
#       subtitle = paste0("Simulation : ", sim_name),
#       x = "Longitude", y = "Latitude"
#     ) +
#     borders("world", colour = "#403d39", size = 0.2) +
#     theme(
#       plot.title = element_text(hjust = 0.5, face = "bold"),
#       plot.subtitle = element_text(hjust = 0.5),
#       axis.text = element_text(size = 10),
#       axis.ticks = element_line(color = "black", size = 0.3),
#       axis.ticks.length = unit(3, "pt"),
#       panel.grid = element_blank(),
#       panel.border = element_rect(color = "black", fill = NA, size = 0.8),
#       axis.line = element_blank(),
#       plot.background = element_rect(fill = "white", color = NA)
#     )
#   
#   ggsave(
#     filename = save_path,
#     plot = plot,
#     device = "pdf",
#     width = 10, height = 6, units = "in"
#   )
# }
# 
# # Boucle principale
# for (sim in simulations) {
#   
#   # 1. Ouvrir le fichier
#   path_sim <- file.path(chem_EDGAR, sim, "GEOSChem.DryDep.2015_m.nc4")
#   if (!file.exists(path_sim)) next
#   nc_sim <- nc_open(path_sim)
#   
#   # 2. Lecture temps et poids
#   time_raw <- ncvar_get(nc_sim, "time")
#   origin_time <- as.POSIXct("2015-01-01 00:00:00", tz = "UTC")
#   time <- origin_time + time_raw * 60
#   weights <- days_in_month(time) / sum(days_in_month(time))
#   weighted_avg <- function(x) sum(x * weights)
#   
#   # 3. Pour chaque espèce
#   for (species in c("Hg0", "Hg2", "HgP")) {
#     
#     varname <- paste0("DryDep_", species)
#     message("Traitement de ", sim, " - ", varname)
#     
#     dep_raw <- ncvar_get(nc_sim, varname)  # [lon, lat, time]
#     dep_mean <- apply(dep_raw, MARGIN = c(1, 2), FUN = weighted_avg)
#     dep_Mg <- dep_mean * unit_conv_dd * kg_Mg
#     
#     # 4. Sauvegarde du plot
#     save_path <- file.path(
#       "C:/Users/colom/Desktop/STAGE/data/clean_mod_data/plot",
#       paste0("EDGAR_", sim, "_DryDep_", species, "_WORLD.pdf")
#     )
#     
#     plot_drydep_map(dep_Mg, sim, species, save_path)
#   }
#   
#   nc_close(nc_sim)
# }
```

``` r
simulations <- c("V03A03", "V20A03", "V20A20", "V03A20")
lon <- seq(-180, 180, length.out = 144)
lat <- seq(-90, 90, length.out = 91)
df_grid <- expand.grid(x = lon, y = lat)



# --- Paramètres graphiques ---
breaks_dep <- c(0.03, 0.05, 0.1, 0.3, 0.5, 0.8, 1, 2.5, 5, Inf)
labels_dep_chr <- c("0.03 - 0.05", "0.05 – 0.1", "0.1 – 0.3", "0.3 – 0.5", "0.5 – 0.8",
                    "0.8 – 1.0", "1.0 – 2.5", "2.5 – 5.0", "> 5.0")
labels_dep_expr <- c(expression("0.03 - 0.05"), expression("0.05 – 0.1"), expression("0.1 – 0.3"),
                     expression("0.3 – 0.5"), expression("0.5 – 0.8"), expression("0.8 – 1.0"),
                     expression("1.0 – 2.5"), expression("2.5 – 5.0"), expression("">" ~5.0"))
colors_dep <- c("#005F73", "#0A9396", "#94D2BD", "#E9D8A6",
                "#EE9B00", "#CA6702", "#BB3E03", "#AE2012", "#9B2226")

# --- Fonction pour générer une carte ---
plot_drydep_map <- function(matrix_dep, sim_name, save_path) {
  df <- df_grid
  df$dep_Mg <- as.vector(matrix_dep)
  df$dep_class <- cut(df$dep_Mg, breaks = breaks_dep, labels = labels_dep_chr, include.lowest = TRUE)
  df <- na.omit(df)

  plot <- ggplot(df, aes(x = x, y = y, fill = dep_class)) +
    geom_raster() +
    scale_fill_manual(
      name = "Dry Dep Hg⁰ [Mg]",
      values = colors_dep,
      labels = labels_dep_expr,
      na.value = "transparent",
      guide = guide_legend(reverse = TRUE)
    ) +
    coord_equal() +
    theme_minimal() +
    labs(
      title = "Annual Dry Deposition of Hg⁰ in 2015",
      subtitle = paste0("Simulation : ", sim_name),
      x = "Longitude", y = "Latitude"
    ) +
    borders("world", colour = "#403d39", size = 0.2) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8),
      axis.text = element_text(size = 10),
      axis.ticks = element_line(color = "black", size = 0.3),
      axis.ticks.length = unit(3, "pt")
    )

  ggsave(filename = save_path, plot = plot, device = "pdf", width = 10, height = 6, units = "in")
}

# --- Boucle principale ---
for (sim in simulations) {
  path_sim <- file.path(chem_EDGAR, sim, "GEOSChem.DryDep.2015_m.nc4")
  if (!file.exists(path_sim)) next

  nc_sim <- nc_open(path_sim)
  time_raw <- ncvar_get(nc_sim, "time")
  time <- as.POSIXct("2015-01-01 00:00:00", tz = "UTC") + time_raw * 60
  weights <- days_in_month(time) / sum(days_in_month(time))
  weighted_avg <- function(x) sum(x * weights)

  varname <- "DryDep_Hg0"
  message("Traitement de ", sim, " - ", varname)

  dep_raw <- ncvar_get(nc_sim, varname)
  dep_mean <- apply(dep_raw, MARGIN = c(1, 2), FUN = weighted_avg)
  dep_Mg <- dep_mean * unit_conv_dd * kg_Mg

  save_path <- file.path("C:/Users/colom/Desktop/STAGE/data/clean_mod_data/plot",
                         paste0("EDGAR_", sim, "_DryDep_Hg0_WORLD_TEST.pdf"))

  plot_drydep_map(dep_Mg, sim, save_path)
  nc_close(nc_sim)
}
```

``` r
# 
# df_grid <- expand.grid(x = lon, y = lat)
# 
# 
# 
# # --- Paramètres graphiques ---
# breaks_dep <- c(0.01, 0.05, 0.1, 0.3, 0.5, 0.8, 1, 1.5, 2, Inf)
# labels_dep_chr <- c("0.01–0.05", "0.05–0.1", "0.1–0.3", "0.3–0.5", "0.5–0.8",
#                     "0.8–1.0", "1.0–1.5", "1.5–2.0", ">2.0")
# labels_dep_expr <- c(expression("0.01–0.05"), expression("0.05–0.1"), expression("0.1–0.3"),
#                      expression("0.3–0.5"), expression("0.5–0.8"), expression("0.8–1.0"),
#                      expression("1.0–1.5"), expression("1.5–2.0"), expression("">"~2.0"))
# 
# colors_dep <- c("#005F73", "#0A9396", "#94D2BD", "#E9D8A6",
#                 "#EE9B00", "#CA6702", "#BB3E03", "#AE2012", "#9B2226")
# 
# # --- Fonction pour générer une carte ---
# plot_drydep_map <- function(matrix_dep, sim_name, save_path) {
#   df <- df_grid
#   df$dep_Mg <- as.vector(matrix_dep)
#   df$dep_class <- cut(df$dep_Mg, breaks = breaks_dep, labels = labels_dep_chr, include.lowest = TRUE)
#   df <- na.omit(df)
# 
#   plot <- ggplot(df, aes(x = x, y = y, fill = dep_class)) +
#     geom_raster() +
#     scale_fill_manual(
#       name = "Dry Dep Hg(II) [Mg]",
#       values = colors_dep,
#       labels = labels_dep_expr,
#       na.value = "transparent",
#       guide = guide_legend(reverse = TRUE)
#     ) +
#     coord_equal() +
#     theme_minimal() +
#     labs(
#       title = "Annual Dry Deposition of Hg(II) in 2015",
#       subtitle = paste0("Simulation : ", sim_name),
#       x = "Longitude", y = "Latitude"
#     ) +
#     borders("world", colour = "#403d39", size = 0.2) +
#     theme(
#       plot.title = element_text(hjust = 0.5, face = "bold"),
#       plot.subtitle = element_text(hjust = 0.5),
#       panel.grid = element_blank(),
#       panel.border = element_rect(color = "black", fill = NA, size = 0.8),
#       axis.text = element_text(size = 10),
#       axis.ticks = element_line(color = "black", size = 0.3),
#       axis.ticks.length = unit(3, "pt")
#     )
# 
#   ggsave(filename = save_path, plot = plot, device = "pdf", width = 10, height = 6, units = "in")
# }
# 
# # --- Boucle principale ---
# for (sim in simulations) {
#   path_sim <- file.path(chem_EDGAR, sim, "GEOSChem.DryDep.2015_m.nc4")
#   if (!file.exists(path_sim)) next
# 
#   nc_sim <- nc_open(path_sim)
#   time_raw <- ncvar_get(nc_sim, "time")
#   time <- as.POSIXct("2015-01-01 00:00:00", tz = "UTC") + time_raw * 60
#   weights <- days_in_month(time) / sum(days_in_month(time))
#   weighted_avg <- function(x) sum(x * weights)
# 
#   varname <- "DryDep_Hg2"
#   message("Traitement de ", sim, " - ", varname)
# 
#   dep_raw <- ncvar_get(nc_sim, varname)
#   dep_mean <- apply(dep_raw, MARGIN = c(1, 2), FUN = weighted_avg)
#   dep_Mg <- dep_mean * unit_conv_dd * kg_Mg
# 
#   save_path <- file.path("C:/Users/colom/Desktop/STAGE/data/clean_mod_data/plot",
#                          paste0("EDGAR_", sim, "_DryDep_Hg2_WORLD.pdf"))
# 
#   plot_drydep_map(dep_Mg, sim, save_path)
#   nc_close(nc_sim)
# }
```

``` r
# lon <- seq(-180, 180, length.out = 144)
# lat <- seq(-90, 90, length.out = 91)
# df_grid <- expand.grid(x = lon, y = lat)
# 
# 
# 
# # --- Paramètres graphiques ---
# breaks_dep <- c(0.01, 0.05, 0.1, 0.3, 0.5, 0.8, 1, 1.5, 2, Inf)
# labels_dep_chr <- c("0.01–0.05", "0.05–0.1", "0.1–0.3", "0.3–0.5", "0.5–0.8",
#                     "0.8–1.0", "1.0–1.5", "1.5–2.0", ">2.0")
# labels_dep_expr <- c(expression("0.01–0.05"), expression("0.05–0.1"), expression("0.1–0.3"),
#                      expression("0.3–0.5"), expression("0.5–0.8"), expression("0.8–1.0"),
#                      expression("1.0–1.5"), expression("1.5–2.0"), expression("">"~2.0"))
# 
# colors_dep <- c("#005F73", "#0A9396", "#94D2BD", "#E9D8A6",
#                 "#EE9B00", "#CA6702", "#BB3E03", "#AE2012", "#9B2226")
# 
# # --- Fonction pour générer une carte ---
# plot_drydep_map <- function(matrix_dep, sim_name, save_path) {
#   df <- df_grid
#   df$dep_Mg <- as.vector(matrix_dep)
#   df$dep_class <- cut(df$dep_Mg, breaks = breaks_dep, labels = labels_dep_chr, include.lowest = TRUE)
#   df <- na.omit(df)
# 
#   plot <- ggplot(df, aes(x = x, y = y, fill = dep_class)) +
#     geom_raster() +
#     scale_fill_manual(
#       name = "Dry Dep Hg(P) [Mg]",
#       values = colors_dep,
#       labels = labels_dep_expr,
#       na.value = "transparent",
#       guide = guide_legend(reverse = TRUE)
#     ) +
#     coord_equal() +
#     theme_minimal() +
#     labs(
#       title = "Annual Dry Deposition of Hg(P) in 2015",
#       subtitle = paste0("Simulation : ", sim_name),
#       x = "Longitude", y = "Latitude"
#     ) +
#     borders("world", colour = "#403d39", size = 0.2) +
#     theme(
#       plot.title = element_text(hjust = 0.5, face = "bold"),
#       plot.subtitle = element_text(hjust = 0.5),
#       panel.grid = element_blank(),
#       panel.border = element_rect(color = "black", fill = NA, size = 0.8),
#       axis.text = element_text(size = 10),
#       axis.ticks = element_line(color = "black", size = 0.3),
#       axis.ticks.length = unit(3, "pt")
#     )
# 
#   ggsave(filename = save_path, plot = plot, device = "pdf", width = 10, height = 6, units = "in")
# }
# 
# # --- Boucle principale ---
# for (sim in simulations) {
#   path_sim <- file.path(chem_EDGAR, sim, "GEOSChem.DryDep.2015_m.nc4")
#   if (!file.exists(path_sim)) next
# 
#   nc_sim <- nc_open(path_sim)
#   time_raw <- ncvar_get(nc_sim, "time")
#   time <- as.POSIXct("2015-01-01 00:00:00", tz = "UTC") + time_raw * 60
#   weights <- days_in_month(time) / sum(days_in_month(time))
#   weighted_avg <- function(x) sum(x * weights)
# 
#   varname <- "DryDep_HgP"
#   message("Traitement de ", sim, " - ", varname)
# 
#   dep_raw <- ncvar_get(nc_sim, varname)
#   dep_mean <- apply(dep_raw, MARGIN = c(1, 2), FUN = weighted_avg)
#   dep_Mg <- dep_mean * unit_conv_dd * kg_Mg
# 
#   save_path <- file.path("C:/Users/colom/Desktop/STAGE/data/clean_mod_data/plot",
#                          paste0("EDGAR_", sim, "_DryDep_HgP_WORLD.pdf"))
# 
#   plot_drydep_map(dep_Mg, sim, save_path)
#   nc_close(nc_sim)
# }
```

``` r
df_grid <- expand.grid(x = lon, y = lat)


# --- Paramètres graphiques ---
breaks_dep <- c(0.01, 0.05, 0.1, 0.3, 0.5, 0.8, 1, 1.5, 2, Inf)
labels_dep_chr <- c("0.01–0.05", "0.05–0.1", "0.1–0.3", "0.3–0.5", "0.5–0.8",
                    "0.8–1.0", "1.0–1.5", "1.5–2.0", ">2.0")
labels_dep_expr <- c(expression("0.01–0.05"), expression("0.05–0.1"), expression("0.1–0.3"),
                     expression("0.3–0.5"), expression("0.5–0.8"), expression("0.8–1.0"),
                     expression("1.0–1.5"), expression("1.5–2.0"), expression("">"~2.0"))
colors_dep <- c("#005F73", "#0A9396", "#94D2BD", "#E9D8A6",
                "#EE9B00", "#CA6702", "#BB3E03", "#AE2012", "#9B2226")

# --- Fonction de carte ---
plot_drydep_map <- function(matrix_dep, sim_name, save_path) {
  df <- df_grid
  df$dep_Mg <- as.vector(matrix_dep)
  df$dep_class <- cut(df$dep_Mg, breaks = breaks_dep, labels = labels_dep_chr, include.lowest = TRUE)
  df <- na.omit(df)

  plot <- ggplot(df, aes(x = x, y = y, fill = dep_class)) +
    geom_raster() +
    scale_fill_manual(
      name = "Dry Dep Hg(II) + Hg(P) [Mg]",
      values = colors_dep,
      labels = labels_dep_expr,
      na.value = "transparent",
      guide = guide_legend(reverse = TRUE)
    ) +
    coord_equal() +
    theme_minimal() +
    labs(
      title = "Annual Dry Deposition of Hg(II) + Hg(P) in 2015",
      subtitle = paste0("Simulation : ", sim_name),
      x = "Longitude", y = "Latitude"
    ) +
    borders("world", colour = "#403d39", size = 0.2) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8),
      axis.text = element_text(size = 10),
      axis.ticks = element_line(color = "black", size = 0.3),
      axis.ticks.length = unit(3, "pt")
    )

  ggsave(filename = save_path, plot = plot, device = "pdf", width = 10, height = 6, units = "in")
}

# --- Boucle principale ---
for (sim in simulations) {
  path_sim <- file.path(chem_EDGAR, sim, "GEOSChem.DryDep.2015_m.nc4")
  if (!file.exists(path_sim)) next

  message("Traitement de ", sim)

  nc_sim <- nc_open(path_sim)
  time_raw <- ncvar_get(nc_sim, "time")
  time <- as.POSIXct("2015-01-01 00:00:00", tz = "UTC") + time_raw * 60
  weights <- days_in_month(time) / sum(days_in_month(time))
  weighted_avg <- function(x) sum(x * weights)

  # Lecture Hg2 et HgP
  dep_Hg2_raw <- ncvar_get(nc_sim, "DryDep_Hg2")
  dep_HgP_raw <- ncvar_get(nc_sim, "DryDep_HgP")

  dep_Hg2_mean <- apply(dep_Hg2_raw, MARGIN = c(1, 2), FUN = weighted_avg)
  dep_HgP_mean <- apply(dep_HgP_raw, MARGIN = c(1, 2), FUN = weighted_avg)

  # Somme et conversion
  dep_total_Mg <- (dep_Hg2_mean + dep_HgP_mean) * unit_conv_dd * kg_Mg

  save_path <- file.path("C:/Users/colom/Desktop/STAGE/data/clean_mod_data/plot",
                         paste0("EDGAR_", sim, "_DryDep_Hg2HgP_WORLD.pdf"))

  plot_drydep_map(dep_total_Mg, sim, save_path)
  nc_close(nc_sim)
}
```

``` r
#WD_V03A03<-nc_open(file.path(chem_EDGAR, "V03A03/GEOSChem.WetLossTot_HgSum.2015_m.nc4"))
#print(WD_V03A03)
# WetLossTot_V03A03<-ncvar_get(WD_V03A03,"WetLossTot_HgSum")
# range(WetLossTot_V03A03)
# weighted_avg <- function(x, weights) sum(x * weights)
#  time_raw <- ncvar_get(WD_V03A03, "time")
#   origin_time <- as.POSIXct("2015-01-01 00:00:00", tz = "UTC")
#   time <- origin_time + time_raw * 60
#   days_in_months <- days_in_month(time)
#   weights <- days_in_months / sum(days_in_months)
#   WetLossTot_V03A03 <- apply(WetLossTot_V03A03, c(1, 2), weighted_avg, weights = weights)
#   WetLossTot_V03A03 <- WetLossTot_V03A03 * s_in_yr * kg_Mg  
# range(WetLossTot_V03A03)  
# dim(WetLossTot_V03A03)
# lon <- seq(-180, 180, length.out = 144)
# lat <- seq(-90, 90, length.out = 91)
# df_grid <- expand.grid(x = lon, y = lat)
# 
# breaks_dep <- c(0.001, 0.005, 0.01,0.1, 0.25, 0.4, 0.5, 0.75, Inf)
# 
# labels_dep_chr <- c(
#   "0.001–0.005", "0.005–0.01", "0.01–0.1",
#   "0.1–0.25", "0.25–0.5", "0.5–0.75", "0.75–1.0", ">1.0"
# )
# 
# labels_dep_expr <- c(
#   expression("1×10"^{-3}*" – 5×10"^{-3}),
#   expression("5×10"^{-3}*" – 1×10"^{-2}),
#   expression("1×10"^{-2}*" – 0.1"),
#   expression("0.1 – 0.25"),
#   expression("0.25 – 0.4"),
#   expression("0.4 – 0.5"),
#   expression("0.5 – 0.75"),
#   expression("">"~0.75")
# )
# 
# colors_dep <- c("#005F73", "#0A9396", "#94D2BD", "#E9D8A6",
#                 "#EE9B00", "#CA6702", "#BB3E03", "#AE2012")
# 
# df <- df_grid
# df$dep_Mg <- as.vector(WetLossTot_V03A03)
# df$dep_class <- cut(df$dep_Mg, breaks = breaks_dep, labels = labels_dep_chr, include.lowest = TRUE)
# df <- na.omit(df)
# 
# 
# wd_V03A03_plot <-ggplot(df, aes(x = x, y = y, fill = dep_class)) +
#   geom_raster() +
#   scale_fill_manual(
#     name = "Wet Dep Total Hg [Mg]",
#     values = colors_dep,
#     labels = labels_dep_expr,
#     na.value = "transparent",
#     guide = guide_legend(reverse = TRUE)
#   ) +
#   coord_equal() +
#   theme_minimal() +
#   labs(
#     title = "Annual Wet Deposition of Hg Total in 2015",
#     subtitle = "Dataset : EDGAR / Simulation : V03A03",
#     x = "Longitude", y = "Latitude"
#   ) +
#   borders("world", colour = "#403d39", size = 0.2) +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold"),
#     plot.subtitle = element_text(hjust = 0.5),
#     panel.grid = element_blank(),
#     panel.border = element_rect(color = "black", fill = NA, size = 0.8),
#     axis.text = element_text(size = 10),
#     axis.ticks = element_line(color = "black", size = 0.3),
#     axis.ticks.length = unit(3, "pt")
#   )
# 
# print(wd_V03A03_plot)
# chemin_pdf <- "C:/Users/colom/Desktop/STAGE/data/clean_mod_data/plot/EDGAR_WD_Total_V03A03.pdf"
# 
# # Enregistrement du plot
# ggsave(
#   filename = chemin_pdf,
#   plot = wd_V03A03_plot,
#   device = "pdf",
#   width = 10, height = 6, units = "in"
# )
```

``` r
lon <- seq(-180, 180, length.out = 144)
lat <- seq(-90, 90, length.out = 91)
df_grid <- expand.grid(x = lon, y = lat)


# --- Paramètres graphiques ---
breaks_dep <- c(0.001, 0.005, 0.01, 0.1, 0.25, 0.4, 0.5, 0.75, Inf)
labels_dep_chr <- c("0.001–0.005", "0.005–0.01", "0.01–0.1",
                    "0.1–0.25", "0.25–0.4", "0.4–0.5", "0.5–0.75", ">0.75")
labels_dep_expr <- c(
  expression("1×10"^{-3}*" – 5×10"^{-3}),
  expression("5×10"^{-3}*" – 1×10"^{-2}),
  expression("1×10"^{-2}*" – 0.1"),
  expression("0.1 – 0.25"),
  expression("0.25 – 0.4"),
  expression("0.4 – 0.5"),
  expression("0.5 – 0.75"),
  expression("">"~0.75")
)
colors_dep <- c("#005F73", "#0A9396", "#94D2BD", "#E9D8A6",
                "#EE9B00", "#CA6702", "#BB3E03", "#AE2012")

# --- Fonction pour générer la carte ---
plot_wetdep_map <- function(matrix_dep, sim_name, save_path) {
  df <- df_grid
  df$dep_Mg <- as.vector(matrix_dep)
  df$dep_class <- cut(df$dep_Mg, breaks = breaks_dep, labels = labels_dep_chr, include.lowest = TRUE)
  df <- na.omit(df)

  plot <- ggplot(df, aes(x = x, y = y, fill = dep_class)) +
    geom_raster() +
    scale_fill_manual(
      name = "Wet Dep Total Hg [Mg]",
      values = colors_dep,
      labels = labels_dep_expr,
      na.value = "transparent",
      guide = guide_legend(reverse = TRUE)
    ) +
    coord_equal() +
    theme_minimal() +
    labs(
      title = "Annual Wet Deposition of Hg Total in 2015",
      subtitle = paste0("Simulation : ", sim_name),
      x = "Longitude", y = "Latitude"
    ) +
    borders("world", colour = "#403d39", size = 0.2) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8),
      axis.text = element_text(size = 10),
      axis.ticks = element_line(color = "black", size = 0.3),
      axis.ticks.length = unit(3, "pt")
    )

  ggsave(filename = save_path, plot = plot, device = "pdf", width = 10, height = 6, units = "in")
}

# --- Boucle principale ---
for (sim in simulations) {
  path_sim <- file.path(chem_EDGAR, sim, "GEOSChem.WetLossTot_HgSum.2015_m.nc4")
  if (!file.exists(path_sim)) next

  nc_sim <- nc_open(path_sim)
  time_raw <- ncvar_get(nc_sim, "time")
  time <- as.POSIXct("2015-01-01 00:00:00", tz = "UTC") + time_raw * 60
  weights <- days_in_month(time) / sum(days_in_month(time))
  weighted_avg <- function(x) sum(x * weights)

  varname <- "WetLossTot_HgSum"
  message("Traitement de ", sim, " - ", varname)

  dep_raw <- ncvar_get(nc_sim, varname)
  dep_mean <- apply(dep_raw, MARGIN = c(1, 2), FUN = weighted_avg)
  dep_Mg <- dep_mean * s_in_yr * kg_Mg  # Conversion kg/s → Mg/an

  save_path <- file.path("C:/Users/colom/Desktop/STAGE/data/clean_mod_data/plot",
                         paste0("EDGAR_", sim, "_WetDep_HgTot_WORLD.pdf"))

  plot_wetdep_map(dep_Mg, sim, save_path)
  nc_close(nc_sim)
}
```

``` r
WD_V03A03<-nc_open(file.path(chem_hist,"EDGAR/V03A03/GEOSChem.WetLossTot_HgSum.2015_m.nc4"))
WD_V20A20<-nc_open(file.path(chem_hist,"EDGAR/V20A20/GEOSChem.WetLossTot_HgSum.2015_m.nc4"))
print(WD_V20A20)
```

    ## File C:/Users/colom/Desktop/STAGE/data/clean_mod_data/14_3_1/EDGAR/V20A20/GEOSChem.WetLossTot_HgSum.2015_m.nc4 (NC_FORMAT_NETCDF4_CLASSIC):
    ## 
    ##      3 variables (excluding dimension variables):
    ##         double lon_bnds[bnds,lon]   (Contiguous storage)  
    ##         double lat_bnds[bnds,lat]   (Contiguous storage)  
    ##         float WetLossTot_HgSum[lon,lat,time]   (Chunking: [144,91,1])  
    ##             _FillValue: -8.99999987309029e+33
    ##             missing_value: -8.99999987309029e+33
    ## 
    ##      4 dimensions:
    ##         time  Size:12   *** is unlimited *** 
    ##             standard_name: time
    ##             long_name: Time
    ##             units: minutes since 2015-01-01 00:00:00
    ##             calendar: gregorian
    ##             axis: T
    ##         lon  Size:144 
    ##             standard_name: longitude
    ##             long_name: Longitude
    ##             units: degrees_east
    ##             axis: X
    ##             bounds: lon_bnds
    ##         bnds  Size:2 (no dimvar)
    ##         lat  Size:91 
    ##             standard_name: latitude
    ##             long_name: Latitude
    ##             units: degrees_north
    ##             axis: Y
    ##             bounds: lat_bnds
    ## 
    ##     12 global attributes:
    ##         CDI: Climate Data Interface version 2.1.1 (https://mpimet.mpg.de/cdi)
    ##         Conventions: CF-1.6
    ##         title: GEOS-Chem diagnostic collection: WetLossLS
    ##         history: Fri Aug 01 10:15:07 2025: cdo expr,WetLossTot_HgSum = WetLossTot_HgLS + WetLossTot_HgConv GEOSChem.WetLossMerge.2015.nc4 GEOSChem.WetLossTot_HgSum.2015_m.nc4
    ## Fri Aug 01 10:15:02 2025: cdo merge GEOSChem.WetLossTotalLS.2015_m.nc4 GEOSChem.WetLossTotal.2015_m.nc4 GEOSChem.WetLossMerge.2015.nc4
    ## Fri Aug 01 10:13:26 2025: cdo selyear,2015 GEOSChem.WetLossTotalLS.alltime_m.nc4 GEOSChem.WetLossTotalLS.2015_m.nc4
    ## Fri Aug 01 10:13:20 2025: cdo vertsum LStemp_allwetdep.nc4 GEOSChem.WetLossTotalLS.alltime_m.nc4
    ## Fri Aug 01 10:13:10 2025: cdo expr,WetLossTot_HgLS = WetLossLS_Hg2ORGP + WetLossLS_Hg2ClP + WetLossLS_HgCl2 + WetLossLS_HgOHOH + WetLossLS_HgOHBrO + WetLossLS_HgOHClO + WetLossLS_HgOHHO2 + WetLossLS_HgOHNO2 + WetLossLS_HgClOH + WetLossLS_HgClBr + WetLossLS_HgClBrO + WetLossLS_HgClClO + WetLossLS_HgClHO2 + WetLossLS_HgClNO2 + WetLossLS_HgBr2 + WetLossLS_HgBrOH + WetLossLS_HgBrClO + WetLossLS_HgBrBrO + WetLossLS_HgBrHO2 + WetLossLS_HgBrNO2 LSfinal.nc4 LStemp_allwetdep.nc4
    ## Fri Aug 01 10:12:58 2025: cdo merge LStemp1.nc4 LStemp2.nc4 LStemp3.nc4 LStemp4.nc4 LSfinal.nc4
    ## Fri Aug 01 10:12:42 2025: cdo selvar,WetLossLS_Hg2ORGP,WetLossLS_Hg2ClP,WetLossLS_HgCl2,WetLossLS_HgOHOH,WetLossLS_HgOHBrO GEOSChem.WetLossLS.alltime.nc4 LStemp1.nc4
    ## Fri Aug 01 10:12:28 2025: cdo mergetime GEOSChem.WetLossLS.20150101_0000z.nc4 GEOSChem.WetLossLS.20150201_0000z.nc4 GEOSChem.WetLossLS.20150301_0000z.nc4 GEOSChem.WetLossLS.20150401_0000z.nc4 GEOSChem.WetLossLS.20150501_0000z.nc4 GEOSChem.WetLossLS.20150601_0000z.nc4 GEOSChem.WetLossLS.20150701_0000z.nc4 GEOSChem.WetLossLS.20150801_0000z.nc4 GEOSChem.WetLossLS.20150901_0000z.nc4 GEOSChem.WetLossLS.20151001_0000z.nc4 GEOSChem.WetLossLS.20151101_0000z.nc4 GEOSChem.WetLossLS.20151201_0000z.nc4 GEOSChem.WetLossLS.alltime.nc4
    ##         format: NetCDF-4
    ##         conventions: COARDS
    ##         ProdDateTime: 
    ##         reference: www.geos-chem.org; wiki.geos-chem.org
    ##         contact: GEOS-Chem Support Team (geos-chem-support@g.harvard.edu)
    ##         simulation_start_date_and_time: 2015-01-01 00:00:00z
    ##         simulation_end_date_and_time: 2016-01-01 00:00:00z
    ##         CDO: Climate Data Operators version 2.1.1 (https://mpimet.mpg.de/cdo)

``` r
WD_V03A03_var <- ncvar_get(WD_V03A03, "WetLossTot_HgSum")
WD_V20A20_var <- ncvar_get(WD_V20A20, "WetLossTot_HgSum")
dim(WD_V03A03_var)
```

    ## [1] 144  91  12

``` r
    weighted_avg <- function(x, weights) sum(x * weights)
    time_raw <- ncvar_get(WD_V03A03, "time")
    origin_time <- as.POSIXct("2015-01-01 00:00:00", tz = "UTC")
    time <- origin_time + time_raw * 60
    days_in_months <- days_in_month(time)
    weights <- days_in_months / sum(days_in_months)
    WD_V03A03_var <- apply(WD_V03A03_var, c(1, 2), weighted_avg, weights = weights)
    WD_V20A20_var <- apply(WD_V20A20_var, c(1, 2), weighted_avg, weights = weights)
range(WD_V03A03_var)
```

    ## [1] 2.033457e-09 3.169220e-05

``` r
range(WD_V20A20_var)
```

    ## [1] 2.138594e-09 4.249527e-05

``` r
dif_wd <- (WD_V20A20_var - WD_V03A03_var) * s_in_yr * kg_Mg  # Conversion kg/s → Mg/an
range(dif_wd)
```

    ## [1] -0.1990711  0.6474713

``` r
mean(dif_wd)
```

    ## [1] 0.0113549

``` r
lon <- seq(-180, 180, length.out = 144)
lat <- seq(-90, 90, length.out = 91)
df_diff <- expand.grid(x = lon, y = lat)
df_diff$diff_conc <- as.vector(dif_wd)

df_diff <- na.omit(df_diff)
colors_div <- c("#2166AC", "#FFFFFF", "#B2182B")  # bleu - blanc - rouge



plot_diff <- ggplot(df_diff, aes(x = x, y = y, fill = diff_conc)) +
  geom_raster() +
  scale_fill_gradient2(
    name = "Différence Hg(0) (ng/m³)\n2020 - 2003",
    low = colors_div[1],
    mid = colors_div[2],
    high = colors_div[3],
    midpoint = 0,
    limits = c(-max(abs(df_diff$diff_conc)), max(abs(df_diff$diff_conc))),
    oob = scales::squish
  ) +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Average annual difference of Hg(0) between 2020 and 2003",
    x = "Longitude", y = "Latitude"
  ) +
  borders("world", colour = "#403d39", size = 0.2) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_text(size = 10),
    axis.ticks = element_line(color = "black", size = 0.3),
    axis.ticks.length = unit(3, "pt"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.line = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  )

plot_diff
```

![](EDGAR_Map_Dep_05_08_2025_files/figure-gfm/diff%20wetdep%202020%202003-1.png)<!-- -->

``` r
lim <- 0.2
plot_diff_zoom <- ggplot(df_diff, aes(x = x, y = y, fill = diff_conc)) +
  geom_raster() +
  scale_fill_gradient2(
    name = "Difference Hg(0) Wetdep [Mg]\n2020 - 2003",
    low = colors_div[1],
    mid = colors_div[2],
    high = colors_div[3],
    midpoint = 0,
    limits = c(-lim, lim),
    oob = scales::squish
  ) +
  coord_equal() +
  theme_minimal() +
  labs(
    title = paste0("Average annual difference of Hg(0) wetdep between 2020 and 2003"),
        subtitle = "Dataset : EDGAR",
    x = "Longitude", y = "Latitude"
  ) +
  borders("world", colour = "#403d39", size = 0.2) +
   theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text = element_text(size = 10),
    axis.ticks = element_line(color = "black", size = 0.3),
    axis.ticks.length = unit(3, "pt"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.line = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  )

plot_diff_zoom
```

![](EDGAR_Map_Dep_05_08_2025_files/figure-gfm/diff%20wetdep%202020%202003-2.png)<!-- -->

``` r
chemin_pdf <- "C:/Users/colom/Desktop/STAGE/data/clean_mod_data/plot/EDGAR_Diff_Wetdep_20_03.pdf"

# Enregistrement du plot
ggsave(
  filename = chemin_pdf,
  plot = plot_diff_zoom,
  device = "pdf",
  width = 10, height = 6, units = "in"
)
```

``` r
dd_V03A03<-nc_open(file.path(chem_hist,"EDGAR/V03A03/GEOSChem.DryDep.2015_m.nc4"))
dd_V20A20<-nc_open(file.path(chem_hist,"EDGAR/V20A20/GEOSChem.DryDep.2015_m.nc4"))
print(dd_V03A03)
```

    ## File C:/Users/colom/Desktop/STAGE/data/clean_mod_data/14_3_1/EDGAR/V03A03/GEOSChem.DryDep.2015_m.nc4 (NC_FORMAT_NETCDF4_CLASSIC):
    ## 
    ##      5 variables (excluding dimension variables):
    ##         double lon_bnds[bnds,lon]   (Contiguous storage)  
    ##         double lat_bnds[bnds,lat]   (Contiguous storage)  
    ##         float DryDep_Hg0[lon,lat,time]   (Chunking: [144,91,1])  
    ##             long_name: Dry deposition flux of species Hg0
    ##             units: molec cm-2 s-1
    ##             averaging_method: time-averaged
    ##         float DryDep_Hg2[lon,lat,time]   (Chunking: [144,91,1])  
    ##             _FillValue: -8.99999987309029e+33
    ##             missing_value: -8.99999987309029e+33
    ##         float DryDep_HgP[lon,lat,time]   (Chunking: [144,91,1])  
    ##             _FillValue: -8.99999987309029e+33
    ##             missing_value: -8.99999987309029e+33
    ## 
    ##      4 dimensions:
    ##         time  Size:12   *** is unlimited *** 
    ##             standard_name: time
    ##             long_name: Time
    ##             units: minutes since 2015-01-01 00:00:00
    ##             calendar: gregorian
    ##             axis: T
    ##         lon  Size:144 
    ##             standard_name: longitude
    ##             long_name: Longitude
    ##             units: degrees_east
    ##             axis: X
    ##             bounds: lon_bnds
    ##         bnds  Size:2 (no dimvar)
    ##         lat  Size:91 
    ##             standard_name: latitude
    ##             long_name: Latitude
    ##             units: degrees_north
    ##             axis: Y
    ##             bounds: lat_bnds
    ## 
    ##     12 global attributes:
    ##         CDI: Climate Data Interface version 2.1.1 (https://mpimet.mpg.de/cdi)
    ##         Conventions: CF-1.6
    ##         title: GEOS-Chem diagnostic collection: DryDep
    ##         history: Fri Aug 01 10:03:21 2025: cdo mergetime DryDep_Hg_all.nc4 GEOSChem.DryDep.2015_m.nc4
    ## Fri Aug 01 10:03:17 2025: cdo merge DryDep_Hg0.nc DryDep_Hg2.nc DryDep_HgP.nc DryDep_Hg_all.nc4
    ## Fri Aug 01 10:03:05 2025: cdo selname,DryDep_Hg0 GEOSChem.DryDep.alltime.nc4 DryDep_Hg0.nc
    ## Fri Aug 01 10:03:02 2025: cdo mergetime GEOSChem.DryDep.20150101_0000z.nc4 GEOSChem.DryDep.20150201_0000z.nc4 GEOSChem.DryDep.20150301_0000z.nc4 GEOSChem.DryDep.20150401_0000z.nc4 GEOSChem.DryDep.20150501_0000z.nc4 GEOSChem.DryDep.20150601_0000z.nc4 GEOSChem.DryDep.20150701_0000z.nc4 GEOSChem.DryDep.20150801_0000z.nc4 GEOSChem.DryDep.20150901_0000z.nc4 GEOSChem.DryDep.20151001_0000z.nc4 GEOSChem.DryDep.20151101_0000z.nc4 GEOSChem.DryDep.20151201_0000z.nc4 GEOSChem.DryDep.alltime.nc4
    ##         format: NetCDF-4
    ##         conventions: COARDS
    ##         ProdDateTime: 
    ##         reference: www.geos-chem.org; wiki.geos-chem.org
    ##         contact: GEOS-Chem Support Team (geos-chem-support@g.harvard.edu)
    ##         simulation_start_date_and_time: 2015-01-01 00:00:00z
    ##         simulation_end_date_and_time: 2016-01-01 00:00:00z
    ##         CDO: Climate Data Operators version 2.1.1 (https://mpimet.mpg.de/cdo)

``` r
dep_Hg2_V03A03 <- ncvar_get(dd_V03A03, "DryDep_Hg2")
dep_HgP_V03A03 <- ncvar_get(dd_V03A03, "DryDep_HgP")
dep_Hg0_V03A03 <- ncvar_get(dd_V03A03, "DryDep_Hg0")

dep_Hg2_V20A20 <- ncvar_get(dd_V20A20, "DryDep_Hg2")
dep_HgP_V20A20 <- ncvar_get(dd_V20A20, "DryDep_HgP")
dep_Hg0_V20A20 <- ncvar_get(dd_V20A20, "DryDep_Hg0")


dim(dep_Hg2_V03A03)
```

    ## [1] 144  91  12

``` r
    weighted_avg <- function(x, weights) sum(x * weights)
    time_raw <- ncvar_get(dd_V03A03, "time")
    origin_time <- as.POSIXct("2015-01-01 00:00:00", tz = "UTC")
    time <- origin_time + time_raw * 60
    days_in_months <- days_in_month(time)
    weights <- days_in_months / sum(days_in_months)
    dep_Hg2_V03A03 <- apply(dep_Hg2_V03A03, c(1, 2), weighted_avg, weights = weights)
    dep_HgP_V03A03 <- apply(dep_HgP_V03A03, c(1, 2), weighted_avg, weights = weights)
    dep_Hg0_V03A03 <- apply(dep_Hg0_V03A03, c(1, 2), weighted_avg, weights = weights)

    dep_Hg2_V20A20 <- apply(dep_Hg2_V20A20, c(1, 2), weighted_avg, weights = weights)
    dep_HgP_V20A20 <- apply(dep_HgP_V20A20, c(1, 2), weighted_avg, weights = weights)
    dep_Hg0_V20A20 <- apply(dep_Hg0_V20A20, c(1, 2), weighted_avg, weights = weights)
    
    
    range(dep_Hg2_V20A20)
```

    ## [1]   1254.517 356902.646

``` r
dif_dd <- ((dep_Hg2_V20A20+dep_HgP_V20A20+dep_Hg0_V20A20)-(dep_Hg2_V03A03+dep_HgP_V03A03+dep_Hg0_V03A03)) * unit_conv_dd * kg_Mg  # Conversion kg/s → Mg/an

range(dif_dd)
```

    ## [1] -1.677868  6.125032

``` r
mean(dif_dd)
```

    ## [1] 0.01454966

``` r
lon <- seq(-180, 180, length.out = 144)
lat <- seq(-90, 90, length.out = 91)
df_diff <- expand.grid(x = lon, y = lat)
df_diff$diff_conc <- as.vector(dif_dd)

df_diff <- na.omit(df_diff)
colors_div <- c("#2166AC", "#FFFFFF", "#B2182B")  # bleu - blanc - rouge



plot_diff <- ggplot(df_diff, aes(x = x, y = y, fill = diff_conc)) +
  geom_raster() +
  scale_fill_gradient2(
    name = "Différence Hg(0) [Mg]\n2020 - 2003",
    low = colors_div[1],
    mid = colors_div[2],
    high = colors_div[3],
    midpoint = 0,
    limits = c(-max(abs(df_diff$diff_conc)), max(abs(df_diff$diff_conc))),
    oob = scales::squish
  ) +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Average annual difference of Hg drydep total between 2020 and 2003",
    x = "Longitude", y = "Latitude"
  ) +
  borders("world", colour = "#403d39", size = 0.2) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_text(size = 10),
    axis.ticks = element_line(color = "black", size = 0.3),
    axis.ticks.length = unit(3, "pt"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.line = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  )

plot_diff
```

![](EDGAR_Map_Dep_05_08_2025_files/figure-gfm/diff%20drydep%202020%202003-1.png)<!-- -->

``` r
lim <- 0.5
plot_diff_zoom <- ggplot(df_diff, aes(x = x, y = y, fill = diff_conc)) +
  geom_raster() +
  scale_fill_gradient2(
    name = "Difference Hg(0) Drydep [Mg]\n2020 - 2003",
    low = colors_div[1],
    mid = colors_div[2],
    high = colors_div[3],
    midpoint = 0,
    limits = c(-lim, lim),
    oob = scales::squish
  ) +
  coord_equal() +
  theme_minimal() +
  labs(
    title = paste0("Average annual difference of Hg drydep total between 2020 and 2003"),
        subtitle = "Dataset : EDGAR",
    x = "Longitude", y = "Latitude"
  ) +
  borders("world", colour = "#403d39", size = 0.2) +
   theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text = element_text(size = 10),
    axis.ticks = element_line(color = "black", size = 0.3),
    axis.ticks.length = unit(3, "pt"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.line = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  )

plot_diff_zoom
```

![](EDGAR_Map_Dep_05_08_2025_files/figure-gfm/diff%20drydep%202020%202003-2.png)<!-- -->

``` r
chemin_pdf <- "C:/Users/colom/Desktop/STAGE/data/clean_mod_data/plot/EDGAR_Diff_Drydep_20_03.pdf"

# Enregistrement du plot
ggsave(
  filename = chemin_pdf,
  plot = plot_diff_zoom,
  device = "pdf",
  width = 10, height = 6, units = "in"
)
```

``` r
simulations <- c("V03A03", "V03A20", "V20A03", "V20A20")

results_amz <- list()

for (sim in simulations) {
  
  # 1. Ouvrir le fichier NetCDF
  path_sim <- file.path(chem_EDGAR, sim, "GEOSChem.DryDep.2015_m.nc4")
  nc_sim <- nc_open(path_sim)
  
  # 2. Lire la variable Hg0 (niveau 1)
dep_Hg2<- ncvar_get(nc_sim, "DryDep_Hg2")
dep_HgP<- ncvar_get(nc_sim, "DryDep_HgP")
dep_Hg0<- ncvar_get(nc_sim, "DryDep_Hg0")
  
  
  
  # 3. Gérer le temps et pondérer
  time_raw <- ncvar_get(nc_sim, "time")
  origin_time <- as.POSIXct("2015-01-01 00:00:00", tz = "UTC")
  time <- origin_time + time_raw * 60
  weights <- days_in_month(time) / sum(days_in_month(time))
  
  # 4. Moyenne pondérée temporelle
  weighted_avg <- function(x) sum(x * weights)
  dep_Hg2_mean <- apply(dep_Hg2, MARGIN = c(1, 2), FUN = weighted_avg)
  dep_HgP_mean <- apply(dep_HgP, MARGIN = c(1, 2), FUN = weighted_avg)
  dep_Hg0_mean <- apply(dep_Hg0, MARGIN = c(1, 2), FUN = weighted_avg)
  
  
  # 5. Conversion en MG/an
  dep_total_mg_contrib <-(dep_Hg2_mean+dep_HgP_mean+dep_Hg0_mean) * unit_conv_dd * kg_Mg  # Conversion kg/s → Mg/an
  
  # 6. Appliquer le masque Amazonie
  dep_total_mg_contrib_amz <- dep_total_mg_contrib * Am_mask
  
  # 7. Stocker le résultat
  results_amz[[sim]] <- dep_total_mg_contrib_amz
  
  # Fermer le fichier
  nc_close(nc_sim)
}


# Résultat : une liste results_amz avec 4 objets nommés :
# results_amz$V03A03, results_amz$V03A20, etc.
V03A03_am <- results_amz[["V03A03"]]
V03A20_am <- results_amz[["V03A20"]]
V20A03_am <- results_amz[["V20A03"]]
V20A20_am <- results_amz[["V20A20"]]

range(V03A03_am)
```

    ## [1] 0.000000 6.819919

``` r
range(V03A20_am)
```

    ## [1]  0.00000 12.72607

``` r
range(V20A03_am)
```

    ## [1] 0.000000 6.816169

``` r
range(V20A20_am)
```

    ## [1]  0.000 12.849

``` r
#1. Décomposition spatiale : SCREEN DANS DOCS GEOSCHEM
Eff_A <- 0.5 * ((V03A20_am - V03A03_am) + (V20A20_am - V20A03_am))
Eff_V <- 0.5 * ((V20A03_am - V03A03_am) + (V20A20_am - V03A20_am))
Interaction <- V20A20_am - V03A03_am - Eff_A - Eff_V

#2. Extraction des valeurs sur l'Amazonie uniquement :

Eff_A_amz <- Eff_A[Am_mask == 1]
Eff_V_amz <- Eff_V[Am_mask == 1]
Interaction_amz <- Interaction[Am_mask == 1]

#3. Moyennes régionales :
mean_A <- mean(Eff_A_amz, na.rm = TRUE)
mean_V <- mean(Eff_V_amz, na.rm = TRUE)
mean_I <- mean(Interaction_amz, na.rm = TRUE)

#4. Camembert avec interaction répartie :
total_A <- mean_A + 0.5 * mean_I
total_V <- mean_V + 0.5 * mean_I
total_sum <- total_A + total_V

contrib <- c(
  Emissions = 100 * total_A / total_sum,
  Végétation  = 100 * total_V / total_sum
)



df <- data.frame(
  source = names(contrib),
  pct = contrib
)

df$label <- paste0(df$source, "\n", sprintf("%.1f", df$pct), "%")

contrib_conc<-ggplot(df, aes(x = "", y = pct, fill = source)) +
  geom_col(width = 1, color = "white") +
  geom_text(aes(label = label), 
            position = position_stack(vjust = 0.5), 
            color = "white", size = 5, fontface = "bold") +
  coord_polar("y") +
  scale_fill_manual(values = c("#D95F02", "#1B9E77")) +
  labs(
    title = "Relative contribution to the change in Hg total drydep",
    subtitle = "Dataset : EDGAR / Amazonia forest focus",
    fill = "Source"
  ) +
  theme_void() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    plot.background = element_rect(fill = "white", color = NA)
  )
print(contrib_conc)
```

![](EDGAR_Map_Dep_05_08_2025_files/figure-gfm/contrib%20drydep%20amz-1.png)<!-- -->

``` r
chemin_pdf <- "C:/Users/colom/Desktop/STAGE/data/clean_mod_data/plot/EDGAR_Contrib_Drydep_amz.pdf"

# Enregistrement du plot
ggsave(
  filename = chemin_pdf,
  plot = contrib_conc,
  device = "pdf",
  width = 10, height = 6, units = "in"
)
```

``` r
simulations <- c("V03A03", "V03A20", "V20A03", "V20A20")

results_amz <- list()

for (sim in simulations) {
  
  path_sim <- file.path(chem_EDGAR, sim, "GEOSChem.DryDep.2015_m.nc4")
  nc_sim <- nc_open(path_sim)
  
  dep_Hg2 <- ncvar_get(nc_sim, "DryDep_Hg2")
  dep_HgP <- ncvar_get(nc_sim, "DryDep_HgP")
  dep_Hg0 <- ncvar_get(nc_sim, "DryDep_Hg0")
  
  time_raw <- ncvar_get(nc_sim, "time")
  origin_time <- as.POSIXct("2015-01-01 00:00:00", tz = "UTC")
  time <- origin_time + time_raw * 60
  weights <- days_in_month(time) / sum(days_in_month(time))
  
  weighted_avg <- function(x) sum(x * weights)
  dep_Hg2_mean <- apply(dep_Hg2, MARGIN = c(1, 2), FUN = weighted_avg)
  dep_HgP_mean <- apply(dep_HgP, MARGIN = c(1, 2), FUN = weighted_avg)
  dep_Hg0_mean <- apply(dep_Hg0, MARGIN = c(1, 2), FUN = weighted_avg)
  
  dep_total_mg_contrib <- (dep_Hg2_mean + dep_HgP_mean + dep_Hg0_mean) * unit_conv_dd * kg_Mg
  
  results_amz[[sim]] <- dep_total_mg_contrib
  
  nc_close(nc_sim)
}

V03A03 <- results_amz[["V03A03"]]
V03A20 <- results_amz[["V03A20"]]
V20A03 <- results_amz[["V20A03"]]
V20A20 <- results_amz[["V20A20"]]

# Analyse croisée
Eff_A <- 0.5 * ((V03A20 - V03A03) + (V20A20 - V20A03))
Eff_V <- 0.5 * ((V20A03 - V03A03) + (V20A20 - V03A20))
Interaction <- V20A20 - V03A03 - Eff_A - Eff_V

# Moyennes globales
mean_A <- mean(Eff_A, na.rm = TRUE)
mean_V <- mean(Eff_V, na.rm = TRUE)
mean_I <- mean(Interaction, na.rm = TRUE)

total_A <- mean_A + 0.5 * mean_I
total_V <- mean_V + 0.5 * mean_I
total_sum <- total_A + total_V

contrib <- c(
  Emissions = 100 * total_A / total_sum,
  Végétation = 100 * total_V / total_sum
)

df <- data.frame(
  source = names(contrib),
  pct = contrib
)

df$label <- paste0(df$source, "\n", sprintf("%.1f", df$pct), "%")

contrib_conc <- ggplot(df, aes(x = "", y = pct, fill = source)) +
  geom_col(width = 1, color = "white") +
  geom_text(aes(label = label), 
            position = position_stack(vjust = 0.5), 
            color = "white", size = 5, fontface = "bold") +
  coord_polar("y") +
  scale_fill_manual(values = c("#D95F02", "#1B9E77")) +
  labs(
    title = "Relative contribution to the change in Hg total drydep",
    subtitle = "Dataset : EDGAR / Global",
    fill = "Source"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    plot.background = element_rect(fill = "white", color = NA)
  )

print(contrib_conc)
```

![](EDGAR_Map_Dep_05_08_2025_files/figure-gfm/contrib%20drydep%20world-1.png)<!-- -->

``` r
chemin_pdf <- "C:/Users/colom/Desktop/STAGE/data/clean_mod_data/plot/EDGAR_Contrib_Drydep.pdf"

# Enregistrement du plot
ggsave(
  filename = chemin_pdf,
  plot = contrib_conc,
  device = "pdf",
  width = 10, height = 6, units = "in"
)
```

``` r
# --- Simulations ---
simulations <- c("V03A03", "V03A20", "V20A03", "V20A20")

results_amz <- list()

for (sim in simulations) {
  
  # 1. Ouvrir le fichier NetCDF
  path_sim <- file.path(chem_EDGAR, sim, "GEOSChem.DryDep.2015_m.nc4")
  nc_sim <- nc_open(path_sim)
  
  # 2. Lire les variables de dépôts
  dep_Hg2 <- ncvar_get(nc_sim, "DryDep_Hg2")
  dep_HgP <- ncvar_get(nc_sim, "DryDep_HgP")
  dep_Hg0 <- ncvar_get(nc_sim, "DryDep_Hg0")
  
  # 3. Temps et pondération par mois
  time_raw <- ncvar_get(nc_sim, "time")
  origin_time <- as.POSIXct("2015-01-01 00:00:00", tz = "UTC")
  time <- origin_time + time_raw * 60
  weights <- days_in_month(time) / sum(days_in_month(time))
  
  weighted_avg <- function(x) sum(x * weights)
  
  dep_Hg2_mean <- apply(dep_Hg2, MARGIN = c(1,2), FUN = weighted_avg)
  dep_HgP_mean <- apply(dep_HgP, MARGIN = c(1,2), FUN = weighted_avg)
  dep_Hg0_mean <- apply(dep_Hg0, MARGIN = c(1,2), FUN = weighted_avg)
  
  # 4. Conversion en mg/m²/jour (ou autre unité déjà définie)
  dep_total_mg <- (dep_Hg2_mean + dep_HgP_mean + dep_Hg0_mean) * unit_conv_dd * kg_Mg
  
  # 5. Appliquer le masque Amazonie
  dep_amz <- dep_total_mg * Am_mask
  
  # 6. Stocker le résultat
  results_amz[[sim]] <- dep_amz
  
  # 7. Fermer le NetCDF
  nc_close(nc_sim)
}

# --- Extraire les simulations ---
V03A03 <- results_amz[["V03A03"]]
V03A20 <- results_amz[["V03A20"]]
V20A03 <- results_amz[["V20A03"]]
V20A20 <- results_amz[["V20A20"]]

# --- Décomposition ANOVA ---
Eff_A <- 0.5 * ((V03A20 - V03A03) + (V20A20 - V20A03))
Eff_V <- 0.5 * ((V20A03 - V03A03) + (V20A20 - V03A20))
Interaction <- V20A20 - V03A03 - Eff_A - Eff_V

# --- Extraire l’Amazonie ---
Eff_A_amz <- Eff_A[Am_mask == 1]
Eff_V_amz <- Eff_V[Am_mask == 1]
Interaction_amz <- Interaction[Am_mask == 1]

# --- Moyennes régionales ---
mean_A <- mean(Eff_A_amz, na.rm = TRUE)
mean_V <- mean(Eff_V_amz, na.rm = TRUE)
mean_I <- mean(Interaction_amz, na.rm = TRUE)

# --- Camembert Emissions / Végétation ---
total_A <- mean_A + 0.5 * mean_I
total_V <- mean_V + 0.5 * mean_I
total_sum <- total_A + total_V

contrib <- c(
  Emissions  = 100 * total_A / total_sum,
  Vegetation = 100 * total_V / total_sum
)

df <- data.frame(
  source = names(contrib),
  pct = contrib
)
df$label <- paste0(df$source, "\n", sprintf("%.1f", df$pct), "%")

contrib_deps <- ggplot(df, aes(x = "", y = pct, fill = source)) +
  geom_col(width = 1, color = "white") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5),
            color = "white", size = 5, fontface = "bold") +
  coord_polar("y") +
  scale_fill_manual(values = c("#D95F02", "#1B9E77")) +
  labs(title = "Relative contribution to the change in Hg total drydep",
       subtitle = "Dataset : EDGAR / Amazonia forest focus",
       fill = "Source") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    plot.background = element_rect(fill = "white", color = NA)
  )

print(contrib_deps)
```

![](EDGAR_Map_Dep_05_08_2025_files/figure-gfm/test%20amazo%20drydep-1.png)<!-- -->

``` r
# --- Sauvegarde ---
chemin_pdf <- "C:/Users/colom/Desktop/STAGE/data/clean_mod_data/plot/TEST_EDGAR_Contrib_Drydep_Amazo.pdf"
ggsave(filename = chemin_pdf, plot = contrib_deps, device = "pdf", width = 10, height = 6, units = "in")
```
