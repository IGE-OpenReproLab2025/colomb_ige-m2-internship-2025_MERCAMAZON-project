Concentrations modélisées de Hg⁰ vs données de référence
================
Martin Colomb
2025-07-04

# Libraries

``` r
library(dplyr)
library(ggplot2)
library(lubridate)
library(ncdf4)
library(viridis)
library(maps)
library(geosphere)
library(patchwork)
```

# Chemins

``` r
chem_martin <- "C:/Users/colom/Desktop/STAGE/data/clean_mod_data/14_3_1/HIST_V3"
chem_feinberg <- "C:/Users/colom/Desktop/STAGE/data/clean_mod_data/14_3_1/HIST_Feinberg"
chem_ref <- "C:/Users/colom/Desktop/STAGE/data/clean_mod_data/14_3_1/ref_wd_con_hg"
```

# Constantes

``` r
M_Hg <- 200.59
R <- 8.3145
P <- 101325
T <- 298.15
```

# Poids mois

``` r
nc_tmp <- nc_open(file.path(chem_martin, "GEOSChem.SpeciesConc.2015_m.nc4"))
time_raw <- ncvar_get(nc_tmp, "time")
nc_close(nc_tmp)
origin_time <- as.POSIXct("2015-01-01 00:00:00", tz = "UTC")
time <- origin_time + time_raw * 60
days_in_months <- days_in_month(time)
weights <- days_in_months / sum(days_in_months)
```

# Fonctions

``` r
load_hg0_conc <- function(filepath, varname, weights) {
  nc <- nc_open(filepath)
  conc <- ncvar_get(nc, varname)[,,1,]
  nc_close(nc)
  conc <- apply(conc, MARGIN = c(1, 2), FUN = function(x) sum(x * weights))
  conc * P * (M_Hg / (R * T)) * 1e9
}

create_df_conc <- function(matrix_conc, lon, lat) {
  expand.grid(lon = lon, lat = lat) %>%
    mutate(conc = as.vector(matrix_conc))
}

plot_concentration_map <- function(df, title, ref_data = NULL, limits = c(0, 5)) {
  p <- ggplot(df, aes(x = lon, y = lat, fill = conc)) +
    geom_raster() +
    scale_fill_viridis_c(option = "magma", name = "Hg⁰ (ng/m³)", direction = -1, limits = limits) +
    coord_fixed() +
    labs(title = title, x = "Longitude", y = "Latitude") +
    borders("world", colour = "black") +
    theme_minimal()
  if (!is.null(ref_data)) {
    p <- p + geom_point(data = ref_data, aes(x = Lon, y = Lat, fill = Hg0),
                        size = 3, shape = 21, inherit.aes = FALSE)
  }
  return(p)
}
```

# Grilles géo

``` r
lon <- seq(-180, 180, length.out = 144)
lat <- seq(-90, 90, length.out = 91)
```

# Plot concentrations Martin

``` r
conc_martin <- load_hg0_conc(file.path(chem_martin, "GEOSChem.SpeciesConc.2015_m.nc4"), "SpeciesConcVV_Hg0", weights)
df_martin <- create_df_conc(conc_martin, lon, lat)
plot_martin <- plot_concentration_map(df_martin, "Concentration de Hg⁰ modélisée (Martin) pour l'année 2015", limits = c(0, 4))
plot_martin
```

![](Map_mod_vs_ref_04_07_2025_files/figure-gfm/concentration_martin-1.png)<!-- -->

# Plot concentrations Feinberg

``` r
conc_feinberg <- load_hg0_conc(file.path(chem_feinberg, "Feinberg_GEOSChem.SpeciesConc.2015_m.nc4"), "SpeciesConc_Hg0", weights)
df_feinberg <- create_df_conc(conc_feinberg, lon, lat)
plot_feinberg <- plot_concentration_map(df_feinberg, "Concentration de Hg⁰ modélisée (Feinberg) pour l'année 2015", limits = c(0, 5))
plot_feinberg
```

![](Map_mod_vs_ref_04_07_2025_files/figure-gfm/concentration_feinberg-1.png)<!-- -->

# Ajout données ref

``` r
ref_hg0 <- read.csv2(file.path(chem_ref, "Hg0_annual_2013-2015_Apr13_2024.csv"), sep = ",", header = TRUE)
colnames(ref_hg0) <- as.character(unlist(ref_hg0[1, ]))
ref_hg0 <- ref_hg0[-1, ] %>%
  mutate(Lat = as.numeric(Lat),
         Lon = as.numeric(Lon),
         Hg0 = as.numeric(Hg0))
```

# Plot concentrations Martin vs ref

``` r
plot_martin_ref <- plot_concentration_map(df_martin,
  "Concentration de Hg⁰ (Martin) vs références pour l'année 2015",
  ref_data = ref_hg0, limits = c(0, 4))
plot_martin_ref
```

![](Map_mod_vs_ref_04_07_2025_files/figure-gfm/martin_vs_ref-1.png)<!-- -->

# Plot concentrations Feinberg vs ref

``` r
plot_feinberg_ref <- plot_concentration_map(df_feinberg,
  "Concentration de Hg⁰ (Feinberg) vs références pour l'année 2015",
  ref_data = ref_hg0, limits = c(0, 5))
plot_feinberg_ref
```

![](Map_mod_vs_ref_04_07_2025_files/figure-gfm/feinberg_vs_ref-1.png)<!-- -->
