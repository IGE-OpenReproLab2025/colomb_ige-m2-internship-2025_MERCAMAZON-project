Map diff conc EDGAR 2003 2020
================
Martin Colomb
2025-08-26

# Chargement des packages

``` r
library(dplyr)
library(ggplot2)
library(lubridate)
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
# Paramètres généraux
simulations <- c("V03A03", "V03A20", "V20A03", "V20A20")
lon <- seq(-180, 180, length.out = 144)
lat <- seq(-90, 90, length.out = 91)
df_grid <- expand.grid(x = lon, y = lat)

# Classes et couleurs pour l'affichage
breaks_conc <- c(0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, Inf)
labels_chr <- c(
  "0.4–0.6", "0.6–0.8", "0.8–1.0", "1.0–1.2",
  "1.2–1.4", "1.4–1.6", "1.6–1.8", "1.8–2.0", ">2.0"
)
labels_expr <- c(
  expression("0.4–0.6"), expression("0.6–0.8"), expression("0.8–1.0"),
  expression("1.0–1.2"), expression("1.2–1.4"), expression("1.4–1.6"),
  expression("1.6–1.8"), expression("1.8–2.0"), expression("">"~2.0")
)
colors_conc <- c(
  "#005F73", "#0A9396", "#94D2BD", "#E9D8A6",
  "#EE9B00", "#CA6702", "#BB3E03", "#AE2012", "#9B2226"
)


plot_simulation_map <- function(conc_matrix, sim_name, save_path) {
  df <- df_grid
  df$conc_ng_m3 <- as.vector(conc_matrix)
  df$conc_class <- cut(df$conc_ng_m3,
                       breaks = breaks_conc,
                       labels = labels_chr,
                       include.lowest = TRUE)
  df <- na.omit(df)
  
  plot <- ggplot(df, aes(x = x, y = y, fill = conc_class)) +
    geom_raster() +
    scale_fill_manual(
      name = "Hg(0) (ng/m³)",
      values = colors_conc,
      labels = labels_expr,
      na.value = "transparent",
      guide = guide_legend(reverse = TRUE)
    ) +
    coord_equal() +
    theme_minimal() +
    labs(
      title = "Annual Mean Atmospheric Hg(0) concentration in 2015",
      subtitle = paste0("Simulation : ", sim_name),
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
  
  ggsave(
    filename = save_path,
    plot = plot,
    device = "pdf",
    width = 10, height = 6, units = "in"
  )
}

# Boucle principale
for (sim in simulations) {
  
  # 1. Ouvrir le fichier
  path_sim <- file.path(chem_EDGAR, sim, "GEOSChem.SpeciesConc.2015_m.nc4")
  nc_sim <- nc_open(path_sim)
  
  # 2. Extraire la concentration Hg0 au niveau 1
  hg0_raw <- ncvar_get(nc_sim, "SpeciesConcVV_Hg0")[, , 1, ]
  
  # 3. Temps et poids mensuels
  time_raw <- ncvar_get(nc_sim, "time")
  origin_time <- as.POSIXct("2015-01-01 00:00:00", tz = "UTC")
  time <- origin_time + time_raw * 60
  weights <- days_in_month(time) / sum(days_in_month(time))
  
  # 4. Moyenne pondérée temporelle
  weighted_avg <- function(x) sum(x * weights)
  hg0_mean <- apply(hg0_raw, MARGIN = c(1, 2), FUN = weighted_avg)
  
  # 5. Conversion en ng/m³
  hg0_ng_m3 <- hg0_mean * (M_Hg / M_air) * mean_Met_AIRDEN * ng_per_kg
  
  # 6. Sauvegarde du plot
  save_path <- file.path(
    "C:/Users/colom/Desktop/STAGE/data/clean_mod_data/plot",
    paste0("EDGAR_", sim, "_Conc_WORLD.pdf")
  )
  
  plot_simulation_map(hg0_ng_m3, sim, save_path)
  
  # 7. Fermer NetCDF
  nc_close(nc_sim)
}
```

``` r
conc_V03A03 <- nc_open(file.path(chem_EDGAR, "V03A03/GEOSChem.SpeciesConc.2015_m.nc4"))
print(conc_V03A03)
```

    ## File C:/Users/colom/Desktop/STAGE/data/clean_mod_data/14_3_1/EDGAR/V03A03/GEOSChem.SpeciesConc.2015_m.nc4 (NC_FORMAT_NETCDF4_CLASSIC):
    ## 
    ##      49 variables (excluding dimension variables):
    ##         double time_bnds[bnds,time]   (Chunking: [2,1])  
    ##         double lon_bnds[bnds,lon]   (Contiguous storage)  
    ##         double lat_bnds[bnds,lat]   (Contiguous storage)  
    ##         double lev_bnds[bnds,lev]   (Contiguous storage)  
    ##             standard_name: atmosphere_hybrid_sigma_pressure_coordinate
    ##             units: level
    ##             formula_terms: ap: ap_bnds b: b_bnds ps: ps
    ##         float AREA[lon,lat]   (Chunking: [144,91])  
    ##             long_name: Surface area
    ##             units: m2
    ##             cell_methods: time: mean
    ##         float SpeciesConcVV_CH4[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species CH4
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_CO[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species CO
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_ClO[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species ClO
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_BrO[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species BrO
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HO2[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HO2
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_O3[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species O3
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_NO[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species NO
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_NO2[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species NO2
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_OH[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species OH
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_Cl[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species Cl
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_Br[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species Br
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgClO[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgClO
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgBrO[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgBrO
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgOHO[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgOHO
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_PHg2Cl[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species PHg2Cl
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_PHg2OH[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species PHg2OH
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_PHg2Br[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species PHg2Br
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_PHg0[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species PHg0
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_PHg2[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species PHg2
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_Hg2STRP[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species Hg2STRP
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_Hg2ORGP[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species Hg2ORGP
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_Hg2ClP[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species Hg2ClP
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgCl2[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgCl2
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgOHOH[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgOHOH
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgOHBrO[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgOHBrO
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgOHClO[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgOHClO
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgOHHO2[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgOHHO2
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgOHNO2[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgOHNO2
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgOH[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgOH
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgClOH[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgClOH
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgClBr[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgClBr
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgClBrO[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgClBrO
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgClClO[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgClClO
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgClHO2[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgClHO2
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgClNO2[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgClNO2
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgCl[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgCl
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgBr2[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgBr2
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgBrOH[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgBrOH
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgBrClO[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgBrClO
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgBrBrO[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgBrBrO
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgBrHO2[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgBrHO2
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgBrNO2[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgBrNO2
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgBr[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgBr
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_Hg0[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species Hg0
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ## 
    ##      5 dimensions:
    ##         time  Size:12   *** is unlimited *** 
    ##             standard_name: time
    ##             long_name: Time
    ##             bounds: time_bnds
    ##             units: minutes since 2015-01-01 00:00:00
    ##             calendar: gregorian
    ##             axis: T
    ##         bnds  Size:2 (no dimvar)
    ##         lon  Size:144 
    ##             standard_name: longitude
    ##             long_name: Longitude
    ##             units: degrees_east
    ##             axis: X
    ##             bounds: lon_bnds
    ##         lat  Size:91 
    ##             standard_name: latitude
    ##             long_name: Latitude
    ##             units: degrees_north
    ##             axis: Y
    ##             bounds: lat_bnds
    ##         lev  Size:47 
    ##             standard_name: atmosphere_hybrid_sigma_pressure_coordinate
    ##             axis: Z
    ##             positive: down
    ##             long_name: hybrid sigma pressure coordinate
    ##             units: level
    ##             formula_terms: ap: ap b: b ps: ps
    ##             bounds: lev_bnds
    ## 
    ##     13 global attributes:
    ##         CDI: Climate Data Interface version 2.1.1 (https://mpimet.mpg.de/cdi)
    ##         Conventions: CF-1.6
    ##         title: GEOS-Chem diagnostic collection: SpeciesConc
    ##         history: Fri Aug 01 10:02:09 2025: cdo monmean GEOSChem.SpeciesConc.2015.nc4 GEOSChem.SpeciesConc.2015_m.nc4
    ## Fri Aug 01 10:01:26 2025: cdo mergetime GEOSChem.SpeciesConc.20150101_0000z.nc4 GEOSChem.SpeciesConc.20150201_0000z.nc4 GEOSChem.SpeciesConc.20150301_0000z.nc4 GEOSChem.SpeciesConc.20150401_0000z.nc4 GEOSChem.SpeciesConc.20150501_0000z.nc4 GEOSChem.SpeciesConc.20150601_0000z.nc4 GEOSChem.SpeciesConc.20150701_0000z.nc4 GEOSChem.SpeciesConc.20150801_0000z.nc4 GEOSChem.SpeciesConc.20150901_0000z.nc4 GEOSChem.SpeciesConc.20151001_0000z.nc4 GEOSChem.SpeciesConc.20151101_0000z.nc4 GEOSChem.SpeciesConc.20151201_0000z.nc4 GEOSChem.SpeciesConc.2015.nc4
    ##         format: NetCDF-4
    ##         conventions: COARDS
    ##         ProdDateTime: 
    ##         reference: www.geos-chem.org; wiki.geos-chem.org
    ##         contact: GEOS-Chem Support Team (geos-chem-support@g.harvard.edu)
    ##         simulation_start_date_and_time: 2015-01-01 00:00:00z
    ##         simulation_end_date_and_time: 2016-01-01 00:00:00z
    ##         frequency: mon
    ##         CDO: Climate Data Operators version 2.1.1 (https://mpimet.mpg.de/cdo)

``` r
hg0_V03A03 <- ncvar_get(conc_V03A03, "SpeciesConcVV_Hg0") #mol mol-1 dry

dim(hg0_V03A03) 
```

    ## [1] 144  91  47  12

``` r
hg0_V03A03 <- hg0_V03A03[, , 1, ]  # Prend le niveau 1, garde lon, lat, time

time_raw <- ncvar_get(conc_V03A03, "time")
origin_time <- as.POSIXct("2015-01-01 00:00:00", tz = "UTC")
time <- origin_time + time_raw * 60
days_in_months <- days_in_month(time)
weights <- days_in_months / sum(days_in_months)
weighted_avg <- function(x) sum(x * weights)
hg0_V03A03 <- apply(hg0_V03A03, MARGIN = c(1, 2), FUN = weighted_avg)
summary(as.vector(hg0_V03A03))
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ## 5.604e-14 9.660e-14 1.011e-13 1.016e-13 1.060e-13 3.003e-13

``` r
dim(hg0_V03A03)#mol mol-1 dry
```

    ## [1] 144  91

``` r
dim(mean_Met_AIRDEN)#Dry air density units: kg m-3
```

    ## [1] 144  91

``` r
M_Hg         # Masse molaire du Hg (kg/mol)
```

    ## [1] 0.20059

``` r
M_air        # Masse molaire de l'air sec (kg/mol)
```

    ## [1] 0.028964

``` r
ng_per_kg    #kg -> ng (1e12)
```

    ## [1] 1e+12

``` r
#Formule : hg0_V03A03 * (M_Hg/M_air) * mean_Met_AIRDEN  = 

V03A03_ng_m3 <-   hg0_V03A03 * (M_Hg/M_air) * mean_Met_AIRDEN * ng_per_kg
range(V03A03_ng_m3)
```

    ## [1] 0.4419279 2.2120839

``` r
dim(V03A03_ng_m3)
```

    ## [1] 144  91

``` r
lon <- seq(-180, 180, length.out = 144)
lat <- seq(-90, 90, length.out = 91)
df_conc <- expand.grid(x = lon, y = lat)
df_conc$conc_ng_m3 <- as.vector(V03A03_ng_m3)
breaks_conc <- c(0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, Inf)
labels_chr <- c(
  "0.4–0.6", "0.6–0.8", "0.8–1.0", "1.0–1.2",
  "1.2–1.4", "1.4–1.6", "1.6–1.8", "1.8–2.0", ">2.0"
)

labels_expr <- c(
  expression("0.4–0.6"), expression("0.6–0.8"), expression("0.8–1.0"),
  expression("1.0–1.2"), expression("1.2–1.4"), expression("1.4–1.6"),
  expression("1.6–1.8"), expression("1.8–2.0"), expression("">"~2.0")
)

colors_conc <- c(
  "#005F73", "#0A9396", "#94D2BD", "#E9D8A6",
  "#EE9B00", "#CA6702", "#BB3E03", "#AE2012", "#9B2226"
)

# Coupe selon les classes
df_conc$conc_class <- cut(
  df_conc$conc_ng_m3,
  breaks = breaks_conc,
  labels = labels_chr,
  include.lowest = TRUE
)

# Supprimer les NA
df_conc <- na.omit(df_conc)
plot_conc <- ggplot(df_conc, aes(x = x, y = y, fill = conc_class)) +
  geom_raster() +
  scale_fill_manual(
    name = "Hg(0) (ng/m³)",
    values = colors_conc,
    labels = labels_expr,
    na.value = "transparent",
    guide = guide_legend(reverse = TRUE)
  ) +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Annual Mean Atmospheric Hg(0) concentration in 2015",
    subtitle = "Simulation : Vegetation 2003 / Emissions EDGAR 2003",
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

plot_conc
```

![](EDGAR_Map_Diff_Conc_01_08_2025_files/figure-gfm/Map%20concentration%20Hg0%202003%20V03A03-1.png)<!-- -->

``` r
# Chemin de sauvegarde (modifie si nécessaire)
chemin_pdf <- "C:/Users/colom/Desktop/STAGE/data/clean_mod_data/plot/EDGAR_V03A03_Conc_WORLD.pdf"

# Enregistrement du plot
ggsave(
  filename = chemin_pdf,
  plot = plot_conc,
  device = "pdf",
  width = 10, height = 6, units = "in"
)
```

``` r
conc_V20A20 <- nc_open(file.path(chem_EDGAR, "V20A20/GEOSChem.SpeciesConc.2015_m.nc4"))
print(conc_V20A20)
```

    ## File C:/Users/colom/Desktop/STAGE/data/clean_mod_data/14_3_1/EDGAR/V20A20/GEOSChem.SpeciesConc.2015_m.nc4 (NC_FORMAT_NETCDF4_CLASSIC):
    ## 
    ##      49 variables (excluding dimension variables):
    ##         double time_bnds[bnds,time]   (Chunking: [2,1])  
    ##         double lon_bnds[bnds,lon]   (Contiguous storage)  
    ##         double lat_bnds[bnds,lat]   (Contiguous storage)  
    ##         double lev_bnds[bnds,lev]   (Contiguous storage)  
    ##             standard_name: atmosphere_hybrid_sigma_pressure_coordinate
    ##             units: level
    ##             formula_terms: ap: ap_bnds b: b_bnds ps: ps
    ##         float AREA[lon,lat]   (Chunking: [144,91])  
    ##             long_name: Surface area
    ##             units: m2
    ##             cell_methods: time: mean
    ##         float SpeciesConcVV_CH4[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species CH4
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_CO[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species CO
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_ClO[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species ClO
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_BrO[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species BrO
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HO2[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HO2
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_O3[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species O3
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_NO[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species NO
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_NO2[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species NO2
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_OH[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species OH
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_Cl[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species Cl
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_Br[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species Br
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgClO[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgClO
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgBrO[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgBrO
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgOHO[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgOHO
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_PHg2Cl[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species PHg2Cl
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_PHg2OH[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species PHg2OH
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_PHg2Br[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species PHg2Br
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_PHg0[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species PHg0
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_PHg2[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species PHg2
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_Hg2STRP[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species Hg2STRP
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_Hg2ORGP[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species Hg2ORGP
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_Hg2ClP[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species Hg2ClP
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgCl2[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgCl2
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgOHOH[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgOHOH
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgOHBrO[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgOHBrO
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgOHClO[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgOHClO
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgOHHO2[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgOHHO2
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgOHNO2[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgOHNO2
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgOH[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgOH
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgClOH[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgClOH
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgClBr[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgClBr
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgClBrO[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgClBrO
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgClClO[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgClClO
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgClHO2[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgClHO2
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgClNO2[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgClNO2
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgCl[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgCl
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgBr2[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgBr2
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgBrOH[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgBrOH
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgBrClO[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgBrClO
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgBrBrO[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgBrBrO
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgBrHO2[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgBrHO2
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgBrNO2[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgBrNO2
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_HgBr[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species HgBr
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ##         float SpeciesConcVV_Hg0[lon,lat,lev,time]   (Chunking: [144,91,1,1])  
    ##             long_name: Concentration of species Hg0
    ##             units: mol mol-1 dry
    ##             cell_methods: time: mean
    ##             averaging_method: time-averaged
    ## 
    ##      5 dimensions:
    ##         time  Size:12   *** is unlimited *** 
    ##             standard_name: time
    ##             long_name: Time
    ##             bounds: time_bnds
    ##             units: minutes since 2015-01-01 00:00:00
    ##             calendar: gregorian
    ##             axis: T
    ##         bnds  Size:2 (no dimvar)
    ##         lon  Size:144 
    ##             standard_name: longitude
    ##             long_name: Longitude
    ##             units: degrees_east
    ##             axis: X
    ##             bounds: lon_bnds
    ##         lat  Size:91 
    ##             standard_name: latitude
    ##             long_name: Latitude
    ##             units: degrees_north
    ##             axis: Y
    ##             bounds: lat_bnds
    ##         lev  Size:47 
    ##             standard_name: atmosphere_hybrid_sigma_pressure_coordinate
    ##             axis: Z
    ##             positive: down
    ##             long_name: hybrid sigma pressure coordinate
    ##             units: level
    ##             formula_terms: ap: ap b: b ps: ps
    ##             bounds: lev_bnds
    ## 
    ##     13 global attributes:
    ##         CDI: Climate Data Interface version 2.1.1 (https://mpimet.mpg.de/cdi)
    ##         Conventions: CF-1.6
    ##         title: GEOS-Chem diagnostic collection: SpeciesConc
    ##         history: Fri Aug 01 10:11:21 2025: cdo monmean GEOSChem.SpeciesConc.2015.nc4 GEOSChem.SpeciesConc.2015_m.nc4
    ## Fri Aug 01 10:10:56 2025: cdo mergetime GEOSChem.SpeciesConc.20150101_0000z.nc4 GEOSChem.SpeciesConc.20150201_0000z.nc4 GEOSChem.SpeciesConc.20150301_0000z.nc4 GEOSChem.SpeciesConc.20150401_0000z.nc4 GEOSChem.SpeciesConc.20150501_0000z.nc4 GEOSChem.SpeciesConc.20150601_0000z.nc4 GEOSChem.SpeciesConc.20150701_0000z.nc4 GEOSChem.SpeciesConc.20150801_0000z.nc4 GEOSChem.SpeciesConc.20150901_0000z.nc4 GEOSChem.SpeciesConc.20151001_0000z.nc4 GEOSChem.SpeciesConc.20151101_0000z.nc4 GEOSChem.SpeciesConc.20151201_0000z.nc4 GEOSChem.SpeciesConc.2015.nc4
    ##         format: NetCDF-4
    ##         conventions: COARDS
    ##         ProdDateTime: 
    ##         reference: www.geos-chem.org; wiki.geos-chem.org
    ##         contact: GEOS-Chem Support Team (geos-chem-support@g.harvard.edu)
    ##         simulation_start_date_and_time: 2015-01-01 00:00:00z
    ##         simulation_end_date_and_time: 2016-01-01 00:00:00z
    ##         frequency: mon
    ##         CDO: Climate Data Operators version 2.1.1 (https://mpimet.mpg.de/cdo)

``` r
hg0_V20A20 <- ncvar_get(conc_V20A20, "SpeciesConcVV_Hg0") #mol mol-1 dry

dim(hg0_V20A20) 
```

    ## [1] 144  91  47  12

``` r
hg0_V20A20 <- hg0_V20A20[, , 1, ]  # Prend le niveau 1, garde lon, lat, time

time_raw <- ncvar_get(conc_V20A20, "time")
origin_time <- as.POSIXct("2015-01-01 00:00:00", tz = "UTC")
time <- origin_time + time_raw * 60
days_in_months <- days_in_month(time)
weights <- days_in_months / sum(days_in_months)
weighted_avg <- function(x) sum(x * weights)
hg0_V20A20 <- apply(hg0_V20A20, MARGIN = c(1, 2), FUN = weighted_avg)
summary(as.vector(hg0_V20A20))
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ## 6.648e-14 1.016e-13 1.068e-13 1.079e-13 1.135e-13 3.650e-13

``` r
dim(hg0_V20A20)#mol mol-1 dry
```

    ## [1] 144  91

``` r
dim(mean_Met_AIRDEN)#Dry air density units: kg m-3
```

    ## [1] 144  91

``` r
M_Hg         # Masse molaire du Hg (kg/mol)
```

    ## [1] 0.20059

``` r
M_air        # Masse molaire de l'air sec (kg/mol)
```

    ## [1] 0.028964

``` r
ng_per_kg    #kg -> ng (1e12)
```

    ## [1] 1e+12

``` r
#Formule : hg0_V20A20 * (M_Hg/M_air) * mean_Met_AIRDEN  = 

V20A20_ng_m3 <-   hg0_V20A20 * (M_Hg/M_air) * mean_Met_AIRDEN * ng_per_kg
range(V20A20_ng_m3)
```

    ## [1] 0.5242534 2.5489646

``` r
dim(V20A20_ng_m3)
```

    ## [1] 144  91

``` r
lon <- seq(-180, 180, length.out = 144)
lat <- seq(-90, 90, length.out = 91)
df_conc <- expand.grid(x = lon, y = lat)
df_conc$conc_ng_m3 <- as.vector(V20A20_ng_m3)
breaks_conc <- c(0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, Inf)
labels_chr <- c(
  "0.4–0.6", "0.6–0.8", "0.8–1.0", "1.0–1.2",
  "1.2–1.4", "1.4–1.6", "1.6–1.8", "1.8–2.0", ">2.0"
)

labels_expr <- c(
  expression("0.4–0.6"), expression("0.6–0.8"), expression("0.8–1.0"),
  expression("1.0–1.2"), expression("1.2–1.4"), expression("1.4–1.6"),
  expression("1.6–1.8"), expression("1.8–2.0"), expression("">"~2.0")
)

colors_conc <- c(
  "#005F73", "#0A9396", "#94D2BD", "#E9D8A6",
  "#EE9B00", "#CA6702", "#BB3E03", "#AE2012", "#9B2226"
)

# Coupe selon les classes
df_conc$conc_class <- cut(
  df_conc$conc_ng_m3,
  breaks = breaks_conc,
  labels = labels_chr,
  include.lowest = TRUE
)

# Supprimer les NA
df_conc <- na.omit(df_conc)
plot_conc <- ggplot(df_conc, aes(x = x, y = y, fill = conc_class)) +
  geom_raster() +
  scale_fill_manual(
    name = "Hg(0) (ng/m³)",
    values = colors_conc,
    labels = labels_expr,
    na.value = "transparent",
    guide = guide_legend(reverse = TRUE)
  ) +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Annual Mean Atmospheric Hg(0) concentration in 2015",
    subtitle = "Simulation : Vegetation 2020 / Emissions EDGAR 2020",
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

plot_conc
```

![](EDGAR_Map_Diff_Conc_01_08_2025_files/figure-gfm/Map%20concentration%202020-1.png)<!-- -->

``` r
# Chemin de sauvegarde (modifie si nécessaire)
chemin_pdf <- "C:/Users/colom/Desktop/STAGE/data/clean_mod_data/plot/EDGAR_V20A20_Conc_WORLD.pdf"

# Enregistrement du plot
ggsave(
  filename = chemin_pdf,
  plot = plot_conc,
  device = "pdf",
  width = 10, height = 6, units = "in"
)
```

``` r
range(V20A20_ng_m3)
```

    ## [1] 0.5242534 2.5489646

``` r
range(V03A03_ng_m3)
```

    ## [1] 0.4419279 2.2120839

``` r
diff03_20<- V20A20_ng_m3 - V03A03_ng_m3 #ng/m3

range(diff03_20)
```

    ## [1] -0.7686564  1.1381597

``` r
lon <- seq(-180, 180, length.out = 144)
lat <- seq(-90, 90, length.out = 91)
df_diff <- expand.grid(x = lon, y = lat)
df_diff$diff_conc <- as.vector(diff03_20)

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
    title = "Différence annuelle moyenne Hg(0) atmosphérique (2020 - 2003)",
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

![](EDGAR_Map_Diff_Conc_01_08_2025_files/figure-gfm/diff%202003%202020%20WORLD-1.png)<!-- -->

``` r
lim <- 0.5
plot_diff_zoom <- ggplot(df_diff, aes(x = x, y = y, fill = diff_conc)) +
  geom_raster() +
  scale_fill_gradient2(
    name = "Hg(0) difference (ng/m³)\n2020 - 2003",
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
    title = paste0("Average annual difference of Hg(0) concentration between 2020 and 2003"),
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

![](EDGAR_Map_Diff_Conc_01_08_2025_files/figure-gfm/diff%202003%202020%20WORLD-2.png)<!-- -->

``` r
# Chemin de sauvegarde (modifie si nécessaire)
chemin_pdf <- "C:/Users/colom/Desktop/STAGE/data/clean_mod_data/plot/EDGAR_diff_20_03_Conc_WORLD.pdf"

# Enregistrement du plot
ggsave(
  filename = chemin_pdf,
  plot = plot_diff_zoom,
  device = "pdf",
  width = 10, height = 6, units = "in"
)
```

<!-- ```{r diff EMI 2003 2020 WORLD} -->
<!-- range(V20A20_ng_m3) -->
<!-- range(V03A03_ng_m3) -->
<!-- diff03_20<- V20A20_ng_m3 - V03A03_ng_m3 #ng/m3 -->
<!-- range(diff03_20) -->
<!-- lon <- seq(-180, 180, length.out = 144) -->
<!-- lat <- seq(-90, 90, length.out = 91) -->
<!-- df_diff <- expand.grid(x = lon, y = lat) -->
<!-- df_diff$diff_conc <- as.vector(diff03_20) -->
<!-- df_diff <- na.omit(df_diff) -->
<!-- colors_div <- c("#2166AC", "#FFFFFF", "#B2182B")  # bleu - blanc - rouge -->
<!-- plot_diff <- ggplot(df_diff, aes(x = x, y = y, fill = diff_conc)) + -->
<!--   geom_raster() + -->
<!--   scale_fill_gradient2( -->
<!--     name = "Différence Hg(0) (ng/m³)\n2020 - 2003", -->
<!--     low = colors_div[1], -->
<!--     mid = colors_div[2], -->
<!--     high = colors_div[3], -->
<!--     midpoint = 0, -->
<!--     limits = c(-max(abs(df_diff$diff_conc)), max(abs(df_diff$diff_conc))), -->
<!--     oob = scales::squish -->
<!--   ) + -->
<!--   coord_equal() + -->
<!--   theme_minimal() + -->
<!--   labs( -->
<!--     title = "Différence annuelle moyenne Hg(0) atmosphérique (2020 - 2003)", -->
<!--     x = "Longitude", y = "Latitude" -->
<!--   ) + -->
<!--   borders("world", colour = "#403d39", size = 0.2) + -->
<!--   theme( -->
<!--     plot.title = element_text(hjust = 0.5, face = "bold"), -->
<!--     axis.text = element_text(size = 10), -->
<!--     axis.ticks = element_line(color = "black", size = 0.3), -->
<!--     axis.ticks.length = unit(3, "pt"), -->
<!--     panel.grid = element_blank(), -->
<!--     panel.border = element_rect(color = "black", fill = NA, size = 0.8), -->
<!--     axis.line = element_blank(), -->
<!--     plot.background = element_rect(fill = "white", color = NA) -->
<!--   ) -->
<!-- plot_diff -->
<!-- lim <- 0.5 -->
<!-- plot_diff_zoom <- ggplot(df_diff, aes(x = x, y = y, fill = diff_conc)) + -->
<!--   geom_raster() + -->
<!--   scale_fill_gradient2( -->
<!--     name = "Hg(0) difference (ng/m³)\n2020 - 2003", -->
<!--     low = colors_div[1], -->
<!--     mid = colors_div[2], -->
<!--     high = colors_div[3], -->
<!--     midpoint = 0, -->
<!--     limits = c(-lim, lim), -->
<!--     oob = scales::squish -->
<!--   ) + -->
<!--   coord_equal() + -->
<!--   theme_minimal() + -->
<!--   labs( -->
<!--     title = paste0("Average annual difference of Hg(0) concentration between 2020 and 2003"), -->
<!--         subtitle = "Dataset : EDGAR", -->
<!--     x = "Longitude", y = "Latitude" -->
<!--   ) + -->
<!--   borders("world", colour = "#403d39", size = 0.2) + -->
<!--    theme( -->
<!--     plot.title = element_text(hjust = 0.5, face = "bold"), -->
<!--     plot.subtitle = element_text(hjust = 0.5), -->
<!--     axis.text = element_text(size = 10), -->
<!--     axis.ticks = element_line(color = "black", size = 0.3), -->
<!--     axis.ticks.length = unit(3, "pt"), -->
<!--     panel.grid = element_blank(), -->
<!--     panel.border = element_rect(color = "black", fill = NA, size = 0.8), -->
<!--     axis.line = element_blank(), -->
<!--     plot.background = element_rect(fill = "white", color = NA) -->
<!--   ) -->
<!-- plot_diff_zoom -->
<!-- # Chemin de sauvegarde (modifie si nécessaire) -->
<!-- chemin_pdf <- "C:/Users/colom/Desktop/STAGE/data/clean_mod_data/plot/EDGAR_diff_20_03_Conc_WORLD.pdf" -->
<!-- # Enregistrement du plot -->
<!-- ggsave( -->
<!--   filename = chemin_pdf, -->
<!--   plot = plot_diff_zoom, -->
<!--   device = "pdf", -->
<!--   width = 10, height = 6, units = "in" -->
<!-- ) -->
<!-- ``` -->

``` r
nc_mask<-nc_open(file.path("C:/Users/colom/Desktop/STAGE/data/clean_mod_data", "Amazon_basin_mask_2x25.nc"))
Am_mask<-ncvar_get(nc_mask, "MASK")

dim(Am_mask)
```

    ## [1] 144  91

``` r
range(Am_mask)
```

    ## [1] 0 1

``` r
df_amazon <- expand.grid(
  x = lon,  # tes longitudes
  y = lat   # tes latitudes
) %>%
  mutate(mask = as.vector(Am_mask))

# Affichage de la carte
ggplot(df_amazon, aes(x = x, y = y, fill = factor(mask))) +
  geom_raster() +
  scale_fill_manual(
    values = c("0" = "transparent", "1" = "#0077BB"),
    labels = c("0" = "Non-Amazonie", "1" = "Amazonie"),
    name = "Masque Amazonie"
  ) +
  borders("world", colour = "black", size = 0.3) +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Masque Amazonien",
    x = "Longitude", y = "Latitude"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "right",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.6),
    plot.background = element_rect(fill = "white", color = NA)
  )
```

![](EDGAR_Map_Diff_Conc_01_08_2025_files/figure-gfm/Amazon%20mask-1.png)<!-- -->

``` r
plot_mask<-ggplot(df_amazon, aes(x = x, y = y, fill = factor(mask))) +
  geom_raster(na.rm = TRUE) +
  scale_fill_manual(
    values = c("0" = "transparent", "1" = "#0077BB"),
    labels = c("1" = "Amazon"),
    name = "Mask",
    na.value = "transparent"
  ) +
  borders("world", colour = "black", size = 0.3) +
  coord_quickmap(xlim = c(-85, -30), ylim = c(-60, 15)) +
  theme_minimal() +
  labs(
    title = "Amazon region mask",
    x = "Longitude", y = "Latitude"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "right",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.6),
    plot.background = element_rect(fill = "white", color = NA)
  )
print(plot_mask)
```

![](EDGAR_Map_Diff_Conc_01_08_2025_files/figure-gfm/Amazon%20mask-2.png)<!-- -->

``` r
# Chemin de sauvegarde (modifie si nécessaire)
chemin_pdf <- "C:/Users/colom/Desktop/STAGE/data/clean_mod_data/plot/Amazon_mask.pdf"

# Enregistrement du plot
ggsave(
  filename = chemin_pdf,
  plot = plot_mask,
  device = "pdf",
  width = 10, height = 6, units = "in"
)
```

``` r
lim <- 0.2
plot_diff_zoom <- ggplot(df_diff, aes(x = x, y = y, fill = diff_conc)) +
  geom_raster() +
  scale_fill_gradient2(
    name = "Difference Hg(0) (ng/m³)\n2020 - 2003",
    low = colors_div[1],
    mid = colors_div[2],
    high = colors_div[3],
    midpoint = 0,
    limits = c(-lim, lim),
    oob = scales::squish
  ) +
  coord_equal() +
  coord_quickmap(xlim = c(-85, -30), ylim = c(-60, 15)) +

  theme_minimal() +
  labs(
    title = paste0("Average annual difference of Hg(0) concentration between 2020 and 2003"),
        subtitle = "Dataset : EDGAR - Latin America focus",
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

![](EDGAR_Map_Diff_Conc_01_08_2025_files/figure-gfm/Focus%20amazonie-1.png)<!-- -->

``` r
# Chemin de sauvegarde (modifie si nécessaire)
chemin_pdf <- "C:/Users/colom/Desktop/STAGE/data/clean_mod_data/plot/EDGAR_Mean_annual_difference_amazonia.pdf"

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
  path_sim <- file.path(chem_EDGAR, sim, "GEOSChem.SpeciesConc.2015_m.nc4")
  nc_sim <- nc_open(path_sim)
  
  # 2. Lire la variable Hg0 (niveau 1)
  hg0 <- ncvar_get(nc_sim, "SpeciesConcVV_Hg0")[, , 1, ]
  
  # 3. Gérer le temps et pondérer
  time_raw <- ncvar_get(nc_sim, "time")
  origin_time <- as.POSIXct("2015-01-01 00:00:00", tz = "UTC")
  time <- origin_time + time_raw * 60
  weights <- days_in_month(time) / sum(days_in_month(time))
  
  # 4. Moyenne pondérée temporelle
  weighted_avg <- function(x) sum(x * weights)
  hg0_mean <- apply(hg0, MARGIN = c(1, 2), FUN = weighted_avg)
  
  # 5. Conversion en ng/m³
  hg0_ng_m3 <- hg0_mean * (M_Hg / M_air) * mean_Met_AIRDEN * ng_per_kg
  
  # 6. Appliquer le masque Amazonie
  hg0_amz <- hg0_ng_m3 * Am_mask
  
  # 7. Stocker le résultat
  results_amz[[sim]] <- hg0_amz
  
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

    ## [1] 0.000000 1.001531

``` r
range(V03A20_am)
```

    ## [1] 0.000000 1.637516

``` r
range(V20A03_am)
```

    ## [1] 0.000000 1.053479

``` r
range(V20A20_am)
```

    ## [1] 0.000000 1.652427

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
  Vegetation  = 100 * total_V / total_sum
)



df <- data.frame(
  source = names(contrib),
  pct = contrib
)

df$label <- paste0(df$source, "\n", sprintf("%.1f", df$pct), "%")

contrib_conc_ama<-ggplot(df, aes(x = "", y = pct, fill = source)) +
  geom_col(width = 1, color = "white") +
  geom_text(aes(label = label), 
            position = position_stack(vjust = 0.5), 
            color = "white", size = 5, fontface = "bold") +
  coord_polar("y") +
  scale_fill_manual(values = c("#D95F02", "#1B9E77")) +
  labs(
    title = "Relative contribution to the change in Hg(0) concentration",
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
print(contrib_conc_ama)
```

![](EDGAR_Map_Diff_Conc_01_08_2025_files/figure-gfm/calcul%20contrib%20vegetation%20amazo-1.png)<!-- -->

``` r
# Chemin de sauvegarde (modifie si nécessaire)
chemin_pdf <- "C:/Users/colom/Desktop/STAGE/data/clean_mod_data/plot/EDGAR_Contrib_Conc_Amazonia.pdf"

# Enregistrement du plot
ggsave(
  filename = chemin_pdf,
  plot = contrib_conc_ama,
  device = "pdf",
  width = 10, height = 6, units = "in"
)
```

``` r
simulations <- c("V03A03", "V03A20", "V20A03", "V20A20")

results <- list()

for (sim in simulations) {
  
  # 1. Ouvrir le fichier NetCDF
  path_sim <- file.path(chem_EDGAR, sim, "GEOSChem.SpeciesConc.2015_m.nc4")
  nc_sim <- nc_open(path_sim)
  
  # 2. Lire la variable Hg0 (niveau 1)
  hg0 <- ncvar_get(nc_sim, "SpeciesConcVV_Hg0")[, , 1, ]
  
  # 3. Gérer le temps et pondérer
  time_raw <- ncvar_get(nc_sim, "time")
  origin_time <- as.POSIXct("2015-01-01 00:00:00", tz = "UTC")
  time <- origin_time + time_raw * 60
  weights <- days_in_month(time) / sum(days_in_month(time))
  
  # 4. Moyenne pondérée temporelle
  weighted_avg <- function(x) sum(x * weights)
  hg0_mean <- apply(hg0, MARGIN = c(1, 2), FUN = weighted_avg)
  
  # 5. Conversion en ng/m³
  hg0_ng_m3 <- hg0_mean * (M_Hg / M_air) * mean_Met_AIRDEN * ng_per_kg
  
 
  hg0 <- hg0_ng_m3
  
  # 7. Stocker le résultat
  results[[sim]] <- hg0
  
  # Fermer le fichier
  nc_close(nc_sim)
}

# Résultat : une liste results_amz avec 4 objets nommés :
# results_amz$V03A03, results_amz$V03A20, etc.
V03A03 <- results[["V03A03"]]
V03A20 <- results[["V03A20"]]
V20A03 <- results[["V20A03"]]
V20A20 <- results[["V20A20"]]

range(V03A03)
```

    ## [1] 0.4419279 2.2120839

``` r
mean(V03A03)
```

    ## [1] 0.8406278

``` r
sd(V03A03)
```

    ## [1] 0.09633196

``` r
range(V03A20)
```

    ## [1] 0.4952566 2.5424424

``` r
range(V20A03)
```

    ## [1] 0.4679975 2.2168357

``` r
range(V20A20)
```

    ## [1] 0.5242534 2.5489646

``` r
mean(V20A20)
```

    ## [1] 0.8916833

``` r
sd(V20A20)
```

    ## [1] 0.1082034

``` r
#1. Décomposition spatiale : SCREEN DANS DOCS GEOSCHEM
Eff_A <- 0.5 * ((V03A20 - V03A03) + (V20A20 - V20A03))
Eff_V <- 0.5 * ((V20A03 - V03A03) + (V20A20 - V03A20))
Interaction <- V20A20 - V03A03 - Eff_A - Eff_V



#3. Moyennes régionales :
mean_A <- mean(Eff_A, na.rm = TRUE)
mean_V <- mean(Eff_V, na.rm = TRUE)
mean_I <- mean(Interaction, na.rm = TRUE)

#4. Camembert avec interaction répartie :
total_A <- mean_A + 0.5 * mean_I
total_V <- mean_V + 0.5 * mean_I
total_sum <- total_A + total_V

contrib <- c(
  Emissions = 100 * total_A / total_sum,
  Vegetation  = 100 * total_V / total_sum
)



df <- data.frame(
  source = names(contrib),
  pct = contrib
)

df$label <- paste0(df$source, "\n", sprintf("%.1f", df$pct), "%")

contrib_conc_tot<-ggplot(df, aes(x = "", y = pct, fill = source)) +
  geom_col(width = 1, color = "white") +
  geom_text(aes(label = label), 
            position = position_stack(vjust = 0.5), 
            color = "white", size = 5, fontface = "bold") +
  coord_polar("y") +
  scale_fill_manual(values = c("#D95F02", "#1B9E77")) +
  labs(
    title = "Relative contribution to the change in Hg(0) concentration",
    subtitle = "Dataset : EDGAR / Global",
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
print(contrib_conc_tot)
```

![](EDGAR_Map_Diff_Conc_01_08_2025_files/figure-gfm/calcul%20contrib%20vegetation%20total-1.png)<!-- -->

``` r
# Chemin de sauvegarde (modifie si nécessaire)
chemin_pdf <- "C:/Users/colom/Desktop/STAGE/data/clean_mod_data/plot/EDGAR_Contrib_Conc_GLOBAL.pdf"

# Enregistrement du plot
ggsave(
  filename = chemin_pdf,
  plot = contrib_conc_tot,
  device = "pdf",
  width = 10, height = 6, units = "in"
)
```

``` r
plots <- contrib_conc_tot + contrib_conc_ama +   plot_layout(ncol = 2)
plots  
```

![](EDGAR_Map_Diff_Conc_01_08_2025_files/figure-gfm/contrib%20vegetation%20total%202%20plots-1.png)<!-- -->
