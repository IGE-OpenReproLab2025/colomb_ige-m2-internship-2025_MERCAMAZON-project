Mercury mass in atmosphere
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

# Lecture des données

``` r
chem_hist_B100<-"C:/Users/colom/Desktop/STAGE/data/clean_mod_data/14_3_1/HIST_B100"
chem_hist<-"C:/Users/colom/Desktop/STAGE/data/clean_mod_data/14_3_1"
con_hg_B100<-nc_open(file.path(chem_hist_B100, "GEOSChem.SpeciesConc.2015_m.nc4"))
statemet_B100<-nc_open(file.path(chem_hist, "GEOSChem.StateMet.2014_2015_m.nc4"))

conc_hg0<-ncvar_get(con_hg_B100, "SpeciesConcVV_Hg0") #Concentration of species Hg0 (mol mol-1 dry)

Met_AIRVOL<-ncvar_get(statemet_B100, "Met_AIRVOL")  #Met_AIRVOL (Volume of grid box m3)
Met_AIRDEN<-ncvar_get(statemet_B100, "Met_AIRDEN")  #Met_AIRDEN (Dry air density kg m-3)
Met_AD<-ncvar_get(statemet_B100, "Met_AD")  #Met_AD (Dry air mass kg)
mean_Met_AD <- apply(Met_AD, c(1,2,3), mean, na.rm = TRUE)


#Amazon basin mmask
nc_mask<-nc_open(file.path("C:/Users/colom/Desktop/STAGE/data/clean_mod_data", "Amazon_basin_mask_2x25.nc"))
Am_mask<-ncvar_get(nc_mask, "MASK")

area<-ncvar_get(con_hg_B100, "AREA")
s_in_yr = 365.2425 * 24 * 3600
unit_conv = s_in_yr
kg_Mg = 1e-3
MW_Hg = 200.59
avo = 6.02e23
g_kg = 1e-3
cm2_m2 = 1e4
unit_conv_dd = MW_Hg / avo * g_kg * cm2_m2 * s_in_yr * area
M_Hg <- 0.20059
M_air <- 0.028964
```

# The conversion formula is:

        mol species                g/mol species
        ----------- * kg dry air * -------------
        mol dry air                g/mol dry air

        = volume mixing ratio * dry air mass * MW species / MW dry air

        = kg species
        
            
        MW dry air is defined as 28.9644 g/mol
        https://github.com/geoschem/geos-chem/blob/4722f288e90291ba904222f4bbe4fc216d17c34a/Headers/physconstants.F90#L26

        
        

# Calcul de la masse de mercure dans l’atmosphère sur l’année

``` r
subdirs <- list.dirs(chem_hist, full.names = TRUE, recursive = FALSE)
subdirs <- subdirs[grepl("HIST_", basename(subdirs))]

results_mass_Hg_total <- tibble()

for (subdir in subdirs) {
  file_path <- file.path(subdir, "GEOSChem.SpeciesConc.2015_m.nc4")
  if (!file.exists(file_path)) next
  
  message("Traitement de : ", basename(subdir))
  
  nc <- nc_open(file_path)
  
  # Récupérer les noms des variables contenant "Hg"
  hg_vars <- names(nc$var)
  hg_vars <- hg_vars[grepl("Hg", hg_vars)]
  
  if (length(hg_vars) == 0) {
    warning(paste("Aucune variable Hg trouvée dans", file_path))
    nc_close(nc)
    next
  }
  
  # Initialiser Hg_total avec la première variable Hg
  Hg_total <- ncvar_get(nc, hg_vars[1])
  
  # Ajouter les autres variables Hg
  if (length(hg_vars) > 1) {
    for (varname in hg_vars[-1]) {
      Hg_total <- Hg_total + ncvar_get(nc, varname)
    }
  }
  
  nc_close(nc)
  
  # Moyenne temporelle sur les 12 mois
  mean_Hg_total <- apply(Hg_total, c(1, 2, 3), mean, na.rm = TRUE)
  mean_Hg_total <- mean_Hg_total[, , 1:dim(mean_Hg_total)[3]]  # si besoin

  # Calcul de la masse Hg (kg) pour chaque niveau vertical
  mass_Hg_total <- mean_Hg_total * mean_Met_AD * (M_Hg / M_air)

  # Somme verticale (colonne atmosphérique)
  mass_Hg_total <- apply(mass_Hg_total, c(1, 2), sum, na.rm = TRUE)

  # Appliquer le masque Amazonie
  mass_Hg_total <- mass_Hg_total * Am_mask

  # Masse totale Hg en Amazonie (tonnes)
  total_mass_Hg_tonnes <- sum(mass_Hg_total, na.rm = TRUE) / 1000

  # Stocker le résultat
  results_mass_Hg_total <- bind_rows(results_mass_Hg_total, tibble(
    case = basename(subdir),
    total_Hg_mass_tonnes = total_mass_Hg_tonnes
  ))
}

print(results_mass_Hg_total)
```

    ## # A tibble: 10 × 2
    ##    case          total_Hg_mass_tonnes
    ##    <chr>                        <dbl>
    ##  1 HIST_B10                      55.9
    ##  2 HIST_B100                     57.5
    ##  3 HIST_B150                     58.0
    ##  4 HIST_B20                      56.1
    ##  5 HIST_B200                     58.4
    ##  6 HIST_B30                      56.4
    ##  7 HIST_B40                      56.6
    ##  8 HIST_Feinberg                 77.5
    ##  9 HIST_V3                       55.6
    ## 10 HIST_V4                       55.6

``` r
subdirs <- list.dirs(chem_hist, full.names = TRUE, recursive = FALSE)
subdirs <- subdirs[grepl("HIST_", basename(subdirs))]

results_mass_Hg0 <- tibble()

for (subdir in subdirs) {
  file_path <- file.path(subdir, "GEOSChem.SpeciesConc.2015_m.nc4")
  if (!file.exists(file_path)) next
  
  message("Traitement de : ", basename(subdir))
  
  nc <- nc_open(file_path)
  
  # Lire la variable Hg0
  conc_hg0 <- ncvar_get(nc, "SpeciesConcVV_Hg0")
  
  nc_close(nc)
  
  # Moyenne temporelle sur les 12 mois
  mean_conc_hg0 <- apply(conc_hg0, c(1, 2, 3), mean, na.rm = TRUE)
  
  # Calcul de la masse Hg0 (kg) par niveau vertical
  mass_hg0 <- mean_conc_hg0 * mean_Met_AD * (M_Hg / M_air)
  
  # Somme verticale (colonne atmosphérique)
  mass_hg0 <- apply(mass_hg0, c(1, 2), sum, na.rm = TRUE)
  
  # Application du masque Amazonie
  mass_hg0 <- mass_hg0 * Am_mask
  
  # Masse totale Hg0 en Amazonie (tonnes)
  total_mass_Hg0_tonnes <- sum(mass_hg0, na.rm = TRUE) / 1000
  
  # Stockage du résultat
  results_mass_Hg0 <- bind_rows(results_mass_Hg0, tibble(
    case = basename(subdir),
    total_Hg0_mass_tonnes = total_mass_Hg0_tonnes
  ))
}

print(results_mass_Hg0)
```

    ## # A tibble: 10 × 2
    ##    case          total_Hg0_mass_tonnes
    ##    <chr>                         <dbl>
    ##  1 HIST_B10                       54.2
    ##  2 HIST_B100                      55.8
    ##  3 HIST_B150                      56.3
    ##  4 HIST_B20                       54.4
    ##  5 HIST_B200                      56.8
    ##  6 HIST_B30                       54.6
    ##  7 HIST_B40                       54.9
    ##  8 HIST_Feinberg                  75.3
    ##  9 HIST_V3                        53.9
    ## 10 HIST_V4                        53.9

``` r
# Dossier principal
chem_hist <- "C:/Users/colom/Desktop/STAGE/data/clean_mod_data/14_3_1"

# Lister les sous-dossiers HIST_*
subdirs <- list.dirs(chem_hist, full.names = TRUE, recursive = FALSE)
subdirs <- subdirs[grepl("HIST_", basename(subdirs))]

# Résultat final
results <- tibble()

for (subdir in subdirs) {
  file_path <- file.path(subdir, "HEMCO_diagnostic_2015_m.nc4")
  if (!file.exists(file_path)) next
  
  message("Traitement de : ", basename(subdir))
  
  hemco_diagn <- nc_open(file_path)
  
  # Temps brut et pondération mensuelle
  time_raw <- ncvar_get(hemco_diagn, "time")
  origin_time <- as.POSIXct("2015-01-01 00:00:00", tz = "UTC")
  time <- origin_time + time_raw * 60
  weights <- days_in_month(time) / sum(days_in_month(time))
  
  # Aire des cellules
  area <- ncvar_get(hemco_diagn, "AREA")  # [lon, lat]
  
  # Variables à traiter
  vars <- list(
    biomass      = "EmisHg0_BioBurn",
    hg0_anthro   = "EmisHg0_Anthro",
    hg2_anthro   = "EmisHg2ClP_Anthro",
    hgp_anthro   = "EmisHgCl2_Anthro",
    geo          = "EmisHg0_Natural",
    asgm         = "EmisHg0_ASGM"
  )
  
  row <- list(case = basename(subdir))
  
  for (name in names(vars)) {
    varname <- vars[[name]]
    emis <- ncvar_get(hemco_diagn, varname)  # [lon, lat, time]
    
    # Moyenne pondérée temporelle
    emis_avg <- apply(emis, c(1, 2), function(x) sum(x * weights))
    
    # Appliquer le masque Amazonie
    emis_avg <- emis_avg * Am_mask  # [lon, lat]
    
    # Total des émissions (tonnes) en Amazonie
    total <- sum(emis_avg * area * s_in_yr * kg_Mg, na.rm = TRUE)
    row[[name]] <- total
  }
  
  # Combinaisons supplémentaires
  row$total_hg0 <- row$hg0_anthro + row$asgm
  row$total_hg2_hgp <- row$hg2_anthro + row$hgp_anthro
  
  # Ajout au tableau final
  results <- bind_rows(results, as_tibble(row))
  
  nc_close(hemco_diagn)
}

print(results)
```

    ## # A tibble: 8 × 9
    ##   case          biomass hg0_anthro hg2_anthro hgp_anthro   geo  asgm total_hg0 total_hg2_hgp
    ##   <chr>           <dbl>      <dbl>      <dbl>      <dbl> <dbl> <dbl>     <dbl>         <dbl>
    ## 1 HIST_B10         25.8       8.31      0.549       1.75  11.5  195.      203.          2.30
    ## 2 HIST_B100        25.8       8.31      0.549       1.75  11.5  195.      203.          2.30
    ## 3 HIST_B150        25.8       8.31      0.549       1.75  11.5  195.      203.          2.30
    ## 4 HIST_B20         25.8       8.31      0.549       1.75  11.5  195.      203.          2.30
    ## 5 HIST_B200        25.8       8.31      0.549       1.75  11.5  195.      203.          2.30
    ## 6 HIST_B30         25.8       8.31      0.549       1.75  11.5  195.      203.          2.30
    ## 7 HIST_B40         25.8       8.31      0.549       1.75  11.5  195.      203.          2.30
    ## 8 HIST_Feinberg    25.8       8.31      0.549       1.75  11.5  195.      203.          2.30

``` r
results_mercuryemis <- tibble()

for (subdir in subdirs) {
  file_path <- file.path(subdir, "GEOSChem.MercuryEmis.alltime_m.nc4")
  
  if (!file.exists(file_path)) next
  
  message("Traitement de : ", basename(subdir))
  
  merc_emiss <- nc_open(file_path)
  
  # Temps et pondération mensuelle
  time_raw <- ncvar_get(merc_emiss, "time")
  origin_time <- as.POSIXct("2015-01-01 00:00:00", tz = "UTC")
  time <- origin_time + time_raw * 60
  weights <- days_in_month(time) / sum(days_in_month(time))
  weighted_avg <- function(x) sum(x * weights)
  
  # Emissions Hg0
  hg0land <- ncvar_get(merc_emiss, "EmisHg0land")  # [lon, lat, time]
  hg0soil <- ncvar_get(merc_emiss, "EmisHg0soil")
  
  # Moyenne temporelle pondérée
  hg0land_mean <- apply(hg0land, c(1, 2), weighted_avg)
  hg0soil_mean <- apply(hg0soil, c(1, 2), weighted_avg)
  
  # Application du masque Amazonie
  hg0land_mean <- hg0land_mean * Am_mask
  hg0soil_mean <- hg0soil_mean * Am_mask
  
  # Totaux (en tonnes)
  total_hg0land <- sum(hg0land_mean, na.rm = TRUE) * s_in_yr * kg_Mg
  total_hg0soil <- sum(hg0soil_mean, na.rm = TRUE) * s_in_yr * kg_Mg
  
  results_mercuryemis <- bind_rows(results_mercuryemis, tibble(
    case = basename(subdir),
    hg0land = total_hg0land,
    hg0soil = total_hg0soil
  ))
  
  nc_close(merc_emiss)
}

print(results_mercuryemis)
```

    ## # A tibble: 7 × 3
    ##   case      hg0land hg0soil
    ##   <chr>       <dbl>   <dbl>
    ## 1 HIST_B10     3.58    84.5
    ## 2 HIST_B100    3.16    84.5
    ## 3 HIST_B150    3.02    84.5
    ## 4 HIST_B20     3.52    84.5
    ## 5 HIST_B200    2.90    84.5
    ## 6 HIST_B30     3.46    84.5
    ## 7 HIST_B40     3.41    84.5

``` r
subdirs <- list.dirs(chem_hist, full.names = TRUE, recursive = FALSE)
subdirs <- subdirs[grepl("HIST_", basename(subdirs))]

# Initialiser les résultats
results_drydep <- tibble()

# Fonction de moyenne pondérée
weighted_avg <- function(x, weights) sum(x * weights)

# Boucle sur les cas
for (subdir in subdirs) {
  file_path <- file.path(subdir, "GEOSChem.DryDep.2015_m.nc4")
  if (!file.exists(file_path)) next
  
  message("Traitement de : ", basename(subdir))
  
  drydep_nc <- nc_open(file_path)
  
  # Temps
  time_raw <- ncvar_get(drydep_nc, "time")
  origin_time <- as.POSIXct("2015-01-01 00:00:00", tz = "UTC")
  time <- origin_time + time_raw * 60
  days_in_months <- days_in_month(time)
  weights <- days_in_months / sum(days_in_months)
  
  # Extraction et traitement des 3 espèces avec masque Amazonie
  get_drydep_sum <- function(varname) {
    var <- ncvar_get(drydep_nc, varname)                     # [lon, lat, time]
    var <- apply(var, c(1, 2), weighted_avg, weights = weights)  # [lon, lat]
    var <- var * Am_mask                                     # masque Amazonie
    var <- var * unit_conv_dd * kg_Mg                        # conversion en tonnes
    sum(var, na.rm = TRUE)
  }
  
  DryDep_Hg0_hist <- get_drydep_sum("DryDep_Hg0")
  DryDep_Hg2_hist <- get_drydep_sum("DryDep_Hg2")
  DryDep_HgP_hist <- get_drydep_sum("DryDep_HgP")

  # Stocker les résultats
  results_drydep <- bind_rows(results_drydep, tibble(
    case = basename(subdir),
    drydep_Hg0_tonnes = DryDep_Hg0_hist,
    drydep_Hg2_tonnes = DryDep_Hg2_hist,
    drydep_HgP_tonnes = DryDep_HgP_hist
  ))
  
  nc_close(drydep_nc)
}

# Afficher les résultats finaux
print(results_drydep)
```

    ## # A tibble: 9 × 4
    ##   case          drydep_Hg0_tonnes drydep_Hg2_tonnes drydep_HgP_tonnes
    ##   <chr>                     <dbl>             <dbl>             <dbl>
    ## 1 HIST_B10                   549.             11.0             0.0669
    ## 2 HIST_B100                  563.              9.30            0.0690
    ## 3 HIST_B150                  568.              8.73            0.0698
    ## 4 HIST_B20                   551.             10.7             0.0672
    ## 5 HIST_B200                  572.              8.28            0.0703
    ## 6 HIST_B30                   553.             10.5             0.0675
    ## 7 HIST_B40                   555.             10.3             0.0678
    ## 8 HIST_Feinberg              570.             13.5             0.0116
    ## 9 HIST_V3                    547.             11.2             0.0665

``` r
# Sous-dossiers HIST_*
subdirs <- list.dirs(chem_hist, full.names = TRUE, recursive = FALSE)
subdirs <- subdirs[grepl("HIST_", basename(subdirs))]

# Initialiser le tableau des résultats
results_wetloss <- tibble()

# Fonction de moyenne pondérée
weighted_avg <- function(x, weights) sum(x * weights)

# Boucle sur les cas
for (subdir in subdirs) {
  file_path <- file.path(subdir, "GEOSChem.WetLossTot_HgSum.2015_m.nc4")
  if (!file.exists(file_path)) next

  message("Traitement de : ", basename(subdir))
  
  wetloss_nc <- nc_open(file_path)
  
  # Temps
  time_raw <- ncvar_get(wetloss_nc, "time")
  origin_time <- as.POSIXct("2015-01-01 00:00:00", tz = "UTC")
  time <- origin_time + time_raw * 60
  days_in_months <- days_in_month(time)
  weights <- days_in_months / sum(days_in_months)
  
  # Extraction et traitement
  wetloss_array <- ncvar_get(wetloss_nc, "WetLossTot_HgSum")  # [lon, lat, time]
  wetloss_array <- apply(wetloss_array, c(1, 2), weighted_avg, weights = weights)  # [lon, lat]
  wetloss_array <- wetloss_array * Am_mask   # Appliquer le masque Amazonie
  wetloss_total_tonnes <- sum(wetloss_array * s_in_yr * kg_Mg, na.rm = TRUE)
  
  # Stocker
  results_wetloss <- bind_rows(results_wetloss, tibble(
    case = basename(subdir),
    wetloss_Hg_total_tonnes = wetloss_total_tonnes
  ))
  
  nc_close(wetloss_nc)
}

# Afficher les résultats
print(results_wetloss)
```

    ## # A tibble: 9 × 2
    ##   case          wetloss_Hg_total_tonnes
    ##   <chr>                           <dbl>
    ## 1 HIST_B10                         34.4
    ## 2 HIST_B100                        33.0
    ## 3 HIST_B150                        32.6
    ## 4 HIST_B20                         34.2
    ## 5 HIST_B200                        32.3
    ## 6 HIST_B30                         34.0
    ## 7 HIST_B40                         33.8
    ## 8 HIST_Feinberg                    45.3
    ## 9 HIST_V3                          34.7

``` r
# Fusionner les tables selon la colonne commune "case"
results_all <- results_wetloss %>%
  full_join(results_drydep, by = "case") %>%
  full_join(results_mercuryemis, by = "case") %>%
  full_join(results, by = "case")

# Écrire dans un fichier Excel
write_xlsx(results_all, "C:/Users/colom/Desktop/STAGE/data/clean_mod_data/14_3_1/resultats_mercury_amazon.xlsx")


# Fusionner les tables selon la colonne commune "case"
results_mass_all <- results_mass_Hg_total  %>%
  full_join(results_mass_Hg0, by = "case")
# Écrire dans un fichier Excel
write_xlsx(results_mass_all, "C:/Users/colom/Desktop/STAGE/data/clean_mod_data/14_3_1/resultats_mass_mercury_amazon.xlsx")
```
