Mercury mass in atmosphere - QIU INVENTORY
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
chem_QIU <- "C:/Users/colom/Desktop/STAGE/data/clean_mod_data/14_3_1/QIU"
con_hg_B100<-nc_open(file.path(chem_hist_B100, "GEOSChem.SpeciesConc.2015_m.nc4"))
statemet_B100<-nc_open(file.path(chem_hist_B100, "GEOSChem.StateMet.2015_m.nc4"))
area<-ncvar_get(con_hg_B100, "AREA")

conc_hg0<-ncvar_get(con_hg_B100, "SpeciesConcVV_Hg0") #Concentration of species Hg0 (mol mol-1 dry)

Met_AD<-ncvar_get(statemet_B100, "Met_AD")  #Met_AD (Dry air mass kg)

area<-ncvar_get(con_hg_B100, "AREA")

s_in_yr = 365.2425 * 24 * 3600
unit_conv = s_in_yr
kg_Mg = 1e-3
MW_Hg = 200.59
avo = 6.02e23
g_kg = 1e-3
cm2_m2 = 1e4
unit_conv_dd = MW_Hg / avo * g_kg * cm2_m2 * s_in_yr * area
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
# Moyenne temporelle sur les 12 mois
mean_Met_AD <- apply(Met_AD, c(1,2,3), mean, na.rm = TRUE)



subdirs <- list.dirs(chem_QIU, full.names = TRUE, recursive = FALSE)
subdirs <- subdirs[grepl("^V", basename(subdirs))]

# Résultat
results_mass_Hg0 <- tibble()

for (subdir in subdirs) {
  file_path <- file.path(subdir, "GEOSChem.SpeciesConc.2015_m.nc4")
  if (!file.exists(file_path)) next
  
  nc <- nc_open(file_path)
  
  # Concentration de Hg0
  conc_hg0 <- ncvar_get(nc, "SpeciesConcVV_Hg0")  # [mol/mol]
  nc_close(nc)
  
  # Moyenne temporelle sur les 12 mois
  mean_conc_hg0 <- apply(conc_hg0, c(1,2,3), mean, na.rm = TRUE)
  
  # Masse (kg) avec mean_Met_AD déjà disponible
  mass_hg0 <- mean_conc_hg0 * mean_Met_AD * (M_Hg / M_air)
  total_mass_hg0_tonnes <- sum(mass_hg0, na.rm = TRUE) / 1000
  
  results_mass_Hg0 <- bind_rows(results_mass_Hg0, tibble(
    case = basename(subdir),
    total_Hg0_mass_tonnes = total_mass_hg0_tonnes
  ))
}

# Résultat final
print(results_mass_Hg0)
```

    ## # A tibble: 2 × 2
    ##   case   total_Hg0_mass_tonnes
    ##   <chr>                  <dbl>
    ## 1 V03A03                 2642.
    ## 2 V20A03                 2660.

``` r
subdirs <- list.dirs(chem_QIU, full.names = TRUE, recursive = FALSE)
subdirs <- subdirs[grepl("^V", basename(subdirs))]

results_mass_Hg_total <- tibble()

for (subdir in subdirs) {
  file_path <- file.path(subdir, "GEOSChem.SpeciesConc.2015_m.nc4")
  if (!file.exists(file_path)) next
  
  nc <- nc_open(file_path)
  
  # Récupérer tous les noms de variables contenant "Hg"
  hg_vars <- names(nc$var)
  hg_vars <- hg_vars[grepl("Hg", hg_vars)]
  
  # Initialiser la somme avec la première variable Hg
  Hg_total <- ncvar_get(nc, hg_vars[1])
  
  # Ajouter les autres variables Hg
  if (length(hg_vars) > 1) {
    for (varname in hg_vars[-1]) {
      Hg_total <- Hg_total + ncvar_get(nc, varname)
    }
  }
  
  nc_close(nc)
  
  # Moyenne temporelle sur les 12 mois
  mean_Hg_total <- apply(Hg_total, c(1,2,3), mean, na.rm = TRUE)
  
  # Calcul de la masse (kg) avec mean_Met_AD déjà disponible
  mass_Hg_total <- mean_Hg_total * mean_Met_AD * (M_Hg / M_air)
  total_mass_Hg_tonnes <- sum(mass_Hg_total, na.rm = TRUE) / 1000
  
  results_mass_Hg_total <- bind_rows(results_mass_Hg_total, tibble(
    case = basename(subdir),
    total_Hg_mass_tonnes = total_mass_Hg_tonnes
  ))
}

print(results_mass_Hg_total)
```

    ## # A tibble: 2 × 2
    ##   case   total_Hg_mass_tonnes
    ##   <chr>                 <dbl>
    ## 1 V03A03                2742.
    ## 2 V20A03                2760.

``` r
# Lister les sous-dossiers HIST_*
subdirs <- list.dirs(chem_QIU, full.names = TRUE, recursive = FALSE)
subdirs <- subdirs[grepl("^V", basename(subdirs))]

# Initialiser un tableau pour stocker les résultats
results <- tibble()

for (subdir in subdirs) {
  file_path <- file.path(subdir, "HEMCO_diagnostic_2015_m.nc4")
  
  if (!file.exists(file_path)) next  # skip si fichier manquant
  
  hemco_diagn <- nc_open(file_path)
  
  # Extraire temps et poids mensuels
  time_raw <- ncvar_get(hemco_diagn, "time")
  origin_time <- as.POSIXct("2015-01-01 00:00:00", tz = "UTC")
  time <- origin_time + time_raw * 60
  weights <- days_in_month(time) / sum(days_in_month(time))
  
  # Fonction de moyenne pondérée
  weighted_avg <- function(x) sum(x * weights)
  
  # Aire des cellules
  area <- ncvar_get(hemco_diagn, "AREA")  # en m²

  # Liste des variables d'intérêt
  vars <- list(
    biomass      = "EmisHg0_BioBurn",
    hg0_anthro   = "EmisHg0_Anthro",
    hg2_anthro   = "EmisHg2ClP_Anthro",
    hgp_anthro   = "EmisHgCl2_Anthro",
    geo          = "EmisHg0_Natural",
    asgm         = "EmisHg0_ASGM"
  )
  
  # Stockage temporaire
  row <- list(case = basename(subdir))
  
  for (name in names(vars)) {
    varname <- vars[[name]]
    emis <- ncvar_get(hemco_diagn, varname)  # [lon, lat, time]
    emis_avg <- apply(emis, c(1, 2), weighted_avg)
    total <- sum(emis_avg * area * s_in_yr * kg_Mg, na.rm = TRUE)
    row[[name]] <- total
  }
  
  # Ajouter les calculs additionnels
  row$total_hg0      <- row$hg0_anthro + row$asgm
  row$total_hg2_hgp  <- row$hg2_anthro + row$hgp_anthro
  
  # Ajouter au tableau final
  results <- bind_rows(results, as_tibble(row))
  
  nc_close(hemco_diagn)
}

print(results)
```

    ## # A tibble: 2 × 9
    ##   case   biomass hg0_anthro hg2_anthro hgp_anthro   geo  asgm total_hg0 total_hg2_hgp
    ##   <chr>    <dbl>      <dbl>      <dbl>      <dbl> <dbl> <dbl>     <dbl>         <dbl>
    ## 1 V03A03    275.          0          0          0  250.     0         0             0
    ## 2 V20A03    275.          0          0          0  250.     0         0             0

``` r
results_mercuryemis <- tibble()

for (subdir in subdirs) {
  file_path <- file.path(subdir, "GEOSChem.MercuryEmis.alltime_m.nc4")
  
  if (!file.exists(file_path)) next
  
  merc_emiss <- nc_open(file_path)
  
  # Temps et poids mensuels
  time_raw <- ncvar_get(merc_emiss, "time")
  origin_time <- as.POSIXct("2015-01-01 00:00:00", tz = "UTC")
  time <- origin_time + time_raw * 60
  weights <- days_in_month(time) / sum(days_in_month(time))
  weighted_avg <- function(x) sum(x * weights)
  
  # Emissions
  hg0land <- ncvar_get(merc_emiss, "EmisHg0land")  # [lon, lat, time]
  hg0soil <- ncvar_get(merc_emiss, "EmisHg0soil")

  # Calculs
  hg0land_mean <- apply(hg0land, c(1, 2), weighted_avg)
  hg0soil_mean <- apply(hg0soil, c(1, 2), weighted_avg)
  
  total_hg0land <- sum(hg0land_mean) * s_in_yr * kg_Mg
  total_hg0soil <- sum(hg0soil_mean) * s_in_yr * kg_Mg
  
  # Ajouter au tableau
  results_mercuryemis <- bind_rows(results_mercuryemis, tibble(
    case = basename(subdir),
    hg0land = total_hg0land,
    hg0soil = total_hg0soil
  ))
  
  nc_close(merc_emiss)
}

print(results_mercuryemis)
```

    ## # A tibble: 2 × 3
    ##   case   hg0land hg0soil
    ##   <chr>    <dbl>   <dbl>
    ## 1 V03A03    32.0    840.
    ## 2 V20A03    32.2    840.

``` r
subdirs <- list.dirs(chem_QIU, full.names = TRUE, recursive = FALSE)
subdirs <- subdirs[grepl("^V", basename(subdirs))]

# Initialiser les résultats
results_drydep <- tibble()

# Fonction de moyenne pondérée
weighted_avg <- function(x, weights) sum(x * weights)

# Boucle sur les cas
for (subdir in subdirs) {
  file_path <- file.path(subdir, "GEOSChem.DryDep.2015_m.nc4")
  if (!file.exists(file_path)) next
  
  drydep_nc <- nc_open(file_path)
  
  # Temps
  time_raw <- ncvar_get(drydep_nc, "time")
  origin_time <- as.POSIXct("2015-01-01 00:00:00", tz = "UTC")
  time <- origin_time + time_raw * 60
  days_in_months <- days_in_month(time)
  weights <- days_in_months / sum(days_in_months)
  
  # Extraction et traitement des 3 espèces
  get_drydep_sum <- function(varname) {
    var <- ncvar_get(drydep_nc, varname)
    var <- apply(var, c(1, 2), weighted_avg, weights = weights)
    var <- var * unit_conv_dd * kg_Mg
    sum(var, na.rm = TRUE)
  }
  
  DryDep_Hg0_hist <- get_drydep_sum("DryDep_Hg0")
  DryDep_Hg2_hist <- get_drydep_sum("DryDep_Hg2")
  DryDep_HgP_hist <- get_drydep_sum("DryDep_HgP")
  
  # Stocker
  results_drydep <- bind_rows(results_drydep, tibble(
    case = basename(subdir),
    drydep_Hg0_tonnes = DryDep_Hg0_hist,
    drydep_Hg2_tonnes = DryDep_Hg2_hist,
    drydep_HgP_tonnes = DryDep_HgP_hist
  ))
  
  nc_close(drydep_nc)
}

# Afficher les résultats
print(results_drydep)
```

    ## # A tibble: 2 × 4
    ##   case   drydep_Hg0_tonnes drydep_Hg2_tonnes drydep_HgP_tonnes
    ##   <chr>              <dbl>             <dbl>             <dbl>
    ## 1 V03A03             1448.              539.              51.1
    ## 2 V20A03             1417.              542.              51.4

``` r
subdirs <- list.dirs(chem_QIU, full.names = TRUE, recursive = FALSE)
subdirs <- subdirs[grepl("^V", basename(subdirs))]

# Initialiser le tableau des résultats
results_wetloss <- tibble()

# Fonction de moyenne pondérée
weighted_avg <- function(x, weights) sum(x * weights)

# Boucle sur les cas
for (subdir in subdirs) {
  file_path <- file.path(subdir, "GEOSChem.WetLossTot_HgSum.2015_m.nc4")
  if (!file.exists(file_path)) next
  
  wetloss_nc <- nc_open(file_path)
  
  # Temps
  time_raw <- ncvar_get(wetloss_nc, "time")
  origin_time <- as.POSIXct("2015-01-01 00:00:00", tz = "UTC")
  time <- origin_time + time_raw * 60
  days_in_months <- days_in_month(time)
  weights <- days_in_months / sum(days_in_months)
  
  # Extraction et traitement
  wetloss_array <- ncvar_get(wetloss_nc, "WetLossTot_HgSum")
  wetloss_array <- apply(wetloss_array, c(1, 2), weighted_avg, weights = weights)
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

    ## # A tibble: 2 × 2
    ##   case   wetloss_Hg_total_tonnes
    ##   <chr>                    <dbl>
    ## 1 V03A03                   1372.
    ## 2 V20A03                   1381.

``` r
# Fusionner les tables selon la colonne commune "case"
results_all <- results_wetloss %>%
  full_join(results_drydep, by = "case") %>%
  full_join(results_mercuryemis, by = "case") %>%
  full_join(results, by = "case")

# Écrire dans un fichier Excel
write_xlsx(results_all, "C:/Users/colom/Desktop/STAGE/data/clean_mod_data/14_3_1/QIU_resultats_mercury.xlsx")


# Fusionner les tables selon la colonne commune "case"
results_mass_all <- results_mass_Hg_total  %>%
  full_join(results_mass_Hg0, by = "case")
# Écrire dans un fichier Excel
write_xlsx(results_mass_all, "C:/Users/colom/Desktop/STAGE/data/clean_mod_data/14_3_1/QIU_resultats_mass_mercury.xlsx")
```
