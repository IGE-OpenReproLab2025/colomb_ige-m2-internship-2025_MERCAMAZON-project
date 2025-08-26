Emis map 2003 2020
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
library(terra)
library(writexl)
library(rnaturalearth)
library(viridis)
library(patchwork)
library(ggtext)
library(sf)
library(rnaturalearthdata)
```

``` r
chem_emis<-"C:/Users/colom/Desktop/STAGE/data/clean_mod_data/14_3_1/Emis_EDGAR"
edgar_2003<-nc_open(file.path(chem_emis, "v8.1_FT2022_TOX_Hg_G_2003_TOTALS_flx.nc"))
edgar_2020<-nc_open(file.path(chem_emis, "v8.1_FT2022_TOX_Hg_G_2020_TOTALS_flx.nc"))
print(edgar_2003)
```

    ## File C:/Users/colom/Desktop/STAGE/data/clean_mod_data/14_3_1/Emis_EDGAR/v8.1_FT2022_TOX_Hg_G_2003_TOTALS_flx.nc (NC_FORMAT_NETCDF4):
    ## 
    ##      1 variables (excluding dimension variables):
    ##         float emi_hg_g[lon,lat,time]   (Chunking: [3600,1800,1])  
    ##             long_name: TOTALS
    ##             units: kg m-2 s-1
    ##             _FillValue: NaN
    ##             missing_value: NaN
    ##             global_total: 1.00kt
    ##             substance: Hg_G
    ##             year: 2003
    ##             release: EDGARv8.1
    ##             description: TOTALS
    ## 
    ##      3 dimensions:
    ##         time  Size:1   *** is unlimited *** 
    ##             axis: T
    ##             calendar: proleptic_gregorian
    ##             standard_name: time
    ##             units: days since 2003-1-15 00:00:00
    ##         lon  Size:3600 
    ##             standard_name: longitude
    ##             long_name: longitude
    ##             units: degrees_east
    ##             axis: X
    ##         lat  Size:1800 
    ##             standard_name: latitude
    ##             long_name: latitude
    ##             units: degrees_north
    ##             axis: Y
    ## 
    ##     13 global attributes:
    ##         CDI: Climate Data Interface version 2.1.1 (https://mpimet.mpg.de/cdi)
    ##         Conventions: CF-1.6
    ##         source: https://edgar.jrc.ec.europa.eu/dataset_tox81
    ##         institution: European Commission, Joint Research Centre
    ##         description: TOTALS
    ##         how_to_cite: https://edgar.jrc.ec.europa.eu/dataset_tox81#howtocite
    ##         copyright_notice: https://edgar.jrc.ec.europa.eu/dataset_tox81#conditions
    ##         contacts: https://edgar.jrc.ec.europa.eu/dataset_tox81#info JRC-EDGAR@ec.europa.eu
    ##         units: kg m-2 s-1
    ##         title: EDGARv8.1 Emission Dataset
    ##         NCO: netCDF Operators version 5.1.4 (Homepage = http://nco.sf.net, Code = http://github.com/nco/nco)
    ##         CDO: Climate Data Operators version 2.1.1 (https://mpimet.mpg.de/cdo)
    ##         history: Tue Jul 22 12:35:09 2025: ncrename -v fluxes,emi_hg_g v8.1_FT2022_TOX_Hg_G_2003_TOTALS_flx.nc
    ## Fri Sep  6 15:20:25 2024: ncap2 -s time=time*-1*0 temp0_v8.1_FT2022_TOX_Hg_G_2003_TOTALS_flx.nc temp1_v8.1_FT2022_TOX_Hg_G_2003_TOTALS_flx.nc
    ## Fri Sep 06 15:20:25 2024: cdo setreftime,2003-01-15,00:00:00 v8.1_FT2022_TOX_Hg_G_2003_TOTALS_flx.nc temp0_v8.1_FT2022_TOX_Hg_G_2003_TOTALS_flx.nc
    ## Fri Sep  6 13:55:46 2024: ncatted -a title,global,o,c,EDGARv8.1 Emission Dataset /workdir/chianti/angoth/ExtData/HEMCO/MERCURY/v2022-06/EDGAR/v8.1_FT2022_TOX_Hg_G_2003_TOTALS_flx.nc0
    ## Fri Sep 06 13:55:46 2024: cdo setreftime,1970-01-15,00:00:00 /workdir/chianti/angoth/ExtData/HEMCO/MERCURY/v2022-06/EDGAR/v8.1_FT2022_TOX_Hg_G_2003_TOTALS_flx.nc /workdir/chianti/angoth/ExtData/HEMCO/MERCURY/v2022-06/EDGAR/v8.1_FT2022_TOX_Hg_G_2003_TOTALS_flx.nc0

``` r
# === 2. Extraire les variables ===
lon <- ncvar_get(edgar_2003, "lon")
lat <- ncvar_get(edgar_2003, "lat")
emi <- ncvar_get(edgar_2003, "emi_hg_g")
mean(emi)
```

    ## [1] 4.568227e-17

``` r
emi<-emi * 31536000 #kg/m2/an
emi_matrix <- t(emi)  # maintenant [lat, lon]
r <- rast(emi_matrix,
          extent = c(min(lon), max(lon), min(lat), max(lat)),
          crs = "EPSG:4326")  # système de coordonnées géographiques
r <- flip(r, direction = "vertical")
df_2003 <- as.data.frame(r, xy = TRUE)
names(df_2003)[3] <- "emi_hg"





df_2003$emi_mg <- df_2003$emi_hg * 1e6
df_2003$emi_mg[df_2003$emi_mg < 2e-9] <- NA  # on supprime visuellement tout < 2e-9

# --- Nouvelles classes (sans la première) ---
breaks  <- c(2e-8, 2e-7, 2e-6, 2e-5, 5e-5, 2e-4, 5e-4, 2e-3, 5e-2, Inf)
labels  <- c(
  "2e-8 – 2e-7", "2e-7 – 2e-6", "2e-6 – 2e-5", "2e-5 – 5e-5",
  "5e-5 – 2e-4", "2e-4 – 5e-4", "5e-4 – 2e-3", "2e-3 – 5e-2", "> 5e-2"
)
colors  <- c(
  "#005F73", "#0A9396", "#94D2BD", "#E9D8A6",
  "#EE9B00", "#CA6702", "#BB3E03", "#AE2012", "#9B2226"
)

df_2003$emi_class <- cut(df_2003$emi_mg, breaks = breaks, labels = labels, include.lowest = TRUE)
df_2003<-na.omit(df_2003)
# --- Plot ---
plot_2003<- ggplot(df_2003, aes(x = x, y = y, fill = emi_class)) +
  geom_raster() +
  scale_fill_manual(
    name = "Hg (mg/m²/an)",
    values = colors,
    na.value = "transparent",
    guide = guide_legend(reverse = TRUE)  

  ) +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Émissions de mercure EDGAR v8.1 - 2003",
    x = "Longitude", y = "Latitude"
  )+
    borders("world", colour = "#403d39") +
    theme_minimal()

print(plot_2003)
```

![](Map_Emis_2003_2020_EDGAR_files/figure-gfm/plot%202003-1.png)<!-- -->

``` r
# === 2. Extraire les variables ===
lon <- ncvar_get(edgar_2020, "lon")
lat <- ncvar_get(edgar_2020, "lat")
emi <- ncvar_get(edgar_2020, "emi_hg_g")
dim(emi)
```

    ## [1] 3600 1800

``` r
emi<-emi * 31536000
emi_matrix <- t(emi)  # maintenant [lat, lon]
r <- rast(emi_matrix,
          extent = c(min(lon), max(lon), min(lat), max(lat)),
          crs = "EPSG:4326")  # système de coordonnées géographiques
r <- flip(r, direction = "vertical")
df_2020 <- as.data.frame(r, xy = TRUE)
names(df_2020)[3] <- "emi_hg"



df_2020$emi_mg <- df_2020$emi_hg * 1e6
df_2020$emi_mg[df_2020$emi_mg < 2e-9] <- NA  # on supprime visuellement tout < 2e-9

# --- Nouvelles classes (sans la première) ---
breaks  <- c(2e-8, 2e-7, 2e-6, 2e-5, 5e-5, 2e-4, 5e-4, 2e-3, 5e-2, Inf)
labels  <- c(
  "2e-8 – 2e-7", "2e-7 – 2e-6", "2e-6 – 2e-5", "2e-5 – 5e-5",
  "5e-5 – 2e-4", "2e-4 – 5e-4", "5e-4 – 2e-3", "2e-3 – 5e-2", "> 5e-2"
)
colors  <- c(
  "#005F73", "#0A9396", "#94D2BD", "#E9D8A6",
  "#EE9B00", "#CA6702", "#BB3E03", "#AE2012", "#9B2226"
)

df_2020$emi_class <- cut(df_2020$emi_mg, breaks = breaks, labels = labels, include.lowest = TRUE)
df_2020<-na.omit(df_2020)



labels <- c(
  expression("2" %*% 10^-8 ~ "-" ~ "2" %*% 10^-7),
  expression("2" %*% 10^-7 ~ "-" ~ "2" %*% 10^-6),
  expression("2" %*% 10^-6 ~ "-" ~ "2" %*% 10^-5),
  expression("2" %*% 10^-5 ~ "-" ~ "5" %*% 10^-5),
  expression("5" %*% 10^-5 ~ "-" ~ "2" %*% 10^-4),
  expression("2" %*% 10^-4 ~ "-" ~ "5" %*% 10^-4),
  expression("5" %*% 10^-4 ~ "-" ~ "2" %*% 10^-3),
  expression("2" %*% 10^-3 ~ "-" ~ "5" %*% 10^-2),
  expression(">" ~ "5" %*% 10^-2)
)

df_2020$emi_class <- cut(df_2020$emi_mg, breaks = breaks, labels = labels, include.lowest = TRUE)
df_2020<-na.omit(df_2020)
# --- Plot ---
plot_2020<-ggplot(df_2020, aes(x = x, y = y, fill = emi_class)) +
  geom_raster() +
 scale_fill_manual(
  name = "Hg(0) (mg/m²/year)",
  values = colors,
  labels = labels,  # ← ici
  na.value = "transparent",
  guide = guide_legend(reverse = TRUE)
) +
  coord_equal() +
  theme_minimal() +
  labs(
  title = "Annual Hg(0) Emissions in 2020",
  subtitle = "Dataset : EDGAR v8.1",
  x = "Longitude", y = "Latitude"
)+
    borders("world", colour = "#403d39", size = 0.2) +
    theme_minimal() +
   theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text = element_text(size = 10),
    axis.ticks = element_line(color = "black", size = 0.3),  # Ajout des petits traits
    axis.ticks.length = unit(3, "pt"),  # Longueur des petits traits
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.line = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  )

plot_2020
```

![](Map_Emis_2003_2020_EDGAR_files/figure-gfm/plot%202020-1.png)<!-- -->

``` r
# Chemin de sauvegarde (modifie si nécessaire)
chemin_pdf <- "C:/Users/colom/Desktop/STAGE/data/clean_mod_data/plot/EDGAR_emi_2020.pdf"

# Enregistrement du plot
ggsave(
  filename = chemin_pdf,
  plot = plot_2020,
  device = "pdf",
  width = 10, height = 6, units = "in"
)
```

``` r
plot_2003 + plot_2020 + plot_layout(ncol=2)
```

![](Map_Emis_2003_2020_EDGAR_files/figure-gfm/plot%20total-1.png)<!-- -->

\#markdown_chemin
\<-“C:/Users/colom/Desktop/STAGE/Markdown/map_wetdep_21_07_2025”
rmarkdown::render(file.path(markdown_chemin , “map_wetdep.Rmd”),
output_format = “github_document”)
