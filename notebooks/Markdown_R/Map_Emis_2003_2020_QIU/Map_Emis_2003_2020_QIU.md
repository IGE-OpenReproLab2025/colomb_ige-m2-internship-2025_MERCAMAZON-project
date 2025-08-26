Map_Emis_2003_2020_QIU
================
Martin
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
```

``` r
chem_emis<-"C:/Users/colom/Desktop/STAGE/data/clean_mod_data/emis/qiu/Hg0_nc"
qiu_2003<-nc_open(file.path(chem_emis, "Hg0_emis_2003.nc"))
qiu_2020<-nc_open(file.path(chem_emis, "Hg0_emis_2020.nc"))
print(qiu_2003)
```

    ## File C:/Users/colom/Desktop/STAGE/data/clean_mod_data/emis/qiu/Hg0_nc/Hg0_emis_2003.nc (NC_FORMAT_NETCDF4):
    ## 
    ##      2 variables (excluding dimension variables):
    ##         float emi_hg_g[longitude,latitude]   (Contiguous storage)  
    ##             units: kg/m2/s
    ##             _FillValue: -1.17549402418441e+38
    ##             grid_mapping: crs
    ##         int crs[]   (Contiguous storage)  
    ##             crs_wkt: GEOGCRS["WGS 84",    ENSEMBLE["World Geodetic System 1984 ensemble",        MEMBER["World Geodetic System 1984 (Transit)"],        MEMBER["World Geodetic System 1984 (G730)"],        MEMBER["World Geodetic System 1984 (G873)"],        MEMBER["World Geodetic System 1984 (G1150)"],        MEMBER["World Geodetic System 1984 (G1674)"],        MEMBER["World Geodetic System 1984 (G1762)"],        MEMBER["World Geodetic System 1984 (G2139)"],        MEMBER["World Geodetic System 1984 (G2296)"],        ELLIPSOID["WGS 84",6378137,298.257223563,            LENGTHUNIT["metre",1]],        ENSEMBLEACCURACY[2.0]],    PRIMEM["Greenwich",0,        ANGLEUNIT["degree",0.0174532925199433]],    CS[ellipsoidal,2],        AXIS["geodetic latitude (Lat)",north,            ORDER[1],            ANGLEUNIT["degree",0.0174532925199433]],        AXIS["geodetic longitude (Lon)",east,            ORDER[2],            ANGLEUNIT["degree",0.0174532925199433]],    USAGE[        SCOPE["Horizontal component of 3D system."],        AREA["World."],        BBOX[-90,-180,90,180]],    ID["EPSG",4326]]
    ##             spatial_ref: GEOGCRS["WGS 84",    ENSEMBLE["World Geodetic System 1984 ensemble",        MEMBER["World Geodetic System 1984 (Transit)"],        MEMBER["World Geodetic System 1984 (G730)"],        MEMBER["World Geodetic System 1984 (G873)"],        MEMBER["World Geodetic System 1984 (G1150)"],        MEMBER["World Geodetic System 1984 (G1674)"],        MEMBER["World Geodetic System 1984 (G1762)"],        MEMBER["World Geodetic System 1984 (G2139)"],        MEMBER["World Geodetic System 1984 (G2296)"],        ELLIPSOID["WGS 84",6378137,298.257223563,            LENGTHUNIT["metre",1]],        ENSEMBLEACCURACY[2.0]],    PRIMEM["Greenwich",0,        ANGLEUNIT["degree",0.0174532925199433]],    CS[ellipsoidal,2],        AXIS["geodetic latitude (Lat)",north,            ORDER[1],            ANGLEUNIT["degree",0.0174532925199433]],        AXIS["geodetic longitude (Lon)",east,            ORDER[2],            ANGLEUNIT["degree",0.0174532925199433]],    USAGE[        SCOPE["Horizontal component of 3D system."],        AREA["World."],        BBOX[-90,-180,90,180]],    ID["EPSG",4326]]
    ##             proj4: +proj=longlat +datum=WGS84 +no_defs
    ##             code: EPSG:4326
    ##             geotransform: -180 0.1000000000000000055511 0 90 0 -0.1000000000000000055511
    ## 
    ##      2 dimensions:
    ##         longitude  Size:3600 
    ##             units: degrees_east
    ##             long_name: longitude
    ##         latitude  Size:1800 
    ##             units: degrees_north
    ##             long_name: latitude
    ## 
    ##     3 global attributes:
    ##         Conventions: CF-1.4
    ##         created_by: R packages ncdf4 and terra (version 1.8-29)
    ##         created_date: 2025-08-01 12:33:06

``` r
# === 2. Extraire les variables ===
lon <- ncvar_get(qiu_2003, "longitude")
lat <- ncvar_get(qiu_2003, "latitude")
emi <- ncvar_get(qiu_2003, "emi_hg_g") 
mean(emi)
```

    ## [1] 4.035383e-17

``` r
dim(emi)
```

    ## [1] 3600 1800

``` r
emi<-emi * 31536000 #kg/m2/an
emi_matrix <- t(emi)  # maintenant [lat, lon]
r <- rast(emi_matrix,
          extent = c(min(lon), max(lon), min(lat), max(lat)),
          crs = "EPSG:4326")  # système de coordonnées géographiques
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
max(df_2003$emi_mg)
```

    ## [1] 124.2736

``` r
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
    title = "Émissions de mercure qiu - 2003",
    x = "Longitude", y = "Latitude"
  )+
    borders("world", colour = "#403d39") +
    theme_minimal()

print(plot_2003)
```

![](Map_Emis_2003_2020_QIU_files/figure-gfm/plot%202003-1.png)<!-- -->

``` r
# === 2. Extraire les variables ===
lon <- ncvar_get(qiu_2020, "longitude")
lat <- ncvar_get(qiu_2020, "latitude")
emi <- ncvar_get(qiu_2020, "emi_hg_g")
dim(emi)
```

    ## [1] 3600 1800

``` r
emi<-emi * 31536000
emi_matrix <- t(emi)  # maintenant [lat, lon]
r <- rast(emi_matrix,
          extent = c(min(lon), max(lon), min(lat), max(lat)),
          crs = "EPSG:4326")  # système de coordonnées géographiques
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
# --- Plot ---
plot_2020<-ggplot(df_2020, aes(x = x, y = y, fill = emi_class)) +
  geom_raster() +
  scale_fill_manual(
    name = "Hg⁰ (mg/m²/an)",
    values = colors,
    na.value = "transparent",
    guide = guide_legend(reverse = TRUE)  

  ) +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Émissions de Hg⁰ qiu v8.1 - 2020",
    x = "Longitude", y = "Latitude"
  )+
    borders("world", colour = "#403d39", size = 0.2) +
    theme_minimal()

print(plot_2020)
```

![](Map_Emis_2003_2020_QIU_files/figure-gfm/plot%202020-1.png)<!-- -->

``` r
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
  subtitle = "Dataset : QIU ",
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


# Chemin de sauvegarde (modifie si nécessaire)
chemin_pdf <- "C:/Users/colom/Desktop/STAGE/data/clean_mod_data/plot/qiu_emi_2020.pdf"

# Enregistrement du plot
ggsave(
  filename = chemin_pdf,
  plot = plot_2020,
  device = "pdf",
  width = 10, height = 6, units = "in"
)
```

\#markdown_chemin
\<-“C:/Users/colom/Desktop/STAGE/Markdown/map_wetdep_21_07_2025”
rmarkdown::render(file.path(markdown_chemin , “map_wetdep.Rmd”),
output_format = “github_document”)
