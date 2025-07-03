# Evol_mass_mercury_2014_2015_27_06_2025

**Auteur :** Martin Colomb  
**Ref :** Hélène Angot  
**Projet :** MERCAMAZON

## Objectifs

Ce code a pour objectifs de :
- Charger les données climatiques et chimiques issues de fichiers NetCDF (.nc4) produits par GEOS-Chem
- Calculer la masse totale de mercure (en tonnes) pour :
  - l’ensemble des espèces Hg (Hg_total)
  - Hg⁰ seul
  - Hg² (calculé comme Hg_total - Hg⁰)
- Visualiser l’évolution mensuelle des masses de Hg, Hg⁰, et Hg² dans l’atmosphère pour toutes les simulations

## Inputs

Ce code a comme entrées : 
- 	Concentrations chimiques : `GEOSChem.SpeciesConc.YYYYMMDD.nc4`
- 	Données météo : `GEOSChem.StateMet.2014_2015_m.nc4`
- 	Répertoires : `HIST_*`

## Outputs

Ce code a comme sorties : 
- Plots de l'évolution de Hg, Hg⁰, et Hg² dans l’atmosphère




