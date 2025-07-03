# Mercury_mass_23_06_2025

**Auteur :** Martin Colomb  
**Ref :** Hélène Angot  
**Projet :** MERCAMAZON

## Objectifs

Ce code a pour objectifs de :
- Charger les données climatiques et chimiques issues de fichiers NetCDF (.nc4) produits par GEOS-Chem
- Calcul de la masse de Hg⁰ dans l’atmosphère
- Traitement des fichiers mensuels pour extraire la masse mensuelle de Hg⁰ sur l’année 2015
- Extraction des émissions mensuelles
- Extraction des émissions de surfaces
- Calcul des dépôts secs totaux
- Calcul des dépôts humimdes totaux
- Export des résultats dans un fichier Excel

  
## Inputs

Ce code a comme entrées : 
- 	Concentrations chimiques : `GEOSChem.SpeciesConc.2015_m.nc4`
- 	Fichiers d'émissions : `HEMCO_diagnostic_2015_m.nc4` & `GEOSChem.MercuryEmis.alltime_m.nc4`
- 	Données météo : `GEOSChem.StateMet.2015_m.nc4`

## Outputs

Ce code a comme sorties : 
- Plots de l'évolution mensuelle de la masse Hg dans l’atmosphère
- Fichier récapitulatif des émissions : `resultats_mercury.xlsx`
- Fichier récapitulatif des masses, dépots secs, dépots humides : `resultats_mass_mercury.xlsx`





