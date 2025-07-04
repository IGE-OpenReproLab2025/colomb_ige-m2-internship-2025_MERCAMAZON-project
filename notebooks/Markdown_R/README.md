# Informations concernant les codes relatifs aux projet MERCAMAZON Rôle des forêts dans la régulation du mercure atmosphérique.
## Peut-on compenser les émissions par la prévention de la déforestation ?

**Auteur :**  Martin Colomb  
**Ref :**  Hélène Angot  
---
## Contexte

Tout comme pour le CO₂, la réduction des émissions de mercure à la source reste la priorité, mais des stratégies complémentaires comme la préservation des forêts jouent un rôle clé. La végétation agissant comme un puits majeur pour ces deux composés : les plantes captent le Hg⁰ via les feuilles et le stockent dans les sols, un mécanisme analogue à la séquestration du carbone [(Feinberg et al., 2022)](https://pubs.rsc.org/en/content/articlelanding/2022/em/d2em00032f). Or, la déforestation libère le mercure préalablement accumulé et réduit cette capacité d’accumulation. Dans ce contexte, mon stage vise à évaluer dans quelle mesure la prévention de la déforestation pourrait compenser les émissions globales de mercure, en maintenant ce réservoir naturel.

---
## Objectif de ce répertoire

Ces codes utilisent des outputs issus du modèle GEOS-Chem (V14.3.1) afin d'observer et calculer des concentrations, masses et dépôts de mercure.

---
## Instructions particulières

Il est nécessaire de remplacer tous les chemins dans chaque code pour pouvoir les faire tourner. Chaque output de modèle est stocké dans un dossier correspondant à la simulation


---
## Scénarios

Les différents scénarios simulés contiennent les mêmes données d'entrées.
- Scénario de référence : `HIST_V3`
- Scénario de référence avec une année (ou 2) de spinup supplémentaire : `HIST_V3_3y`,  `HIST_V3_4y`

Pour les autres simulations, un paramètre a été changé par rapport aux précédentes, `β` (dans les précédentes, `β = 0,01`).

La fréquence de photolyse de cette réaction est paramétrée comme suit :

**J<sub>HgIIP(org)</sub> = β × J<sub>NO₂</sub>**
où  *J<sub>NO₂</sub>* est la fréquence locale de photolyse du NO₂, et le facteur d’échelle **β** est ajusté pour correspondre à la moyenne globale des observations de surface de Hg⁰.
[(Shah et al. 2021)](https://pubs.acs.org/doi/pdf/10.1021/acs.est.1c03160)

Ainsi, B10 représente une augmentation de 10% telle que `β=0,011`, B20 de 20% `β=0,012` etc.

---


## 1. Comp_conc_beta_ref_23_06_2025
[Comp_conc_beta_ref_23_06_2025.md](Comp_conc_beta_ref_23_06_2025/Comp_conc_beta_ref_23_06_2025.md) : 

Ce code a pour objectifs de :
- Charger les données chimiques issues de fichiers NetCDF (.nc4) produits par GEOS-Chem
- Comparer données modélisées et de référence (séparation hémisphère nord - hémisphère sud) pour le scénario B20 (exemple)
- Comparer données modélisées et de référence pour tous les scénarios


---

## 2. Evol_mass_mercury_2014_2015_27_06_2025
[Evol_mass_mercury_2014_2015_27_06_2025.md](Evol_mass_mercury_2014_2015_27_06_2025/Evol_mass_mercury_2014_2015_27_06_2025.md) : 

Ce code a pour objectifs de :
- Charger les données climatiques et chimiques issues de fichiers NetCDF (.nc4) produits par GEOS-Chem
- Calculer la masse totale de mercure (en tonnes) pour :
  - l’ensemble des espèces Hg (Hg_total)
  - Hg⁰ seul
  - Hg² (calculé comme Hg_total - Hg⁰)
- Visualiser l’évolution mensuelle des masses de Hg, Hg⁰, et Hg² dans l’atmosphère pour toutes les simulations


---

## 3. Mercury_mass_16_06_2025.md
[Mercury_mass_23_06_2025.md](Mercury_mass_23_06_2025/Mercury_mass_23_06_2025.md) : 

Ce code a pour objectifs de :
- Charger les données climatiques et chimiques issues de fichiers NetCDF (.nc4) produits par GEOS-Chem
- Calcul de la masse de Hg⁰ dans l’atmosphère
- Traitement des fichiers mensuels pour extraire la masse mensuelle de Hg⁰ sur l’année 2015
- Extraction des émissions mensuelles
- Extraction des émissions de surfaces
- Calcul des dépôts secs totaux
- Calcul des dépôts humimdes totaux
- Export des résultats dans un fichier Excel 


---

## 4. Map_mod_vs_ref_04_07_2025
[Map_mod_vs_ref_04_07_2025.md](Map_mod_vs_ref_04_07_2025/Map_mod_vs_ref_04_07_2025.md) : 

Ce code a pour objectifs de :
- Lire les fichiers NetCDF et calculer les moyennes pondérées annuelles
- Visualiser les concentrations modélisées
- Comparer visuellement les simulations aux observations


