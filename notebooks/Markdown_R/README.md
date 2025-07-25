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

L'objectif de ces changements de valeur de β ont pour objectif d'affiner aux mieux les valeurs de concentrations en Hg<sup>0</sup> et Hg<sup>II</sup>

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

---
## Versions du logiciel et des librairies

Logiciel R / Rstudio : R-4.4.2

| Package     | Version  |
|-------------|----------|
| dplyr       | 1.1.4    |
| geosphere   | 1.5.20   |
| ggplot2     | 3.5.1    |
| lubridate   | 1.9.4    |
| maps        | 3.4.2.1  |
| ncdf4       | 1.23     |
| patchwork   | 1.3.0    |
| stringr     | 1.5.1    |
| tibble      | 3.2.1    |
| viridis     | 0.6.5    |
| writexl     | 1.5.4    |


---
### Stratégies de conservation du code

Dans le cadre de ce stage, une attention particulière a été portée à la conservation du code, tant à court qu’à long terme.
À court terme, le code est organisé, commenté et sauvegardé régulièrement afin de faciliter les tests, les ajustements et les retours pendant la période de développement.
À long terme, l’ensemble du travail est documenté (README, commentaires, structure des fichiers) et conservé dans un dépôt versionné (par exemple Git), permettant sa reprise ou son évolution par d’autres membres de l’équipe ou pour de futurs travaux.

---
### Stratégies de conservation des données

À court terme les données utilisées et produites durant le stage sont conservées de manière organisée et sécurisée sur Zenodo (Actuellement sur [Google Drive](https://drive.google.com/file/d/1nQpKIO5ppbKZYkHsfV--BPHWDcm6tgxA/view?usp=sharing) car problème pour les importer sur Zenodo.).
À long terme, les fichiers importants sont stockés dans un espace dédié (répertoire partagé, serveur de l’équipe, etc.) avec une structure claire et une documentation associée, afin d’en permettre la réutilisation ou la vérification dans le cadre de travaux futurs.
