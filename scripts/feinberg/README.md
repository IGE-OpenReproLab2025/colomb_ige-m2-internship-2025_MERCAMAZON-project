A. Feinberg, Sep 2023
arifeinberg@gmail.com

Analysis and plotting scripts for the publication: Feinberg et al. : Deforestation as an anthropogenic driver of mercury pollution


The scripts related to plotting figures include:
Fig1_deforest_EFs_obs_constraints.py - Plot reproducing Figure 1: comparing emission factors modelled by GEOS-Chem with observational derived values
Fig2_FigS7_LUC_global_factors_emiss_v3.py - Plot reproducing Figure 2 and Figure S7: global deforestation emissions on the country-level under different time horizons
Fig3a_b_Amazon_terr_flux_pres.py - Plot reproducing Figure 3a and b: surface-atmosphere exchange in Amazon
Fig3c_Amazon_mass_balance_sims_v3.py - Plot reproducing Figure 3c: mass balance of Hg in Amazon
Fig4_Reforestation_terr_flux_pres_v2 - Plot reproducing Figure 4: surface-atmosphere exchange differences after reforestation 
Fig5_paper_emisspolicy_defo.py - Plot reproducing Figure 5: comparing land use policies to primary emissions reductions policies
FigS1_GC_soil_emission_param_tune.py - Plot reproducing Figure S1: tuning soil emissions parametrization to match observational constraints
FigS2_AmazonSI_soil_emissions_diff.py - Plot reproducing Figure S2: difference between new and old soil emissions parametrizations
FigS3_paper_soilconc_Am.py - Plot reproducing Figure S3: soil concentrations in the literature from Amazon forested and deforested sites
FigS4_Area_deforest_flux_linear_Amazon.py - Plot reproducing Figure S4: assessing the linearity of deforestation fluxes to deforested area
FigS5_plot_masks_DFR_v2.py - Plot reproducing Figure S5: plotting the masks of regions used to calculate emission factors
FigS6_bar_defo_time_horizon.py - Plot reproducing Figure S6: bar plot comparing global and country-level deforestation emissions depending on considered time horizon
FigS8_AmazonSI_scen.py - Plot reproducing Figure S8: plotting Soares-Filho06 scenario maps 
FigS9_Amazon_LAI_maps.py - Plot reproducing Figure S9: Amazon LAI maps for scenarios HIST, BAU, GOV, and SAV
FigS10_RFR_LAI_maps.py - Plot reproducing Figure S10: Global LAI maps for scenarios HIST, RFR
FigS11_Hg0_diff_maps.py - Plot reproducing Figure S11: surface Hg0 concentrations in simulations

The scripts related to producing input data for simulations include:
deforest_global_nc.py - Script producing input data (LAI + land cover) for deforestation simulations for EF calculations (DFR), distinguishing 8 global regions
Griscom_reforest_world_nc.py - Script producing input data (LAI + land cover) for RFR simulation, based on Griscom et al. (2017) reforestation potential map
savannize_landtype_2003_nc.py - Script producing input data (LAI + land cover) for SAV simulation
Soares_Amazon_bb_nc.py - Script producing biomass burning emissions data for HIST, BAU, and GOV simulations
Soares_Amazon_landtype_nc.py - Script producing LAI + land cover data for HIST, BAU, and GOV simulations

The scripts related to producing the uncertainty analysis of air-land exchange fluxes using offline models are:
unc_Amazon_bb.py - Assessing uncertainty in biomass burning emissions from the Amazon in different scenarios
unc_drydep_GEOSChem_offline_m_all_land_HIST.py - Uncertainty analysis of Hg0 dry deposition from all land grid cells in HIST scenario (to compare with RFR)
unc_drydep_GEOSChem_offline_m_Am_BAU.py - Uncertainty analysis of Hg0 dry deposition from Amazon grid cells in BAU scenario
unc_drydep_GEOSChem_offline_m_Am_GOV.py - Uncertainty analysis of Hg0 dry deposition from Amazon grid cells in GOV scenario
unc_drydep_GEOSChem_offline_m_Am_HIST.py - Uncertainty analysis of Hg0 dry deposition from Amazon grid cells in HIST scenario (to compare with BAU, GOV, SAV)
unc_drydep_GEOSChem_offline_m_Am_SAV.py - Uncertainty analysis of Hg0 dry deposition from Amazon grid cells in SAV scenario
unc_drydep_GEOSChem_offline_m_DFR.py - Uncertainty analysis of Hg0 dry deposition from all region-level grid cells in global deforestation (DFR) scenarios
unc_drydep_GEOSChem_offline_m_realms_HIST.py - Uncertainty analysis of Hg0 dry deposition from all region-level grid cells in HIST scenario (to compare with DFR)
unc_drydep_GEOSChem_offline_m_RFR.py - Uncertainty analysis of Hg0 dry deposition from all land grid cells in HIST scenario
unc_GC_soil_emission_m_BAU_v2.py - Uncertainty analysis of Hg0 soil emissions from all grid cells in BAU scenario
unc_GC_soil_emission_m_DFR.py - Uncertainty analysis of Hg0 soil emissions from all grid cells in DFR scenario
unc_GC_soil_emission_m_GOV.py - Uncertainty analysis of Hg0 soil emissions from all grid cells in GOV scenario
unc_GC_soil_emission_m_HIST.py - Uncertainty analysis of Hg0 soil emissions from all grid cells in HIST scenario
unc_GC_soil_emission_m_RFR.py - Uncertainty analysis of Hg0 soil emissions from all grid cells in RFR scenario
unc_GC_soil_emission_m_SAV.py - Uncertainty analysis of Hg0 soil emissions from all grid cells in SAV scenario
unc_GC_soil_emission_param_tune_auto_v4.py - Producing 100 soil emissions parametrizations that fall within bounds of observed constraints
unc_LHS_parameter_space.py - Producing the Latin Hypercube Sampling of parameter values to be used in uncertainty analysis
unc_RFR_LAI_nc.py - Producing 100 LAI maps for the uncertainty simulations in RFR, to be run before unc_GC_soil_emission_m_RFR.py and unc_drydep_GEOSChem_offline_m_RFR.py

Other scripts include:
GC_soil_emission_h.py - Script to calculate offline soil Hg0 emissions for a full year using hourly meteorological input data
drydep_functions.py - Script containing functions useful for analyzing land cover datasets
helper_functions.py - Script containing miscellaneous functions for analyzing GEOS-Chem output data

Processed data required for running the scripts can be found here in the directory: misc_Data/

Several datasets are already available online (noted in header of scripts) and are therefore not included in this upload. These include:
LUH2 harmonized land use data (v2h release) for use in CMIP6 simulations can be found here: https://luh.umd.edu/data.shtml 
Soares-Filho et al. (2006) Amazon scenarios can be found here: https://doi.org/10.3334/ORNLDAAC/1153
GEOS-Chem meteorological input data can be found here: http://geoschemdata.wustl.edu/ExtData/
Previous GEOS-Chem simulation data from Feinberg et al. (2022) can be found here: https://doi.org/10.7910/DVN/R7NRNK
