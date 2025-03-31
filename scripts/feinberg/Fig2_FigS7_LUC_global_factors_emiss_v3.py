#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 09:15:28 2023
Calculate regional emission factors for Hg Land Use Change, use them to calculate global estimate for deforestation impacts
@author: arifeinberg
"""
#%% import packages
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy.feature as cf
from helper_functions import ds_sel_yr, annual_avg, open_Hg
import xesmf as xe
import pandas as pd
import regionmask
import geopandas as gpd
from mpl_toolkits.axes_grid1 import make_axes_locatable
#%% create masking function for country masks
def mask_square(mask, country, lt_0, lt_1, ln_0, ln_1):
    """Add true values for square in mask, to fix coastline issues
    
    Parameters
    ----------
    lt_0 : float
         min lat
    lt_1 : float
         max lat
    ln_0 : float
         min lon
    ln_1 : float
         max lon
    country : string
         country name
    mask : xarray DataArray
        country mask to fix
    """
    # index of country, will return error if not found
    region_in = np.where(regions_names==country)[0][0]
    
    lat = country_mask.lat.values
    lon = country_mask.lon.values
    
    # indices of coordinates
    ilt_0 = np.where(lat == lt_0)[0][0]
    ilt_1 = np.where(lat == lt_1)[0][0]
    iln_0 = np.where(lon == ln_0)[0][0]
    iln_1 = np.where(lon == ln_1)[0][0]
    
    # set square in country mask to True
    mask[region_in, ilt_0:ilt_1, iln_0:iln_1] =  True
    
    return mask

#%% Get area deforested for each run
realm_names = ['Amazon', 'Afrotrop', 'Indomalay', 'China', 'Neotrop','Palearc','Austral','Nearc']
sim_names = ['0315','0321','0322','0323','0324','0325','0326','0327']
replace_id = [58, 31, 31, 31, 35, 31, 31, 35]# indices of replaced land type
n_realms = len(realm_names)

# load gbox areas
fn_gb = 'misc_Data/gbox_025x025_Olson.nc'
ds_gb = xr.open_dataset(fn_gb)
gbox_area = ds_gb.cell_area # grid_area
gbox_area_v = gbox_area.values # grid_area

# load reference land types 
fn_ols2 = 'misc_Data/Olson_2001_Land_Type_Masks.025x025.generic.nc'
ds_ols2 = xr.open_dataset(fn_ols2)
# save as one np array, first dim is land types
Olson_REF = ds_ols2.to_array().squeeze()

lon = ds_ols2.lon # longitude
lat = ds_ols2.lat # latitude

n_landtypes = 73 # number of Olson land types

# at 2 x 2.5 resolution
fn_ols2_2x25 = 'misc_Data/Olson_2001_Land_Type_Masks.2_25.generic.nc'
ds_ols2_2x25 = xr.open_dataset(fn_ols2_2x25)
# save as one np array, first dim is land types
Olson_REF_2x25 = ds_ols2_2x25.to_array().squeeze()

#%% Calculate land-air exchange in HIST case
sim_i = '0311' # HIST

# Load grid cell area (2x25) for unit conversion of model
fn_gbox = 'misc_Data/GEOSChem_2x25_gboxarea.nc'
ds_gbox = xr.open_dataset(fn_gbox)
gbox_2x25 = ds_gbox.cell_area.values #m2
lat_2x25 = ds_gbox.lat.values
lon_2x25 = ds_gbox.lon.values

fn_emis = '../GC_data/run'+sim_i+'/GEOSChem.MercuryEmis.2015_m.nc4'
ds_emis = xr.open_dataset(fn_emis)

soil_emis_yr = ds_sel_yr(ds_emis, 'EmisHg0soil', 2015) 
land_emis_yr = ds_sel_yr(ds_emis, 'EmisHg0land', 2015) 
geogen_emis_yr = ds_sel_yr(ds_emis, 'EmisHg0geogenic', 2015) 
bb_emis_yr = ds_sel_yr(ds_emis, 'EmisHg0biomass', 2015) 
ant0_emis_yr = ds_sel_yr(ds_emis, 'EmisHg0anthro', 2015) 
ant2_emis_yr = ds_sel_yr(ds_emis, 'EmisHg2HgPanthro', 2015) 

# Convert model data from kg/s to kg/yr for annual average   
s_in_yr = 365.2425 * 24 * 3600 # s in one year

unit_conv = s_in_yr
    
soil_emis = annual_avg(soil_emis_yr) * unit_conv # kg/yr
land_emis = annual_avg(land_emis_yr) * unit_conv # kg/yr
geogen_emis = annual_avg(geogen_emis_yr) * unit_conv # kg/yr
bb_emis = annual_avg(bb_emis_yr) * unit_conv # kg/yr
ant0_emis = annual_avg(ant0_emis_yr) * unit_conv # kg/yr
ant2_emis = annual_avg(ant2_emis_yr) * unit_conv # kg/yr

# Load deposition fluxes
# dry dep
fn_dd = '../GC_data/run'+sim_i+'/GEOSChem.DryDep.2015_m.nc4'
ds_dd = xr.open_dataset(fn_dd)

dd_Hg0_yr = ds_sel_yr(ds_dd, 'DryDep_Hg0', 2015) 
dd_Hg2_yr = ds_sel_yr(ds_dd, 'DryDep_Hg2', 2015) 
dd_HgP_yr = ds_sel_yr(ds_dd, 'DryDep_HgP', 2015) 

# Convert model data from molec/cm2/s to kg/yr for annual average   
s_in_yr = 365.2425 * 24 * 3600 # s in one year
g_kg = 1e-3 # g in kg
cm2_m2 = 1e4 # cm^2 in m^2
MW_Hg = 200.59 # g mol^-1
avo = 6.02e23 # avogadro number molec mol^-1

unit_conv = MW_Hg / avo * g_kg * cm2_m2 * s_in_yr * gbox_2x25 # constant to convert units

dd_Hg0 = annual_avg(dd_Hg0_yr) * unit_conv # kg/yr
dd_Hg2 = annual_avg(dd_Hg2_yr) * unit_conv # kg/yr
dd_HgP = annual_avg(dd_HgP_yr) * unit_conv # kg/yr

# wet dep
fn_dd = '../GC_data/run'+sim_i+'/GEOSChem.WetLossTotal.2015_m.nc4'
ds_dd = xr.open_dataset(fn_dd)

wd_Hg_yr = ds_sel_yr(ds_dd, 'WetLossTot_Hg', 2015) 

# Convert model data from kg/s to kg/yr for annual average   
s_in_yr = 365.2425 * 24 * 3600 # s in one year

unit_conv = s_in_yr

wd_Hg = annual_avg(wd_Hg_yr) * unit_conv # kg/yr

# Calculate totals for fluxes
kg_Mg = 1e3

# emissions
Emiss_tot_HIST = ((soil_emis + land_emis + geogen_emis + bb_emis + ant0_emis + ant2_emis) / kg_Mg).values # Mg/yr
Soil_emiss_HIST = ((soil_emis) / kg_Mg).values # Mg/yr

# deposition
Dep_tot_HIST = ((dd_Hg0 + dd_Hg2 + dd_HgP + wd_Hg ) / kg_Mg).values # Mg/yr

# needed for uncertainty fluxes - change due to other fluxes
Other_HIST = ((land_emis + geogen_emis + bb_emis + ant0_emis + ant2_emis - dd_Hg2 - dd_HgP - wd_Hg) / kg_Mg).values # Mg/yr
#%% loop through land types, calculate area deforested and flux differences
# initialize matrices
Emiss_diff_defo = np.zeros(n_realms)
Dep_diff_defo = np.zeros(n_realms)
LUC_EF_defo = np.zeros(n_realms)
Soil_EF_defo = np.zeros(n_realms)
Other_flux_defo = np.zeros(n_realms)
area_change = np.zeros(n_realms)

for i in range(n_realms):
    realm_i = realm_names[i]
    print(realm_i)
    if (realm_i=='Amazon'):
        fn_realm = '../input_data/SAV/Olson_2001_Land_Type_Masks.025x025.generic_sav58.nc'
        fn_realm_2x25 = '../input_data/SAV/Olson_2001_Land_Type_Masks.2x25.generic_sav58.nc'
    else:
        fn_realm = '../input_data/DFR/Olson_Land_Type_Masks.025x025.DFR_'+realm_i+'.nc'
        fn_realm_2x25 = '../input_data/DFR/Olson_Land_Type_Masks.2x25.DFR_'+realm_i+'.nc'
    # at 0.25 x 0.25 resolution, more accurate area calculation
    ds_realm = xr.open_dataset(fn_realm) 
    Olson_new = ds_realm.to_array().squeeze()
    
    mask_change = Olson_new[replace_id[i], :, :] - Olson_REF[replace_id[i], :, :] # changed area
    
    area_change[i] = (mask_change * gbox_area).sum().values / 1e6 # km2
    print(area_change[i])
    
    # at 2 x 2.5 resolution, for calculated averaged model fluxes    
    ds_realm_2x25 = xr.open_dataset(fn_realm_2x25)
    Olson_new_2x25 = ds_realm_2x25.to_array().squeeze()
    
    mask_change_2x25 = Olson_new_2x25[replace_id[i], :, :] - Olson_REF_2x25[replace_id[i], :, :] # changed area
    # adjust mask so that accounts for all changed area, otherwise missing edge effects
    mask_change_2x25 = xr.where(mask_change_2x25>0, 1, 0) 

    # load simulation data
    sim_i = sim_names[i] #

    fn_emis = '../GC_data/run'+sim_i+'/GEOSChem.MercuryEmis.2015_m.nc4'
    ds_emis = xr.open_dataset(fn_emis)
    
    soil_emis_yr = ds_sel_yr(ds_emis, 'EmisHg0soil', 2015) 
    land_emis_yr = ds_sel_yr(ds_emis, 'EmisHg0land', 2015) 
    geogen_emis_yr = ds_sel_yr(ds_emis, 'EmisHg0geogenic', 2015) 
    bb_emis_yr = ds_sel_yr(ds_emis, 'EmisHg0biomass', 2015) 
    ant0_emis_yr = ds_sel_yr(ds_emis, 'EmisHg0anthro', 2015) 
    ant2_emis_yr = ds_sel_yr(ds_emis, 'EmisHg2HgPanthro', 2015) 
    
    # Convert model data from kg/s to kg/yr for annual average   
    s_in_yr = 365.2425 * 24 * 3600 # s in one year
    
    unit_conv = s_in_yr
        
    soil_emis = annual_avg(soil_emis_yr) * unit_conv # kg/yr
    land_emis = annual_avg(land_emis_yr) * unit_conv # kg/yr
    geogen_emis = annual_avg(geogen_emis_yr) * unit_conv # kg/yr
    bb_emis = annual_avg(bb_emis_yr) * unit_conv # kg/yr
    ant0_emis = annual_avg(ant0_emis_yr) * unit_conv # kg/yr
    ant2_emis = annual_avg(ant2_emis_yr) * unit_conv # kg/yr
    
    # Load deposition fluxes
    # dry dep
    fn_dd = '../GC_data/run'+sim_i+'/GEOSChem.DryDep.2015_m.nc4'
    ds_dd = xr.open_dataset(fn_dd)
    
    dd_Hg0_yr = ds_sel_yr(ds_dd, 'DryDep_Hg0', 2015) 
    dd_Hg2_yr = ds_sel_yr(ds_dd, 'DryDep_Hg2', 2015) 
    dd_HgP_yr = ds_sel_yr(ds_dd, 'DryDep_HgP', 2015) 
    
    # Convert model data from molec/cm2/s to kg/yr for annual average   
    s_in_yr = 365.2425 * 24 * 3600 # s in one year
    g_kg = 1e-3 # g in kg
    cm2_m2 = 1e4 # cm^2 in m^2
    MW_Hg = 200.59 # g mol^-1
    avo = 6.02e23 # avogadro number molec mol^-1
    
    unit_conv = MW_Hg / avo * g_kg * cm2_m2 * s_in_yr * gbox_2x25 # constant to convert units
    
    dd_Hg0 = annual_avg(dd_Hg0_yr) * unit_conv # kg/yr
    dd_Hg2 = annual_avg(dd_Hg2_yr) * unit_conv # kg/yr
    dd_HgP = annual_avg(dd_HgP_yr) * unit_conv # kg/yr
    
    # wet dep
    fn_dd = '../GC_data/run'+sim_i+'/GEOSChem.WetLossTotal.2015_m.nc4'
    ds_dd = xr.open_dataset(fn_dd)
    
    wd_Hg_yr = ds_sel_yr(ds_dd, 'WetLossTot_Hg', 2015) 
    
    # Convert model data from kg/s to kg/yr for annual average   
    s_in_yr = 365.2425 * 24 * 3600 # s in one year
    
    unit_conv = s_in_yr
    
    wd_Hg = annual_avg(wd_Hg_yr) * unit_conv # kg/yr
    
    # Calculate totals for fluxes
    kg_Mg = 1e3
    
    # emissions
    Emiss_tot_i = ((soil_emis + land_emis + geogen_emis + bb_emis + ant0_emis + ant2_emis) / kg_Mg).values # Mg/yr
    Soil_emiss_i = ((soil_emis) / kg_Mg).values # Mg/yr

    # deposition
    Dep_tot_i = ((dd_Hg0 + dd_Hg2 + dd_HgP + wd_Hg ) / kg_Mg).values # Mg/yr
    
    # calculate differences from REF
    Emiss_diff = Emiss_tot_i - Emiss_tot_HIST
    Dep_diff = Dep_tot_i - Dep_tot_HIST

    
    # calculate average only over deforested grid cells
    Emiss_diff_defo[i] = np.sum(mask_change_2x25 * Emiss_diff)
    Dep_diff_defo[i] = np.sum(mask_change_2x25 * Dep_diff)
    Balance = Emiss_diff_defo[i] - Dep_diff_defo[i]
    
    LUC_EF_defo[i] = Balance / area_change[i] # overall Emission factor Mg/yr/km2

    # Factors needed for uncertainty analysis
    Soil_emiss_diff = Soil_emiss_i - Soil_emiss_HIST   
    Soil_emiss_diff_defo = np.sum(mask_change_2x25 * Soil_emiss_diff)
    Soil_EF_defo[i] = (Soil_emiss_diff_defo) / area_change[i] # soil only Emission factor Mg/yr/km2
    
    Other_i = ((land_emis + geogen_emis + bb_emis + ant0_emis + ant2_emis - dd_Hg2 - dd_HgP - wd_Hg) / kg_Mg).values # Mg/yr
    Other_diff = Other_i - Other_HIST
    Other_flux_defo[i] = np.sum(mask_change_2x25 * Other_diff) # total, Mg /yr
    print("Emiss diff (Mg/yr), defo area: " + str(Emiss_diff_defo[i]))
    print("Dep diff (Mg/yr), defo area: " + str(Dep_diff_defo[i]))
    print("Emission factor (Mg/m2/yr), defo area: " + str(LUC_EF_defo[i]))
    print("Soil Emission factor (Mg/m2/yr), defo area: " + str(Soil_EF_defo[i]))
    print("Flux change due to other factors (Mg/yr), defo area: " + str(Other_flux_defo[i]))

# realm_names = ['Amazon', 'Afrotrop', 'Indomalay', 'China', 'Neotrop','Palearc','Austral','Nearc']

# LUC_EF_defo = [7.235329241913561e-05, 1.071806689794155e-05, 2.2819489051113326e-05, 
#                 2.3404428517737673e-05, 7.404479045967743e-06,
#             2.395577187857226e-06, 6.095471571713798e-06, 1.0623755118946325e-05]
# Soil_EF_defo = [1.4104624583103027e-05, 7.965166356949215e-06, 1.2031028933539869e-05,
#                 1.453324446692011e-05,5.359641752377664e-06, 1.368712872976517e-06,
#                 1.7143365857668656e-06, 6.8414553063659995e-06]
# Other_flux_defo = [-9.244362975036687, -0.318027013331827, -1.0131701411212426,
#                     -0.1055723389020468, -0.3064472486730863, -0.12138694100524064, 
#                     -0.10781996357446529, -0.6838366847567511]

#%% Get EFs uncertainty values
# Loop over the realms
# Loop over the realms
Soil_EF_defo_err = np.zeros((2,n_realms)) # error bars for soil EF
LUC_EF_defo_err = np.zeros((2,n_realms)) # error bars for total EF
diff_dep0 = [0] * n_realms # deposition diff
diff_emiss = [0] * n_realms # emissions diff

HIST_emiss =  np.loadtxt("misc_Data/unc_global_emiss_HIST_v2.csv", delimiter=",")

for i in range(n_realms):
    realm = realm_names[i]
    
    # Soil emissions uncertainty 
    if realm=='Amazon':
        realm_emiss = np.loadtxt("misc_Data/unc_global_emiss_SAV.csv", delimiter=",")
    else:
        realm_emiss = np.loadtxt("misc_Data/unc_global_emiss_DFR_" +realm + ".csv", delimiter=",")
    
    # calculate emissions difference due to deforestation perturbation
    diff_emiss[i] = realm_emiss - HIST_emiss
    
    # divide by deforested area to get emission factor
    soil_EF_unc = diff_emiss[i] / area_change[i]
    
    # Get error bounds for soil
    low_soil = np.percentile(soil_EF_unc, 2.5)
    high_soil = np.percentile(soil_EF_unc, 97.5)
    
    # lower bound
    Soil_EF_defo_err[0,i] = Soil_EF_defo[i] - low_soil
    # upper bound
    Soil_EF_defo_err[1,i] = high_soil - Soil_EF_defo[i]

    # Deposition uncertainty 
    realm_dep_hist = np.loadtxt("misc_Data/unc_realm_dep_HIST_" + realm + ".csv", delimiter=",")
    realm_dep_dfr = np.loadtxt("misc_Data/unc_realm_dep_DFR_" + realm + ".csv", delimiter=",")

    # calculate Hg0 deposition difference due to deforestation perturbation
    diff_dep0[i] = realm_dep_dfr - realm_dep_hist

    # divide by deforested area to get emission factor for less deposition (negative sign due to deposition being negative emissions)
    dep_EF_unc = - diff_dep0[i] / area_change[i]

    # Get error bounds for deposition 
    low_dep = np.percentile(dep_EF_unc, 2.5)
    high_dep = np.percentile(dep_EF_unc, 97.5)

    # sum up all differences to get overall low and high flux differences due to deforestation perturbation
    total_diff_low = low_soil + low_dep + Other_flux_defo[i] / area_change[i]
    total_diff_high = high_soil + high_dep + Other_flux_defo[i] / area_change[i]
    
    # lower bound
    LUC_EF_defo_err[0,i] = LUC_EF_defo[i] - total_diff_low
    # upper bound
    LUC_EF_defo_err[1,i] = total_diff_high - LUC_EF_defo[i]

print(LUC_EF_defo_err)
# LUC_EF_defo_err = np.array([[2.77764247e-05, 7.87637585e-06, 7.34510428e-06, 6.03522953e-06,
#                               2.58673852e-06, 2.31978907e-06, 5.26901446e-06, 3.51837723e-06],
#                             [  1.27326075e-04, 1.04295126e-04, 1.90985035e-04, 2.03004306e-04,
#                               5.00345520e-05, 2.10230703e-05, 4.76827582e-05, 5.14893001e-05]])

#%% Prepare region masks, to be used for deforestation flux calculations
# Load masks for biomes and country
realms_masks = ['Nearc', 'Palearc','Afrotrop','Indomalay','Austral','Neotrop', 'China', 'Amazon']
mask_all = [0]*7 # save masks here
# note that mask files have latitude order opposite to GC files 
fn_nea = "misc_Data/Nearctic_mask_025x025.nc"
ds_nea = xr.open_dataset(fn_nea)
mask_all[0] = ds_nea.Nearctic

fn_pal = "misc_Data/Palearctic_mask_025x025.nc"
ds_pal = xr.open_dataset(fn_pal)
mask_all[1] = ds_pal.Palearctic

fn_afr = "misc_Data/Afrotropic_mask_025x025.nc"
ds_afr = xr.open_dataset(fn_afr)
mask_all[2] = ds_afr.Afrotropic

fn_ind = "misc_Data/Indomalayan_mask_025x025.nc"
ds_ind = xr.open_dataset(fn_ind)
mask_all[3] = ds_ind.Indomalayan

fn_aus = "misc_Data/Australasia_mask_025x025.nc"
ds_aus = xr.open_dataset(fn_aus)
fn_oce = "misc_Data/Oceania_mask_025x025.nc"
ds_oce = xr.open_dataset(fn_oce)
mask_all[4] = (ds_aus.Australasia + ds_oce.Oceania)

fn_neo = "misc_Data/Neotropic_mask_025x025.nc"
ds_neo = xr.open_dataset(fn_neo)
mask_all[5] = ds_neo.Neotropic

# load country code
fn_gb = 'misc_Data/staticData_quarterdeg.nc'
ds_gb = xr.open_dataset(fn_gb)
ccode = ds_gb.ccode # country code
mask_all[6] = xr.where(ccode==156., 1, 0).reindex(lat=list(reversed(ds_gb.lat))) # create mask for china

#%% Fix masks for overlap and maximum above 1
masks_dict = {} # dictionary for masks

for i in range(len(mask_all)):
    temp = mask_all[i]
    temp = xr.where(temp>=1, 1, 0) # max at 1
    if i < 6:
         temp = xr.where(mask_all[6]==0, temp, 0) # ensure that no China overlap
    if i < 5:
         temp = xr.where(mask_all[5]==0, temp, 0) # ensure that no overlap in coverage
    if i < 4:
         temp = xr.where(mask_all[4]==0, temp, 0) # ensure that no overlap in coverage
    if i < 3:
         temp = xr.where(mask_all[3]==0, temp, 0) # ensure that no overlap in coverage
    if i < 2:
         temp = xr.where(mask_all[2]==0, temp, 0) # ensure that no overlap in coverage
    if i < 1:
         temp = xr.where(mask_all[1]==0, temp, 0) # ensure that no overlap in coverage
   
    masks_dict[realms_masks[i]] = temp.values
#%% Create Amazon + Neotropic_Non-Amazon rainforest masks
# Amazon square mask used by GEOS-Chem
lat_A_min = -34.125 # southernmost latitude
lat_A_max = 14.125 # northernmost latitude 
lon_A_min = -82.125 # westernmost latitude
lon_A_max = -32.875 # easternmost latitude 

lat_ind_min = np.asarray(np.where(lat==lat_A_min)).flatten()[0]
lat_ind_max = np.asarray(np.where(lat==lat_A_max)).flatten()[0]
lon_ind_min = np.asarray(np.where(lon==lon_A_min)).flatten()[0]
lon_ind_max = np.asarray(np.where(lon==lon_A_max)).flatten()[0]

# find rainforest areas
fn_dd = 'misc_Data/Olson_2001_lt_gridcells.nc'
ds_dd = xr.open_dataset(fn_dd)

# load variables
dd_lt = ds_dd.dd_lt
dd_lt_v = dd_lt.values

is_rainforest = dd_lt==6 # indices of rainforest

mask_other_neo = masks_dict['Neotrop']
mask_Amazon = np.zeros((len(lat), len(lon)))

for ii in range(lat_ind_min, lat_ind_max+1):
    for jj in range(lon_ind_min, lon_ind_max+1):
        if is_rainforest[ii,jj] == True:
            mask_other_neo[ii,jj] = 0 # remove rainforest points
            mask_Amazon[ii,jj] = 1 # add rainforest points


masks_dict['Neotrop'] = mask_other_neo # save as mask
masks_dict['Amazon'] = mask_Amazon # save as mask

#%% load deforestation areas, LUH2 dataset
# Loss of forest
fn_lt = 'misc_Data/transitions_1950_2015_total_for_loss.nc'
ds_lt = xr.open_dataset(fn_lt, decode_times=False)
total_area = ds_lt.total_area.reindex(lat=list(reversed(ds_lt.lat))).values

# total_area = ds_lt.primf_to_secdn+ds_lt.primf_to_urban+ds_lt.primf_to_c3ann + \
#     ds_lt.primf_to_c4ann+ds_lt.primf_to_c3per+ds_lt.primf_to_c4per+ \
#     ds_lt.primf_to_c3nfx+ds_lt.primf_to_pastr+ds_lt.primf_to_range+ \
#     ds_lt.secdf_to_secdn+ds_lt.secdf_to_urban+ds_lt.secdf_to_c3ann+ \
#     ds_lt.secdf_to_c4ann+ds_lt.secdf_to_c3per+ds_lt.secdf_to_c4per+ \
#     ds_lt.secdf_to_c3nfx+ds_lt.secdf_to_pastr+ds_lt.secdf_to_range
#%% Create map of deforestation emission factor, based on mask
LUC_EF_defo_map = np.zeros((len(lat), len(lon))) # best estimate
LUC_EF_defo_map_low = np.zeros((len(lat), len(lon))) # low estimate
LUC_EF_defo_map_high = np.zeros((len(lat), len(lon))) # high estimate

# loop through realms, replace emission factor in corresponding grid cell
for i in range(len(realm_names)):
    LUC_EF_defo_map = np.where(masks_dict[realm_names[i]], LUC_EF_defo[i], LUC_EF_defo_map)
    LUC_EF_defo_map_low = np.where(masks_dict[realm_names[i]], LUC_EF_defo[i]-LUC_EF_defo_err[0,i], LUC_EF_defo_map_low)
    LUC_EF_defo_map_high = np.where(masks_dict[realm_names[i]], LUC_EF_defo[i]+LUC_EF_defo_err[1,i], LUC_EF_defo_map_high)
  
#%% calculate deforestation Hg flux map, for all years
years_all =np.arange(1950,2015)
defo_Hg_flux_map = np.zeros((len(years_all)))
# multiply for all years to calculate deforestation area
total_area_def = total_area * np.broadcast_to(gbox_area_v, total_area.shape) / 1e6 # km2

defo_Hg_flux_map = total_area_def * np.broadcast_to(LUC_EF_defo_map, total_area.shape) # Mg/yr
defo_Hg_flux_map_low = total_area_def * np.broadcast_to(LUC_EF_defo_map_low, total_area.shape) # Mg/yr
defo_Hg_flux_map_high = total_area_def * np.broadcast_to(LUC_EF_defo_map_high, total_area.shape) # Mg/yr

# set nan values to zero
defo_Hg_flux_map = np.nan_to_num(defo_Hg_flux_map)
defo_Hg_flux_map_low = np.nan_to_num(defo_Hg_flux_map_low)
defo_Hg_flux_map_high = np.nan_to_num(defo_Hg_flux_map_high)

#%% accumulate over different time periods
defo_Hg_flux_map_15 = np.sum(defo_Hg_flux_map[50:,:,:], axis=0) / gbox_area_v # convert to Mg/m2
defo_Hg_flux_map_15_low = np.sum(defo_Hg_flux_map_low[50:,:,:], axis=0) / gbox_area_v # convert to Mg/m2
defo_Hg_flux_map_15_high = np.sum(defo_Hg_flux_map_high[50:,:,:], axis=0) / gbox_area_v # convert to Mg/m2

defo_Hg_flux_map_30 = np.sum(defo_Hg_flux_map[35:,:,:], axis=0) / gbox_area_v # convert to Mg/m2
defo_Hg_flux_map_30_low = np.sum(defo_Hg_flux_map_low[35:,:,:], axis=0) / gbox_area_v # convert to Mg/m2
defo_Hg_flux_map_30_high = np.sum(defo_Hg_flux_map_high[35:,:,:], axis=0) / gbox_area_v # convert to Mg/m2

defo_Hg_flux_map_45 = np.sum(defo_Hg_flux_map[20:,:,:], axis=0) / gbox_area_v # convert to Mg/m2
defo_Hg_flux_map_45_low = np.sum(defo_Hg_flux_map_low[20:,:,:], axis=0) / gbox_area_v # convert to Mg/m2
defo_Hg_flux_map_45_high = np.sum(defo_Hg_flux_map_high[20:,:,:], axis=0) / gbox_area_v # convert to Mg/m2

defo_Hg_flux_map_60 = np.sum(defo_Hg_flux_map[5:,:,:], axis=0) / gbox_area_v # convert to Mg/m2
defo_Hg_flux_map_60_low = np.sum(defo_Hg_flux_map_low[5:,:,:], axis=0) / gbox_area_v # convert to Mg/m2
defo_Hg_flux_map_60_high = np.sum(defo_Hg_flux_map_high[5:,:,:], axis=0) / gbox_area_v # convert to Mg/m2

ds_defo_Hg_flux_map_15 = xr.DataArray(data=defo_Hg_flux_map_15, 
                                      coords={'lat':lat, 'lon':lon},
                                      dims=["lat","lon"])
ds_defo_Hg_flux_map_30 = xr.DataArray(data=defo_Hg_flux_map_30, 
                                      coords={'lat':lat, 'lon':lon},
                                      dims=["lat","lon"])
ds_defo_Hg_flux_map_45 = xr.DataArray(data=defo_Hg_flux_map_45, 
                                      coords={'lat':lat, 'lon':lon},
                                      dims=["lat","lon"])
ds_defo_Hg_flux_map_60 = xr.DataArray(data=defo_Hg_flux_map_60, 
                                      coords={'lat':lat, 'lon':lon},
                                      dims=["lat","lon"])

ds_defo_Hg_flux_map_15_low = xr.DataArray(data=defo_Hg_flux_map_15_low, 
                                      coords={'lat':lat, 'lon':lon},
                                      dims=["lat","lon"])
ds_defo_Hg_flux_map_30_low = xr.DataArray(data=defo_Hg_flux_map_30_low, 
                                      coords={'lat':lat, 'lon':lon},
                                      dims=["lat","lon"])
ds_defo_Hg_flux_map_45_low = xr.DataArray(data=defo_Hg_flux_map_45_low, 
                                      coords={'lat':lat, 'lon':lon},
                                      dims=["lat","lon"])
ds_defo_Hg_flux_map_60_low = xr.DataArray(data=defo_Hg_flux_map_60_low, 
                                      coords={'lat':lat, 'lon':lon},
                                      dims=["lat","lon"])

ds_defo_Hg_flux_map_15_high = xr.DataArray(data=defo_Hg_flux_map_15_high, 
                                      coords={'lat':lat, 'lon':lon},
                                      dims=["lat","lon"])
ds_defo_Hg_flux_map_30_high = xr.DataArray(data=defo_Hg_flux_map_30_high, 
                                      coords={'lat':lat, 'lon':lon},
                                      dims=["lat","lon"])
ds_defo_Hg_flux_map_45_high = xr.DataArray(data=defo_Hg_flux_map_45_high, 
                                      coords={'lat':lat, 'lon':lon},
                                      dims=["lat","lon"])
ds_defo_Hg_flux_map_60_high = xr.DataArray(data=defo_Hg_flux_map_60_high, 
                                      coords={'lat':lat, 'lon':lon},
                                      dims=["lat","lon"])

#%% Plot map of deforestation Hg
# extend longitude by 2
lon_extend = np.zeros(lon.size+2)
# fill in internal values
lon_extend[1:-1] = lon # fill up with original values
# fill in extra endpoints
lon_extend[0] = lon[0]-np.diff(lon)[0]
lon_extend[-1] = lon[-1]+np.diff(lon)[-1]
# calculate the midpoints
lon_pcolormesh_midpoints = lon_extend[:-1]+0.5*(np.diff(lon_extend))

# extend latitude by 2
lat_extend = np.zeros(lat.size+2)
# fill in internal values
lat_extend[1:-1] = lat
# fill in extra endpoints
lat_extend[0] = lat[0]-np.diff(lat)[0]
lat_extend[-1] = lat[-1]+np.diff(lat)[-1]
# calculate the midpoints
lat_pcolormesh_midpoints = lat_extend[:-1]+0.5*(np.diff(lat_extend))

# Plot map plots
lat_min = -60
lat_max = 85
# pick the desired colormap, sensible levels, and define a normalization
# instance which takes data values and translates those into levels.
cmap = plt.cm.YlOrRd
f, ax1 = plt.subplots(1, 1, figsize=[16,6],subplot_kw=dict(projection=ccrs.PlateCarree()),
                              gridspec_kw=dict(hspace=0., wspace=0.1))
# grid = plt.GridSpec(1, 4, hspace=1, wspace=1)
# ax1 = plt.subplot(grid[:-1], projection=ccrs.PlateCarree())

h = ax1.pcolormesh(lon_pcolormesh_midpoints, lat_pcolormesh_midpoints, #vmin=0, vmax=70
                   defo_Hg_flux_map_45*1e12, cmap =cmap, norm=mpl.colors.LogNorm(vmin=1e-1),
                   rasterized = True)
#cbar.tick_params(labelsize=14)
ax1.set_ylim([lat_min, lat_max])
# ax1.set_title('Impact of reforestation on surface-atmosphere exchange', fontweight='bold', fontsize=17)

ax1.coastlines(linewidth=0.5, color='k')
#ax1.add_feature(cf.BORDERS, linewidth=0.5, edgecolor='w')
#f.subplots_adjust(bottom=0.4)

cbar = plt.colorbar(h, extend='both',orientation='horizontal',
                    ax=ax1, aspect=40, pad=0.03)
cbar.set_label('Deforestation emissions (\u03bcg m$^{-2}$ yr$^{-1}$)', fontsize = 22)
cbar.ax.tick_params(labelsize=18) 
f.savefig('Figures/Global_deforest_45_noborder_log.pdf',bbox_inches = 'tight')

#%% better to regrid this file to 0.1 x 0.1 resolution for the country calculations
# prepare regridding using xESMF
fn_countries = 'misc_Data/countrymask_0.1x0.1.nc'
ds_countries = xr.open_dataset(fn_countries)
CountryID = ds_countries.CountryID.squeeze().values # load country locations

regridder = xe.Regridder(ds_dd, ds_countries, "conservative")

# regrid datasets to 0.1 x 0.1 for country mask
defo_Hg_flux_map_15_01x01  = regridder(ds_defo_Hg_flux_map_15)
defo_Hg_flux_map_30_01x01  = regridder(ds_defo_Hg_flux_map_30)
defo_Hg_flux_map_45_01x01  = regridder(ds_defo_Hg_flux_map_45)
defo_Hg_flux_map_60_01x01  = regridder(ds_defo_Hg_flux_map_60)

defo_Hg_flux_map_15_01x01_low  = regridder(ds_defo_Hg_flux_map_15_low)
defo_Hg_flux_map_30_01x01_low  = regridder(ds_defo_Hg_flux_map_30_low)
defo_Hg_flux_map_45_01x01_low  = regridder(ds_defo_Hg_flux_map_45_low)
defo_Hg_flux_map_60_01x01_low  = regridder(ds_defo_Hg_flux_map_60_low)

defo_Hg_flux_map_15_01x01_high  = regridder(ds_defo_Hg_flux_map_15_high)
defo_Hg_flux_map_30_01x01_high  = regridder(ds_defo_Hg_flux_map_30_high)
defo_Hg_flux_map_45_01x01_high  = regridder(ds_defo_Hg_flux_map_45_high)
defo_Hg_flux_map_60_01x01_high  = regridder(ds_defo_Hg_flux_map_60_high)

# load gbox areas for 0.1 x 0.1 grid
fn_gb = 'misc_Data/gbox_01x01.nc'
ds_gb = xr.open_dataset(fn_gb)
gbox_area_01x01 = ds_gb.cell_area # grid_area
lon_01x01 = ds_gb.lon
lat_01x01 = ds_gb.lat


# multiply by grid area to get units of Mg
defo_Hg_flux_map_15_01x01_u = defo_Hg_flux_map_15_01x01 * gbox_area_01x01
defo_Hg_flux_map_30_01x01_u = defo_Hg_flux_map_30_01x01 * gbox_area_01x01
defo_Hg_flux_map_45_01x01_u = defo_Hg_flux_map_45_01x01 * gbox_area_01x01
defo_Hg_flux_map_60_01x01_u = defo_Hg_flux_map_60_01x01 * gbox_area_01x01

defo_Hg_flux_map_15_01x01_low_u = defo_Hg_flux_map_15_01x01_low * gbox_area_01x01
defo_Hg_flux_map_30_01x01_low_u = defo_Hg_flux_map_30_01x01_low * gbox_area_01x01
defo_Hg_flux_map_45_01x01_low_u = defo_Hg_flux_map_45_01x01_low * gbox_area_01x01
defo_Hg_flux_map_60_01x01_low_u = defo_Hg_flux_map_60_01x01_low * gbox_area_01x01

defo_Hg_flux_map_15_01x01_high_u = defo_Hg_flux_map_15_01x01_high * gbox_area_01x01
defo_Hg_flux_map_30_01x01_high_u = defo_Hg_flux_map_30_01x01_high * gbox_area_01x01
defo_Hg_flux_map_45_01x01_high_u = defo_Hg_flux_map_45_01x01_high * gbox_area_01x01
defo_Hg_flux_map_60_01x01_high_u = defo_Hg_flux_map_60_01x01_high * gbox_area_01x01

#%% get shape data from regionmask package
countries = regionmask.defined_regions.natural_earth_v5_0_0.countries_110
# create country mask
country_mask = countries.mask_3D(ds_countries)
regions_names = country_mask.region.names.values

#%% Find areas where have issues with mask
country_mask_all = country_mask.sum("region") 

#%% Make edits to mask - to fix some coastline island issues
country_mask = mask_square(country_mask, 'Philippines', 7.65, 20.15, 116.05, 127.05)
country_mask = mask_square(country_mask, 'Philippines', 4.25, 9.75, 121.05, 129.05)
country_mask = mask_square(country_mask, 'Malaysia', 4.55, 7.55, 115.55, 119.45)
country_mask = mask_square(country_mask, 'Malaysia', 2.95, 5.55, 100.45, 102.35)
country_mask = mask_square(country_mask, 'Malaysia', 1.25, 3.25, 110.15, 112.05)
country_mask = mask_square(country_mask, 'Malaysia', 5.05, 6.45, 99.55, 100.75)
country_mask = mask_square(country_mask, 'Malaysia', 1.95, 2.75, 101.95, 103.05)
country_mask = mask_square(country_mask, 'Indonesia', -6.05, 2.25, 115.15, 131.55)
country_mask = mask_square(country_mask, 'Indonesia', -8.05, 1.15, 98.45, 109.25)
country_mask = mask_square(country_mask, 'Indonesia', 1.05, 6.05, 94.35, 99.55)
country_mask = mask_square(country_mask, 'Indonesia', 1.05, 6.05, 94.35, 99.55)
country_mask = mask_square(country_mask, 'Indonesia', -10.85, 0.05, 106.25, 122.75)
country_mask = mask_square(country_mask, 'Thailand',6.75, 9.85, 97.55, 101.05)
country_mask = mask_square(country_mask, 'Thailand',8.95, 10.75, 98.95, 100.45)
country_mask = mask_square(country_mask, 'Vietnam',20.55, 21.55, 106.45, 108.05)
country_mask = mask_square(country_mask, 'Vietnam',9.85, 10.55, 103.85, 104.25)
country_mask = mask_square(country_mask, 'Vietnam',9.35, 10.45, 104.35, 105.55)
country_mask = mask_square(country_mask, 'Vietnam',19.25,21.55,105.25, 108.05)
country_mask = mask_square(country_mask, 'Vietnam',9.85,12.25,107.15, 110.05)
country_mask = mask_square(country_mask, 'Myanmar',18.05, 21.15, 92.45, 94.55)
country_mask = mask_square(country_mask, 'Myanmar',14.45,17.25, 95.95, 98.15)
country_mask = mask_square(country_mask, 'Bangladesh',21.05, 23.55, 89.15, 91.25)
country_mask = mask_square(country_mask, 'Bangladesh',21.05, 22.95, 90.95, 92.35)
country_mask = mask_square(country_mask, 'India',21.05, 22.35, 87.05, 89.05)
country_mask = mask_square(country_mask, 'Sri Lanka',5.55, 9.05, 79.45, 82.75)
country_mask = mask_square(country_mask, 'Madagascar',-26.15, -14.05, 42.55, 51.05)
country_mask = mask_square(country_mask, 'El Salvador',12.75, 13.85, -89.85, -87.95)
country_mask = mask_square(country_mask, 'Taiwan', 21.75,25.65, 119.75, 122.35)
country_mask = mask_square(country_mask, 'India', 7.55, 22.95,70.85, 78.05)
country_mask = mask_square(country_mask, 'China', 27.85, 38.45,116.05, 123.65)
country_mask = mask_square(country_mask, 'Brazil', -31.95, 1.75,-52.85, -33.75)
country_mask = mask_square(country_mask, 'Haiti', 17.75, 20.25, -73.95, -72.15)
country_mask = mask_square(country_mask, 'Honduras', 13.15, 13.85, -87.65, -87.15)
country_mask = mask_square(country_mask, 'Panama', 6.85, 10.55, -82.35, -78.25)
country_mask = mask_square(country_mask, 'Guyana', 5.75, 8.35, -59.45, -57.35)

#%% check for errors 
country_mask_all = country_mask.sum("region")
print("If no errors in mask, returns [0 1]: ")
print(np.unique(country_mask_all)) # make sure don't have two countries assigned to same grid box

#%%
# calculate total deforestation Hg flux by country - first for total of 45 year time horizon
country_Hg_flux_45 = defo_Hg_flux_map_45_01x01_u.weighted(country_mask).sum(dim=("lat", "lon"))
print(sum(country_Hg_flux_45).values)
print(defo_Hg_flux_map_45_01x01_u.sum() - sum(country_Hg_flux_45).values)

# calculate for other fields
country_Hg_flux_15 = defo_Hg_flux_map_15_01x01_u.weighted(country_mask).sum(dim=("lat", "lon")).values
country_Hg_flux_15_low = defo_Hg_flux_map_15_01x01_low_u.weighted(country_mask).sum(dim=("lat", "lon")).values
country_Hg_flux_15_high = defo_Hg_flux_map_15_01x01_high_u.weighted(country_mask).sum(dim=("lat", "lon")).values
country_Hg_flux_30 = defo_Hg_flux_map_30_01x01_u.weighted(country_mask).sum(dim=("lat", "lon")).values
country_Hg_flux_30_low = defo_Hg_flux_map_30_01x01_low_u.weighted(country_mask).sum(dim=("lat", "lon")).values
country_Hg_flux_30_high = defo_Hg_flux_map_30_01x01_high_u.weighted(country_mask).sum(dim=("lat", "lon")).values
country_Hg_flux_45_low = defo_Hg_flux_map_45_01x01_low_u.weighted(country_mask).sum(dim=("lat", "lon")).values
country_Hg_flux_45_high = defo_Hg_flux_map_45_01x01_high_u.weighted(country_mask).sum(dim=("lat", "lon")).values
country_Hg_flux_60 = defo_Hg_flux_map_60_01x01_u.weighted(country_mask).sum(dim=("lat", "lon")).values
country_Hg_flux_60_low = defo_Hg_flux_map_60_01x01_low_u.weighted(country_mask).sum(dim=("lat", "lon")).values
country_Hg_flux_60_high = defo_Hg_flux_map_60_01x01_high_u.weighted(country_mask).sum(dim=("lat", "lon")).values


print("Total emiss 45 low")
print(np.sum(defo_Hg_flux_map_45_01x01_low_u))
print("Total emiss 45 best")
print(np.sum(defo_Hg_flux_map_45_01x01_u))
print("Total emiss 45 high")
print(np.sum(defo_Hg_flux_map_45_01x01_high_u))

# save data files for future comparisons
np.savetxt("misc_Data/country_DFR_15.csv", country_Hg_flux_15, delimiter=",")
np.savetxt("misc_Data/country_DFR_30.csv", country_Hg_flux_30, delimiter=",")
np.savetxt("misc_Data/country_DFR_45.csv", country_Hg_flux_45, delimiter=",")
np.savetxt("misc_Data/country_DFR_60.csv", country_Hg_flux_60, delimiter=",")

np.savetxt("misc_Data/country_DFR_15_low.csv", country_Hg_flux_15_low, delimiter=",")
np.savetxt("misc_Data/country_DFR_30_low.csv", country_Hg_flux_30_low, delimiter=",")
np.savetxt("misc_Data/country_DFR_45_low.csv", country_Hg_flux_45_low, delimiter=",")
np.savetxt("misc_Data/country_DFR_60_low.csv", country_Hg_flux_60_low, delimiter=",")

np.savetxt("misc_Data/country_DFR_15_high.csv", country_Hg_flux_15_high, delimiter=",")
np.savetxt("misc_Data/country_DFR_30_high.csv", country_Hg_flux_30_high, delimiter=",")
np.savetxt("misc_Data/country_DFR_45_high.csv", country_Hg_flux_45_high, delimiter=",")
np.savetxt("misc_Data/country_DFR_60_high.csv", country_Hg_flux_60_high, delimiter=",")

#%% Plot where have missing data from countries due to country mask
lon_extend = np.zeros(lon_01x01.size+2)
# fill in internal values
lon_extend[1:-1] = lon_01x01 # fill up with original values
# fill in extra endpoints
lon_extend[0] = lon_01x01[0]-np.diff(lon_01x01)[0]
lon_extend[-1] = lon_01x01[-1]+np.diff(lon_01x01)[-1]
# calculate the midpoints
lon_pcolormesh_midpoints = lon_extend[:-1]+0.5*(np.diff(lon_extend))

# extend latitude by 2
lat_extend = np.zeros(lat_01x01.size+2)
# fill in internal values
lat_extend[1:-1] = lat_01x01
# fill in extra endpoints
lat_extend[0] = lat_01x01[0]-np.diff(lat_01x01)[0]
lat_extend[-1] = lat_01x01[-1]+np.diff(lat_01x01)[-1]
# calculate the midpoints
lat_pcolormesh_midpoints = lat_extend[:-1]+0.5*(np.diff(lat_extend))

cmap = plt.cm.magma
f, ax1 = plt.subplots(1, 1, figsize=[16,6],subplot_kw=dict(projection=ccrs.PlateCarree()),
                              gridspec_kw=dict(hspace=0., wspace=0.1))
missing_areas = np.where(country_mask_all==True, 0, defo_Hg_flux_map_45_01x01_u)
# Check missing grid boxes
missing_areas_vec = missing_areas.flatten()
missing_areas_max = missing_areas_vec[np.argsort(missing_areas_vec)[-20:]]
# for i in range(20):
#     inds = np.where(missing_areas == missing_areas_max[i])
#     lat_i = lat[inds[0][0]].values
#     lon_i = lon[inds[1][0]].values
#     print(i)
#     print(lat_i)
#     print(lon_i)


h = ax1.pcolormesh(lon_pcolormesh_midpoints, lat_pcolormesh_midpoints,
                   missing_areas*1e3, cmap =cmap,
                   rasterized = True, vmin=0, vmax=2)
#cbar.tick_params(labelsize=14)
lat_min=-60
lat_max=80
ax1.set_ylim([lat_min, lat_max])
# ax1.set_title('Impact of reforestation on surface-atmosphere exchange', fontweight='bold', fontsize=17)

ax1.coastlines(linewidth=0.01, color='w')
#ax1.add_feature(cf.BORDERS, linewidth=0.5, edgecolor='w')
#f.subplots_adjust(bottom=0.4)

cbar = plt.colorbar(h, extend='max',orientation='horizontal',
                    ax=ax1, aspect=40, pad=0.03)
cbar.set_label('Deforestation emissions (kg yr$^{-1}$)', fontsize = 22)
cbar.ax.tick_params(labelsize=18) 
f.savefig('Figures/Missing_country_mask_DFR_nocoast2.pdf',bbox_inches = 'tight')
#%% Assign country level emissions to geopandas dataframe
# deforestation emissions
df_countries = regionmask.defined_regions.natural_earth_v5_0_0.countries_110.to_geodataframe()
df_countries['defo'] = country_Hg_flux_45

# primary emissions
df_primary = pd.read_csv('misc_Data/emissions_countries_GMA_regionnmask.csv')
primary_countries_names = df_primary.iloc[:,1]

primary_emiss = np.zeros(len(df_countries)) # array to say country emiss vals
# compare country names in both arrays to find matching indices
for i in range(len(df_countries)):
    region_i = regions_names[i]
    if region_i == "Côte d'Ivoire": # fix accent issue
        ind_i = list(primary_countries_names).index('Côte dIvoire')
    elif region_i == "N. Cyprus": # not in AMAP, set name to Cyprus
        ind_i = list(primary_countries_names).index('Cyprus')
    elif region_i == "Somaliland": # not in AMAP, set name to Somalia
        ind_i = list(primary_countries_names).index('Somalia')        
    else: 
        ind_i = list(primary_countries_names).index(region_i)        
    primary_emiss[i] = df_primary.iloc[ind_i,3] # emissions in Mg/yr
    
    if region_i == 'Antarctica' or region_i == 'Fr. S. Antarctic Lands': # set emissions to zero for Antarctica
        primary_emiss[i] = 0

df_countries['primary'] = primary_emiss
df_countries['ratio'] = df_countries['defo']/(df_countries['primary'])
ratio_countries = df_countries['defo']/(df_countries['primary'])
np.savetxt("misc_Data/country_ratios_45.csv", ratio_countries.values, delimiter=",")
np.savetxt("misc_Data/country_AMAP_emiss.csv", primary_emiss, delimiter=",")

#%% Plot country emissions as a map
fig, ax = plt.subplots(1, 2, figsize=[16,6],
                       subplot_kw=dict(projection=ccrs.PlateCarree()),
                       gridspec_kw=dict(hspace=0., wspace=0.1))
ax = ax.flatten()
# deforestation emiss plot
df_countries.plot(column='defo', cmap = "YlOrRd", #legend=True, #cax=cax,
                  ax=ax[0], scheme="User_Defined", #vmin=0, vmax=20)
                  classification_kwds=dict(bins=[0.1, 0.5, 1, 2, 3, 5, 10, 15,30, 40]))
ax[0].coastlines(linewidth=0.5)
ax[0].add_feature(cf.BORDERS, linewidth=0.5)
ax[0].set_ylim([-60,85])
ax[0].set_xlim([-180,180])
ax[0].set_xticks([]) 
ax[0].set_yticks([]) 

# ratio with primary emissions plot
df_countries.plot(column='ratio', cmap = "Greens", #legend=True, #cax=cax,
                  ax=ax[1], scheme="User_Defined", #vmin=0, vmax=20)
                  classification_kwds=dict(bins=[0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.]))
ax[1].coastlines(linewidth=0.5)
ax[1].add_feature(cf.BORDERS, linewidth=0.5)
ax[1].set_ylim([-60,85])
ax[1].set_xlim([-180,180])
ax[1].set_xticks([]) 
ax[1].set_yticks([])

fig.savefig('Figures/Map_country_DFR_ratio.pdf',bbox_inches = 'tight')


