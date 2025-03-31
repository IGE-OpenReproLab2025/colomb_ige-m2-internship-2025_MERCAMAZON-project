#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 2023
Running offline version of dry deposition of GEOS-Chem on monthly timescale, hourly average diurnal cycle
Only for DFR grid boxes in HIST scenario, to calculate reference totals
GEOS-Chem input data can be downloaded at: http://geoschemdata.wustl.edu/ExtData/

@author: arifeinberg
"""
#import os
#os.chdir('/Users/arifeinberg/target2/fs03/d0/arifein/python/')
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import scipy.io as sio
from drydep_functions import Compute_Olson_landmap, METERO, DEPVEL_UNC

# realm parameters
realm = 'Amazon'
print(realm)
# overall options for realm names
realms_masks = ['Nearc', 'Palearc','Afrotrop','Indomalay','Austral','Neotrop', 'China', 'Amazon'] 
realm_id = realms_masks.index(realm)

#%% Load necessary variables for running dry deposition code

# required Olson land variables
fn_ols1 = 'misc_Data/Olson_2001_Drydep_Inputs.nc'
ds_ols1 = xr.open_dataset(fn_ols1)

DRYCOEFF = ds_ols1.DRYCOEFF.values # Baldocchi dry deposition polynomial coefficients
IOLSON = ds_ols1.IOLSON.values # Olson land type indices (+1)
IDRYDEP = ds_ols1.IDRYDEP.values # Dry deposition land types

IDEP = ds_ols1.IDEP.values # Mapping index: Olson land type ID to drydep ID
IZO = ds_ols1.IZO.values # Default roughness heights for each Olson land type
IRI = ds_ols1.IRI.values # RI resistance for each dry deposition land type
IRI[2] = 200. # Change RI resistance of coniferous forests to match deciduous, as in GEOS-Chem
IRLU = ds_ols1.IRLU.values # RLU resistance for each dry deposition land type
IRAC = ds_ols1.IRAC.values # RAC resistance for each dry deposition land type
IRGSS = ds_ols1.IRGSS.values # RGSS resistance for each dry deposition land type
IRGSO = ds_ols1.IRGSO.values # RGSO resistance for each dry deposition land type
IRCLS = ds_ols1.IRCLS.values # RCLS resistance for each dry deposition land type
IRCLO = ds_ols1.IRCLO.values # RCLO resistance for each dry deposition land type
IVSMAX = ds_ols1.IVSMAX.values # Max drydep velocity (for aerosol) for each dry deposition land type

# Land type areas - REF
fn_ols2 = 'misc_Data/Olson_2001_Land_Type_Masks.2_25.generic.nc'
ds_ols2 = xr.open_dataset(fn_ols2)
# save as one np array, first dim is land types
Olson_landtype = ds_ols2.to_array().squeeze()
lat = ds_ols2.lat.values # latitude
lon = ds_ols2.lon.values # longitude
lat_l = len(lat)
lon_l = len(lon)

# Land type LAI
fn_xlai  = 'misc_Data/Yuan_proc_MODIS_XLAI.2x25.BAU_2002_m.nc'
ds_xlai = xr.open_dataset(fn_xlai)

# save as one data array, first dim is land types
XLAI_m = ds_xlai.drop('time_bnds').to_array().values
#%% calculate mask for within Amazon area, to know which f0 reactivity to apply
in_AMZN = np.full((lat_l, lon_l), False)
# Amazon square mask used by GEOS-Chem
lat_A_min = -32. # southernmost latitude
lat_A_max = 12. # northernmost latitude 
lon_A_min = -80. # westernmost latitude
lon_A_max = -35. # easternmost latitude 

lat_ind_min = np.asarray(np.where(lat==lat_A_min)).flatten()[0]
lat_ind_max = np.asarray(np.where(lat==lat_A_max)).flatten()[0]
lon_ind_min = np.asarray(np.where(lon==lon_A_min)).flatten()[0]
lon_ind_max = np.asarray(np.where(lon==lon_A_max)).flatten()[0]

in_AMZN[lat_ind_min:lat_ind_max+1, lon_ind_min:lon_ind_max+1] = True
#%% 
# load HIST Hg0 species concentration, for use in dry deposition flux calculations
fn_Hg0 = '../GEOS-Chem_runs/run0311/OutputDir/GEOSChem.SpeciesConc.alltime_m.nc4'
ds_Hg0 = xr.open_dataset(fn_Hg0)
Hg0_conc = ds_Hg0.SpeciesConc_Hg0.isel(lev=0, time=slice(12, 24)).values # Hg0 in mole/mole
#%% load area and masks
# Load grid cell area for unit conversion of model
fn_gbox = '../gbox_areas/GEOSChem_2x25_gboxarea.nc'
ds_gbox = xr.open_dataset(fn_gbox)
gbox_GC = ds_gbox.cell_area.values #m2
Mg_ug = 1e12 # Mg in ug 
#%% Set land type that replaces forest land type
if realm in ['Afrotrop', 'Indomalay','Palearc','Austral','China']: 
    replace_id = 31 # crops and town index that is replaced
elif realm in ['Neotrop','Nearc']:
    replace_id = 35 # Corn and Beans Croplands index that is replaced
elif realm in ['Amazon']:
    replace_id = 58 # Fields and savannah index that is replaced
else:
    print('error with realm name, not found')

#%% Load realm mask
fn_realm_2x25 = 'misc_Data/DFR/Olson_Land_Type_Masks.2x25.DFR_'+realm+'.nc'
# at 2 x 2.5 resolution, for mask calculation    
ds_realm_2x25 = xr.open_dataset(fn_realm_2x25)
Olson_new_2x25 = ds_realm_2x25.to_array().squeeze()

mask_change_2x25 = Olson_new_2x25[replace_id, :, :] - Olson_landtype[replace_id, :, :] # changed area
# adjust mask so that accounts for all changed area, otherwise missing edge effects
mask_change_2x25 = xr.where(mask_change_2x25>0, 1, 0) 

# Get the array indices that need to calculate for realm
re_calc_bool = np.where(mask_change_2x25>0, True, False)

#%% Load land masks
fn_fr = '../misc_datasets/MERRA2_proc/MERRA2_FROCEAN_FRLANDIC_2x25.nc4'
ds_fr = xr.open_dataset(fn_fr)
FRLANDIC = ds_fr.FRLANDIC.squeeze().values # land glaciated areas
FROCEAN = ds_fr.FROCEAN.squeeze().values # ocean areas
# load snow cover, for use in dry deposition flux calculations
fn_sno = '../misc_datasets/MERRA2_proc/MERRA2_FRSNO_2015_2x25_m.nc4'
ds_sno = xr.open_dataset(fn_sno)
FRSNO = ds_sno.FRSNO.values # snow coverage on land, monthly
#%% load uncertainty parameter values
# columns: 1) soil parametrization 2) LAI_pct_replace 3)- f0 Amazon RF 
#          4) f0 other RF 5) f0 elsewhere 6) biomass burning emissions
params_unc = np.loadtxt("misc_Data/LHS_sampling_dep_emis.csv",
                 delimiter=",")
num = len(params_unc) # number of uncertainty calculations
#%% run dep velocity function over all months, diurnal cycle
import time
t0 = time.time()
# constants
MW_air = 28.97 # g/mol, molar mass of dry air
avo = 6.02e23 # avogadro number molec mol^-1
g_kg = 1e3 # g in kg
cm3_m3 = 1e6 # cm3 in m3
ug_g = 1e6 # ug in g
cm2_m2 = 1e4 # cm^2 in m^2
MW_Hg = 200.59 # g mol^-1
s_yr = 60*60*24*365.2425 # s in yr

#  parameters
XMW = 201.0 * 1e-3 # Hg0 molar mass kg/mol
F0 = 3e-5 # Hg0 reactivity
F0_RF = 3e-5 # Hg0 reactivity in other rainforests
F0_AMZN = 0.2 # F0 for Amazon rainforest type

HSTAR = 0.11 # Hg0 

# get additional information about Olson land types
FRCLND, IREG, ILAND, IUSE = Compute_Olson_landmap(Olson_landtype)
IREG_v = IREG.values
ILAND_v = ILAND.values
IUSE_v = IUSE.values

# subset overall non-monthly arrays
gbox_GC_re = gbox_GC[re_calc_bool]
n_re = len(gbox_GC_re) # number of grid cells to calculate
IREG_v_re = IREG_v[re_calc_bool]
FROCEAN_re = FROCEAN[re_calc_bool]
FRLANDIC_re = FRLANDIC[re_calc_bool]
in_AMZN_re = in_AMZN[re_calc_bool]

#3D vars - need to be 2D array in the output, so reshape
re_calc_bool3D = np.broadcast_to(re_calc_bool, Olson_landtype.shape)# need 3D mask for Olson
ILAND_v_re = np.reshape(ILAND_v[re_calc_bool3D], (73, n_re))
IUSE_v_re = np.reshape(IUSE_v[re_calc_bool3D], (73, n_re))

#time array - need to be 2D array in the output, so reshape
re_calc_boolt = np.broadcast_to(re_calc_bool, FRSNO.shape)# need 3D mask for Olson
FRSNO_re = np.reshape(FRSNO[re_calc_boolt], (12, n_re))

# initialize
DV_m = np.zeros((num, 12,n_re)) # array of deposition velocity maps
DF_m = np.zeros((num, 12,n_re)) # array of deposition flux maps

# run monthly loop
for mm in range(12): # monthly loop
    print(mm)
    # load required met variables
    
    fn_met = '../GC_data/run0311/GEOSChem.StateMet.2015' +str(mm+1).zfill(2) + '01_diurnal.nc4'
    ds_met = xr.open_dataset(fn_met)
    
    bxheight = ds_met.Met_BXHEIGHT.isel(lev=0).values # grid box height
    surf_pres = ds_met.Met_PSC2WET.values * 100.0 # surface pressure in Pa
    Z0 = ds_met.Met_Z0.values # surface roughness
    CLDFRC = ds_met.Met_CLDFRC.values # column cloud fraction
    albedo = ds_met.Met_ALBD.values # visible surface albedo
    airden = ds_met.Met_AIRDEN.isel(lev=0).values # air density at surface
    air_molec_cm3 = airden / (MW_air / g_kg) * avo / cm3_m3 # airden in molec/cm3
    hflux = ds_met.Met_HFLUX.values # sensible heat flux
    sw_grnd = ds_met.Met_SWGDN.values # incident shortwave at ground
    surf_t = ds_met.Met_TS.values # surface temperature 
    ustar = ds_met.Met_USTAR.values # friction velocity
    U10M = ds_met.Met_U10M.values # zonal wind at 10 m height
    V10M = ds_met.Met_V10M.values # meridional wind at 10 m height
    SUNCOSmid = ds_met.Met_SUNCOSmid.values # cosine of solar zenith angle, at midpoint of chemistry timestep
    
    Hg0_conc_mm = Hg0_conc[mm, :,:] * np.mean(air_molec_cm3, axis=0) # Hg0 conc in molec/cm3

    # Select only real values
    re_calc_bool_24H = np.broadcast_to(re_calc_bool, bxheight.shape)# need 3D mask for 24H arrays

    bxheight_re = np.reshape(bxheight[re_calc_bool_24H], (24, n_re))
    surf_pres_re = np.reshape(surf_pres[re_calc_bool_24H], (24, n_re))
    Z0_re = np.reshape(Z0[re_calc_bool_24H], (24, n_re)) # surface roughness
    CLDFRC_re = np.reshape(CLDFRC[re_calc_bool_24H], (24, n_re)) # column cloud fraction
    albedo_re = np.reshape(albedo[re_calc_bool_24H], (24, n_re)) # visible surface albedo
    airden_re = np.reshape(airden[re_calc_bool_24H], (24, n_re)) # air density at surface
    hflux_re = np.reshape(hflux[re_calc_bool_24H], (24, n_re)) # sensible heat flux
    sw_grnd_re = np.reshape(sw_grnd[re_calc_bool_24H], (24, n_re)) # incident shortwave at ground
    surf_t_re = np.reshape(surf_t[re_calc_bool_24H], (24, n_re)) # surface temperature 
    ustar_re = np.reshape(ustar[re_calc_bool_24H], (24, n_re)) # friction velocity
    U10M_re = np.reshape(U10M[re_calc_bool_24H], (24, n_re)) # zonal wind at 10 m height
    V10M_re = np.reshape(V10M[re_calc_bool_24H], (24, n_re)) # meridional wind at 10 m height
    SUNCOSmid_re = np.reshape(SUNCOSmid[re_calc_bool_24H], (24, n_re)) # cosine of solar zenith angle, at midpoint of chemistry timestep
    
    Hg0_conc_mm_re = Hg0_conc_mm[re_calc_bool] # Hg0 conc in molec/cm3

    # XLAI, 4 D variable
    XLAI_mm = XLAI_m[:,mm,:,:]
    XLAI_mm_re = np.reshape(XLAI_mm[re_calc_bool3D], (73, n_re))
    # run loop over parameter values
    for pp in range(num): # loop over parameters
        # set F0 parameters for this run
        F0 = params_unc[pp,4] # f0 elsewhere
        F0_AMZN = params_unc[pp,2] # f0 in Amazon rainforest
        F0_RF = params_unc[pp,3] # f0 in other rainforests
        
        DV_h = np.zeros((24,n_re)) # array of hourly deposition velocity maps

        for i in range(24): # loop over hours
            for jj in range(n_re):
                if (FROCEAN_re[jj] + FRLANDIC_re[jj] <1): # don't calculate over water or ice
                    CZ1, LSNOW, OBK, W10 = METERO(bxheight_re[i,jj],\
                                                       albedo_re[i,jj], \
                                                       surf_t_re[i,jj], \
                                                       ustar_re[i,jj], \
                                                       airden_re[i,jj], \
                                                       hflux_re[i,jj], \
                                                       U10M_re[i,jj], \
                                                       V10M_re[i,jj]) 
                            
                    DV_h[i,jj] = DEPVEL_UNC(DRYCOEFF, IOLSON, IDEP, IRI, \
                               IRLU, IRAC, IRGSS, IRGSO, IRCLS, IRCLO, \
                               IREG_v_re[jj],\
                               ILAND_v_re[:,jj],\
                               IUSE_v_re[:,jj], \
                               surf_t_re[i,jj], \
                               XLAI_mm_re[:,jj], \
                               LSNOW, sw_grnd_re[i,jj], CLDFRC_re[i,jj], \
                               SUNCOSmid_re[i,jj], surf_pres_re[i,jj], \
                               ustar_re[i,jj], Z0_re[i,jj], \
                                   CZ1, OBK, XMW, F0, F0_AMZN, F0_RF, in_AMZN_re[jj], HSTAR)
        DV_m[pp,mm,:] = np.mean(DV_h, axis = 0) # take monthly average of diurnal cycles
   
        # missing coverage due to snow + ice + water
        Frac_no_dep = FRSNO_re[mm, :]+ FROCEAN_re + FRLANDIC_re

        DF_m[pp,mm,:] = DV_m[pp,mm, :] * Hg0_conc_mm_re * (1. - Frac_no_dep) # calculate deposition flux
    

DF_ann = np.mean(DF_m, axis=1) # take annual average
unit_conv = MW_Hg / avo * ug_g * cm2_m2 * s_yr
DF_ann = DF_ann * unit_conv # convert to ug/m2/yr

realm_dep = np.zeros(num)  # for listing Amazon totals
for pp in range(num):
    realm_dep[pp] = np.sum(DF_ann[pp,:] * gbox_GC_re ) / Mg_ug # Mg/yr
    print('Realm dry deposition: p' + str(pp))
    print(realm_dep[pp])

# save data to csv, for access
np.savetxt("misc_Data/unc_realm_dep_HIST_" +realm+".csv", realm_dep, delimiter=",")

t1 = time.time()

total_t = t1-t0
print('Execution time in seconds: ' + str(total_t))
