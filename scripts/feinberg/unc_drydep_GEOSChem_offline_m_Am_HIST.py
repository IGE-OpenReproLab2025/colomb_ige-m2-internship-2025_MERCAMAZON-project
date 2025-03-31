#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 6 2023
Running offline version of dry deposition of GEOS-Chem on monthly timescale, hourly average diurnal cycle
Only for Amazon grid boxes
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

# Land type areas
fn_ols2 = 'misc_Data/Olson_Land_Type_Mask.2_25.BAU_2002.nc'
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

# Load Amazon mask
fn_Am_mask = '../masks_2_25/Amazon_basin_mask_2x25.nc'
ds_Am_mask = xr.open_dataset(fn_Am_mask)
Am_mask = ds_Am_mask.MASK.values #unitless
Mg_ug = 1e12 # Mg in ug 

# Get the array indices that need to calculate for Amazon
Am_calc_bool = np.where(Am_mask>0, True, False)

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
HSTAR = 0.11 # Hg0 

# get additional information about Olson land types
FRCLND, IREG, ILAND, IUSE = Compute_Olson_landmap(Olson_landtype)
IREG_v = IREG.values
ILAND_v = ILAND.values
IUSE_v = IUSE.values

# subset overall non-monthly arrays
Am_mask_Am = Am_mask[Am_calc_bool]# weights for Amazon averages
n_Am = len(Am_mask_Am) # number of grid cells to calculate
gbox_GC_Am = gbox_GC[Am_calc_bool]
IREG_v_Am = IREG_v[Am_calc_bool]
FROCEAN_Am = FROCEAN[Am_calc_bool]
FRLANDIC_Am = FRLANDIC[Am_calc_bool]
in_AMZN_Am = in_AMZN[Am_calc_bool]

#3D vars - need to be 2D array in the output, so reshape
Am_calc_bool3D = np.broadcast_to(Am_calc_bool, Olson_landtype.shape)# need 3D mask for Olson
ILAND_v_Am = np.reshape(ILAND_v[Am_calc_bool3D], (73, n_Am))
IUSE_v_Am = np.reshape(IUSE_v[Am_calc_bool3D], (73, n_Am))

#time array - need to be 2D array in the output, so reshape
Am_calc_boolt = np.broadcast_to(Am_calc_bool, FRSNO.shape)# need 3D mask for Olson
FRSNO_Am = np.reshape(FRSNO[Am_calc_boolt], (12, n_Am))

# initialize
DV_m = np.zeros((num, 12,n_Am)) # array of deposition velocity maps
DF_m = np.zeros((num, 12,n_Am)) # array of deposition flux maps

# run monthly loop
for mm in range(12): # monthly loop
    print(mm) # print monthly index
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

    # Select only Amazon basin values
    Am_calc_bool_24H = np.broadcast_to(Am_calc_bool, bxheight.shape)# need 3D mask for 24H arrays

    bxheight_Am = np.reshape(bxheight[Am_calc_bool_24H], (24, n_Am))
    surf_pres_Am = np.reshape(surf_pres[Am_calc_bool_24H], (24, n_Am))
    Z0_Am = np.reshape(Z0[Am_calc_bool_24H], (24, n_Am)) # surface roughness
    CLDFRC_Am = np.reshape(CLDFRC[Am_calc_bool_24H], (24, n_Am)) # column cloud fraction
    albedo_Am = np.reshape(albedo[Am_calc_bool_24H], (24, n_Am)) # visible surface albedo
    airden_Am = np.reshape(airden[Am_calc_bool_24H], (24, n_Am)) # air density at surface
    hflux_Am = np.reshape(hflux[Am_calc_bool_24H], (24, n_Am)) # sensible heat flux
    sw_grnd_Am = np.reshape(sw_grnd[Am_calc_bool_24H], (24, n_Am)) # incident shortwave at ground
    surf_t_Am = np.reshape(surf_t[Am_calc_bool_24H], (24, n_Am)) # surface temperature 
    ustar_Am = np.reshape(ustar[Am_calc_bool_24H], (24, n_Am)) # friction velocity
    U10M_Am = np.reshape(U10M[Am_calc_bool_24H], (24, n_Am)) # zonal wind at 10 m height
    V10M_Am = np.reshape(V10M[Am_calc_bool_24H], (24, n_Am)) # meridional wind at 10 m height
    SUNCOSmid_Am = np.reshape(SUNCOSmid[Am_calc_bool_24H], (24, n_Am)) # cosine of solar zenith angle, at midpoint of chemistry timestep
    
    Hg0_conc_mm_Am = Hg0_conc_mm[Am_calc_bool] # Hg0 conc in molec/cm3

    # XLAI, 4 D variable
    XLAI_mm = XLAI_m[:,mm,:,:]
    XLAI_mm_Am = np.reshape(XLAI_mm[Am_calc_bool3D], (73, n_Am))
    
    # run loop over parameter values
    for pp in range(num): # loop over parameters
        # set F0 parameters for this run
        F0 = params_unc[pp,4] # f0 elsewhere
        F0_AMZN = params_unc[pp,2] # f0 in Amazon rainforest
        F0_RF = params_unc[pp,3] # f0 in other rainforests
        
        # initialize parameters
        DV_h = np.zeros((24,n_Am)) # array of hourly deposition velocity maps
        
        for i in range(24): # loop over hours
            for jj in range(n_Am):
                if (FROCEAN_Am[jj] + FRLANDIC_Am[jj] <1): # don't calculate over water or ice
                    CZ1, LSNOW, OBK, W10 = METERO(bxheight_Am[i,jj],\
                                                       albedo_Am[i,jj], \
                                                       surf_t_Am[i,jj], \
                                                       ustar_Am[i,jj], \
                                                       airden_Am[i,jj], \
                                                       hflux_Am[i,jj], \
                                                       U10M_Am[i,jj], \
                                                       V10M_Am[i,jj]) 
                            
                    DV_h[i,jj] = DEPVEL_UNC(DRYCOEFF, IOLSON, IDEP, IRI, \
                               IRLU, IRAC, IRGSS, IRGSO, IRCLS, IRCLO, \
                               IREG_v_Am[jj],\
                               ILAND_v_Am[:,jj],\
                               IUSE_v_Am[:,jj], \
                               surf_t_Am[i,jj], \
                               XLAI_mm_Am[:,jj], \
                               LSNOW, sw_grnd_Am[i,jj], CLDFRC_Am[i,jj], \
                               SUNCOSmid_Am[i,jj], surf_pres_Am[i,jj], \
                               ustar_Am[i,jj], Z0_Am[i,jj], \
                               CZ1, OBK, XMW, F0, F0_AMZN, F0_RF, in_AMZN_Am[jj], HSTAR)
        DV_m[pp,mm,:] = np.mean(DV_h, axis = 0) # take monthly average of diurnal cycles
       
        # missing coverage due to snow + ice + water
        Frac_no_dep = FRSNO_Am[mm, :]+ FROCEAN_Am + FRLANDIC_Am
    
        DF_m[pp, mm,:] = DV_m[pp, mm, :] * Hg0_conc_mm_Am * (1. - Frac_no_dep) # calculate deposition flux

DF_ann = np.mean(DF_m, axis=1) # take annual average
unit_conv = MW_Hg / avo * ug_g * cm2_m2 * s_yr
DF_ann = DF_ann * unit_conv # convert to ug/m2/yr

Amazon_dep = np.zeros(num)  # for listing Amazon totals
for pp in range(num):
    Amazon_dep[pp] = np.sum(DF_ann[pp,:] * gbox_GC_Am * Am_mask_Am) / Mg_ug # Mg/yr
    print('Amazon dry deposition: p' + str(pp))
    print(Amazon_dep[pp])

# save data to csv, for access
np.savetxt("misc_Data/unc_Amazon_dep_HIST.csv", Amazon_dep, delimiter=",")

t1 = time.time()

total_t = t1-t0
print('Execution time in seconds: ' + str(total_t))
# save data to mat, for easier access
sio.savemat('misc_Data/unc_drydepflux_m_HIST_Am.mat', {"DF": DF_ann})
