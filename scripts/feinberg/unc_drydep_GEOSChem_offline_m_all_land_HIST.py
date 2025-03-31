#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 2023
Running offline version of dry deposition of GEOS-Chem on monthly timescale, hourly average diurnal cycle
For all land boxes HIST scenario, to calculate reference totals to compare with reforestation scenario
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

#%% Load land masks
fn_fr = '../misc_datasets/MERRA2_proc/MERRA2_FROCEAN_FRLANDIC_2x25.nc4'
ds_fr = xr.open_dataset(fn_fr)
FRLANDIC = ds_fr.FRLANDIC.squeeze().values # land glaciated areas
FROCEAN = ds_fr.FROCEAN.squeeze().values # ocean areas
# load snow cover, for use in dry deposition flux calculations
fn_sno = '../misc_datasets/MERRA2_proc/MERRA2_FRSNO_2015_2x25_m.nc4'
ds_sno = xr.open_dataset(fn_sno)
FRSNO = ds_sno.FRSNO.values # snow coverage on land, monthly

# make mask of all land areas
no_dep = FROCEAN + FRLANDIC

mask_land = xr.where(no_dep<1, 1, 0) 

# Get the array indices that need to calculate for land
la_calc_bool = np.where(mask_land>0, True, False)

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
gbox_GC_la = gbox_GC[la_calc_bool]
n_la = len(gbox_GC_la) # number of grid cells to calculate
IREG_v_la = IREG_v[la_calc_bool]
FROCEAN_la = FROCEAN[la_calc_bool]
FRLANDIC_la = FRLANDIC[la_calc_bool]
in_AMZN_la = in_AMZN[la_calc_bool]

#3D vars - need to be 2D array in the output, so reshape
la_calc_bool3D = np.broadcast_to(la_calc_bool, Olson_landtype.shape)# need 3D mask for Olson
ILAND_v_la = np.reshape(ILAND_v[la_calc_bool3D], (73, n_la))
IUSE_v_la = np.reshape(IUSE_v[la_calc_bool3D], (73, n_la))

#time array - need to be 2D array in the output, so reshape
la_calc_boolt = np.broadcast_to(la_calc_bool, FRSNO.shape)# need 3D mask for Olson
FRSNO_la = np.reshape(FRSNO[la_calc_boolt], (12, n_la))

# initialize
DV_m = np.zeros((num, 12,n_la)) # array of deposition velocity maps
DF_m = np.zeros((num, 12,n_la)) # array of deposition flux maps

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
    la_calc_bool_24H = np.broadcast_to(la_calc_bool, bxheight.shape)# need 3D mask for 24H arrays

    bxheight_la = np.reshape(bxheight[la_calc_bool_24H], (24, n_la))
    surf_pres_la = np.reshape(surf_pres[la_calc_bool_24H], (24, n_la))
    Z0_la = np.reshape(Z0[la_calc_bool_24H], (24, n_la)) # surface roughness
    CLDFRC_la = np.reshape(CLDFRC[la_calc_bool_24H], (24, n_la)) # column cloud fraction
    albedo_la = np.reshape(albedo[la_calc_bool_24H], (24, n_la)) # visible surface albedo
    airden_la = np.reshape(airden[la_calc_bool_24H], (24, n_la)) # air density at surface
    hflux_la = np.reshape(hflux[la_calc_bool_24H], (24, n_la)) # sensible heat flux
    sw_grnd_la = np.reshape(sw_grnd[la_calc_bool_24H], (24, n_la)) # incident shortwave at ground
    surf_t_la = np.reshape(surf_t[la_calc_bool_24H], (24, n_la)) # surface temperature 
    ustar_la = np.reshape(ustar[la_calc_bool_24H], (24, n_la)) # friction velocity
    U10M_la = np.reshape(U10M[la_calc_bool_24H], (24, n_la)) # zonal wind at 10 m height
    V10M_la = np.reshape(V10M[la_calc_bool_24H], (24, n_la)) # meridional wind at 10 m height
    SUNCOSmid_la = np.reshape(SUNCOSmid[la_calc_bool_24H], (24, n_la)) # cosine of solar zenith angle, at midpoint of chemistry timestep
    
    Hg0_conc_mm_la = Hg0_conc_mm[la_calc_bool] # Hg0 conc in molec/cm3

    # XLAI, 4 D variable
    XLAI_mm = XLAI_m[:,mm,:,:]
    XLAI_mm_la = np.reshape(XLAI_mm[la_calc_bool3D], (73, n_la))
    # run loop over parameter values
    for pp in range(num): # loop over parameters
        # set F0 parameters for this run
        F0 = params_unc[pp,4] # f0 elsewhere
        F0_AMZN = params_unc[pp,2] # f0 in Amazon rainforest
        F0_RF = params_unc[pp,3] # f0 in other rainforests
        
        DV_h = np.zeros((24,n_la)) # array of hourly deposition velocity maps

        for i in range(24): # loop over hours
            for jj in range(n_la):
                if (FROCEAN_la[jj] + FRLANDIC_la[jj] <1): # don't calculate over water or ice
                    CZ1, LSNOW, OBK, W10 = METERO(bxheight_la[i,jj],\
                                                       albedo_la[i,jj], \
                                                       surf_t_la[i,jj], \
                                                       ustar_la[i,jj], \
                                                       airden_la[i,jj], \
                                                       hflux_la[i,jj], \
                                                       U10M_la[i,jj], \
                                                       V10M_la[i,jj]) 
                            
                    DV_h[i,jj] = DEPVEL_UNC(DRYCOEFF, IOLSON, IDEP, IRI, \
                               IRLU, IRAC, IRGSS, IRGSO, IRCLS, IRCLO, \
                               IREG_v_la[jj],\
                               ILAND_v_la[:,jj],\
                               IUSE_v_la[:,jj], \
                               surf_t_la[i,jj], \
                               XLAI_mm_la[:,jj], \
                               LSNOW, sw_grnd_la[i,jj], CLDFRC_la[i,jj], \
                               SUNCOSmid_la[i,jj], surf_pres_la[i,jj], \
                               ustar_la[i,jj], Z0_la[i,jj], \
                                   CZ1, OBK, XMW, F0, F0_AMZN, F0_RF, in_AMZN_la[jj], HSTAR)
        DV_m[pp,mm,:] = np.mean(DV_h, axis = 0) # take monthly average of diurnal cycles
   
        # missing coverage due to snow + ice + water
        Frac_no_dep = FRSNO_la[mm, :]+ FROCEAN_la + FRLANDIC_la

        DF_m[pp,mm,:] = DV_m[pp,mm, :] * Hg0_conc_mm_la * (1. - Frac_no_dep) # calculate deposition flux
    

DF_ann = np.mean(DF_m, axis=1) # take annual average
unit_conv = MW_Hg / avo * ug_g * cm2_m2 * s_yr
DF_ann = DF_ann * unit_conv # convert to ug/m2/yr

DF_ann_u = np.float32(DF_ann) # convert to smaller units for storage
land_dep = np.zeros(num)  # for listing global totals
for pp in range(num):
    land_dep[pp] = np.sum(DF_ann[pp,:] * gbox_GC_la ) / Mg_ug # Mg/yr
    print('Global dry deposition: p' + str(pp))
    print(land_dep[pp])

# save data to csv, for access
np.savetxt("misc_Data/unc_dep_HIST_global_all_land.csv", land_dep, delimiter=",")
# save data to mat, for easier access
sio.savemat('misc_Data/unc_drydepflux_m_HIST_all_land.mat', {"DF": DF_ann_u})

t1 = time.time()

total_t = t1-t0
print('Execution time in seconds: ' + str(total_t))
