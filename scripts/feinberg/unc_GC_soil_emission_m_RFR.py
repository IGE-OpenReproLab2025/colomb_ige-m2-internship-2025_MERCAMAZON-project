#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 5
Run soil emissions at monthly time resolution for parametrization, but average diurnal cycle
Run RFR scenario
GEOS-Chem input data can be downloaded at: http://geoschemdata.wustl.edu/ExtData/
@author: arifeinberg
"""
import numpy as np
import datetime
import xarray as xr
import numbers
import pandas as pd
from drydep_functions import Compute_Olson_landmap
import time
from scipy.io import savemat

def weird_division1(n, d):
    return n / d if d else np.float32(0) # no change if denominator 0

# create ufunc from weird division1
udiv1 = np.frompyfunc(weird_division1, 2, 1)


#%% functions for soil emissions
def Khan_soil_emiss(Hg_soil, SWGDN, nosnow, SunCos,  LAI, params):
    """Function from Khan et al. (2019) paper for soil emiss
    
    Parameters
    ----------
    Hg_soil : numpy array or single number
        Concentration of Hg in soils (ng Hg /g)
                
    SWGDN : numpy array or single number
        Solar radiation at ground (W/m2)

    nosnow : numpy array or single number
        fraction of snow coverage

    SunCos : numpy array or single number
        cos of solar zenith angle
        
    LAI : numpy array or single number
        leaf area index (m2/m2)
        
    params : numpy array or list
        parameters used in Khan soil emission routine
    """
    soil_emiss_fac = params[0] # prefactor for soil emissions
    mu = params[1] # parameter for LAI shading
    a = params[2] # exponent for soil conc
    b = params[3] # exponent for radiation


    # attenuate solar radiation based on function of leaf area index
    tauz = LAI * mu 
    
    # ensure that suncos is not negative/zero
    suncosvalue = SunCos
    if isinstance(suncosvalue, np.ndarray): # if it is a numpy array
        suncosvalue[suncosvalue<0.09] = 0.09
    elif isinstance(suncosvalue, numbers.Number): # if it is a single number
        suncosvalue = max(suncosvalue,0.09)
    else:
        print("Error1")
    
    # fraction of light reaching the surface is attenuated based on LAI
    lightfrac = np.exp(-tauz / suncosvalue)
    
    # adjust solar radiation by canopy attenuation
    R_g_adj = SWGDN * lightfrac
    
    # change units of Hg soil to ug/g
    Csoil = Hg_soil / 1000.0
    
    # calculate soil emissions in ng/m2/h
    soil_emis_Hg = soil_emiss_fac * (Csoil**a) * (R_g_adj**b)
    
    # calculate soil emissions in ug/m2/yr
    ng_ug = 1000. # ng in ug
    h_yr = 24. * 365.2425 # h in yr
    unit_conv = h_yr /ng_ug
    soil_emiss_Hg = soil_emis_Hg * unit_conv
    
    # only consider flux over snow free land
    soil_emiss_Hg = soil_emiss_Hg * nosnow
    return soil_emiss_Hg

#%% load soil concentrations of Hg
fn_soil = '../emissions/SOIL/soilHg.presentday.v11-01.geosfp.2x25.nc'
ds_soil = xr.open_dataset(fn_soil)

# coords
lon = ds_soil.lon # longitude
lat = ds_soil.lat # latitude

lon_l = len(lon)
lat_l = len(lat)

# scaling factor map
Hg_factor = ds_soil.HG0_DIST # unitless

# unit conversion
preind_Hg = 45.0 # ng/g preindustrial Hg concentration in soils
Hg_soilconc = Hg_factor.squeeze() * preind_Hg


# load snow cover, for use in emission flux calculations
# load land, non glaciated
fn_fr = '../misc_datasets/MERRA2_proc/MERRA2.20150101.CN.2x25.nc4'
ds_fr = xr.open_dataset(fn_fr)

FRLAND = ds_fr.FRLAND.squeeze() # land areas excluding lakes

# load snow coverage
fn_sno = '../misc_datasets/MERRA2_proc/MERRA2_FRSNO_2015_2x25_m.nc4'
ds_sno = xr.open_dataset(fn_sno)

FRSNO = ds_sno.FRSNO# snow coverage on land, monthly
time_t = FRSNO.time
time_l = len(time_t)

# calculate snow free land
f_nosno = FRLAND - FRSNO
f_nosno = xr.where(f_nosno>0, f_nosno, 0)


#%% Convert necessary variables to np arrays
Hg_soil_v = Hg_soilconc.values
f_nosno_v = f_nosno.values

#%% load uncertainty parameter values
# columns: 1) soil parametrization 2) LAI_pct_replace 3)- f0 Amazon RF 
#          4) f0 other RF 5) f0 elsewhere 6) biomass burning emissions
params_unc = np.loadtxt("misc_Data/LHS_sampling_dep_emis.csv",
                 delimiter=",")
num = len(params_unc) # number of uncertainty calculations

# Load parameters for soil parametrization
# columns: 1- prefactor 2- e-soil 3- e-rad
params_unc_soil = np.loadtxt("misc_Data/unc_soil_params_100_v6.csv",
                 delimiter=",")
num_soil = len(params_unc_soil) # number of soil parametrizations
# columns: 1- prefactor 2- e-soil 3- e-rad

#%% load area and masks
# Load grid cell area for unit conversion of model
fn_gbox = '../gbox_areas/GEOSChem_2x25_gboxarea.nc'
ds_gbox = xr.open_dataset(fn_gbox)
gbox_GC = ds_gbox.cell_area.values #m2
Mg_ug = 1e12 # Mg in ug 

#%% loop over parameter combinations, hours, and land types, calculate soil emissions
startTime = time.time()


days_in_month = np.asarray([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]) # days in each month
Amazon_emiss = np.zeros(num) # Amazon soil emissions for the simulation case
global_emiss = np.zeros(num) # Global soil emissions for the simulation case
Hg_soil_emiss_ann = np.zeros((num, lat_l,lon_l)) # Global soil emissions map for the simulation case

for j in range(num):
    print("j" + str(j))
    Hg_soil_emiss_m = np.zeros((time_l, lat_l, lon_l),  np.float32)
    # select correct parametrization for soil
    i_soil = int(params_unc[j,0]) - 1 # for python indexing
    params = [params_unc_soil[i_soil,0], 0.5, params_unc_soil[i_soil,1], params_unc_soil[i_soil,2]] # parameters used
    
    # load LAI
   # Land type LAI - precalculated for each RFR uncertainty run here: unc_RFR_LAI_nc.py
    fn_xlai  = 'misc_Data/RFR/unc_LAI/Yuan_proc_MODIS_XLAI.2x25.RFR_unc_' + str(j) + '_m.nc'
    ds_xlai = xr.open_dataset(fn_xlai)
    
    # saved as one data array, first dim is land types
    XLAI_m = ds_xlai.XLAI.values
    
    # load Olson land types
    # Land type areas
    fn_ols2 = 'misc_Data/RFR/unc_LAI/Olson_Land_Type_Masks.2x25.RFR_unc_' + str(j) + '.nc'
    ds_ols2 = xr.open_dataset(fn_ols2)
    # save as one np array, first dim is land types
    Olson_landtype = ds_ols2.Olson_LT.squeeze()
    # get additional information about Olson land types
    ltype = Olson_landtype["variable"]
    lt_l = len(ltype) # number of surface types
    
    Olson_landtype_v = Olson_landtype.values

    # loop over monthly timesteps        
    for k in range(time_l): 
        #print(k)
        # load solar radiation from GEOS Chem
        # monthly diurnal data
        fn_met_m = '../GC_data/run0311/Soil_emiss_2015_' +str(k+1).zfill(2) + '_diurnal.nc4'
        ds_met_m = xr.open_dataset(fn_met_m)
        SWGDN_m = ds_met_m.Met_SWGDN.values
        SunCos_m = ds_met_m.Met_SUNCOSmid.values
        
        # LAI divided by areal fraction to get true LAI 
        LAI_m_l = np.array(udiv1(XLAI_m[:,k,:,:], Olson_landtype_v),dtype=np.float32)
        Hg_soil_emiss_h = np.zeros((24, lat_l, lon_l),  np.float32)

        for i in range(24): # diurnal cycle per month
            for l in range(lt_l): # loop over land types
                # calculate map of soil emiss for that hour
                Hg_soil_emiss_h_l = Khan_soil_emiss(Hg_soil_v[:,:], SWGDN_m[i,:,:],\
                                                         f_nosno_v[:,:,k], SunCos_m[i,:,:], \
                                                         LAI_m_l[l,:,:], params)
                # sum emiss over all land types, weighted by their area and the days in month
                Hg_soil_emiss_h [i,:,:] =  Hg_soil_emiss_h [i,:,:] + \
                    Hg_soil_emiss_h_l * Olson_landtype_v[l,:,:]
        # take average of diurnal cycle for monthly average, weight by days in month
        Hg_soil_emiss_m[k,:,:] = np.mean(Hg_soil_emiss_h, axis=0) * days_in_month[k]

    # calculate global and Amazon emissions
    Hg_soil_emiss_ann[j,:,:] = np.sum(Hg_soil_emiss_m, axis=0) / sum(days_in_month) # take annual mean, weighted by days in month
    global_emiss[j] = np.sum(Hg_soil_emiss_ann[j,:,:] * gbox_GC) / Mg_ug # Mg/yr
    print(global_emiss[j])


# save data to csv, for access
np.savetxt("misc_Data/unc_global_emiss_RFR.csv", global_emiss, delimiter=",")

# save data to .mat
Hg_soil_emiss_ann_u = np.float32(Hg_soil_emiss_ann) # convert to smaller units for storage
savemat("misc_Data/unc_emiss_map_RFR.mat", {"Hg_soil_emiss_ann": Hg_soil_emiss_ann_u})

executionTime = (time.time() - startTime)
print('Execution time in seconds: ' + str(executionTime))
