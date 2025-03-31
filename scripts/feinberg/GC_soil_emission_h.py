#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  13 16:02:55 2022
Calculate soil emissions offline at hourly time resolution over a year 
GEOS-Chem input data can be downloaded at: http://geoschemdata.wustl.edu/ExtData/
@author: arifeinberg
"""
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import scipy.io as sio
import datetime
import xarray as xr
import numbers
import pandas as pd
from drydep_functions import Compute_Olson_landmap

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

# load solar radiation from GEOS Chem
# hourly data
fn_met_h = '../GEOS-Chem_runs/run0012/OutputDir/Soil_emiss_2015_inputs.nc4'
ds_met_h = xr.open_dataset(fn_met_h)
SWGDN_h = ds_met_h.Met_SWGDN
SunCos_h = ds_met_h.Met_SUNCOSmid

# load snow cover, for use in emission flux calculations
# load land, non glaciated
fn_fr = '../misc_datasets/MERRA2_proc/MERRA2.20150101.CN.2x25.nc4'
ds_fr = xr.open_dataset(fn_fr)

FRLAND = ds_fr.FRLAND.squeeze() # land areas excluding lakes

# load snow coverage
fn_sno = '../misc_datasets/MERRA2_proc/MERRA2_FRSNO_2015_2x25_h.nc4'
ds_sno = xr.open_dataset(fn_sno)

FRSNO = ds_sno.FRSNO # snow coverage on land, monthly

# calculate snow free land
f_nosno = FRLAND - FRSNO
f_nosno = xr.where(f_nosno>0, f_nosno, 0)

# load LAI
# Land type LAI
fn_xlai  = 'misc_Data/Yuan_proc_MODIS_XLAI.2x25.BAU_2002.nc'
ds_xlai = xr.open_dataset(fn_xlai)

# save as one data array, first dim is land types
XLAI_w = ds_xlai.to_array()

# calculate daily values of XLAI 
temp = XLAI_w.reindex(time=pd.date_range("1/1/2002", "12/31/2002")) 
XLAI_d = temp.ffill("time").bfill("time").values # fill forwards, extrapolate first day by filling backwards

time = SunCos_h.time
time_l = len(time)
#%% Convert necessary variables to np arrays
Hg_soil_v = Hg_soilconc.values
SWGDN_v = SWGDN_h.values
f_nosno_v = f_nosno.values
SunCos_v = SunCos_h.values

#%% load Olson land types
# Land type areas
fn_ols2 = 'misc_Data/Olson_Land_Type_Mask.2_25.BAU_2002.nc'
ds_ols2 = xr.open_dataset(fn_ols2)
# save as one np array, first dim is land types
Olson_landtype = ds_ols2.to_array().squeeze()
# get additional information about Olson land types
FRCLND, IREG, ILAND, IUSE = Compute_Olson_landmap(Olson_landtype)
IREG_v = IREG.values #  number of types in grid box
ILAND_v = ILAND.values # index of types actually in grid box
IUSE_v = IUSE.values # fractional coverage of types actually in grid box in per mil
ltype = Olson_landtype["variable"]
lt_l = len(ltype) # number of surface types

Olson_landtype_v = Olson_landtype.values

#%% loop over hours, and land types, calculate soil emissions
Hg_soil_emiss_h = np.zeros((time_l, lat_l, lon_l),  np.float32)
Hg_soil_emiss_h_l = np.zeros((lt_l, time_l, lat_l, lon_l),  np.float32)
params = [71., 0.5, 2.5, 0.76] # final tuning

days_in_month = np.asarray([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]) # days in each month

# loop over hourly timesteps        
for k in range(time_l): 
    if (k % 50 ==0):
        print(k)

    # calculate days of year
    day_mo = time[k].dt.day.values.item() - 1  # for python indexing, date of month
    mo_no = time[k].dt.month.values.item() - 1 # for python indexing, month number
    day_no = sum(days_in_month[:mo_no]) + day_mo  # day of year 

    # select correct LAI day    
    XLAI_h =  XLAI_d[:,day_no,:, :] #1.5
    LAI_h_l = np.array(udiv1(XLAI_h, Olson_landtype_v),dtype=np.float32) # LAI divided by areal fraction to get true LAI 
    for l in range(lt_l): # loop over land types
      # calculate map of soil emiss for that hour
      Hg_soil_emiss_h_l[l,k,:,:] = Khan_soil_emiss(Hg_soil_v[:,:], SWGDN_v[k,:,:],\
                                               f_nosno_v[:,:,k], SunCos_v[k,:,:], \
                                               LAI_h_l[l,:,:], params)
      # sum emiss over all land types, weighted by their area    
      Hg_soil_emiss_h [k,:,:] =  Hg_soil_emiss_h [k,:,:] + \
          Hg_soil_emiss_h_l[l,k,:,:] * Olson_landtype_v[l,:,:]

# save data to mat, for easier access
fn_output = 'misc_Data/soil_emiss_h_2002LAI.mat'
sio.savemat(fn_output, {"Hg_soil_emiss_h": Hg_soil_emiss_h,"params": params})
