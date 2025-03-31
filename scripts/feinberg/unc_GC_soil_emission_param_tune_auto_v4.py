#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Apr 17 2023
Automatically find soil emission tuning parameters for certain set of constraints - use LAI data for inputs, monthly averages with diurnality
GEOS-Chem input data can be downloaded at: http://geoschemdata.wustl.edu/ExtData/
@author: arifeinberg
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
import datetime
import xarray as xr
import numbers
from smt.sampling_methods import LHS
import time

#%% functions for soil emissions
def weird_division1(n, d):
    return n / d if d else np.float32(0) # no change if denominator 0

# create ufunc from weird division1
udiv1 = np.frompyfunc(weird_division1, 2, 1)

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
    
    soil_emiss_fac = params[0] # (10**0.709) / 1.5 # another factor for soil emissions
    mu = params[1] # 0.5 # constant for LAI shading
    a = params[2] # 0.119 # exponent for soil conc
    b = params[3] # 0.137 # exponent for radiation

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

def run_monthly_diurnal(Hg_soil, nosnow, XLAI, Olson_landtype, params):
    """calculating annually the emissions, each month average diurnal cycle
    
    Parameters
    ----------
    Hg_soil : numpy array or single number
        Concentration of Hg in soils (ng Hg /g)
                
    nosnow : numpy array or single number
        fraction of snow coverage
        
    XLAI : numpy array 
        leaf area index for each land type (m2/m2)
        
    Olson_landtype : numpy array 
        Olson landtypes for each grid cell (fraction of grid cell)
    
    params : numpy array or list
        parameters used in Khan soil emission routine
    """
    days_in_month = np.asarray([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]) # days in each month
    lt_l = np.shape(Olson_landtype)[0] # number of land types
    lat_l = np.shape(Olson_landtype)[1] # number of lats
    lon_l = np.shape(Olson_landtype)[2] # number of lons

    Hg_soil_emiss_m = np.zeros((12, lat_l, lon_l),  np.float32) # initialize

    # loop over monthly timesteps        
    for k in range(12): 
        #print(k)
        # load solar radiation from GEOS Chem
        # monthly diurnal data
        fn_met_m = '../GC_data/run0311/Soil_emiss_2015_' +str(k+1).zfill(2) + '_diurnal.nc4'
        ds_met_m = xr.open_dataset(fn_met_m)
        SWGDN_m = ds_met_m.Met_SWGDN.values # Solar radiation at ground (W/m2)
        SunCos_m = ds_met_m.Met_SUNCOSmid.values # cos of solar zenith angle
        
        # LAI divided by areal fraction to get true LAI 
        LAI_m_l = np.array(udiv1(XLAI[:,k,:,:], Olson_landtype),dtype=np.float32)
        Hg_soil_emiss_h = np.zeros((24, lat_l, lon_l),  np.float32)

        for i in range(24): # diurnal cycle per month
            for l in range(lt_l): # loop over land types
                # calculate map of soil emiss for that hour
                Hg_soil_emiss_h_l = Khan_soil_emiss(Hg_soil, SWGDN_m[i,:,:],\
                                                         nosnow[k,:,:], SunCos_m[i,:,:], \
                                                         LAI_m_l[l,:,:], params)
                # sum emiss over all land types, weighted by their area and the days in month
                Hg_soil_emiss_h [i,:,:] =  Hg_soil_emiss_h [i,:,:] + \
                    Hg_soil_emiss_h_l * Olson_landtype[l,:,:]
        # take average of diurnal cycle for monthly average, weight by days in month
        Hg_soil_emiss_m[k,:,:] = np.mean(Hg_soil_emiss_h, axis=0) * days_in_month[k]
    # return annual emissions map
    Hg_soil_emiss_ann = np.sum(Hg_soil_emiss_m, axis=0) / sum(days_in_month) # take annual mean, weighted by days in month
    return Hg_soil_emiss_ann
#%% load soil concentrations of Hg
fn_soil = '../emissions/SOIL/soilHg.presentday.v11-01.geosfp.2x25.nc'
ds_soil = xr.open_dataset(fn_soil)

# coords
lon = ds_soil.lon # longitude
lat = ds_soil.lat # latitude

# scaling factor map
Hg_factor = ds_soil.HG0_DIST.values # unitless

# unit conversion
preind_Hg = 45.0 # ng/g preindustrial Hg concentration in soils
Hg_soilconc = Hg_factor.squeeze() * preind_Hg

# load solar radiation from GEOS Chem

# monthly data
fn_met = '../GEOS-Chem_runs/run0012/OutputDir/Soil_emiss_2015_inputs_m.nc4'
ds_met = xr.open_dataset(fn_met)

SWGDN = ds_met.Met_SWGDN.values
SWGDN_ann = SWGDN # monthly mean
SunCos = ds_met.Met_SUNCOSmid.values
SunCos_ann = SunCos # monthly mean

# load snow cover, for use in emission flux calculations
# load land, non glaciated
fn_fr = '../misc_datasets/MERRA2_proc/MERRA2.20150101.CN.2x25.nc4'
ds_fr = xr.open_dataset(fn_fr)

FRLAND = ds_fr.FRLAND.squeeze().values # land areas excluding lakes

# load snow coverage
fn_sno = '../misc_datasets/MERRA2_proc/MERRA2_FRSNO_2015_2x25_m.nc4'
ds_sno = xr.open_dataset(fn_sno)

FRSNO = ds_sno.FRSNO.values # snow coverage on land, monthly
FRSNO_ann = FRSNO # monthly mean

# calculate snow free land
f_nosno = FRLAND - FRSNO_ann
f_nosno = np.where(f_nosno>0, f_nosno, 0)

# load LAI for REF
fn_xlai  = 'misc_Data/Yuan_proc_MODIS_XLAI.2x25.BAU_2002_m.nc'
ds_xlai = xr.open_dataset(fn_xlai)

# save as one data array, first dim is land types
XLAI_HIST = ds_xlai.drop_vars('time_bnds').to_array().values # monthly mean

# load LAI for SAV
fn_xlai  = 'misc_Data/Yuan_proc_MODIS_XLAI.2x25.SAV_58_m.nc'
ds_xlai = xr.open_dataset(fn_xlai)

# save as one data array, first dim is land types
XLAI_SAV = ds_xlai.drop_vars('time_bnds').to_array().values # monthly mean
#%% load Olson land types
# Land type areas
fn_ols2 = 'misc_Data/Olson_Land_Type_Mask.2_25.BAU_2002.nc'
ds_ols2 = xr.open_dataset(fn_ols2)
# save as one np array, first dim is land types
Olson_landtype_HIST = ds_ols2.to_array().squeeze().values

fn_ols2 = 'misc_Data/Olson_Land_Type_Mask.2_25.SAV.nc'
ds_ols2 = xr.open_dataset(fn_ols2)
# save as one np array, first dim is land types
Olson_landtype_SAV = ds_ols2.to_array().squeeze().values

#%% Limits of constraints, setting up experimental design
# limits of observable constraints
rat_obs_df_f_lims = [1.8, 31] # 6.7 media ratio of deforest to forest
rat_obs_AM_ET_lims = [3.5, 8] # 5.3 ratio of Amazon to extratropical grasslands
obs_ET_df_lims = [3.5, 11.4] # 4.3 # value emissions Extratropical grasslands

# 23 median
obs_AM_df_max = 79 # max value emissions Amazon deforested
obs_AM_df_min = 9.8 # min value emissions Amazon deforested

# sample using Latin Hypercube Sampling
xlimits = np.array([rat_obs_df_f_lims, rat_obs_AM_ET_lims, obs_ET_df_lims])
sampling = LHS(xlimits=xlimits) # already sampled, just to get params)

num = 100
x = sampling(num) # already sampled, just to get params

#np.savetxt("misc_Data/LHS_sampling_x.csv", x, delimiter=",") # already sampled, just to get params
#x = np.loadtxt("misc_Data/LHS_sampling_x.csv",
#                 delimiter=",")

# initialize array
params_all = np.zeros((3,num))
Amazon_Hg_v_a = np.zeros(num)

#%% run loop to calculate parameters for soil emissions
startTime = time.time()

for i in range(num):
    print(i)
    rat_obs_df_f = x[i,0] # ratio of deforest to forest
    rat_obs_AM_ET = x[i,1]# ratio of Amazon to extratropical grasslands
    obs_ET_df = x[i,2] # value emissions Extratropical grasslands
    
    shading = 0.5 # keep this fixed
    
    # STEP 1: solve for e_rad, since that is just sensitive to ratio of deforest to forest
    n = 5 # number of parameter options
    e_rad = np.linspace(0.1, 3, n) # range of reasonable options
    prefactor_c = 71 # can keep constant for this step
    e_soil_c = 2.5 # can keep constant for this step
    # initialize arrays
    Amazon_Hg_v = np.zeros(n)
    NH_ml_Hg_v = np.zeros(n)
    Amazon_Hg_RF_v = np.zeros(n)
    # amazon limited bounding box
    ilt_l_A = 39 # -12. 
    ilt_h_A = 45 # 0.
    iln_l_A = 43 #-72.
    iln_h_A = 52 #-50.
    
    #  used for constraining Northern Hemisphere midlatitudes
    ilt_l_N = 60 # 30. 
    ilt_h_N = 70 # 50.
    # create loops
    for l in range(n): # e_rad
        params = [prefactor_c, shading, e_soil_c, e_rad[l]]
        # calculate soil map - with deforestation
        s_map = run_monthly_diurnal(Hg_soilconc, f_nosno, XLAI_SAV, Olson_landtype_SAV, params)
         # calculate soil map -  with normal LAI conditions
        s_map_RF = run_monthly_diurnal(Hg_soilconc, f_nosno, XLAI_HIST, Olson_landtype_HIST, params)
       
        # restrict to Amazon
        # (savannah, cropland, grassland etc)
        Amazon_Hg = s_map[ilt_l_A:ilt_h_A+1, iln_l_A:iln_h_A+1]
        Amazon_Hg_v[l] = np.median(Amazon_Hg)
        # (rainforest)
        Amazon_Hg_RF = s_map_RF[ilt_l_A:ilt_h_A+1, iln_l_A:iln_h_A+1]
        Amazon_Hg_RF_v[l] = np.median(Amazon_Hg_RF)
    
    # calculate ratio of deforested vs. forested emissions
    ratio_RF_forest = Amazon_Hg_v / Amazon_Hg_RF_v
    
    # Find value of e_rad that fits observed constraint
    e_rad_tune = np.interp(rat_obs_df_f, ratio_RF_forest, e_rad)
    
    # STEP 2: solve for e_soil, since this is sensitive to ratio of NH to Amazon
    e_soil = np.linspace(1., 6., n) # range of reasonable options
    for k in range(n): # e_rad
        params = [prefactor_c, shading, e_soil[k], e_rad_tune]
        # calculate soil map - with deforestation
        s_map = run_monthly_diurnal(Hg_soilconc, f_nosno, XLAI_SAV, Olson_landtype_SAV, params)
       
        # restrict to Amazon
        # LAI=(savannah, cropland, grassland etc)
        Amazon_Hg = s_map[ilt_l_A:ilt_h_A+1, iln_l_A:iln_h_A+1]
        Amazon_Hg_v[k] = np.median(Amazon_Hg)
        
        # restrict to NH midlatitudes
        NH_Hg = s_map[ilt_l_N:ilt_h_N+1, :]
        NH_ml_Hg_v[k]= np.median(NH_Hg[NH_Hg>0.])
    
    # calculate ratio of Amazon vs. extratropical emissions
    ratio_Amazon_NH = Amazon_Hg_v / NH_ml_Hg_v # Amazon deforested to NH midlats
    # Find value of e_rad that fits observed constraint
    e_soil_tune = np.interp(rat_obs_AM_ET, ratio_Amazon_NH, e_soil)
    #print(e_soil_tune)
    #  STEP 3: calculate the prefactor 
    params = [prefactor_c, shading, e_soil_tune, e_rad_tune]
    s_map = run_monthly_diurnal(Hg_soilconc, f_nosno, XLAI_SAV, Olson_landtype_SAV, params)
      
    # restrict to Amazon
    # LAI (savannah, cropland, grassland etc)
    Amazon_Hg = s_map[ilt_l_A:ilt_h_A+1, iln_l_A:iln_h_A+1]
    Amazon_Hg_v = np.median(Amazon_Hg)
    
    # restrict to NH midlatitudes
    NH_Hg = s_map[ilt_l_N:ilt_h_N+1, :]
    NH_ml_Hg_v= np.median(NH_Hg[NH_Hg>0.])
    
    # tune prefactor by ratio of obs constraint and current NH extratrop value
    prefactor_tune = obs_ET_df / NH_ml_Hg_v * prefactor_c
    #print(prefactor_tune)
    
    # STEP 4: check that the prefactor doesn't yield overly large Amazon emissions
    params = [prefactor_tune, shading, e_soil_tune, e_rad_tune]
    s_map = run_monthly_diurnal(Hg_soilconc, f_nosno, XLAI_SAV, Olson_landtype_SAV, params)
    
    Amazon_Hg = s_map[ilt_l_A:ilt_h_A+1, iln_l_A:iln_h_A+1]
    Amazon_Hg_v = np.median(Amazon_Hg)
    
    Amazon_Hg_v_a[i] = Amazon_Hg_v
    print([prefactor_tune, e_soil_tune, e_rad_tune])
    # check if reasonable value
    if (Amazon_Hg_v < obs_AM_df_max and Amazon_Hg_v < obs_AM_df_max):
        params_all[:,i] = [prefactor_tune, e_soil_tune, e_rad_tune]
    else:
        params_all[:,i] = [np.nan, np.nan, np.nan]
        print("Fail")
executionTime = (time.time() - startTime)
print('Execution time in seconds: ' + str(executionTime))

#%% Count number of fails
# count zeros in array
# n_zeros = np.sum(np.isnan(params_all))
# display the count of zeros
# print("Number of fails:")
# print(n_zeros/3)
# print("Number Total")
# print(num)
# print("Param1 bounds")
# print(np.nanmin(params_all[0,:]))
# print(np.nanmax(params_all[0,:]))
# print("Param2 bounds")
# print(np.nanmin(params_all[1,:]))
# print(np.nanmax(params_all[1,:]))
# print("Param3 bounds")
# print(np.nanmin(params_all[2,:]))
# print(np.nanmax(params_all[2,:]))
#%% save non-nan fails as csv
rows_nonnan = np.argwhere(~np.isnan(params_all[2,:])).squeeze()# find nonnan rows
params_all_nonan = np.transpose(params_all[:,rows_nonnan])
np.savetxt("misc_Data/unc_soil_params_100_v6.csv", params_all_nonan, delimiter=",")
