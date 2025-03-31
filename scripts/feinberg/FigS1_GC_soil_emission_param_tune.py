#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  9 16:02:55 2022
Tune soil emission parametrization from Khan19 to fit with observations
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

def GC_soil_emiss(Hg_soil, SWGDN, nosnow, SunCos,  LAI):
    """Function currently implemented in GEOS-Chem for soil emissions
    
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

    """
    soil_emiss_fac = 1.6e-2 * 0.9688 # another factor for soil emissions
    mu = 0.5 # constant for LAI shading
    gamma = 0.0011 # constant for emissions

    # attenuate solar radiation based on function of leaf area index
    tauz = LAI * mu 
    
    # ensure that suncos is not negative/zero
    suncosvalue = SunCos
    if isinstance(SunCos, np.ndarray): # if it is a numpy array
        suncosvalue[suncosvalue<0.09] = 0.09
    elif isinstance(SunCos, numbers.Number): # if it is a single number
        suncosvalue = max(suncosvalue,0.09)
    else:
        print("Error1")
    
    # fraction of light reaching the surface is attenuated based on LAI
    lightfrac = np.exp(-tauz / suncosvalue)
    
    # calculate soil emissions in ng/m2/h
    soil_emis_Hg = np.exp(gamma * SWGDN * lightfrac) * Hg_soil * soil_emiss_fac
    
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

# scaling factor map
Hg_factor = ds_soil.HG0_DIST # unitless

# unit conversion
preind_Hg = 45.0 # ng/g preindustrial Hg concentration in soils
Hg_soilconc = Hg_factor.squeeze() * preind_Hg

# load solar radiation from GEOS Chem
# daily data
fn_met_d = '../GEOS-Chem_runs/run0012/OutputDir/Soil_emiss_2015_inputs.nc4'
ds_met_d = xr.open_dataset(fn_met_d)
SWGDN_d = ds_met_d.Met_SWGDN.isel(time=slice(0,24))
SunCos_d = ds_met_d.Met_SUNCOSmid.isel(time=slice(0,24))

# monthly data
fn_met = '../GEOS-Chem_runs/run0012/OutputDir/Soil_emiss_2015_inputs_m.nc4'
ds_met = xr.open_dataset(fn_met)

SWGDN = ds_met.Met_SWGDN
SWGDN_ann = SWGDN.mean("time")
SunCos = ds_met.Met_SUNCOSmid
SunCos_ann = SunCos.mean("time")

# load snow cover, for use in emission flux calculations
# load land, non glaciated
fn_fr = '../misc_datasets/MERRA2_proc/MERRA2.20150101.CN.2x25.nc4'
ds_fr = xr.open_dataset(fn_fr)

FRLAND = ds_fr.FRLAND.squeeze() # land areas excluding lakes

# load snow coverage
fn_sno = '../misc_datasets/MERRA2_proc/MERRA2_FRSNO_2015_2x25_m.nc4'
ds_sno = xr.open_dataset(fn_sno)

FRSNO = ds_sno.FRSNO # snow coverage on land, monthly
FRSNO_ann = FRSNO.mean("time") # annual mean

# calculate snow free land
f_nosno = FRLAND - FRSNO_ann
f_nosno = xr.where(f_nosno>0, f_nosno, 0)


#%% loop over many parameter options, see how affects key constraints
n = 10 # number of parameter options
prefactor = np.linspace(10, 250., n)
shading = 0.5
e_soil = np.linspace(1, 3, n)
e_rad = np.linspace(0.1, 1.2, n)

# initialize arrays
Amazon_Hg_v = np.zeros((n,n,n))
NH_ml_Hg_v = np.zeros((n,n,n))
Amazon_Hg_RF_v = np.zeros((n,n,n))

# amazon limited bounding box
ilt_l_A = 39 # -12. 
ilt_h_A = 45 # 0.
iln_l_A = 43 #-72.
iln_h_A = 52 #-50.

# median NH midlat emiss
ilt_l_N = 60 # 30. 
ilt_h_N = 70 # 50.

# create loops
for i in range(n): # prefactor
    for k in range(n): # e_soil
        for l in range(n): # e_rad
            params = [prefactor[i], shading, e_soil[k], e_rad[l]]
            # calculate soil map - with deforestation
            s_map = Khan_soil_emiss(Hg_soilconc.values, SWGDN_ann.values, f_nosno.values, 0.5, 1.5, params)
             # calculate soil map -  with normal LAI conditions
            s_map_RF = Khan_soil_emiss(Hg_soilconc.values, SWGDN_ann.values, f_nosno.values, 0.5, 4, params)
           
            # restrict to Amazon
            # LAI=1.5 (savannah, cropland, grassland etc)
            Amazon_Hg = s_map[ilt_l_A:ilt_h_A+1, iln_l_A:iln_h_A+1]
            Amazon_Hg_v[i,k,l] = np.median(Amazon_Hg)
            # LAI=4 (rainforest)
            Amazon_Hg_RF = s_map_RF[ilt_l_A:ilt_h_A+1, iln_l_A:iln_h_A+1]
            Amazon_Hg_RF_v[i,k,l] = np.median(Amazon_Hg_RF)
            
            # restrict to NH midlatitudes
            NH_Hg = s_map[ilt_l_N:ilt_h_N+1, :]
            NH_ml_Hg_v[i,k,l]= np.median(NH_Hg[NH_Hg>0.])
                
# calculate ratios - better constraint
ratio_Amazon_NH = Amazon_Hg_v / NH_ml_Hg_v # Amazon deforested to NH midlats
# calculate ratio of deforested vs. forested emissions
ratio_RF_forest = Amazon_Hg_v / Amazon_Hg_RF_v
#%% Fix parameters at certain level, calculate prefactor
n = 10
shading_f = 0.5 # 0.5
e_rad_f = 0.76
e_soil_f = 2.5
prefactor_f = np.linspace(1, 300., n) #673 # 310 #220

Amazon_Hg_f = np.zeros(n)
Amazon_Hg_RF_f = np.zeros(n)
NH_ml_Hg_f = np.zeros(n)
# create loops
for i in range(n): # prefactor
    params = [prefactor_f[i], shading_f, e_soil_f, e_rad_f]
    # calculate soil map
    s_map = Khan_soil_emiss(Hg_soilconc.values, SWGDN_ann.values, f_nosno.values, 0.5, 1.5, params)
     # calculate soil map - LAI increased
    s_map_RF = Khan_soil_emiss(Hg_soilconc.values, SWGDN_ann.values, f_nosno.values, 0.5, 4, params)
   
    # restrict to Amazon
    # LAI=1.5 (savannah, cropland, grassland etc)
    Amazon_Hg = s_map[ilt_l_A:ilt_h_A+1, iln_l_A:iln_h_A+1]
    Amazon_Hg_f[i] = np.median(Amazon_Hg)
    # LAI=4 (rainforest)
    Amazon_Hg_RF = s_map_RF[ilt_l_A:ilt_h_A+1, iln_l_A:iln_h_A+1]
    Amazon_Hg_RF_f[i] = np.median(Amazon_Hg_RF)
    
    # restrict to NH midlatitudes
    NH_Hg = s_map[ilt_l_N:ilt_h_N+1, :]
    NH_ml_Hg_f[i]= np.median(NH_Hg[NH_Hg>0.])

#%% Plot clean figures
# plots for Amazon emiss
f, axes = plt.subplots(2, 2, figsize=[8,8],
                              gridspec_kw=dict(hspace=0.3, wspace=0.3))

axes = axes.flatten()

axes[0].plot(e_rad, np.mean(ratio_RF_forest, axis=(0,1)))
axes[0].set_ylabel('Ratio')
axes[0].set_xlabel('$c$ (radiation exponent)')
axes[0].axhline(6.7,c='k',ls='dashed')

axes[0].set_title('Ratio rainforest:deforest')

axes[1].plot(e_soil, np.mean(ratio_Amazon_NH, axis=(0,2)))
axes[1].set_ylabel('Ratio')
axes[1].set_xlabel('$b$ (soil exponent)')
axes[1].set_title('Ratio Amazon:Northern Hemisphere extratropics')
axes[1].axhline(5.3,c='k',ls='dashed')

axes[2].plot(prefactor_f, Amazon_Hg_f)
axes[2].axhline(23,c='k',ls='dashed')
axes[2].set_ylabel('Emiss (\u03bcg m$^{-2}$ yr$^{-1}$)')
axes[2].set_xlabel('$a$ (prefactor)')
axes[2].set_title('Amazon (deforested)')
  
axes[3].plot(prefactor_f, NH_ml_Hg_f)
axes[3].axhline(4.3,c='k',ls='dashed')
axes[3].set_ylabel('Emiss (\u03bcg m$^{-2}$ yr$^{-1}$)')
axes[3].set_xlabel('$a$ (prefactor)')
axes[3].set_title('Northern Hemisphere extratropics (grassland)')

axes[0].legend(['Model','Observational constraint'])
f.savefig('Figures/SI_param_tune_soil_emiss.pdf',bbox_inches = 'tight')

