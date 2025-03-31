#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu April 14 2022
Calculating uncertainty in biomass burning emissions for different scenarios, emission factors
@author: arifeinberg
"""
import xarray as xr
import numpy as np
import netCDF4
from datetime import date
import matplotlib.pyplot as plt
from calendar import monthrange

#%% Setup variables
# Choose scenario
sce = 'BAU' # scenario: 'HIST','BAU', or 'GOV'
sce_l = sce.lower() # lower case scenario

# load parameters for uncertainty analysis
params_unc = np.loadtxt("misc_Data/LHS_sampling_dep_emis.csv",
                 delimiter=",")
num = len(params_unc) # number of uncertainty calculations
bb_emiss_fac = params_unc[:,5] # biomass burning emis factors ug/m2/yr

bb_emiss_ref = 380 # reference biomass burning emis factor ug/m2/yr

# Scaling factors for emissions
sav_scale = 7.23E-08 # SAV - savannah burning
borf_scale = 1.50E-07 # BORF - boreal forest burning
temp_scale = 1.50E-07 # TEMP - temperate forest burning
defo_scale = 5.85E-08 # DEFO - deforestation burning
peat_scale = 7.56E-08 # PEAT - peatland burning
agri_scale = 4.48E-08 # AGRI - agricultural burning

s_yr = 365.2425*24*60*60 # yr to s

kg_ton = 1000 # kg per Mg

# Load Amazon mask
fn_Am_mask = '../masks/Amazon_basin_mask_025x025.nc'
ds_Am_mask = xr.open_dataset(fn_Am_mask)
Am_mask = ds_Am_mask.MASK.values #unitless

# Get the array indices that need to calculate for Amazon
Am_calc_bool = np.where(Am_mask>0, True, False)

#%% Load necessary variables for Olson land map, LAI
# load gbox areas for South America
fn_gb = '../gbox_areas/gbox_025x025_Olson.nc'
ds_gb = xr.open_dataset(fn_gb)
lon = ds_gb.lon # longitude
lat = ds_gb.lat # latitude

gbox_area = ds_gb.cell_area.values # grid_area
gbox_area_A = gbox_area[Am_calc_bool]

#%% list file names
if sce == 'HIST':
    fn_BB = 'misc_Data/BAU/BB_slow/2002/GFED4_gen.025x025.BAU_2002'
    yr = 2002
elif sce == 'BAU':
    fn_BB = 'misc_Data/BAU/BB_slow/2050/GFED4_gen.025x025.BAU_2050'
    yr = 2050    
elif sce == 'GOV':
    fn_BB = 'misc_Data/GOV/BB_slow/2050/GFED4_gen.025x025.GOV_2050'
    yr = 2050    
else:
    print('error, scenario name is incorrect')

#%% calculate the amount of emissions in the reference EF value scenario
c_mt = 0 # counter

BB_tot_AMZN = np.zeros(12) # total BB in Amazon
BB_DEFO_AMZN = np.zeros(12) # DEFO BB in Amazon

for k in range(1, 13): # looping over months
    mt = f'{k:02d}' # add leading zeros
    num_days = monthrange(yr, k)[1] # num_days = 28
    s_in_mt = num_days * 24 * 60 * 60

    # input file name
    fn_BB_tot = fn_BB + mt + '.nc'
        
    # load emissions
    ds_BB = xr.open_dataset(fn_BB_tot)
    DM_sav = ds_BB.DM_SAVA.values.squeeze() # savannah burning
    DM_borf = ds_BB.DM_BORF.values.squeeze() # boreal forest burning
    DM_temp = ds_BB.DM_TEMP.values.squeeze() # temperate forest burning
    DM_defo = ds_BB.DM_DEFO.values.squeeze() # deforestation burning
    DM_peat = ds_BB.DM_PEAT.values.squeeze() # peatland burning
    DM_agri = ds_BB.DM_AGRI.values.squeeze() # agricultural burning
    
    # calculate scaled total Hg in kg/m^2/s
    Hg_total = sav_scale * DM_sav + borf_scale * DM_borf + \
            temp_scale * DM_temp + defo_scale * DM_defo + \
                peat_scale * DM_peat + agri_scale * DM_agri
    Hg_DEFO = defo_scale * DM_defo
    
    Hg_total_A = Hg_total[Am_calc_bool]
    Hg_DEFO_A = Hg_DEFO[Am_calc_bool]

    BB_tot_AMZN[c_mt] = np.sum(Hg_total_A * gbox_area_A) / kg_ton * s_in_mt # ton/mt
    BB_DEFO_AMZN[c_mt] = np.sum(Hg_DEFO_A * gbox_area_A) / kg_ton * s_in_mt # ton/mt
            
    # increment counter
    c_mt = c_mt + 1

# calculate annual total
BB_tot_AMZN_y = np.sum(BB_tot_AMZN)
BB_DEFO_AMZN_y = np.sum(BB_DEFO_AMZN)

# calculate BB without deforestation
BB_oth_AMZN_y = BB_tot_AMZN_y - BB_DEFO_AMZN_y

#%% Loop through parameter values, calculate BB emissions
Amazon_BB = np.zeros(num)  # for listing Amazon totals

# looping through param options
for i in range(num):
    scale = bb_emiss_fac[i] / bb_emiss_ref # scaling factor for deforestation emission
    BB_DEFO_i = BB_DEFO_AMZN_y * scale # scale deforestation emissions
    Amazon_BB[i] = BB_DEFO_i + BB_oth_AMZN_y # sum with other to get total emissions

# save data to csv, for access
ofn_bb = "misc_Data/unc_Amazon_bb_" + sce + ".csv"
np.savetxt(ofn_bb, Amazon_BB, delimiter=",")
