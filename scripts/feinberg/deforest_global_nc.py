#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on April 3 2023
Creating Olson land types and XLAI files for global deforestation perturbation runs
@author: arifeinberg
"""
import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import time
import netCDF4
#%% Load necessary variables for Olson land map, LAI
# load gbox areas 
fn_gb = 'misc_Data/gbox_025x025_Olson.nc'
ds_gb = xr.open_dataset(fn_gb)
gbox_area = ds_gb.cell_area # grid_area
gbox_area_v = gbox_area.values # grid_area

# required Olson land variables
fn_ols1 = 'misc_Data/Olson_2001_Drydep_Inputs.nc'
ds_ols1 = xr.open_dataset(fn_ols1)

DRYCOEFF = ds_ols1.DRYCOEFF.values # Baldocchi dry deposition polynomial coefficients
IOLSON = ds_ols1.IOLSON.values # Olson land type indices (+1)
IDRYDEP = ds_ols1.IDRYDEP.values # Dry deposition land types

IDEP = ds_ols1.IDEP.values # Mapping index: Olson land type ID to drydep ID

IOLSON_6 = np.asarray(np.where(IDEP==6)).flatten() # Indices of rainforest
IOLSON_not6 = np.asarray(np.where(IDEP!=6)).flatten() # Indices of all other landtypes
IOLSON_not6 = IOLSON_not6[:-1] # don't include last entry in array

# Land type areas
fn_ols2 = 'misc_Data/Olson_2001_Land_Type_Masks.025x025.generic.nc'
ds_ols2 = xr.open_dataset(fn_ols2)
# save as one np array, first dim is land types
Olson_landtype = ds_ols2.to_array().squeeze()

lon = ds_ols2.lon # longitude
lat = ds_ols2.lat # latitude

n_landtypes = 73 # number of Olson land types

# Land type LAI
fn_xlai  = 'misc_Data/Yuan_proc_MODIS_XLAI.025x025.2003.nc' # Use 2003 as reference year
ds_xlai = xr.open_dataset(fn_xlai)

# save as one data array, first dim is land types
XLAI_w = ds_xlai.to_array()
# time variable
time_w = XLAI_w.time
time_l = len(time_w)
# calculate mean lai for old conditions
lai_old = XLAI_w.sum('variable').mean('time').values
#%% Load masks for biomes and country
realm_names = ['Nearctic', 'Palearctic','Afrotropic','Indomalayan','Australasia+Oceania','Neotropic', 'China']
mask_all = [0]*len(realm_names) # save masks here

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

#%% Load deforestation areas dataset
# Loss of forest
fn_lt = 'misc_Data/transitions_t_total_forest.nc'
ds_lt = xr.open_dataset(fn_lt, decode_times=False)
area_df = xr.where(ds_lt.total_area.squeeze()>0,1,0).reindex(lat=list(reversed(ds_lt.lat)))
area_df = area_df.values# maximize area deforested for perturbation experiments

#%% Fix masks for overlap and maximum above 1
mask_all_np = np.zeros((7,720, 1440))
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
   
    mask_all_np[i,:,:] = temp.values
#%% For each of the realms/regions, calculate the average LAI over replacement areas
replace_lt = [0]*7 # potential replacement land types
replace_lt[0] = 35 # Nearctic
replace_lt[1] = 31 # Palearctic, No China
replace_lt[2] = 31 # Afrotropic
replace_lt[3] = 31 # Indomalaya
replace_lt[4] = 31 # Australasia + Oceania
replace_lt[5] = 35 # Neotropic
replace_lt[6] = 31 # China
XLAI_dict = np.zeros((7, time_l)) # array for XLAI values

# create LAI threshold to average over
LAI_thres = 0.1 # if below this, then not reliable amount of vegetation
for j in range(len(replace_lt)):
    # select raster of realm area
    mask_bio = mask_all_np[j,:,:] # load region mask
    gbox_area_mask = gbox_area.where(mask_bio)
    IOLSON_id = replace_lt[j] # land_type
    # select XLAI_w_dd
    XLAI_w_dd = XLAI_w[IOLSON_id,:,:,:]
    # mask the other variables
    XLAI_mask = XLAI_w_dd.where(mask_bio)
    # calculate timeseries of average LAI
    XLAI_avg = np.zeros(time_l)
    # calculate XLAI average in region/biome everywhere as loop:
    for k in range(time_l):
        # Calculate average LAI over domain
        XLAI_sel = XLAI_mask[k, :, :]
        # take weighted-area average of all values
        XLAI_avg[k] = np.mean(XLAI_sel.where(XLAI_sel>LAI_thres) * gbox_area_mask.where(XLAI_sel>LAI_thres))\
            / np.mean(gbox_area_mask.where(XLAI_sel>LAI_thres)) 
    # save in dictionary
    XLAI_dict[j,:] = XLAI_avg
#%% Create Neotropic, Non-Amazon rainforest mask for replacement
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

mask_other_neo = mask_all_np[5,:,:]

for ii in range(lat_ind_min, lat_ind_max+1):
    for jj in range(lon_ind_min, lon_ind_max+1):
        if is_rainforest[ii,jj] == True:
            mask_other_neo[ii,jj] = 0 # remove rainforest points
#%% Check whether deforestation occurs in a forest grid box, calculate new LAI
# input filenames
ifn_LT = fn_ols2
ifn_LAI = fn_xlai
realm_names_short = ['Nearc', 'Palearc','Afrotrop','Indomalay','Austral','Neotrop', 'China']

# loop over all regions
for ll in range(len(realm_names)): # region index
    print(realm_names[ll])
    print(ll)
    
    # # Output file names - LAI and Olson land types
    ofn_LT = 'misc_Data/DFR/Olson_Land_Type_Masks.025x025.DFR_' +realm_names_short[ll] + '.nc'
    ofn_LAI = 'misc_Data/DFR/Yuan_proc_MODIS_XLAI.025x025.DFR_' +realm_names_short[ll] + '.nc'
    # select only deforested landtypes within region
    mask_bio = mask_all_np[ll,:,:] # load region mask
    if ll ==5: # if Neotropic, only do the non-Amazon rainforest area
        mask_bio = mask_other_neo
        
    area_df_ll = area_df * mask_bio # select only areas within biome
    # find indices where have deforestation
    inds_DF = np.asarray(np.where((area_df_ll>0. ))).T 
    n_df = len(inds_DF)
    # create new matrices for edited land use data
    XLAI_w_def = XLAI_w.copy()
    Olson_landtype_def = Olson_landtype.copy()
    
    # initialize matrices 
    if_deforested = np.zeros(len(inds_DF)) # check if actually deforested (original forest area, higher LAI than replace)
    
    for dd in range(n_df): # loop over all gridcells
        lat_i = inds_DF[dd,0]
        lon_i = inds_DF[dd,1]
        # LAI to reforest is the average for that realm
        XLAI_replace = XLAI_dict[ll, :]
        # check whether the current land type is a forest land type:
        if (dd_lt_v[lat_i, lon_i] in [2, 3, 6]): # forest land types
            if (XLAI_replace.mean() < lai_old[lat_i, lon_i]): # deforest in this case
                Olson_landtype_def[:, lat_i, lon_i] = 0. # set all to 0
                Olson_landtype_def[replace_lt[ll], lat_i, lon_i] = 1. # deforested land type occupies whole box
                XLAI_w_def[:, :, lat_i, lon_i] =0 # remove existing vegetation
                XLAI_w_def[replace_lt[ll], :, lat_i, lon_i] = XLAI_replace # replace with deforested land type
                if_deforested[dd] = 1 # Actually deforested
    
    inds_DF_final = inds_DF[np.asarray(np.where((if_deforested==1))).flatten(),:]
    total_deforested_area = np.sum(gbox_area_v[inds_DF_final[:,0],inds_DF_final[:,1]])/ 1e6 # km^2
    print(total_deforested_area)
    print('Check land type max and min')
    print(Olson_landtype_def.sum(axis=0).max())
    print(Olson_landtype_def.sum(axis=0).min())
    
    os.system('cp ' + ifn_LT + ' ' + ofn_LT) # copy files to new file names
    os.system('cp ' + ifn_LAI + ' ' + ofn_LAI) # copy files to new file names
    
    dset = netCDF4.Dataset(ofn_LT, 'r+')
    dset2 = netCDF4.Dataset(ofn_LAI, 'r+')
    
    # loop through, replace  land type
    for i in range(72):
        varname = 'LANDTYPE' + str(i).zfill(2) # variable name in netcdf
        dset[varname][:] = Olson_landtype_def[i,:,:][:]
        varname2 = 'XLAI' + str(i).zfill(2) # variable name in netcdf
        dset2[varname2][:] = XLAI_w_def[i,:,:,:][:]
        
    dset.close()
    dset2.close()
