#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 2022
Creating Olson land types and XLAI files for savanna land types, with 2003 as base year
@author: arifeinberg
"""
import os
import xarray as xr
import numpy as np
import netCDF4

#%% Load necessary variables 
# load gbox areas for South America
fn_gb = 'misc_Data/gbox_025x025_Olson.nc'
ds_gb = xr.open_dataset(fn_gb)
gbox_area = ds_gb.cell_area # grid_area

# required Olson land variables
fn_ols1 = 'misc_Data/Olson_2001_Drydep_Inputs.nc'
ds_ols1 = xr.open_dataset(fn_ols1)

DRYCOEFF = ds_ols1.DRYCOEFF.values # Baldocchi dry deposition polynomial coefficients
IOLSON = ds_ols1.IOLSON.values # Olson land type indices (+1)
IDRYDEP = ds_ols1.IDRYDEP.values # Dry deposition land types

IDEP = ds_ols1.IDEP.values # Mapping index: Olson land type ID to drydep ID

# Land type areas
fn_ols2 = 'misc_Data/Olson_2001_Land_Type_Masks.025x025.generic.nc'
ds_ols2 = xr.open_dataset(fn_ols2)
# save as one np array, first dim is land types
Olson_landtype = ds_ols2.to_array().squeeze()

lon = ds_ols2.lon # longitude
lat = ds_ols2.lat # latitude

n_landtypes = 73 # number of Olson land types

# Land type LAI
fn_xlai  = 'misc_Data/Yuan_proc_MODIS_XLAI.025x025.2003.nc'
ds_xlai = xr.open_dataset(fn_xlai)

# save as one data array, first dim is land types
XLAI_w = ds_xlai.to_array()

#%% select Amazon coordinate bounds
lat_min = -19.375
lat_max = 8.375
lon_min = -79.625
lon_max = -41.875

# restrict lat and lon
lat_A = lat[(lat<=lat_max) & (lat>=lat_min)]
lon_A = lon[(lon<=lon_max) & (lon>=lon_min)]

lon_l = len(lon_A)
lat_l = len(lat_A)

Olson_landtype_A = Olson_landtype[:,(lat<=lat_max)&(lat>=lat_min),(lon<=lon_max)&(lon>=lon_min)]
XLAI_w_A = XLAI_w[:,:, (lat<=lat_max)&(lat>=lat_min), (lon<=lon_max)&(lon>=lon_min)]

#%% Change tropical rainforest landclasses to savanna in Amazon
savanna_id = 58 # Olson id: 58 = Fields and Woody Savanna

gbox_area_A = gbox_area[(lat<=lat_max) & (lat>=lat_min),(lon<=lon_max) & (lon>=lon_min)]

# Indices of rainforest
IOLSON_6 = np.asarray(np.where(IDEP==6)).flatten()

# create new matrices for edited land use data
XLAI_w_sav = XLAI_w_A.copy()
time_w = XLAI_w_sav.time
time_l = len(time_w)

Olson_landtype_sav = Olson_landtype_A.copy()


# create LAI threshold to average over
LAI_thres = 0.1 # if below this, then not reliable amount of savanna vegetation

# calculate timeseries of South American savanna average LAI
XLAI_avg = np.zeros(time_l)
for i in range(time_l):
    # Calculate average savanna LAI over Amazon domain
    XLAI_sel = XLAI_w_A[savanna_id, i, :, :]
    
    XLAI_avg[i] =  np.mean(XLAI_sel.where(XLAI_sel>LAI_thres) * gbox_area_A.where(XLAI_sel>LAI_thres)) \
        / np.mean(gbox_area_A.where(XLAI_sel>LAI_thres)) # take average of all mixed values

# find indices where have tropical rainforest
inds_RF = np.asarray(np.where(Olson_landtype_A[IOLSON_6,:,:] > 0)).T
a = 0
# loop over these indices, correct to savanna
for i in range(len(inds_RF)):
    lat_i = inds_RF[i,1]
    lon_i = inds_RF[i,2]
    # previously occupied RF area goes to savanna
    total_RF_area = sum(Olson_landtype_A[IOLSON_6,lat_i,lon_i])
    Olson_landtype_sav[savanna_id,lat_i,lon_i] = Olson_landtype_A[savanna_id,lat_i,lon_i] + total_RF_area
    Olson_landtype_sav[IOLSON_6,lat_i,lon_i] = 0.
    # replace by savanna type vegetation
    if (np.mean(XLAI_w_sav[savanna_id, :, lat_i,lon_i]) <= LAI_thres): # no savanna vegetation already
        XLAI_w_sav[savanna_id, :, lat_i,lon_i] = XLAI_avg
    else:
        a +=1  
    XLAI_w_sav[IOLSON_6,:,lat_i,lon_i] = 0.  

# Replace within global matrix
Olson_landtype[:,(lat<=lat_max)&(lat>=lat_min),(lon<=lon_max)&(lon>=lon_min)] = Olson_landtype_sav

XLAI_w[:,:, (lat<=lat_max)&(lat>=lat_min), (lon<=lon_max)&(lon>=lon_min)] = XLAI_w_sav
#%% Rewrite netcdf
fname = 'misc_Data/SAV/Olson_2001_Land_Type_Masks.025x025.generic_sav' + str(savanna_id) + '.nc'
fname2 = 'misc_Data/SAV/Yuan_proc_MODIS_XLAI.025x025.SAV_' + str(savanna_id) + '.nc'

os.system('cp misc_Data/Olson_2001_Land_Type_Masks.025x025.generic.nc ' + fname)
os.system('cp misc_Data/Yuan_proc_MODIS_XLAI.025x025.2003.nc ' + fname2)

dset = netCDF4.Dataset(fname, 'r+')
dset2 = netCDF4.Dataset(fname2, 'r+')

# loop through, replace rainforest land type
for i in IOLSON_6:
    varname = 'LANDTYPE' + str(i) # variable name in netcdf
    dset[varname][:] = Olson_landtype[i,:,:][:]
    varname2 = 'XLAI' + str(i) # variable name in netcdf
    dset2[varname2][:] = XLAI_w[i,:,:,:][:]

# replace savanna land type
varname = 'LANDTYPE' + str(savanna_id) # variable name in netcdf
dset[varname][:] = Olson_landtype[savanna_id,:,:][:]
dset.close()

varname2 = 'XLAI' + str(savanna_id) # variable name in netcdf
dset2[varname2][:] = XLAI_w[savanna_id,:,:,:][:]
dset2.close()



