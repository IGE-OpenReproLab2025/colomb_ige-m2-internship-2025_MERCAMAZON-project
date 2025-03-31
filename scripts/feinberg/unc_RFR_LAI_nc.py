#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 2 2023
Creating Olson land types and XLAI files for Griscom reforest scenario
uncertainty assessment - create 100 nc files, one for each uncertainty scenario
@author: arifeinberg
"""
import matplotlib
matplotlib.use('Agg')
import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import time
import netCDF4
import xesmf as xe
import warnings
from scipy.io import savemat, loadmat

def most_common(lst): # need this function for finding neighbour
    lst_edit = [i for i in lst if i != 0] # remove zeros from list
    if not lst_edit:
        return 0
    else:
        return max(set(lst_edit), key=lst_edit.count)

#%% Load reforestation potential dataset
# required Olson land variables
fn_lt = 'misc_Data/Olson_2001_lt_gridcells.nc'
ds_lt = xr.open_dataset(fn_lt)

# reforestation potential
fn_rf = 'misc_Data/warped_mask_griscom_025x025.nc'
ds_rf = xr.open_dataset(fn_rf)
lat = ds_rf.lat
lon = ds_rf.lon

# load variables
lt = ds_lt.lt
dd_lt = ds_lt.dd_lt
rf = ds_rf.x_reforest
#%% Load necessary variables for Olson land map, LAI
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
fn_xlai  = 'misc_Data/Yuan_proc_MODIS_XLAI.025x025.2003.nc' # Use 2003 as reference year
ds_xlai = xr.open_dataset(fn_xlai)

# save as one data array, first dim is land types
XLAI_w = ds_xlai.to_array()
# time variable
time_w = XLAI_w.time
time_l = len(time_w)

#%% Load biomes raster file
fn_bio = "misc_Data/biome_ecozones_025x025.nc"
ds_bio = xr.open_dataset(fn_bio)
#%% prepare regridding data from 025x025 to 2x25
# example 2 x 2.5 dataset
fn_ols2x25 = 'misc_Data/Olson_2001_Land_Type_Masks.2_25.generic.nc'
ds_ols2x25 = xr.open_dataset(fn_ols2x25)

# prepare regridding using xESMF
regridder = xe.Regridder(ds_xlai, ds_ols2x25, "conservative")

#%% Check how often reforestation is in a forest grid box
bio_to_dd = [6, 2, 3, 2, 3, 3, 2, 2, 2, 2, np.nan, 2, 2, 9] # the dry dep category for each biome number
realm_names = ['Nearctic', 'Palearctic','Afrotropic','Indomalayan','Australasia',
               'Oceania','Neotropic']
# input filenames
ifn_LT = fn_ols2
ifn_LAI = fn_xlai


# find indices where have deforestation
inds_RF = np.asarray(np.where((rf>0. ))).T 
n_rf = len(inds_RF)

XLAI_w_avg = XLAI_w.mean("time").sum("variable")

# initialize matrices 
old_lt = np.zeros(n_rf)
new_lt = np.zeros(n_rf)
lai_old = np.zeros(n_rf)
rf_val = np.zeros(n_rf)

realm_list = [0] * n_rf

break_again = False # for skipping out of 2 loops
start = time.time() # for timing loop
# loop over these indices, add correct land type
for rr in range(n_rf):
    lat_i = inds_RF[rr,0]
    lon_i = inds_RF[rr,1]
    rf_val[rr] = rf.isel(lat=lat_i, lon=lon_i) # reforestation potential
    old_lt[rr] = dd_lt.isel(lat=lat_i, lon=lon_i) # old dd land type
    lai_old[rr] = XLAI_w_avg.isel(lat=lat_i, lon=lon_i) # old mean LAI
    # check what the land type is in this grid cell
    for j in realm_names:
        for i in range(1, 15):
            if (not np.isnan(bio_to_dd[i-1])): #don't calculate LAI for non-forest biomes
                varname = j + "_" + str(i)
                if (not varname in ds_bio.keys()): # combination doesn't exist
                    continue # skip iteration
                # check if it is within area:
                if (ds_bio[varname].isel(lat=lat_i, lon=lon_i) == 1):
                    break_again = True # for skippping realm loop
                    realm_list[rr] = varname
                    new_lt[rr] = bio_to_dd[i-1]
                    break
        if (break_again == True): # check if already found the new forest category
            break_again = False # reset this variable
            break
    if (new_lt[rr] == 0):
        print("Issue with rr: " + str(rr))
        # in this case, investigate 9 neighbouring grid cells, find most common land type and use this
        new_lt_neigh = np.zeros(9)
        realm_list_neigh  = [0] * 9
        count = 0 # counter variable
        for ll in range(-1,2): # lat iterator
            for kk in range(-1,2): # lon iterator
                break_again = False # for skipping out of 2 loops
                lat_neigh = inds_RF[rr,0] + ll
                lon_neigh = inds_RF[rr,1] + kk
                for j in realm_names: # find the realm of neighbour
                    for i in range(1, 15):
                        if (not np.isnan(bio_to_dd[i-1])): #don't calculate LAI for non-forest biomes
                            varname = j + "_" + str(i)
                            if (not varname in ds_bio.keys()): # combination doesn't exist
                                continue # skip iteration
                            # check if it is within area:
                            if (ds_bio[varname].isel(lat=lat_neigh, lon=lon_neigh) == 1):
                                break_again = True # for skippping realm loop
                                realm_list_neigh[count] = varname
                                new_lt_neigh[count] = bio_to_dd[i-1]
                                break
                    if (break_again == True): # check if already found the new forest category
                        break_again = False # reset this variable
                        break
                count = count + 1 # increase counter
        # figure out which was the most common value
        realm_list[rr] = most_common(realm_list_neigh)
        # find index of most common value
        idx_neigh = realm_list_neigh.index(most_common(realm_list_neigh))
        # set the new land type to this common value
        new_lt[rr] = new_lt_neigh[idx_neigh]
        # account for cases where have no neighbouring grid cells
        if realm_list[rr] == 0: # will just skip these rare cases
            print("No biome found for rr: " + str(rr))
# save temporary data for easiness of debugging and running
savemat("misc_Data/temp_data_RFR.mat", {"old_lt": old_lt, "new_lt": new_lt, 
                                        "lai_old": lai_old, "rf_val": rf_val,
                                        "realm_list": realm_list})

time loop
now = time.time()
elapsed = now - start                 
print('elapsed: ' + str(elapsed))


# load previously calculated data
data_RFR = loadmat("misc_Data/temp_data_RFR.mat")
old_lt = data_RFR['old_lt'].squeeze()
new_lt = data_RFR['new_lt'].squeeze()
lai_old = data_RFR['lai_old'].squeeze()
rf_val = data_RFR['rf_val'].squeeze()
realm_list = data_RFR['realm_list']

#%% load uncertainty parameter values
# columns: 0) soil parametrization 1) LAI_pct_replace 2)- f0 Amazon RF 
#          3) f0 other RF 4) f0 elsewhere 5) biomass burning emissions
params_unc = np.loadtxt("misc_Data/LHS_sampling_dep_emis.csv",
                 delimiter=",")
num = len(params_unc) # number of uncertainty calculations

# loop through, create nc files for each LAI percentile scenario
for ll in range(num):
    print('LL is ' + str(ll))
    # calculate replaced LAI map
    LAI_pct = params_unc[ll,1]
    
    start = time.time() # for timing loop
    
    # For each of the realms/regions, calculate the average LAI over forested areas
    XLAI_dict = {} # dictionary for XLAI values
    # create LAI threshold to average over
    LAI_thres = 0.1 # if below this, then not reliable amount of vegetation
    
    # create new matrices for edited land use data
    XLAI_w_ref = XLAI_w.copy()
    Olson_landtype_ref = Olson_landtype.copy()

    for j in realm_names:
        for i in range(1, 15):
            if (not np.isnan(bio_to_dd[i-1])): #don't calculate LAI for non-forest biomes
                varname = j + "_" + str(i)
                if (not varname in ds_bio.keys()): # combination doesn't exist
                    continue # skip iteration
                dd_id = bio_to_dd[i-1] # dry deposition (dd) land type
                IOLSON_id = np.asarray(np.where(IDEP==dd_id)).flatten() # Olson indices of dd type
                if (dd_id==9): # assume mangrove area, not all wetlands
                    IOLSON_id = 72 # mangrove land type
                    # select XLAI_w_dd
                    XLAI_w_dd = XLAI_w[IOLSON_id,:,:,:]
                else: # all other land types
                    # sum XLAI over the dd type
                    XLAI_w_dd = XLAI_w[IOLSON_id,:,:,:].sum("variable")
                # select raster of realm area
                mask_bio = ds_bio[varname] # load region mask
                # mask the other variables
                XLAI_mask = XLAI_w_dd.where(mask_bio)
                # calculate timeseries of average LAI
                XLAI_avg = np.zeros(time_l)
                # calculate XLAI average in region/biome everywhere as loop:
                with warnings.catch_warnings(): # ignore warnings for empty arrays
                    warnings.filterwarnings(action='ignore', message='All-NaN slice encountered') 
                    for k in range(time_l):
                        # Calculate percentile LAI over domain
                        XLAI_sel = XLAI_mask[k, :, :]
                        XLAI_avg[k] = np.nanpercentile(XLAI_sel.where(XLAI_sel>LAI_thres),LAI_pct)
                        #ISSUE with nans?
                # save in dictionary
                XLAI_dict[varname] = XLAI_avg
    
    
    # Do the reforestation of LAI
    count_no_LAI = 0
    # indices of land types to reforest
    reforest_indices = [999, 999,  26, 27, 999, 999, 33, # note extra index for python starting at 0
                        999, 999, 72, 999, 999] 
    lai_new = np.zeros(n_rf)
    
    for rr in range(n_rf): # loop over all gridcells
        lat_i = inds_RF[rr,0]
        lon_i = inds_RF[rr,1]
        # realm used to replace, found in last loop
        realm_u = realm_list[rr].replace(" ", "")
        # LAI to reforest is the average for that realm
        XLAI_replace = XLAI_dict[realm_u]
        # check if the replaced LAI is nan
        if (np.isnan(XLAI_replace.sum())):
            count_no_LAI = count_no_LAI + 1
            lai_new[rr] = lai_old[rr]
            # print(realm_u)
            # around 150 are missing, or around 0.2%
            continue
        # check whether the current land type is the reforested land type:
        if (old_lt[rr] == new_lt[rr]): # no change in land type
           lt_use = lt.isel(lat=lat_i, lon=lon_i)
           # make forest denser by adding average LAI
           XLAI_w_ref[lt_use, :, lat_i,lon_i] = XLAI_w_ref[lt_use, :, lat_i,lon_i] + rf_val[rr] * XLAI_replace
        else: # change in land type
           # only reforest if it will increase LAI of overall grid cell
           if (XLAI_replace.mean() < lai_old[rr]):
              # print("don't reforest with lower LAI")
              # print(lat[lat_i])
              # print(lon[lon_i])
              # print(realm_u)
              # print(new_lt[rr])
              # print(old_lt[rr])
              continue # don't reforest
           # add new forest area
           ind_O = reforest_indices[int(new_lt[rr])] # Olson land type
           Olson_landtype_ref[ind_O, lat_i, lon_i] = rf_val[rr]
           XLAI_w_ref[ind_O, :, lat_i, lon_i] = XLAI_replace * rf_val[rr] # need to normalize by area
           # old land area needs to be adjusted
           lt_use = lt.isel(lat=lat_i, lon=lon_i)
           Olson_landtype_ref[lt_use, lat_i, lon_i] = 1. - rf_val[rr]
           # old land type LAI doesn't change, but needs to be normalized by new area
           XLAI_w_ref[lt_use, :, lat_i, lon_i] = XLAI_w_ref[lt_use, :, lat_i, lon_i] * (1. - rf_val[rr])
    
    # Final outputs
    
    # make check 
    print("Min Coverage")
    print(np.sum(Olson_landtype_ref, axis=0).min())
    print("Max Coverage")
    print(np.sum(Olson_landtype_ref, axis=0).max())
    
    # Resample to monthly values
    XLAI_w_ref_m = XLAI_w_ref.resample(time="1MS").mean(dim="time")
    
    # regrid using xESMF
    XLAI_w_ref_rm  = regridder(XLAI_w_ref_m, keep_attrs=True)
    Olson_landtype_ref_rm  = regridder(Olson_landtype_ref, keep_attrs=True)
    
    XLAI_w_ref_rm.name = 'XLAI'
    Olson_landtype_ref_rm.name = 'Olson_LT'
    
    # Output file names - LAI and Olson land types
    ofn_LT = 'misc_Data/RFR/unc_LAI/Olson_Land_Type_Masks.2x25.RFR_unc_' + str(ll) + '.nc'
    ofn_LAI = 'misc_Data/RFR/unc_LAI/Yuan_proc_MODIS_XLAI.2x25.RFR_unc_' + str(ll) + '_m.nc'
    
    # save files
    Olson_landtype_ref_rm.to_netcdf(ofn_LT)
    XLAI_w_ref_rm.to_netcdf(ofn_LAI)
    
    # delete variable for safety
    del Olson_landtype_ref, XLAI_w_ref_m
    
# time loop
now = time.time()
elapsed = now - start                 
print('elapsed: ' + str(elapsed))
