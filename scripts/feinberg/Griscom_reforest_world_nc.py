#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 19 2022
Creating Olson land types and XLAI files for Griscom reforest scenario
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

#%% Load biomes raster file
fn_bio = "misc_Data/biome_ecozones_025x025.nc"
ds_bio = xr.open_dataset(fn_bio)
# # %% calculate how many overlaps have in ecozones
# sum_a = xr.zeros_like(ds_bio.Neotropic_14)
# for ii in ds_bio.data_vars:
#     sum_a = sum_a + ds_bio[ii]
# sum_a.plot.pcolormesh()
#%% For each of the realms/regions, calculate the average LAI over forested areas
bio_to_dd = [6, 2, 3, 2, 3, 3, 2, 2, 2, 2, np.nan, 2, 2, 9] # the dry dep category for each biome number
realm_names = ['Nearctic', 'Palearctic','Afrotropic','Indomalayan','Australasia',
               'Oceania','Neotropic']
XLAI_dict = {} # dictionary for XLAI values
# create LAI threshold to average over
LAI_thres = 0.1 # if below this, then not reliable amount of vegetation

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
            gbox_area_mask = gbox_area.where(mask_bio)
            # calculate timeseries of average LAI
            XLAI_avg = np.zeros(time_l)
            # calculate XLAI average in region/biome everywhere as loop:
            for k in range(time_l):
                # Calculate average rainforest LAI over Amazon domain
                XLAI_sel = XLAI_mask[k, :, :]
                # take weighted-area average of all values
                XLAI_avg[k] = np.mean(XLAI_sel.where(XLAI_sel>LAI_thres) * gbox_area_mask.where(XLAI_sel>LAI_thres)) \
                    / np.mean(gbox_area_mask.where(XLAI_sel>LAI_thres)) 
            # save in dictionary
            XLAI_dict[varname] = XLAI_avg

#%% Create plots
f, axes = plt.subplots(2, 3, figsize=[15,8],
                              gridspec_kw=dict(hspace=0.4, wspace=0.3))
axes = axes.flatten()
realm_names_plot = ['Nearctic', 'Palearctic','Afrotropic','Indomalayan','Australasia',
               'Neotropic']
for j, realm in enumerate(realm_names_plot):
    for i in range(1, 15):
        if (not np.isnan(bio_to_dd[i-1])): #don't calculate LAI for non-forest biomes
            varname = realm + "_" + str(i)
            if (not varname in ds_bio.keys()): # combination doesn't exist
                continue # skip iteration
            XLAI_realm = XLAI_dict[varname]
            axes[j].plot(time_w, XLAI_realm,linewidth=2, label=varname)
    axes[j].set_ylabel('Leaf Area Index (m$^2$ / m$^2$)'  )
    axes[j].set_title(realm, fontsize=13)
    axes[j].legend(loc='best')
    
    # fix ticks
    t = axes[j].get_xticks()
    t_new = t[0:7:2]
    axes[j].set_xticks(t_new)

f.savefig('Figures/LAI_averages_realms.pdf',bbox_inches = 'tight')
#%% Check how often reforestation is in a forest grid box
# input filenames
ifn_LT = fn_ols2
ifn_LAI = fn_xlai

# Output file names - LAI and Olson land types
ofn_LT = 'misc_Data/RFR/Olson_Land_Type_Masks.025x025.RFR_Griscom.nc'
ofn_LAI = 'misc_Data/RFR/Yuan_proc_MODIS_XLAI.025x025.RFR_Griscom.nc'

# find indices where have deforestation
inds_RF = np.asarray(np.where((rf>0. ))).T 
n_rf = len(inds_RF)
# create new matrices for edited land use data
XLAI_w_ref = XLAI_w.copy()
Olson_landtype_ref = Olson_landtype.copy()
XLAI_w_avg = XLAI_w.mean("time").sum("variable")

# initialize matrices 
old_lt = np.zeros(n_rf)
new_lt = np.zeros(n_rf)
lai_old = np.zeros(n_rf)
lai_new = np.zeros(n_rf)
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

# time loop
now = time.time()
elapsed = now - start                 
print('elapsed: ' + str(elapsed))

#%% Do the reforestation of LAI
count_no_LAI = 0
# indices of land types to reforest
reforest_indices = [999, 999,  26, 27, 999, 999, 33, # note extra index for python starting at 0
                    999, 999, 72, 999, 999] 
for rr in range(n_rf): # loop over all gridcells
    lat_i = inds_RF[rr,0]
    lon_i = inds_RF[rr,1]
    # realm used to replace, found in last loop
    realm_u = realm_list[rr]
    # LAI to reforest is the average for that realm
    XLAI_replace = XLAI_dict[realm_u]
    # check if the replaced LAI is nan
    if (np.isnan(XLAI_replace.sum())):
        count_no_LAI = count_no_LAI + 1
        lai_new[rr] = lai_old[rr]
        print(realm_u)
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
          print("don't reforest with lower LAI")
          print(lat[lat_i])
          print(lon[lon_i])
          print(realm_u)
          print(new_lt[rr])
          print(old_lt[rr])
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
    
# check how LAI changes in Jan and July. make sure that don't have decreases, looks normal

#%% Reforest areas which were removed pre- 2002

# make check 
print("Min Coverage")
print(np.sum(Olson_landtype_ref, axis=0).min())
print("Max Coverage")
print(np.sum(Olson_landtype_ref, axis=0).max())
print("Max average LAI")
print(XLAI_w_ref.sum("variable").mean("time").max())


os.system('cp ' + ifn_LT + ' ' + ofn_LT) # copy files to new file names
os.system('cp ' + ifn_LAI + ' ' + ofn_LAI) # copy files to new file names

dset = netCDF4.Dataset(ofn_LT, 'r+')
dset2 = netCDF4.Dataset(ofn_LAI, 'r+')

# loop through, replace  land type
for i in range(72):
    varname = 'LANDTYPE' + str(i).zfill(2) # variable name in netcdf
    dset[varname][:] = Olson_landtype_ref[i,:,:][:]
    varname2 = 'XLAI' + str(i).zfill(2) # variable name in netcdf
    dset2[varname2][:] = XLAI_w_ref[i,:,:,:][:]


#%% Compare Olson land map with 2002 reference map from Soares
# Indices of rainforest
f, axes = plt.subplots(3, 3, figsize=[15,8],subplot_kw=dict(projection=ccrs.PlateCarree()),
                              gridspec_kw=dict(hspace=0.4, wspace=0.2))
axes = axes.flatten()
XLAI_new = XLAI_w_ref.sum("variable").mean("time")
XLAI_old = XLAI_w.sum("variable").mean("time")
diff_plot = XLAI_new - XLAI_old

h = axes[0].pcolormesh(lon,lat,XLAI_old,vmin=0, vmax=7,
                       rasterized = True)
cbar = plt.colorbar(h, ax=axes[0], fraction=0.046, pad=0.04,extend='max',
                    orientation='horizontal')
cbar.set_label('LAI')   
axes[0].set_title('Original LAI average',fontsize=16)
axes[0].coastlines(linewidth=2)

h = axes[1].pcolormesh(lon,lat,XLAI_new,vmin=0, vmax=7,
                       rasterized = True)
cbar = plt.colorbar(h, ax=axes[1], fraction=0.046, pad=0.04,extend='max',
                    orientation='horizontal')
cbar.set_label('LAI')   
axes[1].set_title('Griscom LAI average',fontsize=16)
axes[1].coastlines(linewidth=2)

h = axes[2].pcolormesh(lon,lat,diff_plot,cmap='bwr',vmin=-4, vmax=4,
                       rasterized = True)
cbar = plt.colorbar(h, ax=axes[2], fraction=0.046, pad=0.04,extend='both',
                    orientation='horizontal')
cbar.set_label('LAI')
axes[2].set_title('Griscom LAI average',fontsize=16)
axes[2].coastlines(linewidth=2)

XLAI_new = XLAI_w_ref.sum("variable").isel(time=2) # January
XLAI_old = XLAI_w.sum("variable").isel(time=2) # January
diff_plot = XLAI_new - XLAI_old

h = axes[3].pcolormesh(lon,lat,XLAI_old,vmin=0, vmax=7,
                       rasterized = True)
cbar = plt.colorbar(h, ax=axes[3], fraction=0.046, pad=0.04,extend='max',
                    orientation='horizontal')
cbar.set_label('LAI')   
axes[3].set_title('Original LAI Jan',fontsize=16)
axes[3].coastlines(linewidth=2)

h = axes[4].pcolormesh(lon,lat,XLAI_new,vmin=0, vmax=7,
                       rasterized = True)
cbar = plt.colorbar(h, ax=axes[4], fraction=0.046, pad=0.04,extend='max',
                    orientation='horizontal')
cbar.set_label('LAI')   
axes[4].set_title('Griscom LAI Jan',fontsize=16)
axes[4].coastlines(linewidth=2)

h = axes[5].pcolormesh(lon,lat,diff_plot,cmap='bwr',vmin=-4, vmax=4,
                       rasterized = True)
cbar = plt.colorbar(h, ax=axes[5], fraction=0.046, pad=0.04,extend='both',
                    orientation='horizontal')
cbar.set_label('LAI')
axes[5].set_title('Griscom LAI Jan',fontsize=16)
axes[5].coastlines(linewidth=2)

XLAI_new = XLAI_w_ref.sum("variable").isel(time=24) # July
XLAI_old = XLAI_w.sum("variable").isel(time=24) # July
diff_plot = XLAI_new - XLAI_old

h = axes[6].pcolormesh(lon,lat,XLAI_old,vmin=0, vmax=7,
                       rasterized = True)
cbar = plt.colorbar(h, ax=axes[6], fraction=0.046, pad=0.04,extend='max',
                    orientation='horizontal')
cbar.set_label('LAI')   
axes[6].set_title('Original LAI July',fontsize=16)
axes[6].coastlines(linewidth=2)

h = axes[7].pcolormesh(lon,lat,XLAI_new,vmin=0, vmax=7,
                       rasterized = True)
cbar = plt.colorbar(h, ax=axes[7], fraction=0.046, pad=0.04,extend='max',
                    orientation='horizontal')
cbar.set_label('LAI')   
axes[7].set_title('Griscom LAI July',fontsize=16)
axes[7].coastlines(linewidth=2)

h = axes[8].pcolormesh(lon,lat,diff_plot,cmap='bwr',vmin=-4, vmax=4,
                       rasterized = True)
cbar = plt.colorbar(h, ax=axes[8], fraction=0.046, pad=0.04,extend='both',
                    orientation='horizontal')
cbar.set_label('LAI')
axes[8].set_title('Griscom LAI July',fontsize=16)
axes[8].coastlines(linewidth=2)
f.savefig('Figures/RFR_LAI_Comparison_Griscom.pdf',bbox_inches = 'tight')
