#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on May 6 2023
Make plots of masks for DFR runs
@author: arifeinberg
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import time
#%% Load necessary variables for Olson land map, LAI
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
#%% Load masks for biomes and country
realm_names = ['Nearctic', 'Palearctic','Afrotropic','Indomalayan','Australasia+Oceania','Neotropic', 'China']
mask_all = [0]*len(realm_names) # save masks here
# ISSUE WITH MASKS ORDER OF LATITUDE
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
mask_Am_neo = np.zeros(np.shape(mask_all_np[5,:,:]))

for ii in range(lat_ind_min, lat_ind_max+1):
    for jj in range(lon_ind_min, lon_ind_max+1):
        if is_rainforest[ii,jj] == True:
            mask_other_neo[ii,jj] = 0 # remove rainforest points
            mask_Am_neo[ii,jj] = 1 # add rainforest points
#%% Create plots
mask_all_plot = np.zeros((720, 1440))
mask_names = ['Nearctic', 'Palearctic','Afrotropic','Indomalayan','Australasia','Neotropic-Amazon', 'Neotropic-other', 'China']
mask_all_plot = np.where(mask_all_np[0,:,:]==1.,1, mask_all_plot)  # Nearctic
mask_all_plot = np.where(mask_all_np[1,:,:]==1.,2, mask_all_plot)  # Palearctic
mask_all_plot = np.where(mask_all_np[2,:,:]==1.,3, mask_all_plot)  # Afrotropic
mask_all_plot = np.where(mask_all_np[3,:,:]==1.,5, mask_all_plot)  # Indomalayan
mask_all_plot = np.where(mask_all_np[4,:,:]==1.,6, mask_all_plot)  # Australasia
mask_all_plot = np.where(mask_Am_neo==1.,4, mask_all_plot)  # Neotropic Amazon
mask_all_plot = np.where(mask_other_neo==1.,7, mask_all_plot)  # Neotropic other
mask_all_plot = np.where(mask_all_np[6,:,:]==1.,8, mask_all_plot)  # China
from matplotlib.colors import ListedColormap
colors = [(1, 1, 1), (0.2, 0.133, 0.5725), (0.53333333, 0.8, 0.93333333), \
          (0.26666667, 0.66666667, 0.6), (0.06666667, 0.46666667, 0.2),\
          (0.6, 0.6, 0.2), (0.86666667, 0.8, 0.46666667), \
          (0.8, 0.4, 0.46666667), (0.53333333, 0.13333333, 0.33333333)]
cm_colorblind = ListedColormap(colors, name='clb')
f, axes = plt.subplots(1, 1, figsize=[12,8],subplot_kw=dict(projection=ccrs.PlateCarree()),
                              gridspec_kw=dict(hspace=0.3, wspace=0.2))
h = axes.pcolormesh(lon,lat,mask_all_plot,vmin=-0.5, vmax=8.5,
                   rasterized = True,cmap = cm_colorblind)# cmap='Pastel1_r')
#cbar = plt.colorbar(h, ax=axes, fraction=0.046, pad=0.04,
#                orientation='horizontal')
#axes.set_title(mask_names[i],fontsize=14)
axes.coastlines(linewidth=1)
axes.set_ylim([-60, 85])
f.savefig('Figures/realm_masks_DFR_v2.pdf',bbox_inches = 'tight')


