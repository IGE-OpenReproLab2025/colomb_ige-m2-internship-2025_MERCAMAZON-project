#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 16:17:44 2022
Plot LAI maps needed for SI of paper- Amazon
@author: arifeinberg
"""

import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import scipy.io as sio
import datetime
import xarray as xr
import matplotlib.colors as colors


#%% load data file
# 2003 LAI 
fn_HIST_lai = '../input_data/HIST/Yuan_proc_MODIS_XLAI.025x025.BAU_2002_2015.nc'
ds_HIST_lai = xr.open_dataset(fn_HIST_lai)

lon = ds_HIST_lai.lon # longitude
lat = ds_HIST_lai.lat # latitude

# save as one data array, first dim is land types
XLAI_w = ds_HIST_lai.to_array()
LAI_2003 = XLAI_w.mean("time").sum(axis=0) # take annual mean

# BAU LAI 
fn_BAU_lai = '../input_data/BAU/Yuan_proc_MODIS_XLAI.025x025.BAU_2050_2015.nc'
ds_BAU_lai = xr.open_dataset(fn_BAU_lai)

# save as one data array, first dim is land types
XLAI_w = ds_BAU_lai.to_array()
LAI_BAU = XLAI_w.mean("time").sum(axis=0) # take annual mean

# GOV LAI 
fn_GOV_lai = '../input_data/GOV/Yuan_proc_MODIS_XLAI.025x025.GOV_2050_2015.nc'
ds_GOV_lai = xr.open_dataset(fn_GOV_lai)

lon = ds_GOV_lai.lon # longitude
lat = ds_GOV_lai.lat # latitude

# save as one data array, first dim is land types
XLAI_w = ds_GOV_lai.to_array()
LAI_GOV = XLAI_w.mean("time").sum(axis=0) # take annual mean

# SAV LAI 
fn_SAV_lai = '../input_data/SAV/Yuan_proc_MODIS_XLAI.025x025.SAV_58_2015.nc'
ds_SAV_lai = xr.open_dataset(fn_SAV_lai)

# save as one data array, first dim is land types
XLAI_w = ds_SAV_lai.to_array()
LAI_SAV = XLAI_w.mean("time").sum(axis=0) # take annual mean

#%% Shift pcolormesh lat and lon
# extend longitude by 2
lon_extend = np.zeros(lon.size+2)
# fill in internal values
lon_extend[1:-1] = lon # fill up with original values
# fill in extra endpoints
lon_extend[0] = lon[0]-np.diff(lon)[0]
lon_extend[-1] = lon[-1]+np.diff(lon)[-1]
# calculate the midpoints
lon_pcolormesh_midpoints = lon_extend[:-1]+0.5*(np.diff(lon_extend))

# extend latitude by 2
lat_extend = np.zeros(lat.size+2)
# fill in internal values
lat_extend[1:-1] = lat
# fill in extra endpoints
lat_extend[0] = lat[0]-np.diff(lat)[0]
lat_extend[-1] = lat[-1]+np.diff(lat)[-1]
# calculate the midpoints
lat_pcolormesh_midpoints = lat_extend[:-1]+0.5*(np.diff(lat_extend))

#%% select Amazon coordinate bounds
lat_min = -19.375
lat_max = 13.375
lon_min = -82.625
lon_max = -31.875

#%% Plot data
f, axes = plt.subplots(2, 2, figsize=[16,12],subplot_kw=dict(projection=ccrs.PlateCarree()),
                              gridspec_kw=dict(hspace=0.2, wspace=0.1))
axes = axes.flatten()
LAI_to_plot = [LAI_2003, LAI_BAU, LAI_GOV, LAI_SAV]
titles = ['HIST','BAU','GOV','SAV']
# Create colormap
# Have lowest color as white
upper = plt.cm.Greens(np.arange(256))
lower = np.ones((1,4))
cmap_new = np.vstack((lower, upper))
# convert to matplotlib colormap
cmap_new = colors.ListedColormap(cmap_new, name='Greens_n', N=cmap_new.shape[0])
for i in range(4):
    h = axes[i].pcolormesh(lon_pcolormesh_midpoints, lat_pcolormesh_midpoints,
                   LAI_to_plot[i], cmap =cmap_new,
                   vmin=0, vmax=7,
                   rasterized = True)
    axes[i].set_title(titles[i], fontsize=22)
    axes[i].set_xlim([lon_min, lon_max])
    axes[i].set_ylim([lat_min, lat_max])
    axes[i].coastlines(linewidth=1.5)

cbar = f.colorbar(h, extend='max',orientation='horizontal',
                  ax=axes[:].ravel().tolist(), pad=0.03,aspect=50)
cbar.set_label('LAI (m$^{2}$ m$^{-2}$)', fontsize=20)   
cbar.ax.tick_params(labelsize=18)

f.savefig('Figures/SI_LAI_Amazon.pdf',bbox_inches = 'tight')
