#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 16:17:44 2022
Plot LAI maps needed for SI of paper - reforestation
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

# RFR LAI 
fn_RFR_lai = '../input_data/RFR/Yuan_proc_MODIS_XLAI.025x025.RFR_Griscom_2015.nc'
ds_RFR_lai = xr.open_dataset(fn_RFR_lai)

# save as one data array, first dim is land types
XLAI_w = ds_RFR_lai.to_array()
LAI_RFR = XLAI_w.mean("time").sum(axis=0) # take annual mean

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

#%% select global coordinate bounds
lat_min = -60
lat_max = 85

#%% Plot data
f, axes = plt.subplots(3, 1, figsize=[12,14],subplot_kw=dict(projection=ccrs.PlateCarree()),
                              gridspec_kw=dict(hspace=0.2, wspace=0.1))
axes = axes.flatten()
LAI_to_plot = [LAI_2003, LAI_RFR]
titles = ['HIST','RFR']
# Create colormap
# Have lowest color as white
upper = plt.cm.Greens(np.arange(256))
lower = np.ones((1,4))
cmap_new = np.vstack((lower, upper))
# convert to matplotlib colormap
cmap_new = colors.ListedColormap(cmap_new, name='Greens_n', N=cmap_new.shape[0])
letter =['a','b']
for i in range(2):
    h = axes[i].pcolormesh(lon_pcolormesh_midpoints, lat_pcolormesh_midpoints,
                   LAI_to_plot[i], cmap =cmap_new,
                   vmin=0, vmax=4,
                   rasterized = True)
    axes[i].set_title(titles[i], fontsize=22)
    axes[i].set_ylim([lat_min, lat_max])
    axes[i].coastlines(linewidth=1)

    cbar = plt.colorbar(h, extend='max',orientation='vertical',
                   ax=axes[i],pad=0.03)
    cbar.set_label('LAI (m$^{2}$ m$^{-2}$)', fontsize=20)   
    cbar.ax.tick_params(labelsize=18)
    axes[i].text(0.0, 1.04, letter[i], fontweight='bold',fontsize=17,
             horizontalalignment='left', verticalalignment='center',
             transform=axes[i].transAxes)

upper = plt.cm.Reds(np.arange(256))
lower = np.ones((1,4))
cmap_new2 = np.vstack((lower, upper))
# convert to matplotlib colormap
cmap_new2 = colors.ListedColormap(cmap_new2, name='Reds_n', N=cmap_new2.shape[0])

h = axes[2].pcolormesh(lon_pcolormesh_midpoints, lat_pcolormesh_midpoints,
                LAI_RFR-LAI_2003, cmap =cmap_new2,
                vmin=0, vmax=3,
                rasterized = True)
axes[2].set_title('RFR – HIST', fontsize=22)
axes[2].set_ylim([lat_min, lat_max])
axes[2].coastlines(linewidth=1)
axes[2].text(0.0, 1.04, 'c', fontweight='bold',fontsize=17,
         horizontalalignment='left', verticalalignment='center',
         transform=axes[2].transAxes)

cbar = plt.colorbar(h, extend='max',orientation='vertical',
               ax=axes[2],pad=0.03)
cbar.set_label('Difference (m$^{2}$ m$^{-2}$)', fontsize=20)   
cbar.ax.tick_params(labelsize=18)

f.savefig('Figures/SI_LAI_RFR.pdf',bbox_inches = 'tight')
