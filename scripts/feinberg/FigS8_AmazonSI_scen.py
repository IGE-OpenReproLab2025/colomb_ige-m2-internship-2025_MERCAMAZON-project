#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu April 14 2022
Creating maps for scenarios of Soares-Filho06
The .tif data for Soares-Filho06 can be downloaded here: https://doi.org/10.3334/ORNLDAAC/1153
@author: arifeinberg
"""
import os
import xarray as xr
import numpy as np
import netCDF4
import rioxarray
import cartopy.crs as ccrs

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

#%% load file
        
# Load GEOTiff files to calculate deforested area in current year
# past year
fn1 = '../Amazon_scenarios/BAU_amazonia/amazonbau2002.tif'
ds1 = xr.open_dataset(fn1, engine="rasterio")
landtypes_2002 = ds1.band_data.squeeze()

fn1 = '../Amazon_scenarios/BAU_amazonia/amazonbau2050.tif'
ds1 = xr.open_dataset(fn1, engine="rasterio")
landtypes_2050 = ds1.band_data.squeeze()

fn1 = '../Amazon_scenarios/GOV_amazonia/amazongov2050.tif'
ds1 = xr.open_dataset(fn1, engine="rasterio")
landtypes_GOV = ds1.band_data.squeeze()

#%% plot
# create colormap
colors = ['#d95f02', '#1b9e77', '#e6ab02']  # R -> G -> B
cmap = LinearSegmentedColormap.from_list('my_list', colors, N=3)

f, axes = plt.subplots(3, 1, figsize=[12,14],subplot_kw=dict(projection=ccrs.PlateCarree()),
                              gridspec_kw=dict(hspace=0.2, wspace=0.1))
axes = axes.flatten()
landtypes_2002.plot.pcolormesh(cmap=cmap, ax=axes[0],
                               add_colorbar=False, rasterized = True)
axes[0].set_title('HIST', fontsize=22)
axes[0].text(0.0, 1.04, 'a', fontweight='bold',fontsize=17,
          horizontalalignment='left', verticalalignment='center',
          transform=axes[0].transAxes)

landtypes_2050.plot.pcolormesh(cmap=cmap, ax=axes[1],
                               add_colorbar=False, rasterized = True)
axes[1].set_title('BAU', fontsize=22)
axes[1].text(0.0, 1.04, 'b', fontweight='bold',fontsize=17,
          horizontalalignment='left', verticalalignment='center',
          transform=axes[1].transAxes)

landtypes_GOV.plot.pcolormesh(cmap=cmap, ax=axes[2],
                               add_colorbar=False, rasterized = True)
axes[2].set_title('GOV', fontsize=22)
axes[2].text(0.0, 1.04, 'c', fontweight='bold',fontsize=17,
          horizontalalignment='left', verticalalignment='center',
          transform=axes[2].transAxes)

f.savefig('Figures/SI_scen_S1.pdf',bbox_inches = 'tight')
