#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 16:17:44 2022
Plot differences in soil emissions
GC run0105 results (old soil emissions parametrization) can be downloaded here: https://doi.org/10.7910/DVN/R7NRNK
@author: arifeinberg
"""

import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import xarray as xr
import sys
sys.path.append('pythonHgBenchmark')
from helper_functions import ds_sel_yr, annual_avg, open_Hg
import matplotlib.colors as colors

#%% load area and masks
# Load grid cell area for unit conversion of model
fn_gbox = '../gbox_areas/GEOSChem_2x25_gboxarea.nc'
ds_gbox = xr.open_dataset(fn_gbox)
gbox_GC = ds_gbox.cell_area #m2

#%% Load emissions fluxes
sim = ['0311','0105']
sim_name = ['newsoil', 'oldsoil']

# initialize matrices
soil_emiss_a = []

# run loop
for i in range(len(sim)):
    sim_i = sim[i]
    fn_emis = '../GEOS-Chem_runs/run'+sim_i+'/OutputDir/GEOSChem.MercuryEmis.2015_m.nc4'
    ds_emis = xr.open_dataset(fn_emis)
    lon = ds_emis.lon
    lat = ds_emis.lat
    soil_emis_yr = ds_sel_yr(ds_emis, 'EmisHg0soil', 2015) 
    
    
    # Convert model data from kg/s to ug/m2/yr for annual average   
    s_in_yr = 365.2425 * 24 * 3600 # s in one year
    ug_kg = 1e9 # ug in kg
    
    unit_conv = s_in_yr / gbox_GC * ug_kg
        
    soil_emis = annual_avg(soil_emis_yr) * unit_conv # kg/yr
    
    soil_emiss_a.append(soil_emis) # kg/yr
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
#%% Plot map plots
lat_min = -60
lat_max = 85
Change_soil = soil_emiss_a[0]-soil_emiss_a[1]
# pick the desired colormap, sensible levels, and define a normalization
# instance which takes data values and translates those into levels.
f, axes = plt.subplots(2, 1, figsize=[16,12],subplot_kw=dict(projection=ccrs.PlateCarree()),
                              gridspec_kw=dict(hspace=0.1, wspace=0.1))
axes = axes.flatten()
ax1 = axes[0]
# Have lowest color as white
upper = plt.cm.Oranges(np.arange(256))
lower = np.ones((1,4))
cmap_new = np.vstack((lower, upper))
# convert to matplotlib colormap
cmap_new = colors.ListedColormap(cmap_new, name='Oranges_n', N=cmap_new.shape[0])

h = ax1.pcolormesh(lon_pcolormesh_midpoints, lat_pcolormesh_midpoints,
                   soil_emiss_a[0], cmap =cmap_new,
                   vmin=0, vmax=35,
                   rasterized = True)
h.cmap.set_under('w')
#cbar.tick_params(labelsize=14)
ax1.set_ylim([lat_min, lat_max])
# ax1.set_title('Impact of reforestation on surface-atmosphere exchange', fontweight='bold', fontsize=17)
ax1.text(0.0, 1.04, 'a', fontweight='bold',fontsize=17,
             horizontalalignment='left', verticalalignment='center',
             transform=ax1.transAxes)

ax1.coastlines(linewidth=1.5)

cbar = plt.colorbar(h, extend='max',orientation='vertical',
                    ax=ax1, aspect=30, pad=0.03)
cbar.set_label('Soil emissions (\u03bcg m$^{-2}$ yr$^{-1}$)', fontsize = 20)
cbar.ax.tick_params(labelsize=18) 
ax2 = axes[1]
h = ax2.pcolormesh(lon_pcolormesh_midpoints, lat_pcolormesh_midpoints,
                   Change_soil, cmap ='bwr',
                   vmin=-20, vmax=20,
                   rasterized = True)
#cbar.tick_params(labelsize=14)
ax2.set_ylim([lat_min, lat_max])
# ax1.set_title('Impact of reforestation on surface-atmosphere exchange', fontweight='bold', fontsize=17)
ax2.text(0.0, 1.04, 'b', fontweight='bold',fontsize=17,
             horizontalalignment='left', verticalalignment='center',
             transform=ax2.transAxes)

ax2.coastlines(linewidth=1.5)

cbar = plt.colorbar(h, extend='both',orientation='vertical',
                    ax=ax2, aspect=30, pad=0.03)
cbar.set_label('Difference (\u03bcg m$^{-2}$ yr$^{-1}$)', fontsize = 20)
cbar.ax.tick_params(labelsize=18) 

f.savefig('Figures/Soil_emissions_map.pdf',bbox_inches = 'tight')
