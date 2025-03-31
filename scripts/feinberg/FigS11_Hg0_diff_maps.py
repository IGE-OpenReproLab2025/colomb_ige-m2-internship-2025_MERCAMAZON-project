#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 16:17:44 2022
Plot differences in Hg0 concentrations
@author: arifeinberg
"""

import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import xarray as xr
from helper_functions import ds_sel_yr, annual_avg, open_Hg
import matplotlib.colors as colors

#%% Load concentrations
sim = ['0311','0312','0313','0315','0314']
sim_name = ['BAU', 'GOV','SAV','RFR']

# initialize matrices
Hg0_conc_a = []

# run loop
for i in range(len(sim)):
    sim_i = sim[i]
    fn_conc = '../GC_data/run'+sim_i+'/GEOSChem.SpeciesConc.2015_m.nc4'
    ds_conc = xr.open_dataset(fn_conc)
    lon = ds_conc.lon
    lat = ds_conc.lat
    Hg0_yr = ds_sel_yr(ds_conc, 'SpeciesConc_Hg0', 2015) 
    
    # Now more traceable
    R = 8.314462 # m^3 Pa K^-1 mol ^-1
    MW_Hg = 200.59 # g mol^-1
    ng_g = 1e9 # ng/g
    pg_ng = 1e3 # pg/ng
    
    stdpressure = 101325 # Pascals
    stdtemp = 273.15 # Kelvins
    
    unit_conv = stdpressure / R / stdtemp * MW_Hg * ng_g # converter from vmr to ng m^-3
    
        
    Hg0 = annual_avg(Hg0_yr) * unit_conv # ng m^-3
    
    Hg0_conc_a.append(Hg0) # ng m^-3
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
lat_min = -75
lat_max = 30
lon_min = -100
lon_max = -25
f, axes = plt.subplots(2, 2, figsize=[8,12],subplot_kw=dict(projection=ccrs.PlateCarree()),
                              gridspec_kw=dict(hspace=0.4, wspace=0.00))
axes = axes.flatten()
titles = ['BAU', 'GOV','SAV','RFR']

# Create colormap
for i in range(4):
    diff = Hg0_conc_a[i+1] - Hg0_conc_a[0]
    diff_max = abs(diff).max().values
    h = axes[i].pcolormesh(lon_pcolormesh_midpoints, lat_pcolormesh_midpoints,
                   diff.isel(lev=0), cmap ='bwr',
                   vmin=-diff_max, vmax=diff_max,
                   rasterized = True)
    axes[i].set_title(titles[i], fontsize=22, fontweight='bold')
    axes[i].set_xlim([lon_min, lon_max])
    axes[i].set_ylim([lat_min, lat_max])
    axes[i].coastlines(linewidth=1.5)

    cbar = plt.colorbar(h, extend='neither',orientation='horizontal',
                      ax=axes[i], pad=0.03, fraction=0.04)
    cbar.set_label('[Hg$^0$] (ng m$^{-3}$)', fontsize=20)   
    cbar.ax.tick_params(labelsize=18)

f.savefig('Figures/SI_Hg0_Amazon.pdf',bbox_inches = 'tight')
