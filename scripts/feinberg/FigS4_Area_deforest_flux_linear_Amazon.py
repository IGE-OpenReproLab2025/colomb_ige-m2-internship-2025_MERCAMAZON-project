#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 2023
Check if area deforested in Amazon is linear with flux
@author: arifeinberg
"""

import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import xarray as xr
from sklearn.linear_model import LinearRegression
#%% relate area deforested to flux change
scenario = ['HIST','BAU', 'GOV', 'SAV']
# FOR TOMORROW, CONSIDER HG0 DEPOSITION ALONE
# erosion_flux = np.array([63.6,84.5, 72.6, 124.7])
emiss_flux = np.array([54.55071828346036,89.56596268079161, 70.41627539662036, 143.78284107423264]) 
dep_flux = np.array([570.513986133593,465.6423721764879, 523.6492170874038, 211.15800609700054])
air_balance = np.array([-332.4664582334524,-179.74791888431753, -271.93805744966653, 108.93042228883377])
#%% calculate area deforested due to flux change
area_deforest= np.zeros(4) # store area changed here in km2
# load gbox areas
fn_gb = 'misc_Data/gbox_025x025_Olson.nc'
ds_gb = xr.open_dataset(fn_gb)
gbox_area = ds_gb.cell_area # grid_area
gbox_area_v = gbox_area.values # grid_area

# calculate areas
ds_HIST = xr.open_dataset('../input_data/HIST/Olson_Land_Type_Masks.025x025.BAU_200212.nc')
HIST_Olson = ds_HIST.LANDTYPE58.squeeze()

ds_BAU = xr.open_dataset('../input_data/BAU/Olson_Land_Type_Masks.025x025.BAU_205012.nc')
BAU_Olson = ds_BAU.LANDTYPE58.squeeze()
BAU_mask_change = BAU_Olson - HIST_Olson # changed area
area_deforest[1] = (BAU_mask_change * gbox_area).sum().values / 1e6 # km2

ds_GOV = xr.open_dataset('../input_data/GOV/Olson_Land_Type_Masks.025x025.GOV_205012.nc')
GOV_Olson = ds_GOV.LANDTYPE58.squeeze()
GOV_mask_change = GOV_Olson - HIST_Olson # changed area
area_deforest[2] = (GOV_mask_change * gbox_area).sum().values / 1e6 # km2

ds_SAV = xr.open_dataset('../input_data/SAV/Olson_2001_Land_Type_Masks.025x025.generic_sav58.nc')
SAV_Olson = ds_SAV.LANDTYPE58.squeeze()
SAV_mask_change = SAV_Olson - HIST_Olson # changed area
area_deforest[3] = (SAV_mask_change * gbox_area).sum().values / 1e6 # km2
#%% create linear regression models
model = LinearRegression()
model.fit(area_deforest.reshape(-1, 1), air_balance)
model_balance_rsq = model.score(area_deforest.reshape(-1, 1), air_balance)
model_balance_ypred = model.predict(area_deforest.reshape(-1, 1))

model.fit(area_deforest.reshape(-1, 1), emiss_flux)
model_emis_rsq = model.score(area_deforest.reshape(-1, 1), emiss_flux)
model_emis_ypred = model.predict(area_deforest.reshape(-1, 1))

model.fit(area_deforest.reshape(-1, 1), dep_flux)
model_dep_rsq = model.score(area_deforest.reshape(-1, 1), dep_flux)
model_dep_ypred = model.predict(area_deforest.reshape(-1, 1))

#%%
f, ax2 = plt.subplots(1, 1, figsize=[9,6],
                              gridspec_kw=dict(hspace=0.2, wspace=0.1))
ax2.plot(area_deforest/1e6, air_balance,'o', label='Land-air balance', c='#1f77b4')
ax2.plot(area_deforest/1e6, model_balance_ypred,'-', c='#1f77b4')
ax2.plot(area_deforest/1e6, emiss_flux,'o', label='Soil Hg$^0$ emission',  c='#2ca02c')
ax2.plot(area_deforest/1e6, model_emis_ypred,'-', c='#2ca02c')
ax2.plot(area_deforest/1e6, dep_flux,'o', label='Hg$^0$ dry deposition',  c='#d62728')
ax2.plot(area_deforest/1e6, model_dep_ypred,'-', c='#d62728')

ax2.set_ylabel('Flux (Mg yr$^{-1}$)', fontsize=15)
ax2.set_xlabel('Area deforested (10$^6$ km$^2$)', fontsize=15)
ax2.tick_params(axis='both', which='major', labelsize=13)
ax2.legend(fontsize=13, loc='lower right')
ax2.set_xlim([-0.1,7])
ax2.text(5, -50, 'R$^2$ = ' + str(round(model_balance_rsq,3)), c='#1f77b4',
         fontweight='bold',fontsize=13, horizontalalignment='left', verticalalignment='center')
ax2.text(5, 160, 'R$^2$ = ' + str(round(model_emis_rsq,3)), c='#2ca02c',
         fontweight='bold',fontsize=13, horizontalalignment='left', verticalalignment='center')
ax2.text(5, 340, 'R$^2$ = ' + str(round(model_dep_rsq,3)), c='#d62728',
         fontweight='bold',fontsize=13, horizontalalignment='left', verticalalignment='center')
f.savefig('Figures/lin_relationship_defo_flux.pdf',bbox_inches = 'tight')
