#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  9 08:29:33 2023
Plot bar plot of country-level and total emissions from deforestation
@author: arifeinberg
"""
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

#%% load country-level emissions
country_names = np.loadtxt("misc_Data/country_names_regionmask.csv", delimiter=',', dtype='str')

country_Hg_flux_15 = np.loadtxt("misc_Data/country_DFR_15.csv")
country_Hg_flux_30 = np.loadtxt("misc_Data/country_DFR_30.csv")
country_Hg_flux_45 = np.loadtxt("misc_Data/country_DFR_45.csv")
country_Hg_flux_60 = np.loadtxt("misc_Data/country_DFR_60.csv")

country_Hg_flux_15_low = np.loadtxt("misc_Data/country_DFR_15_low.csv")
country_Hg_flux_30_low = np.loadtxt("misc_Data/country_DFR_30_low.csv")
country_Hg_flux_45_low = np.loadtxt("misc_Data/country_DFR_45_low.csv")
country_Hg_flux_60_low = np.loadtxt("misc_Data/country_DFR_60_low.csv")

country_Hg_flux_15_high = np.loadtxt("misc_Data/country_DFR_15_high.csv")
country_Hg_flux_30_high = np.loadtxt("misc_Data/country_DFR_30_high.csv")
country_Hg_flux_45_high = np.loadtxt("misc_Data/country_DFR_45_high.csv")
country_Hg_flux_60_high = np.loadtxt("misc_Data/country_DFR_60_high.csv")
#%% sort countries and look at top ten - descending
s_idx_15 = np.flip(np.argsort(country_Hg_flux_15)) # indices for sort
s_idx_30 = np.flip(np.argsort(country_Hg_flux_30)) # indices for sort
s_idx_45 = np.flip(np.argsort(country_Hg_flux_45)) # indices for sort
s_idx_60 = np.flip(np.argsort(country_Hg_flux_60)) # indices for sort

# top 15 countries for 45 yr horizon
n = 15
country_45 = country_names[s_idx_45[:n]]

print(country_45)
# simplify country names
country_45[country_45=='Dem. Rep. Congo'] = 'DR Congo'
#%% calculate errors
# global total
glob_Hg_flux_60_err = np.zeros((2,1)) # for error bars
glob_Hg_flux_60_err[0,0] = np.sum(country_Hg_flux_60) - np.sum(country_Hg_flux_60_low)# low
glob_Hg_flux_60_err[1,0] = np.sum(country_Hg_flux_60_high) - np.sum(country_Hg_flux_60) # high
glob_Hg_flux_45_err = np.zeros((2,1)) # for error bars
glob_Hg_flux_45_err[0,0] = np.sum(country_Hg_flux_45) - np.sum(country_Hg_flux_45_low)# low
glob_Hg_flux_45_err[1,0] = np.sum(country_Hg_flux_45_high) - np.sum(country_Hg_flux_45) # high
glob_Hg_flux_30_err = np.zeros((2,1)) # for error bars
glob_Hg_flux_30_err[0,0] = np.sum(country_Hg_flux_30) - np.sum(country_Hg_flux_30_low)# low
glob_Hg_flux_30_err[1,0] = np.sum(country_Hg_flux_30_high) - np.sum(country_Hg_flux_30) # high
glob_Hg_flux_15_err = np.zeros((2,1)) # for error bars
glob_Hg_flux_15_err[0,0] = np.sum(country_Hg_flux_15) - np.sum(country_Hg_flux_15_low)# low
glob_Hg_flux_15_err[1,0] = np.sum(country_Hg_flux_15_high) - np.sum(country_Hg_flux_15) # high

# country total
country_Hg_flux_60_err = np.zeros((2,len(country_names))) # for error bars
country_Hg_flux_60_err[0,:] = country_Hg_flux_60 - country_Hg_flux_60_low # low
country_Hg_flux_60_err[1,:] = country_Hg_flux_60_high - country_Hg_flux_60 # high
country_Hg_flux_45_err = np.zeros((2,len(country_names))) # for error bars
country_Hg_flux_45_err[0,:] = country_Hg_flux_45 - country_Hg_flux_45_low # low
country_Hg_flux_45_err[1,:] = country_Hg_flux_45_high - country_Hg_flux_45 # high
country_Hg_flux_30_err = np.zeros((2,len(country_names))) # for error bars
country_Hg_flux_30_err[0,:] = country_Hg_flux_30 - country_Hg_flux_30_low # low
country_Hg_flux_30_err[1,:] = country_Hg_flux_30_high - country_Hg_flux_30 # high
country_Hg_flux_15_err = np.zeros((2,len(country_names))) # for error bars
country_Hg_flux_15_err[0,:] = country_Hg_flux_15 - country_Hg_flux_15_low # low
country_Hg_flux_15_err[1,:] = country_Hg_flux_15_high - country_Hg_flux_15 # high

#%% Plot bar graph
fig, ax = plt.subplots(1, 2, figsize=[15,7],
                       gridspec_kw=dict(hspace=0.1, wspace=0.2, width_ratios=[1,6]))
fig.subplots_adjust(bottom=0.2, top=0.9, left=0.12, right=0.98) # <-- Change the 0.02 to work for your plot.
ax = ax.flatten()
# global
ax[0].bar(1, np.sum(country_Hg_flux_60), width=0.15, color='#006d2c', edgecolor='k')
ax[0].bar(1.2, np.sum(country_Hg_flux_45), width=0.15, color='#238b45', edgecolor='k')
ax[0].bar(1.4, np.sum(country_Hg_flux_30), width=0.15, color='#74c476', edgecolor='k')
ax[0].bar(1.6, np.sum(country_Hg_flux_15), width=0.15, color='#bae4b3', edgecolor='k')
# add errors
ax[0].errorbar(1,np.sum(country_Hg_flux_60),yerr=glob_Hg_flux_60_err, capthick=1, capsize=4, elinewidth=1, 
            ecolor='k')
ax[0].errorbar(1.2,np.sum(country_Hg_flux_45),yerr=glob_Hg_flux_45_err, capthick=1, capsize=4, elinewidth=1, 
            ecolor='k')
ax[0].errorbar(1.4,np.sum(country_Hg_flux_30),yerr=glob_Hg_flux_30_err, capthick=1, capsize=4, elinewidth=1, 
            ecolor='k')
ax[0].errorbar(1.6,np.sum(country_Hg_flux_15),yerr=glob_Hg_flux_15_err, capthick=1, capsize=4, elinewidth=1, 
            ecolor='k')

x_labs  = ['Global']
ax[0].set_xticks([1.3])
ax[0].set_xticklabels(x_labs, fontsize=17)
ax[0].tick_params(axis='y', which='major', labelsize=17)
ax[0].set_ylabel('Deforestation emissions (Mg Hg yr$^{-1}$)',  fontsize=18)
ax[0].set_xlim([0.6,2])
ax[0].text(0., 1.04, 'a', fontweight='bold',fontsize=17,
             horizontalalignment='left', verticalalignment='center',
             transform=ax[0].transAxes)
ax[0].set_yscale('log')

# country-level
ax[1].bar(np.arange(1.0,n+0.9), country_Hg_flux_60[s_idx_45[:n]], width=0.15, color='#006d2c', edgecolor='k')
ax[1].bar(np.arange(1.2,n+1.1), country_Hg_flux_45[s_idx_45[:n]], width=0.15, color='#238b45', edgecolor='k')
ax[1].bar(np.arange(1.4,n+1.3), country_Hg_flux_30[s_idx_45[:n]], width=0.15, color='#74c476', edgecolor='k')
ax[1].bar(np.arange(1.6,n+1.5), country_Hg_flux_15[s_idx_45[:n]], width=0.15, color='#bae4b3', edgecolor='k')
# add errors
ax[1].errorbar(np.arange(1.0,n+0.9),country_Hg_flux_60[s_idx_45[:n]],
               yerr=country_Hg_flux_60_err[:,s_idx_45[:n]], ls='none',
               capthick=1, capsize=2, elinewidth=1, 
            ecolor='k')
ax[1].errorbar(np.arange(1.2,n+1.1),country_Hg_flux_45[s_idx_45[:n]],
               yerr=country_Hg_flux_45_err[:,s_idx_45[:n]], ls='none',
               capthick=1, capsize=2, elinewidth=1, 
            ecolor='k')
ax[1].errorbar(np.arange(1.4,n+1.3),country_Hg_flux_30[s_idx_45[:n]],
               yerr=country_Hg_flux_30_err[:,s_idx_45[:n]], ls='none',
               capthick=1, capsize=2, elinewidth=1, 
            ecolor='k')
ax[1].errorbar(np.arange(1.6,n+1.5),country_Hg_flux_15[s_idx_45[:n]],
               yerr=country_Hg_flux_15_err[:,s_idx_45[:n]], ls='none',
               capthick=1, capsize=2, elinewidth=1, 
            ecolor='k')

ax[1].set_yscale('log')

x_vals = np.arange(1.3,n+1.3)
#x_vals = np.concatenate((np.arange(1,17), np.arange(18,20), np.arange(21,23)), axis=None)
x_labs  = country_45
ax[1].set_xticks(x_vals)
ax[1].set_xticklabels(x_labs, rotation = 45, fontsize=17, ha='right', rotation_mode="anchor")
#ax.axvline(17, ls='solid', lw=1, color='k')
#ax.axvline(20, ls='dashed', lw=1, color='k', ymax=0.88)
ax[1].tick_params(axis='y', which='major', labelsize=17)
ax[1].set_ylabel('Deforestation emissions (Mg Hg yr$^{-1}$)',  fontsize=18)
#ax[1].set_title('Emissions considering different deforestation time horizons',  fontsize=18)

ax[1].set_xlim([0.6,16])
ax[1].legend(['60 yr (1955–2014)', '45 yr (1970–2014)', '30 yr (1985–2014)',
            '15 yr (2000–2014)'], fontsize=18)
ax[1].text(0., 1.04, 'b', fontweight='bold',fontsize=17,
             horizontalalignment='left', verticalalignment='center',
             transform=ax[1].transAxes)

fig.savefig('Figures/Bar_country_DFR_v2.pdf',bbox_inches = 'tight')
