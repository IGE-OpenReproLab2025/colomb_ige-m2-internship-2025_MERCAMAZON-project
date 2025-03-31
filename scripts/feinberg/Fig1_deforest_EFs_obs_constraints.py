#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  24 2023
Compare deforestation Hg EFs calculated in GEOS-Chem with observed values
@author: arifeinberg
"""
#%% import packages
import numpy as np
import matplotlib.pyplot as plt
#%% EFs from model calculations
realm_names = ['Amazon', 'China', 'Nearctic', 'Afrotropic', 'Indomalayan',  'Neotropic','Palearctic','Australasia']
# realm_names = ['Amazon', 'Afrotrop', 'Indomalay', 'China', 'Neotrop','Palearc','Austral','Nearc']

# EF data and uncertainties calculated in LUC_global_factors_emiss_v2.py
LUC_EF_defo = [7.235329241913561e-05,2.3404428517737673e-05, 1.0623755118946325e-05,
               1.071806689794155e-05, 2.2819489051113326e-05, 7.404479045967743e-06,
          2.395577187857226e-06, 6.095471571713798e-06]
Soil_EF_defo = [1.4104624583103027e-05, 1.453324446692011e-05, 6.8414553063659995e-06, 
                7.965166356949215e-06, 1.2031028933539869e-05, 5.359641752377664e-06, 
                1.368712872976517e-06, 1.7143365857668656e-06]
LUC_EF_defo_err = np.array([[2.77764247e-05, 6.03522956e-06, 3.51837723e-06, 7.87637585e-06,
        7.34510428e-06, 2.58673852e-06, 2.31978907e-06, 5.26901446e-06],
       [1.27326075e-04, 2.03004305e-04, 5.14893001e-05, 1.04295126e-04,
        1.90985035e-04, 5.00345520e-05, 2.10230703e-05, 4.76827583e-05]])

Soil_EF_defo_err = np.array([[1.45863498e-06, 2.56310088e-06, 1.42419025e-06, 4.15465678e-06,
        1.23593938e-06, 1.38929819e-06, 1.18763283e-06, 9.32424964e-07],
       [9.86179216e-05, 1.86797382e-04, 4.66670003e-05, 6.30033105e-05,
        1.28840152e-04, 4.17886721e-05, 1.37715245e-05, 5.01275249e-06]])

n = len(realm_names)
#%% EFs from observations
# Total estimates
# from Amazon soil measurements of paired forested-deforested sites
amazon_mb_u = 3.18e-4 # mean
amazon_mb_err = np.zeros((2,1)) # for error bars
amazon_mb_err[0,0] = amazon_mb_u - 2.50e-5 # low
amazon_mb_err[1,0] = 7.21e-4 - amazon_mb_u # high

# China - Wang et al. (2016)
wang_16_u = 3.86e-5
wang_16_err = np.zeros((2,1)) # for error bars
wang_16_err[0,0] = wang_16_u - 5.26e-6 # low
wang_16_err[1,0] = 6.75e-5 - wang_16_u # high

# Nearctic - Homann et al. (2015)
homann_15_u = 1.49e-5
homann_15_err = np.zeros((2,1)) # for error bars
homann_15_err[0,0] = homann_15_u - 8.57e-6 # low
homann_15_err[1,0] = 2.22e-5 - homann_15_u # high

# Nearctic - Gamby et al. (2015)
gamby_15_u = 2.93e-5

# Soil estimates
# Amazon - Carpi et al. (2014)
carpi_14_am_u = (1.65e-5+2.24e-4)/2 # not actually shown, just range
carpi_14_am_err = np.zeros((2,1)) # for error bars
carpi_14_am_err[0,0] = carpi_14_am_u - 1.65e-5 # low
carpi_14_am_err[1,0] = 2.24e-4 - carpi_14_am_u # high

# Amazon - Magarelli and Fostier (2005)
magarelli_05_u = 1.82e-5
magarelli_05_err = np.zeros((2,1)) # for error bars
magarelli_05_err[0,0] = magarelli_05_u - 9.20e-6 # low
magarelli_05_err[1,0] = 2.80e-5 - magarelli_05_u # high

# Amazon - Almeida et al. (2009)
almeida_09_u = 3.5e-5

# Nearctic - Mazur et al. (2014)
mazur_14_u = (1.65e-5+3.90e-5)/2 # not actually shown, just range
mazur_14_err = np.zeros((2,1)) # for error bars
mazur_14_err[0,0] = mazur_14_u - 1.65e-5 # low
mazur_14_err[1,0] = 3.90e-5 - mazur_14_u # high

# Nearctic - Eckley et al. (2021)
eckley_21_u = 7.7e-6

# Nearctic - Carpi et al. (2014)
carpi_14_ny_u = 5.4e-5

# China - Ma et al. (2013)
ma_13_u = 5.71e-5

#%% Plot bar graph
fig, ax = plt.subplots(2, 1, figsize=[12,9], sharex=True,
                       gridspec_kw=dict(hspace=0.1, wspace=0.1))
fig.subplots_adjust(bottom=0.15, top=0.95, left=0.12, right=0.98) # <-- Change the 0.02 to work for your plot.
ax = ax.flatten()
# Total emiss subplot
# Model estimates
ax[0].errorbar(np.arange(1,n+1), LUC_EF_defo, yerr=LUC_EF_defo_err,
              fmt='o',capthick=2, capsize=4, elinewidth=2,
              ecolor='k', mfc='k', mec='k', markersize=8,
        label='GEOS-Chem model')

# # add blank line in legend plot
# l = ax[0].plot([1],[1e-5],color="none", label=' ')

# Observed estimates
ax[0].errorbar(1.15, amazon_mb_u,yerr=amazon_mb_err,
            fmt='s', capthick=2, capsize=4, elinewidth=2, 
            ecolor='grey', mfc='grey', mec='k', markersize=8, label='Amazon refs (see caption)')
ax[0].errorbar(2.15, wang_16_u,yerr=wang_16_err,
            fmt='h', capthick=2, capsize=4, elinewidth=2, 
            ecolor='grey', mfc='grey', mec='k', markersize=8, label='Wang et al. 2016')
ax[0].errorbar(3.15, homann_15_u,yerr=homann_15_err,
            fmt='d', capthick=2, capsize=4, elinewidth=2, 
            ecolor='grey', mfc='grey', mec='k', markersize=8, label='Homann et al. 2015')
ax[0].errorbar(3.3, gamby_15_u,fmt='^', mfc='grey', mec='k', ms=8,elinewidth=0,
        label='Gamby et al. 2015')

ax[0].set_yscale('log')
#ax[0].set_xticks([])
#ax.axvline(20, ls='dashed', lw=1, color='k', ymax=0.88)
# add lines separating realms
ax[0].axvline(1.725, ls='dashed', lw=0.5, color='#bdbdbd')
ax[0].axvline(2.725, ls='dashed', lw=0.5, color='#bdbdbd')
ax[0].axvline(3.725, ls='dashed', lw=0.5, color='#bdbdbd')
ax[0].axvline(4.725, ls='dashed', lw=0.5, color='#bdbdbd')
ax[0].axvline(5.725, ls='dashed', lw=0.5, color='#bdbdbd')
ax[0].axvline(6.725, ls='dashed', lw=0.5, color='#bdbdbd')
ax[0].axvline(7.725, ls='dashed', lw=0.5, color='#bdbdbd')

ax[0].tick_params(axis='y', which='major', labelsize=17)
ax[0].set_ylabel('Emission Factor (Mg m$^{-2}$ yr$^{-1}$)',  fontsize=16)
#ax.set_ylim([0,100])
ax[0].legend(fontsize=13)
ax[0].set_title('Total', fontsize=18, fontweight='bold')
ax[0].xaxis.set_tick_params(bottom=False)

# Soil emiss subplot
# model estimates
ax[1].errorbar(np.arange(1,n+1), Soil_EF_defo, yerr=Soil_EF_defo_err,
              fmt='o',capthick=2, capsize=4, elinewidth=2,
              ecolor='k', mfc='k', mec='k', markersize=8,
        label='GEOS-Chem model')

# Observed estimates
ax[1].errorbar(1.15, carpi_14_am_u,yerr=carpi_14_am_err,
            fmt='None', capthick=2, capsize=4, elinewidth=2, 
            ecolor='#525252', mfc='#525252', mec='k', markersize=8, label='Carpi et al. 2014')
ax[1].errorbar(1.3, magarelli_05_u,yerr=magarelli_05_err,
            fmt='v', capthick=2, capsize=4, elinewidth=2, 
            ecolor='grey', mfc='grey', mec='k', markersize=8, label='Magarelli and Fostier 2005')
ax[1].errorbar(1.45, almeida_09_u,fmt='X', mfc='grey', mec='k', ms=8,elinewidth=0,
        label='Almeida et al. 2009')
ax[1].errorbar(2.15, ma_13_u,fmt='p', mfc='grey', mec='k', ms=9,elinewidth=0,
        label='Ma et al. 2013')
ax[1].errorbar(3.15, mazur_14_u,yerr=mazur_14_err,
            fmt='None', capthick=2, capsize=4, elinewidth=2, 
            ecolor='grey', mfc='grey', mec='k', markersize=8, label='Mazur et al. 2014')
ax[1].errorbar(3.30, carpi_14_ny_u,fmt='P', mfc='grey', mec='k', ms=8,elinewidth=0,
        label='Carpi et al. 2014')
ax[1].errorbar(3.45, eckley_21_u,fmt='*', mfc='grey', mec='k', ms=10,elinewidth=0,
        label='Eckley et al. 2021')

ax[1].set_yscale('log')
x_vals = np.arange(1.225,n+1.225)
x_labs  = realm_names
ax[1].set_xticks(x_vals)
ax[1].set_xticklabels(x_labs, rotation = 45, fontsize=17, ha='right', rotation_mode="anchor")
# add lines separating realms
ax[1].axvline(1.725, ls='dashed', lw=0.5, color='#bdbdbd')
ax[1].axvline(2.725, ls='dashed', lw=0.5, color='#bdbdbd')
ax[1].axvline(3.725, ls='dashed', lw=0.5, color='#bdbdbd')
ax[1].axvline(4.725, ls='dashed', lw=0.5, color='#bdbdbd')
ax[1].axvline(5.725, ls='dashed', lw=0.5, color='#bdbdbd')
ax[1].axvline(6.725, ls='dashed', lw=0.5, color='#bdbdbd')
ax[1].axvline(7.725, ls='dashed', lw=0.5, color='#bdbdbd')

ax[1].tick_params(axis='y', which='major', labelsize=17)
ax[1].set_ylabel('Emission Factor (Mg m$^{-2}$ yr$^{-1}$)',  fontsize=16)
ax[1].set_xlim([0.725,8.725])
ax[1].legend(fontsize=14, ncol=2)
ax[1].set_title('Soil', fontsize=18, fontweight='bold')

fig.savefig('Figures/deforest_EFs_obs_compare.pdf',bbox_inches = 'tight')
