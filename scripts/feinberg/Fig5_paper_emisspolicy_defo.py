#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 3 2023
Compare fluxes from deforestation to Hg emissions policies
@author: arifeinberg
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib
#matplotlib.use('Agg')

#%% global emissions
emiss_GMA = 2222.55
emiss_defo_45 = 216.77614895916028 # deforestation (1970-2014) emissions in 2015
# Uncertainties
emiss_GMA_err = np.zeros((2,1)) # for error bars
emiss_GMA_err[0,0] = emiss_GMA - 2000 # low
emiss_GMA_err[1,0] = 2820 - emiss_GMA # high

emiss_defo_45_err = np.zeros((2,1)) # for error bars
emiss_defo_45_err[0,0] = emiss_defo_45 - 133.67848762 # low
emiss_defo_45_err[1,0] = 1649.66520133 - emiss_defo_45 # high


#%% policies - defo
Am_cons = 92.190138565349
Gb_rfr = 97.867
# Uncertainties
Am_cons_err = np.zeros((2,1)) # for error bars
Am_cons_err[0,0] = Am_cons - 59.17835055195607 # low
Am_cons_err[1,0] = 234.19363142909015 - Am_cons # high

Gb_rfr_err = np.zeros((2,1)) # for error bars
Gb_rfr_err[0,0] = Gb_rfr -64.06478412780498 # low
Gb_rfr_err[1,0] = 449.3995726760292 - Am_cons # high

#%% policies - primary past
# MATS standards (25 tons - 86% reduction to 4 tons Table 4 in EPA 2019 - https://www.govinfo.gov/content/pkg/FR-2019-02-07/pdf/2019-00936.pdf)
MATS = 25.
# China air pollution control 2013 to 2017 (127 tons - Liu et al. 2019)
China_APC = 127

# Canada reduced Hg emissions by 4.81 tons /yr over 2007–2017 due to Risk Management Strategy for Mercury. Table 1 in ECCC 2020:
#https://www.canada.ca/en/environment-climate-change/services/management-toxic-substances/evaluation-effectiveness-risk-management-measures-mercury.html. 
Canada_RMS = 4.81
#%% policies - primary projected
# Pacyna et al. 2016 - 940 tons/yr - new policy vs. Current policy 2035.
# Global_NP_CP = 940
#ASGM - Bruno et al. 2023 group paper (262 tons - reduction of 50% of emissions in 50% of ASGM sites), 
ASGM_50_50 = 262
# China Minamata + Stringent Climate Policy (69 tons - Mulvaney et al. 2020)
China_Minamata = 69
# Rafaj et al. 2014 - 44 tons for European decarbonization
Europe_decarb = 44.
#%% Plot stacked bar graph
fig, ax = plt.subplots(2, 1, figsize=[13,7],
                       gridspec_kw=dict(hspace=0.5, wspace=0.1, height_ratios=[1,4]))
fig.subplots_adjust(bottom=0.1, top=0.9, left=0.5, right=0.98) # <-- Change the 0.02 to work for your plot.

ax = ax.flatten()
ax[0].barh(1, emiss_GMA, height=0.15, color='grey', edgecolor='k') # Global Hg emissions
ax[0].errorbar(emiss_GMA, 1, xerr=emiss_GMA_err, capthick=1, capsize=4, elinewidth=1, 
            ecolor='k')
ax[0].barh(0.8, emiss_defo_45, height=0.15, color='#238b45', edgecolor='k')
ax[0].errorbar(emiss_defo_45, 0.8, xerr=emiss_defo_45_err, capthick=1, capsize=4, elinewidth=1, 
            ecolor='k')

#x_vals = np.concatenate((np.arange(1,17), np.arange(18,20), np.arange(21,23)), axis=None)
y_labs  = ['Net emissions in 2015 from 1970–2014 deforestation', # \n(this study)',
           'Primary anthropogenic emissions for 2015'] #\n(GMA 2018)']
ax[0].set_yticks([0.8, 1])
ax[0].set_yticklabels(y_labs, fontsize=13, ha='right')
#ax.axvline(17, ls='solid', lw=1, color='k')
#ax.axvline(20, ls='dashed', lw=1, color='k', ymax=0.88)
ax[0].tick_params(axis='x', which='major', labelsize=15)
ax[0].set_xlabel('Emissions (Mg yr$^{-1}$)',  fontsize=15)
ax[0].set_title('Global Hg emissions',  fontsize=15, fontweight='bold')
ax[0].text(0., 1.14, 'a', fontweight='bold',fontsize=17,
             horizontalalignment='left', verticalalignment='center',
             transform=ax[0].transAxes)

# Policy plot
ax[1].barh(1, Am_cons, height=0.15, color='#238b45', edgecolor='k') # Global Hg emissions
ax[1].errorbar(Am_cons, 1, xerr=Am_cons_err, capthick=1, capsize=4, elinewidth=1, 
            ecolor='k')
ax[1].barh(0.8, Gb_rfr, height=0.15, color='#238b45', edgecolor='k')
ax[1].errorbar(Gb_rfr, 0.8, xerr=Gb_rfr_err, capthick=1, capsize=4, elinewidth=1, 
            ecolor='k')

ax[1].barh(0.6, China_APC, height=0.15, color='#525252', edgecolor='k') # Global Hg emissions
ax[1].barh(0.4, MATS, height=0.15, color='#525252', edgecolor='k') # Global Hg emissions
ax[1].barh(0.2, Canada_RMS, height=0.15, color='#525252', edgecolor='k') # Global Hg emissions

ax[1].barh(0.0, ASGM_50_50, height=0.15, color='#969696', edgecolor='k') # Global Hg emissions
ax[1].barh(-0.2, China_Minamata, height=0.15, color='#969696', edgecolor='k') # Global Hg emissions
ax[1].barh(-0.4, Europe_decarb, height=0.15, color='#969696', edgecolor='k') # Global Hg emissions

#x_vals = np.concatenate((np.arange(1,17), np.arange(18,20), np.arange(21,23)), axis=None)
y_labs  = ['European decarbonization of energy system by 2050', #\n(Rafaj et al. 2014)',
           'China implementation of Minamata + stringent climate policy by 2030', #\n(Mulvaney et al. 2020)',
           'Global reduction of 50% of emissions at 50% ASGM sites', #\n(Bruno et al. 2023)',
           'Canada Risk Management Strategy for Mercury 2007-2017', #\n(ECCC 2020)',
           'US Mercury and Air Toxics Standards (MATS) 2010-2017', #\n(EPA 2019)',
            'China air pollution control 2013–2017', #\n(Liu et al. 2019)',
           'Global reforestation scenario', #\n(this study)',
           'Amazon conservation by 2050'] #\n(this study)']
ax[1].set_yticks([-.4,-.2,0, 0.2, 0.4, 0.6, 0.8, 1])
ax[1].set_yticklabels(y_labs, fontsize=13, ha='right')
#ax.axvline(17, ls='solid', lw=1, color='k')
#ax.axvline(20, ls='dashed', lw=1, color='k', ymax=0.88)
ax[1].tick_params(axis='x', which='major', labelsize=15)
ax[1].set_xlabel('Emissions reductions (Mg yr$^{-1}$)',  fontsize=15)
ax[1].set_title('Policy impacts on net atmospheric Hg fluxes',  fontsize=15, fontweight='bold')
ax[1].text(425,0.91, 'Land use scenarios', horizontalalignment='right',
      verticalalignment='center', fontsize=15, fontweight='bold', color='#238b45')
ax[1].text(425,0.35, 'Past policy impacts', horizontalalignment='right',
      verticalalignment='center', fontsize=15, fontweight='bold', color='#525252')
ax[1].text(425,-0.25, 'Other policy scenarios', horizontalalignment='right',
      verticalalignment='center', fontsize=15, fontweight='bold', color='#969696')
ax[1].text(0., 1.04, 'b', fontweight='bold',fontsize=17,
             horizontalalignment='left', verticalalignment='center',
             transform=ax[1].transAxes)
fig.savefig('Figures/paper_defo_policy_emiss_bar.pdf',bbox_inches = 'tight')
