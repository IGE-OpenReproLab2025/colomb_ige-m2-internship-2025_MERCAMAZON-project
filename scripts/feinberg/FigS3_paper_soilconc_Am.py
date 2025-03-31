#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 10:34:33 2023
Plot box plots of forested and deforested concentrations
@author: arifeinberg
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.patches as patches
from scipy.stats import wilcoxon
#%% load deforestation impact data
fn_Hg= 'misc_Data/Hg_soil_loss_data.csv' # Hg conc data
df_Hg = pd.read_csv(fn_Hg)
Hg_for = df_Hg.iloc[:,2].values
Hg_df = df_Hg.iloc[:,3].values
Hg_diff = Hg_for-Hg_df
#%% Make boxplot plot
#color definitions
frc1 = '#33a02c' # dark forest
frc2 = '#b2df8a' # light forest
dfc1 = '#ff7f00' # dark deforest
dfc2 = '#fdbf6f' # light deforest

# stats definitions
IQR_for = np.percentile(Hg_for, 75) - np.percentile(Hg_for, 25) # IQR
P25_for = np.percentile(Hg_for, 25) # P25
P75_for = np.percentile(Hg_for, 75) # P75
med_for = np.percentile(Hg_for, 50) # median
P90_for = np.percentile(Hg_for, 90) # P90
P10_for = np.percentile(Hg_for, 10) # P10

IQR_df = np.percentile(Hg_df, 75) - np.percentile(Hg_df, 25) # IQR
P25_df = np.percentile(Hg_df, 25) # P25
P75_df = np.percentile(Hg_df, 75) # P75
med_df = np.percentile(Hg_df, 50) # median
P90_df = np.percentile(Hg_df, 90) # P90
P10_df = np.percentile(Hg_df, 10) # P10

# conduct Wilcoxon signed rank test for paired sample (non-parametric version of paired T test)
res = wilcoxon(Hg_for, Hg_df)

fig, ax = plt.subplots(figsize=[3.5,8])
fig.subplots_adjust(left=0.26, right=0.98) # <-- Change the 0.02 to work for your plot.

# manual boxplots
# forest
# median line
ax.plot([0.86, 1.14], [med_for, med_for], lw = 3, c = frc1, zorder=20)
# add line for P75-P90
ax.plot([1., 1.], [P75_for, P90_for], lw = 2, c = frc1 , zorder=0)
# add line for P10-P25
ax.plot([1., 1.], [P10_for, P25_for], lw = 2, c = frc1 , zorder=0)
# add whiskers for P10/P90
ax.plot([0.95, 1.05], [P90_for, P90_for], lw = 2, c = frc1 , zorder=0)
ax.plot([0.95, 1.05], [P10_for, P10_for], lw = 2, c = frc1 , zorder=0)
# add box for IQR
rect_for = patches.Rectangle((0.85, P25_for), 0.3, IQR_for, zorder=10,
                         linewidth=1, edgecolor=frc1, facecolor=frc2)
ax.add_patch(rect_for)

# deforest
# median line
ax.plot([1.46, 1.74], [med_df, med_df], lw = 3, c = dfc1, zorder=20)
# add line df P75-P90
ax.plot([1.6, 1.6], [P75_df, P90_df], lw = 2, c = dfc1 , zorder=0)
# add line df P10-P25
ax.plot([1.6, 1.6], [P10_df, P25_df], lw = 2, c = dfc1 , zorder=0)
# add whiskers df P10/P90
ax.plot([1.55, 1.65], [P90_df, P90_df], lw = 2, c = dfc1 , zorder=0)
ax.plot([1.55, 1.65], [P10_df, P10_df], lw = 2, c = dfc1 , zorder=0)
# add box df IQR
rect_df = patches.Rectangle((1.45, P25_df), 0.3, IQR_df, zorder=10,
                         linewidth=1, edgecolor=dfc1, facecolor=dfc2)
ax.add_patch(rect_df)

# draw connecting lines
for i in range(len(Hg_for)):
    ax.plot([1,1.6], [Hg_for[i], Hg_df[i]],'-', c= [0.7,0.7,0.7], zorder=25)
# add points for all samples
ax.plot([1]*len(Hg_for), Hg_for, '.',mec='k', c=frc1, zorder=30)
ax.plot([1.6]*len(Hg_for), Hg_df, '.',mec='k', c=dfc1, zorder=30)

ax.set_xticks([])

# ax.set_xticks([1.,1.6])
# ax.set_xticklabels(['Forest', 'Deforested'])
# ax.tick_params(axis='x', which='major', labelsize=17)
ax.tick_params(axis='y', which='major', labelsize=15)
ax.set_ylabel('Hg concentration (ng g$^{-1}$)', fontsize = 17)
ax.set_xlim([0.75,1.85])
ax.set_title('Amazon soil \n measurements \n (literature)', fontsize=16, fontweight='bold')
title_s = 'Amazon soil \n measurements \n (literature)'
#ax.text(1.45, 295, title_s, horizontalalignment='center',
#      verticalalignment='center', fontsize=15, fontweight = 'bold',color='black')

text1 = "Median diff:\n 25 ng g$^{-1}$"
text2 = "Wilcoxon:\n p < 0.001"
ax.text(1.5, 290, text1, horizontalalignment='center',
      verticalalignment='center', fontsize=14, color='black')
ax.text(1.5, 260, text2, horizontalalignment='center',
      verticalalignment='center', fontsize=14, color='black')
ax.text(1., 13, 'Forest', horizontalalignment='center',
      verticalalignment='center', fontsize=16, color=frc1)
ax.text(1.58, 13, 'Deforested', horizontalalignment='center',
      verticalalignment='center', fontsize=16, color=dfc1)
fig.savefig('Figures/paper_soil_conc_lit.pdf',bbox_inches = 'tight')
