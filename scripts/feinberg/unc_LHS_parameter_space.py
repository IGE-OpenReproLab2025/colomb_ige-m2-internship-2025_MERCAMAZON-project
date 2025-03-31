#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 2023
Create Latin Hypercube design of parameter space, to be used in uncertainty analysis for land use change paper
@author: arifeinberg
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
import numbers
from smt.sampling_methods import LHS

#%% Limits of parameters, setting up experimental design
# limits of parameters
soil_param = [1,101] # soil emission parametrization number
LAI_replace= [10,90] # percentile of LAI with which to replace data
f0_Amazon= [np.log10(0.01), np.log10(0.5)] # reactivity for dry deposition (Amazon)
f0_tropics_RF = [np.log10(1e-5), np.log10(0.2)] # reactivity for dry deposition (other tropical rainforests)
f0_else = [1e-5, 5e-5] # reactivity for dry deposition (non-rainforest)
bb_EF = [233.3, 410] # Biomass burning emission factor for Hg0,
                       # lower and upper limits adjusted down 50% for Carpi et al. (2014),
                       # soil emissions in first year after fire
# sample using Latin Hypercube Sampling
xlimits = np.array([soil_param, LAI_replace, f0_Amazon, f0_tropics_RF, f0_else, bb_EF])
sampling = LHS(xlimits=xlimits) # do LHS sampling

num = 100 # take 100 samples
x = sampling(num) #  get params
# adjust x to have proper values for each variable
x_adj = x
# switch from log-values to linear
x_adj[:,2] = 10**x_adj[:,2]
x_adj[:,3] = 10**x_adj[:,3]
# take integer for first column of soil param number
x_adj[:,0] = np.floor(x_adj[:,0])
print(len(np.unique(x_adj[:,0]))) # check that all soil params are used
np.savetxt("misc_Data/LHS_sampling_dep_emis.csv", x_adj, delimiter=",") # already sampled, just to get params
