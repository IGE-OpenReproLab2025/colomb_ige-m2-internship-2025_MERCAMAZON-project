#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 9 2022
Creating Olson land types and XLAI files for Soares scenarios
The .tif data for Soares-Filho06 can be downloaded here: https://doi.org/10.3334/ORNLDAAC/1153
@author: arifeinberg
"""
import os
import xarray as xr
import numpy as np
import netCDF4
import rioxarray
from cdo import *
cdo = Cdo()
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf

def weird_division1(n, d):
    return n / d if d else 1 # no change if denominator 0
# create ufunc from weird division1
udiv1 = np.frompyfunc(weird_division1, 2, 1)
#%% Setup loop
# since now perturb from June, need to run the full set of years
startyear = 2002
endyear = 2050
year_r = np.arange(startyear, endyear+1)

sce = 'BAU' # scenario: 'BAU' or 'GOV'
sce_l = sce.lower() # scenario: 'BAU' or 'GOV'

#%% First load 2002 geotif, used as reference for future changes
# Load GEOTiff files 
fn = '../Amazon_scenarios/' + sce + '_amazonia/amazon' + sce_l +'2002.tif'
ds = xr.open_dataset(fn, engine="rasterio")
landtypes = ds.band_data.squeeze()
# xr.plot.pcolormesh(landtypes)

# Get forest area for 2002 as reference
# Before remapping, create new variable that assigns area fraction for each land type
x_forest = xr.where(landtypes==2,1.0,0.0) # intact forest

# save dataArrays as dataset
ds_new = x_forest.to_dataset(name='x_forest')

# clean up data array so that readable by cdo
ds_new = ds_new.reset_coords(['band', 'spatial_ref'] ,drop=True)
ds_new = ds_new.rename({'x': 'lon','y': 'lat'})
ds_new = ds_new.reindex(lat=ds_new.lat[::-1]) # reverse lat dimension so that increasing

# add attributes for lat/lon
ds_new['lon'].attrs['units']="degrees east"
ds_new['lat'].attrs['units']="degrees north"

# Do regridding on 0.25 x 0.25 degree grid to compare with Olson
ds_remap = cdo.remapcon("misc_Data/geos.025x025_SA.grid",input=ds_new, returnXDataset=True)

# Get variables needed
x_forest_ref = ds_remap.x_forest.values # forested area in 2002
lat_rm = ds_remap.lat.values
lon_rm = ds_remap.lon.values

#%% Load necessary variables for Olson land map, LAI
# load gbox areas for South America
fn_gb = 'misc_Data/gbox_025x025_Olson.nc'
ds_gb = xr.open_dataset(fn_gb)
gbox_area = ds_gb.cell_area # grid_area

# required Olson land variables
fn_ols1 = 'misc_Data/Olson_2001_Drydep_Inputs.nc'
ds_ols1 = xr.open_dataset(fn_ols1)

DRYCOEFF = ds_ols1.DRYCOEFF.values # Baldocchi dry deposition polynomial coefficients
IOLSON = ds_ols1.IOLSON.values # Olson land type indices (+1)
IDRYDEP = ds_ols1.IDRYDEP.values # Dry deposition land types

IDEP = ds_ols1.IDEP.values # Mapping index: Olson land type ID to drydep ID

IOLSON_6 = np.asarray(np.where(IDEP==6)).flatten() # Indices of rainforest

# Land type areas
fn_ols2 = 'misc_Data/Olson_2001_Land_Type_Masks.025x025.generic.nc'
ds_ols2 = xr.open_dataset(fn_ols2)
# save as one np array, first dim is land types
Olson_landtype = ds_ols2.to_array().squeeze()

lon = ds_ols2.lon # longitude
lat = ds_ols2.lat # latitude

n_landtypes = 73 # number of Olson land types

# Land type LAI
fn_xlai  = 'misc_Data/Yuan_proc_MODIS_XLAI.025x025.2003.nc' # Use 2003 as reference year
ds_xlai = xr.open_dataset(fn_xlai)

# save as one data array, first dim is land types
XLAI_w = ds_xlai.to_array()

# restrict variables to area of South America
# select Amazon coordinate bounds
lat_min = lat_rm.min()
lat_max = lat_rm.max()
lon_min = lon_rm.min()
lon_max = lon_rm.max()

lon_l = len(lon_rm)
lat_l = len(lat_rm)

Olson_landtype_A = Olson_landtype[:,(lat<=lat_max)&(lat>=lat_min),(lon<=lon_max)&(lon>=lon_min)].values
XLAI_w_A = XLAI_w[:,:, (lat<=lat_max)&(lat>=lat_min), (lon<=lon_max)&(lon>=lon_min)].values
gbox_area_A = gbox_area[(lat<=lat_max) & (lat>=lat_min),(lon<=lon_max) & (lon>=lon_min)].values
#%% Prepare land change adjustment - calculate average savanna LAI to replace rainforest with
# Change tropical rainforest landclasses to savanna in Amazon
savanna_id = 58 # Olson id: 58 = Fields and Woody Savanna

# time variable
time_w = XLAI_w.time
time_l = len(time_w)

# create LAI threshold to average over
LAI_thres = 0.1 # if below this, then not reliable amount of savanna vegetation

# calculate timeseries of South American savanna average LAI
XLAI_avg = np.zeros(time_l)

# calculate XLAI average in savannah everywhere as loop:
for i in range(time_l):
    # Calculate average savanna LAI over Amazon domain
    XLAI_sel = XLAI_w_A[savanna_id, i, :, :]
    
    XLAI_avg[i] =  np.mean(XLAI_sel[XLAI_sel>LAI_thres] * gbox_area_A[XLAI_sel>LAI_thres]) \
        / np.mean(gbox_area_A[XLAI_sel>LAI_thres]) # take average of all mixed values

#%% Run loop
for i, yr in enumerate(year_r): # looping over years to create input files
    
    yr_s = str(yr) # year as string
    print(yr_s)
    # Input file names - LAI and Olson land types
    ifn_LT = fn_ols2
    ifn_LAI = fn_xlai

    # Output file names - LAI and Olson land types
    ofn_LT = 'misc_Data/' + sce + '/LAI/Olson_Land_Type_Masks.025x025.' + sce + '_' + yr_s +'.nc'
    ofn_LAI = 'misc_Data/' + sce + '/LAI/Yuan_proc_MODIS_XLAI.025x025.' + sce + '_' + yr_s +'.nc'

    if yr == 2002: # for first year, take reference values for 2003
                   # note use 2003 instead of 2002 because MODIS-LAI shows 
                   # measurement artifact before 2002
        os.system('cp ' + ifn_LT + ' ' + ofn_LT) # copy files to new file names
        os.system('cp ' + ifn_LAI + ' ' + ofn_LAI) # copy files to new file names
        
        # edit metadata so that time is correct for the desired year
        # units
        atr = "minutes since " + yr_s + "-01-01 00:00:00"
        cmd = 'ncatted -h -O -a  units,time,m,c,"' + atr +'" ' + ofn_LAI
        os.system(cmd)
        # begin_date
        atr = yr_s + "0101"
        cmd = 'ncatted -h -O -a  begin_date,time,m,c,"' + atr +'" ' + ofn_LAI
        os.system(cmd)
        
    else: # After yr 2002, need to take into account deforestation changes
        # Load GEOTiff files
        fn = '../Amazon_scenarios/' + sce + '_amazonia/amazon' + sce_l + yr_s +'.tif'
        ds = xr.open_dataset(fn, engine="rasterio")
        landtypes = ds.band_data.squeeze()
        # xr.plot.pcolormesh(landtypes)
        # Before remapping, create new variable that assigns area fraction for forest type
        x_forest = xr.where(landtypes==2,1.0,0.0) # intact forest
        # save dataArray as dataset
        ds_new = x_forest.to_dataset(name='x_forest')
        # clean up data array so that readable by cdo
        ds_new = ds_new.reset_coords(['band', 'spatial_ref'] ,drop=True)
        ds_new = ds_new.rename({'x': 'lon','y': 'lat'})
        ds_new = ds_new.reindex(lat=ds_new.lat[::-1]) # reverse lat dimension so that increasing
        # add attributes for lat/lon
        ds_new['lon'].attrs['units']="degrees east"
        ds_new['lat'].attrs['units']="degrees north"
        # Do regridding on 0.25 x 0.25 degree grid to compare with Olson
        ds_remap = cdo.remapcon("misc_Data/geos.025x025_SA.grid",input=ds_new, returnXDataset=True)
        # Get variables needed
        x_forest_rm = ds_remap.x_forest.values
        
        # Calculate change indices - how much was deforested relative to forest
        # in 2002 Soares dataset
        x_change_v = udiv1(x_forest_rm, x_forest_ref) # division, return 1 if 0 denominator

        # create new matrices for edited land use data
        XLAI_w_def = XLAI_w_A.copy()
        Olson_landtype_def = Olson_landtype_A.copy()

        # Before June (index 19)
        if yr > 2003: # For 2003 keep until June the 2002 LAI
            # Land change adjustment    
            # find indices where have deforestation
            inds_DF = np.asarray(np.where((x_change_v_prev<1 ))).T
            counter_no_RF = 0
            total_areas = np.zeros(73)
            # loop over these indices, correct to savanna
            for i in range(len(inds_DF)):
                lat_i = inds_DF[i,0]
                lon_i = inds_DF[i,1]
                # check if actually have rainforest in this grid cell, otherwise do nothing
                total_RF_area = sum(Olson_landtype_A[IOLSON_6,lat_i,lon_i])
                if total_RF_area > 0:
                    change_factor = x_change_v_prev[lat_i,lon_i]
                    # alter the LAI of the rainforest land types to adjust for area scaling of LAI
                    # index 19 is for start of June
                    XLAI_w_def[IOLSON_6,:19,lat_i,lon_i] = XLAI_w_def[IOLSON_6,:19,lat_i,lon_i] * change_factor
                    # replace lost forest by savanna type vegetation
                    # index 19 is for start of June
                    XLAI_w_def[savanna_id,:19, lat_i,lon_i] = XLAI_avg[:19] * Olson_landtype_def_prev[savanna_id,lat_i,lon_i]

        # After June (index 19)
        
        # Land change adjustment    
        # find indices where have deforestation
        inds_DF = np.asarray(np.where((x_change_v<1 ))).T
        counter_no_RF = 0
        total_areas = np.zeros(73)
    
        # loop over these indices, correct to savanna
        for i in range(len(inds_DF)):
            lat_i = inds_DF[i,0]
            lon_i = inds_DF[i,1]
            # check if actually have rainforest in this grid cell, otherwise do nothing
            total_RF_area = sum(Olson_landtype_A[IOLSON_6,lat_i,lon_i])
            if total_RF_area > 0:
                change_factor = x_change_v[lat_i,lon_i]
                # decrease rainforest landtypes by change factor
                Olson_landtype_def[IOLSON_6,lat_i,lon_i] = Olson_landtype_A[IOLSON_6,lat_i,lon_i] * change_factor
                # calculate total rainforest area lost
                RF_lost = total_RF_area * (1. - change_factor)
                # increase savanna landtypes by area lost
                Olson_landtype_def[savanna_id,lat_i,lon_i] = Olson_landtype_A[savanna_id,lat_i,lon_i] + RF_lost
                # alter the LAI of the rainforest land types to adjust for area scaling of LAI
                # index 19 is for start of June
                XLAI_w_def[IOLSON_6,19:,lat_i,lon_i] = XLAI_w_def[IOLSON_6,19:,lat_i,lon_i] * change_factor
                # replace lost forest by savanna type vegetation
                # index 19 is for start of June
                XLAI_w_def[savanna_id, 19:, lat_i,lon_i] = XLAI_avg[19:] * Olson_landtype_def[savanna_id,lat_i,lon_i]
            else:
                # calculate what is mean landtype in these grid cells
                total_areas = total_areas + Olson_landtype_A[:,lat_i,lon_i]
                counter_no_RF += 1
                
        # make previous year results accessible from the next year loop
        x_change_v_prev = x_change_v
        Olson_landtype_def_prev = Olson_landtype_def
        
        # Replace within global matrix
        # create copies of variables
        Olson_landtype_c = Olson_landtype.copy() 
        XLAI_w_c = XLAI_w.copy()
        
        Olson_landtype_c[:,(lat<=lat_max)&(lat>=lat_min),(lon<=lon_max)&(lon>=lon_min)] = Olson_landtype_def
    
        XLAI_w_c[:,:, (lat<=lat_max)&(lat>=lat_min), (lon<=lon_max)&(lon>=lon_min)] = XLAI_w_def
        
        os.system('cp ' + ifn_LT + ' ' + ofn_LT) # copy files to new file names
        os.system('cp ' + ifn_LAI + ' ' + ofn_LAI) # copy files to new file names

        dset = netCDF4.Dataset(ofn_LT, 'r+')
        dset2 = netCDF4.Dataset(ofn_LAI, 'r+')
        
        # loop through, replace rainforest land type
        for i in IOLSON_6:
            varname = 'LANDTYPE' + str(i) # variable name in netcdf
            dset[varname][:] = Olson_landtype_c[i,:,:][:]
            varname2 = 'XLAI' + str(i) # variable name in netcdf
            dset2[varname2][:] = XLAI_w_c[i,:,:,:][:]
        
        # replace savanna land type
        varname = 'LANDTYPE' + str(savanna_id) # variable name in netcdf
        dset[varname][:] = Olson_landtype_c[savanna_id,:,:][:]
        dset.close()
        
        varname2 = 'XLAI' + str(savanna_id) # variable name in netcdf
        dset2[varname2][:] = XLAI_w_c[savanna_id,:,:,:][:]
        dset2.close()
        
        # edit metadata so that time is correct for the desired year
        # units
        atr = "minutes since " + yr_s + "-01-01 00:00:00"
        cmd = 'ncatted -h -O -a  units,time,m,c,"' + atr +'" ' + ofn_LAI
        os.system(cmd)
        # begin_date
        atr = yr_s + "0101"
        cmd = 'ncatted -h -O -a  begin_date,time,m,c,"' + atr +'" ' + ofn_LAI
        os.system(cmd)
#%% Run a new loop to create Olson land map symbolic links
cmd = 'rm -r misc_Data/' + sce + '/Olson_LT/'
os.system(cmd)

cmd = 'mkdir misc_Data/' + sce + '/Olson_LT/'
os.system(cmd)

# Creating monthly symbolic links to correct land map
for i, yr in enumerate(year_r): # looping over years to create input files
    yr_s = str(yr) # year as string
    yr_s1 = str(yr - 1) # year minus 1 as string
    for k in range(1, 13): # looping over months
        mt = f'{k:02d}' # add leading zeros
        if k<6 and yr>2002: # use previous year values
            # Inputt file names - annual Olson land type
            ifn_LT = '/net/fs03/d0/arifein/python/misc_Data/' + sce + '/LAI/Olson_Land_Type_Masks.025x025.' + sce + '_' + yr_s1 +'.nc'
        else:
            # Input file names - annual Olson land type
            ifn_LT = '/net/fs03/d0/arifein/python/misc_Data/' + sce + '/LAI/Olson_Land_Type_Masks.025x025.' + sce + '_' + yr_s +'.nc'
            
        ofn_LT = 'misc_Data/' + sce + '/Olson_LT/Olson_Land_Type_Masks.025x025.' + sce + '_' + yr_s + mt + '.nc'
                
        cmd = "ln -s " + ifn_LT + " " + ofn_LT
        os.system(cmd)

