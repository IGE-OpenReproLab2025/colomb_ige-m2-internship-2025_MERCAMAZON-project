#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu April 14 2022
Creating Biomass burning emission files for Soares-Filho06 scenarios
The .tif data for Soares-Filho06 can be downloaded here: https://doi.org/10.3334/ORNLDAAC/1153
GEOS-Chem input data (i.e., biomass burning emissions from GFED) can be downloaded at: http://geoschemdata.wustl.edu/ExtData/
@author: arifeinberg
"""
import os
import xarray as xr
import numpy as np
import netCDF4
import rioxarray
from cdo import *
cdo = Cdo()
from datetime import date

def weird_division1(n, d):
    return n / d if d else 1 # no change if denominator 0
# create ufunc from weird division1
udiv1 = np.frompyfunc(weird_division1, 2, 1)
#%% Setup loop
startyear = 2002
endyear = 2050
year_r = np.arange(startyear, endyear+1)

sce = 'GOV' # scenario: 'BAU' or 'GOV'
sce_l = sce.lower() # scenario: 'BAU' or 'GOV'

DEFO_Hg_rate = 3.8e-7 # kg Hg/m^2 over year
defo_scale = 5.85E-08 # scaling factor kg Hg/kg DM

time_denom = (30 + 31) * 24 * 60 * 60 # 61 days in Aug and Sept, converted to s 

DEFO_Hg_u = DEFO_Hg_rate / defo_scale / time_denom # kg DM/m^2(deforested)/s

# Slow soil emissions from burning
slow_Hg_rate = DEFO_Hg_rate * 0.50 # kg Hg/m^2 over year, 50% extra emissions
s_yr = 365.2425*24*60*60 # yr to s
slow_Hg_u = slow_Hg_rate / defo_scale / s_yr
#%% Load necessary variables for biomass burning emissions
# load gbox areas for South America
fn_gb = '../gbox_areas/gbox_025x025_Olson.nc'
ds_gb = xr.open_dataset(fn_gb)
gbox_area = ds_gb.cell_area # grid_area
lat = ds_gb.lat # lat
lon = ds_gb.lon # lon


# restrict variables to area of South America
# select Amazon coordinate bounds
lat_min = -19.375
lat_max = 8.375
lon_min = -79.625
lon_max = -41.875

gbox_area_A = gbox_area[(lat<=lat_max) & (lat>=lat_min),(lon<=lon_max) & (lon>=lon_min)].values

#%% Run loop
for i, yr in enumerate(year_r): # looping over years to create input files
    
    yr_s = str(yr) # year as string
    yr_s1 = str(yr - 1) # year minus 1 as string
    print(yr_s)
    if (yr > 2002): # for year 2002, keep at reference BB emissions
        
        # Load GEOTiff files to calculate deforested area in current year
        # past year
        fn1 = '../Amazon_scenarios/' + sce + '_amazonia/amazon' + sce_l + yr_s1 +'.tif'
        ds1 = xr.open_dataset(fn1, engine="rasterio")
        landtypes1 = ds1.band_data.squeeze()
        # Before remapping, create new variable that assigns area fraction for forest type
        x_forest1 = xr.where(landtypes1==2,1.0,0.0) # intact forest
    
        # current year
        fn = '../Amazon_scenarios/' + sce + '_amazonia/amazon' + sce_l + yr_s +'.tif'
        ds = xr.open_dataset(fn, engine="rasterio")
        landtypes = ds.band_data.squeeze()
        # Before remapping, create new variable that assigns area fraction for forest type
        x_forest = xr.where(landtypes==2,1.0,0.0) # intact forest
        
        # Calculate where deforested pixels in past year are
        x_deforest = x_forest1 - x_forest
        # save dataArray as dataset
        ds_new = x_deforest.to_dataset(name='x_deforest')
        # clean up data array so that readable by cdo
        ds_new = ds_new.reset_coords(['band', 'spatial_ref'] ,drop=True)
        ds_new = ds_new.rename({'x': 'lon','y': 'lat'})
        ds_new = ds_new.reindex(lat=ds_new.lat[::-1]) # reverse lat dimension so that increasing
        # add attributes for lat/lon
        ds_new['lon'].attrs['units']="degrees east"
        ds_new['lat'].attrs['units']="degrees north"
        # Do regridding on 0.25 x 0.25 degree grid to compare with BB emissions
        ds_remap = cdo.remapcon("misc_Data/geos.025x025_SA.grid",input=ds_new, returnXDataset=True)
        # Get variables needed
        x_deforest_rm = ds_remap.x_deforest.values

        # x_deforest_rm is fraction of grid box that is deforested in that year
        # multiplying it by DEFO_Hg_u gives additional biomass burning emissions
        BB_deforest = x_deforest_rm * DEFO_Hg_u # immediate burning emissions 
        BB_slow = x_deforest_rm * slow_Hg_u # slow burning emissions from soils
        
    for k in range(1, 13): # looping over months
        mt = f'{k:02d}' # add leading zeros
        
        # Input file names - BB 
        ifn_BB = '../emissions/GFED4/v2020-02/2003/GFED4_gen.025x025.2003' + mt + '.nc'
        
        # make year directory if doesn't already exist
        os.system('mkdir -p misc_Data/' + sce + '/BB_slow/' + yr_s +'/')
        
        # Output file names - BB
        ofn_BB = 'misc_Data/' + sce + '/BB_slow/' + yr_s + '/GFED4_gen.025x025.' + sce + '_' + yr_s + mt + '.nc'
        
        # Copy reference 2003 emissions
        os.system('cp ' + ifn_BB + ' ' + ofn_BB) # copy files to new file names
 
        # open dataset for editing
        dset = netCDF4.Dataset(ofn_BB, 'r+')
        
        # Only add additional deforestation emissions after 2002
        if ((k == 8) or (k == 9)) and (yr > 2002): # August and September, burning events
            # load old deforestation emiss for 2003
            DEFO_old = xr.open_dataset(ofn_BB).DM_DEFO
            # create new version
            DEFO_new = DEFO_old.copy()
            # use new deforestation emissions for new values
            DEFO_new[0,(lat<=lat_max)&(lat>=lat_min),(lon<=lon_max)&(lon>=lon_min)] =\
                BB_deforest + BB_slow
            
            dset['DM_DEFO'][:] = DEFO_new
            
        if ((k < 8 )) and (yr > 2003): # add slow emissions from previous year
            # load old deforestation emiss for 2003
            DEFO_old = xr.open_dataset(ofn_BB).DM_DEFO
            # create new version
            DEFO_new = DEFO_old.copy()
            # use new deforestation emissions for new values
            DEFO_new[0,(lat<=lat_max)&(lat>=lat_min),(lon<=lon_max)&(lon>=lon_min)] +=\
                BB_slow_prev
            
            dset['DM_DEFO'][:] = DEFO_new
        
        if ((k > 9 )) and (yr > 2002): # add slow emissions from this year
            # load old deforestation emiss for 2003
            DEFO_old = xr.open_dataset(ofn_BB).DM_DEFO
            # create new version
            DEFO_new = DEFO_old.copy()
            # use new deforestation emissions for new values
            DEFO_new[0,(lat<=lat_max)&(lat>=lat_min),(lon<=lon_max)&(lon>=lon_min)] +=\
                BB_slow
            
            dset['DM_DEFO'][:] = DEFO_new
        
            
        # edit time so that correct for the desired year
        c_date = date(1985,1,1) # ref date
        n_date = date(yr, k, 1) # new date
        t_diff = (n_date - c_date).days * 24 # difference in hours
        
        n_time = t_diff # new time to edit in
        
        dset['time'][:] = n_time
        dset.close()

    if (yr > 2002): # save slow emissions for the next year
        BB_slow_prev = BB_slow 
        
