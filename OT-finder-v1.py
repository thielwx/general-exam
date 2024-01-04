#!/usr/bin/env python
# coding: utf-8

# This script identifies potential OTs using ERA5 data along with GOES-16 ABI CMI(14) CONUS data. Tropopause relative IR brightness temperatures are calculated and peaks are identified, then saved into a pandas dataframe

# In[2]:


import netCDF4 as nc
import pandas as pd
import numpy as np
from datetime import datetime
from datetime import timedelta
from pyproj import Proj
from pyresample import SwathDefinition, kd_tree
from skimage.feature import peak_local_max

from glob import glob
import os
import sys

sys.path.insert(1, '/localdata/General-Exam/')
#sys.path.insert(1, '')
from case_bounds import *


# In[3]:


#========Super important variables==========GLOBAL========
tropopause_relative_IR_temp_threshold = 0 #Threshold between the ERA5 Ttrop1 and CMIP14 that
max_window_length = 6 # Translates to no OTs within 6 km of another (moving window is 12 km wide)

anvil_radius = 5 #Means we're selecting points in a 10 by 10 box surrounding the potential OT
anvil_temp_threshold = 225 #Only sampling points that we beleive are in the anvil to identify contrast with pOT points


# In[20]:


#Reading in arguments on a case by case basis
args = sys.argv
#args = [0, '20180514-Derecho'] #For local testing/debugging

case = args[1]
case_loc = '/localdata/General-Exam/cases/'+case+'/'

#Setting the variables on a case by case basis
if case=='20180514-Derecho':
    start_time = datetime(2018,5,14,15,40,0)
    #start_time = datetime(2018,5,14,20,0,0) #For local testing/debugging
    end_time = datetime(2018,5,15,1,0,0)
    #end_time = datetime(2018,5,14,22,0,0) #For local testing/debugging
    #case_loc = '../Test-Data/' #For local testing/debugging
    
elif case=='20180628-Derecho':
    start_time = datetime(2018,6,28,22,40,0)
    end_time = datetime(2018,6,29,16,55,0)
elif case=='20190719-Derecho':
    start_time = datetime(2019,7,19,21,10,0)
    end_time = datetime(2019,7,20,7,50,0)
elif case=='20180514-MCS':
    start_time = datetime(2018,5,14,20,5,0)
    end_time = datetime(2018,5,15,0,35,0)
elif case=='20180624-MCS':
    start_time = datetime(2018,6,23,21,25,0)
    end_time = datetime(2018,6,24,9,55,0)
elif case=='20190717-MCS':
    start_time = datetime(2019,7,17,8,15,0)
    end_time = datetime(2019,7,17,22,20,0)
else:
    print ('---ERROR: CASE UNDEFINED---')
    sys.exit()
    
    
    
#Creating the list of times to pull from
time_list = pd.date_range(start=start_time, end=end_time, freq='5T').to_list()


# In[21]:


def time_string(time):
    year = time.strftime('%Y')
    doy = time.strftime("%j") #day of year
    month = time.strftime("%m")
    day = time.strftime("%d") #day of month
    
    hr = time.strftime("%H")
    minute = time.strftime("%M")
    
    return year, doy, month, day, hr, minute


# # Function Land (Grabbing the data)

# In[22]:


def ERA_data_grabber(time,case_loc):
    #Getting the file name
    year, doy, month, day, hr, minute = time_string(time)
    era_file_str = year + month + day + '_trop.nc'
    
    era_data_loc = case_loc + 'ERA5-data/' + era_file_str #NOTE, move ERA5 data to correct positions and make directories in case file
    
    
    #Reading in the netCDF file
    era = nc.Dataset(era_data_loc,'r')
    
    #Pulling the ERA variables available
    era_lon = era.variables['Longitude'][:]
    era_lat = era.variables['Latitude'][:]

    Ttrop1 = era.variables['Ttrop1'][:,:,:]
    Ttrop2 = era.variables['Ttrop2'][:,:,:]
    
    ztrop1 = era.variables['ztrop1'][:,:,:]
    ztrop2 = era.variables['ztrop2'][:,:,:]
    
    ptrop1 = era.variables['ptrop1'][:,:,:]
    ptrop2 = era.variables['ptrop2'][:,:,:]
    
    return era_lon, era_lat, Ttrop1, Ttrop2, ztrop1, ztrop2, ptrop1, ptrop2


# In[23]:


#Takes in file locations and names and returns any files found from the previous five mintues
#Adapted from https://github.com/thielwx/composites/blob/master/make_df.py

def fivechecker(floc,flist):
    
    files = ()
    #Checking the files over the previous five mintues
    for i in range(len(flist)):
        a = glob(floc+flist[i])
        files = np.append(files,a)
    return files


# In[24]:


#A function that takes in a netCDF Dataset and returns its x/y coordinates as as lon/lat
def latlon(data):    
    sat_h = data.variables['goes_imager_projection'].perspective_point_height
    sat_lon = data.variables['goes_imager_projection'].longitude_of_projection_origin
    sat_sweep = data.variables['goes_imager_projection'].sweep_angle_axis
    
    X = data.variables['x'][:] * sat_h
    Y = data.variables['y'][:] * sat_h

    p = Proj(proj='geos', h=sat_h, lon_0=sat_lon, sweep=sat_sweep)
    YY, XX = np.meshgrid(Y, X)
    lons, lats = p(XX, YY, inverse=True)
    lons[lons > 10000] = np.nan
    lats[lats > 10000] = np.nan
    
    return lons.T, lats.T


# In[25]:


def cmi_data_grabber(time, case_loc):
    dt = timedelta(seconds = 60)
    
    cmi_str = 'OR_ABI-L2-CMIPC-M?C??_G16_s'
    
    #Creating a list of files over the previous 5 minutes, so we can search for them later
    cmi_list = np.empty(0)
    for i in range(5):
        cur_time = time - timedelta(seconds = 60*i) #Going backwards in time from the current time
        year, doy, month, day, hr, minute = time_string(cur_time)
        cmi_time_str = year+doy+hr+minute+'*.nc'
        
        cmi_list = np.append(cmi_list,cmi_str+cmi_time_str)
    
    cmi08_base = case_loc+'CMIPC08/'
    cmi13_base = case_loc+'CMIPC13/'
    cmi14_base = case_loc+'CMIPC14/'
    
    #Searching for files
    cmi08_files = fivechecker(cmi08_base,cmi_list)
    cmi13_files = fivechecker(cmi13_base,cmi_list)
    cmi14_files = fivechecker(cmi14_base,cmi_list)
    
    if len(cmi08_files)==1:
        cmi08_file = nc.Dataset(cmi08_files[0])
        cmi08 = cmi08_file.variables['CMI'][:,:]
    else:
        cmi08=[-999]
    
    if len(cmi13_files)==1:
        cmi13_file = nc.Dataset(cmi13_files[0])
        cmi13 = cmi13_file.variables['CMI'][:,:]
    else:
        cmi13=[-999]
        
    if len(cmi14_files)==1:
        cmi14_file = nc.Dataset(cmi14_files[0])
        cmi14 = cmi14_file.variables['CMI'][:,:]
        sat_lons, sat_lats = latlon(cmi14_file)
        proj = cmi14_file.variables['goes_imager_projection']
    else:
        cmi14=[-999]
        sat_lons=-999
        sat_lats=-999
        proj = -999
        
    return cmi08, cmi13, cmi14, sat_lons, sat_lats, proj


# # Funciton Land (Processing the data)

# In[26]:


def neg_lons(lon):
    array = np.arange(0,len(lon),1)
    for i in array:
        if lon[i]>180:
            lon[i] = lon[i]-360
    return lon


# In[27]:


'''
Modified from https://github.com/thielwx/utilities/blob/master/fn_master.py
'''

def resample(ERA_field, ERA_lats, ERA_lons, target_lats, target_lons):

    #Adjusting the negative longitudes on the ERA5 Grid
    #Creating 2D grids of the ERA data
    ERA_lons = neg_lons(ERA_lons)
    ERA_lats2d, ERA_lons2d = np.meshgrid(ERA_lats,ERA_lons)

    #Creating the swath definitions
    target_swath = SwathDefinition(lons=target_lons,lats=target_lats)
    ERA_swath = SwathDefinition(lons=ERA_lons2d.T,lats=ERA_lats2d.T)

    #Resampling using a KD-tree to fit the data to a grid
    output = kd_tree.resample_nearest(source_geo_def=ERA_swath,
                            data=ERA_field,
                            target_geo_def=target_swath,
                            radius_of_influence=4e4)
    
    return output


# In[28]:


def max_finder(diff, sat_lats, sat_lons):
    #Finding the local maxima in relative temps
    xy = peak_local_max(diff, min_distance=max_window_length, threshold_abs=tropopause_relative_IR_temp_threshold)
    
    #Caputring the lats and lons of the max points
    max_lats = []
    max_lons = []
    for i in range(xy.shape[0]):
        row_int = xy[i,0]
        col_int = xy[i,1]
        max_lats = np.append(max_lats, sat_lats[row_int,col_int])
        max_lons = np.append(max_lons, sat_lons[row_int,col_int])
        
    return xy, max_lats, max_lons


# In[29]:


#Working backward to select the points that use the satellite scene
def sat_data_finder(var,xy):
    var_max_pt = []
    
    #Looping through each of the max lat points
    for i in range(xy.shape[0]):
        #Getting the integers for each row and column
        row_int = xy[i,0]
        col_int = xy[i,1]
        
        #Finding the variable selected from each row/column
        var_selected = var[row_int, col_int]
            
        #Adding the selected variable to the array that is returned
        var_max_pt = np.append(var_max_pt,var_selected)
        
    return var_max_pt


# In[30]:


#How to pull in data from around the OT

def anvil_finder(var, xy):
    di = anvil_radius #Number of satellite data points to pull from the pOT point
    di_inner = 3 #Number of satellite data points at the center to exclude (removing OT area from calculation)
    
    #Empty arrays for each percentile class you want to have
    var_p25 = []
    var_p50 = []
    var_p75 = []
    var_mean = []
    
    #Ensuring that the anvil pixels are brightness temperatures less than 225 K (Bedka et al 2010, Dworak et al 2012)
    var[var>anvil_temp_threshold] = np.nan
    
    for i in range(xy.shape[0]):
        var_temp = var.copy() #Setting at temporary copy so we can remove the center pixels from each pOT point
        row_int = xy[i,0]
        col_int = xy[i,1]
        
        #Removing areas directly around the OT
        var_temp[row_int-di_inner:row_int+di_inner,col_int-di_inner:col_int+di_inner] = np.nan
        #Selecting only the points within the dataset 
        var_cut = var_temp[row_int-di:row_int+di, col_int-di:col_int+di]
        
        var_p25 = np.append(var_p25, np.nanpercentile(var_cut, 25))
        var_p50 = np.append(var_p50, np.nanpercentile(var_cut, 50))
        var_p75 = np.append(var_p75, np.nanpercentile(var_cut, 75))
        var_mean = np.append(var_mean, np.nanmean(var_cut))
        
    return var_p25, var_p50, var_p75, var_mean


# # The funciton that drives it all

# In[36]:


#Trying to make it so we don't have to pull new ERA data every time (for speed)
era_check_day = (time_list[0]- timedelta(days=1)).strftime("%d")  
era_check_hr = (time_list[0] - timedelta(hours=1)).strftime("%H") 

output_df = pd.DataFrame()

for i in range(len(time_list)):
    print (time_list[i])
    
    #Getting strings of the current time
    year, doy, month, day, hr, minute = time_string(time_list[i])
    
    #Checking if the day has changed and if we need to pull from a different ERA file
    if era_check_day != day:
        #Grabbing the ERA data for the current day
        era_lon, era_lat, Ttrop1, Ttrop2, ztrop1, ztrop2, ptrop1, ptrop2 = ERA_data_grabber(time_list[i],case_loc)
        #Updating the check day so it doesn't hit it again until the next hour
        era_check_day = day
    print ('ERA5 data checked')
    
    #Grabbing the GOES-16 CONUS, ABI data from bands 13 and 14
    cmi08, cmi13, cmi14, sat_lons, sat_lats, proj = cmi_data_grabber(time_list[i],case_loc)
    
    #Checking that we have the CMI data (no dropouts)
    if (len(cmi08)==1) or (len(cmi13)==1) or (len(cmi14)==1):
        print ('***CMI DATA MISSING***')
        continue
    else:
        print('CMI data checked')
    
    #Resampling the ERA data to match the GOES-16 CMI data
    if era_check_hr != hr:
        Ttrop1_res = resample(Ttrop1[int(hr),:,:], era_lat, era_lon, sat_lats, sat_lons)
        Ttrop2_res = resample(Ttrop2[int(hr),:,:], era_lat, era_lon, sat_lats, sat_lons)
        ptrop1_res = resample(ptrop1[int(hr),:,:], era_lat, era_lon, sat_lats, sat_lons)
        ptrop2_res = resample(ptrop2[int(hr),:,:], era_lat, era_lon, sat_lats, sat_lons)
        ztrop1_res = resample(ztrop1[int(hr),:,:], era_lat, era_lon, sat_lats, sat_lons)
        ztrop2_res = resample(ztrop2[int(hr),:,:], era_lat, era_lon, sat_lats, sat_lons)
        era_check_hr = hr
    print ('ERA5 data resampled')
    #Calculating the tropopause relative IR temps
    #Note: IR T_B colder than the tropopause are POSITIVE
    diff = Ttrop1_res-cmi14    
    
    #Creating the bounds for our pOT ID tests
    if case=='20180514-Derecho':
        diff = bounds_20180514_Derecho(time_list[i], diff, sat_lons, sat_lats)
    elif case=='20180628-Derecho':
        diff = bounds_20180628_Derecho(time_list[i], diff, sat_lons, sat_lats)
    elif case=='20190719-Derecho':
        diff = bounds_20190719_Derecho(time_list[i], diff, sat_lons, sat_lats)
    elif case=='20180514-MCS':
        diff = bounds_20180514_MCS(time_list[i], diff, sat_lons, sat_lats)
    elif case=='20180624-MCS':
        diff = bounds_20180624_MCS(time_list[i], diff, sat_lons, sat_lats)
    elif case=='20190717-MCS':
        diff = bounds_20190717_MCS(time_list[i], diff, sat_lons, sat_lats)
    
    #Finding the maximum points
    xy, max_lats, max_lons = max_finder(diff, sat_lats, sat_lons)
    print ('Potential OTs found:'+str(len(max_lats)))
    
    #Getting the ERA5 and CMI data from the selected points
    Ttrop1_final = sat_data_finder(Ttrop1_res, xy)
    Ttrop2_final = sat_data_finder(Ttrop2_res, xy)
    ztrop1_final = sat_data_finder(ztrop1_res, xy)
    ztrop2_final = sat_data_finder(ztrop2_res, xy)
    ptrop1_final = sat_data_finder(ptrop1_res, xy)
    ptrop2_final = sat_data_finder(ptrop2_res, xy)
    cmi08_final = sat_data_finder(cmi08,xy)
    cmi13_final = sat_data_finder(cmi13,xy)
    cmi14_final = sat_data_finder(cmi14,xy)
    
    
    cmi13_p25, cmi13_p50, cmi13_p75, cmi13_mean = anvil_finder(cmi13,xy)
    cmi14_p25, cmi14_p50, cmi14_p75, cmi14_mean = anvil_finder(cmi14,xy)
    
    #Preparing data to be put into a dataframe
    time_array = np.full(len(max_lats),time_list[i])
    
    data = {
        'time': time_array,
        'lat': max_lats,
        'lon': max_lons,
        'scene_row':xy[:,0],
        'scene_col':xy[:,1],
        'era5_Ttrop1':Ttrop1_final,
        'era5_Ttrop2':Ttrop2_final,
        'era5_ptrop1':ptrop1_final,
        'era5_ptrop2':ptrop2_final,
        'era5_ztrop1':ztrop1_final,
        'era5_ztrop2':ztrop2_final,
        'cmi08': cmi08_final,
        'cmi13': cmi13_final,
        'cmi14': cmi14_final,
        'cmi13_p25':cmi13_p25,
        'cmi13_p50':cmi13_p50,
        'cmi13_p75':cmi13_p75,
        'cmi13_mean':cmi13_mean,
        'cmi14_p25':cmi14_p25,
        'cmi14_p50':cmi14_p50,
        'cmi14_mean':cmi14_mean,
        'cmi14_p75':cmi14_p75,
    }
    
    #Creating dataframe and appending the output
    curr_df = pd.DataFrame(data=data)
    output_df = pd.concat((output_df,curr_df), axis=0)


# In[40]:


output_str = '/localdata/General-Exam/OTs/'+case+'-unfiltered-trop_rel_IR_temp_'+str(tropopause_relative_IR_temp_threshold)+'-v1.pkl'
print (output_str)
output_df.to_pickle(output_str)

