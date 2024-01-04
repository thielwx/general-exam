#!/usr/bin/env python
# coding: utf-8

# In[115]:


import pandas as pd
import numpy as np
import netCDF4 as nc
import sys
from datetime import datetime
from datetime import timedelta
from glob import glob
import os
import satpy.modifiers.parallax as plax
from pyresample import SwathDefinition, kd_tree

import matplotlib.pyplot as plt


# In[31]:


def time_string(time):
    year = time.strftime('%Y')
    doy = time.strftime("%j") #day of year
    month = time.strftime("%m")
    day = time.strftime("%d") #day of month
    
    hr = time.strftime("%H")
    minute = time.strftime("%M")
    
    return year, doy, month, day, hr, minute


# In[121]:


def MRMS_grabber(file_date, file_loc, hr, minute):
    
    #Ensuring we get the right tenths place for the minute (only grabbing on file)
    tens_place = str(int(int(minute)/10))
    #If we're on the 0 minute, use the MRMS file on the 0 minute
    if (int(minute)%10)==0:
        file = glob(file_loc+file_date+'-'+hr+tens_place+'0*.netcdf')
    
    #If we're on the 5 minute, use the MRMS file on the 0 minute
    elif (int(minute)%10)==5:
        file = glob(file_loc+file_date+'-'+hr+tens_place+'4*.netcdf')
    
    else:
        file = ['NOFILE']
        print ('NO MRMS FILE FOUND')
    
    return file


# In[118]:


def MRMS_data_formatter(file, var):
    #loading in the data from the MRMS netcdf file
    dset = nc.Dataset(file,'r')
    x_pix = dset.variables['pixel_x'][:] #Pixel locations (indicies) for LATITUDE
    y_pix = dset.variables['pixel_y'][:] #Pixel locations (indicies) for LONGITUDE
    data = dset.variables[var][:]
    
    u_lat = dset.Latitude #Upper-most latitude
    l_lon = dset.Longitude #Left-most longitude
    
    #Creating the arrays for the lat and lon coordinates
    y = dset.dimensions['Lat'].size #3500
    x = dset.dimensions['Lon'].size #7000
    lat = np.arange(u_lat, u_lat-(y*0.01),-0.01) #Going from upper to lower
    lon = np.arange(l_lon, l_lon+(x*0.01),0.01) #Going from left to right
    
    #Using the pixel indicides to get the pixel latitudes and longitudes
    lat_data = lat[x_pix] #Remember x_pixel represents LATITUDE
    lon_data = lon[y_pix] #Remember y_pixel represent LONGITUDE
    
    #Removing all data west of 103W and also any false data (less data to process)
    locs = np.where((data>0))[0]
    lon_data = lon_data[locs]
    lat_data = lat_data[locs]
    data = data[locs]
    
    #Creating the lon/lat grid
    lon_grid, lat_grid = np.meshgrid(lon, lat)
    
    #Defining the swaths of our listed and gridded data lat/lons
    MRMS_grid_swath = SwathDefinition(lons=lon_grid, lats=lat_grid)
    MRMS_point_swath = SwathDefinition(lons=lon_data, lats=lat_data)
    
    #Putting the data into a grid
    output = kd_tree.resample_nearest(source_geo_def=MRMS_point_swath,
                            data=data,
                            target_geo_def=MRMS_grid_swath,
                            radius_of_influence=1e2)
    
    
    return lat_grid, lon_grid, output, lat, lon


# In[141]:


def MRMS_max_finder(sel_OTs, data, lat_grid, lon_grid, lat_list, lon_list):
    
    OT_time = np.array(sel_OTs['time'])
    print (OT_time)
    
    OT_lats = sel_OTs['lat'].values
    OT_lons = sel_OTs['lon'].values
    
    #Going to just use geopotential height for geometric height becasue they are very similar in the troposphere
    #and lower stratosphere
    OT_zTrop = np.array(sel_OTs['era5_ztrop1'].values) * 1000
    
    #Using parallax correction to get more accurate values
    lon_search, lat_search = plax.get_parallax_corrected_lonlats(sat_lon=-75.0, sat_lat=0.0, sat_alt=35786023.0,
                                            lon=OT_lons, lat=OT_lats, height=OT_zTrop)
    
    for i in range(sel_OTs.shape[0]):
        print (i)
        #Creating a boolean grid of truth values to pick from
        selects = (lat_grid>=lat_search[i]-0.05) & (lat_grid<=lat_search[i]+0.05) & (lon_grid>=lon_search[i]-0.05) & (lon_grid<=lon_search[i]+0.05)
        
        select_data = data.copy()
        #Getting rid of the data outside outside of the area of interest
        select_data[~selects] = np.nan
        
        #Getting the maximum of what's left
        var_max = np.nanmax(select_data)
        
        print (var_max)
        
        OTs.loc[((OTs['lat']==OT_lats[i]) & (OTs['lon']==OT_lons[i]) & (OTs['time']==OT_time[i])), 'MergedReflectivityQCComposite'] = var_max
        
        


# In[142]:


#Reading in the processed OT file to use
OTs = pd.read_pickle('processed-OTs/OT-dataset-TTrop1_diff_0-cmi13_anvil_diff_6point5.pkl')

#Adding a new columns for our MRMS data
OTs['MergedReflectivityQCComposite'] = np.full(OTs.shape[0],np.nan)

case_list = ['20180514-Derecho', '20180628-Derecho', '20190719-Derecho', '20180514-MCS', '20180624-MCS', '20190717-MCS']

for case in case_list:
    
    mrms_loc = '/localdata/General-Exam/cases/'+case+'/MergedReflectivityQCComposite/'  

    #Setting the variables on a case by case basis
    if case=='20180514-Derecho':
        start_time = datetime(2018,5,14,15,40,0)
        #start_time = datetime(2018,5,14,20,0,0) #For local testing/debugging
        end_time = datetime(2018,5,15,1,0,0)
        #end_time = datetime(2018,5,14,20,20,0) #For local testing/debugging
        #mrms_loc = '../Test-Data/MergedReflectivityQCComposite/' #For local testing/debugging

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



    #Creating the list of times to pull from
    time_list = pd.date_range(start=start_time, end=end_time, freq='5T').to_list()
    
    
    for i in range(len(time_list)):
        #Selecting OTs from the current time
        sel_OTs = OTs.loc[OTs['time']==time_list[i]]

        #Checking that we have OTs that need MRMS data
        if sel_OTs.shape[0]>0:
            print (time_list[i])

            #Getting the time
            year, doy, month, day, hr, minute = time_string(time_list[i])
            file_date = year + month + day

            #Grabbing the file name of the MRMS data
            comp_dbz_file = MRMS_grabber(file_date, mrms_loc, hr, minute)
            
	    if comp_dbz_file[0]=='NOFILE':
	    	continue
		

            #Getting the MRMS data
            lat_grid, lon_grid, comp_dbz_grid, lat_list, lon_list = MRMS_data_formatter(comp_dbz_file[0], 'MergedReflectivityQCComposite')

            MRMS_max_finder(sel_OTs, comp_dbz_grid, lat_grid, lon_grid, lat_list, lon_list)


        else:
            continue
            

    OTs.to_pickle('processed-OTs/OT-dataset-TTrop1_diff_0-cmi13_anvil_diff_6point5-MRMSadded.pkl')

