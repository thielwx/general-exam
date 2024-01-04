#!/usr/bin/env python
# coding: utf-8

# A collections of functions that set the bounds for the MCSs and Derechos

# In[1]:


from datetime import datetime
import numpy as np


# In[ ]:


def bounds_20180514_Derecho(time, diff, sat_lons, sat_lats):
    if (time>=datetime(2018,5,14,15,40,0)) & (time<=datetime(2018,5,15,1,0,0)):
        latmax = 42.00
        latmin = 36.50
        lonmax = -75.00
        lonmin = -85.00
    else:
        print ('Outside of range, no bounds applied')
        latmax = 90
        latmin = 0
        lonmax = 0
        lonmin = -180
    
    #Finding the points that are in our time-relative bounds
    picker = (sat_lons>lonmin) & (sat_lons<lonmax) & (sat_lats>latmin) & (sat_lats<latmax)
    diff[~picker] = np.nan #Setting points outside those boudns to nan
    
    return diff


# In[ ]:


def bounds_20180628_Derecho(time, diff, sat_lons, sat_lats):
    if (time>=datetime(2018, 6, 28, 22, 40, 0)) & (time<datetime(2018, 6, 29, 2, 30, 0)):
        latmax = 52.00
        latmin = 45.00
        lonmax = -103.00
        lonmin = -109.00
    elif (time>=datetime(2018, 6, 29, 2, 30, 0)) & (time<datetime(2018, 6, 29, 9, 30, 0)):
        latmax = 52.00
        latmin = 45.00
        lonmax = -95.00
        lonmin = -109.00    
    elif (time>=datetime(2018, 6, 28, 9, 30, 0)) & (time<datetime(2018, 6, 29, 13, 0, 0)):
        latmax = 52.00
        latmin = 45.00
        lonmax = -85.00
        lonmin = -102.00
    elif (time>=datetime(2018, 6, 29, 13, 0, 0)) & (time<=datetime(2018, 6, 29, 16, 55, 0)):
        latmax = 52.00
        latmin = 45.00
        lonmax = -85.00
        lonmin = -99.00
    else:
        print ('Outside of range, no bounds applied')
        latmax = 90
        latmin = 0
        lonmax = 0
        lonmin = -180
    
    #Finding the points that are in our time-relative bounds
    picker = (sat_lons>lonmin) & (sat_lons<lonmax) & (sat_lats>latmin) & (sat_lats<latmax)
    diff[~picker] = np.nan #Setting points outside those boudns to nan
    return diff


# In[ ]:


def bounds_20190719_Derecho(time, diff, sat_lons, sat_lats):
    if (time>=datetime(2019,7,19,21,10,0)) & (time<datetime(2019,7,20,2,30,0)):
        latmax = 48.00
        latmin = 43.50
        lonmax = -86.00
        lonmin = -97.00
    elif (time>=datetime(2019,7,20,2,30,0)) & (time<datetime(2019,7,20,6,0,0)):
        latmax = 48.00
        latmin = 42.00
        lonmax = -83.00
        lonmin = -92.00
    elif (time>=datetime(2019,7,20,6,0,0)) & (time<=datetime(2019,7,20,7,50,0)):
        latmax = 48.00
        latmin = 41.00
        lonmax = -83.00
        lonmin = -90.60
    else:
        print ('Outside of range, no bounds applied')
        latmax = 90
        latmin = 0
        lonmax = 0
        lonmin = -180
    
    #Finding the points that are in our time-relative bounds
    picker = (sat_lons>lonmin) & (sat_lons<lonmax) & (sat_lats>latmin) & (sat_lats<latmax)
    diff[~picker] = np.nan #Setting points outside those boudns to nan
    
    return diff


# In[ ]:


def bounds_20180514_MCS(time, diff, sat_lons, sat_lats):
    if (time>=datetime(2018,5,14,20,5,0)) & (time<datetime(2018,5,14,21,30,0)):
        latmax = 40.50
        latmin = 38.00
        lonmax = -88.00
        lonmin = -94.00
    elif (time>=datetime(2018,5,14,21,30,0)) & (time<datetime(2018,5,14,23,0,0)):
        latmax = 40.50
        latmin = 38.00
        lonmax = -88.00
        lonmin = -93.00
    elif (time>=datetime(2018,5,14,23,0,0)) & (time<=datetime(2018,5,15,0,35,0)):
        latmax = 40.50
        latmin = 38.00
        lonmax = -88.00
        lonmin = -92.00
    else:
        print ('Outside of range, no bounds applied')
        latmax = 90
        latmin = 0
        lonmax = 0
        lonmin = -180
    
    #Finding the points that are in our time-relative bounds
    picker = (sat_lons>lonmin) & (sat_lons<lonmax) & (sat_lats>latmin) & (sat_lats<latmax)
    diff[~picker] = np.nan #Setting points outside those boudns to nan
    
    return diff


# In[ ]:


def bounds_20180624_MCS(time, diff, sat_lons, sat_lats):
    if (time>=datetime(2018,6,23,21,25,0)) & (time<datetime(2018,6,23,23,30,0)):
        latmax = 40.50
        latmin = 38.00
        lonmax = -101.00
        lonmin = -106.00
    elif (time>=datetime(2018,6,23,23,30,0)) & (time<datetime(2018,6,24,4,0,0)):
        latmax = 39.50
        latmin = 36.00
        lonmax = -99.00
        lonmin = -106.00
    elif (time>=datetime(2018,6,24,4,0,0)) & (time<datetime(2018,6,24,7,0,0)):
        latmax = 39.00
        latmin = 35.50
        lonmax = -97.50
        lonmin = -102.00
    elif (time>=datetime(2018,6,24,7,0,0)) & (time<datetime(2018,6,24,8,0,0)):
        latmax = 38.00
        latmin = 34.50
        lonmax = -97.00
        lonmin = -101.00
    elif (time>=datetime(2018,6,24,8,0,0)) & (time<datetime(2018,6,24,9,0,0)):
        latmax = 38.00
        latmin = 34.50
        lonmax = -96.50
        lonmin = -100.00
    elif (time>=datetime(2018,6,24,9,0,0)) & (time<=datetime(2018,6,24,9,55,0)):
        latmax = 38.00
        latmin = 34.50
        lonmax = -96.00
        lonmin = -99.00
    else:
        print ('Outside of range, no bounds applied')
        latmax = 90
        latmin = 0
        lonmax = 0
        lonmin = -180
    
    #Finding the points that are in our time-relative bounds
    picker = (sat_lons>lonmin) & (sat_lons<lonmax) & (sat_lats>latmin) & (sat_lats<latmax)
    diff[~picker] = np.nan #Setting points outside those boudns to nan
    
    return diff


# In[ ]:


def bounds_20190717_MCS(time, diff, sat_lons, sat_lats):
    if (time>=datetime(2019,7,17,8,15,0)) & (time<datetime(2019,7,17,13,0,0)):
        latmax = 45.00
        latmin = 41.50
        lonmax = -94.00
        lonmin = -102.00
    elif (time>=datetime(2019,7,17,13,0,0)) & (time<datetime(2019,7,17,14,30,0)):
        latmax = 44.00
        latmin = 41.00
        lonmax = -92.00
        lonmin = -98.00
    elif (time>=datetime(2019,7,17,14,30,0)) & (time<datetime(2019,7,17,17,0,0)):
        latmax = 44.00
        latmin = 41.00
        lonmax = -91.00
        lonmin = -96.50
    elif (time>=datetime(2019,7,17,17,0,0)) & (time<=datetime(2019,7,17,22,20,0)):
        latmax = 44.00
        latmin = 38.50
        lonmax = -90.00
        lonmin = -96.00
    else:
        print ('Outside of range, no bounds applied')
        latmax = 90
        latmin = 0
        lonmax = 0
        lonmin = -180
    
    #Finding the points that are in our time-relative bounds
    picker = (sat_lons>lonmin) & (sat_lons<lonmax) & (sat_lats>latmin) & (sat_lats<latmax)
    diff[~picker] = np.nan #Setting points outside those boudns to nan
    
    return diff

