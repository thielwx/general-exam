#!/usr/bin/env python
# coding: utf-8

# This script pulls GOES data from the Google Cloud storage on an hourly basis

# In[1]:


import subprocess as sp
import pandas as pd
import numpy as np
import os
import sys
from datetime import timedelta
from datetime import datetime


# In[46]:


#Setting up the input arguments/constants

args = sys.argv
#args = [0, '20180514-Derecho', 'ABI', 'ACHAC', 17, '201805141200', '201805141600', 3]

case = args[1] #YYYYMMDD-Derecho (or MCS)
sensor = args[2] #ABI or GLM
product = args[3] #CMIPC, ACHAC
channel_num = args[4] #ABI channel number (if unnneeded, set to 17)
t_start = args[5] #Start time of the data pull (YYYYMMDDHHmm)
t_end = args[6] #End time of the data pull (YYYYMMDDHHmm)
meso_num = args[7] #Mesoscale domain number (if unneeded, set to 3)




# Constants
if int(channel_num)<17:
    save_loc = '/localdata/General-Exam/cases/'+case+'/'+product+str(channel_num).zfill(2)+'/'
else:
    save_loc = '/localdata/General-Exam/cases/'+case+'/'+product+'/'

print (save_loc)    
    
if not os.path.exists(save_loc):
    os.makedirs(save_loc) 


# In[47]:


t_start_dt = datetime.strptime(t_start, '%Y%m%d%H%M%S') #Converting into datetimes
t_end_dt = datetime.strptime(t_end, '%Y%m%d%H%M%S') #Converting into datetimes

tdelta = timedelta(minutes=60) #Creating a time delta, so the data doesn't run over

#Creating a list of dates and times that we can pull our hourly analysis from
time_list = pd.date_range(start=t_start_dt, end=t_end_dt-tdelta, freq='H').to_list()


# In[48]:


def command_maker(t, sensor, product, channel_num, meso_num, save_loc):
    year = t.strftime('%Y')
    doy = t.strftime("%j")
    hr = t.strftime("%H")
    #CONUS CMIP
    if (product=='CMIPC'):
        half2 = 'OR_'+sensor+'-L2-'+product+'-M?C'+str(channel_num).zfill(2)+'_G16_s'+year+doy+hr+'* '+save_loc
        command = 'gsutil -m cp gs://gcp-public-data-goes-16/'+sensor+'-L2-'+product+'/'+year+'/'+doy+'/'+hr+'/'+half2
    #Mesosector CMIP    
    elif (product=='CMIPM'):
        half2 = 'OR_'+sensor+'-L2-'+product+str(meso_num)+'-M?C'+str(channel_num).zfill(2)+'_G16_s'+year+doy+hr+'* '+save_loc
        command = 'gsutil -m cp gs://gcp-public-data-goes-16/'+sensor+'-L2-'+product+'/'+year+'/'+doy+'/'+hr+'/'+half2
    #Mesosector L2    
    elif (int(meso_num)!=3):
        half2 = 'OR_'+sensor+'-L2-'+product+str(meso_num)+'-M?_G16_s'+year+doy+hr+'* '+save_loc
        command = 'gsutil -m cp gs://gcp-public-data-goes-16/'+sensor+'-L2-'+product+'/'+year+'/'+doy+'/'+hr+'/'+half2
    #GLM    
    elif (product=='GLM'):
        command = 'gsutil -m cp gs://gcp-public-data-goes-16/GLM-L2-LCFA/'+year+'/'+doy+'/'+hr+'/OR_GLM-L2-LCFA_G16_s'+year+doy+hr+'* '+save_loc
    #CONUS L2    
    else: #All L2 CONUS products
        half2 = 'OR_'+sensor+'-L2-'+product+'-M?_G16_s'+year+doy+hr+'* '+save_loc
        command = 'gsutil -m cp gs://gcp-public-data-goes-16/'+sensor+'-L2-'+product+'/'+year+'/'+doy+'/'+hr+'/'+half2
        
    print (command)
    
    return command  


# In[50]:


for t in time_list:
    command = command_maker(t, sensor, product, channel_num, meso_num, save_loc)
    
    p = sp.Popen(command,shell=True)
    p.wait()


# In[ ]:




