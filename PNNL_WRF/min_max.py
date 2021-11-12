"""
Atlertaive way for findings LLJs. Using argrelextrema. DOESN'T WORK!!!! Appears close but not sure what is off
"""

from scipy.signal import argrelextrema
from scipy.ndimage import gaussian_filter

import wrf
from netCDF4 import Dataset
from matplotlib.cm import get_cmap
import cartopy.crs as crs
import cartopy
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
#import xarray

import os
import pandas as pd
import numpy as np
from numpy import exp,where,ma,cos,sin,pi,amax,amin,asarray
import datetime
from datetime import timedelta, date
from scipy import stats
import re
from scipy.signal import argrelextrema

file_in_dt_fmt = '%Y%m%d%H%M'

Rs_da=287.05          # Specific gas const for dry air, J kg^{-1} K^{-1}
Rs_v=461.51           # Specific gas const for water vapour, J kg^{-1} K^{-1}

## vapor pressure calculation function

def VaporPressure(tempc,phase="liquid"):
    """Water vapor pressure over liquid water or ice.

    INPUTS: 
    tempc: (C) OR dwpt (C), if SATURATION vapour pressure is desired.
    phase: ['liquid'],'ice'. If 'liquid', do simple dew point. If 'ice',
    return saturation vapour pressure as follows:

    Tc>=0: es = es_liquid
    Tc <0: es = es_ice

   
    RETURNS: e_sat  (Pa)
    
    SOURCE: http://cires.colorado.edu/~voemel/vp.html (#2:
    CIMO guide (WMO 2008), modified to return values in Pa)
    
    This formulation is chosen because of its appealing simplicity, 
    but it performs very well with respect to the reference forms
    at temperatures above -40 C. At some point I'll implement Goff-Gratch
    (from the same resource).
    """

    over_liquid=6.112*exp(17.67*tempc/(tempc+243.12))*100.
    over_ice=6.112*exp(22.46*tempc/(tempc+272.62))*100.
    # return where(tempc<0,over_ice,over_liquid)

    if phase=="liquid":
        # return 6.112*exp(17.67*tempc/(tempc+243.12))*100.
        return over_liquid
    elif phase=="ice":
        # return 6.112*exp(22.46*tempc/(tempc+272.62))*100.
        return where(tempc<0,over_ice,over_liquid)
    else:
        raise NotImplementedError


matplotlib.rcParams.update({'font.size': 25})


#change the start and end times to be plotted
x_start_day = 20181101
x_end_day = 20181218

path = '/home/disk/monsoon/relampago/qc/sounding/eol_composite/High_Res_v1.1'

#pressure to get
min_pres = 450
max_pres = 970

outpath = '/home/disk/meso-home/crs326/Documents/Research/PNNL_WRF/'

Station = 'CordobaAero'

found_time_list = []

total_count_obs = 0

#criteria
# 1: original criteria from Oliviera et al 2018
crit1 = [1, 3000, 4000, 23.326, 11.663]
# 2: high - allows for jets centered higher in the atmosphere, and a bit weaker compared to mean flow
crit2 = [2 ,3200, 5700, 23.326, 11.663]

crit = crit2

# ----------------------------------------------------------------------------------------

# get dates from times to only run for correct day folders
input_start_date = pd.to_datetime(str(x_start_day)[:8], format='%Y%m%d', errors='ignore')
input_end_date = pd.to_datetime(str(x_end_day)[:8], format='%Y%m%d', errors='ignore')


crit_num = crit[0]
max_search_hgt = crit[1]
min_search_hgt = crit[2]
max_wind_threshold = crit[3]
decrease_to_min = crit[4]
print crit

time_list = []
sounding_filepaths = []
list_folder_paths = [] # list that will contain folders (path+name)
total_count_obs = 0
N_max_wind_time_list = []
N_drct_mode_list = []
N_speed_mean_list = []
N_max_wind_pres_list = []
S_max_wind_time_list = []
S_drct_mode_list = []
S_speed_mean_list = []
S_max_wind_pres_list = []
EW_max_wind_time_list = []
EW_drct_mode_list = []
EW_speed_mean_list = []
EW_max_wind_pres_list = []

# get list of folders and sounding file paths
print 'Getting files'

for item in os.listdir(path):
    if item.endswith('.cls'):
        file_date = pd.to_datetime(item.split('Res_')[1].split('.cls')[0], format='%Y%m%d', errors='ignore')

        # read in file if within specified time
        if (file_date >= input_start_date) and (file_date <= input_end_date):

            file_path = path+'/'+item

                #print 'Processing %s'%item

            # -------- READ IN FILE ---------
            break_lines = []
            header_lines = []

            file_data = []
            format_list = []
            time_list = []
            lat = []
            lon = []
            header_found = False

            htcol_list = []
            prescol_list = []
            tempcol_list = []
            dewcol_list = []
            spdcol_list = []
            drctcol_list = []

            with open(file_path, 'r') as rawdata:

                # read in sounding data for sounding 
                for il, line in enumerate(rawdata):

                    if 'data type' in line.lower():
                        break_lines.append(il)

                    if 'site id' in line.lower():
                        site = line.split('Site ID:')[1].strip(' ').split(',')[0].strip(' ')
                        format_list.append(re.sub('[^A-Za-z0-9]+', '', site))

                    if 'release location' in line.lower():
                        lat.append(line.split('S,')[1].split(',')[1].strip(' '))
                        lon.append(line.split('S,')[1].split(',')[0].strip(' '))

                    if header_found == True:
                        header_split = line.lower().split()[:]

                        htcol_list.append(header_split.index('alt'))
                        prescol_list.append(header_split.index('press'))
                        tempcol_list.append(header_split.index('temp'))
                        dewcol_list.append(header_split.index('dewpt'))
                        spdcol_list.append(header_split.index('spd'))
                        drctcol_list.append(header_split.index('dir'))

                        header_found = False

                    if 'nominal' in line.lower():
                        header_lines.append(il+1)
                        header_found = True

    #                        full_time = line.split('(y,m,d,h,m,s):')[1].strip(' ')
    #                        year = full_time[:4]
    #                        month = full_time[6:8]
    #                        day = full_time[10:12]
    #                        time = full_time[14:16]+full_time[17:19]+full_time[20:22]
    #                        time_list.append(year+month+day+time)

                    if 'utc' in line.lower(): 
                        full_time = line.split('(y,m,d,h,m,s):')[1].strip(' ')
                        year = full_time[:4]
                        month = full_time[6:8]
                        day = full_time[10:12]
                        hour = full_time[14:16]
                        minute = full_time[17:19]
                        second = full_time[20:22]
                        time = hour+minute+second
    #                        print site
    #                        print hour
    #                        print minute
    #                        print year+month+day+time
    #                        if int(minute) < 30:
    #                            minute = '00'
    #                        if int(minute) >= 30:
    #                            if int(hour) < 10:
    #                                hour = int(hour) + 1
    #                                hour = '0' + str(hour)
    #                            else:
    #                                hour = str(int(hour) + 1)
    #                            minute = '00'
    #                        
    #                        second = '00' 
    #                        
    #                        if int(hour) == 24:
    #                            hour = '00'
    #                            if int(day) < 10:
    #                                day = int(day) + 1
    #                                day = '0' + str(day)
    #                            else:
    #                                day = str(int(day) + 1)
    #                            
    #                        print hour
    #                        print minute
    #                            
    #                        time = hour+minute+second
    #                        print year+month+day+time
                        time_list.append(year+month+day+time)

                    length=il

                break_lines.append(il)

            # go through and find soundings from correct station
                for zz in range(len(header_lines)):

                    # take only soundings from station specifed 
                    if format_list[zz] == Station: 

                        if zz == len(header_lines)-1:
                            file_data = np.genfromtxt(file_path, skip_header=header_lines[zz]+3)
                        else:
                            file_data = np.genfromtxt(file_path, skip_header=header_lines[zz]+3, skip_footer=length+1-break_lines[zz+1])

                        height = file_data[:,htcol_list[zz]]
                        pres = file_data[:,prescol_list[zz]]
                        temp = file_data[:,tempcol_list[zz]]
                        dew = file_data[:,dewcol_list[zz]]
                        wspd = file_data[:,spdcol_list[zz]]*1.94 # convert to knots
                        drct = file_data[:,drctcol_list[zz]]

                        drct[drct > 360] = 0
                        wspd[wspd > 998] = 0

                        height[height == 99999.0] = np.nan
                        temp[temp == 999.00] = np.nan
                        dew[dew == 999.00] = np.nan

                        data=dict(zip(('hght','pres','temp','dwpt','sknt','drct'),(height,pres,temp,dew,wspd,drct)))

                        # get the indices to use to get data between user defined pressures
                        ind_pres = np.nonzero(np.logical_and(data['pres']>=min_pres, data['pres']<=max_pres))

                        #get the correct data
                        pres = gaussian_filter(data['pres'][ind_pres], sigma=1.5)
                        #print 'hereeeeeeeeeeeeeeeeeeeeeeeeee'
                        #print len(pres)
                        hght = gaussian_filter(data['hght'][ind_pres], sigma=1.5)
                        
                        #print 'hght', hght
                        
                        #hght_station = hght[0]
                        #hght = hght - hght_station
                        
                        #print 'hght_station', hght_station
                        #print 'hght', hght
                        
                        temp = data['temp'][ind_pres]
                        speed = gaussian_filter(data['sknt'][ind_pres], sigma=1.5)
                        drct = data['drct'][ind_pres]
                        dwpt = data['dwpt'][ind_pres]

                        #calculate some new variables                              
                        e_sat = VaporPressure(temp)

                        e = VaporPressure(dwpt)

                        w = (e*Rs_da) / (Rs_v*(pres*100. - e))
                        q = w / (w+1) #kg/kg

                        RH = (e / e_sat)*100.

                        pres_s = data['pres'][0]

                        print time_list[zz][:12]

                        # get the current sounding file name and use name to get a time
                        file_time = datetime.datetime.strptime(time_list[zz][:12], file_in_dt_fmt)
                        print file_time    

                        # round up hour if greater than 10 mins
                        if file_time.minute/10 >= 1:
                            file_time = file_time.replace(second=0, microsecond=0, minute=0, hour=file_time.hour)+timedelta(hours=1)

                        elif file_time.minute/10 < 1:
                            file_time = file_time.replace(second=0, microsecond=0, minute=0, hour=file_time.hour)

                        else:
                            print 'DATETIME ERROR'

                        # if 10/31 23Z doesn't round up, round it up anyway because it was meant for 0Z
                        if file_time == datetime.datetime(2018, 10, 31, 23):
                            file_time = datetime.datetime(2018, 11, 1, 0)
                        print file_time

                        # convert time to datetime to use as x coordinate
                        file_time = pd.to_datetime(file_time, format='%Y%m%d%H', errors='ignore')
                        length_pres = len(pres)
                        time = [file_time]*length_pres

                        # get variable from low pressure to 3000m (used to find max wind)
                        height_cutoff_max = where(hght <= max_search_hgt)
                        #print 'hght', hght
                        #print 'max_search_hgt', max_search_hgt
                        #print 'height_cutoff_max', height_cutoff_max
                        hght_below_max = hght[height_cutoff_max]
                        speed_below_max = speed[height_cutoff_max]
                        pres_below_max = pres[height_cutoff_max]
                        drct_below_max = drct[height_cutoff_max]

                        found_index = np.nan

                        # find local max(s) and min(s)
                        max_vals_index = argrelextrema(speed_below_max,np.greater)[0]
                        min_vals_index = argrelextrema(speed_below_max,np.less)[0]

                        max_vals = speed_below_max[max_vals_index]
                        min_vals = speed_below_max[min_vals_index]
                        
                        max_vals_heights = hght_below_max[max_vals_index]
                        min_vals_heights = hght_below_max[min_vals_index]

                        max_vals_pres = pres_below_max[max_vals_index]
                        min_vals_pres = pres_below_max[min_vals_index]
                        

                        # add the top search as doesn't include end mins
                        min_vals_index = np.append(min_vals_index, [len(speed_below_max)-1])
                        min_vals = np.append(min_vals, speed_below_max[-1])
                        
                        min_vals_heights = np.append(min_vals_heights, hght_below_max[-1])
                        
                        max_vals_pres = np.append(max_vals_pres,pres_below_max[-1])

                        print 'max_vals_index', max_vals_index
                        print 'max_vals', max_vals
                        print 'min_vals_index', min_vals_index
                        print 'min_vals', min_vals
                        print '

                        for z in range(len(max_vals)):

                            if(max_vals[z] > max_wind_threshold): # make sure max wind is above threashold

                                if(max_vals[z]-min_vals[z]) > decrease_to_min:
                                    found_index = z
                                    total_count_obs = total_count_obs + 1
                                    print total_count_obs
                                    break
                        #max_wind = max_vals[z]
                        #min_wind = min_vals[z]
                        #max_wind_height = max_vals_heights[z]
                        #min_wind_height = min_vals_heights[z]
                        
print total_count_obs
                        