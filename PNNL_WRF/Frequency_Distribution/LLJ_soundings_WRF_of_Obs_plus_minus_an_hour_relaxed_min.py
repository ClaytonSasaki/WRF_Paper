'''
First finds LLJ soundings that match criteria from observations. Then finds if those are found in the PNNL WRF data.

Also, get distribution of height and speed
'''

import wrf
from netCDF4 import Dataset
from matplotlib.cm import get_cmap
import cartopy.crs as crs
import cartopy
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import math
from scipy.interpolate import interp1d
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
from matplotlib.ticker import MaxNLocator

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

outpath = '/home/disk/meso-home/crs326/Documents/Research/PNNL_WRF/Frequency_Distribution/'

Station = 'VillaDeMariaDelRioSeco'

obs_found_time_list = []

total_count_obs = 0

#criteria
# 1: original criteria from Oliviera et al 2018
crit1 = [1, 3000, 4000, 23.326, 7.775]
# 2: high - allows for jets centered higher in the atmosphere, and a bit weaker compared to mean flow
crit2 = [2 ,3600, 6000, 23.326, 7.775]

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
                        pres = data['pres'][ind_pres]
                        #print 'hereeeeeeeeeeeeeeeeeeeeeeeeee'
                        #print len(pres)
                        hght = data['hght'][ind_pres]
                        
                        #print 'hght', hght
                        
                        #hght_station = hght[0]
                        #hght = hght - hght_station
                        
                        #print 'hght_station', hght_station
                        #print 'hght', hght
                        
                        temp = data['temp'][ind_pres]
                        speed = data['sknt'][ind_pres]
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
                        #print 'speed', speed
                        #print 'speed_below_max', speed_below_max
                        # find max wind
                        max_wind = amax(speed_below_max)
                        #print 'max_wind', max_wind

                        # find index of max wind
                        max_wind_index_list = where(speed_below_max==max_wind)
                        #print 'max_wind_index_list', max_wind_index_list

                        # find max wind heights and at lowest height as defined in Oliviera et al 2018
                        max_wind_height_list = hght_below_max[max_wind_index_list]
                        max_wind_height = amin(max_wind_height_list)
                        #print 'max_wind_height_list', max_wind_height_list

                    # find lowest height (highest pressure) with maximum wind as laid out in Oliviera et al 2018
                        if max_wind > max_wind_threshold and max_wind_height < max_search_hgt:

                            # find index of min height for max wind
                            max_wind_height_index = where(hght == max_wind_height)

                            # find pres for max wind
                            max_wind_pres = pres_below_max[max_wind_height_index][0]

                            # find most common direction around max wind speed
                            # NOTE: If you make the layer too wide and peak is too close 
                            drct_near_max_wind_height = drct_below_max[int(max_wind_height_index[0])-12:int(max_wind_height_index[0])+12]
                            drct_mode, count = stats.mode(drct_near_max_wind_height)

                            # avg wind around max wind speed
                            speed_near_max_wind_height = speed_below_max[int(max_wind_height_index[0])-10:int(max_wind_height_index[0])+10]
                            speed_mean = np.nanmean(speed_near_max_wind_height)
                            #print 'hght', hght
                            #print 'min_search_hgt', min_search_hgt
                            #print 'max_wind_height', max_wind_height
                            # find variables from max wind  up to 4000m (use to find if it is a peak)

                            # function which returns True when constraints are satisfied.
                            height_cutoff_funct = lambda hght: hght <=min_search_hgt and hght >= max_wind_height 

                            # Apply constraints element-wise to the dists array.
                            height_cutoff_min = np.vectorize(height_cutoff_funct)(hght) 
                            height_cutoff_min = np.where(height_cutoff_min)[0] # Get output. 
                            #print 'max height for min wind', hght[height_cutoff_min][-1]
                            speed_min = speed[height_cutoff_min]
                            #print 'speed_min', speed_min

                            min_found = False

                            while(min_found == False):

                                try: 
                                    ind_first_min = argrelextrema(speed_min, np.less)[0][0]
                                    #print 'ind_first_min', ind_first_min
                                    current_min = speed_min[ind_first_min]
                                    #print 'current_min', current_min
                                    speed_min = speed_min[ind_first_min:] #get rid of everything under current min
                                    #print 'speed_min', speed_min
                                    ind_next_max = np.argmax(speed_min_)

                                    if(speed_min[ind_next_max] - current_min) > (max_wind - current_min) or len(speed_min) == 1:

                                        min_found = True
                                        #print 'min_found' 
                                except: 
                                    ind_first_min = np.argmin(speed_min) # used if last value is min so there is no minima
                                    current_min = speed_min[ind_first_min]
                                    #print 'current_min', current_min
                                    break


                            # find index for min wind
                            min_wind_index = where(speed == current_min)

                             # find pres for min wind
                            min_wind_pres = pres[min_wind_index][0]
                            print 'min_wind_pres', min_wind_pres

                            #print 'decrease_to_min', max_wind - current_min

                            # see if difference between max wind and min wind is large enough
                            if (max_wind - current_min) > decrease_to_min:

                                # If yes a LLJ is considered to be found
                                #print 'LLJ'

                                if len(drct_mode) != 0:
                                    total_count_obs = total_count_obs + 1

                                    obs_found_time_list.append(file_time)

        
        
print 'total_count_obs', total_count_obs
print 'obs_found_time_list', obs_found_time_list
        
#added_times = ['201810311106', '201811092322', '201811171130', '201812130231', '201812130830']
#rounded
#added_times = ['201810311100', '201811100000', '201811171200', '201812130300', '201812130900']

#print type(obs_found_time_list[0])

#for a in range(len(added_times)):
#    added_time_dt = datetime.datetime.strptime(added_times[a][:10], '%Y%m%d%H')
#    obs_found_time_list.append(added_time_dt)

#print 'obs_found_time_list', obs_found_time_list
        
#print type(obs_found_time_list[-1])    
        
        
        
        
        
        
        

        
        
        
        
        
path = '/home/disk/monsoon/relampago/raw/wrf/'  

#pressure to plot
pressure_top = 450
pressure_bottom = 967

#change the start and end days to be plotted
start_date = 2018110100
end_date = 2018121821

#################################################################

# get dates from times to only run for correct day folders
input_start_day = pd.to_datetime(str(start_date)[:8], format='%Y%m%d', errors='ignore')
input_end_day = pd.to_datetime(str(end_date)[:8], format='%Y%m%d', errors='ignore')

time_plot = []

total_count_WRF = 0

WRF_found_time_list = []
WRF_speed_mean_list = []
WRF_max_wind_pres_list = []

full_mb_p = np.arange(pressure_bottom,pressure_top-1,-1)
#print 'full_mb_p', full_mb_p
crit = [2 ,3600, 6000, 23.326, 11.663]

#crit reduced: 5 knots max wind, 3 knots speed shear
#crit = [2 ,3500, 6000, 18.326, 5.831]

crit_num = crit[0]
max_search_hgt = crit[1]
min_search_hgt = crit[2]
max_wind_threshold = crit[3]
decrease_to_min = crit[4]

time_list = []

# read in LLJ list file
# For COR
if 'CordobaAero' in Station:
    lat_s = -31.298
    lon_s = -64.212
    area = 'COR' # used for file name and text
    Station_short = 'COR'
# For VMRS
if 'VillaDeMariaDelRioSeco' in Station:
    lat_s = -29.906
    lon_s = -63.726
    area = 'VMRS' # used for file name and text
    Station_short = 'VMRS'

    WRF_time_list = []

for i, LLJ_list in enumerate([obs_found_time_list]):

    total_count = 0
    total_ERA_LLJ_count = 0
    
    ERA_max_wind_time_list = []
    ERA_drct_mode_list = []
    ERA_speed_mean_list = []
    ERA_max_wind_pres_list = []
    
    ERA_LLJ_max_wind_time_list = []
    ERA_LLJ_drct_mode_list = []
    ERA_LLJ_speed_mean_list = []
    ERA_LLJ_max_wind_pres_list = []
    
    OBS_ERA_LLJ_speed = [] #list of mean speeds from obs for times where LLJ found in obs AND in ERA5
    
    for zz, LLJ_date_dt in enumerate(LLJ_list):

        #print 'LLJ_date_dt', LLJ_date_dt
        #LLJ_date_dt = LLJ_date_dt+timedelta(hours=1)
        hour = int(LLJ_date_dt.hour) # get hour from LLJ sounding

        # add all files to a list
        for folder_date in os.listdir(path):
            
            folder_path = os.path.join(path, folder_date)
            
            folder_date_dt = pd.to_datetime(folder_date, format='%Y%m%d', errors='ignore')
            
            LLJ_date_dt = LLJ_date_dt.replace(hour=0) # replace hour to 0 so that date can be checked without consideration of hour
            
            if folder_date_dt == LLJ_date_dt:
                
                print 'found day'
                
                LLJfound = False
                
                # start with search_hour equal to hour from observations, if a LLJ is found the stop, if not try the hour before and then hour after
                for num in range(3):
                    
                    print 'LLJ_date_dt', LLJ_date_dt
                    
                    print 'num', num
                    
                    if(LLJfound == False):

                        print 'num use', num

                        if num == 0:
                            search_hour = hour
                        elif num == 1:
                            search_hour = hour-1
                        elif num == 2:
                            search_hour = hour+1

                        print 'search_hour', search_hour

                        for hourly_file in os.listdir(folder_path):

                            file_path = os.path.join(folder_path, hourly_file)

                            #print 'hourly_file', hourly_file

                            file_hour = int(hourly_file[22:24])

                            #print 'hour', hour
                            #print type(hour)
                            #print 'file_hour', file_hour
                            #print type(file_hour)
                            #print ' '


                            if(file_hour == search_hour):

                                print 'found hour'

                                nc = Dataset(file_path,'r')

                                u = wrf.getvar(nc, 'ua', units='kt')[:,:,:-1] # in kts
                                v = wrf.getvar(nc, 'va', units='kt')[:,:-1,:] # in kts
                                hght = wrf.getvar(nc, 'z')

                                #speed = np.sqrt(np.power(u,2) + np.power(v, 2))
                                #wind_dir_trig_to = np.arctan2(u/speed, v/speed) 
                                #wind_dir_trig_to_degrees = wind_dir_trig_to * 180/math.pi
                                #wind_dir_trig_from_degrees = wind_dir_trig_to_degrees + 180

                                speed, drct = wrf.getvar(nc, 'wspd_wdir', units='kt')

                                lats = wrf.getvar(nc, 'XLAT', meta=False)
                                lons = wrf.getvar(nc, 'XLONG', meta=False)

                                #print 'lats', lats.shape
                                #print 'lats', lats.shape
                                #print 'pres', pres.shape

                                #PB = wrf.getvar(nc, 'PB')
                                #P = wrf.getvar(nc, 'P')
                                #pres = (PB + P)/100 # convect to hPa

                                pres = wrf.getvar(nc, 'pressure')

                                #print pres

                                Station_xy = wrf.ll_to_xy(nc, lat_s, lon_s)

                                # get the indices to use to get data between user defined levels
                                pres_top_WRF_i = (np.abs(pres[:,Station_xy[0],Station_xy[1]] - pressure_top)).argmin()
                                pres_bottom_WRF_i = (np.abs(pres[:,Station_xy[0],Station_xy[1]] - pressure_bottom)).argmin()

                                #print 'pres_top_WRF_i', pres_top_WRF_i
                                #print 'pres_bottom_WRF_i', pres_bottom_WRF_i

        #                        interp_hght = wrf.vinterp(nc, field=hght, vert_coord='pres', interp_levels=full_mb_p)
        #                        interp_speed = wrf.vinterp(nc, field=speed, vert_coord='pres', interp_levels=full_mb_p)
        #                        interp_pres = wrf.vinterp(nc, field=pres, vert_coord='pres', interp_levels=full_mb_p)
        #                        interp_drct = wrf.vinterp(nc, field=drct, vert_coord='pres', interp_levels=full_mb_p)
        #
        #                        #print 'interp_drct', interp_drct.shape
        #                        
        #                        # narrow down speed, pres, and drct to one location and set pres boundaries
        #                        hght = interp_hght[:,Station_xy[0],Station_xy[1]]
        #                        speed = interp_speed[:,Station_xy[0],Station_xy[1]]
        #                        pres = interp_pres[:,Station_xy[0],Station_xy[1]]
        #                        drct = interp_drct[:,Station_xy[0],Station_xy[1]]


        #                        # narrow down speed, pres, and drct to one location and set pres boundaries
        #                        hght = hght[pres_bottom_WRF_i:pres_top_WRF_i,Station_xy[0],Station_xy[1]]
        #                        speed = speed[pres_bottom_WRF_i:pres_top_WRF_i,Station_xy[0],Station_xy[1]]
        #                        pres = pres[pres_bottom_WRF_i:pres_top_WRF_i,Station_xy[0],Station_xy[1]]
        #                        print 'pres', pres
        #                        drct = drct[pres_bottom_WRF_i:pres_top_WRF_i,Station_xy[0],Station_xy[1]]
        #                    

                                # narrow down speed, pres, and drct to one location and set pres boundaries
                                hght = hght[:,Station_xy[0],Station_xy[1]]
                                speed = speed[:,Station_xy[0],Station_xy[1]]
                                pres = pres[:,Station_xy[0],Station_xy[1]]
                                drct = drct[:,Station_xy[0],Station_xy[1]]

                                #print 'pres', pres

                                #print 'full_mb_p', full_mb_p

                                f_hght = interp1d(pres,hght) 
                                f_speed = interp1d(pres,speed)
                                f_pres = interp1d(pres,pres)
                                f_drct = interp1d(pres,drct)

                                hght = f_hght(full_mb_p)
                                speed = f_speed(full_mb_p)
                                pres = f_pres(full_mb_p)
                                drct = f_drct(full_mb_p)

                                WRF_file_time = folder_date_dt.replace(hour=hour)

                                # get variable from low pressure to 3000m (used to find max wind)
                                height_cutoff_max = where(hght <= max_search_hgt)
                                #print 'hght', hght
                                #print 'max_search_hgt', max_search_hgt
                                #print 'height_cutoff_max', height_cutoff_max
                                hght_below_max = hght[height_cutoff_max]
                                speed_below_max = speed[height_cutoff_max]
                                pres_below_max = pres[height_cutoff_max]
                                drct_below_max = drct[height_cutoff_max]
                                #print 'speed', speed
                                #print 'speed_below_max', speed_below_max
                                # find max wind
                                max_wind = amax(speed_below_max)
                                #print 'max_wind', max_wind

                                # find index of max wind
                                max_wind_index_list = where(speed_below_max==max_wind)
                                #print 'max_wind_index_list', max_wind_index_list

                                # find max wind heights and at lowest height as defined in Oliviera et al 2018
                                max_wind_height_list = hght_below_max[max_wind_index_list]
                                max_wind_height = amin(max_wind_height_list)
                                #print 'max_wind_height_list', max_wind_height_list

                            # find lowest height (highest pressure) with maximum wind as laid out in Oliviera et al 2018
                                if max_wind > max_wind_threshold and max_wind_height < max_search_hgt:

                                    # find index of min height for max wind
                                    max_wind_height_index = where(hght == max_wind_height)

                                    # find pres for max wind
                                    max_wind_pres = pres_below_max[max_wind_height_index][0]

                                    # find most common direction around max wind speed
                                    # NOTE: If you make the layer too wide and peak is too close 
                                    drct_near_max_wind_height = drct_below_max[int(max_wind_height_index[0])-12:int(max_wind_height_index[0])+12]
                                    drct_mode, count = stats.mode(drct_near_max_wind_height)

                                    # avg wind around max wind speed
                                    speed_near_max_wind_height = speed_below_max[int(max_wind_height_index[0])-10:int(max_wind_height_index[0])+10]
                                    speed_mean = np.nanmean(speed_near_max_wind_height)
                                    #print 'hght', hght
                                    #print 'min_search_hgt', min_search_hgt
                                    #print 'max_wind_height', max_wind_height
                                    # find variables from max wind  up to 4000m (use to find if it is a peak)

                                    # function which returns True when constraints are satisfied.
                                    height_cutoff_funct = lambda hght: hght <=min_search_hgt and hght >= max_wind_height 

                                    # Apply constraints element-wise to the dists array.
                                    height_cutoff_min = np.vectorize(height_cutoff_funct)(hght) 
                                    height_cutoff_min = np.where(height_cutoff_min)[0] # Get output. 
                                    #print 'max height for min wind', hght[height_cutoff_min][-1]
                                    speed_min = speed[height_cutoff_min]
                                    #print 'speed_min', speed_min

                                    min_found = False

                                    while(min_found == False):

                                        try: 
                                            ind_first_min = argrelextrema(speed_min, np.less)[0][0]
                                            #print 'ind_first_min', ind_first_min
                                            current_min = speed_min[ind_first_min]
                                            #print 'current_min', current_min
                                            speed_min = speed_min[ind_first_min:] #get rid of everything under current min
                                            #print 'speed_min', speed_min
                                            ind_next_max = np.argmax(speed_min_)

                                            if(speed_min[ind_next_max] - current_min) > (max_wind - current_min) or len(speed_min) == 1:

                                                min_found = True
                                                #print 'min_found' 
                                        except: 
                                            ind_first_min = np.argmin(speed_min) # used if last value is min so there is no minima
                                            current_min = speed_min[ind_first_min]
                                            #print 'current_min', current_min
                                            break


                                    # find index for min wind
                                    min_wind_index = where(speed == current_min)

                                     # find pres for min wind
                                    min_wind_pres = pres[min_wind_index][0]
                                    #print 'min_wind_pres', min_wind_pres

                                    #print 'decrease_to_min', max_wind - current_min

                                    #print 'max - min', max_wind - current_min

                                    # see if difference between max wind and min wind is large enough
                                    if (max_wind - current_min) > decrease_to_min:

                                        # If yes a LLJ is considered to be found
                                        #print 'LLJ'

                                        if len(drct_mode) != 0:

                                            LLJfound = True

                                            total_count_WRF = total_count_WRF + 1

                                            print 'WRF LLJ count', total_count_WRF, 'num', num

                                            WRF_found_time_list.append(file_time)

                                            WRF_speed_mean_list.append(speed_mean)

                                            WRF_max_wind_pres_list.append(max_wind_pres)


fig2, axs2 = plt.subplots(1, 2, sharey=True, tight_layout=True, figsize=(20,20))
        
# We can set the bins with the `bins` kwarg
counts_pres, bins_pres, patches_pres = axs2[0].hist(WRF_max_wind_pres_list, bins=np.arange(650,975,25), weights=np.ones(len(WRF_max_wind_pres_list)) / len(WRF_max_wind_pres_list) * 100)
counts_speed, bins_speed, patches_speed = axs2[1].hist(WRF_speed_mean_list, bins=np.arange(20,48,2), weights=np.ones(len(WRF_max_wind_pres_list)) / len(WRF_max_wind_pres_list) * 100)

speed_mean_median = np.median(WRF_speed_mean_list)
print 'speed_mean_median', speed_mean_median

axs2[1].axvline(x=speed_mean_median, color='black', linewidth=3)

fig2.text(0.001,0.98, 'Station: %s'%(Station_short))
fig2.text(0.75, 0.98, 'Number of Soundings: %s'%(total_count_WRF))

axs2[0].set_xlabel('Pressure Level of Max Wind (hPa)')
#axs2[0].xaxis.set_ticks([700,750,800,850,900])
for tick in axs2[0].get_xticklabels():
    tick.set_rotation(45)
axs2[1].set_xlabel('Max Wind Speeds (kts)')
axs2[0].set_ylabel('Percentage of Soundings (%)')

axs2[0].set_xticks(bins_pres[::2])
axs2[1].set_xticks(bins_speed[::2])

axs2[0].set_ylim(0,20)

axs2[0].yaxis.set_major_locator(MaxNLocator(integer=True))
axs2[1].yaxis.set_major_locator(MaxNLocator(integer=True))

plt.savefig(outpath + 'hour_leeway_relaxed_min_f000000_%s_max_wind_pres_pdf_%d.png'%(Station, crit_num))

plt.close()



