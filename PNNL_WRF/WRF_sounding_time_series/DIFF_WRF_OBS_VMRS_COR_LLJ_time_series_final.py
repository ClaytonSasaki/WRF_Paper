#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 16:50:00 2021

@author: crs326
"""

import matplotlib
matplotlib.use('Agg')

import wrf
from netCDF4 import Dataset
import os
import pandas as pd
import numpy as np
from numpy import exp,where,ma,cos,sin,pi,amax,amin
import matplotlib.dates as dates
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.interpolate import interp1d
from scipy import stats
import re
import math
import datetime
from datetime import timedelta, date
import matplotlib.dates as mdates
from scipy.signal import argrelextrema
plt.rcParams.update({'font.size': 28})

file_in_dt_fmt_obs = '%Y%m%d%H%M'
file_in_dt_fmt_wrf = '%Y-%m-%d_%H:%M:%S'

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

# function to truncate color map to better fit data, used right now with jet and for winds
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

cmap = plt.get_cmap('jet')
new_cmap = truncate_colormap(cmap, 0.2, 1)

# --------------------- START OF MAIN BODY OF FILE ------------------------

##########CHANGE#################################################

# stations to plot time series
station_list = ['CordobaAero', 'VillaDeMariaDelRioSeco']

#criteria

# 2: high - allows for jets centered higher in the atmosphere, and a bit weaker compared to mean flow
crit = [2 ,3600, 6000, 23.326, 11.663]

colorN = ['black', 'gray']
colorS = ['black', 'gray']

markerN = ['s', 'X']
markerS = ['s', 'X']

sizeN = [200, 500]
sizeS = [200, 500]

# path where files are saved
outpath = '/home/disk/meso-home/crs326/Documents/Research/PNNL_WRF/'

# min and max pressure to plot
min_pres = 450
max_pres = 975

#change the start and end times to be plotted
x_start_day = 20181203
x_end_day = 20181218

#################################################################

# path to get data from
path_obs = '/home/disk/monsoon/relampago/qc/sounding/eol_composite/High_Res_v1.1'
path_wrf = '/home/disk/monsoon/relampago/raw/wrf/'

# get dates from times to only run for correct day folders
input_start_day = pd.to_datetime(str(x_start_day), format='%Y%m%d', errors='ignore')
input_end_day = pd.to_datetime(str(x_end_day), format='%Y%m%d', errors='ignore')

#num_days = int(str(x_end_time)[6:8])-int(str(x_start_time)[6:8])
#print num_days

#set for now because of issues with cross months
num_days = 16

# go through code for each Station in station_list
for Station in station_list:
    
    TimeList1 = []
    TimeList2 = []
    
    if Station == 'CordobaAero':
        
        full_station = 'Cordoba'
        short_station = 'COR'
        lat_s = -31.298
        lon_s = -64.212
        
    if Station == 'VillaDeMariaDelRioSeco':
    
        full_station = 'VMRS'
        short_station = 'VMRS'
        lat_s = -29.906
        lon_s = -63.726
        
    else:
        
        full_station = Station
        short_station = Station
    
    print full_station

# FOR 3 SEPERATE PLOTS
    # initalize figures for plots and change size of figure depending on number of days plotting
    fig, ax = plt.subplots(figsize=(2*num_days+13, 16))
    
    fig2, ax2 = plt.subplots(figsize=(2*num_days+13, 16))
    
    fig3, ax3 = plt.subplots(figsize=(2*num_days+13, 16))
    
    ##for GOES full time series use this##

# FOR 1 LONG PLOT
#    # initalize figures for plots and change size of figure depending on number of days plotting
#    fig, ax = plt.subplots(figsize=(2*48+4, 16))
#    
#    fig2, ax2 = plt.subplots(figsize=(2*48+4, 16))
#    
#    fig3, ax3 = plt.subplots(figsize=(2*48+4, 16))
    
    #############

    crit_num = crit[0]
    max_search_hgt = crit[1]
    min_search_hgt = crit[2]
    max_wind_threshold = crit[3]
    decrease_to_min = crit[4]
    
    full_mb_p = np.arange(max_pres,min_pres-1,-1)
    
    # ------------------------ OBS -----------------------------
    
    time_list = []
    sounding_filepaths = []
    list_folder_paths = [] # list that will contain folders (path+name)
    total_count = 0
    calm_count = 0
    max_wind_time_list = []
    drct_mode_list = []
    speed_mean_list = []
    max_wind_pres_list = []
    
    all_header_lines = []
   
    # get list of folders and sounding file paths
    print 'Getting files'
    
    sorted_files = sorted(os.listdir(path_obs))
    
    for item in sorted_files:
        if item.endswith('.cls'):
            file_date = pd.to_datetime(item.split('Res_')[1].split('.cls')[0], format='%Y%m%d', errors='ignore')
              
            # read in file if within specified time
            if (file_date >= input_start_day) and (file_date <= input_end_day):
    
                file_path = path_obs+'/'+item
            
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
                            
                            print 'file_path', file_path
                        
                            if zz == len(header_lines)-1:
                                file_data = np.genfromtxt(file_path, skip_header=header_lines[zz]+3)
                            else:
                                file_data = np.genfromtxt(file_path, skip_header=header_lines[zz]+3, skip_footer=length+1-break_lines[zz+1])
                
                            height_obs = file_data[:,htcol_list[zz]]
                            pres_obs = file_data[:,prescol_list[zz]]
                            temp_obs = file_data[:,tempcol_list[zz]]
                            dew_obs = file_data[:,dewcol_list[zz]]
                            wspd_obs = file_data[:,spdcol_list[zz]]*1.94 # convert to knots
                            drct_obs = file_data[:,drctcol_list[zz]]
                               
                            drct_obs[drct_obs > 360] = 0
                            wspd_obs[wspd_obs > 998] = 0
                            
                            height_obs[height_obs == 99999.0] = np.nan
                            temp_obs[temp_obs == 999.00] = np.nan
                            dew_obs[dew_obs == 999.00] = np.nan
                    
                            data=dict(zip(('hght','pres','temp','dwpt','sknt','drct'),(height_obs,pres_obs,temp_obs,dew_obs,wspd_obs,drct_obs)))
                                            
                            # get the indices to use to get data between user defined pressures
                            #ind_pres = np.nonzero(np.logical_and(data['pres']>=min_pres, data['pres']<=max_pres))
                        
                            #get the correct data
                            pres_obs = data['pres']
                            #print 'pres_obs', pres_obs
                            #print 'hereeeeeeeeeeeeeeeeeeeeeeeeee'
                            #print len(pres)
                            
                            pres_s_obs = pres_obs[0]
                            
                            # height MSL (Mean sea level)
                            hght_obs_MSL = data['hght']

                            #print 'hght_MSL', hght_MSL

                            # elevation
                            elev = hght_obs_MSL[0]

                            #print 'elev', elev

                            # height AGL (Above ground level)
                            #hght = np.subtract(data['hght'],elev)
                            hght_obs = data['hght']

                            #print 'hght', hght
                            
                            temp_obs = data['temp']
                            speed_obs = data['sknt']
                            drct_obs = data['drct']
                            dwpt_obs = data['dwpt']
                            
                            f_hght_obs = interp1d(pres_obs,hght_obs)
                            f_temp_obs = interp1d(pres_obs,temp_obs)
                            f_speed_obs = interp1d(pres_obs,speed_obs)
                            f_pres_obs = interp1d(pres_obs,pres_obs)
                            f_drct_obs = interp1d(pres_obs,drct_obs)
                            f_dwpt_obs = interp1d(pres_obs,dwpt_obs)
                            
                            if(pres_s_obs<=max_pres):
                                
                                max_plot_pres = math.floor(pres_s_obs)
                                
                                if(pres_obs[-1]>=min_pres):
                                    min_plot_pres = math.ceil(pres_obs[-1])
                                    part_mb_p_obs = np.arange(max_plot_pres,min_plot_pres-1,-1)
                                else:   
                                    part_mb_p_obs = np.arange(max_plot_pres,min_pres-1,-1)
                                    
                            elif(pres_obs[-1]>=min_pres):
                                min_plot_pres = math.ceil(pres_obs[-1])
                                part_mb_p_obs = np.arange(max_pres,min_plot_pres-1,-1)
                                
                            else:
                                part_mb_p_obs = full_mb_p
                                max_plot_pres = max_pres
                            
                            #print 'pres_obs', pres_obs
                            #print 'len(pres_obs)', len(pres_obs)
                            
                            #print 'part_mb_p_obs', part_mb_p_obs
                            #print 'len(part_mb_p_obs)', len(part_mb_p_obs)
                            
                            hght_obs = f_hght_obs(part_mb_p_obs)
                            temp_obs = f_temp_obs(part_mb_p_obs)
                            speed_obs = f_speed_obs(part_mb_p_obs)
                            pres_obs = f_pres_obs(part_mb_p_obs)
                            drct_obs = f_drct_obs(part_mb_p_obs)
                            dwpt_obs = f_dwpt_obs(part_mb_p_obs)
                            
                            #calculate some new variables                              
                            e_sat_obs = VaporPressure(temp_obs)
                                
                            e_obs = VaporPressure(dwpt_obs)
                            
                            w_obs = (e_obs*Rs_da) / (Rs_v*(pres_obs*100. - e_obs))
                            q_obs = w_obs / (w_obs+1) #kg/kg
                                    
                            RH_obs = (e_obs / e_sat_obs)*100.
                            
                            # do calcuations for inputs for wind barbs
                            rdir_obs = (270.-drct_obs)*(pi/180.)
                            uu_obs = ma.masked_invalid(speed_obs*cos(rdir_obs))
                            vv_obs = ma.masked_invalid(speed_obs*sin(rdir_obs))
                            
                            vp_obs = pres_obs <= pres_s_obs # valid pressures to plot
            
                            nbarbs = (~uu_obs.mask).sum()
                        
                            barb_levels = np.logspace(np.log10(min_pres), np.log10(max_pres), 30)
                        
                            pres_barbs = np.unique(np.array( [np.argmin(np.abs(_ - pres_obs)) for _ in barb_levels] ))
                            pres_barbs = pres_barbs[pres_barbs<pres_obs[vp_obs].shape[0]]  
                        
                            # get the current sounding file name and use name to get a time
                            #print 'str format', time_list[zz]
                            file_time_obs = datetime.datetime.strptime(time_list[zz][:12], file_in_dt_fmt_obs)
                            #print 'dt format, non-rounded', file_time_obs
                            
                            # round up hour if greater than 30 mins
                            file_time_obs = file_time_obs.replace(second=0, microsecond=0, minute=0, hour=file_time_obs.hour)+timedelta(hours=file_time_obs.minute//30)
                            #print 'dt format, rounded', file_time_obs
             
                            # if 10/31 23Z doesn't round up, round it up anyway because it was meant for 0Z
                            if file_time_obs == datetime.datetime(2018, 10, 31, 23):
                                file_time_obs = datetime.datetime(2018, 11, 1, 0)
                            #print 'dt format, rounded', file_time_obs
                            
                            # convert time to datetime to use as x coordinate
                            file_time_obs = pd.to_datetime(file_time_obs, format='%Y%m%d%H', errors='ignore')
                            length_pres = len(pres_obs)
                            time = [file_time_obs]*length_pres
                            
                            #print 'dt format, rounded', file_time_obs


                            ######## find jet core from observations ###########


                            # get variable from low pressure to 3000m (used to find max wind)
                            height_cutoff_max = where(hght_obs <= max_search_hgt)
                            #print 'hght', hght
                            #print 'max_search_hgt', max_search_hgt
                            #print 'height_cutoff_max', height_cutoff_max
                            hght_below_max = hght_obs[height_cutoff_max]
                            speed_below_max = speed_obs[height_cutoff_max]
                            pres_below_max = pres_obs[height_cutoff_max]
                            drct_below_max = drct_obs[height_cutoff_max]
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
                                max_wind_height_index = where(hght_obs == max_wind_height)

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
                                height_cutoff_funct = lambda hght_obs: hght_obs <=min_search_hgt and hght_obs >= max_wind_height 

                                # Apply constraints element-wise to the dists array.
                                height_cutoff_min = np.vectorize(height_cutoff_funct)(hght_obs) 
                                height_cutoff_min = np.where(height_cutoff_min)[0] # Get output. 
                                #print 'max height for min wind', hght[height_cutoff_min][-1]
                                speed_min = speed_obs[height_cutoff_min]
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
                                min_wind_index = where(speed_obs == current_min)

                                 # find pres for min wind
                                min_wind_pres = pres_obs[min_wind_index][0]
                                #print 'min_wind_pres', min_wind_pres

                                #print 'decrease_to_min', max_wind - current_min

                                # see if difference between max wind and min wind is large enough
                                if (max_wind - current_min) > decrease_to_min:

                                    # If yes a LLJ is considered to be found
                                    #print 'LLJ'

                                    if len(drct_mode) != 0:
                                        cc='gray'
                                        ms='X'
                                        ss=500
                                        ax.scatter(dates.date2num(file_time_obs), max_wind_pres, color=cc, marker=ms, s=ss, zorder=3)



                            ######## get corresponding WRF file ################
                            
                            #get obs hour
                            obs_hour = file_time_obs.hour
                            
                            #Find corresponding WRF file and read it in
                            sorted_folders = sorted(os.listdir(path_wrf))

                            for folder_date in sorted_folders:
                                
                                folder_date_dt = pd.to_datetime(folder_date, format='%Y%m%d', errors='ignore')
                                
                                #print 'folder_date_dt', folder_date_dt
                                #print 'file_time_obs', file_time_obs
                                
                                #obs date without hour, min, sec
                                file_time_obs = file_time_obs.replace(hour=0, minute=0, second=0)
                                

                                if(folder_date_dt == file_time_obs):
                                    
                                    #print 'found day'

                                    folder_path = os.path.join(path_wrf, folder_date)

                                    sorted_hour_files = sorted(os.listdir(folder_path))

                                    for hourly_file in sorted_hour_files:
                                        
                                        file_hour = int(hourly_file[22:24])
                                        
                                        #print 'file_hour', file_hour
                                        #print 'obs_hour', obs_hour
                                        
                                        if(file_hour == obs_hour):
                                            #print 'found hourly file'
                                            
                                            print 'hourly_file', hourly_file
                                            
                                            file_path_wrf = os.path.join(folder_path, hourly_file)

                                            # get the current sounding file name and use name to get a time
                                            file_time_wrf = datetime.datetime.strptime(hourly_file[11:], file_in_dt_fmt_wrf)
                                            print file_time_wrf 

                                            nc = Dataset(file_path_wrf,'r')

                                            u_wrf = wrf.getvar(nc, 'ua', units='kt')[:,:,:-1] # in kts
                                            v_wrf = wrf.getvar(nc, 'va', units='kt')[:,:-1,:] # in kts
                                            hght_wrf = wrf.getvar(nc, 'z')

                                            #speed = np.sqrt(np.power(u,2) + np.power(v, 2))
                                            #wind_dir_trig_to = np.arctan2(u/speed, v/speed) 
                                            #wind_dir_trig_to_degrees = wind_dir_trig_to * 180/math.pi
                                            #wind_dir_trig_from_degrees = wind_dir_trig_to_degrees + 180

                                            speed_wrf, drct_wrf = wrf.getvar(nc, 'wspd_wdir', units='kt')

                                            lats_wrf = wrf.getvar(nc, 'XLAT', meta=False)
                                            lons_wrf = wrf.getvar(nc, 'XLONG', meta=False)

                                            #print 'lats', lats.shape
                                            #print 'lats', lats.shape
                                            #print 'pres', pres.shape

                                            #PB = wrf.getvar(nc, 'PB')
                                            #P = wrf.getvar(nc, 'P')
                                            #pres = (PB + P)/100 # convect to hPa

                                            pres_wrf = wrf.getvar(nc, 'pressure')

                                            #print pres

                                            Station_xy = wrf.ll_to_xy(nc, lat_s, lon_s)

                                            # get the indices to use to get data between user defined levels
                                            #pres_top_WRF_i = (np.abs(pres[:,Station_xy[0],Station_xy[1]] - pressure_top)).argmin()
                                            #pres_bottom_WRF_i = (np.abs(pres[:,Station_xy[0],Station_xy[1]] - pressure_bottom)).argmin()

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
                                            hght_wrf = hght_wrf[:,Station_xy[0],Station_xy[1]]
                                            speed_wrf = speed_wrf[:,Station_xy[0],Station_xy[1]]
                                            v_wrf = v_wrf[:,Station_xy[0],Station_xy[1]]
                                            pres_wrf = pres_wrf[:,Station_xy[0],Station_xy[1]]
                                            drct_wrf = drct_wrf[:,Station_xy[0],Station_xy[1]]

                                            #print 'pres', pres_wrf[0]

                                            #print 'pres', pres

                                            #print 'full_mb_p', full_mb_p

                                            f_hght_wrf = interp1d(pres_wrf,hght_wrf) 
                                            f_speed_wrf = interp1d(pres_wrf,speed_wrf)
                                            f_v_wrf = interp1d(pres_wrf,v_wrf)
                                            f_pres_wrf = interp1d(pres_wrf,pres_wrf)
                                            f_drct_wrf = interp1d(pres_wrf,drct_wrf)
                                            
                                            

                                            part_mb_p_wrf = part_mb_p_obs
                                            
                                            
                                            

                                            hght_wrf = f_hght_wrf(part_mb_p_wrf)
                                            speed_wrf = f_speed_wrf(part_mb_p_wrf)
                                            v_wrf = f_v_wrf(part_mb_p_wrf)
                                            pres_wrf = f_pres_wrf(part_mb_p_wrf)
                                            drct_wrf = f_drct_wrf(part_mb_p_wrf)

#                                            # convert time to datetime to use as x coordinate
#                                            length_pres = len(pres_wrf)
#                                            time = [file_time_wrf]*length_pres

                                            #if lowest pressure at wrf is greater than obs than cut off wrf to highest obs pressure
                                            if pres_obs[0] < pres_wrf[0]:
                                                obs_start_i = np.where(pres_wrf == pres_obs[0])[0][0]
                
                                                #print 'obs_start_i', obs_start_i
                                                hght_wrf = hght_wrf[obs_start_i:]
                                                speed_wrf = speed_wrf[obs_start_i:]
                                                v_wrf = v_wrf[obs_start_i:]
                                                pres_wrf = pres_wrf[obs_start_i:]
                                                drct_wrf = drct_wrf[obs_start_i:]
    
    
                                            diff_obs_wrf = vv_obs - v_wrf
        
                                            #print vv_obs
                                            #print v_wrf
                                            #print diff_obs_wrf

                                            # plot meridional wind speeds on ax where postive is from north
                                            diff_v_wind = ax.scatter(time, pres_obs, c=diff_obs_wrf, vmin=-40, vmax=40, cmap='seismic_r', zorder=1)



    
    
    
#
#    time_list = []
#
#    #GO THROUGH ALL FOLDERS
#
#    total_LLJ_count_WRF = 0
#    total_count_WRF = 0
#
#    WRF_file_time_list = []
#
#    WRF_found_time_list = []
#    WRF_speed_mean_list = []
#    WRF_max_wind_pres_list = []
#    v_plot = np.zeros((len(full_mb_p),num_days*24))
#
#
#    OBS_ERA_LLJ_speed = [] #list of mean speeds from obs for times where LLJ found in obs AND in ERA5
#
#    # add all files to a list
#    #print '1'
#
#    sorted_folders = sorted(os.listdir(path))
#
#    for folder_date in sorted_folders:
#
#        if folder_date.startswith('20'):
#
#            #print '2'
#
#            folder_date_dt = pd.to_datetime(folder_date, format='%Y%m%d', errors='ignore')
#
#            #print 'folder_date_dt', folder_date_dt
#            #print 'input_start_day', input_start_day
#            #print 'input_end_day', input_end_day
#
#            if (folder_date_dt >= input_start_day) and (folder_date_dt <= input_end_day):
#
#                #print '3'
#
#                folder_path = os.path.join(path, folder_date)
#
#                sorted_hour_files = sorted(os.listdir(folder_path))
#
#                for hourly_file in sorted_hour_files:
#
#                    #print '4'
#
#                    file_path = os.path.join(folder_path, hourly_file)
#
#                    print 'hourly_file', hourly_file
#
#                    file_hour = int(hourly_file[22:24])
#
#                    print 'file_hour', file_hour
#                    #print type(file_hour)
#                    #print ' '
#
#                    # get the current sounding file name and use name to get a time
#                    file_time = datetime.datetime.strptime(hourly_file[11:], file_in_dt_fmt)
#                    print file_time 
#
#                    nc = Dataset(file_path,'r')
#
#                    u = wrf.getvar(nc, 'ua', units='kt')[:,:,:-1] # in kts
#                    v = wrf.getvar(nc, 'va', units='kt')[:,:-1,:] # in kts
#                    hght = wrf.getvar(nc, 'z')
#
#                    #speed = np.sqrt(np.power(u,2) + np.power(v, 2))
#                    #wind_dir_trig_to = np.arctan2(u/speed, v/speed) 
#                    #wind_dir_trig_to_degrees = wind_dir_trig_to * 180/math.pi
#                    #wind_dir_trig_from_degrees = wind_dir_trig_to_degrees + 180
#
#                    speed, drct = wrf.getvar(nc, 'wspd_wdir', units='kt')
#
#                    lats = wrf.getvar(nc, 'XLAT', meta=False)
#                    lons = wrf.getvar(nc, 'XLONG', meta=False)
#
#                    #print 'lats', lats.shape
#                    #print 'lats', lats.shape
#                    #print 'pres', pres.shape
#
#                    #PB = wrf.getvar(nc, 'PB')
#                    #P = wrf.getvar(nc, 'P')
#                    #pres = (PB + P)/100 # convect to hPa
#
#                    pres = wrf.getvar(nc, 'pressure')
#
#                    #print pres
#
#                    Station_xy = wrf.ll_to_xy(nc, lat_s, lon_s)
#
#                    # get the indices to use to get data between user defined levels
#                    #pres_top_WRF_i = (np.abs(pres[:,Station_xy[0],Station_xy[1]] - pressure_top)).argmin()
#                    #pres_bottom_WRF_i = (np.abs(pres[:,Station_xy[0],Station_xy[1]] - pressure_bottom)).argmin()
#
#                    #print 'pres_top_WRF_i', pres_top_WRF_i
#                    #print 'pres_bottom_WRF_i', pres_bottom_WRF_i
#
#            #                        interp_hght = wrf.vinterp(nc, field=hght, vert_coord='pres', interp_levels=full_mb_p)
#            #                        interp_speed = wrf.vinterp(nc, field=speed, vert_coord='pres', interp_levels=full_mb_p)
#            #                        interp_pres = wrf.vinterp(nc, field=pres, vert_coord='pres', interp_levels=full_mb_p)
#            #                        interp_drct = wrf.vinterp(nc, field=drct, vert_coord='pres', interp_levels=full_mb_p)
#            #
#            #                        #print 'interp_drct', interp_drct.shape
#            #                        
#            #                        # narrow down speed, pres, and drct to one location and set pres boundaries
#            #                        hght = interp_hght[:,Station_xy[0],Station_xy[1]]
#            #                        speed = interp_speed[:,Station_xy[0],Station_xy[1]]
#            #                        pres = interp_pres[:,Station_xy[0],Station_xy[1]]
#            #                        drct = interp_drct[:,Station_xy[0],Station_xy[1]]
#
#
#            #                        # narrow down speed, pres, and drct to one location and set pres boundaries
#            #                        hght = hght[pres_bottom_WRF_i:pres_top_WRF_i,Station_xy[0],Station_xy[1]]
#            #                        speed = speed[pres_bottom_WRF_i:pres_top_WRF_i,Station_xy[0],Station_xy[1]]
#            #                        pres = pres[pres_bottom_WRF_i:pres_top_WRF_i,Station_xy[0],Station_xy[1]]
#            #                        print 'pres', pres
#            #                        drct = drct[pres_bottom_WRF_i:pres_top_WRF_i,Station_xy[0],Station_xy[1]]
#            #                    
#
#                    # narrow down speed, pres, and drct to one location and set pres boundaries
#                    hght = hght[:,Station_xy[0],Station_xy[1]]
#                    speed = speed[:,Station_xy[0],Station_xy[1]]
#                    v = v[:,Station_xy[0],Station_xy[1]]
#                    pres = pres[:,Station_xy[0],Station_xy[1]]
#                    drct = drct[:,Station_xy[0],Station_xy[1]]
#                    
#                    print 'pres', pres[0]
#
#                    #print 'pres', pres
#
#                    #print 'full_mb_p', full_mb_p
#
#                    f_hght = interp1d(pres,hght) 
#                    f_speed = interp1d(pres,speed)
#                    f_v = interp1d(pres,v)
#                    f_pres = interp1d(pres,pres)
#                    f_drct = interp1d(pres,drct)
#
#                    hght = f_hght(full_mb_p)
#                    speed = f_speed(full_mb_p)
#                    v = f_v(full_mb_p)
#                    pres = f_pres(full_mb_p)
#                    drct = f_drct(full_mb_p)
#
#                    # convert time to datetime to use as x coordinate
#                    length_pres = len(pres)
#                    time = [file_time]*length_pres
#
#                    # plot meridional wind speeds on ax where negaive is from north
#                    v_plot[:,total_count_WRF] = v
#
#                    #print 'v_plot shape 1', v_plot.shape
#
#                    total_count_WRF = total_count_WRF + 1
#
#                    #print 'total_count_WRF', total_count_WRF
#
#                    WRF_file_time_list.append(file_time)
#
#                    #print 'WRF_file_time_list', WRF_file_time_list
#
#                    #wind_speed = ax.scatter(time, pres, c=v, vmin=-40, vmax=40, cmap='seismic_r', zorder=1)
#
#                    # get variable from low pressure to 3000m (used to find max wind)
#                    height_cutoff_max = where(hght <= max_search_hgt)
#                    #print 'hght', hght
#                    #print 'max_search_hgt', max_search_hgt
#                    #print 'height_cutoff_max', height_cutoff_max
#                    hght_below_max = hght[height_cutoff_max]
#                    speed_below_max = speed[height_cutoff_max]
#                    pres_below_max = pres[height_cutoff_max]
#                    drct_below_max = drct[height_cutoff_max]
#                    #print 'speed', speed
#                    #print 'speed_below_max', speed_below_max
#                    # find max wind
#                    max_wind = amax(speed_below_max)
#                    #print 'max_wind', max_wind
#
#                    # find index of max wind
#                    max_wind_index_list = where(speed_below_max==max_wind)
#                    #print 'max_wind_index_list', max_wind_index_list
#
#                    # find max wind heights and at lowest height as defined in Oliviera et al 2018
#                    max_wind_height_list = hght_below_max[max_wind_index_list]
#                    max_wind_height = amin(max_wind_height_list)
#                    #print 'max_wind_height_list', max_wind_height_list
#
#                # find lowest height (highest pressure) with maximum wind as laid out in Oliviera et al 2018
#                    if max_wind > max_wind_threshold and max_wind_height < max_search_hgt:
#
#                        # find index of min height for max wind
#                        max_wind_height_index = where(hght == max_wind_height)
#
#                        # find pres for max wind
#                        max_wind_pres = pres_below_max[max_wind_height_index][0]
#
#                        # find most common direction around max wind speed
#                        # NOTE: If you make the layer too wide and peak is too close 
#                        drct_near_max_wind_height = drct_below_max[int(max_wind_height_index[0])-12:int(max_wind_height_index[0])+12]
#                        drct_mode, count = stats.mode(drct_near_max_wind_height)
#
#                        # avg wind around max wind speed
#                        speed_near_max_wind_height = speed_below_max[int(max_wind_height_index[0])-10:int(max_wind_height_index[0])+10]
#                        speed_mean = np.nanmean(speed_near_max_wind_height)
#                        #print 'hght', hght
#                        #print 'min_search_hgt', min_search_hgt
#                        #print 'max_wind_height', max_wind_height
#                        # find variables from max wind  up to 4000m (use to find if it is a peak)
#
#                        # function which returns True when constraints are satisfied.
#                        height_cutoff_funct = lambda hght: hght <=min_search_hgt and hght >= max_wind_height 
#
#                        # Apply constraints element-wise to the dists array.
#                        height_cutoff_min = np.vectorize(height_cutoff_funct)(hght) 
#                        height_cutoff_min = np.where(height_cutoff_min)[0] # Get output. 
#                        #print 'max height for min wind', hght[height_cutoff_min][-1]
#                        speed_min = speed[height_cutoff_min]
#                        #print 'speed_min', speed_min
#
#                        min_found = False
#
#                        while(min_found == False):
#
#                            try: 
#                                ind_first_min = argrelextrema(speed_min, np.less)[0][0]
#                                #print 'ind_first_min', ind_first_min
#                                current_min = speed_min[ind_first_min]
#                                #print 'current_min', current_min
#                                speed_min = speed_min[ind_first_min:] #get rid of everything under current min
#                                #print 'speed_min', speed_min
#                                ind_next_max = np.argmax(speed_min_)
#
#                                if(speed_min[ind_next_max] - current_min) > (max_wind - current_min) or len(speed_min) == 1:
#
#                                    min_found = True
#                                    #print 'min_found' 
#                            except: 
#                                ind_first_min = np.argmin(speed_min) # used if last value is min so there is no minima
#                                current_min = speed_min[ind_first_min]
#                                #print 'current_min', current_min
#                                break
#
#
#                        # find index for min wind
#                        min_wind_index = where(speed == current_min)
#
#                         # find pres for min wind
#                        min_wind_pres = pres[min_wind_index][0]
#                        #print 'min_wind_pres', min_wind_pres
#
#                        #print 'decrease_to_min', max_wind - current_min
#
#                        #print 'max - min', max_wind - current_min
#
#                        # see if difference between max wind and min wind is large enough
#                        if (max_wind - current_min) > decrease_to_min:
#
#                            # If yes a LLJ is considered to be found
#                            #print 'LLJ'
#
#                            if len(drct_mode) != 0:
#
#
#                                NE = 45
#                                NW = 315
#                                SE = 135
#                                SW = 225
#
#                                total_LLJ_count_WRF = total_LLJ_count_WRF + 1
#
#                                #print 'WRF LLJ count', total_LLJ_count_WRF
#
#                                WRF_found_time_list.append(file_time)
#
#                                WRF_speed_mean_list.append(speed_mean)
#
#                                WRF_max_wind_pres_list.append(max_wind_pres) 
#
#                                if drct_mode<=NE or drct_mode>=NW:
#                                    cc=colorN[crit_num-1]
#                                    ms=markerN[crit_num-1]
#                                    ss=sizeN[crit_num-1]
#                                    # plot arrow to indicate LLJ and pressure of max wind speed         
#                                    ax.scatter(dates.date2num(file_time), max_wind_pres, color=cc, marker=ms, s=ss, zorder=3)
#                                    ax2.scatter(dates.date2num(file_time), max_wind_pres, color=cc, marker=ms, s=ss, zorder=3)
#
#                                elif drct_mode>=SE and drct_mode<=SW:
#                                    cc=colorS[crit_num-1]
#                                    ms=markerS[crit_num-1]
#                                    ss=sizeS[crit_num-1]
#                                    # plot arrow to indicate LLJ and pressure of max wind speed         
#                                    ax.scatter(dates.date2num(file_time), max_wind_pres, color=cc, marker=ms, s=ss, zorder=3)
#                                    ax2.scatter(dates.date2num(file_time), max_wind_pres, color=cc, marker=ms, s=ss, zorder=3)
#
#                                else:
#                                    cc='black'
#
#                                    #total_count = total_count + 1
#                                    #print file_time
#                                    #max_wind_time_list.append(file_time)
#                                    #drct_mode_list.append(float(drct_mode[0]))
#                                    #speed_mean_list.append(speed_mean)
#                                    #max_wind_pres_list.append(max_wind_pres)
#                    else:
#                        # does not meet peak LLJ criteria 
#                        cc='gray'
#                else:
#                    # does not meet max wind LLJ criteria 
#                    cc='gray'
#
#    print 'v_plot shape 2', v_plot.shape  
#    print 'v_plot shape 2', v_plot.T.shape 
#
#    print 'WRF_file_time_list', WRF_file_time_list
#    print 'WRF_file_time_list shape', len(WRF_file_time_list)
#
#    print 'full_mb_p', full_mb_p
#    print 'full_mb_p shape', np.squeeze(full_mb_p.shape)
#
#    #plot v-wind
#    v_wind = ax.contourf(WRF_file_time_list,np.squeeze(full_mb_p),v_plot, cmap='seismic_r', levels=np.arange(-40,45,5), extend='both', zorder=1)
#
#    #v_wind = ax.contourf(WRF_file_time_list,np.squeeze(full_mb_p),c=v_plot, cmap='seismic_r', zorder=1)

    # set axes format

    ax.set_xlim([input_start_day, input_end_day + datetime.timedelta(days=1)])
    ax.set_ylim(ax.get_ylim()[::-1])
    #ax.set_xticks(WRF_file_time_list[::24], minor=False)
    #ax.set_xticks(WRF_file_time_list[::2], minor=True)
    ax.set_yticks(np.arange(max_pres,min_pres,-50), minor=False)
    ax.set_yticks(np.arange(max_pres,min_pres,-25), minor=True)
    ax.tick_params(axis='x', which='minor', bottom=True)
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=2))
    ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    ax.xaxis.set_minor_formatter(mdates.DateFormatter(""))
    ax.tick_params(axis='x', which='both', bottom=True, labelbottom=True, length=8, width=2)
    ax.tick_params(axis='y', which='both', length=8, width=2)
    #ax.set_xticklabels(time_list, rotation=20)

    # add color bars to wind speed and specific humidity plots
    cbar = fig.colorbar(diff_v_wind, ax=ax, orientation='vertical')

    cbar.set_label('obs-wrf v-wind diff (kts)', fontsize=30)
    cbar.ax.tick_params(labelsize=30) 

    # add labels to plot axes
    ax.set_ylabel('Pressure (hPa)', fontsize=30)   
    ax.set_xlabel('Time (MM-DD)', fontsize=30)


    #### used when adding GOES data ###

    #Nov 1-15
    #fig.savefig('%s/%s'%(outpath, 'f0GOES_' + short_station + 'vv_wspd_Nov01-Nov15_time_series.png'), bbox_inches='tight', dpi=120)

    #Nov 16-30
    #fig.savefig('%s/%s'%(outpath, 'f0GOES_' + short_station + 'vv_wspd_Nov16-Nov30_time_series.png'), bbox_inches='tight', dpi=120)

    #Dec 1-18
    #fig.savefig('%s/%s'%(outpath, 'f0GOES_' + short_station + 'vv_wspd_Dec01-Dec18_time_series.png'), bbox_inches='tight', dpi=120)

    #Nov 1 - Dec 18
    #fig.savefig('%s/%s'%(outpath, 'f000000_' + short_station + 'vv_wspd_Nov01-Dec18_time_series.png'), bbox_inches='tight', dpi=120)

    ##############################

    #Nov 1-16
    #fig.savefig('%s/%s'%(outpath, 'Diff_w_jet_core_' + short_station + 'vv_wspd_Nov01-Nov16_time_series.png'), bbox_inches='tight', dpi=120)

    #Nov 17- Dec 2
    #fig.savefig('%s/%s'%(outpath, 'Diff_w_jet_core_' + short_station + 'vv_wspd_Nov17-Dec02_time_series.png'), bbox_inches='tight', dpi=120)

    #Dec 3-18
    fig.savefig('%s/%s'%(outpath, 'Diff_w_jet_core_' + short_station + 'vv_wspd_Dec03-Dec18_time_series.png'), bbox_inches='tight', dpi=120)

    #Nov 1-15
    #fig.savefig('%s/%s'%(outpath, 'f0000_' + short_station + 'vv_wspd_Nov01-Nov15_time_series.png'), bbox_inches='tight', dpi=120)
    #fig2.savefig('%s/%s'%(outpath, 'f0_' + short_station + 'vv_q_Nov01-Nov15_time_series.png'), bbox_inches='tight', dpi=120)
    #fig3.savefig('%s/%s'%(outpath, 'f0_' + short_station + 'vv_barb_Nov01-Nov15_time_series.png'), bbox_inches='tight', dpi=120)

    #Nov 16-30
    #fig.savefig('%s/%s'%(outpath, 'f0000_' + short_station + 'vv_wspd_Nov16-Nov30_time_series.png'), bbox_inches='tight', dpi=120)
    #fig2.savefig('%s/%s'%(outpath, 'f0_' + short_station + 'vv_q_Nov16-Nov30_time_series.png'), bbox_inches='tight', dpi=120)
    #fig3.savefig('%s/%s'%(outpath, 'f0_' + short_station + 'vv_barb_Nov16-Nov30_time_series.png'), bbox_inches='tight', dpi=120)

    #Dec 1-18
    #fig.savefig('%s/%s'%(outpath, 'f0000_' + short_station + 'vv_wspd_Dec01-Dec18_time_series.png'), bbox_inches='tight', dpi=120)
    #fig2.savefig('%s/%s'%(outpath, 'f0_' + short_station + 'vv_q_Dec01-Dec18_time_series.png'), bbox_inches='tight', dpi=120)
    #fig3.savefig('%s/%s'%(outpath, 'f0_' + short_station + 'vv_barb_Dec01-Dec18_time_series.png'), bbox_inches='tight', dpi=120)      

    ##Dec 5-14
    #fig.savefig('%s/%s'%(outpath, short_station + 'vv_wspd_Dec05-Dec14_time_series.png'), bbox_inches='tight', dpi=320)
    #fig2.savefig('%s/%s'%(outpath, short_station + 'vv_q_Dec05-Dec14_time_series.png'), bbox_inches='tight', dpi=120)
    #fig3.savefig('%s/%s'%(outpath, short_station + 'vv_barb_Dec05-Dec14_time_series.png'), bbox_inches='tight', dpi=120)

    ##Dec 11-14
    #fig.savefig('%s/%s'%(outpath, short_station + 'vv_wspd_Dec11-Dec14_time_series.png'), bbox_inches='tight', dpi=320)
    #fig2.savefig('%s/%s'%(outpath, short_station + 'vv_q_Dec11-Dec14_time_series.png'), bbox_inches='tight', dpi=120)
    #fig3.savefig('%s/%s'%(outpath, short_station + 'vv_barb_Dec11-Dec14_time_series.png'), bbox_inches='tight', dpi=120)

    ##Nov 9-13
    #fig.savefig('%s/%s'%(outpath, short_station + 'vv_wspd_Nov09-Nov13_time_series.png'), bbox_inches='tight', dpi=320)
    #fig2.savefig('%s/%s'%(outpath, short_station + 'vv_q_Nov09-Nov13_time_series.png'), bbox_inches='tight', dpi=120)
    #fig3.savefig('%s/%s'%(outpath, short_station + 'v_barb_Nov09-Nov13_time_series.png'), bbox_inches='tight', dpi=120)

    ##Nov 14-24
    #fig.savefig('%s/%s'%(outpath, short_station + 'vv_wspd_Nov14-Nov24_time_series.png'), bbox_inches='tight', dpi=320)
    #fig2.savefig('%s/%s'%(outpath, short_station + 'vv_q_Nov14-Nov24_time_series.png'), bbox_inches='tight', dpi=120)
    #fig3.savefig('%s/%s'%(outpath, short_station + 'vv_barb_Nov14-Nov24_time_series.png'), bbox_inches='tight', dpi=120)

    ##Nov 25 - Nov 30
    #fig.savefig('%s/%s'%(outpath, short_station + 'vv_wspd_Nov25-Nov30_time_series.png'), bbox_inches='tight', dpi=320)
    #fig2.savefig('%s/%s'%(outpath, short_station + 'vv_q_Nov25-Nov30_time_series.png'), bbox_inches='tight', dpi=120)
    #fig3.savefig('%s/%s'%(outpath, short_station + 'vv_barb_Nov25-Nov30_time_series.png'), bbox_inches='tight', dpi=120)

    ##Dec 1 - Dec 4
    #fig.savefig('%s/%s'%(outpath, short_station + 'vv_wspd_Dec01-Dec04_time_series.png'), bbox_inches='tight', dpi=320)
    #fig2.savefig('%s/%s'%(outpath, short_station + 'vv_q_Dec01-Dec04_time_series.png'), bbox_inches='tight', dpi=120)
    #fig3.savefig('%s/%s'%(outpath, short_station + 'vv_barb_Dec01-Dec04_time_series.png'), bbox_inches='tight', dpi=120)

    ##Dec 14- Dec 20
    #fig.savefig('%s/%s'%(outpath, short_station + 'vv_wspd_Dec14-Dec20_time_series.png'), bbox_inches='tight', dpi=320)
    #fig2.savefig('%s/%s'%(outpath, short_station + 'vv_q_Dec14-Dec20_time_series.png'), bbox_inches='tight', dpi=120)
    #fig3.savefig('%s/%s'%(outpath, short_station + 'vv_barb_Dec14-Dec20_time_series.png'), bbox_inches='tight', dpi=120)

    plt.close('all')

    #print '# of all soundings: ', len(all_header_lines)
    #print '# of LLJ soundings: ', len(TimeList1)
    #print 'crit 1', TimeList1
    #print '# of LLJ soundings: ', len(TimeList2)
    #print 'crit 2', TimeList2
