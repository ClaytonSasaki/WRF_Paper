import matplotlib
matplotlib.use('Agg') 

from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as crs
import numpy as np
import numpy.ma as ma
from numpy import exp,where,ma,cos,sin,pi,amax,amin
import pandas as pd
import os
from scipy import stats
from scipy.interpolate import interp1d
import math
from cartopy.feature import NaturalEarthFeature, LAND, OCEAN, COASTLINE, BORDERS, LAKES, RIVERS
import time

# ------------------------- cmaps ------------------------------

from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, interplevel, ll_to_xy_proj, get_basemap, xy_to_ll)

import matplotlib.colors as colors

start_time = time.time()

# used for contourf of terrain 
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval), cmap(np.linspace(minval, maxval, n)))
    return new_cmap

arr = np.linspace(0, 50, 100).reshape((10, 10))
fig, ax = plt.subplots(ncols=2)

cmap = plt.get_cmap('terrain')
new_cmap = truncate_colormap(cmap, 0.21, 1)
ax[0].imshow(arr, interpolation='nearest', cmap=cmap)
ax[1].imshow(arr, interpolation='nearest', cmap=new_cmap)

# used for pcolormesh of wind
#def truncate_colormap2(cmap, minval=0.0, maxval=1.0, n=100):
#    new_cmap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval), cmap(np.linspace(minval, maxval, n)))
#    return new_cmap
#
#arr = np.linspace(0, 50, 100).reshape((10, 10))
#fig, ax = plt.subplots(ncols=2)
#
#cmap2 = plt.get_cmap('gray')
#new_cmap2 = truncate_colormap(cmap2, 0.2, 0.8)
#ax[0].imshow(arr, interpolation='nearest', cmap=cmap2)
#ax[1].imshow(arr, interpolation='nearest', cmap=new_cmap2)

# used for contour of wind
def truncate_colormap2(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval), cmap(np.linspace(minval, maxval, n)))
    return new_cmap

arr = np.linspace(0, 50, 100).reshape((10, 10))
fig, ax = plt.subplots(ncols=2)

cmap2 = plt.get_cmap('gnuplot2')
new_cmap2 = truncate_colormap(cmap2, 0.2, 1.0)
ax[0].imshow(arr, interpolation='nearest', cmap=cmap2)
ax[1].imshow(arr, interpolation='nearest', cmap=new_cmap2)

# ----------------- reading in file/plotting terrain -----------------

# Open the NetCDF file
path_wrf = '/home/disk/monsoon/relampago/raw/wrf/'

#chosen_date_time_list = ['2018-11-02_07_00_00']

#chosen_date_time_list = ['2018-11-02_03_00_00', '2018-11-02_04_00_00', '2018-11-02_05_00_00', '2018-11-02_06_00_00', '2018-11-02_07_00_00', '2018-11-02_08_00_00', '2018-11-02_09_00_00', '2018-11-02_10_00_00', '2018-11-02_11_00_00', '2018-11-02_12_00_00', '2018-11-02_13_00_00', '2018-11-02_14_00_00', '2018-11-02_15_00_00', '2018-11-02_16_00_00', '2018-11-02_17_00_00', '2018-11-02_18_00_00', '2018-11-02_19_00_00', '2018-11-02_20_00_00', '2018-11-02_21_00_00', '2018-11-02_22_00_00','2018-11-02_23_00_00', '2018-11-03_00_00_00']

#chosen_date_time_list = ['2018-11-02_00_00_00', '2018-11-02_01_00_00', '2018-11-02_02_00_00', '2018-11-02_03_00_00', '2018-11-02_04_00_00', '2018-11-02_05_00_00', '2018-11-02_06_00_00', '2018-11-02_07_00_00', '2018-11-02_08_00_00', '2018-11-02_09_00_00', '2018-11-02_10_00_00', '2018-11-02_11_00_00', '2018-11-02_12_00_00', '2018-11-02_13_00_00', '2018-11-02_14_00_00', '2018-11-02_15_00_00', '2018-11-02_16_00_00', '2018-11-02_17_00_00', '2018-11-02_18_00_00', '2018-11-02_19_00_00', '2018-11-02_20_00_00', '2018-11-02_21_00_00', '2018-11-02_22_00_00','2018-11-02_23_00_00', '2018-11-03_00_00_00']

#chosen_date_time_list = ['2018-11-08_00_00_00', '2018-11-08_01_00_00', '2018-11-08_02_00_00', '2018-11-08_03_00_00', '2018-11-08_04_00_00', '2018-11-08_05_00_00', '2018-11-08_06_00_00', '2018-11-08_07_00_00', '2018-11-08_08_00_00', '2018-11-08_09_00_00', '2018-11-08_10_00_00', '2018-11-08_11_00_00', '2018-11-08_12_00_00', '2018-11-08_13_00_00', '2018-11-08_14_00_00', '2018-11-08_15_00_00', '2018-11-08_16_00_00', '2018-11-08_17_00_00', '2018-11-08_18_00_00', '2018-11-08_19_00_00', '2018-11-08_20_00_00', '2018-11-08_21_00_00', '2018-11-08_22_00_00','2018-11-08_23_00_00', '2018-11-09_00_00_00']

#chosen_date_time_list = ['2018-11-09_00_00_00', '2018-11-09_01_00_00', '2018-11-09_02_00_00', '2018-11-09_03_00_00', '2018-11-09_04_00_00', '2018-11-09_05_00_00', '2018-11-09_06_00_00', '2018-11-09_07_00_00', '2018-11-09_08_00_00', '2018-11-09_09_00_00', '2018-11-09_10_00_00', '2018-11-09_11_00_00', '2018-11-09_12_00_00', '2018-11-09_13_00_00', '2018-11-09_14_00_00', '2018-11-09_15_00_00', '2018-11-09_16_00_00', '2018-11-09_17_00_00', '2018-11-09_18_00_00', '2018-11-09_19_00_00', '2018-11-09_20_00_00', '2018-11-09_21_00_00', '2018-11-09_22_00_00','2018-11-09_23_00_00', '2018-11-10_00_00_00']

#chosen_date_time_list = ['2018-11-10_00_00_00', '2018-11-10_01_00_00', '2018-11-10_02_00_00', '2018-11-10_03_00_00', '2018-11-10_04_00_00', '2018-11-10_05_00_00', '2018-11-10_06_00_00', '2018-11-10_07_00_00', '2018-11-10_08_00_00', '2018-11-10_09_00_00', '2018-11-10_10_00_00', '2018-11-10_11_00_00', '2018-11-10_12_00_00', '2018-11-10_13_00_00', '2018-11-10_14_00_00', '2018-11-10_15_00_00', '2018-11-10_16_00_00', '2018-11-10_17_00_00', '2018-11-10_18_00_00', '2018-11-10_19_00_00', '2018-11-10_20_00_00', '2018-11-10_21_00_00', '2018-11-10_22_00_00','2018-11-10_23_00_00', '2018-11-11_00_00_00','2018-11-11_01_00_00', '2018-11-11_02_00_00', '2018-11-11_03_00_00', '2018-11-11_04_00_00', '2018-11-11_05_00_00', '2018-11-11_06_00_00', '2018-11-11_07_00_00', '2018-11-11_08_00_00', '2018-11-11_09_00_00', '2018-11-11_10_00_00', '2018-11-11_11_00_00', '2018-11-11_12_00_00', '2018-11-11_13_00_00', '2018-11-11_14_00_00', '2018-11-11_15_00_00', '2018-11-11_16_00_00', '2018-11-11_17_00_00', '2018-11-11_18_00_00', '2018-11-11_19_00_00', '2018-11-11_20_00_00', '2018-11-11_21_00_00', '2018-11-11_22_00_00','2018-11-11_23_00_00', '2018-11-12_00_00_00', '2018-11-12_01_00_00', '2018-11-12_02_00_00', '2018-11-12_03_00_00', '2018-11-12_04_00_00', '2018-11-12_05_00_00', '2018-11-12_06_00_00', '2018-11-12_07_00_00', '2018-11-12_08_00_00', '2018-11-12_09_00_00', '2018-11-12_10_00_00', '2018-11-12_11_00_00', '2018-11-12_12_00_00', '2018-11-12_13_00_00', '2018-11-12_14_00_00', '2018-11-12_15_00_00', '2018-11-12_16_00_00', '2018-11-12_17_00_00', '2018-11-12_18_00_00', '2018-11-12_19_00_00', '2018-11-12_20_00_00', '2018-11-12_21_00_00', '2018-11-12_22_00_00','2018-11-12_23_00_00', '2018-11-13_00_00_00']

#chosen_date_time_list = ['2018-11-11_00_00_00', '2018-11-11_01_00_00', '2018-11-11_02_00_00', '2018-11-11_03_00_00', '2018-11-11_04_00_00', '2018-11-11_05_00_00', '2018-11-11_06_00_00', '2018-11-11_07_00_00', '2018-11-11_08_00_00', '2018-11-11_09_00_00', '2018-11-11_10_00_00', '2018-11-11_11_00_00', '2018-11-11_12_00_00', '2018-11-11_13_00_00', '2018-11-11_14_00_00', '2018-11-11_15_00_00', '2018-11-11_16_00_00', '2018-11-11_17_00_00', '2018-11-11_18_00_00', '2018-11-11_19_00_00', '2018-11-11_20_00_00', '2018-11-11_21_00_00', '2018-11-11_22_00_00','2018-11-11_23_00_00', '2018-11-12_00_00_00']

#chosen_date_time_list = ['2018-11-12_00_00_00', '2018-11-12_01_00_00', '2018-11-12_02_00_00', '2018-11-12_03_00_00', '2018-11-12_04_00_00', '2018-11-12_05_00_00', '2018-11-12_06_00_00', '2018-11-12_07_00_00', '2018-11-12_08_00_00', '2018-11-12_09_00_00', '2018-11-12_10_00_00', '2018-11-12_11_00_00', '2018-11-12_12_00_00', '2018-11-12_13_00_00', '2018-11-12_14_00_00', '2018-11-12_15_00_00', '2018-11-12_16_00_00', '2018-11-12_17_00_00', '2018-11-12_18_00_00', '2018-11-12_19_00_00', '2018-11-12_20_00_00', '2018-11-12_21_00_00', '2018-11-12_22_00_00','2018-11-12_23_00_00', '2018-11-13_00_00_00']

#chosen_date_time_list = ['2018-11-20_00_00_00', '2018-11-20_01_00_00', '2018-11-20_02_00_00', '2018-11-20_03_00_00', '2018-11-20_04_00_00', '2018-11-20_05_00_00', '2018-11-20_06_00_00', '2018-11-20_07_00_00', '2018-11-20_08_00_00', '2018-11-20_09_00_00', '2018-11-20_10_00_00', '2018-11-20_11_00_00', '2018-11-20_12_00_00', '2018-11-20_13_00_00', '2018-11-20_14_00_00', '2018-11-20_15_00_00', '2018-11-20_16_00_00', '2018-11-20_17_00_00', '2018-11-20_18_00_00', '2018-11-20_19_00_00', '2018-11-20_20_00_00', '2018-11-20_21_00_00', '2018-11-20_22_00_00', '2018-11-20_23_00_00', '2018-11-21_00_00_00', '2018-11-21_01_00_00', '2018-11-21_02_00_00', '2018-11-21_03_00_00', '2018-11-21_04_00_00', '2018-11-21_05_00_00', '2018-11-21_06_00_00', '2018-11-21_07_00_00', '2018-11-21_08_00_00', '2018-11-21_09_00_00', '2018-11-21_10_00_00', '2018-11-21_11_00_00', '2018-11-21_12_00_00', '2018-11-21_13_00_00', '2018-11-21_14_00_00', '2018-11-21_15_00_00', '2018-11-21_16_00_00', '2018-11-21_17_00_00', '2018-11-21_18_00_00', '2018-11-21_19_00_00', '2018-11-21_20_00_00', '2018-11-21_21_00_00', '2018-11-21_22_00_00', '2018-11-21_23_00_00', '2018-11-22_00_00_00']


#chosen_date_time_list = ['2018-12-08_00_00_00', '2018-12-08_01_00_00', '2018-12-08_02_00_00', '2018-12-08_03_00_00', '2018-12-08_04_00_00', '2018-12-08_05_00_00', '2018-12-08_06_00_00', '2018-12-08_07_00_00', '2018-12-08_08_00_00', '2018-12-08_09_00_00', '2018-12-08_10_00_00', '2018-12-08_11_00_00', '2018-12-08_12_00_00', '2018-12-08_13_00_00', '2018-12-08_14_00_00', '2018-12-08_15_00_00', '2018-12-08_16_00_00', '2018-12-08_17_00_00', '2018-12-08_18_00_00', '2018-12-08_19_00_00', '2018-12-08_20_00_00', '2018-12-08_21_00_00', '2018-12-08_22_00_00','2018-12-08_23_00_00', '2018-12-09_00_00_00']

#chosen_date_time_list = ['2018-12-13_00_00_00', '2018-12-13_01_00_00', '2018-12-13_02_00_00', '2018-12-13_03_00_00', '2018-12-13_04_00_00', '2018-12-13_05_00_00', '2018-12-13_06_00_00', '2018-12-13_07_00_00', '2018-12-13_08_00_00', '2018-12-13_09_00_00', '2018-12-13_10_00_00', '2018-12-13_11_00_00', '2018-12-13_12_00_00', '2018-12-13_13_00_00', '2018-12-13_14_00_00', '2018-12-13_15_00_00', '2018-12-13_16_00_00', '2018-12-13_17_00_00', '2018-12-13_18_00_00', '2018-12-13_19_00_00', '2018-12-13_20_00_00', '2018-12-13_21_00_00', '2018-12-13_22_00_00','2018-12-13_23_00_00', '2018-12-14_00_00_00']

chosen_date_time_list = ['2018-12-13_02_00_00', '2018-12-13_03_00_00', '2018-12-13_04_00_00', '2018-12-13_05_00_00', '2018-12-13_06_00_00', '2018-12-13_07_00_00', '2018-12-13_08_00_00', '2018-12-13_09_00_00', '2018-12-13_10_00_00', '2018-12-13_11_00_00', '2018-12-13_12_00_00', '2018-12-13_13_00_00', '2018-12-13_14_00_00', '2018-12-13_15_00_00', '2018-12-13_16_00_00', '2018-12-13_17_00_00', '2018-12-13_18_00_00', '2018-12-13_19_00_00', '2018-12-13_20_00_00', '2018-12-13_21_00_00', '2018-12-13_22_00_00','2018-12-13_23_00_00', '2018-12-14_00_00_00', '2018-12-12_00_00_00', '2018-12-12_01_00_00', '2018-12-12_02_00_00', '2018-12-12_03_00_00', '2018-12-12_04_00_00', '2018-12-12_05_00_00', '2018-12-12_06_00_00', '2018-12-12_07_00_00', '2018-12-12_08_00_00', '2018-12-12_09_00_00', '2018-12-12_10_00_00', '2018-12-12_11_00_00', '2018-12-12_12_00_00', '2018-12-12_13_00_00', '2018-12-12_14_00_00', '2018-12-12_15_00_00', '2018-12-12_16_00_00', '2018-12-12_17_00_00', '2018-12-12_18_00_00', '2018-12-12_19_00_00', '2018-12-12_20_00_00', '2018-12-12_21_00_00', '2018-12-12_22_00_00','2018-12-12_23_00_00', '2018-12-11_00_00_00', '2018-12-11_01_00_00', '2018-12-11_02_00_00', '2018-12-11_03_00_00', '2018-12-11_04_00_00', '2018-12-11_05_00_00', '2018-12-11_06_00_00', '2018-12-11_07_00_00', '2018-12-11_08_00_00', '2018-12-11_09_00_00', '2018-12-11_10_00_00', '2018-12-11_11_00_00', '2018-12-11_12_00_00', '2018-12-11_13_00_00', '2018-12-11_14_00_00', '2018-12-11_15_00_00', '2018-12-11_16_00_00', '2018-12-11_17_00_00', '2018-12-11_18_00_00', '2018-12-11_19_00_00', '2018-12-11_20_00_00', '2018-12-11_21_00_00', '2018-12-11_22_00_00','2018-12-11_23_00_00']
########### for finding LLJs ###########

# 2: high - allows for jets centered higher in the atmosphere, and a bit weaker compared to mean flow
crit = [2 ,3600, 6000, 23.326, 11.663]

# min and max pressure to plot
min_pres = 450
max_pres = 975

crit_num = crit[0]
max_search_hgt = crit[1]
min_search_hgt = crit[2]
max_wind_threshold = crit[3]
decrease_to_min = crit[4]

full_mb_p = np.arange(max_pres,min_pres-1,-1)

######################################

#print '---header: %s seconds ---' %(time.time() - start_time)
#start_time2 = time.time()

for chosen_date_time in chosen_date_time_list:
    
    chosen_date_dt = pd.to_datetime(chosen_date_time[:10], format='%Y-%m-%d', errors='ignore')
    
    chosen_hour = int(chosen_date_time[11:13])
    
    sorted_folders = sorted(os.listdir(path_wrf))
    
    #print 'sorted_folders', sorted_folders
    
    # get the correct file to match the given date_time
    for folder_date in sorted_folders:

        #print 'folder_date', folder_date

        if folder_date.startswith('20'):

            #print 'yes'

            folder_path = os.path.join(path_wrf, folder_date)

            folder_date_dt = pd.to_datetime(folder_date, format='%Y%m%d', errors='ignore')

            #print 'folder_date_dt', folder_date_dt
            #print 'chosen_date_dt', chosen_date_dt

            if folder_date_dt == chosen_date_dt:

                #print 'found'

                sorted_hour_files = sorted(os.listdir(folder_path))

                for hourly_file in sorted_hour_files:

                    file_path = os.path.join(folder_path, hourly_file)

                    #print 'hourly_file', hourly_file

                    file_hour = int(hourly_file[22:24])

                    #print 'chosen_hour', chosen_hour
                    #print type(chosen_hour)
                    #print 'file_hour', file_hour
                    #print type(file_hour)
                    #print ' '

                    if file_hour == chosen_hour:

                        #print 'found'

                        ncfile = Dataset(file_path,'r')
                        chosen_file_name = hourly_file

    #print '---get file: %s seconds ---' %(time.time() - start_time2)
    #start_time3 = time.time()
    # Get the terrain height
    terrain = getvar(ncfile, 'HGT')

    # Get the latitude and longitude points
    lats, lons = latlon_coords(terrain)

    # Get the cartopy mapping object
    cart_proj = get_cartopy(terrain)

    # Create a figure
    fig = plt.figure(figsize=(12,6))
    # Set the GeoAxes to the projection used by WRF
    ax = plt.axes(projection=cart_proj)

    #ax.add_feature(LAND)
    ax.add_feature(OCEAN, zorder=2)
    #ax.add_feature(COASTLINE)
    ax.add_feature(LAKES, alpha=0.5, zorder=2)
    ax.add_feature(RIVERS, zorder=2)
    ax.add_feature(BORDERS, edgecolor='gray', zorder=2)

    # Make  filled contours for the terrain
    terrain_plot = plt.contourf(to_np(lons), to_np(lats), to_np(terrain), levels=[0,250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3250,3500,3750,4000,4250,4500,4750,5000,5250,5500,5750,6000,6250], extend='max', transform=crs.PlateCarree(), cmap=new_cmap, zorder=1)

    # Add a color bar
    cbar = plt.colorbar(terrain_plot, ax=ax, shrink=.98)
    cbar.set_label('Elevation (m)')

    # Set the map bounds
    ax.set_xlim(cartopy_xlim(terrain))
    ax.set_ylim(cartopy_ylim(terrain))

    # Add the gridlines
    #ax.gridlines(color="gray", linestyle="dotted", draw_labels=True)
    ax.gridlines(crs=crs.PlateCarree(), linewidth=1, color='gray', alpha=0.5, linestyle='--')
    
    ############################# read in some other variables #############################
    
    hght = getvar(ncfile, 'height', units='m') # in m NOTE: this is MSL not AGL!!!
    pres = getvar(ncfile, 'pressure')

    u = getvar(ncfile, 'ua', units='m s-1') # in m/s
    v = getvar(ncfile, 'va', units='m s-1') # in m/s
    speed, drct = getvar(ncfile, 'wspd_wdir', units='m s-1') # in m/s
    
    # ------------------------- plot winds and height ------------------------------
    level = 930
    
    #mark out areas above level specifed, can be used to see if area is included in analysis at a certain height
    #sfc_pres = getvar(ncfile, 'PSFC')/100 #in hPa
    #sfc_pres.values[sfc_pres.values <= level] = -99
    #sfc_pres.values[sfc_pres.values > level] = np.nan
    #
    #invalid = plt.contourf(to_np(lons), to_np(lats), to_np(sfc_pres), transform=crs.PlateCarree(), cmap='gray', zorder=2)
    #
    #plt.colorbar(invalid, ax=ax)

    u_level = np.squeeze(interplevel(u, pres, level))
    v_level = np.squeeze(interplevel(v, pres, level))
#    speed_level = np.squeeze(interplevel(speed, pres, level))

    #speed_level.values[speed_level.values < 12] = np.nan
#    smooth_speed_level = smooth2d(speed_level, 10)

    ##get area of v-wind >12 m/s 
    #v_level_signFlip = -np.squeeze(v_level)
#    levels=[12,13,14,15,16,17,18,19,20,21,22,23,24,25,26]
    #norm = colors.BoundaryNorm(levels, ncolors=new_cmap2.N, clip=True)
    #grt12 = plt.pcolormesh(to_np(lons), to_np(lats), to_np(smooth_speed_level), transform=crs.PlateCarree(), cmap=new_cmap2, norm=norm, zorder=3)
    
#    levels=[12]
#    grt12 = plt.contour(to_np(lons), to_np(lats), to_np(smooth_speed_level), levels=levels, transform=crs.PlateCarree(), cmap=new_cmap2, zorder=3)
#    ax.clabel(grt12, inline=1, fontsize=10, fmt='%.0f')

    # Add a color bar
    #cbar2 = plt.colorbar(grt12, ax=ax)
    #cbar2.set_label('full wind (m/s)')
    
    ####################### get points that meet LLJ criteria #######################
    
    #print 'terrain shape', terrain.shape
    
    df_barb_data = pd.DataFrame(columns = ['lat', 'lon', 'u_max', 'v_max', 'color'])

    u_max = to_np(terrain)
    u_max[:,:] = 0
    v_max = to_np(terrain)
    v_max[:,:] = 0

    print u_max
    print to_np(u_max).shape
    print v_max
    print to_np(v_max).shape

    #print '---get data, plot terrain: %s seconds ---' %(time.time() - start_time3)
    
    # go through each point on the grid (skips the western most 100 columns, over Chile and the Andes, and the southern most 40 rows to expediate code)
    for y in range(40,len(terrain[:,-1])):
        for x in range(100,len(terrain[-1,:])):

            #start_time4 = time.time()

            LLJ_found = False
            
            print 'x', x
            print 'y', y
            
            # narrow down speed, pres, and drct to one location and set pres boundaries
            hght_i = hght[:,y,x]
            speed_i = speed[:,y,x]
            v_i = v[:,y,x]
            pres_i = pres[:,y,x]
            drct_i = drct[:,y,x]

            #print '---get point data: %s seconds ---' %(time.time() - start_time4)
            #start_time5 = time.time()

            #print 'pres', pres[0]

            #print 'pres', pres

            #print 'full_mb_p', full_mb_p

            f_hght = interp1d(pres_i,hght_i) 
            f_speed = interp1d(pres_i,speed_i)
            f_v = interp1d(pres_i,v_i)
            f_pres = interp1d(pres_i,pres_i)
            f_drct = interp1d(pres_i,drct_i)

            # the highest pressure is less than max_pres we need cannot interpolate down to max_pres
            if(pres_i[0]<=max_pres):
                part_mb_p = np.arange(math.floor(pres_i[0]),min_pres-1,-1)
                max_plot_pres = math.floor(pres_i[0])
            else:
                part_mb_p = full_mb_p
                max_plot_pres = max_pres

            #print 'math.floor(pres[0])', math.floor(pres[0])
            #print 'max_plot_pres', max_plot_pres

            #print 'len(full_mb_p)', len(full_mb_p)
            #print 'len(part_mb_p)', len(part_mb_p)

            #print 'full_mb_p[0]', full_mb_p[0]
            #print 'part_mb_p[0]', part_mb_p[0]

            hght_i = f_hght(part_mb_p)
            speed_i = f_speed(part_mb_p)
            v_i = f_v(part_mb_p)
            pres_i = f_pres(part_mb_p)
            drct_i = f_drct(part_mb_p)
            
            #print '---interpolate point data: %s seconds ---' %(time.time() - start_time5)
            #start_time6 = time.time()

            #print 'to_np(hght_i[0])', to_np(hght_i)[0]
            #print 'max_search_hgt-10', max_search_hgt-10
            
            # check if below max search hght is above the terrain? Sea level presure?
            if(to_np(hght_i)[0] <= max_search_hgt-10):

                # get variable from low pressure to 3000m (used to find max wind)
                height_cutoff_max = where(hght_i <= max_search_hgt)
                #print 'hght_i', hght_i
                #print 'max_search_hgt', max_search_hgt
                #print 'height_cutoff_max', height_cutoff_max
                hght_below_max = hght_i[height_cutoff_max]
                speed_below_max = speed_i[height_cutoff_max]
                pres_below_max = pres_i[height_cutoff_max]
                drct_below_max = drct_i[height_cutoff_max]
                #print 'speed', speed
                #print 'speed_i'
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

                #print '---find max: %s seconds ---' %(time.time() - start_time6)
                #start_time7 = time.time()

            # find lowest height (highest pressure) with maximum wind as laid out in Oliviera et al 2018
                if max_wind > max_wind_threshold and max_wind_height < max_search_hgt:

                    # find index of min height for max wind
                    max_wind_height_index = where(hght_i == max_wind_height)

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
                    height_cutoff_funct = lambda hght_i: hght_i <=min_search_hgt and hght_i >= max_wind_height 

                    # Apply constraints element-wise to the dists array.
                    height_cutoff_min = np.vectorize(height_cutoff_funct)(hght_i) 
                    height_cutoff_min = np.where(height_cutoff_min)[0] # Get output. 
                    #print 'max height for min wind', hght[height_cutoff_min][-1]
                    speed_min = speed_i[height_cutoff_min]
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

                    #print '---find min: %s seconds ---' %(time.time() - start_time7)
                    #start_time8 = time.time()

                    # find index for min wind
                    min_wind_index = where(speed_i == current_min)

                     # find pres for min wind
                    min_wind_pres = pres_i[min_wind_index][0]
                    #print 'min_wind_pres', min_wind_pres

                    #print 'decrease_to_min', max_wind - current_min

                    #print 'max - min', max_wind - current_min

                    # see if difference between max wind and min wind is large enough
                    if (max_wind - current_min) > decrease_to_min:

                        # If yes a LLJ is considered to be found
                        #print 'LLJ'

                        if len(drct_mode) != 0:

                            NE = 45
                            NW = 315
                            SE = 135
                            SW = 225
                            
                            # plot dot at point that meets LLJ criteria, red for northerly, blue for southerly
                            if drct_mode<=NE or drct_mode>=NW:
                                cc='red'

                            elif drct_mode>=SE and drct_mode<=SW:
                                cc='blue'

                            else:
                                cc='black'

                            lon_lat = to_np(xy_to_ll(ncfile, x, y))
                            
                            lat = lon_lat[0]
                            lon = lon_lat[1]
                            
                            if max_wind_pres <= 700:
                                height_index = 'gold'
                            elif max_wind_pres > 700 and max_wind_pres <=750:
                                height_index = 'orange'
                            elif max_wind_pres > 750 and max_wind_pres <=800:
                                height_index = 'orangered'
                            elif max_wind_pres > 800 and max_wind_pres <=850:
                                height_index = 'mediumvioletred'
                            elif max_wind_pres > 850 and max_wind_pres <=900:
                                height_index = 'darkviolet'
                            elif max_wind_pres > 900:
                                height_index = 'mediumblue'
                            
                            plt.scatter(lon, lat, s=15, color=height_index, zorder=4, transform=crs.PlateCarree())
                            
                            # get winds at max speed and plot barbs for those winds
                            u_max_level = np.squeeze(interplevel(u, pres, max_wind_pres))
                            v_max_level = np.squeeze(interplevel(v, pres, max_wind_pres))
                            
                            #print lon
                            #print lat
                            #print to_np(u_max_level[y,x])
                            #print to_np(v_max_level[y,x])
                            u_max[y,x] = to_np(u_max_level[y,x])
                            #print 'u_max[y,x]', u_max[y,x]
                            v_max[y,x] = to_np(v_max_level[y,x])
                            
                            df_barb_data = df_barb_data.append({'lat' : lon, 'lon' : lat, 'u_max' : to_np(u_max_level[y,x]), 'v_max' : to_np(v_max_level[y,x]), 'color' : cc}, ignore_index = True)
                            #print 'df_barb_data', df_barb_data
                            #plt.barbs(lon, lat, to_np(u_max_level[y,x]), to_np(v_max_level[y,x]), transform=crs.PlateCarree(), length=5, zorder=6, color=cc)
                            #print 'met criteria'
                            
                            LLJ_found = True

                            #print '---plot LLJ: %s seconds ---' %(time.time() - start_time8)

                    #else:
                        # does not meet peak LLJ criteria 

                #else:
                    # does not meet max wind LLJ criteria
            if LLJ_found == False:
                lon_lat = to_np(xy_to_ll(ncfile, x, y))
                            
                lat = lon_lat[0]
                lon = lon_lat[1]
                
                df_barb_data = df_barb_data.append({'lat' : lon, 'lon' : lat, 'u_max' : 0, 'v_max' : 0, 'color' : 'black'}, ignore_index = True)
                #print 'df_barb_data', df_barb_data

                #print '---plot non-LLJ: %s seconds ---' %(time.time() - start_time6)

#    print 'df_barb_data', df_barb_data
#    print 'df_barb_data.u_max', df_barb_data.u_max
#    print 'type(df_barb_data.u_max)', type(df_barb_data.u_max)
#    print 'np.array(df_barb_data.u_max)', np.array(df_barb_data.u_max)
#    print 'type(np.array(df_barb_data.u_max))', type(np.array(df_barb_data.u_max))
#    print '[ma.getdata(x) for x in (np.array(df_barb_data.u_max))]', [ma.getdata(x) for x in (np.array(df_barb_data.u_max))]
#    print 'type(ma.getdata(x) for x in (np.array(df_barb_data.u_max)))', type([ma.getdata(x) for x in (np.array(df_barb_data.u_max))]) 
#    
#    print 'u_max', np.array([ma.getdata(x) for x in (np.array(df_barb_data.u_max))])
#    print 'u_max type', type(np.array([ma.getdata(x) for x in (np.array(df_barb_data.u_max))]))
#    print 'u_max shape', np.array([ma.getdata(x) for x in (np.array(df_barb_data.u_max))]).shape
#    
#    print 'v_max', np.array([ma.getdata(x) for x in (np.array(df_barb_data.v_max))])
#    print 'v_max type', type(np.array([ma.getdata(x) for x in (np.array(df_barb_data.v_max))]))
#    print 'v_max shape', np.array([ma.getdata(x) for x in (np.array(df_barb_data.v_max))]).shape
#    
#    print 'lon', np.array(df_barb_data.lon)
#    print 'lon type', type(np.array(df_barb_data.lon))
#    print 'lon shape', np.array(df_barb_data.lon).shape
#    
#    print 'lat', np.array(df_barb_data.lat)
#    print 'u_max type', type(np.array(df_barb_data.lat))
#    print 'u_max shape', np.array(df_barb_data.lat).shape
    
#    plt.scatter(np.array(df_barb_data.lon), np.array(df_barb_data.lat), s=10, transform=crs.PlateCarree(), zorder=6)
#    plt.barbs(np.array(df_barb_data.lon[::10]), np.array(df_barb_data.lat[::10]), np.array([ma.getdata(x) for x in (np.array(df_barb_data.u_max[::10]))]), np.array([ma.getdata(x) for x in (np.array(df_barb_data.v_max[::10]))]), transform=crs.PlateCarree(), length=5, zorder=6, color=np.array(df_barb_data.color))
    plt.barbs(to_np(lons[::20,::20]), to_np(lats[::20,::20]),  to_np(u_max[::20, ::20]), to_np(v_max[::20, ::20]), transform=crs.PlateCarree(), length=5, zorder=5)
    #plt.barbs(to_np(lons[::25,::25]), to_np(lats[::25,::25]),  to_np(u_level[::25, ::25]), to_np(v_level[::25, ::25]), transform=crs.PlateCarree(), length=5, zorder=6)
    #plt.barbs(to_np(lons[::10,::10]), to_np(lats[::10,::10]),  to_np(u_level[::10, ::10]), to_np(v_level[::10, ::10]), transform=crs.PlateCarree(), length=5)

    plt.title("%shPa Winds: %sZ" % (str(level), chosen_file_name[11:24]))

    #COR_xy = ll_to_xy_proj(ncfile, -31.298, -64.212, map_proj=1)
    #VMRS_xy = ll_to_xy_proj(ncfile, -29.906, -63.726, map_proj=1)
    
    #print 'COR_xy', COR_xy
    #print 'VMRS_xy', VMRS_xy
    
    #print lats[-140000,50000]
    #print lons[-140000,50000]
    
    # Get the basemap object
    bm = get_basemap(terrain)

    # Convert the lat/lon points in to x/y points in the projection space
    x, y = bm(to_np(lons), to_np(lats))
    
    #print 'x shape', x.shape
    #print 'y shape', y.shape

    #ax.scatter(-140000.,50000., s=70, color='blue', zorder=5)
    #ax.scatter(-60.,-30., s=70, color='yellow', zorder=5)
    #ax.scatter(-50.,-30., s=70, color='red', zorder=5)

    ax.scatter(-64.212, -31.298, s=30, color='gray', transform=crs.PlateCarree(), zorder=6)
    ax.scatter(-63.726, -29.906, s=30, color='gray', transform=crs.PlateCarree(), zorder=6)

    ## Set the map bounds
    #ax.set_xlim([-66,-60])
    #ax.set_ylim([-34,-28])

    #print 'saving'

    plt.savefig('/home/disk/meso-home/crs326/Documents/Research/PNNL_WRF/WRF_MAPS/222LLJ_smooth_PNNL_WRF_%sZ_Wind_%s_full-wind_grt12ms.png' %(chosen_file_name[11:24], str(level)), dpi=200)

    #print 'saved'
    print '---total time: %s seconds ---' %(time.time() - start_time)

