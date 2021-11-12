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
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, interplevel, ll_to_xy_proj, get_basemap, xy_to_ll, CoordPair, vertcross)

# Open the NetCDF file
path_wrf = '/home/disk/monsoon/relampago/raw/wrf/'

#chosen_date_time_list = ['2018-11-02_07_00_00']

chosen_date_time_list = ['2018-11-02_03_00_00', '2018-11-02_04_00_00', '2018-11-02_05_00_00', '2018-11-02_06_00_00', '2018-11-02_07_00_00', '2018-11-02_08_00_00', '2018-11-02_09_00_00', '2018-11-02_10_00_00', '2018-11-02_11_00_00', '2018-11-02_12_00_00', '2018-11-02_13_00_00', '2018-11-02_14_00_00', '2018-11-02_15_00_00', '2018-11-02_16_00_00', '2018-11-02_17_00_00', '2018-11-02_18_00_00', '2018-11-02_19_00_00', '2018-11-02_20_00_00', '2018-11-02_21_00_00', '2018-11-02_22_00_00','2018-11-02_23_00_00', '2018-11-03_00_00_00']


######################################

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

    # Create a figure
    fig = plt.figure(figsize=(12,6))

    # Define the cross section start and end points
    cross_start = CoordPair(lat=25, lon=-75)
    cross_end = CoordPair(lat=31, lon=-65)
    
    print cross_start
    
    # Get the WRF variables
    ht = getvar(ncfile, 'z', units='m') # in meters
    v = getvar(ncfile, 'va', units='m s-1') # in m/s
    terrain = getvar(ncfile, 'HGT')
    
    # Get the lat/lon points
    lats, lons = latlon_coords(v)

    # Get the cartopy projection object
    cart_proj = get_cartopy(v)

    # Compute the vertical cross-section interpolation.  Also, include the
    # lat/lon points along the cross-section in the metadata by setting latlon
    # to True.
    v_cross = vertcross(v, ht, wrfin=ncfile, projection=cart_proj, start_point=cross_start, end_point=cross_end, latlon=True, meta=True)

    # To remove the slight gap between the dbz contours and terrain due to the
    # contouring of gridded data, a new vertical grid spacing, and model grid
    # staggering, fill in the lower grid cells with the first non-missing value
    # for each column.

    # Make a copy of the z cross data. Let's use regular numpy arrays for this.
    v_cross_filled = np.ma.copy(to_np(v_cross))

    # For each cross section column, find the first index with non-missing
    # values and copy these to the missing elements below.
    for i in range(v_cross_filled.shape[-1]):
        column_vals = v_cross_filled[:,i]
        # Let's find the lowest index that isn't filled. The nonzero function
        # finds all unmasked values greater than 0. Since 0 is a valid value
        # for dBZ, let's change that threshold to be -200 dBZ instead.
        first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
        v_cross_filled[0:first_idx, i] = v_cross_filled[first_idx, i]

    # Get the terrain heights along the cross section line
    ter_line = interpline(ter, wrfin=ncfile, start_point=cross_start,
                          end_point=cross_end)

    # Create the figure
    fig = pyplot.figure(figsize=(8,6))
    ax_cross = pyplot.axes()

    dbz_levels = np.arange(5., 50., 5.)

    # Make the cross section plot for v
    v_levels = np.arange(5.,75.,5.)
    xs = np.arange(0, v_cross.shape[-1], 1)
    ys = to_np(v_cross.coords["vertical"])
    v_contours = v_cross.contourf(xs, ys, to_np(v_cross_filled), levels=v_levels, cmap='YlOrRd', extend='max')

    # Add the color bar
    cb_v = fig.colorbar(v_contours, ax=ax_cross)
    cb_v.ax.tick_params(labelsize=8)

    # Fill in the mountain area
    ht_fill = ax_cross.fill_between(xs, 0, to_np(ter_line),
                                    facecolor="saddlebrown")

    # Set the x-ticks to use latitude and longitude labels
    coord_pairs = to_np(v_cross.coords["xy_loc"])
    x_ticks = np.arange(coord_pairs.shape[0])
    x_labels = [pair.latlon_str() for pair in to_np(coord_pairs)]

    # Set the desired number of x ticks below
    num_ticks = 5
    thin = int((len(x_ticks) / num_ticks) + .5)
    ax_cross.set_xticks(x_ticks[::thin])
    ax_cross.set_xticklabels(x_labels[::thin], rotation=45, fontsize=8)

    # Set the x-axis and  y-axis labels
    ax_cross.set_xlabel("Latitude, Longitude", fontsize=12)
    ax_cross.set_ylabel("Height (m)", fontsize=12)

    # Add a title
    ax_cross.set_title("Cross-Section of V-Wind (m/s)", {"fontsize" : 14})

    plt.savefig('/home/disk/meso-home/crs326/Documents/Research/PNNL_WRF/WRF_Cross_section/test_%s.png' %(chosen_file_name[11:24]), dpi=200)
    print 'saved'
