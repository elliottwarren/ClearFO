"""
Script to plot the MURK aerosol and RH for each ceil location
(should try to generalise it so it can plot the rain/fog/cloud later on)

Created by Elliott 10/10/2016
"""

import matplotlib.pyplot as plt
from matplotlib.dates import date2num
from matplotlib.dates import DateFormatter

from cycler import cycler
import numpy as np
import scipy.io as sp
import math
import datetime as dt

import ceilUtils as ceil
import cristinaCeilUtils as cutils
from ellUtils import add_at
import ellUtils as eu

import cartopy.crs as ccrs
import iris
import iris.fileformats.pp as pp
from netCDF4 import Dataset

def main():


    # -------------------------------------
    # Setup
    # -------------------------------------

    maindir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/'
    datadir = maindir + 'data/'


    # days to loop between
    dayStart = 125
    dayEnd = 127

    # year date for ceil data
    year = 2016

    # which modelled data to read in and directories
    model_type = 'UKV'
    savedir = maindir + 'figures/' + model_type + '/'

    # model resolution
    if model_type == 'UKV':
        res = '1p5km'
    elif model_type == '333m':
        res = '0p3'


    # -------------------------------------
    # Read
    # -------------------------------------

    # read in the ceil locations
    loc_filename = 'CeilsCSV.csv'
    loc_fname = datadir + loc_filename

    ceil_rawmeta = eu.csv_read(loc_fname)

    # convert into dictionary - [lon, lat]
    # skip the headers
    ceil_metadata = {}
    for i in ceil_rawmeta[1:]:
        ceil_metadata[i[1]] = [float(i[3]), float(i[4])]

    for day in dayStart:

        # temp filename for the day
        mod_filename = 'extract_prodm_op_ukv_20160504_21_full.nc'

        # concatenate the names
        mod_fname = datadir + mod_filename

        # Read in the modelled data for London
        mod_data = eu.netCDF_read(mod_fname)

        # rotate lats and lons to find location in rotated space

        # create rotated pole from metadata

        # rot_pole1 = temp_full.coord('grid_latitude').coord_system.as_cartopy_crs()
        # datafile = Dataset(mod_fname, 'r')
        # rot_pole1 = datafile.create_grid_obj

        for site, loc in ceil_metadata.iteritems():

            # give full names to play safe
            lat = loc[0]
            lon = loc[1]

            # ToDo - projection may well change if not UKV
            if model_type == 'UKV':
                # create the rotated pole system from the UKV metadata
                rot_pole1 = iris.coord_systems.RotatedGeogCS(37.5, 177.5, ellipsoid=iris.coord_systems.GeogCS(6371229.0)).as_cartopy_crs()

            ll = ccrs.Geodetic()
            target_xy1 = rot_pole1.transform_point(lat, lon, ll)

            if (res == '2p2km') or (res == '0p5km') or (res == '1p5km'):
                glon = target_xy1[0] + 360
            else:
                glon = target_xy1[0]

            glat = target_xy1[1]

            # define half grid spacing
            if res == '2p2km':
                delta = 0.01
            elif res == '0p5km':
                delta = 0.0022
            elif res == '0p2km':
                delta = 0.0009
            elif res == '0p2km':  # guess using 0p2km and that 1p5 is just OVER 1/3 of 0p5
                delta = 0.00044
            elif res == '1p5km':
                delta = 0.5 * 0.0135

            # get idx location in rotated space
            # define rotated lat and lon
            # don't know why some have 360 and others do not
            # for some reason unknown to science, 100 m does not want + 360...
            if (res == '2p2km') or (res == '0p5km') or (res == '1p5km'):
                glon = target_xy1[0] + 360
            else:
                glon = target_xy1[0]

            glat = target_xy1[1]

            # idx match with the lat and lon from the model
            mod_glon, idx_lon, diff_lon = eu.nearest(mod_data['longitude'], glon)
            mod_glat, idx_lat, diff_lat = eu.nearest(mod_data['latitude'], glat)


            # extract the variables for that location
            mod_aer = mod_data['aerosol_for_visibility'][:, :, idx_lat, idx_lon]
            mod_q = mod_data['specific_humidity'][:, :, idx_lat, idx_lon]
            mod_p = mod_data['air_pressure'][:, :, idx_lat, idx_lon]
            mod_t = mod_data['air_temperature'][:, :, idx_lat, idx_lon]

            # convert temperature to degree C
            mod_t_celsius = mod_t-273.15

            # calculate relative humidity
            # Thermal Physics of the Atmosphere - Maarten's book.
            # -----------
            # saturated vapour pressure (hPa, then Pa) - Clausius-Clapeyron eq 5.18, pp. 108
            e_s_hpa = 6.112*(np.exp((17.67*mod_t_celsius)/(mod_t_celsius+243.5)))
            e_s = e_s_hpa*100

            # mass mixing ratio of water vapour pp. 100
            r_v = mod_q / (1-mod_q)
            # mass mixing ratio of water vapour at saturation eq 5.22, pp. 100
            r_vs = 0.622 * (e_s / mod_p)

            # relative humidity (variant of eq 5.24, pp 101)
            rh = r_v/ r_vs


            # ------------------------------------
            # Plot
            # ------------------------------------

            # fontsize
            fs = 12

            # RH
            fig, ax = plt.subplots(2, 1, figsize=(8, 5))
            im_rh = ax[0].pcolormesh(mod_data['time'], mod_data['level_height'], np.transpose(100.0 * rh), vmin=0, vmax=100)
            im_aer = ax[1].pcolormesh(mod_data['time'], mod_data['level_height'], np.transpose(mod_aer), vmin=0, vmax=150)

            # adjust subplots to fit colour bars
            plt.subplots_adjust(right=0.8)

            # # set up colour bars
            # pos = ax[0].get_position()._get_extents()
            # rh_cax = fig.add_axes([pos[2]+0.01, pos[1], 0.05, pos[3]-0.54])
            # fig.colorbar(im_rh, cax=rh_cax)
            #
            # pos = ax[1].get_position()._get_extents()
            # aer_cax = fig.add_axes([pos[2]+0.01, pos[1], 0.05, pos[3]-0.1])
            # fig.colorbar(im_aer, cax=aer_cax)

            # prettify
            ax0 = eu.fig_majorAxis(fig)
            ax0.set_xlabel('Time [HH:MM]', fontsize=fs)
            ax0.set_ylabel('Height [m]', fontsize=fs)

            ax[0].set_ylim([0, 8000])
            ax[0].axis('tight')
            ax[1].set_ylim([0, 8000])
            ax[1].set_yticklabels('')
            ax[1].axis('tight')
            plt.tight_layout()

            # adjust label position, needs to be done after tight_layout() because it doesn't work well with
            # separate overlaying axis
            ax0.yaxis.set_label_coords(-0.08, 0.5)
            ax0.xaxis.set_label_coords(0.4, -0.15)

            # alphabet label subplots
            add_at(ax[0], 'a)', loc=2)
            add_at(ax[1], 'b)', loc=2)

            # colourbar
            # cax, kw = mpl.colorbar.make_axes([ax for ax in ax.flat], aspect=12)
            # plt.colorbar(im, cax=cax, **kw)

            # save
            plt.savefig(savedir + 'dailyPlots/' +
                        model_type + '_' + site + '_RH_MMR_' + str(day) + 'raw.png')  # filename


            plt.close(fig)

    return



if __name__ == '__main__':
    main()