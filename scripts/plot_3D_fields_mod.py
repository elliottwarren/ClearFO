"""
Plots cross sections and 3D vector fields of model variables to help diagnose stuff

Created by Elliott 06/06/17
"""

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
from matplotlib.dates import DateFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable

import numpy as np
import datetime as dt

from copy import deepcopy
import colorsys

import ellUtils as eu
import ceilUtils as ceil
from forward_operator import FOUtils as FO
from forward_operator import FOconstants as FOcon

def dateList_to_datetime(dayList):

    """ Convert list of string dates into datetimes """

    datetimeDays = []

    for d in dayList:

        datetimeDays += [dt.datetime(int(d[0:4]), int(d[4:6]), int(d[6:8]))]

    return datetimeDays

def plot_time_height_cross(site, mod_all_data, CSAT3, model_type, Z, day, idx_lon, idx_lat, savedir):

    """ Plots cross section (time x height) of u, v and w for a site"""

    # fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(7, 8))
    fig, (ax3, ax4) = plt.subplots(2, 1, figsize=(7, 4))
    # # height[11] = 555 m; .  '4D variable'[time, height, lat, lon]
    # mesh1 = ax1.pcolormesh(mod_all_data['time'], mod_all_data['level_height'],
    #                      np.transpose(mod_all_data['x_wind'][:, :, idx_lat, idx_lon]),
    #                      vmin=-4.0, vmax=4.0, cmap=cm.get_cmap('jet'))
    #
    # mesh2 = ax2.pcolormesh(mod_all_data['time'], mod_all_data['level_height'],
    #                      np.transpose(mod_all_data['y_wind'][:, :, idx_lat, idx_lon]),
    #                      vmin=-4.0, vmax=4.0, cmap=cm.get_cmap('jet'))

    mesh3 = ax3.pcolormesh(mod_all_data['time'], mod_all_data['level_height'],
                         np.transpose(mod_all_data['upward_air_velocity'][:, :, idx_lat, idx_lon]),
                         vmin=-0.05, vmax=0.05, cmap=cm.get_cmap('jet'))

    line = ax4.plot(CSAT3['time'], CSAT3['w'], linestyle='-', marker='o', markersize=4)


    # prettify
    # for mesh, ax in zip((mesh1, mesh2, mesh3), (ax1, ax2, ax3)):
    #
    #     ax.xaxis.set_major_formatter(DateFormatter('%H:%M'))
    #     ax.yaxis.label.set_size(10)
    #     ax.xaxis.label.set_size(10)
    #     ax.set_xlim([day, day + dt.timedelta(days=1)])
    #     ax.set_ylim([0, 1500.0])
    #
    #     divider = make_axes_locatable(ax)
    #     cax = divider.append_axes("right", size="5%", pad=0.05)
    #     plt.colorbar(mesh, cax=cax)

    for ax in (ax3, ax4):

        ax.xaxis.set_major_formatter(DateFormatter('%H:%M'))
        ax.yaxis.label.set_size(10)
        ax.xaxis.label.set_size(10)
        ax.set_xlim([day, day + dt.timedelta(days=1)])

    ax3.set_ylim([0, 1500.0])

    divider = make_axes_locatable(ax3)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(mesh3, cax=cax)

    divider = make_axes_locatable(ax4)
    cax2 = divider.append_axes("right", size="5%", pad=0.05)
    #cb = plt.colorbar(mesh3, cax=cax2)
    #cb.remove()
    cax2.axis('off')

    # plt.suptitle(site + ': ' + str(day) + ' units: m s-1')

    # eu.add_at(ax1, r'$u$', loc=2)
    # eu.add_at(ax2, r'$v$', loc=2)
    # eu.add_at(ax3, r'$w$', loc=2)

    ax3.set_ylabel(r'$Height \/\/[m]$')
    ax4.set_ylabel(r'$w_{o} \/\/[m \/s^{-1}]$')
    ax4.set_xlabel(r'$Time \/\/[HH:MM]$')
    ax4.set_ylim([-0.4, 0.4])

    # plt.tight_layout(h_pad=0.1)

    plt.savefig(savedir + model_type + '_' + site + '_wind_CSAT3_time_height_cross_' + day.strftime('%Y%m%d') +
                '_' + str(Z) + 'Z.png')  # filename

    plt.close(fig)

    return

def plot_2D_wind(mod_all_data, ceil_metadata, model_type, res, savedir, Z, t, h_levels):

    for h in np.arange(h_levels):
        # plot 2D fields at fixed time and height of area over London
        fig, ax = plt.subplots(1, 1, figsize=(6, 6))

        # ----------- NOTE! ----------------------
        # current issue is that x_wind has shape (14,12) whereas longitude is (15) not (14)...
        # would need to trim off of y wind and longitude... but which side to trim; do I trim the start or end?

        # make ure aerosol is transposed, as it needs shape (latitude, longitude) because it is (rows, columns), (y, x)
        # it WILL plot if rows and columns are the wrong way round AND if the wrong shape!
        # height[11] = 555 m; time[15] = 15:00.  '4D variable'[time, height, lat, lon]
        quiv = ax.quiver(mod_all_data['longitude'][:-1], mod_all_data['latitude'],
                         mod_all_data['x_wind'][t, h, :, :], mod_all_data['y_wind'][t, h, :, :-1],
                         np.transpose(mod_all_data['aerosol_for_visibility'][t, h, :, :-1]),
                         units='width', clim=[10.0, 80.0])

        plt.colorbar(quiv)

        qk = plt.quiverkey(quiv, 0.9, 0.9, 2, r'$2 \frac{m}{s}$', labelpos='E',
                           coordinates='figure')

        # plot each ceilometer location
        for site, loc in ceil_metadata.iteritems():
            idx_lon, idx_lat, glon, glat = FO.get_site_loc_idx_in_mod(mod_all_data, loc, model_type, res)
            plt.scatter(glon, glat, color='black')
            plt.annotate(site, (glon, glat))

            # prettify
            # ax.xaxis.set_major_formatter(DateFormatter('%H:%M'))
            # ax.yaxis.label.set_size(10)
            # ax.xaxis.label.set_size(10)
            # ax.set_xlim([day, day + dt.timedelta(days=1)])
            # ax.set_ylim([0, 1500.0])

            # divider = make_axes_locatable(ax)
            # cax = divider.append_axes("right", size="5%", pad=0.05)
            # plt.colorbar(mesh, cax=cax)

        plt.suptitle(
            site + ': ' + str(mod_all_data['time'][t]) + ' height: ' + str(mod_all_data['level_height'][h]) + ' m')

        plt.savefig(savedir + '2DwindFields/' + model_type + '_2Dfield_' +
                    mod_all_data['time'][t].strftime('%Y%m%d_%Hhr') + '_' + str(Z) + 'Z' + '_' +
                    str(mod_all_data['level_height'][h]) + 'm.png')  # filename


    return

def plot_cloud_time_height_cross(mod_all_data, cloud_obs_i, idx_lat_main, idx_lon_main, main_site, savedir, model_type, Z, day):

    """ single subplot of cloud fraction, time and height cross section over the ceilometer site"""

    fig, ax = plt.subplots(1, 1, figsize=(7, 3))

    # height[11] = 555 m; .  '4D variable'[time, height, lat, lon]
    mesh = ax.pcolormesh(mod_all_data['time'], mod_all_data['level_height'],
                         np.transpose(mod_all_data['cloud_volume_fraction_in_atmosphere_layer'][:, :, idx_lat_main,
                                      idx_lon_main]), vmin=0.0, vmax=1.0, cmap=cm.get_cmap('Blues'))

    # overylay plot of CBH
    line = ax.plot_date(cloud_obs_i['time'],cloud_obs_i['CLD_Height_L1'], marker='x', color='red')

    # # top ceilometer height
    # ax.plot_date([cloud_obs_i['time'][0],cloud_obs_i['time'][-1]],
    #              [cloud_obs_i['height']+7500.0,[cloud_obs_i['height']+7500.0]],
    #               fmt='-', color='black')

    # prettify
    ax.xaxis.set_major_formatter(DateFormatter('%H:%M'))
    ax.yaxis.label.set_size(10)
    ax.xaxis.label.set_size(10)
    ax.set_xlim([day, day + dt.timedelta(days=1)])
    ax.set_ylim([0, 12000.0])

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(mesh, cax=cax)

    ax.set_ylabel(r'$Height \/\/[m]$')
    ax.set_xlabel(r'$Time \/\/[HH:MM]$')

    plt.tight_layout()

    # eu.add_at(ax, r'$cloud \/fract$', loc=1)

    plt.savefig(savedir + model_type + '_' + main_site + '_cloud_time_height_cross_' + day.strftime('%Y%m%d') +
                '_' + str(Z) + 'Z.png')  # filename

    plt.close(fig)


    return


def main():

    # ==============================================================================
    # Setup
    # ==============================================================================

    # which modelled data to read in
    model_type = 'UKV'
    res = FOcon.model_resolution[model_type]

    # directories
    maindir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/'
    datadir = maindir + 'data/'
    savedir = maindir + 'figures/' + model_type + '/CrossSecs_3Dfields/'
    # savedir = maindir + 'figures/' + model_type + '/highPmCase/'

    # data
    ceilDatadir = datadir + 'L0/'
    modDatadir = datadir + model_type + '/'
    rhDatadir = datadir + 'L1/'
    aerDatadir = datadir + 'LAQN/'

    # statistics to run
    pm10_stats = False
    rh_stats = True


    site_bsc = FOcon.site_bsc
    site_bsc_i = {'CL31-C_MR': 4.5}

    site_bsc_colours = FOcon.site_bsc_colours

    # daystrList = ['20150414', '20150415', '20150421', '20150611', '20160504', '20160823', '20160911', '20161125',
    #               '20161129', '20161130', '20161204']

    # days_iterate = dateList_to_datetime(daystrList)

    corr_max_height = 2000.0

    # forecast data start time
    Z='21'

    # main site for the u, v and w plot.
    main_site = 'MR'
    ceil_id = 'CL31-C'
    ceil_id_full = ceil_id + '_' + main_site

    days_iterate = [dt.datetime(2016,01,19)]# new PM10 case study day
    # day = [dt.datetime(2016, 05, 04)] # one of my old case study days

    # plot ?
    windfield = True # 2D wind field
    windcross = False # u, v, w at site in height and time
    cloudcross = False # 2D cloud field

    # ==============================================================================
    # Read data
    # ==============================================================================

    # Read Ceilometer metadata

    # ceilometer list to use
    ceilsitefile = 'CeilsCSVfull.csv'
    ceil_metadata = FO.read_ceil_metadata(datadir, ceilsitefile)

    # get metadata just for single site
    ceil_data_i = {main_site: ceil_metadata[main_site]}
    loc = ceil_metadata[main_site]

    for day in days_iterate:

        print 'day = ' + day.strftime('%Y-%m-%d')

        # Read all UKV forecast and automatically run the FO
        mod_all_data = FO.read_all_mod_data(modDatadir, day, Z)

        # get time idx for day in mod_all_data and trim it
        range_time = FO.get_time_idx_forecast(mod_all_data, day)
        mod_all_data['x_wind'] = mod_all_data['x_wind'][range_time, :, :, :]
        mod_all_data['y_wind'] = mod_all_data['y_wind'][range_time, :, :, :]
        mod_all_data['upward_air_velocity'] = mod_all_data['upward_air_velocity'][range_time, :, :, :]
        mod_all_data['aerosol_for_visibility'] = mod_all_data['aerosol_for_visibility'][range_time, :, :, :]
        mod_all_data['cloud_volume_fraction_in_atmosphere_layer'] = \
            mod_all_data['cloud_volume_fraction_in_atmosphere_layer'][range_time, :, :, :]
        mod_all_data['time'] = np.array(mod_all_data['time'])[range_time]

        mod_site_data = {main_site: mod_all_data}

        # read in cloud base heights from obs
        cloud_obs = FO.read_ceil_CLD_obs(day, site_bsc_i, ceilDatadir, mod_site_data, dayExtract=True)

        # read in CSAT w obs
        datapath = datadir + 'L1/CSAT3_KSSW_2016019_30min.nc'
        CSAT3 = eu.netCDF_read(datapath, vars='')


        # get idx location for instrument
        idx_lon_main, idx_lat_main, _, _ = FO.get_site_loc_idx_in_mod(mod_all_data, loc, model_type, res)

        # plots time x height cross section of vertical wind for an instrument
        # plot_time_height_cross(site, mod_all_data, model_type, Z, day, idx_lon_main, idx_lat_main, savedir)

        if windcross == True:

            plot_time_height_cross(main_site, mod_all_data, CSAT3, model_type, Z, day, idx_lon_main, idx_lat_main, savedir)

        if windfield == True:

            # plot 2D wind fields
            t = 19 # time slice
            h_levels = 12 # number of height levels (will a plot for every level between 0 and h_levels)
            plot_2D_wind(mod_all_data, ceil_metadata, model_type, res, savedir, Z, t, h_levels)

        if cloudcross == True:

            # cloud fraction time height cross section over main ceilometer site
            plot_cloud_time_height_cross(mod_all_data, cloud_obs['CL31-C_MR'], idx_lat_main, idx_lon_main, main_site, savedir,
                                         model_type, Z, day)

    return

if __name__ == '__main__':
    main()


