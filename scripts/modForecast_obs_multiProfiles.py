"""
Plot a profile each hour of the forward modelled backscatter and ceil obs. As well as the whole period of RH and PM10.
Improved upon the original version to take netCDF of UKV forecast.

Created by Elliott Tues 18 Oct 2016

"""


import matplotlib.pyplot as plt
from matplotlib.dates import date2num
from matplotlib.dates import DateFormatter

import numpy as np
import datetime as dt
import cartopy.crs as ccrs
import iris

import ceilUtils as ceil
import ellUtils as eu
import FOUtils as FO

# Read

def read_ceil_metadata(datadir):

    """
    Read in ceil metadata (lon, lat) into a dictionary
    :param datadir:
    :return:
    """

    # read in the ceil locations
    loc_filename = 'CeilsCSV.csv'
    loc_fname = datadir + loc_filename

    ceil_rawmeta = eu.csv_read(loc_fname)

    # convert into dictionary - [lon, lat]
    # skip the headers
    ceil_metadata = {}
    for i in ceil_rawmeta[1:]:
        ceil_metadata[i[1]] = [float(i[3]), float(i[4])]

    return ceil_metadata

def mod_site_extract_calc(ceil_metadata, mod_all_data, model_type, res, day):

    # define mod_data array
    mod_data = {}

    for site, loc in ceil_metadata.iteritems():

        # define dictionary for the site
        mod_data[site] = {}

        # give full names to play safe
        lat = loc[0]
        lon = loc[1]

        # ToDo - projection may well change if not UKV
        if model_type == 'UKV':
            # create the rotated pole system from the UKV metadata
            rot_pole1 = iris.coord_systems.RotatedGeogCS(37.5, 177.5, ellipsoid=iris.coord_systems.GeogCS(
                6371229.0)).as_cartopy_crs()

        ll = ccrs.Geodetic()
        target_xy1 = rot_pole1.transform_point(lat, lon, ll)

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
        mod_glon, idx_lon, diff_lon = eu.nearest(mod_all_data['longitude'], glon)
        mod_glat, idx_lat, diff_lat = eu.nearest(mod_all_data['latitude'], glat)

        # print ""
        # print site + ' - idx_lon = ' + str(idx_lon) + '; idx_lat = ' + str(idx_lat)

        # only extract data for the main day
        _, startIdx_time, _ = eu.nearest(mod_all_data['time'], (day + dt.timedelta(hours=24)))
        _, endIdx_time, _ = eu.nearest(mod_all_data['time'], (day + dt.timedelta(hours=48)))
        range_time = np.arange(startIdx_time, endIdx_time + 1)

        # extract the variables for that location
        mod_aer = mod_all_data['aerosol_for_visibility'][range_time, :, idx_lat, idx_lon]
        mod_q = mod_all_data['specific_humidity'][range_time, :, idx_lat, idx_lon]
        mod_p = mod_all_data['air_pressure'][range_time, :, idx_lat, idx_lon]
        mod_t = mod_all_data['air_temperature'][range_time, :, idx_lat, idx_lon]
        mod_h = mod_all_data['level_height']
        mod_time = [mod_all_data['time'][i] for i in range_time]

        # convert temperature to degree C
        mod_t_celsius = mod_t - 273.15

        # calculate relative humidity
        # Thermal Physics of the Atmosphere - Maarten's book.
        # -----------
        # saturated vapour pressure (hPa, then Pa) - Clausius-Clapeyron eq 5.18, pp. 108
        e_s_hpa = 6.112 * (np.exp((17.67 * mod_t_celsius) / (mod_t_celsius + 243.5)))
        e_s = e_s_hpa * 100

        # mass mixing ratio of water vapour pp. 100
        r_v = mod_q / (1 - mod_q)
        # mass mixing ratio of water vapour at saturation eq 5.22, pp. 100
        r_vs = 0.622 * (e_s / mod_p)

        # relative humidity (variant of eq 5.24, pp 101)
        mod_rh = r_v / r_vs

        # prcoess forward modelled backscatter for each site
        mod_data[site]['backscatter'] = FO.forward_operator(mod_aer, mod_rh, mod_h)

        # store MURK aerosol, RH and heights in mod_data dictionary
        mod_data[site]['RH'] = mod_rh
        mod_data[site]['aerosol_for_visibility'] = mod_aer
        mod_data[site]['level_height'] = mod_h
        mod_data[site]['time'] = mod_time

    return mod_data

def read_all_rh_obs(day, site_rh, rhDatadir, mod_data):

    """
    Read in day and following day's data, for all rh obs.

    :param day: day string
    :param site_rh: all rh sites
    :param rhDatadir: data directory for rh
    :return: rh obs: dictionary
    """

    # define array
    rh_obs = {}

    # get date string for obs of the main and following days
    doyStr = (day + dt.timedelta(hours=24)).strftime('%Y%j')
    doyStr2 = (day + dt.timedelta(hours=48)).strftime('%Y%j')

    for site, height in site_rh.iteritems():

        rh_obs[site] = {}

        rh_fnames = [rhDatadir + site + '_' + doyStr + '_1min.nc',
                     rhDatadir + site + '_' + doyStr2 + '_1min.nc']

        # read in all data
        data_obs = eu.netCDF_read(rh_fnames, vars=['RH', 'time'])
        data_obs['height'] = height

        # find nearest time in rh time
        # pull out ALL the nearest time idxs and differences
        t_idx = np.array([eu.nearest(data_obs['time'], t)[1] for t in mod_data['NK']['time']])
        t_diff = np.array([eu.nearest(data_obs['time'], t)[2] for t in mod_data['NK']['time']])

        # extract hours
        rh_obs[site]['RH'] = data_obs['RH'][t_idx]
        rh_obs[site]['height'] = data_obs['height']
        rh_obs[site]['time'] = [data_obs['time'][i] for i in t_idx]

        # overwrite t_idx locations where t_diff is too high with nans
        # only keep t_idx values where the difference is below 5 minutes
        bad = np.array([abs(i.days * 86400 + i.seconds) > 10 * 60 for i in t_diff])
        rh_obs[site]['RH'][bad] = np.nan

        # change flags to nans
        rh_obs[site]['RH'][np.where(rh_obs[site]['RH'] < 0)] = np.nan

    return rh_obs

def read_all_pm10_obs(dayStart, dayEnd, site_aer, aerDatadir, mod_data):

    """
    Read in the LAQN data.
    It comes in a single file and can contain mempty strings for missing data, but the code
    will replace them with nans. Extracts pm10 and time, but inflexable if data columns would change.
    :param dayStart:
    :param dayEnd:
    :param site_aer:
    :param aerDatadir:
    :return:
    """

    # define array
    pm10_obs = {}

    # get date string for obs of the main and following days
    dateStr = (dayStart + dt.timedelta(hours=24)).strftime('%Y%m%d') + '-' + \
              (dayEnd + dt.timedelta(hours=48)).strftime('%Y%m%d')

    for site, height in site_aer.iteritems():

        pm10_obs[site] = {}

        aer_fname = aerDatadir + site + '_' + dateStr + '.csv'

        raw_aer = np.genfromtxt(aer_fname, delimiter=',', skip_header=1, dtype="|S20")

        # convert missing sections ('' in LAQN data) to NaNs
        missing_idx = np.where(raw_aer[:, 3] == '')
        raw_aer[missing_idx, 3] = np.nan

        # extract obs and time together as a dictionary entry for the site.
        data_obs = \
            {'pm_10': np.array(raw_aer[:, 3], dtype=float),
             'time': np.array(
                 [dt.datetime.strptime(raw_aer[i, 2], '%d/%m/%Y %H:%M')
                  for i in np.arange(len(raw_aer))]),
             'height': height}


        # find nearest time in rh time
        # pull out ALL the nearest time idxs and differences
        t_idx = np.array([eu.nearest(data_obs['time'], t)[1] for t in mod_data['NK']['time']])
        t_diff = np.array([eu.nearest(data_obs['time'], t)[2] for t in mod_data['NK']['time']])

        # extract ALL nearest hours data, regardless of time difference
        pm10_obs[site]['pm_10'] = data_obs['pm_10'][t_idx]
        pm10_obs[site]['height'] = data_obs['height']
        pm10_obs[site]['time'] = [data_obs['time'][i] for i in t_idx]

        # overwrite t_idx locations where t_diff is too high with nans
        # only keep t_idx values where the difference is below 5 minutes
        bad = np.array([abs(i.days * 86400 + i.seconds) > 10 * 60 for i in t_diff])
        pm10_obs[site]['pm_10'][bad] = np.nan

    return pm10_obs

def read_ceil_obs(day, site_bsc, ceilDatadir, mod_data):

    """
    Read in ceilometer backscatter, time, height and SNR data and strip the hours out of it.

    :param day:
    :param ceilDatadir:
    :return: data_obs: dictionary
    """

    # contains all the sites time-upscaled data
    bsc_obs = {}

    for site in site_bsc.keys():

        # this sites time-upscaled data
        bsc_obs[site] = {}

        # date for the main day
        doyStr = (day + dt.timedelta(hours=24)).strftime('%Y%j')

        # get filename
        bsc_fname = ceilDatadir + site + '_' + doyStr + '_15sec.nc'

        # read backscatter data
        data_obs = ceil.netCDF_read_BSC(bsc_fname)

        # find nearest time in ceil time
        # pull out ALL the nearest time idxs and differences
        t_idx = np.array([eu.nearest(data_obs['time'], t)[1] for t in mod_data['NK']['time']])
        t_diff = np.array([eu.nearest(data_obs['time'], t)[2] for t in mod_data['NK']['time']])

        # extract data
        # for var, data in data_obs.iteritems():
        bsc_obs[site]['SNR'] = data_obs['SNR'][t_idx, :]
        bsc_obs[site]['backscatter'] = data_obs['backscatter'][t_idx, :]
        bsc_obs[site]['height'] = data_obs['height']
        bsc_obs[site]['time'] = [data_obs['time'][i] for i in t_idx]

        # overwrite t_idx locations where t_diff is too high with nans
        # only keep t_idx values where the difference is below 5 minutes
        bad = np.array([abs(i.days * 86400 + i.seconds) > 10 * 60 for i in t_diff])

        bsc_obs[site]['SNR'][bad, :] = np.nan
        bsc_obs[site]['backscatter'][bad, :] = np.nan

    return bsc_obs


# Plotting

def bsc_profile_plot(ax, mod_data, bsc_obs, t):

    """
    plot the modelled and observed backscatter profiles

    :param ax:
    :param mod_data:
    :param bsc_obs:
    :param t:
    :return: ax
    """


    # ToDo - if the legend gets messy. Remove label from obs and just use the one from mod. Then set the colour code
    # to do a set rotation. As dictionary can be plotted in any order, the mod site will need to 'get' the colour
    # of the obs at the same site, to know they match.


    # plot obs
    for site_obs, data_obs in bsc_obs.iteritems():

        ax.plot(data_obs['backscatter'][t, :], data_obs['height'], label=site_obs, ls='-')

    # plot model
    for site_mod, data_mod in mod_data.iteritems():

        ax.plot(data_mod['backscatter'][t, :], data_mod['level_height'], label=site_mod, linewidth=2, ls='--')


    # Prettify figure at the end
    ax.set_xlabel(r'$log_{10}(\beta) \/\/\mathrm{[m^{-1} \/sr^{-1}]}$')
    ax.set_ylabel('Height [m]')

    ax.set_ylim([0, 2000])
    ax.set_xlim([-7, -5])
    ax.xaxis.set_ticks(np.arange(-7, -4, 0.5))
    # ax.set_prop_cycle(cycler('color', ['b', 'c', 'm', 'y', 'k']))
    eu.add_at(ax, 'a)', loc=2)

    return ax

def rh_profile_plot(ax, site_rh, rh_obs, mod_data, t):

    """
    Plot the RH profile and obs

    :param ax:
    :param site_rh:
    :param rh_obs:
    :param mod_data:
    :param t:
    :return: ax
    """

    # plot obs rh
    for site, height in site_rh.iteritems():

        # currently no label, so it can be added manually later
        ax.scatter(rh_obs[site]['RH'][t], height, color='green', edgecolors='black')


    # plot model rh
    for site_mod, data_mod in mod_data.iteritems():

        ax.plot(data_mod['RH'][t, :] * 100, data_mod['level_height'], 'b-', linewidth=2)
        # ax.plot(data_mod['model'][t] * 100, data_mod['level_height'], 'b-', label='RH'+r'$_{100}$', linewidth=2)

    # plot reference line
    ax.plot([38.0, 38.0], [0, 10000], color='black', ls='--')

    # prettify
    ax.set_ylim([0, 2000])
    ax.set_xlim([20, 100])
    ax.set_xlabel(r'$\mathrm{RH\/\/[\%]}$')
    ax.set_yticklabels('')
    # ax.set_prop_cycle(cycler('color', ['b', 'c', 'm', 'y', 'k']))
    eu.add_at(ax, 'b)', loc=2)

    return ax

def aer_profile_plot(ax, site_aer, pm10_obs, mod_data, t):

    """
    Plot the pm10 obs from LAQN and the modelled murk aerosol

    :param ax3:
    :param site_aer:
    :param pm10_obs:
    :param mod_data:
    :param t:
    :return: ax
    """

    # plot obs pm10
    for site, height in site_aer.iteritems():

        ax.scatter(pm10_obs[site]['pm_10'][t], height, color='red', edgecolors='black')

    # plot model pm10
    for site_mod, data_mod in mod_data.iteritems():

        ax.plot(data_mod['aerosol_for_visibility'][t, :], data_mod['level_height'],
                color='black', linewidth=2, ls='--')

    # prettify
    ax.set_ylim([0, 2000])
    ax.set_xlim([0, 100])
    ax.set_yticklabels('')
    # ax.set_xlabel(r'$m$' + ' [' + r'$\mu$' + 'g kg' + r'$-1$' +']')
    ax.set_xlabel(r'$m \/\mathrm{[\mu g\/ kg^{-1}]}$')
    # ax.set_prop_cycle(cycler('color', ['b', 'c', 'm', 'y', 'k']))
    eu.add_at(ax, 'c)', loc=1)

    return ax

def rh_ts_plot(ax, rh_obs, mod_data, t):

    """
    Plot the rh time series for the period. Uses average rh model input calculated
    outside this function and the unchanged RH observations.

    :param ax4:
    :param site_rh:
    :param rh_obs:
    :param mod_data:
    :param t:
    :return: ax
    """


    # plot rh obs profile plots
    for site_obs, data_obs in rh_obs.iteritems():
        ax.plot_date(date2num(data_obs['time']),data_obs['RH'], label=site_obs, fmt='-')

    # plot model for a single height near the surface (5 m)
    # height index to plot
    h_idx = 0

    for site_mod, data_mod in mod_data.iteritems():
        # ax.plot_date(date2num(data_mod['time']), data_mod['RH'][t, h_idx] * 100,
        #              label=str(data_mod['level_height'][h_idx]) + ' m model', linewidth=2, fmt='-', ls='-')
        ax.plot_date(date2num(data_mod['time']), data_mod['RH'][:, h_idx] * 100,
                     linewidth=2, fmt='-', ls='-')

    # plot reference line to show where profile lies
    ax.plot_date([t, t], [0, 100], color='black', ls='--', fmt='-')

    # prettify
    ax.set_ylabel(r'$\mathrm{RH\/\/[\%]}$')
    ax.set_ylim([20, 80])
    ax.set_xlim([data_mod['time'][0], data_mod['time'][-1]])
    # ax.set_prop_cycle(cycler('color', ['b', 'c', 'm', 'y', 'k']))
    ax.xaxis.set_major_formatter(DateFormatter('%H:%M'))
    ax.set_xlabel('Time [HH:MM]')
    eu.add_at(ax, 'd)', loc=2)

    return ax

def aer_ts_plot(ax, pm10_obs, mod_data, t):


    # plot pm10 obs profile plots

    for site_obs, data_obs in pm10_obs.iteritems():
        ax.plot_date(date2num(data_obs['time']),data_obs['pm_10'], label=site_obs, fmt='-')

    # plot model for a single height near the surface (5 m)
    # height index to plot
    h_idx = 0

    for site_mod, data_mod in mod_data.iteritems():
        # ax.plot_date(date2num(data_mod['time']), data_mod['RH'][t, h_idx] * 100,
        #              label=str(data_mod['level_height'][h_idx]) + ' m model', linewidth=2, fmt='-', ls='-')
        ax.plot_date(date2num(data_mod['time']), data_mod['aerosol_for_visibility'][:, h_idx],
                     linewidth=2, fmt='-', ls='-')

    # plot reference line to show where profile lies
    ax.plot_date([t, t], [0, 100], color='black', ls='--', fmt='-')

    # prettify
    ax.set_ylabel(r'$m \/\mathrm{[\mu g\/ kg^{-1}]}$')
    ax.set_ylim([0, 100])
    ax.set_xlim([data_mod['time'][0], data_mod['time'][-1]])
    # ax.set_prop_cycle(cycler('color', ['b', 'c', 'm', 'y', 'k']))
    ax.xaxis.set_major_formatter(DateFormatter('%H:%M'))
    ax.set_xlabel('Time [HH:MM]')
    eu.add_at(ax, 'e)', loc=2)

    return ax



def main():

    # ==============================================================================
    # Setup
    # ==============================================================================

    # which modelled data to read in
    model_type = 'UKV'

    # model resolution
    if model_type == 'UKV':
        res = '1p5km'
    elif model_type == '333m':
        res = '0p3km'

    # directories
    maindir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/'
    datadir = maindir + 'data/'
    savedir = maindir + 'figures/' + model_type + '/'

    # data
    ceilDatadir = datadir + 'L1/'
    modDatadir = datadir + model_type + '/'
    rhDatadir = datadir + 'L1/'
    aerDatadir = datadir + 'LAQN/'

    # sites with height pairs
    # height above sea level (LUMA metadata)- height above ground (DTM) = height above surface
    # surface height taken from DTM (Grimmond and ...) - Digital Terrain Model (surface only, no buildings)
    # ToDo IMU needs updating!!!!
    # ToDo RH heights needs filling in too!
    site_bsc = {'CL31-A_BSC_IMU': 79.0 - 14.7, 'CL31-B_BSC_RGS': 28.1 - 19.4, 'CL31-C_BSC_MR': 32.0 - 27.5,
                'CL31-D_BSC_NK': 27.0 - 23.2}
    site_rh = {'WXT_KSSW': 0, 'Davis_BCT': 0, 'Davis_BFCL': 0, 'Davis_BGH': 0, 'Davis_IMU': 0, 'Davis_IML': 0}
    site_aer = {'PM10_RGS': 23.0 - 19.4, 'PM10_MR': 32.0 - 27.5, 'PM10_NK': 26.0 - 23.2}
    # sites = site_bsc + site_rh + site_aer

    # ceil_height_offset = [79.0 - 14.7, 28.1 - 19.4, 32.0 - 27.5, 27.0 - 23.2]
    aer_heights = [23.0 - 19.4, 32.0 - 27.5, 26.0 - 23.2]
    rh_heights = [0, 0, 0, 0, 0, 0]

    # days to loop between
    dayStart = dt.datetime(2016, 05, 03)
    dayEnd = dt.datetime(2016, 05, 05)


    # ==============================================================================
    # Read and process modelled data
    # ==============================================================================

    # 1. Read Ceilometer metadata
    # ----------------------------
    ceil_metadata = read_ceil_metadata(datadir)


    # datetime range to iterate over
    if dayStart != dayEnd:
        days_iterate = eu.date_range(dayStart, dayEnd, 1, 'days')
    else:
        days_iterate = [dayStart]

    for day in days_iterate:


        # 2. Read UKV forecast in
        # -----------------------

        # date string
        dateStr = day.strftime('%Y%m%d')

        # temp filename for the day
        mod_filename = 'extract_prodm_op_ukv_' + dateStr + '_21_full.nc'

        # concatenate the names
        mod_fname = modDatadir + mod_filename

        # Read in the modelled data for London
        mod_all_data = eu.netCDF_read(mod_fname)

        # extract MURK aerosol and calculate RH for each of the sites in the ceil metadata
        # (can be different locations to sites_bsc)
        mod_data = mod_site_extract_calc(ceil_metadata, mod_all_data, model_type, res, day)


        # 3. Read WXT and Davis
        # ------------------------

        # ToDo only get the nearest idx to the hour, like the bsc obs
        # all RH obs for the two main days, for all sites, read in at once
        rh_obs = read_all_rh_obs(day, site_rh, rhDatadir, mod_data)


        # 4. Read PM10
        # ------------------------

        # all the pm10 obs for ALL THREE DAYS,
        # ToDo only get the nearest idx to the hour, like the bsc obs
        pm10_obs = read_all_pm10_obs(dayStart, dayEnd, site_aer, aerDatadir, mod_data)


        # set up plot
        # fig = plt.figure(figsize=(10, 7))


        # for site in site_bsc:

        # 5. Read ceilometer backscatter
        # --------------------------------

        # ToDo transpose so time = x, height = y

        bsc_obs = read_ceil_obs(day, site_bsc, ceilDatadir, mod_data)


        # do any stats on time reduced data here...


        # ==============================================================================
        # Plotting
        # ==============================================================================

        for t in np.arange(bsc_obs['CL31-A_BSC_IMU']['time']):

            # current hour
            t_hr_str = bsc_obs['CL31-A_BSC_IMU']['time'][t].strftime("%Y-%m-%d %H:%M:%S")

            # set up plot
            fig = plt.figure(figsize=(10, 7))

            # set up subplots
            ax1 = plt.subplot2grid((3, 3), (0, 0))  # bsc
            ax2 = plt.subplot2grid((3, 3), (0, 1))  # aer_mmr
            ax3 = plt.subplot2grid((3, 3), (0, 2))  # rh
            ax4 = plt.subplot2grid((3, 3), (1, 0), colspan=3)  # rh time series
            ax5 = plt.subplot2grid((3, 3), (2, 0), colspan=3)  # aer + pm10 time series

            # ceilometer profile
            ax1 = bsc_profile_plot(ax1, mod_data, bsc_obs, t)

            # RH profile plot
            ax2 = rh_profile_plot(ax2, site_rh, rh_obs, mod_data, t)

            # PM10 profile plot
            ax3 = aer_profile_plot(ax3, site_aer, pm10_obs, mod_data, t)

            # RH time series
            ax4 = rh_ts_plot(ax4, rh_obs, mod_data, t)

            # PM10 time series
            ax5 = aer_ts_plot(ax5, pm10_obs, mod_data, t)

            # final prettifying for figure
            fig.suptitle('BSC: ' + t.strftime("%Y-%m-%d %H:%M:%S"), fontsize=12)
            plt.tight_layout(h_pad=0)

            # adjust axis so the ledgends can be placed on the side
            fig.subplots_adjust(top=0.95, right=0.8)

            # legends
            han1, lab1 = ax1.get_legend_handles_labels()
            han2, lab2 = ax2.get_legend_handles_labels()
            han3, lab3 = ax3.get_legend_handles_labels()

            # lab2[1] = 'RH observations'
            # lab3 = ['PM' + r'$_{10}$' + ' observations']
            # han3 = [han3[3]]  # ignore the other lines there

            # ax3.legend(han1 + han2 + han3, lab1 + lab2 + lab3, fontsize=8, bbox_to_anchor=(1.07, 1), loc=2,
            #            borderaxespad=0.0)
            # ax4.legend(fontsize=8, bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.0)
            # ax5.legend(fontsize=8, bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.0)

            # save the figure
            plt.savefig(savedir + 'profiles/' + 'BscRhAer_' + model_type + '_' + t.strftime("%Y%m%d_%H%M") + '.png')


    print 'END PROGRAM'

    plt.close('all')

if __name__ == '__main__':
    main()