"""
Calibration smoothing to test what the best way to use the calibration data is

Created on 15/05/17 by EW
"""

import matplotlib.pyplot as plt
from matplotlib.dates import date2num
from matplotlib.dates import DateFormatter
from pm_RH_vs_mod_obs_backscatter import dateList_to_datetime

import numpy as np
import datetime as dt

import ellUtils as eu
import ceilUtils as ceil
from forward_operator import FOUtils as FO
from forward_operator import FOconstants as FOcon
import calibrations as cal

def read_cal():

    """Read calibration data"""

    sites = {'KSS45W':cal.KSS45W,'RGS':cal.RGS,'MR':cal.MR}

    calib = {}
    for key, item in sites.iteritems():
        # read in calibration

        calib[key] = {'Dates': item.Dates, 'wv_cal': np.array(item.C_modes_wv),
                      'samples': item.profile_total}

    return calib

def read_periods(perioddatadir, site_id):

    """ Read in ceilometer period data, e.g. firmware, hardware, window transmission and pulse changes"""

    datapath = perioddatadir + 'ceilometer_periods_' + site_id + '.csv'
    # read in ceilometer periods
    # periods = np.genfromtxt(datapath, delimiter=',')
    periods_raw = eu.csv_read(datapath)

    # sort data out into dictionary
    periods = {}
    for i in np.arange(len(periods_raw[0])):
        periods[periods_raw[0][i]] = [j[i] for j in periods_raw[1:]]

    # convert date strings into datetimes
    periods['Date'] = [dt.datetime.strptime(d, '%d/%m/%Y') for d in periods['Date']]

    return periods

def read_window_trans(site, ceildatadir):

    # site id (short)
    site_id = site.split('_')[-1]
    ceil_id = site.split('_')[0]

    # get filename
    fname1 = ceildatadir + ceil_id + '_CCW30_' + site_id + '_2015_15min.nc'
    fname2 = ceildatadir + ceil_id + '_CCW30_' + site_id + '_2016_15min.nc'

    # read window transmission
    transmission1 = ceil.netCDF_read_CCW30(fname1, var_type='transmission')
    transmission2 = ceil.netCDF_read_CCW30(fname2, var_type='transmission')

    # merge dictionaries
    # transmission = eu.merge_dicts(transmission1, transmission1)

    transmission = {}
    for var, var_data in transmission1.iteritems():

        if var == 'time':
            transmission[var] = var_data + transmission2[var]
        elif (var == 'level_height') | (var == 'height'):
            transmission[var] = var_data
        else:
            # transmission[var] = np.ma.hstack((transmission1[var], np.ma.masked_array(transmission2[var],mask=np.zeros(len(transmission2[var])))))
            transmission[var] = np.hstack((transmission1[var], transmission2[var]))



    return transmission

def read_pulse(site, ceildatadir):

    # site id (short)
    site_id = site.split('_')[-1]
    ceil_id = site.split('_')[0]

    # get filename
    fname1 = ceildatadir + ceil_id + '_CCW30_' + site_id + '_2015_15min.nc'
    fname2 = ceildatadir + ceil_id + '_CCW30_' + site_id + '_2016_15min.nc'

    # read window transmission
    transmission1 = ceil.netCDF_read_CCW30(fname1, var_type='pulse')
    transmission2 = ceil.netCDF_read_CCW30(fname2, var_type='pulse')

    # merge dictionaries
    # transmission = eu.merge_dicts(transmission1, transmission1)

    transmission = {}
    for var, var_data in transmission1.iteritems():

        if var == 'time':
            transmission[var] = var_data + transmission2[var]
        elif (var == 'level_height') | (var == 'height'):
            transmission[var] = var_data
        else:
            # transmission[var] = np.ma.hstack((transmission1[var], np.ma.masked_array(transmission2[var],mask=np.zeros(len(transmission2[var])))))
            transmission[var] = np.hstack((transmission1[var], transmission2[var]))


    return transmission

def remove_events(periods, window_trans_daily):

    # remove days in window_trans_daily when events happened (firmware change, window cleaning etc)
    period_event_days = np.array([i.date() for i in periods['Date']])

    for day in period_event_days:
        idx = np.where(np.array(window_trans_daily['dates']) == day)

        # turn daily data into nan as it is unreliable
        for key in window_trans_daily.iterkeys():
            if key != 'dates':
                window_trans_daily[key][idx] = np.nan

    return window_trans_daily

def calc_daily_window_trans(window_trans, calib_dates_days, calib_raw):

    """
    Calculate the daily maximum window transmission, and pair it with the calibration coefficient that factors in
    water vapour.
    :param window_trans:
    :param calib_dates_days:
    :param calib_raw:
    :return:
    """

    # find all unique days in window_trans time
    # loop through days
    # find all dates with that day
    # max(c) and store

    dayMin = window_trans['time'][0].date()
    dayMax = window_trans['time'][-1].date()
    # np.array([i.date() for i in window_trans['time']])

    window_trans_dates = np.array([i.date() for i in window_trans['time']])

    daysRange = eu.date_range(dayMin, dayMax, 1, 'days')
    window_trans_daily = {'dates': daysRange, 'max_window': np.empty(len(daysRange)),
                          'c_wv': np.empty(len(daysRange)), 'samples': np.empty(len(daysRange))}
    window_trans_daily['c_wv'][:] = np.nan
    window_trans_daily['max_window'][:] = np.nan
    window_trans_daily['samples'][:] = np.nan

    for day, dayIdx in zip(daysRange, np.arange(len(daysRange))):

        # idx for all days on this day
        winIdx = np.where(window_trans_dates == day)

        # store maximum window transmission for the day
        if winIdx[0].size != 0:
            window_trans_daily['max_window'][dayIdx] = np.nanmax(window_trans['transmission'][winIdx])

        # get c
        cIdx = np.where(calib_dates_days == day)

        # store c for the day
        if cIdx[0].size != 0:
            window_trans_daily['c_wv'][dayIdx] = calib_raw[cIdx]


        # get c
        sIdx = np.where(calib_dates_days == day)

        # store c for the day
        if sIdx[0].size != 0:
            window_trans_daily['samples'][dayIdx] = calib_raw[sIdx]


    return window_trans_daily

def plot_smooth(dates, calibration, ma_7, ma_10, ma_30):

    """
    Plot the smoothed curves...
    :param dates:
    :param calibration:
    :param ma_7:
    :param ma_10:
    :param ma_30:
    :return:
    """

    # plot
    plt.figure()
    plt.plot_date(dates, calibration, label='linear interp.', fmt='-', ls='-')
    # plt.plot_date(dates, raw, label='raw', fmt='-', ls='-')
    plt.plot_date(dates, ma_7, label='7 day', fmt='-', ls='-')
    plt.plot_date(dates, ma_10, label='10 day', fmt='-', ls='-')
    plt.plot_date(dates, ma_30, label='30 day', fmt='-', ls='-')

    ax = plt.gca()
    ax.xaxis.set_major_formatter(DateFormatter('%m/%y'))
    ax.set_xlim([dt.datetime(2015, 2, 1), dt.datetime(2016, 6, 1)])
    ax.set_xlabel('Time [mm/YY]')
    ax.set_ylabel('water vapour C_modes')
    plt.suptitle(site)
    plt.legend()
    plt.savefig(savedir + 'calibration_wv_' + site + '.png')

    return

def plot_cal_wind_pulse_timeseries(calib_dates, calib_raw, window_trans, pulse, savedir, site, clear_days):

        """ Plot time series of points, of raw calibration, window transmission and pulse energy"""

        fig = plt.figure()
        ax1 = fig.add_subplot(1, 1, 1)
        ax2 = ax1.twinx()
        marker = 4
        ax1.plot_date(calib_dates, calib_raw, label='raw', color='g', markersize=marker)
        ax2.plot_date(window_trans['time'], window_trans['transmission'], label='window_trans', markersize=marker)
        ax2.plot_date(pulse['time'], pulse['pulse'], label='pulse', markersize=marker)

        # put vertical lines in for each of the clear days
        for day in clear_days:

            ax2.plot_date([day, day], [-100, 200], fmt='-', ls='--', color='black')

        # Prettify
        ax2.set_ylim([0.0, 100.0])
        ax = plt.gca()
        ax.xaxis.set_major_formatter(DateFormatter('%m/%y'))
        ax.set_xlim([window_trans['time'][0], window_trans['time'][-1]])
        ax.set_xlabel('Time [mm/YY]',labelpad=0)
        ax1.set_ylabel('water vapour C')
        ax2.set_ylabel('Percentage')
        plt.suptitle(site)
        plt.legend()
        plt.savefig(savedir + 'rawCalib_trans_pulse_' + site + '.png')

def scatter_window_c_all(window_trans_daily, site, savedir, period=False):

    # scatter plot for the calibration data. daily max window trans vs calibration factor.
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.scatter(window_trans_daily['max_window'], window_trans_daily['c_wv'], color='g')

    # add a linear fit to the plot
    m, b = eu.linear_fit_plot(window_trans_daily['max_window'], window_trans_daily['c_wv'], ax, ls='-', color='black')

    # string for equation
    eqstr = 'y = %.2fx + %.2f' % (m, b)

    # Prettify
    ax.set_ylim([np.nanmin(window_trans_daily['c_wv']), np.nanmax(window_trans_daily['c_wv'])])
    ax.set_xlim([np.nanmin(window_trans_daily['max_window']), np.nanmax(window_trans_daily['max_window'])])
    # ax.xaxis.set_major_formatter(DateFormatter('%m/%y'))
    # ax.set_xlim([dt.datetime(2015, 2, 1), dt.datetime(2016, 6, 1)])
    ax.set_xlabel('Window transmission (daily max) [%]')
    ax.set_ylabel('water vapour C')
    plt.suptitle(site + '; ' + eqstr)
    plt.legend()
    plt.savefig(savedir + 'rawCalib_vs_maxWindow_daily_' + site + '.png')


    return

def scatter_window_c_period(window_trans_daily, site, savedir, startDay, endDay):

    # scatter plot for the calibration data. daily max window trans vs calibration factor.
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.scatter(window_trans_daily['max_window'], window_trans_daily['c_wv'], c = window_trans_daily['sample'], color='g')

    # add a linear fit to the plot
    m, b = eu.linear_fit_plot(window_trans_daily['max_window'], window_trans_daily['c_wv'], ax, ls='-', color='black')

    # string for equation
    eqstr = 'y = %.2fx + %.2f' % (m, b)

    # period range as a string for saving
    period_range_str = startDay.strftime('%Y%m%d') + '-' + endDay.strftime('%Y%m%d')

    # Prettify
    ylim = [np.nanmin(window_trans_daily['c_wv']), np.nanmax(window_trans_daily['c_wv'])]
    xlim = [np.nanmin(window_trans_daily['max_window']), np.nanmax(window_trans_daily['max_window'])]

    # some periods only have 100 % window transmission, leading to weird limits being set on the xlim.
    # hence this deals with that edge case
    if xlim[0] == xlim[1]:
        xlim[0] = 90.0
        xlim[1] = 100.0

    ax.set_ylim(ylim)
    ax.set_xlim(xlim)
    # ax.xaxis.set_major_formatter(DateFormatter('%m/%y'))
    # ax.set_xlim([dt.datetime(2015, 2, 1), dt.datetime(2016, 6, 1)])
    ax.set_xlabel('Window transmission (daily max) [%]')
    ax.set_ylabel('water vapour C')
    plt.suptitle(site + '; ' + eqstr + '\n' + period_range_str)
    plt.legend()

    plt.savefig(savedir + '/' + site + '/rawCalib_vs_maxWindow_daily_' + period_range_str + '_' + site + '.png')

    return

def dateList_to_datetime_calib_format(dayList):

    """ Convert list of string dates into datetimes """

    datetimeDays = []

    for d in dayList:

        datetimeDays += [dt.datetime(int(d[0:4]), int(d[5:7]), int(d[8:10]))]

    return datetimeDays

if __name__ == '__main__':

    # ==============================================================================
    # Setup
    # ==============================================================================

    # directories
    maindir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/'
    datadir = maindir + 'data/'
    ceildatadir = datadir + 'L1/'
    perioddatadir = datadir + 'Calibrations_for_LUMO_Ceilometers/'
    savedir = maindir + 'figures/calibration_test/'

    site_bsc = {'CL31-B_RGS': 28.1 - 19.4, 'CL31-C_MR': 32.0 - 27.5, 'CL31-A_KSS45W': 64.3}

    # clear sky days
    daystrList = ['20150414', '20150415', '20150421', '20150611', '20160504', '20160823', '20160911', '20161125',
                  '20161129', '20161130', '20161204']

    clear_days = dateList_to_datetime(daystrList)

    # read calibration data
    calib_all = read_cal()


    # for each site
    for site in site_bsc:

        # get site id
        site_id = site.split('_')[-1]

        # this site calibration data
        calib = calib_all[site_id]

        # read periods
        periods = read_periods(perioddatadir, site_id)

        # read transmission
        window_trans = read_window_trans(site, ceildatadir)


        # read pulse energy
        pulse = read_pulse(site, ceildatadir)

        # turn dates into datetimes
        calib_dates = dateList_to_datetime_calib_format(calib['Dates'])
        calib_dates_days = np.array([i.date() for i in calib_dates])

        # linearly interpolat the nans
        # calibration = eu.linear_interpolation(data['wv_cal'])

        # mask raw data
        calib_raw = np.ma.masked_where(np.isnan(calib['wv_cal']), calib['wv_cal'])

        # create  time series of dailymax window transmission
        window_trans_daily = calc_daily_window_trans(window_trans, calib_dates_days, calib_raw)

        # remove event days from window transmission (firmware, hardware, cleaing etc)
        window_trans_daily = remove_events(periods, window_trans_daily)


        # remove transmission enteries in periods
        # find all NONE transmission enteries
        idx = np.where(np.array(periods['Type']) != 'Transmission')[0]

        if idx[0].size != 0:
            for key in periods.iterkeys():
                periods[key] = [periods[key][i] for i in idx]



        # find all window, c and sample in current date range
        periods_days = np.array([i.date() for i in periods['Date']])

        # loop through start and end day pairs
        for startDay, endDay in zip(periods_days[:-1], periods_days[1:]):

            # get idx position for this period
            # for some reason the idx array is within another array... is this to do with having two arguments?
            idx = np.where((np.array(window_trans_daily['dates']) > startDay) &
                                     (np.array(window_trans_daily['dates']) < endDay))

            # if window_trans_daily data exists for this date, extract and plot the data
            if idx[0].size != 0:

                # extract data for this period
                window_trans_daily_period = {}
                for key in window_trans_daily.iterkeys():
                    window_trans_daily_period[key] = [window_trans_daily[key][i] for i in idx[0]]

                # if there is data for window transmission and c then plot
                if np.any(np.isfinite(window_trans_daily_period['c_wv'])):
                    # use the scatter function to make the graphs
                    scatter_window_c_period(window_trans_daily_period, site, savedir, startDay, endDay)

        # plot window vs c, coloured with sample size



        # ----------------
        # moving average the data
        #ma_7 = eu.moving_average(calibration, 7)
        #ma_10 = eu.moving_average(calibration, 10)
        #ma_30 = eu.moving_average(calibration, 30)
        #
        # plot the smoothed data
        # plot_smooth(dates, calibration, ma_7, ma_10, ma_30)
        # ----------------

        # time-series point plot of calibration, window trans and pulse energy.
        # plot_cal_wind_pulse_timeseries(calib_dates, calib_raw, window_trans, pulse, savedir, site, clear_days)

        # scatter all the data for a site. Max daily window trans verses calibration coeff
        # scatter_window_c_all(window_trans_daily, site, savedir)



    print 'END PROGRAM'