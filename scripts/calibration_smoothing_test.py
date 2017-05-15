"""
Calibration smoothing to test what the best way to use the calibration data is

Created on 15/05/17 by EW
"""

import matplotlib.pyplot as plt
from matplotlib.dates import date2num
from matplotlib.dates import DateFormatter

import numpy as np
import datetime as dt

import ellUtils as eu
from forward_operator import FOUtils as FO
from forward_operator import FOconstants as FOcon
import calibrations as cal

def read_cal():

    """Read calibration data"""

    sites = {'KSS45W':cal.KSS45W,'RGS':cal.RGS,'MR':cal.MR}

    calib = {}
    for key, item in sites.iteritems():
        # read in calibration

        calib[key] = {'Dates': item.Dates, 'wv_cal': np.array(item.C_modes_wv)}

    return calib

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
    savedir = maindir + 'figures/calibration_test/'

    # read calibration data
    calib = read_cal()


    for site, data in calib.iteritems():

        # turn dates into datetimes
        dates = dateList_to_datetime_calib_format(data['Dates'])

        # linearly interpolat the nans
        calibration = eu.linear_interpolation(data['wv_cal'])

        # mask raw data
        raw = np.ma.masked_where(np.isnan(data['wv_cal']), data['wv_cal'])

        # moving average the data
        ma_7 = eu.moving_average(calibration, 7)
        ma_10 = eu.moving_average(calibration, 10)
        ma_30 = eu.moving_average(calibration, 30)

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


































    print 'END PROGRAM'