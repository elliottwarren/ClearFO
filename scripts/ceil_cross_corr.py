"""
Cross correlate the ceilometer observations with and without the calibration values

Created by Elliott Mon 22/05/2017
"""


import matplotlib.pyplot as plt
from matplotlib.dates import date2num
from matplotlib.dates import DateFormatter

import numpy as np
import datetime as dt
from scipy.stats import spearmanr

import ellUtils as eu
from forward_operator import FOUtils as FO
from forward_operator import FOconstants as FOcon

def dateList_to_datetime(dayList):

    """ Convert list of string dates into datetimes """

    datetimeDays = []

    for d in dayList:

        datetimeDays += [dt.datetime(int(d[0:4]), int(d[4:6]), int(d[6:8]))]

    return datetimeDays

def plot_day(bsc_avg_day, bsc_obs, day, savedir, calib):

    """Plot the average vertical profile for each sensor, for the current day"""

    # plot each day separately
    fig = plt.figure(figsize=(6, 3.5))
    ax = plt.subplot2grid((1, 1), (0, 0))

    for site, avg in bsc_avg_day.iteritems():
        ax.plot(avg, bsc_obs[site]['height'], label=site)

    ax.set_xlabel('beta')
    ax.set_ylabel('Height [m]')
    ax.set_xlim([0.0, 5.0e-6])
    plt.legend()
    dayStr = day.strftime('%Y%m%d')
    plt.suptitle(dayStr)

    if calib == True:
        plt.savefig(savedir + 'avg_Vert_profiles_calib_' + dayStr + '.png')  # filename
    else:
        plt.savefig(savedir + 'avg_Vert_profiles_uncalib_' + dayStr + '.png')  # filename
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
    savedir = maindir + 'figures/calibration_test/calib_check/'

    # data
    ceilDatadir = datadir + 'L1/'
    modDatadir = datadir + model_type + '/'
    rhDatadir = datadir + 'L1/'
    aerDatadir = datadir + 'LAQN/'

    # instruments and other settings
    site_bsc = FOcon.site_bsc
    site_rh = FOcon.site_rh
    site_aer = FOcon.site_aer
    site_bsc_colours = FOcon.site_bsc_colours

    # day start and end of the MAIN DAYS, inclusively(forecast day + 1)
    # KSS45W days
    daystrList = ['20150414', '20150415', '20150421', '20150611']

    bsc_avg_day = {}

    calib=True

    # ==============================================================================
    # Read data
    # ==============================================================================

    # 1. Read Ceilometer metadata
    # ----------------------------
    ceil_metadata = FO.read_ceil_metadata(datadir)

    # datetime range to iterate over
    days_iterate = dateList_to_datetime(daystrList)

    for day in days_iterate:

        print 'day = ' + day.strftime('%Y-%m-%d')

        # Read ceilometer backscatter
        # --------------------------------

        bsc_obs = FO.read_ceil_obs_only(day, site_bsc, ceilDatadir, calib=calib)

        # set up the average dicts
        if day == days_iterate[0]:

            bsc_avg_day_all = {}
            for site in bsc_obs.iterkeys():
                bsc_avg_day_all[site] = []

            bsc_avg_total = {}
            for site in bsc_obs.iterkeys():
                bsc_avg_total[site] = []

        # average in the vertical
        # store data in new array.

        for site, bsc_site_obs in bsc_obs.iteritems():

            # short site id that matches the model id
            site_id = site.split('_')[-1]

            # average ceil data in vertical
            bsc_avg_day[site] = np.nanmean(bsc_obs[site]['backscatter'], axis=0)
            bsc_avg_day_all[site] += [bsc_avg_day[site]]

        # plot avergae profile for the current day
        plot_day(bsc_avg_day, bsc_obs, day, savedir, calib)


    # create the average of the average profiles
    # first turn list of lists into a single numpy array, then average across an axis
    for site in bsc_obs.iterkeys():
        bsc_avg_total[site] = np.array([np.array(xi) for xi in bsc_avg_day_all[site]])

        bsc_avg_total[site] = np.nanmean(bsc_avg_total[site], axis=0)

    # plot the average of average profiles
    fig = plt.figure(figsize=(6, 3.5))
    ax = plt.subplot2grid((1, 1), (0, 0))

    for site, avg_total in bsc_avg_total.iteritems():
        ax.plot(avg_total, bsc_obs[site]['height'], label=site)

    ax.set_xlabel('beta')
    ax.set_ylabel('Height [m]')
    ax.set_xlim([0.0, 5.0e-6])
    plt.legend()
    dayStr = day.strftime('%Y%m%d')

    if calib == True:
        plt.savefig(savedir + 'avg_ofAvg_Vert_profiles_calib.png')  # filename
    else:
        plt.savefig(savedir + 'avg_ofAvg_Vert_profiles_uncalib.png')  # filename

    plt.close(fig)


    return

if __name__ == '__main__':
    main()