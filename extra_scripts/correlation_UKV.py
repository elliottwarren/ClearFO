"""
Correlated UKV variables between sites

Created by Elliott Tues 29/06/17
"""

import matplotlib.pyplot as plt

import numpy as np
import datetime as dt
from scipy.stats import spearmanr
from matplotlib.dates import DateFormatter

from forward_operator import FOUtils as FO
from forward_operator import FOconstants as FOcon

# def create_stats_entry(site_id, statistics={}):
#
#     """
#     Define or expand the almighty statistics array
#
#     :param site_bsc:
#     :param mbe_limit_max:
#     :param mbe_limit_step:
#     :return: statistics (dict)
#
#     statistics[site]['r'] = [...]
#     statistics[site]['MBE'] = {'0-500': ..., '500-1000': ...}
#     statistics[site]['time'] = [...]
#     """
#
#     # statistics will be grouped based on hour, so create a simple hourly array [0 ... 23]
#     hrs = np.arange(0, 24)
#
#     # Structure of statistics:
#     # statistics[site]['r'] = [...]
#     # statistics[site]['MBE'] = {'0-500': ..., '500-1000': ...}
#     # statistics[site]['time'] = [...]
#
#     if site_id not in statistics:
#
#         # define site based lists to store the correlation results in
#         statistics[site_id] = {'r': [], 'p': [],
#                                'diff': [],
#                                'aer_diff': [],
#                                'aer_mod': [],
#                                'aer_obs': [],
#                                'rh_diff': [],
#                                'rh_mod': [],
#                                'rh_obs': [],
#                                'back_point_diff': [],
#                                'RMSE': [],
#                                'MBE': [],
#                                'hr': []}
#
#         # for key in statistics[site_id].iterkeys():
#         #     for hr in hrs:
#         #         statistics[site_id][key][str(hr)] = []
#
#     return statistics

def dateList_to_datetime(dayList):

    """ Convert list of string dates into datetimes """

    datetimeDays = []

    for d in dayList:

        datetimeDays += [dt.datetime(int(d[0:4]), int(d[4:6]), int(d[6:8]))]

    return datetimeDays

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
    savedir = maindir + 'figures/' + model_type + '/clearSkyPeriod/'

    # data
    ceilDatadir = datadir + 'L1/'
    modDatadir = datadir + model_type + '/'
    rhDatadir = datadir + 'L1/'
    aerDatadir = datadir + 'LAQN/'

    # statistics to run
    pm10_stats = True
    rh_stats = False

    # # instruments and other settings
    #site_rh = FOcon.site_rh
    #site_rh = {'WXT_KSSW': 50.3}
    #rh_instrument = site_rh.keys()[0]

    #site = 'NK'
    #ceil_id = 'CL31-D'
    #ceil = ceil_id + '_' + site

    # instruments and other settings
    #site_rh = FOcon.site_rh
    #site_rh = {'Davis_IMU': 72.8}
    #rh_instrument = site_rh.keys()[0]


    site = 'MR'
    ceil_id = 'CL31-C'
    # ceil = ceil_id + '_BSC_' + site
    ceil = ceil_id + '_' + site

    site_bsc = {ceil: FOcon.site_bsc[ceil]}
    # site_bsc = {ceil: FOcon.site_bsc[ceil], 'CL31-E_BSC_NK': 27.0 - 23.2}

    if pm10_stats == True:
        site_aer = {'PM10_'+site: FOcon.site_aer['PM10_'+site]}

    if rh_stats == True:
        site_rh = {'WXT_KSSW': 50.3}
        rh_instrument = site_rh.keys()[0]

    site_bsc_colours = FOcon.site_bsc_colours


    # MR calib days
    #daystrList = ['20150414', '20150415', '20150421', '20150611', '20160504', '20160823', '20160911', '20161125',
    #              '20161129', '20161130', '20161204']

    # high pollution case study day
    daystrList = ['20160119']


    days_iterate = dateList_to_datetime(daystrList)

    # variable to compare (whatever the mod_site_extract_calc function has named them)
    variable = 'Q_H'

    # define statistics dictionary
    statistics = {}
    sampleSize = 0 # add to this

    # store RH(z = z_i) for each
    var_MR = []
    var_KSS45W = []

    time_MR = []
    time_KSS45W = []

    # ==============================================================================
    # Read data
    # ==============================================================================

    # Read Ceilometer metadata

    # ceilometer list to use
    ceilsitefile = 'UKV_correlatesites.csv'
    ceil_metadata = FO.read_ceil_metadata(datadir, ceilsitefile)

    for day in days_iterate:

        print 'day = ' + day.strftime('%Y-%m-%d')

        # Read UKV forecast and automatically run the FO

        # extract MURK aerosol and calculate RH for each of the sites in the ceil metadata
        # reads all london model data, extracts site data, stores in single dictionary
        mod_data = FO.mod_site_extract_calc(day, ceil_metadata, modDatadir, model_type, res, 910, version=0.2, allvars=True)

        # idx heights to pull z
        zidx_mr = 1
        zidx_kss45w = 3

        # actual heights of both
        z_mr = mod_data['MR']['level_height'][zidx_mr]
        z_kss45w = mod_data['MR']['level_height'][zidx_kss45w]

        # store these for correlating later
        var_MR += [mod_data['MR'][variable][:, zidx_mr]]
        var_KSS45W += [mod_data['KSS45W'][variable][:, zidx_kss45w]]

        # store their times
        time_MR += [mod_data['MR']['time']]
        time_KSS45W += [mod_data['KSS45W']['time']]

    var_MR = np.hstack(var_MR)
    var_KSS45W = np.hstack(var_KSS45W)

    time_MR = np.hstack(time_MR)
    time_KSS45W = np.hstack(time_KSS45W)

    # do correlation
    corr = {}
    corr['r'], corr['p'] = spearmanr(var_MR, var_KSS45W, nan_policy='omit')
    n = len(var_MR)

    print 'r = ' + str(corr['r'])
    print 'p = ' + str(corr['p'])

    # SCATTER
    # plt.plot([-250.0, 2000.0], [-250.0, 2000.0], color='black')
    # plt.scatter(RH_MR*100.0, RH_KSS45W*100.0)
    # plt.xlabel('MR')
    # plt.ylabel('KSS45W')
    # plt.suptitle(str(len(daystrList)) + ' sample days; n=' + str(n) + '; z_mr=' + str(z_mr) + 'm ; z_kss45w=' + str(z_kss45w) + 'm\n'
    #               'r=' + str(corr['r']) + '; p=' + str(corr['p']))
    #
    # fname = 'corr_'+variable+'_MR' + '_' + str(z_mr) + 'm-KSS45W_' + str(z_kss45w) + 'm.png'
    # plt.savefig(savedir + 'point_diff/picking_heights/' + fname)

    # LINE PLOT
    plt.plot_date(time_MR, var_MR, '-', label='MR')
    plt.plot_date(time_KSS45W, var_KSS45W, '-', label='KSS45W')
    plt.xlabel('time')
    plt.ylabel('Q_H')
    ax = plt.gca()
    ax.xaxis.set_major_formatter(DateFormatter('%H:%M'))
    plt.suptitle(str(len(daystrList)) + ' sample days; n=' + str(n) + '; z_mr=' + str(z_mr) + 'm ; z_kss45w=' + str(z_kss45w) + 'm\n'
                  'r=' + str(corr['r']) + '; p=' + str(corr['p']))
    plt.legend()

    fname = 'corr_'+variable+'_MR' + '_' + str(z_mr) + 'm-KSS45W_' + str(z_kss45w) + 'm_lineplot.png'
    plt.savefig(savedir + 'point_diff/picking_heights/' + fname)

    # plt.close('all')

















    return

if __name__ == '__main__':
    main()























print 'END PROGRAM'
