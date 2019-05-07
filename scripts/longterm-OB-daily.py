"""
Hopefully the start of something good.
Make daily plots with hourly statistics with standard deviations and everything

Created by Elliott Thurs 27th April 2017
"""

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

import numpy as np
import datetime as dt
from scipy.stats import spearmanr

from copy import deepcopy

import ellUtils as eu
from forward_operator import FOUtils as FO
from forward_operator import FOconstants as FOcon


def unique_pairs(obs_idx, diff):

    """
    Find range that excludes duplicate occurances. Keeps the pair with the smallest height difference and removes
    the rest.

    :param obs_idx:
    :param diff:
    :return: unique_pairs_range

    At this point, the two arrays are like:
    obs_idx = [0, 0, 0, 1, 3, 5, .... 769, 769, 769]
    mod_idx = [0, 1, 2, 3, 4, 4, .... 67,  68,  69 ]
    By finding the unique pairs index array for obs_idx, the same array can be used
    on the mod_idx, as they are already paired up and of equal lengths. E.g. from above
    0-0, 0-1, ..., 3-4, 5-4 etc.
    """

    # 1. remove start duplicates
    # -------------------------------
    # find start idx to remove duplicate pairs
    duplicates = np.where(obs_idx == obs_idx[0])[0]  # find duplicates

    if len(duplicates) > 1:
        lowest_diff = np.argmin(abs(diff[duplicates]))  # find which has smallest difference
        pairs_idx_start = duplicates[lowest_diff]  # set start position for pairing at this point
    else:
        pairs_idx_start = 0

    # 2. remove end duplicates
    # -------------------------------
    # find end idx to remove duplicate pairs
    duplicates = np.where(obs_idx == obs_idx[-1])[0]  # find duplicates
    if len(duplicates) > 1:
        lowest_diff = np.argmin(abs(diff[duplicates]))  # find which has smallest difference
        pairs_idx_end = duplicates[lowest_diff]  # set start position for pairing at this point
    else:
        pairs_idx_end = len(obs_idx)

    # create range in order to extract the unique pairs
    unique_pairs_range = np.arange(pairs_idx_start, pairs_idx_end + 1)

    return unique_pairs_range

def create_stats_entry(site_id, statistics={}):

    """
    Define or expand the almighty statistics array

    :param site_bsc:
    :param mbe_limit_max:
    :param mbe_limit_step:
    :return: statistics (dict)

    statistics[site]['r'] = [...]
    statistics[site]['MBE'] = {'0-500': ..., '500-1000': ...}
    statistics[site]['time'] = [...]
    """

    # statistics will be grouped based on hour, so create a simple hourly array [0 ... 23]
    hrs = np.arange(0, 24)

    # Structure of statistics:
    # statistics[site]['r'] = [...]
    # statistics[site]['MBE'] = {'0-500': ..., '500-1000': ...}
    # statistics[site]['time'] = [...]

    if site_id not in statistics:

        # define site based lists to store the correlation results in
        statistics[site_id] = {'r': {}, 'p': {},
                               'diff': {},
                               'RMSE': {},
                               'MBE': {},
                               'MAE': {},
                               'AE': {}}


        for hr in hrs:
            statistics[site_id]['r'][str(hr)] = []
            statistics[site_id]['p'][str(hr)] = []
            statistics[site_id]['diff'][str(hr)] = []
            statistics[site_id]['RMSE'][str(hr)] = []
            statistics[site_id]['MAE'][str(hr)] = []
            statistics[site_id]['AE'][str(hr)] = []

    return statistics

def create_stats_summary_dict(site_id, summary={}):

    """
    Define or expand the almighty summary statistics array
    summary will be a dictionary for a different statistic e.g. corr or rmse

    :param site_id:
    :param summary:
    :return: summary
    """

    # statistics will be grouped based on hour, so create a simple hourly array [0 ... 23]
    hrs = np.arange(0, 24)

    # define nanArray
    # idea is to preestablish the array size, so it can be indexed into. When the statistic is made for that hour,
    # it can be placed in the correct position. If the stat cannot be made, the idx will remain nan.
    nanArray = np.empty(24)
    nanArray[:] = np.nan

    if site_id not in summary:

        # define site based lists to store the correlation results in
        summary[site_id] = {'median': deepcopy(nanArray), 'q25': deepcopy(nanArray),'q75': deepcopy(nanArray),
                               'IQR': deepcopy(nanArray), 'mean': deepcopy(nanArray), 'stdev': deepcopy(nanArray),
                                'hrs': hrs}

    return summary

def dateList_to_datetime(dayList):

    """ Convert list of string dates into datetimes """

    datetimeDays = []

    for d in dayList:

        datetimeDays += [dt.datetime(int(d[0:4]), int(d[4:6]), int(d[6:8]))]

    return datetimeDays

def nearest_heights(mod_height, obs_height, minHeight=-np.inf, maxHeight=np.inf):

    """
    Get the nearest ceilometer height gate to each model level

    :param mod_height:
    :param obs_height:
    :return:

    obs_idx = ALL nearest gate idx
    mod_idx = idx of the model height that each obs_idx are paired to
    """


    a = np.array([eu.nearest(obs_height, i) for i in mod_height])
    values = a[:, 0]
    obs_idx = np.array(a[:, 1], dtype=int)
    diff = a[:, 2]
    mod_idx = np.arange(len(mod_height))  # mod_idx should be paired with obs_idx spots.

    # Trim off the ends of obs_idx, as UKV and obs z0 and zmax are different, leading to the same gate matching multiple ukvs
    # assumes no duplicates in the middle of the arrays, just at the end

    # At this point, variables are like:
    # obs_idx = [0, 0, 0, 1, 3, 5, .... 769, 769, 769]
    # mod_idx = [0, 1, 2, 3, 4, 4, .... 67,  68,  69 ]
    unique_pairs_range = unique_pairs(obs_idx, diff)

    # ALL unique pairs
    # Use these to plot correlations for all possible pairs, regardless of height
    obs_unique_pairs = obs_idx[unique_pairs_range]
    mod_unique_pairs = mod_idx[unique_pairs_range]
    values_unique_pairs = values[unique_pairs_range]
    diff_unique_pairs = diff[unique_pairs_range]

    # ~~~~~~~~~~~~~~~~~~~~ #


    # # hc = height cut - original using corr_max_height
    # hc_unique_pairs_range = np.where(values_unique_pairs <= corr_max_height)[0]

    # cut off below min height and above max height. If undefined, then limits set to +/- inf
    hc_unique_pairs_range = np.logical_and(values_unique_pairs >= minHeight,
                                           values_unique_pairs <= maxHeight)

    # trim off unique pairs that are below the minimum height and above the maximum height
    obs_hc_unique_pairs = obs_unique_pairs[hc_unique_pairs_range]
    mod_hc_unique_pairs = mod_unique_pairs[hc_unique_pairs_range]
    pairs_hc_unique_values = values_unique_pairs[hc_unique_pairs_range]
    pairs_hc_unique_diff = diff_unique_pairs[hc_unique_pairs_range]


    return obs_hc_unique_pairs, mod_hc_unique_pairs, \
           pairs_hc_unique_values, pairs_hc_unique_diff

def summary_statistics(stat_i, site_i, site_stats_i):

    """
    Calculate the summary statistics for this hour.
    :param stat_i:
    :param site_i:
    :param hr:
    :param stat_data_hr:
    :return: stat_i
    """

    # takes the hours in any order (as they are keys in a dictionary), and places the calculated value in the corresponding
    # idx position, as the arrays have already been made and filled with nans
    for hr, stat_data_hr in site_stats_i.iteritems():

        hridx = int(hr)

        stat_i[site_i]['median'][hridx] = deepcopy(np.nanmedian(stat_data_hr))
        stat_i[site_i]['mean'][hridx] = np.nanmean(stat_data_hr)
        stat_i[site_i]['stdev'][hridx] = np.nanstd(stat_data_hr)

        nanFree = np.array(deepcopy(stat_data_hr))[~np.isnan(np.array(deepcopy(stat_data_hr)))]

        # # if there is data to take a percentile from...
        # if nanFree.shape[0] != 0:
        #
        #     stat_i[site_i]['q25'][hridx] = np.percentile(nanFree, 25)
        #     stat_i[site_i]['q75'][hridx] = np.percentile(nanFree, 75)
        #
        #     stat_i[site_i]['IQR'][hridx] = np.percentile(nanFree, 75) - np.percentile(nanFree, 25)
        #
        # else:
        #     stat_i[site_i]['q25'][hridx] = np.nan
        #     stat_i[site_i]['q75'][hridx] = np.nan
        #
        #     stat_i[site_i]['IQR'][hridx] = np.nan

        # if there is data to take a percentile from...

        stat_i[site_i]['q25'][hridx] = np.percentile(nanFree, 25)
        stat_i[site_i]['q75'][hridx] = np.percentile(nanFree, 75)

        stat_i[site_i]['IQR'][hridx] = np.percentile(nanFree, 75) - np.percentile(nanFree, 25)


    return stat_i

def plot_corr(stat, savedir, site_bsc_colours, model_type, extra=''):

    """
    Plot the median and IQR of the correlation
    :return:
    """

    fig = plt.figure(figsize=(6, 3.5))
    ax = plt.subplot2grid((1, 1), (0, 0))

    for site, site_summary in stat.iteritems():

        # median
        ax.plot(site_summary['hrs'], site_summary['median'],
                label=site, linewidth=1, ls='--', color=site_bsc_colours[site])

        # shade IQR (25th - 75th percentile)
        ax.fill_between(site_summary['hrs'], site_summary['q25'], site_summary['q75'],
                        alpha=0.2, facecolor=site_bsc_colours[site])

    # # prettify
    # # fig.suptitle(data['time'][0].strftime("%Y%m%d") + '-' + data['time'][-1].strftime("%Y%m%d"), fontsize=12)
    ax.set_xlim([0.0, 23.0])
    # ax.set_ylim([-0.5, 1])
    ax.set_xlabel('Time [HH]')
    ax.set_ylabel(r'$Spearman  \/\/correlation$')
    ax.legend(loc='best', fontsize=8)
    # fig.suptitle('median and IQR')
    plt.tight_layout()

    plt.savefig(savedir +'correlations/' +
                      model_type + '_SpearCorrTs_' + '4clearDaysSample_med_IQR_' +
                      extra + '.png')  # filename


    return fig

def plot_rmse(stat, savedir, site_bsc_colours, model_type, extra=''):

        """
        Plot the median and IQR of the rmse
        :return:
        """

        fig = plt.figure(figsize=(6, 3.5))
        ax = plt.subplot2grid((1, 1), (0, 0))

        for site, site_summary in stat.iteritems():
            # median
            ax.plot(site_summary['hrs'], site_summary['median'],
                    label=site, linewidth=1, ls='--', color=site_bsc_colours[site])

            # shade IQR (25th - 75th percentile)
            ax.fill_between(site_summary['hrs'], site_summary['q25'], site_summary['q75'],
                            alpha=0.2, facecolor=site_bsc_colours[site])

        # # prettify
        # # fig.suptitle(data['time'][0].strftime("%Y%m%d") + '-' + data['time'][-1].strftime("%Y%m%d"), fontsize=12)
        ax.set_xlim([0.0, 23.0])
        # ax.set_ylim([-0.5, 1])
        ax.set_xlabel('Hour')
        ax.set_ylabel('RMSE')
        ax.legend(loc='best', fontsize=8)
        fig.suptitle('median and IQR')
        plt.tight_layout()
        plt.savefig(savedir + 'rmse/' +
                    model_type + '_rmse_' + 'clearDaysSample_med_IQR' + extra + '.png')  # filename

        return fig

def plot_diff(stat, savedir, site_bsc_colours, model_type, extra=''):
    """
    Plot the median and IQR of the diff
    :return:
    """

    fig = plt.figure(figsize=(6, 3.5))
    ax = plt.subplot2grid((1, 1), (0, 0))

    for site, site_summary in stat.iteritems():
        # median
        ax.semilogy(site_summary['hrs'], site_summary['median'],
                label=site, linewidth=1, ls='--', color=site_bsc_colours[site])

        # shade IQR (25th - 75th percentile)
        ax.fill_between(site_summary['hrs'], site_summary['q25'], site_summary['q75'],
                        alpha=0.2, facecolor=site_bsc_colours[site])

    # # prettify
    # # fig.suptitle(data['time'][0].strftime("%Y%m%d") + '-' + data['time'][-1].strftime("%Y%m%d"), fontsize=12)
    ax.set_xlim([0.0, 23.0])
    ax.set_ylim([0.0, 50.0])
    ax.set_xlabel('Hour')
    # ax.set_ylabel(r'$Difference \/\mathrm{(log_{10}(\beta_m) - log_{10}(\beta_o))}$')
    ax.set_ylabel('obs_x / mod_y')
    ax.legend(loc='best', fontsize=8)
    # fig.suptitle('median and IQR')
    plt.tight_layout()
    plt.savefig(savedir + 'diff/' +
                model_type + '_difference_normal_' + 'clearDaysSample_med_IQR' + extra + '.png')  # filename

    return fig

def plot_MAE(stat, savedir, site_bsc_colours, model_type, extra=''):

    """
    Plot the median and IQR of the MAE
    :return:
    """

    fig = plt.figure(figsize=(6, 3.5))
    ax = plt.subplot2grid((1, 1), (0, 0))

    for site, site_summary in MAE.iteritems():

        # median
        ax.plot(site_summary['hrs'], site_summary['median'],
                label=site, linewidth=1, ls='--', color=site_bsc_colours[site])

        # shade IQR (25th - 75th percentile)
        ax.fill_between(site_summary['hrs'], site_summary['q25'], site_summary['q75'],
                        alpha=0.2, facecolor=site_bsc_colours[site])

    # # prettify
    # # fig.suptitle(data['time'][0].strftime("%Y%m%d") + '-' + data['time'][-1].strftime("%Y%m%d"), fontsize=12)
    ax.set_xlim([0.0, 23.0])
    # ax.set_ylim([-0.5, 1])
    ax.set_xlabel('Time [HH]')
    ax.set_ylabel(r'$Mean \/\/Absolute \/\/Error\/\/\/ \mathrm{(|\beta_m - \beta_o|)}$')
    ax.legend(loc='best', fontsize=8)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
    # fig.suptitle('median and IQR')
    plt.tight_layout()

    plt.savefig(savedir +'mae/' +
                      model_type + '_hourly_composite_' + 'clearDaysSample_med_IQR_' +
                      extra + '.png')  # filename


    return fig

def plot_AE_med(stat, savedir, site_bsc_colours, model_type, extra=''):

    """
    Plot the median and IQR of the MAE
    :return:
    """

    fig = plt.figure(figsize=(6, 3.5))
    ax = plt.subplot2grid((1, 1), (0, 0))

    for site, site_summary in AE.iteritems():

        # median
        ax.plot(site_summary['hrs'], site_summary['median'],
                label=site, linewidth=1, ls='--', color=site_bsc_colours[site])

        # shade IQR (25th - 75th percentile)
        ax.fill_between(site_summary['hrs'], site_summary['q25'], site_summary['q75'],
                        alpha=0.2, facecolor=site_bsc_colours[site])

    # get number of samples in each hour and position where to plot them (the x values!)
    n = [len(statistics['MR']['AE'][str(i)]) for i in range(24)]
    pos = range(24)

    # add sample size at the top of plot for each box and whiskers
    # pos_t = np.arange(numBoxes) + 1
    upperLabels = [str(np.round(n_i, 2)) for n_i in n]
    weights = ['bold', 'semibold']
    y_lim_max = ax.get_ylim()[1]
    for tick in range(len(pos)):
        k = tick % 2
        ax.text(pos[tick], y_lim_max + (y_lim_max * 0.02), upperLabels[tick],
                 horizontalalignment='center', size='x-small')
    ax2 = ax.twiny()
    ax2.set_xlim([0.0, 23.0])
    ax2.xaxis.set_ticklabels([])

    # # prettify
    # # fig.suptitle(data['time'][0].strftime("%Y%m%d") + '-' + data['time'][-1].strftime("%Y%m%d"), fontsize=12)
    ax.set_xlim(ax2.get_xlim())
    ax.set_xlabel('Time [HH]')
    ax.set_ylabel(r'$Absolute \/\/Error\/\/\/ \mathrm{(|\beta_m - \beta_o|)}$')
    ax.legend(loc='best', fontsize=8)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
    # fig.suptitle('median and IQR')
    plt.tight_layout()

    plt.savefig(savedir +'ae/' +
                      model_type + '_hourly_composite_' + '4clearDaysSample_med_IQR_' +
                      extra + '_new.png')  # filename


    return fig

def plot_AE_mean(stat, savedir, site_bsc_colours, model_type, extra=''):

    """
    Plot the median and IQR of the MAE
    :return:
    """

    fig = plt.figure(figsize=(6, 3.5))
    ax = plt.subplot2grid((1, 1), (0, 0))

    for site, site_summary in AE.iteritems():

        # median
        ax.plot(site_summary['hrs'], site_summary['mean'],
                label=site, linewidth=1, ls='--', color=site_bsc_colours[site])

        # shade IQR (25th - 75th percentile)
        ax.fill_between(site_summary['hrs'],
                        site_summary['mean'] - site_summary['stdev'],
                        site_summary['mean'] + site_summary['stdev'],
                        alpha=0.2, facecolor=site_bsc_colours[site])

    # # prettify
    # # fig.suptitle(data['time'][0].strftime("%Y%m%d") + '-' + data['time'][-1].strftime("%Y%m%d"), fontsize=12)
    ax.set_xlim([0.0, 23.0])
    # ax.set_ylim([-0.5, 1])
    ax.set_xlabel('Time [HH]')
    ax.set_ylabel(r'$Absolute \/\/Error\/\/\/ \mathrm{(|\beta_m - \beta_o|)}$')
    ax.legend(loc='best', fontsize=8)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
    # fig.suptitle('median and IQR')
    plt.tight_layout()

    plt.savefig(savedir +'ae/' +
                      model_type + '_hourly_composite_' + '4clearDaysSample_mean_1stdev_' +
                      extra + '.png')  # filename

    return fig

if __name__ == '__main__':

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

    # instruments and other settings
    # site_bsc = FOcon.site_bsc # turn on for correlations with all ceilometers
    site_rh = FOcon.site_rh
    site_aer = FOcon.site_aer
    site_bsc_colours = FOcon.site_bsc_colours

    site_bsc = {'CL31-A_KSS45W': 64.3 - 31.4,
     'CL31-B_RGS': 8.700000000000003,
     'CL31-C_MR': 4.5,
     'CL31-D_NK': 3.8000000000000007}

    # # # true list (5 Feb 2015 - 31 Dec 2016)
    # daystrList = ['20150414', '20150415', '20150421', '20150611', '20160504', '20160823', '20160911', '20161125',
    #               '20161129', '20161130', '20161204']

    # # true list + 19/01/16 case
    # daystrList = ['20150414', '20150415', '20150421', '20150611', '20160119', '20160504', '20160823', '20160911',
    #               '20161125', '20161129', '20161130', '20161204']

    # main list (5 Feb 2015 - 31 Dec 2016) without 11th June (very high region of RH, near surface in the morning)
    #daystrList = ['20150414', '20150415', '20150421', '20160504', '20160823', '20160911', '20161125',
    #              '20161129', '20161130', '20161204']

    # # days when KSS45W, RGS, MR and NK all had data
    daystrList = ['20150414', '20150415', '20150421', '20160119']
    # # daystrList = ['20160119']

    days_iterate = dateList_to_datetime(daystrList)

    # site for MLH data
    site = 'MR'
    ceil_id = 'CL31-C'
    ceil_id_full = ceil_id + '_' + site

    # statistics to run
    stats_corr = True
    stats_diff = True # Currently calib is turned off!
    stats_RMSE = True # Currently calib is turned off!

    # reduce MLH so it is safely within the BL across ALL ceilometers and not just MR? [%]
    reduce_MLH = 10.0

    # define statistics dictionary
    statistics = {}
    corr = {}
    diff = {}
    MAE = {}
    AE = {}
    rmse = {}
    sampleSize = {}

    # ==============================================================================
    # Read data
    # ==============================================================================

    # Read Ceilometer metadata

    # ceilometer list to use
    ceilsitefile = 'CeilsCSVclearFO.csv'
    # ceilsitefile = 'CeilsCSVcalibrated.csv'
    # ceilsitefile = 'UKV_sitesToExtract.csv'
    ceil_metadata = FO.read_ceil_metadata(datadir, ceilsitefile)

    # just get MLH for MR
    ceil_data_i = {ceil_id_full: site_bsc[ceil_id_full]}

    for day in days_iterate:

        print 'day = ' + day.strftime('%Y-%m-%d')

        # Read UKV forecast and automatically run the FO

        # extract MURK aerosol and calculate RH for each of the sites in the ceil metadata
        # reads all london model data, extracts site data, stores in single dictionary
        mod_data = FO.mod_site_extract_calc(day, ceil_metadata, modDatadir, model_type, res, 910, version=0.2)

        # Read ceilometer backscatter

        # will only read in data is the site is there!
        # ToDo Remove the time sampling part and put it into its own function further down.
        # bsc_obs = FO.read_ceil_obs(day, site_bsc, ceilDatadir, mod_data, calib=False)
        bsc_obs = FO.read_all_ceil_BSC(day, site_bsc, ceilDatadir, timeMatch=mod_data, calib=True)

        # mlh_obs = FO.read_all_ceil_obs(day, site_bsc, ceilDatadir, fType= 'MLH', timeMatch=mod_data, calib=True)
        # mlh_obs = FO.read_all_ceil_obs(day, site_bsc, ceilDatadir, timeMatch=mod_data, fType='MLH')
        mlh_obs = FO.read_all_ceil_obs(day, ceil_data_i, ceilDatadir, timeMatch=mod_data, fType='MLH')


        # ==============================================================================
        # Process
        # ==============================================================================

        # requires model data to be at ceilometer location!
        for site, bsc_site_obs in bsc_obs.iteritems():

            # also requires MLH data
            # if site in mlh_obs:
            if site in bsc_obs:

                # short site id that matches the model id
                site_id = site.split('_')[-1]
                print '     Processing for site: ' + site_id

                # Get unique height pairs between obs and model
                # each height is only paired up once
                # no cutting of maximium height limit yet... that is done in the t loop

                # minHeight for all sites was chosen to be 73.0, which is just higher than IMU
                # Maximum height for all sites are the MLHs estimated from ceilometer backscatter (Simone's algorithm).
                obs_hc_unique_pairs, mod_hc_unique_pairs, \
                pairs_hc_unique_values, pairs_hc_unique_diff = \
                    nearest_heights(mod_data[site_id]['level_height'], bsc_site_obs['height'], minHeight=73.0)

                # create entry in the dictionary if one does not exist
                statistics = create_stats_entry(site_id, statistics)

                # # KSS45W has MLH whereas MR and NK have MH.
                # # Currently they're different so put a work around in for now.
                # if 'MLH' in mlh_obs[site]:
                #     mlh = mlh_obs[site]['MLH']
                # elif 'MH' in mlh_obs[site]:
                #     mlh = mlh_obs[site]['MH']

                # KSS45W has MLH whereas MR and NK have MH.
                # Currently they're different so put a work around in for now.
                if 'MLH' in mlh_obs[ceil_id_full]:
                    mlh = mlh_obs[ceil_id_full]['MLH']
                elif 'MH' in mlh_obs[ceil_id_full]:
                    mlh = mlh_obs[ceil_id_full]['MH']

                # for each hour possible in the day
                for t, MLH_i in zip(np.arange(0, 24), mlh):

                    # reduce MLHby a percentage based on reduce_MLH defined in the Setup
                    MLH_i -= ((MLH_i/100.0)*reduce_MLH)

                    hr = str(t)

                    # pairs below the MLH height at this time
                    mlh_unique_pairs_range = np.where(pairs_hc_unique_values <= MLH_i)[0]

                    # get unique pairs under the MLH for this timestep
                    # do not include this in nearest_heights() as MLH_i varies with each time step and it is simpler
                    # this way, as the unique pairs, prior to MLH_i trimming are always the same.
                    obs_mlh_unique_pairs = obs_hc_unique_pairs[mlh_unique_pairs_range]
                    mod_mlh_unique_pairs = mod_hc_unique_pairs[mlh_unique_pairs_range]

                    # extract out all unique pairs below the upper height limit
                    # these are time and height matched now
                    obs_x = bsc_site_obs['backscatter'][t, obs_mlh_unique_pairs]
                    mod_y = mod_data[site_id]['backscatter'][t, mod_mlh_unique_pairs]

                    # STATISTICS
                    # ---------------

                    # store time
                    # statistics[site_id]['time'] += [mod_data[site_id]['time'][t]]

                    # Correlations
                    if stats_corr == True:

                        # correlate and store
                        # if number of remaining pairs is too low, set r and p to nan
                        try:
                            r, p = spearmanr(obs_x, mod_y, nan_policy='omit')
                        except:
                            r = np.nan
                            p = np.nan

                        statistics[site_id]['r'][hr] += [r]
                        statistics[site_id]['p'][hr] += [p]

                    if stats_diff == True:

                        statistics[site_id]['diff'][hr] += list(mod_y - obs_x)
                        # statistics[site_id]['diff'][hr] += [np.nanmedian(obs_x / mod_y)]
                        statistics[site_id]['MAE'][hr] += [np.nanmean([np.abs(mod_y - obs_x)])]
                        statistics[site_id]['AE'][hr] += list(np.abs(mod_y - obs_x))
                        # statistics[site_id]['diff'][hr] += [np.nanmean(np.log10(mod_y) - np.log10(obs_x))]

                    if stats_RMSE == True:

                        # statistics[site_id]['RMSE'][hr] += [eu.rmse(np.log10(mod_y), np.log10(obs_x))]
                        statistics[site_id]['RMSE'][hr] += [eu.rmse(mod_y, obs_x)]

    # # check to see sample size
    # for h in statistics[site_id]['AE'].iterkeys():
    #     print h
    #     print 'n =' + str(len(statistics[site_id]['AE'][h]))
    #     print ''

    # gather up statistics...
    # create a mean and standard deviation for each hour for plotting

    print '\n' + 'Gathering statistics...'

    corr = {}
    diff = {}
    rmse = {}
    MAE={}
    AE={}
    sampleSize = {}

    for site_i, site_stats in statistics.iteritems():

        print 'site_i = ' + site_i

        # setup site within the summary statistics and carry out statistics
        if stats_corr == True:
            corr = create_stats_summary_dict(site_i, corr)
            corr = summary_statistics(corr, site_i, site_stats['r'])
        if stats_RMSE == True:
            rmse = create_stats_summary_dict(site_i, rmse)
            rmse = summary_statistics(rmse, site_i, site_stats['RMSE'])
        if stats_diff == True:
            diff = create_stats_summary_dict(site_i, diff)
            diff = summary_statistics(diff, site_i, site_stats['diff'])

            MAE = create_stats_summary_dict(site_i, MAE)
            MAE = summary_statistics(MAE, site_i, site_stats['MAE'])

            AE = create_stats_summary_dict(site_i, AE)
            AE = summary_statistics(AE, site_i, site_stats['AE'])
    # plot!
    extra='73mtoNewMLH_minus'+str(int(reduce_MLH)) +'pct'

    if stats_corr == True:
        fig = plot_corr(corr, savedir, site_bsc_colours, model_type, extra=extra)
    if stats_RMSE == True:
        fig = plot_rmse(rmse, savedir, site_bsc_colours, model_type, extra=extra)
    if stats_diff == True:

        fig = plot_diff(diff, savedir, site_bsc_colours, model_type, extra=extra)

        # plot MAE in corr style
        fig = plot_MAE(MAE, savedir, site_bsc_colours, model_type, extra=extra)

        # plot AE in corr style
        fig = plot_AE_med(AE, savedir, site_bsc_colours, model_type, extra=extra)
        fig = plot_AE_mean(AE, savedir, site_bsc_colours, model_type, extra=extra)


    plt.close('all')



print 'END PROGRAM'
