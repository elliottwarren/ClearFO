"""
Hopefully the start of something good.
Make daily plots with hourly statistics with standard deviations and everything

Created by Elliott Thurs 27th April 2017
"""

import matplotlib.pyplot as plt
from matplotlib.dates import date2num
from matplotlib.dates import DateFormatter

import numpy as np
import datetime as dt
from scipy.stats import spearmanr

import ellUtils as eu
from mod_obs_stats_plot import unique_pairs
from forward_operator import FOUtils as FO
from forward_operator import FOconstants as FOcon

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
                               'MBE': {}}

        for hr in hrs:
            statistics[site_id]['r'][str(hr)] = []
            statistics[site_id]['p'][str(hr)] = []
            statistics[site_id]['diff'][str(hr)] = []
            statistics[site_id]['RMSE'][str(hr)] = []

    return statistics

def create_stats_summary_dict_mean(site_id, summary={}):

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
        summary[site_id] = {'mean': nanArray, 'mean_plus_stdev': nanArray,
                               'mean_minus_stdev': nanArray,
                               'stdev': nanArray,
                               'hrs': hrs}

    return summary

def create_stats_summary_dict_med(site_id, summary={}):

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
        summary[site_id] = {'median': nanArray, 'q25': nanArray,
                               'q75': nanArray,
                               'IQR': nanArray,
                               'hrs': hrs}

    return summary

def dateList_to_datetime(dayList):

    """ Convert list of string dates into datetimes """

    datetimeDays = []

    for d in dayList:

        datetimeDays += [dt.datetime(int(d[0:4]), int(d[4:6]), int(d[6:8]))]

    return datetimeDays

def nearest_heights(mod_height, obs_height, corr_max_height):

    """
    Get the nearest ceilometer height gate to each model level

    :param mod_height:
    :param obs_height:
    :param corr_max_height:
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

    # Remove pairs where obs is above the max allowed height.
    # hc = height cut
    hc_unique_pairs_range = np.where(values_unique_pairs <= corr_max_height)[0]

    # trim off unique pairs that are above the maximum height
    obs_hc_unique_pairs = obs_unique_pairs[hc_unique_pairs_range]
    mod_hc_unique_pairs = mod_unique_pairs[hc_unique_pairs_range]
    pairs_hc_unique_values = values_unique_pairs[hc_unique_pairs_range]
    pairs_hc_unique_diff = diff_unique_pairs[hc_unique_pairs_range]


    return obs_hc_unique_pairs, mod_hc_unique_pairs, \
           pairs_hc_unique_values, pairs_hc_unique_diff

def summary_statistics_mean(stat_i,site_i, hr, stat_data_hr):

    """
    Calculate the summary statistics for this hour.
    :param stat_i:
    :param site_i:
    :param hr:
    :param stat_data_hr:
    :return: stat_i
    """

    hridx = int(hr)

    stat_i[site_i]['mean'][hridx] = np.nanmean(stat_data_hr)
    stat_i[site_i]['stdev'][hridx] = np.nanstd(stat_data_hr)

    stat_i[site_i]['mean_minus_stdev'][hridx] = stat_i[site_i]['mean'][hridx] - stat_i[site_i]['stdev'][hridx]
    stat_i[site_i]['mean_plus_stdev'][hridx] = stat_i[site_i]['mean'][hridx] + stat_i[site_i]['stdev'][hridx]


    return stat_i

def summary_statistics_med(stat_i,site_i, hr, stat_data_hr):

    """
    Calculate the summary statistics for this hour.
    :param stat_i:
    :param site_i:
    :param hr:
    :param stat_data_hr:
    :return: stat_i
    """

    hridx = int(hr)

    stat_i[site_i]['median'][hridx] = np.nanmedian(stat_data_hr)

    nanFree = np.array(stat_data_hr)[~np.isnan(np.array(stat_data_hr))]

    stat_i[site_i]['q25'][hridx] = np.percentile(nanFree, 25)
    stat_i[site_i]['q75'][hridx] = np.percentile(nanFree, 75)

    stat_i[site_i]['IQR'][hridx] = stat_i[site_i]['q75'] - stat_i[site_i]['q25']


    return stat_i

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
    savedir = maindir + 'figures/' + model_type + '/longterm/'

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

    # day list
    # full list
    #daystrList = ['20150414', '20150415', '20150421', '20150611', '20160504', '20160823', '20160911', '20161125',
    #              '20161129', '20161130', '20161204', '20170120', '20170122', '20170325', '20170408']

    # current list
    daystrList = ['20150414', '20150415', '20150421', '20150611', '20160504', '20160823', '20160911', '20161125',
                  '20161129']

    # daystrList = ['20160504', '20160505']

    days_iterate = dateList_to_datetime(daystrList)

    # statistics to run
    stats_corr = True
    stats_diff = True
    stats_RMSE = True

    # correlation max height
    corr_max_height = 2000

    # define statistics dictionary
    statistics = {}
    corr = {}
    diff = {}
    rmse = {}
    sampleSize = {}


    # ==============================================================================
    # Read data
    # ==============================================================================

    # Read Ceilometer metadata

    # ceilometer list to use
    ceilsitefile = 'CeilsCSVfull.csv'
    ceil_metadata = FO.read_ceil_metadata(datadir, ceilsitefile)

    for day in days_iterate:

        print 'day = ' + day.strftime('%Y-%m-%d')

        # Read UKV forecast and automatically run the FO

        # extract MURK aerosol and calculate RH for each of the sites in the ceil metadata
        # reads all london model data, extracts site data, stores in single dictionary
        mod_data = FO.mod_site_extract_calc(day, ceil_metadata, modDatadir, model_type, res, 910)

        # Read ceilometer backscatter

        # will only read in data is the site is there!
        # ToDo Remove the time sampling part and put it into its own function further down.
        bsc_obs = FO.read_ceil_obs(day, site_bsc, ceilDatadir, mod_data)

        # ==============================================================================
        # Process
        # ==============================================================================

        # requires model data to be at ceilometer location!
        for site, bsc_site_obs in bsc_obs.iteritems():

            # short site id that matches the model id
            site_id = site.split('_')[-1]
            print '     Processing for site: ' + site_id

            # Get unique height pairs between obs and model
            # each height is only paired up once
            # heights above a maximum limit are cut (define by corr_max_height)

            obs_hc_unique_pairs, mod_hc_unique_pairs, \
            pairs_hc_unique_values, pairs_hc_unique_diff = \
                nearest_heights(mod_data[site_id]['level_height'], bsc_site_obs['height'], corr_max_height)

            # create entry in the dictionary if one does not exist
            statistics = create_stats_entry(site_id, statistics)
            #stat_mean  = create_stats_entry(site_id, stat_mean)
            #stat_stdev = create_stats_entry(site_id, stat_stdev)
            #stat_n = create_stats_entry(site_id, stat_n)

            # for each hour possible in the day
            for t in np.arange(0, 24):

                hr = str(t)

                # extract out all unique pairs below the upper height limit
                # these are time and height matched now
                obs_x = bsc_site_obs['backscatter'][t, obs_hc_unique_pairs]
                mod_y = mod_data[site_id]['backscatter'][t, mod_hc_unique_pairs]

                # STATISTICS
                # ---------------

                # store time
                # statistics[site_id]['time'] += [mod_data[site_id]['time'][t]]

                # Correlations
                if stats_corr == True:

                    # correlate and store
                    # if number of remaining pairs is too low, set r and p to nan
                    try:
                        r, p = spearmanr(np.log10(obs_x), np.log10(mod_y), nan_policy='omit')
                    except:
                        r = np.nan
                        p = np.nan

                    statistics[site_id]['r'][hr] += [r]
                    statistics[site_id]['p'][hr] += [p]

                if stats_diff == True:

                    statistics[site_id]['diff'][hr] += [np.nanmedian(np.log10(mod_y) - np.log10(obs_x))]
                    # statistics[site_id]['diff'][hr] += [np.nanmean(np.log10(mod_y) - np.log10(obs_x))]

                if stats_RMSE == True:

                    statistics[site_id]['RMSE'][hr] += [eu.rmse(np.log10(mod_y), np.log10(obs_x))]


    # gather up statistics...
    # create a mean and standard deviation for each hour for plotting

    print '\n' + 'Gathering statistics...'

    for site_i, site_stats in statistics.iteritems():

        # setup site within the summary statistics
        corr = create_stats_summary_dict_med(site_i, corr)
        rmse = create_stats_summary_dict_med(site_i, rmse)
        diff = create_stats_summary_dict_med(site_i, diff)

        for stat, stat_data_all_hrs in site_stats.iteritems():

            for hr, stat_data_hr in stat_data_all_hrs.iteritems():

                corr = summary_statistics_med(corr, site_i, hr, stat_data_hr)
                rmse = summary_statistics_med(rmse, site_i, hr, stat_data_hr)
                diff = summary_statistics_med(diff, site_i, hr, stat_data_hr)

                #corr = summary_statistics_mean(corr, site_i, hr, stat_data_hr)
                #rmse = summary_statistics_mean(rmse, site_i, hr, stat_data_hr)
                #diff = summary_statistics_mean(diff, site_i, hr, stat_data_hr)

    # for site_i, site_stats in statistics.iteritems():
    #
    #     # setup site within the summary statistics
    #     corr = create_stats_summary_dict_med(site_i, corr)
    #     rmse = create_stats_summary_dict_med(site_i, rmse)
    #     diff = create_stats_summary_dict_med(site_i, diff)
    #
    #     for stat, stat_data_all_hrs in site_stats.iteritems():
    #
    #         for hr, stat_data_hr in stat_data_all_hrs.iteritems():
    #
    #             corr = summary_statistics_med(corr, site_i, hr, stat_data_hr)
    #             rmse = summary_statistics_med(rmse, site_i, hr, stat_data_hr)
    #             diff = summary_statistics_med(diff, site_i, hr, stat_data_hr)
    #
    #             #corr = summary_statistics_mean(corr, site_i, hr, stat_data_hr)
    #             #rmse = summary_statistics_mean(rmse, site_i, hr, stat_data_hr)
    #             #diff = summary_statistics_mean(diff, site_i, hr, stat_data_hr)


    # plot!

    #"""
    #plot the correlation statistics (\beta_m, site vs \beta_o, site) and save.#
    #
    #:return: fig
    #"""

    fig = plt.figure(figsize=(6, 3.5))
    ax = plt.subplot2grid((1, 1), (0, 0))

    for site, site_summary in corr.iteritems():

        # # mean
        # ax.plot(site_summary['hrs'], site_summary['mean'],
        #         label='mean', linewidth=1, ls='--', color=site_bsc_colours[site])
        #
        # # shade the range of stdev
        # ax.fill_between(site_summary['hrs'], site_summary['mean_minus_stdev'], site_summary['mean_plus_stdev'],
        #                 alpha=0.5, facecolor=site_bsc_colours[site])

        # mean
        ax.plot(site_summary['hrs'], site_summary['median'],
                label='mean', linewidth=1, ls='--', color=site_bsc_colours[site])

        # shade the range of stdev
        ax.fill_between(site_summary['hrs'], site_summary['q25'], site_summary['q75'],
                        alpha=0.2, facecolor=site_bsc_colours[site])

    # # prettify
    # # fig.suptitle(data['time'][0].strftime("%Y%m%d") + '-' + data['time'][-1].strftime("%Y%m%d"), fontsize=12)
    ax.set_xlim(site_summary['hrs'].min(), site_summary['hrs'].max())
    # ax.set_ylim([-0.5, 1])
    ax.set_xlabel('Hour')
    ax.set_ylabel(r'$Spearman \/\/\rho \/\/correlation$')
    # ax.xaxis.set_major_formatter(DateFormatter('%d/ %H:%M'))
    # ax.legend(loc='best', fontsize=8)

    plt.savefig(savedir +'correlations/' +
                      model_type + '_SpearCorrTs_' + 'clearDaysSample_' +
                      str(corr_max_height) + 'm.png')  # filename





















    return

if __name__ == '__main__':
    main()























print 'END PROGRAM'
