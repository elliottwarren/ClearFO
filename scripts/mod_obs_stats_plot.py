"""
Script to do all the stats to the FO output. Correlations first...

Created by Elliott Thur 27th Oct 2016
"""

import matplotlib.pyplot as plt
from matplotlib.dates import date2num
from matplotlib.dates import DateFormatter

import numpy as np
import datetime as dt
from scipy.stats import spearmanr

import ceilUtils as ceil
import ellUtils as eu
from forward_operator import FOUtils as FO
from forward_operator import FOconstants as FOcon

def stats_dic(site_bsc, mbe_limit_max, mbe_limit_step):

    """
    Define the almighty statistics dictionary.

    :param site_bsc:
    :param mbe_limit_max:
    :param mbe_limit_step:
    :return: statistics (dict)

    statistics[site]['r'] = [...]
    statistics[site]['MBE'] = {'0-500': ..., '500-1000': ...}
    statistics[site]['time'] = [...]
    """

    # calculate the height groups
    mbe_height_limits = np.arange(0, mbe_limit_max + mbe_limit_step, mbe_limit_step)
    # list of strings to match the ranges
    height_groups_order = np.array([str(i) + '-' + str(i + mbe_limit_step) for i in mbe_height_limits[:-1]])

    # define array to hold statistics
    statistics = {}
    # statistics[site]['r'] = [...]
    # statistics[site]['MBE'] = {'0-500': ..., '500-1000': ...}
    # statistics[site]['time'] = [...]

    # define site based lists to store the correlation results in
    for site in site_bsc.iterkeys():

        site_id = site.split('_')[-1]
        statistics[site_id] = {'r': [], 'p': [],
                               'MBE': {},
                               'time': []}

        for hg in height_groups_order:
            statistics[site_id]['MBE'][hg] = []

    return statistics

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

def plot_correlations(savedir, model_type, statistics, corr_max_height):

    """
    plot the correlation statistics (\beta_m, site vs \beta_o, site) and save.

    :return: fig
    """

    fig = plt.figure(figsize=(8, 5))
    ax = plt.subplot2grid((1, 1), (0, 0))

    for site, data in statistics.iteritems():

        plt.plot_date(data['time'], data['r'], label=site, linewidth=1, fmt='-')

    # plot reference line to show where profile lies
    ax.plot_date([statistics['RGS']['time'][24], statistics['RGS']['time'][24]], [-1, 1], color='black', ls='--', fmt='--')
    ax.plot_date([statistics['RGS']['time'][49], statistics['RGS']['time'][49]], [-1, 1], color='black', ls='--', fmt='--')

    # prettify
    # fig.suptitle(data['time'][0].strftime("%Y%m%d") + '-' + data['time'][-1].strftime("%Y%m%d"), fontsize=12)
    ax.set_xlim([data['time'][0], data['time'][-1]])
    ax.set_ylim([-0.5, 1])
    ax.set_xlabel('Time [DD/ HH:MM]')
    ax.set_ylabel(r'$Spearman \/\/\rho \/\/correlation$')
    ax.xaxis.set_major_formatter(DateFormatter('%d/ %H:%M'))
    ax.legend(loc='best', fontsize=8)

    plt.savefig(savedir +'correlations/' +
                      model_type + '_SpearCorrTs_' +
                      data['time'][0].strftime("%Y%m%d") + '-' + data['time'][-1].strftime("%Y%m%d")+'_'+
                      str(corr_max_height) + 'm.png')  # filename

    return fig

def plot_mbe(savedir, model_type, statistics, height_groups_order, site_bsc_colours):

    """
    Plot the mean bias error (MBE) statistics

    :param savedir:
    :param model_type:
    :param statistics:
    :param height_groups_order:
    :param site_bsc_colours:
    :return: fig
    """


    fig, ax = plt.subplots(4, 1, figsize=(10, 5))

    # plot mbe data
    for p, hg in zip(ax, height_groups_order):

        for site, data in statistics.iteritems():

            p.plot_date(date2num(data['time']),
                            data['MBE'][hg],
                            label=site, linewidth=1, color=site_bsc_colours[site], fmt='-')

        # prettify - plot specific and done once the data is plotted
        p.set_xlim([data['time'][0], data['time'][-1]])
        eu.add_at(p, hg + ' m', loc=4)
        p.set_ylim([-1.8, 1.8])
        p.yaxis.set_ticks(np.arange(1.5, -2.5, -1))

    # reference lines
    for p in np.arange(len(ax)):

        time = statistics['RGS']['time']

        ax[p].plot_date(date2num([statistics['RGS']['time'][0], statistics['RGS']['time'][-1]]), [0, 0],
                        color='black', ls='--', fmt='-')

        # turn off the x axis labels so only the bottom one plots them
        if p < len(ax):
            ax[p].set_xticklabels([])





    # prettify - overall

    ax[0].legend(fontsize=8, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)

    # setup a figure wide axis for labels
    ax0 = eu.fig_majorAxis(fig)
    ax0.set_ylabel('Difference (log10(model) - log10(obs))')
    ax0.set_ylabel(r'$Difference \/\mathrm{(log_{10}(\beta_m) - log_{10}(\beta_o))}$')
    ax0.set_xlabel('Time [HH:MM]')
    ax[-1].xaxis.set_major_formatter(DateFormatter('%H:%M'))
    plt.setp(ax[-1].get_xticklabels())

    plt.tight_layout(h_pad=0.2)  # need to modify the padding to improve the plot
    fig.subplots_adjust(right=0.8)

    plt.savefig(savedir + 'mbe/' +
                model_type + '_meanBiasError_sameSites.png')  # filename

    return fig

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
    savedir = maindir + 'figures/' + model_type + '/'

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
    dayStart = dt.datetime(2016, 05, 04)
    dayEnd = dt.datetime(2016, 05, 06)

    # statistics to run
    stats_corr = True
    stats_mbe = False

    # mbe ranges
    mbe_limit_step = 500
    mbe_limit_max = 2000

    # correlation max height
    corr_max_height = 2000

    # calculate the height groups and matching strings
    mbe_height_limits = np.arange(0, mbe_limit_max + mbe_limit_step, mbe_limit_step)
    height_groups_order = np.array([str(i) + '-' + str(i + mbe_limit_step) for i in mbe_height_limits[:-1]])

    # set up statistics dictionary
    statistics = stats_dic(site_bsc, mbe_limit_max, mbe_limit_step)

    # ==============================================================================
    # Read data
    # ==============================================================================

    # 1. Read Ceilometer metadata
    # ----------------------------
    ceil_metadata = FO.read_ceil_metadata(datadir)

    # datetime range to iterate over
    days_iterate = eu.date_range(dayStart, dayEnd, 1, 'days')

    for day in days_iterate:


        # 1. Read UKV forecast in
        # -----------------------

        # extract MURK aerosol and calculate RH for each of the sites in the ceil metadata
        # (can be different locations to sites_bsc)
        # reads all london model data, extracts site data, stores in single dictionary
        mod_data = FO.mod_site_extract_calc(day, ceil_metadata, modDatadir, model_type, res, 910)


        # 2. Read ceilometer backscatter
        # --------------------------------

        bsc_obs = FO.read_ceil_obs(day, site_bsc, ceilDatadir, mod_data)

        # ==============================================================================
        # Process modelled data
        # ==============================================================================


        # requires model data to be at ceilometer location!
        for site, bsc_site_obs in bsc_obs.iteritems():

            # short site id that matches the model id
            site_id = site.split('_')[-1]

            # get the nearest ceilometer height gate to each model level
            # obs_idx = ALL nearest gate idx
            # mod_idx = idx of the model height that each obs_idx are paired to
            a = np.array([eu.nearest(bsc_site_obs['height'], i)for i in mod_data[site_id]['level_height']])
            values = a[:, 0]
            obs_idx = np.array(a[:, 1], dtype=int)
            diff = a[:, 2]
            mod_idx = np.arange(len(mod_data[site_id]['level_height'])) # mod_idx should be paired with obs_idx spots.

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

            # Remove pairs where obs is above 2000 m.
            # hc = height cut
            hc_unique_pairs_range = np.where(values_unique_pairs <= corr_max_height)[0]

            # trim off unique pairs that are above the maximum height
            obs_hc_unique_pairs = obs_unique_pairs[hc_unique_pairs_range]
            mod_hc_unique_pairs = mod_unique_pairs[hc_unique_pairs_range]
            pairs_hc_unique_values = values_unique_pairs[hc_unique_pairs_range]
            pairs_hc_unique_diff = diff_unique_pairs[hc_unique_pairs_range]

            # statistics
            for t in np.arange(len(mod_data[site_id]['time'])):

                # extract out all unique pairs below the upper height limit
                obs_x = bsc_site_obs['backscatter'][t, obs_hc_unique_pairs]
                mod_y = mod_data[site_id]['backscatter'][t, mod_hc_unique_pairs]

                # store time
                statistics[site_id]['time'] += [mod_data[site_id]['time'][t]]

                if stats_corr == True:

                    # correlate and store
                    # if number of remaining pairs is too low, set r and p to nan
                    try:
                        r, p = spearmanr(np.log10(obs_x), np.log10(mod_y), nan_policy='omit')
                    except:
                        r = np.nan
                        p = np.nan

                    statistics[site_id]['r'] += [r]
                    statistics[site_id]['p'] += [p]

                if stats_mbe == True:

                    # arrays of idx for each group
                    # list type to retain plotting order
                    height_groups_idx = [np.where((pairs_hc_unique_values >= i) &
                                            (pairs_hc_unique_values < i + mbe_limit_step))[0]
                                    for i in mbe_height_limits[:-1]]

                    # calculate mbe for each height group
                    for i in np.arange(len(height_groups_idx)):

                        # further divide the data based on the current height group (hg)
                        obs_x_hg = obs_x[height_groups_idx[i]]
                        mod_y_hg = mod_y[height_groups_idx[i]]

                        statistics[site_id]['MBE'][height_groups_order[i]] += \
                        [np.nanmean(np.log10(mod_y_hg) - np.log10(obs_x_hg))]

    # ==============================================================================
    # Plotting
    # ==============================================================================

    # After all day's stats are done

    if stats_corr == True:
        fig = plot_correlations(savedir, model_type, statistics, corr_max_height)

    if stats_mbe == True:

        fig = plot_mbe(savedir, model_type, statistics, height_groups_order, site_bsc_colours)

    print 'END PROGRAM'

    plt.close('all')

    return

if __name__ == '__main__':
    main()