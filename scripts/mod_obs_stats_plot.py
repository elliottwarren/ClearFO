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
    site_bsc = FO.site_bsc
    site_rh = FO.site_rh
    site_aer = FO.site_aer
    site_bsc_colours = FO.site_bsc_colours

    # day start and end
    dayStart = dt.datetime(2016, 05, 04)
    dayEnd = dt.datetime(2016, 05, 06)


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
        mod_data = FO.mod_site_extract_calc(day, ceil_metadata, modDatadir, model_type, res)


        # 2. Read ceilometer backscatter
        # --------------------------------

        bsc_obs = FO.read_ceil_obs(day, site_bsc, ceilDatadir, mod_data)

        # ==============================================================================
        # Process modelled data
        # ==============================================================================

        if day == days_iterate[0]:

            # define array to hold statistics
            statistics = {}
            # statistics['r'] = [...]
            # statistics['MBE'] = [...]
            # statistics['time'] = [...]

            # define site based lists to store the correlation results in
            for site in bsc_obs.iterkeys():

                site_id = site.split('_')[-1]
                statistics[site_id] = {'r': [], 'p': [], 'time': []}

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

            # trim off the ends of obs_idx, as UKV and obs z0 and zmax are different, leading to the same gate matching multiple ukvs

            # assumes no duplicates in the middle of the arrays, just at the end
            # Hence, data within the array is complete.

            # At this point, variables are like:
            # obs_idx = [0, 0, 0, 1, 3, 5, .... 769, 769, 769]
            # mod_idx = [0, 1, 2, 3, 4, 4, .... 67,  68,  69 ]
            # need to remove duplicate pairs where obs are used multiple times and keeping the ones
            # where the differences is smallest -> It will remove some mod_idx locations too e.g. mod_idx = 0
            # will need to be removed if the diff with obs_idx = 0 is larger compared to mod_idx = 1.

            # find start idx to remove duplicate pairs
            duplicates = np.where(obs_idx == obs_idx[0])[0] # find duplicates
            if len(duplicates) > 1:
                lowest_diff = np.argmin(abs(diff[duplicates]))  # find which has smallest difference
                pairs_idx_start = duplicates[lowest_diff]  # set start position for pairing at this point
            else:
                pairs_idx_start = 0

            # find end idx to remove duplicate pairs
            duplicates = np.where(obs_idx == obs_idx[-1])[0] # find duplicates
            if len(duplicates) > 1:
                lowest_diff = np.argmin(abs(diff[duplicates])) # find which has smallest difference
                pairs_idx_end = duplicates[lowest_diff] # set start position for pairing at this point
            else:
                pairs_idx_end = len(obs_idx)

            # create range in order to extract the unique pairs
            unique_pairs_range = np.arange(pairs_idx_start, pairs_idx_end + 1)

            # trim off obs_idx positions that would be used multiple times in correlation
            # keep the values and differences for those pairs
            # use these to plot correlations for all possible pairs, regardless of height
            obs_unique_pairs = obs_idx[unique_pairs_range]
            mod_unique_pairs = mod_idx[unique_pairs_range]
            values_unique_pairs = values[unique_pairs_range]
            diff_unique_pairs = diff[unique_pairs_range]

            # ~~~~~~~~~~~~~~~~~~~~ #

            # Remove pairs where obs is above 2000 m.
            # hc = height cut
            hc_unique_pairs_range = np.where(values_unique_pairs < 2000)[0]

            # trim off unique pairs that are above the maximum height
            obs_hc_unique_pairs = obs_unique_pairs[hc_unique_pairs_range]
            mod_hc_unique_pairs = mod_unique_pairs[hc_unique_pairs_range]
            pairs_hc_unique_values = values_unique_pairs[hc_unique_pairs_range]
            pairs_hc_unique_diff = diff_unique_pairs[hc_unique_pairs_range]

            # correlate \beta_m with \beta_o for the same location
            for t in np.arange(len(mod_data[site_id]['time'])):

                # extract out the pairs to correlate
                obs_x = bsc_site_obs['backscatter'][t, obs_hc_unique_pairs]
                mod_y = mod_data[site_id]['backscatter'][t, mod_hc_unique_pairs]

                # correlate and store
                r, p = spearmanr(obs_x, mod_y)
                statistics[site_id]['r'] += [r]
                statistics[site_id]['p'] += [p]
                statistics[site_id]['time'] += [mod_data[site_id]['time'][t]]

    # ==============================================================================
    # Plotting
    # ==============================================================================

    # After all day's stats are done

    fig = plt.figure(figsize=(7, 5))
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
    ax.set_xlabel('Time [HH:MM]')
    ax.set_ylabel(r'$Spearman \/\/\rho \/\/correlation$')
    ax.xaxis.set_major_formatter(DateFormatter('%d/ %H:%M'))
    ax.legend(loc='best', fontsize=8)


    plt.savefig(savedir +'correlations/' +
                      model_type + '_SpearCorrTs_' +
                      data['time'][0].strftime("%Y%m%d") + '-' + data['time'][-1].strftime("%Y%m%d")+ '.png')  # filename


    print 'END PROGRAM'

    plt.close('all')

    return

if __name__ == '__main__':
    main()