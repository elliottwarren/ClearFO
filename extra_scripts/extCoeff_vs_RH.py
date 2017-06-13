"""
Plot extinction coefficient vs RH

Created by Elliott 13/06/2017
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

    from mod_obs_stats_plot import unique_pairs

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
    savedir = maindir + 'figures/' + model_type + '/sanityCheck/'

    # data
    ceilDatadir = datadir + 'L1/'
    modDatadir = datadir + model_type + '/'
    rhDatadir = datadir + 'L1/'
    aerDatadir = datadir + 'LAQN/'

    # site for model to extract for
    site = 'MR'
    ceil_id = 'CL31-C'
    ceil_id_full = ceil_id + '_' + site

    #daystrList = ['20150414', '20150415', '20150421', '20150611', '20160504', '20160823', '20160911', '20161125',
    #              '20161129', '20161130', '20161204']
    #days_iterate = dateList_to_datetime(daystrList)

    corr_max_height = 2000.0

    # forecast data start time
    Z='21'

    days_iterate = [dt.datetime(2016,01,19)]# new PM10 case study day

    # ==============================================================================
    # Read data
    # ==============================================================================

    # Read Ceilometer metadata

    # ceilometer list to use
    ceilsitefile = 'CeilsCSVfull.csv'
    ceil_metadata = FO.read_ceil_metadata(datadir, ceilsitefile)

    # extract out current site only
    ceil_data_i = {site: ceil_metadata[site]}

    for day in days_iterate:

        print 'day = ' + day.strftime('%Y-%m-%d')

        # extract MURK aerosol and calculate RH for each of the sites in the ceil metadata
        # reads all london model data, extracts site data, stores in single dictionary
        mod_data_old = FO.mod_site_extract_calc(day, ceil_data_i, modDatadir, model_type, res, 910,
                                            Z=Z, allvars=True, version=0.1)

        mod_data_new = FO.mod_site_extract_calc(day, ceil_data_i, modDatadir, model_type, res, 910,
                                            Z=Z, allvars=True, version=0.2)


        # store ext_coeff and RH for lowest ~2000 m
        if day == days_iterate[0]:

            RH_all_old = mod_data_old['MR']['RH'][:, 0:24]
            RH_all_new = mod_data_new['MR']['RH'][:, 0:24]
            ext_old = mod_data_old['MR']['alpha_a'][:, 0:24]
            ext_new = mod_data_new['MR']['alpha_a'][:, 0:24]

        else:

            RH_all_old = np.append(RH_all_old, mod_data_old['MR']['RH'][:, 0:24])
            RH_all_new = np.append(RH_all_new, mod_data_new['MR']['RH'][:, 0:24])
            ext_old = np.append(ext_old, mod_data_old['MR']['alpha_a'][:, 0:24])
            ext_new = np.append(ext_new, mod_data_new['MR']['alpha_a'][:, 0:24])

        # plot alpha_a vs RH (day)
        fig, ax = plt.subplots(1, 1, figsize=(6, 6))
        s=3

        plt.scatter(mod_data_old[site]['RH']*100.0, mod_data_old[site]['alpha_a']*1000.0, label='r_m (swollen)', s=s)
        plt.scatter(mod_data_new[site]['RH']*100.0, mod_data_new[site]['alpha_a']*1000.0, label='r_md (dry)', s=s)

        ax.set_xlabel('RH [%]', fontsize=10, labelpad=10)
        ax.set_ylabel('extinction coefficient [km-1]', fontsize=10, labelpad=2)

        plt.suptitle('day = '+ str(day.date()) +'; lambda = 910nm')
        ax.set_xlim([0.0, 100.0])
        plt.legend()
        plt.savefig(savedir + model_type + '-' + site + '-' + str(day.date()) + '-extCoeff_vs_RH.png')  # filename
        plt.close(fig)

    # plot alpha_a vs RH (all)
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    s = 5

    plt.scatter(RH_all_old * 100.0, ext_old * 1000.0, label='r_m (swollen)', s=s)
    plt.scatter(RH_all_new * 100.0, ext_new * 1000.0, label='r_md (dry)', s=s)

    ax.set_xlabel('RH [%]', fontsize=10, labelpad=10)
    ax.set_ylabel('extinction coefficient [km-1]', fontsize=10, labelpad=2)

    plt.suptitle('11 sample days; lambda = 910nm')
    ax.set_xlim([0.0, 100.0])
    ax.set_ylim([0.0, 1.2])
    plt.legend()
    plt.savefig(savedir + model_type + '-' + site + '-allDays-extCoeff_vs_RH.png')  # filename
    plt.close(fig)





















    print 'END PROGRAM'


if __name__ == '__main__':
    main()