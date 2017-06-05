"""
4 panel plot of beta_m, beta_o, m and RH for a day

Created by Elliott 30/05/17
"""

import matplotlib.pyplot as plt
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
    savedir = maindir + 'figures/' + model_type + '/4panel/'

    # data
    ceilDatadir = datadir + 'L1/'
    modDatadir = datadir + model_type + '/'
    rhDatadir = datadir + 'L1/'
    aerDatadir = datadir + 'LAQN/'

    # statistics to run
    pm10_stats = False
    rh_stats = True

    site = 'MR'
    ceil_id = 'CL31-C'
    ceil_id_full = ceil_id + '_' + site


    site_bsc = {ceil_id_full: FOcon.site_bsc[ceil_id_full]}
    # site_bsc = {ceil: FOcon.site_bsc[ceil], 'CL31-E_BSC_NK': 27.0 - 23.2}

    # site_aer = {'PM10_'+site: FOcon.site_aer['PM10_'+site]}

    site_rh = {'WXT_KSSW': 50.3}
    rh_instrument = site_rh.keys()[0]

    site_bsc_colours = FOcon.site_bsc_colours

    #daystrList = ['20150414', '20150415', '20150421', '20150611', '20160504', '20160823', '20160911', '20161125',
    #              '20161129', '20161130', '20161204']
    #days_iterate = dateList_to_datetime(daystrList)

    corr_max_height = 2000.0


    days_iterate = [dt.datetime(2016,01,19)]# new PM10 case study day
    # day = [dt.datetime(2016, 05, 04)] # one of my old case study days

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

        # Read UKV forecast and automatically run the FO

        # extract MURK aerosol and calculate RH for each of the sites in the ceil metadata
        # reads all london model data, extracts site data, stores in single dictionary
        mod_data = FO.mod_site_extract_calc(day, ceil_data_i, modDatadir, model_type, res, 910)

        # Read ceilometer backscatter

        # will only read in data is the site is there!
        # ToDo Remove the time sampling part and put it into its own function further down.
        # bsc_obs = FO.read_ceil_obs(day, site_bsc, ceilDatadir, mod_data, calib=True)
        bsc_obs = FO.read_ceil_obs_all(day, site_bsc, ceilDatadir)

        # sub sampled in time
        bsc_obs_sub = FO.read_ceil_obs(day, site_bsc, ceilDatadir, mod_data)


        # # read in PM10 data and extract data for the current day
        # pm10 = FO.read_pm10_obs(site_aer, aerDatadir, matchModSample=False)
        #
        # # extract the current day out of pm10
        # # .date() from pm10 dates
        # dates = np.array([i.date() for i in pm10['PM10_'+site]['time']])
        # idx = np.where(dates == day.date())
        #
        # # extract
        # pm10['PM10_'+site]['pm_10'] = pm10['PM10_'+site]['pm_10'][idx]
        # pm10['PM10_'+site]['time'] = [pm10['PM10_'+site]['time'][i] for i in idx]

        # read in RH data
        # rh_obs = FO.read_all_rh_obs(day, site_rh, rhDatadir, mod_data)


        # MBE #

        # get unique pairs of ceilometer obs to model, then do a simple MBE
        obs_hc_unique_pairs, mod_hc_unique_pairs, \
        pairs_hc_unique_values, pairs_hc_unique_diff = \
            nearest_heights(mod_data[site]['level_height'], bsc_obs_sub[ceil_id_full]['height'], corr_max_height)

        # extract out all unique pairs below the upper height limit
        # these are time and height matched now
        obs_x = bsc_obs_sub[ceil_id_full]['backscatter'][:, obs_hc_unique_pairs]
        mod_y = mod_data[site]['backscatter'][:, mod_hc_unique_pairs]

        mbe = np.log10(mod_y) - np.log10(obs_x)


        # plot the data
        # 4 panel, beta_o, beta_m, m with pm10 overlay, rh with rh_obs (KSSW) overlay
        fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, figsize=(7, 6.5))

        site_id = site.split('_')[-1]

        # beta_o
        # mesh1 = ax1.pcolormesh(bsc_obs[ceil_id_full]['time'], bsc_obs[ceil_id_full]['height'], np.transpose(bsc_obs[ceil_id_full]['backscatter']),
        #                                   norm=LogNorm(vmin=1e-7, vmax=1e-5), cmap=cm.get_cmap('jet'))

        mesh1 = ax1.pcolormesh(bsc_obs[ceil_id_full]['time'], bsc_obs[ceil_id_full]['height'],
                               np.transpose(bsc_obs[ceil_id_full]['backscatter']),
                                          norm=LogNorm(vmin=1e-7, vmax=1e-5), cmap=cm.get_cmap('jet'))

        # beta_m
        mesh2 = ax2.pcolormesh(mod_data[site_id]['time'], mod_data[site_id]['level_height'],
                               np.transpose(mod_data[site_id]['backscatter']),
                                          norm=LogNorm(vmin=1e-7, vmax=1e-5), cmap=cm.get_cmap('jet'))

        # ax2.plot([mod_data[site_id]['time'][0],mod_data[site_id]['time'][-1]], [111.67, 111.67], ls='--', color='black')

        mesh3 = ax3.pcolormesh(mod_data[site_id]['time'], pairs_hc_unique_values,
                               np.transpose(mbe), cmap=cm.get_cmap('jet'))

        # m
        mesh4 = ax4.pcolormesh(mod_data[site_id]['time'], mod_data[site_id]['level_height'], np.transpose(mod_data[site_id]['aerosol_for_visibility']),
                                          vmin=0, vmax=30, cmap=cm.get_cmap('OrRd')) #  vmin=0, vmax=100

        # RH
        mesh5 = ax5.pcolormesh(mod_data[site_id]['time'], mod_data[site_id]['level_height'], np.transpose(mod_data[site_id]['RH'])*100,
                                          cmap = cm.get_cmap('Blues'), vmin=20.0, vmax=100.0)


        plt.subplots_adjust(right=0.8)

        # prettify
        for mesh, ax in zip((mesh1, mesh2, mesh3, mesh4, mesh5),(ax1, ax2, ax3, ax4, ax5)):
            ax.xaxis.set_major_formatter(DateFormatter('%H:%M'))
            ax.yaxis.label.set_size(10)
            ax.xaxis.label.set_size(10)
            ax.set_xlim([day, day + dt.timedelta(days=1)])
            ax.set_ylim([0, 1000.0])


        divider = make_axes_locatable(ax1)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(mesh1, cax=cax)

        divider = make_axes_locatable(ax2)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(mesh2, cax=cax)

        divider = make_axes_locatable(ax3)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(mesh3, cax=cax)

        divider = make_axes_locatable(ax4)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(mesh4, cax=cax)

        divider = make_axes_locatable(ax5)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(mesh5, cax=cax)

        ax1.get_xaxis().set_ticks([])
        ax2.get_xaxis().set_ticks([])
        ax3.get_xaxis().set_ticks([])
        ax4.get_xaxis().set_ticks([])

        eu.add_at(ax1, r'$\beta_{o}$', loc=2)
        eu.add_at(ax2, r'$\beta_{m}$', loc=2)
        eu.add_at(ax3, r'$log_{10}(\beta_{m}) - log_{10}(\beta_{o})$', loc=2)
        eu.add_at(ax4, r'$m$', loc=2)
        eu.add_at(ax5, r'$RH$', loc=2)

        ax0 = eu.fig_majorAxis(fig)
        ax0.set_xlabel('Time [HH:MM]', fontsize=10, labelpad=2)
        ax0.set_ylabel('Height [m]', fontsize=10, labelpad=10)

        plt.tight_layout(h_pad=0.1)

        plt.savefig(savedir + model_type + '-' + site + '-beta_o_beta_m_MBE_m_RH_' + day.strftime('%Y%m%d') + '.png')  # filename

        plt.close(fig)

    print 'END PROGRAM'

    return



if __name__ == '__main__':
    main()