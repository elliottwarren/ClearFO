"""
Plot the PM10 or RH near surface difference, against mod - obs attenuated backscatter

Created by Elliott Tues 09/05/17
"""

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm

import numpy as np
import datetime as dt
from scipy.stats import spearmanr

from copy import deepcopy
import colorsys

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
                               'aer_diff': {},
                               'aer_mod': {},
                               'aer_obs': {},
                               'rh_diff': {},
                               'rh_mod': {},
                               'rh_obs': {},
                               'back_point_diff': {},
                               'RMSE': {},
                               'MBE': {}}

        for key in statistics[site_id].iterkeys():
            for hr in hrs:
                statistics[site_id][key][str(hr)] = []

    return statistics

def dateList_to_datetime(dayList):

    """ Convert list of string dates into datetimes """

    datetimeDays = []

    for d in dayList:

        datetimeDays += [dt.datetime(int(d[0:4]), int(d[4:6]), int(d[6:8]))]

    return datetimeDays

def get_nearest_ceil_mod_height_idx(mod_height, obs_height, ceil_gate_num):

    """
    Returns the ceilometer height index for the ceilometer gate number, and the idx for
    the nearest model height.

    :param mod_height:
    :param obs_height:
    :param ceil_gate_num:
    :return:
    """

    # idx and height for ceil
    ceil_gate_idx = ceil_gate_num - 1

    ceil_gate_height = obs_height[ceil_gate_idx]

    # find nearest height and idx for mod
    a = eu.nearest(mod_height, ceil_gate_height)

    mod_pair_height = a[0]
    mod_height_idx = a[1]
    height_diff = a[2]


    return ceil_gate_idx, mod_height_idx

def plot_back_point_diff(var_diff, back_point_diff, savedir, model_type, ceil_gate_num, ceil, sampleSize, corr, var_type):

    """
    Plot the rh or aer difference vs backscatter diff
    :return:
    """

    rgb = colour_range(24)

    fig = plt.figure(figsize=(6, 3.5))
    ax = plt.subplot2grid((1, 1), (0, 0))

    # variable specific labels and names
    if var_type == 'RH':
        xlab = r'$Difference \/\mathrm{(RH_{ukv} - RH_{obs})}$'

    elif var_type == 'aerosol':
        xlab = r'$Difference \/\mathrm{(m_{MURK} - PM_{10})}$'

    # plot each hour
    for t in np.arange(0, 24):

        hr = str(t)
        hr_colour = rgb[t]
        plt.scatter(var_diff[hr], back_point_diff[hr], color=hr_colour, s=6)

    ax.set_xlabel(xlab)
    ax.set_ylabel(r'$Difference \/\mathrm{(log_{10}(\beta_m) - log_{10}(\beta_o))}$')

    # Fake a ScalarMappable so I can display a colormap
    cmap, norm = mcolors.from_levels_and_colors(range(24 + 1), rgb)
    sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm)

    fig.suptitle(ceil + '; n = ' + str(sampleSize) + '; r = ' + '{:1.2f}'.format(corr['r']) +
                 '; p = ' + '%1.2f' % corr['p'])
    plt.tight_layout()
    plt.subplots_adjust(top=0.90)
    plt.savefig(savedir + 'point_diff/' +
                model_type + '_' + var_type + '_diff_' + ceil + '_clearDays_gate' + str(ceil_gate_num) + 'v0.2.png')  # filename

    return fig

def plot_back_point_diff_6hr(var_diff, back_point_diff, savedir, model_type, ceil_gate_num, ceil, sampleSize, corr, var_type):

    """
    Plot the rh or aer difference vs backscatter diff
    :return:
    """

    rgb = colour_range(24)

    fig = plt.figure(figsize=(6, 3.5))
    ax = plt.subplot2grid((1, 1), (0, 0))

    # variable specific labels and names
    if var_type == 'RH':
        xlab = r'$Difference \/\mathrm{(RH_{ukv} - RH_{obs})}$'

    elif var_type == 'aerosol':
        xlab = r'$Difference \/\mathrm{(m_{MURK} - PM_{10})}$'

    for t in [0,6,12,18]:

        t_range = np.arange(t, t+6)

        for t_i in t_range:

            hr = str(t_i)

            # hr_colour = rgb[t]

            plt.scatter(var_diff[hr], back_point_diff[hr], s=4)

        ax.set_xlabel(xlab)
        ax.set_ylabel(r'$Difference \/\mathrm{(log_{10}(\beta_m) - log_{10}(\beta_o))}$')

        # # Fake a ScalarMappable so I can display a colormap
        # cmap, norm = mcolors.from_levels_and_colors(range(6 + 1), rgb)
        # sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
        # sm.set_array([])
        # fig.colorbar(sm)

        fig.suptitle(ceil + '; n = ' + str(sampleSize) + '; r = ' + '{:1.2f}'.format(corr['r']) +
                     '; p = ' + '%1.2f' % corr['p'])
        plt.tight_layout()
        plt.subplots_adjust(top=0.90)
        plt.savefig(savedir + 'point_diff/' +
                    model_type + '_' + var_type + '_diff_' + ceil + '_clearDays_gate' + str(ceil_gate_num) + 't'+str(t)+'.png')  # filename

    return fig

def colour_range(num_colours=24.0):

    """Makes a simple range of colours"""

    for i in range(num_colours):

        rgb = [colorsys.hsv_to_rgb(i / 72., 1.0, 1.0) for i in range(num_colours)]

        # rgb = colorsys.hsv_to_rgb(i / 72.0, 1.0, 1.0)
        # print(i, [round(255 * x) for x in rgb])

        return rgb

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
    pm10_stats = False
    rh_stats = True

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


    site = 'KSS45W'
    ceil_id = 'CL31-A'
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

    # day list
    # clear sky days (5 Feb 2015 - 31 Dec 2016)
    # daystrList = ['20150414', '20150415', '20150421', '20150611', '20160504', '20160823', '20160911', '20161125',
    #               '20161129', '20161130', '20161204']

    # KSS45W days
    daystrList = ['20150414', '20150415', '20150421', '20150611']

    # MR calib days
    # daystrList = ['20150414', '20150415', '20150421', '20150611', '20160504', '20160823', '20160911', '20161125',
    #              '20161129', '20161130', '20161204']

    # NK_D calib days
    # daystrList = ['20150414', '20150415', '20150421', '20150611', '20160504']

    # IMU days
    # daystrList = ['20160504', '20160823', '20160911', '20161125',
    #              '20161129', '20161130', '20161204']

    days_iterate = dateList_to_datetime(daystrList)



    # correlation max height
    corr_max_height = 2000

    # ceilometer gate number to use for backscatter comparison
    # 1 - noisy
    # 2 - more stable
    # see Kotthaus et al (2016) for more.
    ceil_gate_num = 2

    # define statistics dictionary
    statistics = {}
    sampleSize = 0 # add to this


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
        mod_data = FO.mod_site_extract_calc(day, ceil_metadata, modDatadir, model_type, res, 910, version=0.2)

        # Read ceilometer backscatter

        # will only read in data is the site is there!
        # ToDo Remove the time sampling part and put it into its own function further down.
        bsc_obs = FO.read_ceil_obs(day, site_bsc, ceilDatadir, mod_data, calib=True)

        bsc_obs_uncal = FO.read_ceil_obs(day, site_bsc, ceilDatadir, mod_data, calib=False)

        if pm10_stats == True:
            # read in PM10 data and extract data for the current day
            pm10 = FO.read_pm10_obs(site_aer, aerDatadir, mod_data)

        # read in RH data
        if rh_stats == True:
            rh_obs = FO.read_all_rh_obs(day, site_rh, rhDatadir, mod_data)

        # ==============================================================================
        # Process
        # ==============================================================================

        # extract the single site out of bsc_obs
        # if ceil in bsc_obs:
        #     bsc_site_obs = bsc_obs[ceil]
        # else:
        #     bsc_site_obs = bsc_obs['CL31-E_BSC_NK']

        if ceil in bsc_obs:

            bsc_site_obs = bsc_obs[ceil]

            # short site id that matches the model id
            site_id = site.split('_')[-1]
            print '     Processing for site: ' + site_id

            # get the ceilometer and model height index for the define ceilometer range gate
            if rh_stats == True:
                mod_rh_height_idx = eu.nearest(mod_data[site]['level_height'], site_rh[rh_instrument])[1]

            ceil_height_idx, mod_height_idx =\
                 get_nearest_ceil_mod_height_idx(mod_data[site_id]['level_height'], bsc_site_obs['height'], ceil_gate_num)

            # create entry in the dictionary if one does not exist
            statistics = create_stats_entry(site_id, statistics)

            # for each hour possible in the day
            for t in np.arange(0, 24):

                hr = str(t)

                # extract out all unique pairs below the upper height limit
                # these are time and height matched now
                #obs_x = bsc_site_obs['backscatter'][t, obs_hc_unique_pairs]
                #mod_y = mod_data[site_id]['backscatter'][t, mod_hc_unique_pairs]

                # extract pairs of values used in statistics
                if pm10_stats == True:
                    pm10_i = pm10['PM10_'+site]['pm_10'][t]
                    murk_i = mod_data[site]['aerosol_for_visibility'][t, 0] # 0th height = 5 m

                if rh_stats == True:
                    rh_obs_i = rh_obs[rh_instrument]['RH'][t]
                    rh_mod_i = mod_data[site]['RH'][t, mod_rh_height_idx] * 100.0 # convert from [fraction] to [%]

                obs_back_i = bsc_site_obs['backscatter'][t, ceil_height_idx]
                mod_back_i = mod_data[site]['backscatter'][t, mod_height_idx]

                # STATISTICS
                # ---------------

                # length of aer_diff[hr] and ['back_point_diff'] hour should and MUST be the same length
                # such that their idx positions line up
                if pm10_stats == True:
                    statistics[site_id]['aer_diff'][hr] += [murk_i - pm10_i]

                    # if the difference pairs do not posses an NaN (and will therefore be plotted), add 1 to sample size
                    if ~np.isnan(murk_i - pm10_i) & ~np.isnan(np.log10(mod_back_i) - np.log10(obs_back_i)):
                        sampleSize += 1

                if rh_stats == True:
                    statistics[site_id]['rh_diff'][hr] += [rh_mod_i - rh_obs_i]

                    # if the difference pairs do not posses an NaN (and will therefore be plotted), add 1 to sample size
                    if ~np.isnan(rh_mod_i - rh_obs_i) & ~np.isnan(np.log10(mod_back_i) - np.log10(obs_back_i)):
                        sampleSize += 1


                # all extra stats slots
                statistics[site_id]['back_point_diff'][hr] += [np.log10(mod_back_i) - np.log10(obs_back_i)]

                statistics[site_id]['aer_mod'][hr] += murk_i
                statistics[site_id]['aer_obs'][hr] += pm10_i
                statistics[site_id]['rh_mod'][hr]  += rh_mod_i
                statistics[site_id]['rh_obs'][hr]  += rh_obs_i


    # gather up statistics...


    # do correlation
    if rh_stats == True:
        corr = {}
        corr['rh_diff_all'] = [i for sublist in statistics[site_id]['rh_diff'].values() for i in sublist]
        corr['back_diff_all'] = [i for sublist in statistics[site_id]['back_point_diff'].values() for i in sublist]
        corr['r'], corr['p'] = spearmanr(corr['rh_diff_all'], corr['back_diff_all'], nan_policy='omit')

    if pm10_stats == True:
        corr = {}
        corr['aer_diff_all'] = [i for sublist in statistics[site_id]['aer_diff'].values() for i in sublist]
        corr['back_diff_all'] = [i for sublist in statistics[site_id]['back_point_diff'].values() for i in sublist]
        corr['r'], corr['p'] = spearmanr(corr['aer_diff_all'], corr['back_diff_all'], nan_policy='omit')

    # plot!

    #"""
    #plot the correlation statistics (\beta_m, site vs \beta_o, site) and save.#
    #
    #:return: fig
    #"""

#    fig = plot_diff(diff, savedir, site_bsc_colours, model_type)
    if pm10_stats == True:
        fig = plot_back_point_diff(statistics[site_id]['aer_diff'], statistics[site_id]['back_point_diff'],
                                   savedir, model_type, ceil_gate_num, ceil, sampleSize, corr, var_type='aerosol')

    if rh_stats == True:
        fig = plot_back_point_diff(statistics[site_id]['rh_diff'], statistics[site_id]['back_point_diff'],
                                   savedir, model_type, ceil_gate_num, ceil, sampleSize, corr, var_type='RH')


    plt.close('all')

















    return

if __name__ == '__main__':
    main()























print 'END PROGRAM'
