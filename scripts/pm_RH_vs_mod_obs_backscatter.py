"""
Plot the PM10 or RH near surface difference, against mod - obs attenuated backscatter

Also calculate Mean Absolute Error (MAE) binned by delta_m

Created by Elliott Tues 09/05/17
"""

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import FormatStrFormatter

import numpy as np
import datetime as dt
from scipy.stats import spearmanr
from scipy.stats import pearsonr

from copy import deepcopy
import colorsys

import ellUtils as eu
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

    # Structure of statistics:
    # statistics[site]['r'] = [...]
    # statistics[site]['MBE'] = {'0-500': ..., '500-1000': ...}
    # statistics[site]['time'] = [...]

    if site_id not in statistics:

        # define site based lists to store the correlation results in
        statistics[site_id] = {'r': [], 'p': [],
                               'diff': [],
                               'aer_diff': [],
                               'aer_diff_normalsed_aer_obs': [],
                               'aer_mod': [],
                               'aer_obs': [],
                               'rh_diff': [],
                               'rh_mod': [],
                               'rh_obs': [],
                               'back_diff_log': [],
                               'back_diff_norm': [],
                               'back_diff_norm_normalised_back_obs': [],
                               'abs_back_diff_log': [],
                               'abs_back_diff_norm': [],
                               'back_obs': [],
                               'back_mod': [],
                               'RMSE': [],
                               'MBE': [],
                               'hr': [],
                               'datetime':[]}


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

def plot_back_point_diff(stats_site, savedir, model_type, ceil_gate_num, ceil, sampleSize=np.nan, corr=np.nan,
                         var_type=np.nan, c_type='rh_mod', savename='', extra=''):

    """
    Plot the rh or aer difference vs backscatter diff
    :return:
    """

    # stats_site = statistics[site_id]
    c_type = 'rh_mod'
    var_type = 'aerosol'

    # sample size for backscatter, aer and rh - help work out why some points don't show up
    back_bool = ~np.isnan(np.array(stats_site['back_diff_norm']))
    aer_bool = ~np.isnan(np.array(stats_site['aer_diff']))
    rh_bool = ~np.isnan(np.array(stats_site['rh_mod']))

    back_n = len(np.where(back_bool == True)[0])
    aer_n = len(np.where(aer_bool == True)[0])
    rh_n = len(np.where(rh_bool == True)[0])

    bool = back_bool & aer_bool & rh_bool
    sampleSize = len(np.where(bool == True)[0])

    # variable plotting against backscatter point diff
    if var_type == 'aerosol':
        # var_diff = stats_site['aer_diff']
        var_diff = stats_site['aer_diff_normalsed_aer_obs']
    elif var_type == 'RH':
        var_diff = stats_site['rh_diff']

    # backscatter point difference
    # back_point_diff = stats_site['back_diff_norm']
    back_point_diff = stats_site['back_diff_norm_normalised_back_obs']
    # back_point_diff = stats_site['back_diff_log']

    fig = plt.figure(figsize=(6, 3.5))
    ax = plt.subplot2grid((1, 1), (0, 0))

    # variable specific labels and names
    if var_type == 'RH':
        xlab = r'$Difference \/\mathrm{(RH_{ukv} - RH_{obs})}$'

    elif var_type == 'aerosol':
        xlab = r'$\frac{m_{MURK} - PM_{10}}{PM_{10}}}$'
        #xlab = r'$Difference \/\mathrm{(m_{MURK} - PM_{10})}$'

    c_min = 20
    c_max = 100.0

    # define the colormap
    cmap, norm = discrete_colour_map(c_min, c_max, 17)

    # plot data
    # scat = plt.scatter(var_diff, back_point_diff, c=stats_site[c_type], s=6, vmin=40.0, vmax=100.0, cmap=cmap, norm=norm)
    scat = plt.scatter(var_diff, back_point_diff, c=stats_site[c_type], s=6, vmin=c_min, vmax=c_max, cmap=cmap, norm=norm)

    # add 0 lines
    ax.axhline(linestyle='--', color='grey', alpha=0.5)
    ax.axvline(linestyle='--', color='grey', alpha=0.5)

    # ax.set_ylim([-5e-06, 5e-06])

    # add colourbar on the side
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(scat, cax=cax, norm=norm)
    # cbar.set_label(r'$RH\/[\%]$', labelpad=-38, y=1.075, rotation=0)
    cbar.set_label(r'$\mathrm{RH\/[\%]}$', labelpad=-36, y=1.075, rotation=0)
    # cbar.set_label(r'$m\/[/mu g \/kg^{-1}]$', labelpad=-38, y=1.075, rotation=0)


    ax.set_xlabel(xlab)
    # ax.set_ylim([-1e-5, 1e-5])
    # ax.set_ylim([-6.0e-6, 6.0e-6]) # paper 1
    # ax.set_ylim([-6.0e-5, 6.0e-5])
    ax.set_ylim([-1.2, 15.0])
    ax.set_xlim([-1.2, 6.0])
    # ax.set_ylabel(r'$Difference \/\mathrm{(log_{10}(\beta_m) - log_{10}(\beta_o))}$')
    # ax.set_ylabel(r'$Difference \/\mathrm{(\beta_m - \beta_o)}$')
    ax.set_ylabel(r'$\frac{\beta_m - \beta_o}{\beta_o}}$')
    #ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))

    # Fake a ScalarMappable so I can display a colormap
    # cmap, norm = mcolors.from_levels_and_colors(range(24 + 1), rgb)
    # sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
    # sm.set_array([])
    # fig.colorbar(sm)
    # plt.colorbar()

    fig.suptitle('n='+str(sampleSize)+': back_n='+str(back_n)+'; aer='+str(aer_n)+' rh_n='+str(rh_n))
    #fig.suptitle(ceil + '; n = ' + str(len(back_point_diff)) + '; r = ' + '{:1.2f}'.format(corr['r']) +
    #             '; p = ' + '%1.2f' % corr['p'])

    plt.tight_layout()
    plt.subplots_adjust(top=0.90)
    if savename is '':
        plt.savefig(savedir + model_type + '_' + var_type + '_diff_' + ceil + '_clearDays_gate' + str(ceil_gate_num) + '_c' + c_type +
               '_' + extra + '.png')  # filename
    else:
        plt.savefig(savedir + savename + extra)

    plt.close(fig)

    return

def plot_back_point_ratio(stats_site, savedir, model_type, ceil_gate_num, ceil, sampleSize, corr, var_type, c_type='hr', extra=''):

    """
    Plot the rh or aer difference vs backscatter diff
    :return:
    """

    # variable plotting against backscatter point diff
    if var_type == 'aerosol':
        var_ratio = stats_site['aer_diff']
    elif var_type == 'RH':
        var_ratio = stats_site['rh_diff']

    # backscatter point difference
    back_point_diff = stats_site['back_diff_norm']
    # back_point_diff = stats_site['back_diff_log']

    fig = plt.figure(figsize=(6, 3.5))
    ax = plt.subplot2grid((1, 1), (0, 0))

    # variable specific labels and names
    if var_type == 'RH':
        xlab = r'$Difference \/\mathrm{(RH_{ukv} - RH_{obs})}$'

    elif var_type == 'aerosol':
        xlab = r'$Difference \/\mathrm{(m_{MURK} - PM_{10})}$'

    # define the colormap
    cmap, norm = discrete_colour_map(40, 100, 13)

    # plot data
    scat = plt.scatter(var_ratio, back_point_diff, c=stats_site[c_type], s=6, vmin=40.0, vmax=100.0, cmap=cmap, norm=norm)

    # add 0 lines
    ax.axhline(linestyle='--', color='grey', alpha=0.5)
    ax.axvline(linestyle='--', color='grey', alpha=0.5)

    # ax.set_ylim([-5e-06, 5e-06])

    # add colourbar on the side
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(scat, cax=cax, norm=norm)
    cbar.set_label(r'$RH\/[\%]$', labelpad=-38, y=1.075, rotation=0)


    ax.set_xlabel(xlab)
    ax.set_ylim([-1e-5, 1e-5])
    # ax.set_ylabel(r'$Difference \/\mathrm{(log_{10}(\beta_m) - log_{10}(\beta_o))}$')
    ax.set_ylabel(r'$Difference \/\mathrm{(\beta_m - \beta_o)}$')
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))

    # Fake a ScalarMappable so I can display a colormap
    # cmap, norm = mcolors.from_levels_and_colors(range(24 + 1), rgb)
    # sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
    # sm.set_array([])
    # fig.colorbar(sm)
    # plt.colorbar()

    fig.suptitle(ceil + '; n = ' + str(sampleSize) + '; r = ' + '{:1.2f}'.format(corr['r']) +
                 '; p = ' + '%1.2f' % corr['p'])
    plt.tight_layout()
    plt.subplots_adjust(top=0.90)
    plt.savefig(savedir + 'point_diff/' +
                model_type + '_' + var_type + '_diff_' + ceil + '_clearDays_gate' + str(ceil_gate_num) + '_c' + c_type +
                '_' + extra + '.png')  # filename

    plt.close(fig)

    return

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

        rgb = [colorsys.hsv_to_rgb(i / (num_colours*3.0), 1.0, 1.0) for i in range(num_colours)]

        # rgb = colorsys.hsv_to_rgb(i / 72.0, 1.0, 1.0)
        # print(i, [round(255 * x) for x in rgb])

        return rgb

def discrete_colour_map(lower_bound, upper_bound, spacing):

    """Create a discrete colour map"""

    import matplotlib.pyplot as plt

    cmap = plt.cm.jet
    # extract all colors from the .jet map
    cmaplist = [cmap(i) for i in range(cmap.N)]
    # force the first color entry to be grey
    # cmaplist[0] = (.5, .5, .5, 1.0)
    # create the new map
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)

    # define the bins and normalize
    bounds = np.linspace(lower_bound, upper_bound, spacing)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    return cmap, norm

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
    savedir = maindir + 'figures/' + model_type + '/clearSkyPeriod/point_diff/'

    # data
    ceilDatadir = datadir + 'L1/'
    modDatadir = datadir + model_type + '/'
    rhDatadir = datadir + 'L1/'
    aerDatadir = datadir + 'LAQN/'

    # statistics to run
    pm10_stats = True
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

    # # pm10
    site = 'MR'
    ceil_id = 'CL31-C'
    # ceil = ceil_id + '_BSC_' + site
    ceil = ceil_id + '_' + site

    # rh
    #site = 'KSS45W'
    #ceil_id = 'CL31-A'
    ## ceil = ceil_id + '_BSC_' + site
    #ceil = ceil_id + '_' + site

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
    daystrList = ['20150414', '20150415', '20150421', '20150611', '20160119', '20160504', '20160823', '20160911', '20161125',
                  '20161129', '20161130', '20161204']

    if site == 'KSS45W':
    # KSS45W days
        daystrList = ['20150414', '20150415', '20150421', '20150611']

    # MR calib days
    elif site == 'MR':
        daystrList = ['20150414', '20150415', '20150421', '20150611', '20160119', '20160504', '20160823', '20160911',
                      '20161125', '20161129', '20161130', '20161204'] # with 19/1/16

    # daystrList = ['20150414']
    elif site == 'NK':
    # NK_D calib days
        daystrList = ['20150414', '20150415', '20150421', '20150611', '20160504']

    elif site == 'IMU':
    # IMU days
        daystrList = ['20160504', '20160823', '20160911', '20161125',
                      '20161129', '20161130', '20161204']

    days_iterate = dateList_to_datetime(daystrList)


    # ceilometer gate number to use for backscatter comparison
    # 1 - noisy
    # 2 - more stable
    # see Kotthaus et al (2016) for more.
    ceil_gate_num = 2

    # define statistics dictionary
    statistics = {}
    sampleSize = 0 # add to this

    # aerFO run and bsc obs calibrate version number
    version=1.0

    # ==============================================================================
    # Read data
    # ==============================================================================

    # Read Ceilometer metadata

    # ceilometer list to use
    ceilsitefile = 'CeilsCSVfull.csv'
    ceil_metadata = FO.read_ceil_metadata(datadir, ceilsitefile)

    ceil_data_i = {site: ceil_metadata[site]}

    for day in days_iterate:

        print 'day = ' + day.strftime('%Y-%m-%d')

        # Read UKV forecast and automatically run the FO

        # extract MURK aerosol and calculate RH for each of the sites in the ceil metadata
        # reads all london model data, extracts site data, stores in single dictionary
        mod_data = FO.mod_site_extract_calc(day, ceil_data_i, modDatadir, model_type, res, 910, version=version, allvars=True)

        # mod_data_all = FO.mod_site_extract_calc(day, ceil_metadata, modDatadir, model_type, res, 910, version=0.2, allvars=True)


        # Read ceilometer backscatter

        # will only read in data is the site is there!
        # ToDo Remove the time sampling part and put it into its own function further down.
        bsc_obs = FO.read_ceil_obs(day, site_bsc, ceilDatadir, mod_data, calib=True, version=version)


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
            # get the nearest rh_height based on the instrument height, however if rh obs are not being read in,
            # default it to 5 m, near the surface, and the closest gate to the PM10 obs.
            if rh_stats == True:
                mod_rh_height_idx = eu.nearest(mod_data[site]['level_height'], site_rh[rh_instrument])[1]
            else:
                mod_rh_height_idx = 0

            ceil_height_idx, mod_height_idx =\
                 get_nearest_ceil_mod_height_idx(mod_data[site_id]['level_height'], bsc_site_obs['height'], ceil_gate_num)

            # create entry in the dictionary if one does not exist
            statistics = create_stats_entry(site_id, statistics)

            # for each hour possible in the day
            for t in np.arange(0, 24):

                hr = str(t)

                # add the hour to statstics dict
                statistics[site_id]['hr'] += [t]
                statistics[site_id]['datetime'] += [mod_data[site]['time'][t]]

                # extract out all unique pairs below the upper height limit
                # these are time and height matched now
                #obs_x = bsc_site_obs['backscatter'][t, obs_hc_unique_pairs]
                #mod_y = mod_data[site_id]['backscatter'][t, mod_hc_unique_pairs]

                # extract pairs of values used in statistics
                if pm10_stats == True:
                    pm10_i = pm10['PM10_'+site]['pm_10'][t]
                    # murk_i = mod_data[site]['aerosol_for_visibility'][t, 0] # 0th height = 5 m

                if rh_stats == True:
                    rh_obs_i = rh_obs[rh_instrument]['RH'][t]

                # get murk_i and rh_i regardless if pm10_stats and rh_stats is True or not
                murk_i = mod_data[site]['aerosol_concentration_dry_air'][t, 0]  # 0th height = 5 m
                rh_mod_i = mod_data[site]['RH'][t, mod_rh_height_idx] * 100.0  # convert from [fraction] to [%]

                obs_back_i = bsc_site_obs['backscatter'][t, ceil_height_idx]
                mod_back_i = mod_data[site]['backscatter'][t, mod_height_idx]

                statistics[site_id]['back_obs'] += [obs_back_i]
                statistics[site_id]['back_mod'] += [mod_back_i]

                statistics[site_id]['aer_mod'] += [murk_i]
                statistics[site_id]['rh_mod'] += [rh_mod_i]

                # STATISTICS
                # ---------------

                # length of aer_diff[hr] and ['back_point_diff'] hour should and MUST be the same length
                # such that their idx positions line up
                if pm10_stats == True:
                    statistics[site_id]['aer_diff'] += [murk_i - pm10_i]
                    statistics[site_id]['aer_diff_normalsed_aer_obs'] += [(murk_i - pm10_i)/pm10_i]
                    statistics[site_id]['aer_obs'] += [pm10_i]

                    # if the difference pairs do not posses an NaN (and will therefore be plotted), add 1 to sample size
                    if ~np.isnan(murk_i - pm10_i) & ~np.isnan(np.log10(mod_back_i) - np.log10(obs_back_i)):
                        sampleSize += 1

                if rh_stats == True:
                    statistics[site_id]['rh_diff'] += [rh_mod_i - rh_obs_i]
                    statistics[site_id]['rh_obs'] += [rh_obs_i]

                    # if the difference pairs do not posses an NaN (and will therefore be plotted), add 1 to sample size
                    if ~np.isnan(rh_mod_i - rh_obs_i) & ~np.isnan(np.log10(mod_back_i) - np.log10(obs_back_i)):
                        sampleSize += 1


                # all extra stats slots
                statistics[site_id]['back_diff_log'] += [np.log10(mod_back_i) - np.log10(obs_back_i)]
                statistics[site_id]['back_diff_norm'] += [mod_back_i - obs_back_i]
                statistics[site_id]['back_diff_norm_normalised_back_obs'] += [(mod_back_i - obs_back_i)/obs_back_i]
                statistics[site_id]['abs_back_diff_norm'] += [np.abs(mod_back_i - obs_back_i)]
                statistics[site_id]['abs_back_diff_log'] += [np.abs(np.log10(mod_back_i) - np.log10(obs_back_i))]


    # calculate basic statistics for a table in paper 1
    basic_stats = {}
    # for var in ['back', 'rh', 'aer']:
    for var in ['rh', 'aer']:
        var_mod = var+'_mod'
        var_obs = var+'_obs'
        var_diff = var+'_diff'

        basic_stats[var] = {}

        basic_stats[var]['mod_mean'] = np.nanmean(statistics[site_id][var_mod])
        basic_stats[var]['mod_median'] = np.nanmedian(statistics[site_id][var_mod])
        basic_stats[var]['mod_stdev'] = np.nanstd(statistics[site_id][var_mod])
        basic_stats[var]['mod_IQR'] = np.nanpercentile(statistics[site_id][var_mod], 75) - np.nanpercentile(statistics[site_id][var_mod], 25)

        basic_stats[var]['obs_mean'] = np.nanmean(statistics[site_id][var_obs])
        basic_stats[var]['obs_median'] = np.nanmedian(statistics[site_id][var_obs])
        basic_stats[var]['obs_stdev'] = np.nanstd(statistics[site_id][var_obs])
        basic_stats[var]['obs_IQR'] = np.nanpercentile(statistics[site_id][var_obs], 75) - np.nanpercentile(statistics[site_id][var_obs], 25)

        basic_stats[var]['mean_diff'] = np.nanmean(statistics[site_id][var_diff])
        basic_stats[var]['pct_diff'] = []
        #basic_stats[var]['ratio'] = np.nanmean([i / j for i, j in [statistics[site_id][var_mod], statistics[site_id][var_obs]]])


    # # quick remove RH > x
    # for val in np.arange(40, 110, 10):
    #
    #     stats_copy = deepcopy(statistics[site_id])
    #
    #     idx_lt = np.where(np.array(stats_copy['rh_mod']) >= val)
    #
    #     for key in stats_copy.iterkeys():
    #         if len(stats_copy[key]) != 0:
    #             for i in idx_lt[0]:
    #                 stats_copy[key][i] = np.nan
    #     s = np.sum(~np.isnan(stats_copy['aer_diff']))
    #     r, p = spearmanr(stats_copy['aer_diff'], stats_copy['back_diff_norm'],
    #                      nan_policy='omit')
    #
    #     print 'val =' + str(val)
    #     print r
    #     print p
    #     print s
    #     print ''


    # do correlation
    if rh_stats == True:
        corr = {}
        corr['r'], corr['p'] = spearmanr(statistics[site_id]['rh_diff'], statistics[site_id]['back_diff_norm'], nan_policy='omit')
        a1 = np.array(statistics[site_id]['aer_mod']) # extract the rh difference dataset from statistics

    if pm10_stats == True:
        corr = {}
        corr['r'], corr['p'] = spearmanr(statistics[site_id]['aer_diff'], statistics[site_id]['back_diff_norm'], nan_policy='omit')
        a1 = np.array(statistics[site_id]['aer_diff']) # extract the aerosol difference dataset from statistics
    # plot!

    # ----------------

    # pearson correlation
    b1 = np.array(statistics[site_id]['back_diff_norm'])

    a1_idx = np.where(np.isnan(a1))
    b1_idx = np.where(np.isnan(b1))

    a1[b1_idx] = np.nan
    b1[a1_idx] = np.nan

    an_idx = ~np.isnan(a1)

    pearsonr(a1[an_idx],b1[an_idx])

    # --------------

    # calculate mean absolute error, binned by delta_m (n = sample size)
    MAE = {'MAE': [], 'AE': [], 'stddev MAE': [], 'aer_diff': [], 'n': []}

    bin = 25.0
    bin_start = -50
    bin_end = 125

    bin_start_all = np.arange(bin_start, bin_end, bin) # only start of each bin
    bin_end_all = np.arange(bin_start+bin, bin_end+bin, bin) # end of each bin
    pos = np.arange(bin_start + (0.5*bin), bin_end, bin) # centre position of each bin
    bin_widths = np.repeat(bin, len(pos)) # width of every bin

    # -50 to 100
    for bin_s_i, bin_e_i in zip(bin_start_all, bin_end_all):

        # find all abs_back_diff_norm for current bin
        m_bin_bool = np.logical_and(statistics[site_id]['aer_diff'] >= bin_s_i,
                                    statistics[site_id]['aer_diff'] <= bin_e_i)
        m_bin_idx = np.where(m_bin_bool == True)[0]

        # fill MAE statistics
        MAE['MAE'] += [np.nanmean(np.array(statistics[site_id]['abs_back_diff_norm'])[m_bin_idx])]
        # MAE['AE'] += [np.array(statistics[site_id]['abs_back_diff_norm'])[m_bin_idx]]
        # temp
        a = np.array(statistics[site_id]['back_diff_norm'])[m_bin_idx]
        MAE['AE'] += [a[~np.isnan(a)]]

        MAE['stddev MAE'] += [np.nanstd(np.array(statistics[site_id]['abs_back_diff_norm'])[m_bin_idx])]
        MAE['aer_diff'] += [[bin_s_i, bin_e_i]]
        MAE['n'] += [len(m_bin_idx)]
        # MAE['n'] += [1 if type(m_bin_idx)...

    # plot MAE
    # work around needed for end box as it only has two values.
    fig, ax = plt.subplots(1, 1, figsize=(6, 3.5))

    # if the last box of the boxplot only has two points, scatter the points instead of making the boxplot
    if len(MAE['AE'][-1]) < 4:
        plt.boxplot(MAE['AE'][:-1], widths=bin, positions=pos[:-1], whis=[5, 95], sym='x', hold=True)
        # plt.scatter([pos[-1],pos[-1]], MAE['AE'][-1], marker='x', color='black')
        plt.scatter(np.repeat(pos[-1], len(MAE['AE'][-1])), MAE['AE'][-1], marker='x', color='black')
        plt.scatter(pos[:-1], MAE['MAE'][:-1], marker='o', color='blue', edgecolor='blue')

    else:
        plt.boxplot(MAE['AE'], widths=bin, positions=pos, whis=[5, 95], sym='x', hold=True)
        plt.scatter(pos, MAE['MAE'], marker='o', color='blue', edgecolor='blue')
        # raise ValueError('End box n != 2, make sure data is being plotted correctly \n'
        #                  'a specific work around is currently implemented!')


    # prettify
    ax.set_xlim([bin_start, bin_end])
    ax.set_xticks(np.arange(bin_start, bin_end+bin, bin))
    ax.set_xticklabels(np.arange(bin_start, bin_end+bin, bin))
    ax.set_xlabel(r'$Difference \/\mathrm{(m_{MURK} - PM_{10})}$')
    ax.set_ylabel(r'$Absolute \/\/Error\/\/\/ \mathrm{(|\beta_m - \beta_o|)}$')
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
    #ax.set_ylim([0.0, 6.0e-06]) #orig
    ax.set_ylim([-10.0e-06, 10.0e-06])

    # add sample size at the top of plot for each box and whiskers
    # pos_t = np.arange(numBoxes) + 1
    upperLabels = [str(np.round(n, 2)) for n in MAE['n']]
    weights = ['bold', 'semibold']
    for tick, label in zip(range(len(pos)), ax.get_xticklabels()):
        k = tick % 2
        ax.text(pos[tick], 6.0e-06 - (6.0e-06 * 0.05), upperLabels[tick],
                 horizontalalignment='center', size='x-small')

    plt.subplots_adjust(left=0.1)
    plt.tight_layout()
    plt.savefig(savedir + 'mae/' +
                'MAE_clearSkyPeriod_median_IQR_meanDots2.png')
    plt.close(fig)




    # ---------------

    # pass in statistics with site id!
    if pm10_stats == True:
        plot_back_point_diff(statistics[site_id],
                                   savedir, model_type, ceil_gate_num, ceil, sampleSize, corr, var_type='aerosol',
                                   c_type='rh_mod', extra = '_norm_back_newpress')

    if rh_stats == True:
        plot_back_point_diff(statistics[site_id],
                                   savedir, model_type, ceil_gate_num, ceil, sampleSize, corr, var_type='RH',
                                   c_type='rh_mod')


    plt.close('all')

    # # scatter diff in backscatter vs pm10 to see if it should be linear
    # r, p = spearmanr(statistics[site_id]['aer_obs'], statistics[site_id]['back_obs'], nan_policy='omit')
    # plt.scatter(statistics[site_id]['aer_obs'], statistics[site_id]['back_obs'])
    # plt.xlabel('pm10')
    # plt.ylabel('log back diff')
    # plt.ylim([min(statistics[site_id]['back_obs']), max(statistics[site_id]['back_obs'])])
    # plt.suptitle('r=' + str(r) + '; p='+str(p))


    # # plot heights a(i) and height(i)
    # orog = {'NK':28.6743, 'KSS45W': 21.969, 'RGS': 19.9921, 'MR': 43.1841, 'IMU': 29.6641}
    #
    # a = np.array([5.00000,21.6667,45.0000,75.0000,111.667,155.000,205.000,261.667
    #     ,325.000,395.000,471.667,555.000,645.000,741.667,845.000,955.000
    #     ,1071.67,1195.00,1325.00,1461.67,1605.00,1755.00,1911.67,2075.00
    #     ,2245.00,2421.67,2605.00,2795.00,2991.67,3195.00,3405.00,3621.67
    #     ,3845.00,4075.00,4311.67,4555.00,4805.00,5061.67,5325.00,5595.00
    #     ,5871.67,6155.01,6445.15,6742.49,7047.82,7362.36,7687.92,8026.93
    #     ,8382.58,8758.92,9160.94,9594.76,10067.7,10588.3,11166.8,11814.9
    #     ,12546.0,13375.7,14321.3,15402.7,16642.0,18063.9,19696.0,21568.9
    #     ,23716.1,26174.7,28985.5,32192.7,35845.0,40000.0])
    #
    # b = np.array([0.999424,0.997504,0.994820,0.991375,0.987171,0.982215,0.976512,0.970069
    #     ,0.962893,0.954993,0.946377,0.937057,0.927043,0.916346,0.904981,0.892961
    #     ,0.880300,0.867014,0.853118,0.838632,0.823572,0.807957,0.791808,0.775146
    #     ,0.757992,0.740368,0.722298,0.703807,0.684920,0.665662,0.646062,0.626146
    #     ,0.605944,0.585484,0.564799,0.543919,0.522876,0.501704,0.480437,0.459110
    #     ,0.437758,0.416418,0.395119,0.373871,0.352663,0.331463,0.310213,0.288832
    #     ,0.267223,0.245272,0.222861,0.199882,0.176257,0.151965,0.127085,0.101852,
    #      0.0767339,0.0525319,0.0305214, 0.0126308, 0.00167859, 0.00000, 0.00000, 0.00000
    #     , 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000])
    #
    # fig = plt.figure(figsize=(5,5))
    #
    # hmod_abs = {}
    # for key, orog_i in orog.iteritems():
    #
    #     hmod_abs[key] = a + (b * orog_i) - orog_i
    #
    #     plt.plot(a - hmod_abs[key], a, label=key, color = site_bsc_colours[key])
    #
    #     plt.ylim([0.0, 2000.0])
    #     plt.xlim([0.0, 10.0])
    #     plt.ylabel('Height [m]')
    #     plt.xlabel(r'$h_{sea} - h_{surface}$' + ' [m]')
    #     plt.legend()
    #     plt.tight_layout()
    #     plt.savefig(savedir + 'a_vs_height(i).png')


    # ----------------------

    # SCATTER - beta_m vs beta_o, coloured by rh_obs
    # simple plot for paper 1, given the reviewers comments

    # extract x and y data so the plotting code is more readable
    x_data = statistics[site]['back_obs']
    y_data = statistics[site]['back_mod']

    c_key = 'rh_diff'
    c_data = statistics[site][c_key]

    # c_min = -40.0
    # c_max = 100.0

    # c_min = 40.0
    # c_max = 100.0

    c_min = -30.0
    c_max = 30.0

    fig = plt.figure(figsize=(7.5, 4))
    ax = plt.subplot2grid((1, 1), (0, 0))

    # define the colormap
    cmap, norm = discrete_colour_map(c_min, c_max, 13)

    # plot data
    # scat = plt.scatter(x_data, y_data, c=statistics[site]['rh_obs'], s=6, vmin=40.0, vmax=100.0, cmap=cmap, norm=norm)
    scat = plt.scatter(x_data, y_data, c=c_data, s=6, vmin=c_min, vmax=c_max, cmap=cmap, norm=norm)

    # add 1-to-1 line
    plt.plot([np.nanmin([x_data, y_data]), np.nanmin([x_data, y_data])],
             [np.nanmax([x_data, y_data]), np.nanmax([x_data, y_data])],
             linestyle='--', color='grey', alpha=0.5)

    # ax.set_ylim([-5e-06, 5e-06])

    # add colourbar on the side
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(scat, cax=cax, norm=norm)
    if c_key[:2] == 'rh':
        cbar.set_label(r'$\Delta RH\/[\%]$', labelpad=-38, y=1.075, rotation=0)
    elif c_key[:3] == 'aer' :
        cbar.set_label(r'$m\/[\mu g\/ kg^{-1}]$', labelpad=-38, y=1.075, rotation=0)

    ax.set_xlabel(r'$\mathrm{\beta_o}$')
    ax.set_ylabel(r'$\mathrm{\beta_m}$')
    ax.set_ylim([np.nanmin(y_data), np.nanmax(y_data)])
    #ax.set_ylim([np.nanmin(y_data), np.nanmax(y_data)])
    ax.set_xlim([np.nanmin(x_data), np.nanmax(x_data)])
    #ax.set_ylim([-1e-5, 1e-5])
    #ax.set_xlim([-1e-5, 1e-5])
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1e'))

    #fig.suptitle(ceil + '; n = ' + str(sampleSize) + '; r = ' + '{:1.2f}'.format(corr['r']) +
    #             '; p = ' + '%1.2f' % corr['p'])
    plt.tight_layout()
    plt.subplots_adjust(top=0.90)
    plt.savefig(savedir + 'beta_o_vs_beta_m/' +
                model_type + '_scatter_beta_m_vs_beta_o_' + ceil + '_clearDays_gate' + str(ceil_gate_num) + '_c_'+c_key+'_'+
                '.png')  # filename

    plt.close(fig)


    # -----------------------------------------------

    # RATIO SCATTER - beta_m vs beta_o, coloured by rh_obs
    # simple plot for paper 1, given the reviewers comments

    # extract x and y data so the plotting code is more readable
    x_data = statistics[site]['rh_mod']
    y_data = np.array(statistics[site]['back_mod']) / np.array(statistics[site]['back_obs'])
    c_key = 'aer_mod'
    c_data = statistics[site][c_key]

    # # aer_diff
    # c_min = -40.0
    # c_max = 100.0

    # aer mod
    c_min = 0.0
    c_max = 100.0

    # # rh diff
    # c_min = -30.0
    # c_max = 30.0


    fig = plt.figure(figsize=(6, 4.5))
    ax = plt.subplot2grid((1, 1), (0, 0))

    # define the colormap
    # cmap, norm = discrete_colour_map(c_min, c_max, 13)
    cmap, norm = discrete_colour_map(c_min, c_max, 11)

    # plot data
    # scat = plt.scatter(x_data, y_data, c=statistics[site]['rh_obs'], s=6, vmin=40.0, vmax=100.0, cmap=cmap, norm=norm)
    scat = plt.scatter(x_data, y_data, c=c_data, s=6, vmin=c_min, vmax=c_max, cmap=cmap, norm=norm)

    # add 1-to-1 line
    plt.plot([np.nanmin([x_data, y_data]), np.nanmin([x_data, y_data])],
             [np.nanmax([x_data, y_data]), np.nanmax([x_data, y_data])],
             linestyle='--', color='grey', alpha=0.5)

    # # add ratio=1 line
    ax.axhline(1, linestyle='--', color='grey', alpha=0.5)
    # ax.axvline(linestyle='--', color='grey', alpha=0.5)

    # ax.set_ylim([-5e-06, 5e-06])
    ax.set_ylim([0.0, 3.5])
    # ax.set_ylim([np.nanmin([x_data, y_data]), np.nanmax([x_data, y_data])])
    ax.set_xlim([np.nanmin([x_data, y_data]), np.nanmax([x_data, y_data])])
    # ax.set_xlim([-30.0, 30.0])

    # add colourbar on the side
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(scat, cax=cax, norm=norm)
    if c_key[:2] == 'rh':
        cbar.set_label(r'$\Delta RH\/[\%]$', labelpad=-38, y=1.075, rotation=0)
    elif c_key[:3] == 'aer' :
        # cbar.set_label(r'$\Delta\/m\/[\mu g\/ kg^{-1}]$', labelpad=-38, y=1.075, rotation=0)
        cbar.set_label(r'$m_{MURK}\/[\mu g\/ kg^{-1}]$', labelpad=-38, y=1.075, rotation=0)


    ax.set_xlabel(r'$RH$')
    ax.set_ylabel(r'$\mathrm{\beta_m / \beta_o}$')

    plt.tight_layout()
    plt.subplots_adjust(top=0.90)
    plt.savefig(savedir + 'beta_o_vs_beta_m/' +
                model_type + '_scatter_beta_m_beta_o_ratio_vs_aer_diff_' + ceil + '_clearDays_gate' + str(ceil_gate_num) + '_c_'+c_key+'_'+
                '.png')  # filename

    plt.close(fig)



    print 'END PROGRAM'