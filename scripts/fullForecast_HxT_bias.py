"""
Plot the height x time bias profiles of the full forecast profiles together (3 profile sets)
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.dates import DateFormatter

import datetime as dt
import numpy as np
from copy import deepcopy

import ceilUtils as ceil
import ellUtils as eu
from forward_operator import FOUtils as FO
from forward_operator import FOconstants as FOcon

from mod_obs_stats_plot import unique_pairs

def read_ceil_obs(day, site, height, ceilDatadir):

    """
    Crude function to read in several days worth of ceilometer data that correspond to the model forecast time
    :param day:
    :param site:
    :param height:
    :param ceilDatadir:
    :return:
    """

    # read in ALL ceilometer data, without subsampling times to match the model
    bsc_obs_0 = FO.read_ceil_obs_all(day + dt.timedelta(hours=-24), {site: height}, ceilDatadir)  # forecast day
    bsc_obs_1 = FO.read_ceil_obs_all(day, {site: height}, ceilDatadir)  # main day
    bsc_obs_2 = FO.read_ceil_obs_all(day + dt.timedelta(hours=+24), {site: height},
                                     ceilDatadir)  # extra bit up to hr 36

    # merge ceilometer data together
    # ToDo - very crude function
    bsc_obs = merge_dicts(bsc_obs_0, bsc_obs_1, bsc_obs_2)


def merge_dicts(bsc_obs_0, bsc_obs_1, bsc_obs_2):

    """
    Really crude bsc_obs merger. dictionaries need to be ordered before going in
    :param dict_args:
    :return:
    """

    from copy import deepcopy

    bsc_obs = deepcopy(bsc_obs_0)
    dicts = [bsc_obs_1, bsc_obs_2]

    for d in dicts: # you can list as many input dicts as you want here

        for site, site_data in d.iteritems():

            bsc_obs[site]['backscatter'] = np.vstack((bsc_obs[site]['backscatter'], site_data['backscatter']))
            bsc_obs[site]['SNR'] = np.vstack((bsc_obs[site]['SNR'], site_data['SNR']))
            bsc_obs[site]['time'] += site_data['time']

    return bsc_obs

def trim_obs_time(bsc_obs, mod_data):

    """
    Crudely trim the bsc_obs data so the start and end times match the modelled data.
    NO subsampling is occuring
    :param bsc_obs:
    :param mod_data:
    :return:
    """

    for site in bsc_obs.iterkeys():

        # short site id that matches the model id
        site_id = site.split('_')[-1]

        # start and end idxs
        t_start = eu.nearest(bsc_obs[site]['time'], mod_data[site_id]['time'][0])[1]
        t_end = eu.nearest(bsc_obs[site]['time'], mod_data[site_id]['time'][-1])[1]

        # time range to extract data for
        t_idx = np.arange(t_start, t_end + 1)

        # extract data
        # for var, data in data_obs.iteritems():
        bsc_obs[site]['SNR'] = bsc_obs[site]['SNR'][t_idx, :]
        bsc_obs[site]['backscatter'] = bsc_obs[site]['backscatter'][t_idx, :]
        bsc_obs[site]['time'] = [bsc_obs[site]['time'][i] for i in t_idx]

    return bsc_obs

# process

def nearest_heights(mod_height, bsc_height):

    """
    Get nearest heights that are not unique for calculating the bias
    :param mod_height:
    :param bsc_height:
    :return:
    """

    # get nearest mod height level for each ceilometer gate
    a = np.array([eu.nearest(mod_height, i) for i in bsc_height])
    values = a[:, 0]
    diff = a[:, 2]

    mod_idx = np.array(a[:, 1], dtype=int)
    obs_idx = np.arange(len(bsc_height))  # mod_idx should be paired with obs_idx spots.

    # UNIQUE PAIRS NOT TAKEN

    # Remove pairs where mod is above 2000 m.
    # hc = height cut
    hc_pairs_range = np.where(values <= 2000)[0]

    # trim off unique pairs where the model level is above the maximum height
    obs_hc_pairs = obs_idx[hc_pairs_range]
    mod_hc_pairs = mod_idx[hc_pairs_range]
    pairs_hc_values = values[hc_pairs_range]
    pairs_hc_diff = diff[hc_pairs_range]

    return obs_hc_pairs, mod_hc_pairs

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

    # instruments and other settings
    site_bsc = FO.site_bsc
    site_bsc_colours = FO.site_bsc_colours

    # start and end
    dayStart = dt.datetime(2016, 05, 04)
    dayEnd = dt.datetime(2016, 05, 06)

    # ==============================================================================
    # Read data
    # ==============================================================================

    # 1. Read Ceilometer metadata
    # ----------------------------
    ceil_metadata = FO.read_ceil_metadata(datadir)

    # datetime range to iterate over
    # day the forecast started (21UTC)
    days_iterate = eu.date_range(dayStart, dayEnd, 1, 'days')

    for site, height in site_bsc.iteritems():

        # store each set of forecast statistics
        # store the day in a list so it can be recalled in order
        stat_bias = {}
        stat_bias_days = []

        # short site id that matches the model id
        site_id = site.split('_')[-1]

        for day in days_iterate:

            # string for the day
            dayStr = day.strftime('%Y%m%d')

            # 1. Read UKV forecast in
            # -----------------------

            # extract MURK aerosol and calculate RH for each of the sites in the ceil metadata
            # (can be different locations to sites_bsc)
            # reads all london model data, extracts site data, stores in single dictionary
            mod_data = FO.mod_site_extract_calc(day, {site_id: ceil_metadata[site_id]}, modDatadir, model_type, res, fullForecast=True)

            # 2. Read ceilometer
            # -----------------------

            bsc_obs = read_ceil_obs(day, site, height, ceilDatadir)

            # trim times
            # ToDo - very crude function
            bsc_obs = trim_obs_time(bsc_obs, mod_data)

            # 3. Empty variables to store bias
            # -------------------------------------
            stat_bias[dayStr] = {'time': [], 'height': [], 'bias': []}
            stat_bias_days += [dayStr]

            # ==============================================================================
            # Process
            # ==============================================================================

            # nearest heights
            obs_hc_pairs, mod_hc_pairs = nearest_heights(mod_data[site_id]['level_height'], bsc_obs[site['height']])

            # store heights
            stat_bias[dayStr]['bias'] = np.empty((bsc_obs[site]['backscatter'].shape[0],
                                                  obs_hc_pairs.shape[0]))

            for t in np.arange(len(bsc_obs[site]['time'])):

                # find t_idx in mod that match ceil
                # take difference for the whole profile
                # store it... somewhere

                # nearest time in mod to obs
                t_mod_idx = eu.nearest(mod_data[site_id]['time'], bsc_obs[site]['time'][t])[1]

                # extract profiles
                obs_x = bsc_obs[site]['backscatter'][t, obs_hc_pairs]
                mod_y = mod_data[site_id]['backscatter'][t_mod_idx, mod_hc_pairs]

                # bias (log10(obs) - log10(mod))
                stat_bias[dayStr]['bias'][t, :] = np.log10(obs_x) - np.log10(mod_y)

            # store the times
            stat_bias[dayStr]['time'] = deepcopy(bsc_obs[site]['time'])

            # store the heights
            stat_bias[dayStr]['height'] = deepcopy(bsc_obs[site]['height'][obs_hc_pairs])

            # convert bias profiles from list to numpy array
            # ToDo a more robust way of storing the data. It is lucky this creates the right shape tbh...
            # stat_bias[dayStr]['bias'] = np.array(stat_bias[dayStr]['bias'])

        # ==============================================================================
        # Plotting
        # ==============================================================================

        # start once the bia for, all the forecasts, for this site have been calculated
        # create a staggering bias plot
            #----------------
                    #----------------------
                                    #--------------------

        # setup plot
        fig, axs = plt.subplots(len(stat_bias_days), 1, figsize=(12, 7))

        # plot
        for ax, forecastDay in zip(axs.reshape(-1), stat_bias_days):

            # data for this forecast
            # extract using a list way, instead of looped over, so it plots in order
            data = stat_bias[forecastDay]

            mesh = ax.pcolormesh(data['time'], data['height'], np.transpose(data['bias']),
                          vmin=-2.0, vmax=2.0)

            # ax prettify
            ax.set_xlim([stat_bias[stat_bias_days[0]]['time'][0],
                         stat_bias[stat_bias_days[-1]]['time'][-1]])

            ax.set_ylim([np.nanmin(data['height']),
                         np.nanmax(data['height'])])

            ax.xaxis.set_major_formatter(DateFormatter('%d/ %H:%M'))

            if forecastDay != stat_bias_days[-1]:
                ax.axes.xaxis.set_ticklabels([])

        # figure prettify
        plt.tight_layout(h_pad=0.1)

        plt.subplots_adjust(left = 0.05, bottom = 0.05, right=1)

        # figure colourbar
        cax, kw = mpl.colorbar.make_axes(list(axs))
        cbar = plt.colorbar(mesh, cax=cax, **kw)


        ax0 = eu.fig_majorAxis(fig)
        ax0.set_title(site_id + ' log10(obs)-log10(mod): ' + str(stat_bias[stat_bias_days[0]]['time'][0]) + ' to ' +
                                        str(stat_bias[stat_bias_days[-1]]['time'][-1]))

        # save fig
        plt.savefig(savedir + 'dailyPlots/fullForecast/' + 'obs-mod-bias_' + site_id + '_' + model_type + '_' +
                    stat_bias[stat_bias_days[0]]['time'][0].strftime("%Y%m%d") + '-' +
                    stat_bias[stat_bias_days[-1]]['time'][-1].strftime("%Y%m%d") + '.png')





    print 'END PROGRAM'

    return


if __name__ == '__main__':
    main()