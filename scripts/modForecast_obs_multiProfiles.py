"""
Plot a profile each hour of the forward modelled backscatter and ceil obs. As well as the whole period of RH and PM10.
Improved upon the original version to take netCDF of UKV forecast.

Created by Elliott Tues 18 Oct 2016

"""


import matplotlib.pyplot as plt
from matplotlib.dates import date2num
from matplotlib.dates import DateFormatter

import numpy as np
import datetime as dt

import ellUtils as eu
from forward_operator import FOUtils as FO

# Plotting

def bsc_profile_plot(ax, mod_data, bsc_obs, site_bsc_colours, t):

    """
    plot the modelled and observed backscatter profiles

    :param ax:
    :param mod_data:
    :param bsc_obs:
    :param t:
    :return: ax
    """


    # ToDo - if the legend gets messy. Remove label from obs and just use the one from mod. Then set the colour code
    # to do a set rotation. As dictionary can be plotted in any order, the mod site will need to 'get' the colour
    # of the obs at the same site, to know they match.


    # plot obs
    for site_obs, data_obs in bsc_obs.iteritems():

        # find colour for the site
        colour = site_bsc_colours[site_obs.split('_')[-1]]

        ax.plot(np.log10(data_obs['backscatter'][t, :]), data_obs['height'],
                color=colour, ls='-')


    # plot model
    for site_mod, data_mod in mod_data.iteritems():

        # get matching colour for the site
        colour = site_bsc_colours[site_mod]

        ax.plot(np.log10(data_mod['backscatter'][t, :]), data_mod['level_height'],
                color=colour, label=site_mod, linewidth=2, ls='--')


    # Prettify figure at the end
    ax.set_xlabel(r'$log_{10}(\beta) \/\/\mathrm{[m^{-1} \/sr^{-1}]}$')
    ax.set_ylabel('Height [m]')

    ax.set_ylim([0, 2000])
    ax.set_xlim([-7, -5])
    ax.xaxis.set_ticks(np.arange(-7, -4.5, 0.5))
    # ax.set_prop_cycle(cycler('color', ['b', 'c', 'm', 'y', 'k']))
    eu.add_at(ax, 'a)', loc=2)

    return ax

def rh_profile_plot(ax, site_rh, rh_obs, mod_data, site_bsc_colours, t):

    """
    Plot the RH profile and obs

    :param ax:
    :param site_rh:
    :param rh_obs:
    :param mod_data:
    :param t:
    :return: ax
    """

    # plot model rh
    for site_mod, data_mod in mod_data.iteritems():

        # matching colour for the site with corresponding \beta_m estimates
        colour = site_bsc_colours[site_mod]

        ax.plot(data_mod['RH'][t, :] * 100, data_mod['level_height'], color=colour, linewidth=1, ls='--')

    # plot obs rh ontop of lines
    for site, height in site_rh.iteritems():

        # currently no label, so it can be added manually later
        ax.scatter(rh_obs[site]['RH'][t], height, label=site, color='green', edgecolors='black')

    # plot reference line
    ax.plot([38.0, 38.0], [0, 10000], color='black', ls='--')

    # prettify
    ax.set_ylim([0, 2000])
    ax.set_xlim([20, 100])
    ax.set_xlabel(r'$\mathrm{RH\/\/[\%]}$')
    ax.set_yticklabels('')
    # ax.set_prop_cycle(cycler('color', ['b', 'c', 'm', 'y', 'k']))
    eu.add_at(ax, 'b)', loc=2)

    return ax

def aer_profile_plot(ax, site_aer, pm10_obs, mod_data, site_bsc_colours, t):

    """
    Plot the pm10 obs from LAQN and the modelled murk aerosol

    :param ax3:
    :param site_aer:
    :param pm10_obs:
    :param mod_data:
    :param t:
    :return: ax
    """

    # plot model pm10
    for site_mod, data_mod in mod_data.iteritems():

        # get matching colour for the site
        colour = site_bsc_colours[site_mod]

        ax.plot(data_mod['aerosol_for_visibility'][t, :], data_mod['level_height'],
                color=colour, linewidth=2, ls='--')

    # plot obs pm10 ontop of lines
    for site, height in site_aer.iteritems():

        ax.scatter(pm10_obs[site]['pm_10'][t], height, label=site, color='red', edgecolors='black')

    # prettify
    ax.set_ylim([0, 2000])
    ax.set_xlim([0, 100])
    ax.set_yticklabels('')
    # ax.set_xlabel(r'$m$' + ' [' + r'$\mu$' + 'g kg' + r'$-1$' +']')
    ax.set_xlabel(r'$m \/\mathrm{[\mu g\/ kg^{-1}]}$')
    # ax.set_prop_cycle(cycler('color', ['b', 'c', 'm', 'y', 'k']))
    eu.add_at(ax, 'c)', loc=1)

    return ax

def rh_ts_plot(ax, rh_obs, mod_data, site_bsc_colours, t):

    """
    Plot the rh time series for the period. Uses average rh model input calculated
    outside this function and the unchanged RH observations.

    :param ax4:
    :param site_rh:
    :param rh_obs:
    :param mod_data:
    :param t:
    :return: ax
    """


    # plot rh obs profile plots
    for site_obs, data_obs in rh_obs.iteritems():

        ax.plot_date(date2num(data_obs['time']),data_obs['RH'],
                     label=site_obs, fmt='-')

    # plot model for a single height near the surface (5 m)
    # height index to plot
    h_idx = 0

    for site_mod, data_mod in mod_data.iteritems():

        # get matching colour for the site
        colour = site_bsc_colours[site_mod]

        ax.plot_date(date2num(data_mod['time']), data_mod['RH'][:, h_idx] * 100,
                     color=colour, label=site_mod, linewidth=2, fmt='-', ls='--')

    # plot reference line to show where profile lies
    ax.plot_date([t, t], [0, 100], color='black', ls='--', fmt='--')

    # prettify
    ax.set_ylabel(r'$\mathrm{RH\/\/[\%]}$')
    ax.set_ylim([20, 80])
    ax.set_xlim([data_mod['time'][0], data_mod['time'][-1]])
    # ax.set_prop_cycle(cycler('color', ['b', 'c', 'm', 'y', 'k']))
    ax.xaxis.set_major_formatter(DateFormatter('%H:%M'))
    ax.set_xlabel('Time [HH:MM]')
    eu.add_at(ax, 'd)', loc=2)

    return ax

def aer_ts_plot(ax, pm10_obs, mod_data, site_bsc_colours, t):


    # plot pm10 obs profile plots

    for site_obs, data_obs in pm10_obs.iteritems():
        ax.plot_date(date2num(data_obs['time']),data_obs['pm_10'], label=site_obs, fmt='-')

    # plot model for a single height near the surface (5 m)
    # height index to plot
    h_idx = 0

    for site_mod, data_mod in mod_data.iteritems():

        colour = site_bsc_colours[site_mod]

        ax.plot_date(date2num(data_mod['time']), data_mod['aerosol_for_visibility'][:, h_idx],
                     color=colour, label=site_mod, linewidth=2, fmt='-', ls='--')

    # plot reference line to show where profile lies
    ax.plot_date([t, t], [0, 100], color='black', ls='--', fmt='-')

    # prettify
    ax.set_ylabel(r'$m \/\mathrm{[\mu g\/ kg^{-1}]}$')
    ax.set_ylim([0, 100])
    ax.set_xlim([data_mod['time'][0], data_mod['time'][-1]])
    # ax.set_prop_cycle(cycler('color', ['b', 'c', 'm', 'y', 'k']))
    ax.xaxis.set_major_formatter(DateFormatter('%H:%M'))
    ax.set_xlabel('Time [HH:MM]')
    eu.add_at(ax, 'e)', loc=2)

    return ax



def main():

    # ==============================================================================
    # Setup
    # ==============================================================================

    # which modelled data to read in
    model_type = 'UKV'

    # model resolution
    res = FO.model_resolution[model_type]

    # directories
    maindir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/'
    datadir = maindir + 'data/'
    savedir = maindir + 'figures/' + model_type + '/'

    # data
    ceilDatadir = datadir + 'L1/'
    modDatadir = datadir + model_type + '/'
    rhDatadir = datadir + 'L1/'
    aerDatadir = datadir + 'LAQN/'

    # sites with height pairs
    # height above sea level (LUMA metadata)- height above ground (DTM) = height above surface
    # surface height taken from DTM (Grimmond and ...) - Digital Terrain Model (surface only, no buildings)
    # ToDo IMU needs updating!!!!
    # ToDo RH heights needs filling in too!
    site_bsc = {'CL31-A_BSC_IMU': 79.0 - 14.7, 'CL31-B_BSC_RGS': 28.1 - 19.4, 'CL31-C_BSC_MR': 32.0 - 27.5,
                'CL31-D_BSC_NK': 27.0 - 23.2}
    site_rh = {'WXT_KSSW': 0, 'Davis_BCT': 0, 'Davis_BFCL': 0, 'Davis_BGH': 0, 'Davis_IMU': 0, 'Davis_IML': 0}
    site_aer = {'PM10_RGS': 23.0 - 19.4, 'PM10_MR': 32.0 - 27.5, 'PM10_NK': 26.0 - 23.2}
    # sites = site_bsc + site_rh + site_aer

    # colours for plotting
    site_bsc_colours = {'IMU': 'b', 'RGS': 'y', 'MR': 'r', 'NK': 'm', 'SW': 'k', 'KSSW': 'c'}
    # ax.set_prop_cycle(cycler('color', ['b', 'c', 'm', 'y', 'k']))

    # ceil_height_offset = [79.0 - 14.7, 28.1 - 19.4, 32.0 - 27.5, 27.0 - 23.2]
    aer_heights = [23.0 - 19.4, 32.0 - 27.5, 26.0 - 23.2]
    rh_heights = [0, 0, 0, 0, 0, 0]

    # days to loop between
    # dayStart = dt.datetime(2016, 05, 03)
    # dayEnd = dt.datetime(2016, 05, 05)

    dayStart = dt.datetime(2016, 05, 04)
    dayEnd = dt.datetime(2016, 05, 06)


    # ==============================================================================
    # Read and process modelled data
    # ==============================================================================

    # 1. Read Ceilometer metadata
    # ----------------------------
    ceil_metadata = FO.read_ceil_metadata(datadir)


    # datetime range to iterate over
    days_iterate = eu.date_range(dayStart, dayEnd, 1, 'days')

    for day in days_iterate:


        # 2. Read UKV forecast in
        # -----------------------

        # extract MURK aerosol and calculate RH for each of the sites in the ceil metadata
        # (can be different locations to sites_bsc)
        mod_data = FO.mod_site_extract_calc(day, ceil_metadata, modDatadir, model_type, res)


        # 3. Read WXT and Davis
        # ------------------------

        # all RH obs for the main day, for all sites
        rh_obs = FO.read_all_rh_obs(day, site_rh, rhDatadir, mod_data)


        # 4. Read PM10
        # ------------------------

        pm10_obs = FO.read_all_pm10_obs(dayStart, dayEnd, site_aer, aerDatadir, mod_data)



        # 5. Read ceilometer backscatter
        # --------------------------------

        bsc_obs = FO.read_ceil_obs(day, site_bsc, ceilDatadir, mod_data)


        # do any stats on time reduced data here...


        # ==============================================================================
        # Plotting
        # ==============================================================================

        for t in np.arange(len(bsc_obs['CL31-A_BSC_IMU']['time'])):

            # current hour as string
            t_hr = bsc_obs['CL31-A_BSC_IMU']['time'][t]

            # set up plot
            fig = plt.figure(figsize=(10, 7))

            # set up subplots
            ax1 = plt.subplot2grid((3, 3), (0, 0))  # bsc
            ax2 = plt.subplot2grid((3, 3), (0, 1))  # aer_mmr
            ax3 = plt.subplot2grid((3, 3), (0, 2))  # rh
            ax4 = plt.subplot2grid((3, 3), (1, 0), colspan=3)  # rh time series
            ax5 = plt.subplot2grid((3, 3), (2, 0), colspan=3)  # aer + pm10 time series

            # ceilometer profile
            ax1 = bsc_profile_plot(ax1, mod_data, bsc_obs, site_bsc_colours, t)

            # RH profile plot
            ax2 = rh_profile_plot(ax2, site_rh, rh_obs, mod_data, site_bsc_colours, t)

            # PM10 profile plot
            ax3 = aer_profile_plot(ax3, site_aer, pm10_obs, mod_data, site_bsc_colours, t)

            # RH time series
            ax4 = rh_ts_plot(ax4, rh_obs, mod_data, site_bsc_colours, t)

            # PM10 time series
            ax5 = aer_ts_plot(ax5, pm10_obs, mod_data, site_bsc_colours, t)

            # final prettifying for figure
            fig.suptitle('BSC: ' + t_hr.strftime("%Y-%m-%d %H:%M:%S"), fontsize=12)
            plt.tight_layout(h_pad=0)

            # adjust axis so the ledgends can be placed on the side
            fig.subplots_adjust(top=0.95, right=0.8)

            # legends
            han1, lab1 = ax1.get_legend_handles_labels()
            han2, lab2 = ax2.get_legend_handles_labels()
            han3, lab3 = ax3.get_legend_handles_labels()

            lab2 = ['RH observations']
            lab3 = ['PM' + r'$_{10}$' + ' observations']
            han2 = [han2[2]]
            han3 = [han3[2]]  # ignore the other points there

            ax3.legend(han1 + han2 + han3, lab1 + lab2 + lab3, fontsize=8, bbox_to_anchor=(1.07, 1), loc=2,
                       borderaxespad=0.0)
            ax4.legend(fontsize=8, bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.0)
            ax5.legend(fontsize=8, bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.0)

            # save the figure
            plt.savefig(savedir + 'profiles/' + 'BscRhAer_' + model_type + '_' + t_hr.strftime("%Y%m%d_%H%M") + '.png')

            # close figure
            plt.close(fig)

    print 'END PROGRAM'

    plt.close('all')

    return

if __name__ == '__main__':
    main()