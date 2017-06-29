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
import os

import ellUtils as eu
from forward_operator import FOUtils as FO
from forward_operator import FOconstants as FOcon

def dateList_to_datetime(dayList):

    """ Convert list of string dates into datetimes """

    datetimeDays = []

    for d in dayList:

        datetimeDays += [dt.datetime(int(d[0:4]), int(d[4:6]), int(d[6:8]))]

    return datetimeDays

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
    # ax.set_xlim([-7.5, -5.0]) # normal
    # ax.xaxis.set_ticks(np.arange(-7.5, -4.5, 0.5))

    ax.set_xlim([-7.0, -4.5]) # high PM cases
    ax.xaxis.set_ticks(np.arange(-7.0, -4.0, 0.5))

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
        ax.scatter(rh_obs[site]['RH'][t], height, label=site, color='green', edgecolors='black', s=120)

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

        ax.scatter(pm10_obs[site]['pm_10'][t], height, label=site, color='red', edgecolors='black', s=120)

    # prettify
    ax.set_ylim([0, 2000])
    # ax.set_xlim([0, 100]) # normal
    ax.set_xlim([0, 150]) # high PM case

    ax.set_yticklabels('')
    ax.set_xlabel(r'$m \/\mathrm{[\mu g\/ kg^{-1}]}$')
    eu.add_at(ax, 'c)', loc=1)

    return ax

def rh_ts_plot(ax, rh_obs, mod_data, site_bsc_colours, t_hr):

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
    ax.plot_date([t_hr, t_hr], [0, 100], color='black', ls='--', fmt='--')

    # prettify
    ax.set_ylabel(r'$\mathrm{RH\/\/[\%]}$')
    # ax.set_ylim([20, 80]) # normal
    ax.set_ylim([20, 100]) # normal
    ax.set_xlim([data_mod['time'][0], data_mod['time'][-1]])
    # ax.set_prop_cycle(cycler('color', ['b', 'c', 'm', 'y', 'k']))
    ax.xaxis.set_major_formatter(DateFormatter('%H:%M'))
    ax.set_xlabel('Time [HH:MM]')
    eu.add_at(ax, 'd)', loc=2)

    return ax

def aer_ts_plot(ax, pm10_obs, mod_data, site_bsc_colours, t_hr):


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
    ax.plot_date([t_hr, t_hr], [0, 1000], color='black', ls='--', fmt='-')

    # prettify
    ax.set_ylabel(r'$m \/\mathrm{[\mu g\/ kg^{-1}]}$')
    # ax.set_ylim([0, 100]) # normal
    ax.set_ylim([20, 160]) # high PM case
    ax.set_xlim([data_mod['time'][0], data_mod['time'][-1]])
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

    # ceilometer wavelength [nm]
    ceil_lam = 910

    # model resolution
    res = FOcon.model_resolution[model_type]

    # directories
    maindir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/'
    datadir = maindir + 'data/'
    savedir = maindir + 'figures/' + model_type + '/profiles/'

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

    # sites to run for

    # days to loop between - original clear sky case study
    # dayStart = dt.datetime(2016, 05, 04)
    # dayEnd = dt.datetime(2016, 05, 06)

    dayStart = dt.datetime(2016, 05, 04)
    dayEnd = dt.datetime(2016, 05, 06)

    # daystrList = ['20150414', '20150415', '20150421', '20150611', '20160504', '20160823', '20160911', '20161125',
    #               '20161129', '20161130', '20161204']

    # daystrList = ['20160119']

    daystrList =['20161204']

    days_iterate = dateList_to_datetime(daystrList)

    # ==============================================================================
    # Read and process modelled data
    # ==============================================================================

    # 1. Read Ceilometer metadata
    # ----------------------------
    ceil_metadata = FO.read_ceil_metadata(datadir, loc_filename='CeilsCSVclearFO.csv')


    # datetime range to iterate over
    # days_iterate = eu.date_range(dayStart, dayEnd, 1, 'days')

    for day in days_iterate:


        # 2. Read UKV forecast in
        # -----------------------

        # extract MURK aerosol and calculate RH for each of the sites in the ceil metadata
        # (can be different locations to sites_bsc)
        mod_data = FO.mod_site_extract_calc(day, ceil_metadata, modDatadir, model_type, res, ceil_lam)


        # 3. Read WXT and Davis
        # ------------------------

        # all RH obs for the main day, for all sites
        rh_obs = FO.read_all_rh_obs(day, site_rh, rhDatadir, mod_data)


        # 4. Read PM10
        # ------------------------

        pm10_obs = FO.read_pm10_obs(site_aer, aerDatadir, mod_data)



        # 5. Read ceilometer backscatter
        # --------------------------------

        bsc_obs = FO.read_ceil_obs(day, site_bsc, ceilDatadir, mod_data, calib=False)


        # only keep mod data where ceilometer data also exists
        bsc_sites_present = [i.split('_')[-1] for i in bsc_obs.keys()]
        keys = mod_data.keys()
        for key in keys:
            if (key in bsc_sites_present) == False:
                del mod_data[key]

        # ==============================================================================
        # Plotting
        # ==============================================================================

        for t in np.arange(len(mod_data.values()[0]['time'])):

            print t

            # current hour as string
            t_hr = mod_data.values()[0]['time'][t]

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
            ax4 = rh_ts_plot(ax4, rh_obs, mod_data, site_bsc_colours, t_hr)

            # PM10 time series
            ax5 = aer_ts_plot(ax5, pm10_obs, mod_data, site_bsc_colours, t_hr)

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
            if (os.path.exists(savedir + day.strftime("%Y%m%d"))) == False:
                os.mkdir(savedir + day.strftime("%Y%m%d"))
            plt.savefig(savedir + day.strftime("%Y%m%d") + '/BscRhAer_' + model_type + '_' + t_hr.strftime("%Y%m%d_%H%M") + '.png')

            # close figure
            plt.close(fig)

    print 'END PROGRAM'

    return

if __name__ == '__main__':
    main()