"""
4 panel plot of beta_m, beta_o, m and RH for a day

Created by Elliott 30/05/17
"""

import matplotlib.pyplot as plt

import numpy as np
import datetime as dt

from copy import deepcopy
import colorsys

import ellUtils as eu
import ceilUtils as ceil
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
    savedir = maindir + 'figures/' + model_type + '/clearSkyPeriod/'

    # data
    ceilDatadir = datadir + 'L1/'
    modDatadir = datadir + model_type + '/'
    rhDatadir = datadir + 'L1/'
    aerDatadir = datadir + 'LAQN/'

    # statistics to run
    pm10_stats = False
    rh_stats = True

    site = 'NK'
    ceil_id = 'CL31-D'
    ceil_id_full = ceil_id + '_' + site


    site_bsc = {ceil_id_full: FOcon.site_bsc[ceil_id_full]}
    # site_bsc = {ceil: FOcon.site_bsc[ceil], 'CL31-E_BSC_NK': 27.0 - 23.2}

    site_aer = {'PM10_'+site: FOcon.site_aer['PM10_'+site]}

    site_rh = {'WXT_KSSW': 50.3}
    rh_instrument = site_rh.keys()[0]

    site_bsc_colours = FOcon.site_bsc_colours

    # day = dt.datetime(2016, 01, 19) # new PM10 case study day
    day = dt.datetime(2016, 05, 04) # one of my old case study days

    # ==============================================================================
    # Read data
    # ==============================================================================

    # Read Ceilometer metadata

    # ceilometer list to use
    ceilsitefile = 'CeilsCSVfull.csv'
    ceil_metadata = FO.read_ceil_metadata(datadir, ceilsitefile)

    # extract out current site only
    ceil_data_i = {site: ceil_metadata[site]}

    print 'day = ' + day.strftime('%Y-%m-%d')

    # Read UKV forecast and automatically run the FO

    # extract MURK aerosol and calculate RH for each of the sites in the ceil metadata
    # reads all london model data, extracts site data, stores in single dictionary
    mod_data = FO.mod_site_extract_calc(day, ceil_data_i, modDatadir, model_type, res, 910)

    # Read ceilometer backscatter

    # will only read in data is the site is there!
    # ToDo Remove the time sampling part and put it into its own function further down.
    bsc_obs = FO.read_ceil_obs(day, site_bsc, ceilDatadir, mod_data, calib=True)


    # read in PM10 data and extract data for the current day
    pm10 = FO.read_pm10_obs(site_aer, aerDatadir, matchModSample=False)

    # extract the current day out of pm10
    # .date() from pm10 dates
    dates = np.array([i.date() for i in pm10['PM10_NK']['time']])
    idx = np.where(dates == day.date())

    # extract
    pm10['PM10_NK']['pm_10'] = pm10['PM10_NK']['pm_10'][idx]
    pm10['PM10_NK']['time'] = [pm10['PM10_NK']['time'][i] for i in idx]


    # read in RH data
    rh_obs = FO.read_all_rh_obs(day, site_rh, rhDatadir, mod_data)



    # plot the data
    # 4 panel, beta_o, beta_m, m with pm10 overlay, rh with rh_obs (KSSW) overlay
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(8, 4))

    site_id = site.split('_')[-1]

    # beta_o
    mesh1, ax1 = ceil.ceil_plot_to_ax(bsc_obs[site]['time'], bsc_obs[site]['height'], bsc_obs[site]['backscatter'],
                                      ax1, hmax=2000, vmin=1e-7, vmax=5e-6,
                                      tmin=day, tmax=day + dt.timedelta(days=1))

    # beta_m
    mesh2, ax2 = ceil.ceil_plot_to_ax(mod_data[site_id]['time'], mod_data[site_id]['level_height'], mod_data[site_id]['backscatter'],
                                      ax2, hmax=2000, vmin=1e-7, vmax=5e-6,
                                      tmin=day, tmax=day + dt.timedelta(days=1))

    # m
    mesh3, ax3 = ceil.ceil_plot_to_ax(mod_data[site_id]['time'], mod_data[site_id]['level_height'], mod_data[site_id]['aerosol_for_visibility'],
                                      ax3, hmax=2000,
                                      tmin=day, tmax=day + dt.timedelta(days=1))

    mesh4, ax4 = ceil.ceil_plot_to_ax(mod_data[site_id]['time'], mod_data[site_id]['level_height'], mod_data[site_id]['RH'],
                                      ax4, hmax=2000,
                                      tmin=day, tmax=day + dt.timedelta(days=1))



    print 'END PROGRAM'

    return



if __name__ == '__main__':
    main()