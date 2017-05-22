"""
Cross correlate the ceilometer observations with and without the calibration values

Created by Elliott Mon 22/05/2017
"""


import matplotlib.pyplot as plt
from matplotlib.dates import date2num
from matplotlib.dates import DateFormatter

import numpy as np
import datetime as dt
from scipy.stats import spearmanr

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