"""
Script to create sets of height x time profile plots. E.g. All 4 ceilometers for a day, all 4 \beta_m for a day etc.

Created by Fri Elliott 28th Oct 2016
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

def main():


    # ==============================================================================
    # Setup
    # ==============================================================================

    # which modelled data to read in
    model_type = 'UKV'
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

    # instruments and other settings
    site_bsc = FO.site_bsc
    site_rh = FO.site_rh
    site_aer = FO.site_aer
    site_bsc_colours = FO.site_bsc_colours

    # day start and end
    dayStart = dt.datetime(2016, 05, 04)
    dayEnd = dt.datetime(2016, 05, 04)


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

        # 2. Read ceilometer
        # -----------------------

        # read in ALL ceilometer data, without subsampling times to match the model
        bsc_obs = FO.read_ceil_obs_all(day, site_bsc, ceilDatadir, mod_data)

        # ==============================================================================
        # Plotting
        # ==============================================================================











    print 'END PROGRAM'

    return

if __name__ == '__main__':
    main()