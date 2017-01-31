"""
Load in several days worth of data and make some simple summary plots, such as a distribution of radii, RH etc.

"""


import matplotlib.pyplot as plt

import numpy as np
import datetime as dt

import ellUtils as eu
from forward_operator import FOUtils as FO
from forward_operator import FOconstants as FOcon


def main():

    # ==============================================================================
    # Setup
    # ==============================================================================

    # which modelled data to read in
    model_type = 'UKV'

    # model resolution
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

    # days to loop between
    dayStart = dt.datetime(2016, 05, 04)
    dayEnd = dt.datetime(2016, 05, 06)


    # ==============================================================================
    # Read and process modelled data
    # ==============================================================================

    # 1. Read Ceilometer metadata
    # ----------------------------
    ceil_metadata = FO.read_ceil_metadata(datadir)

    # 2. Read calibration data
    # ----------------------------


    # datetime range to iterate over
    days_iterate = eu.date_range(dayStart, dayEnd, 1, 'days')

    for day in days_iterate:


        # 2. Read UKV forecast in
        # -----------------------

        # extract MURK aerosol and calculate RH for each of the sites in the ceil metadata
        # (can be different locations to sites_bsc)
        mod_data = FO.mod_site_extract_calc(day, ceil_metadata, modDatadir, model_type, res)


        # ==============================================================================
        # Plotting
        # ==============================================================================



    print 'END PROGRAM'

    plt.close('all')

    return

if __name__ == '__main__':
    main()