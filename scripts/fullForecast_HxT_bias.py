"""
Plot the height x time bias profiles of the full forecast profiles together (3 profile sets)
"""

import matplotlib.pyplot as plt
import matplotlib as mpl

import datetime as dt

import ceilUtils as ceil
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

    # instruments and other settings
    site_bsc = FO.site_bsc
    site_bsc_colours = FO.site_bsc_colours

    # day start and end
    dayStart = dt.datetime(2016, 05, 04)
    dayEnd = dt.datetime(2016, 05, 06)




























    print 'END PROGRAM'

    return


if __name__ == '__main__':
    main()