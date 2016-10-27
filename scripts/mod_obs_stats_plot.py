"""

Script to do all the stats to the FO output. Correlations first...

"""

import matplotlib.pyplot as plt
from matplotlib.dates import date2num
from matplotlib.dates import DateFormatter

import numpy as np
import datetime as dt

import ceilUtils as ceil
import ellUtils as eu
import FOUtils as FO

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









    print 'END PROGRAM'

    plt.close('all')

    return

if __name__ == '__main__':
    main()