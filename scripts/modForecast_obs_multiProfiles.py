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

import ceilUtils as ceil
import ellUtils as eu
import FOUtils as FO


def main():

    # ==============================================================================
    # Setup
    # ==============================================================================

    maindir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/'
    datadir = maindir + 'data/'

    # data dirs
    ceilDataDir = datadir + 'L1/'
    rhDataDir = datadir + 'L1/'
    aerDataDir = datadir + 'LAQN/'

    # sites
    site_bsc = ['CL31-A_BSC_IMU', 'CL31-B_BSC_RGS', 'CL31-C_BSC_MR', 'CL31-D_BSC_NK']
    site_RH = ['WXT_KSSW', 'Davis_BCT', 'Davis_BFCL', 'Davis_BGH', 'Davis_IMU', 'Davis_IML']
    site_aer = ['PM10_RGS', 'PM10_MR', 'PM10_NK']
    # aer_start = len(site_bsc) + len(site_RH)
    sites = site_bsc + site_RH + site_aer

    # height above sea level (LUMA metadata)- height above ground (DTM) = height above surface
    # surface height taken from DTM (Grimmond and ...) - Digital Terrain Model (surface only, no buildings)
    # ToDo IMU needs updating!!!!
    ceil_height_offset = [79.0 - 14.7, 28.1 - 19.4, 32.0 - 27.5, 27.0 - 23.2]
    # ID string
    # IDs = ceil.createID(site_bsc)
    site = [i.split('_')[-1] for i in sites]
    sensor = [i.split('_')[0] for i in sites]
    # currently uses height above sea level from google earth height
    aer_heights = [23.0 - 19.4, 32.0 - 27.5, 26.0 - 23.2]

    # days to loop between
    dayStart = 125
    dayEnd = 127

    # year date for ceil data
    year = 2016

    # modelled data date
    date_mod = [2016, 05, 04]

    # MURK aerosol ceofficients to artificially alter the mass mixing ratio
    # aer_mod_coeff = [0.2, 0.5, 1.0]

    # which modelled data to read in
    model_type = 'UKV'

    # plot save
    savedir = maindir + 'figures/' + model_type + '/'


    # ==============================================================================
    # Read and process modelled data
    # ==============================================================================

    # read UKV forecast in



    # read WXT and Davis

    # read PM10




































    print 'END PROGRAM'

    plt.close('all')

if __name__ == '__main__':
    main()