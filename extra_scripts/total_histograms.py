"""
Create histograms of UKV variables after loading in a sample of days

Created by Elliott 27/06/17
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

def dateList_to_datetime(dayList):

    """ Convert list of string dates into datetimes """

    datetimeDays = []

    for d in dayList:

        datetimeDays += [dt.datetime(int(d[0:4]), int(d[4:6]), int(d[6:8]))]

    return datetimeDays

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
    savedir = maindir + 'figures/' + model_type + '/highPmCase/'

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
    # dayStart = dt.datetime(2016, 05, 04)
    # dayEnd = dt.datetime(2016, 05, 06)

    daystrList = ['20150414', '20150415', '20150421', '20150611', '20160504', '20160823', '20160911', '20161125',
                  '20161129', '20161130', '20161204']

    days_iterate = dateList_to_datetime(daystrList)

    # statistics to run
    stats_corr = False
    stats_mbe = True

    # mbe ranges
    mbe_limit_step = 500
    mbe_limit_max = 2000

    # correlation max height
    corr_max_height = 2000

    Z='21'

    # ==============================================================================
    # Read data
    # ==============================================================================

    # Read Ceilometer metadata

    # ceilometer list to use
    ceilsitefile = 'CeilsCSVfull.csv'
    ceil_metadata = FO.read_ceil_metadata(datadir, ceilsitefile)

    for day in days_iterate:

        print 'day = ' + day.strftime('%Y-%m-%d')

        # Read UKV forecast and automatically run the FO
        # multiply m by a coeff to modify it
        m_coeff = np.ones((25, 70))
        m_layer_coeff = 1.0
        m_coeff[:, :5] = m_layer_coeff

        # extract MURK aerosol and calculate RH for each of the sites in the ceil metadata
        # reads all london model data, extracts site data, stores in single dictionary
        mod_data = FO.mod_site_extract_calc(day, ceil_metadata, modDatadir, model_type, res, 910,
                                            m_coeff=m_coeff, Z=Z, version=0.2)

        a=1











    return

if __name__ == '__main__':
    main()