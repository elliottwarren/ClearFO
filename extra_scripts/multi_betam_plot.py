"""
Quick plots of beta_m at several locations

Created by Elliott Thurs 1st September 2017
"""

import matplotlib.pyplot as plt

import numpy as np
import datetime as dt
from scipy.stats import spearmanr

from copy import deepcopy

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
    savedir = maindir + 'figures/' + model_type + '/clearSkyPeriod/'

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

    # true list
    daystrList = ['20150414', '20150415', '20150421', '20150611', '20160504', '20160823', '20160911', '20161125',
                  '20161129', '20161130', '20161204']

    # daystrList = ['20150415'] # KSS45W has MLH, MR and NK have MH.

    # daystrList = ['20160504', '20160505']

    days_iterate = dateList_to_datetime(daystrList)

    # ==============================================================================
    # Read data
    # ==============================================================================

    # Read Ceilometer metadata

    # ceilometer list to use
    ceilsitefile = 'CeilsCSVcalibrated.csv'
    # ceilsitefile = 'UKV_sitesToExtract.csv'
    ceil_metadata = FO.read_ceil_metadata(datadir, ceilsitefile)

    for day in days_iterate:

        print 'day = ' + day.strftime('%Y-%m-%d')

        # Read UKV forecast and automatically run the FO

        # extract MURK aerosol and calculate RH for each of the sites in the ceil metadata
        # reads all london model data, extracts site data, stores in single dictionary
        mod_data = FO.mod_site_extract_calc(day, ceil_metadata, modDatadir, model_type, res, 910, version=0.2)
























    return

if __name__ == '__main__':
    main()