"""
Load in several days worth of data and make some simple summary plots, such as a distribution of radii, RH etc.

"""


import matplotlib.pyplot as plt

import numpy as np
import datetime as dt
from copy import deepcopy

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
    savedir = maindir + 'figures/' + model_type + '/summary/'

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

    # Read Ceilometer metadata
    ceil_metadata = FO.read_ceil_metadata(datadir)

    # datetime range to iterate over
    days_iterate = eu.date_range(dayStart, dayEnd, 1, 'days')

    for day in days_iterate:

        # Read and concatonate data
        # ---------------------------

        # extract MURK aerosol and calculate RH for each of the sites in the ceil metadata
        # (can be different locations to sites_bsc)
        mod_day_data = FO.mod_site_extract_calc(day, ceil_metadata, modDatadir, model_type, res, allvars=True)


        # concatonate data
        if day == days_iterate[0]:
            mod_data = deepcopy(mod_day_data)
        else:
            mod_data = eu.merge_dicts(mod_data, mod_day_data)

    # ==============================================================================
    # Plotting
    # ==============================================================================


    # r_m
    fig = plt.figure(figsize=(8, 4.5))

    for site, site_data in mod_data.iteritems():

        r_m = site_data['r_m'] * 1.0e6 # micrometers

        # line style
        y, binEdges = np.histogram(r_m, bins=50)
        bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
        plt.plot(bincenters, y, '-', label=site)

        plt.xlabel(r'$radius \/\mathrm{[\mu m]}]}$')
        plt.ylabel('frequency')
        plt.grid()
        plt.legend(loc='best', fancybox=True, framealpha=0.5)
        plt.tight_layout()
        plt.savefig(savedir + 'r_m_dist_linestyle' + '.png')
        plt.close(fig)

        # stack style
        n, bins, patches = plt.hist(r_m, 50, alpha=0.75, label=site, edgecolor='none')

    print 'END PROGRAM'

    plt.close('all')

    return

if __name__ == '__main__':
    main()