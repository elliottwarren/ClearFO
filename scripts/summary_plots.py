"""
Load in several days worth of data and make some simple summary plots, such as a distribution of radii, RH etc.

EW 21/02/2017
"""


import matplotlib.pyplot as plt

import numpy as np
import datetime as dt
from copy import deepcopy

import ellUtils as eu
from forward_operator import FOUtils as FO
from forward_operator import FOconstants as FOcon

def plot_r_m(mod_data, savedir):

    """
    Plot r_m frequency distribution as a line graph
    EW 22/02/17
    :param mod_data:
    :param savedir:
    :return:
    """

    fig = plt.figure(figsize=(8, 4.5))

    for site, site_data in mod_data.iteritems():

        r_m = site_data['r_m'] * 1.0e6 # micrometers
        # N = site_data['N'] * 1e-06  # micrometers

        # reference for dN/dlog(D) graphs (bottom of page 3)
        # http://www.tsi.com/uploadedFiles/_Site_Root/Products/Literature/Application_Notes
        # /PR-001-RevA_Aerosol-Statistics-AppNote.pdf

        # line style
        y, binEdges = np.histogram(r_m, bins=100)
        bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
        plt.plot(bincenters, y, '-', label=site)

        plt.xlabel(r'$radius \/\mathrm{[\mu m]}}$')
        # plt.xlabel(r'$number conc.\/\mathrm{[cm^{-1}]}}$')
        plt.ylabel('frequency')
        plt.grid()
        plt.legend(loc='best', fancybox=True, framealpha=0.5)
        plt.tight_layout()

    plt.savefig(savedir + 'r_m_dist' + '.png')
    plt.close(fig)

    return

def plot_N(mod_data, savedir):
    """
    Plot r_m frequency distribution as a line graph
    EW 22/02/17
    :param mod_data:
    :param savedir:
    :return:
    """

    fig = plt.figure(figsize=(8, 4.5))

    for site, site_data in mod_data.iteritems():

        N = site_data['N'] * 1e-06  # cm-1

        # reference for dN/dlog(D) graphs (bottom of page 3)
        # http://www.tsi.com/uploadedFiles/_Site_Root/Products/Literature/Application_Notes
        # /PR-001-RevA_Aerosol-Statistics-AppNote.pdf

        # line style
        y, binEdges = np.histogram(N, bins=100)
        bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
        plt.plot(bincenters, y, '-', label=site)

        plt.xlabel(r'$number conc.\/\mathrm{[cm^{-1}]}}$')
        plt.ylabel('frequency')
        plt.grid()
        plt.legend(loc='best', fancybox=True, framealpha=0.5)
        plt.tight_layout()

    plt.savefig(savedir + 'N_dist' + '.png')
    plt.close(fig)

    return

def plot_Q_line(mod_data, savedir):

    """
    Plot line histogram of Q
    EW 23/02/17
    :param mod_data:
    :param savedir:
    :return:
    """

    fig = plt.figure(figsize=(8, 4.5))

    for site, site_data in mod_data.iteritems():

        Q = site_data['Q'] # micrometers
        # N = site_data['N'] * 1e-06  # micrometers

        # reference for dN/dlog(D) graphs (bottom of page 3)
        # http://www.tsi.com/uploadedFiles/_Site_Root/Products/Literature/Application_Notes
        # /PR-001-RevA_Aerosol-Statistics-AppNote.pdf

        # line style
        y, binEdges = np.histogram(Q, bins=100)
        bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
        plt.plot(bincenters, y, '-', label=site)

        plt.xlabel('Q')
        # plt.xlim([0, 4])
        plt.ylabel('Frequency')
        plt.grid()
        plt.legend(loc='best', fancybox=True, framealpha=0.5)
        plt.tight_layout()

    plt.savefig(savedir + 'Q_dist' + '.png')
    plt.close(fig)

    return

def pcolor_Q(mod_data, site, ceil_lam, savedir):
    """
    Plot pcolor of Q for a site
    EW 23/02/17
    :param mod_data:
    :param savedir:
    :return:
    """

    fig = plt.figure(figsize=(8, 4.5))

    Q = mod_data[site]['Q']  # micrometers
    time = mod_data[site]['time']
    height = mod_data[site]['level_height']
    plt.pcolor(np.transpose(Q))
    plt.colorbar()
    plt.axis('tight')
    plt.xlabel('time index')
    plt.ylabel('height index')
    plt.tight_layout()

    plt.savefig(savedir + 'Q_pcolor_' + site + '_' + str(ceil_lam) + 'nm.png')
    plt.close(fig)

    return

def plot_dNdlogr(mod_data):

    """
    NOT YET FUNCTIONAL!!!!
    Create a dN/dlog(r) curve.
    :param mod_data:
    :return:
    """

    # try making the dN/dlog(D) graph
    N = mod_data['IMU']['N'] *1.0e-6
    r_m = mod_data['IMU']['r_m']*1.0e6
    logr_m = np.log10(r_m)

    N_mean = []
    N_cum =[]
    #concentration [micro meter-1 cm-1]
    conc = []
    idx_store =[]
    _, binEdges = np.histogram(r_m, bins=100)
    bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
    # binEdges[1:] - binEdges[:-1]
    dr = binEdges[1] - binEdges[0]


    for i in np.arange(len(binEdges[:-1])):
        idx = np.where((r_m >= binEdges[i]) & (r_m <= binEdges[i+1]))
        idx_store.append(idx)
        N_mean_i = np.nansum(N[idx])/np.sum(np.ones(N.shape)) #?????

        N_mean.append(N_mean_i)
        conc.append(N_mean_i/dr)
    #
    # N_mean = np.array(N_mean)
    # N_mean[np.isnan(N_mean)] = 0
    # N_cum = np.cumsum(N_mean)

    dN = np.array(N_cum[1:]) - np.array(N_cum[:-1])
    dlogr = binEdges[1:] - binEdges[:-1]
    # dN = np.array(N_sum[1:]) - np.array(N_sum[:-1])
    # plt.semilogx(N_sum / dlogr)
    plt.semilogx(dN / dlogr[:-1])


    return

def trim_to_lowest_heights(mod_data, max_height=2000):

    """
    Trim the data, such that only data below the [max_height] is kept. Default = 2000 m
    EW 22/02/17
    :param mod_data:
    :param max_height:
    :return:
    """

    # only keep data below 2000 m
    for site in mod_data.iterkeys():
        height_idx = np.where(mod_data[site]['level_height'] <= max_height)[0]

        for var in mod_data[site].iterkeys():
            if var == 'level_height':
                mod_data[site][var] = mod_data[site][var][height_idx]
            elif var != 'time':
                mod_data[site][var] = mod_data[site][var][:, height_idx]

    return mod_data

def main():

    # ==============================================================================
    # Setup
    # ==============================================================================

    # which modelled data to read in
    model_type = 'UKV'

    # ceilometer resolution [nm]
    ceil_lam = 910

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

    # append savdir to include the date range and mkdir if it doesn't exist
    savedir = savedir + days_iterate[0].strftime('%j') + '-' + days_iterate[-1].strftime('%j%Y') + '/'
    eu.ensure_dir(savedir)

    for day in days_iterate:

        # Read and concatonate data
        # ---------------------------

        # extract MURK aerosol and calculate RH for each of the sites in the ceil metadata
        # (can be different locations to sites_bsc)
        mod_day_data = FO.mod_site_extract_calc(day, ceil_metadata, modDatadir, model_type, res, ceil_lam, allvars=True)


        # concatonate data
        if day == days_iterate[0]:
            mod_data = deepcopy(mod_day_data)
        else:
            mod_data = eu.merge_dicts(mod_data, mod_day_data)

        #

    # trim data, such that only data below the [maximum_height] is kept
    mod_data = trim_to_lowest_heights(mod_data, max_height=2000)


    # ==============================================================================
    # Plotting
    # ==============================================================================

    # plotting
    plot_r_m(mod_data, savedir)
    plot_N(mod_data, savedir)
    plot_Q_line(mod_data, savedir)
    pcolor_Q(mod_data, 'IMU', ceil_lam, savedir)



    print 'END PROGRAM'

    plt.close('all')

    return

if __name__ == '__main__':
    main()