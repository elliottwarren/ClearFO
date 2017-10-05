"""
Plot sensible head flux (Q_H) at KSSW and several location+height combinations from the UKV

Created by Elliott Tues 29/06/17
"""

import matplotlib.pyplot as plt

import numpy as np
import datetime as dt
from scipy.stats import spearmanr
from matplotlib.dates import DateFormatter

from forward_operator import FOUtils as FO
from forward_operator import FOconstants as FOcon
import ellUtils as eu


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

    # # instruments and other settings
    #site_rh = FOcon.site_rh
    #site_rh = {'WXT_KSSW': 50.3}
    #rh_instrument = site_rh.keys()[0]

    #site = 'NK'
    #ceil_id = 'CL31-D'
    #ceil = ceil_id + '_' + site

    # instruments and other settings
    #site_rh = FOcon.site_rh
    #site_rh = {'Davis_IMU': 72.8}
    #rh_instrument = site_rh.keys()[0]

    site = 'MR'
    ceil_id = 'CL31-C'
    # ceil = ceil_id + '_BSC_' + site
    ceil = ceil_id + '_' + site

    site_bsc = {ceil: FOcon.site_bsc[ceil]}
    # site_bsc = {ceil: FOcon.site_bsc[ceil], 'CL31-E_BSC_NK': 27.0 - 23.2}

    site_bsc_colours = FOcon.site_bsc_colours

    # high pollution case study day
    daystrList = ['20160119']

    days_iterate = dateList_to_datetime(daystrList)

    # variable to compare (whatever the mod_site_extract_calc function has named them)
    variable = 'Q_H'


    # ==============================================================================
    # Read data
    # ==============================================================================

    # Read Ceilometer metadata

    # ceilometer list to use
    # ceilsitefile = 'UKV_correlatesites.csv'
    ceilsitefile = 'CeilsCSV_qgis_map.csv'
    ceil_metadata = FO.read_ceil_metadata(datadir, ceilsitefile)

    for day in days_iterate:

        print 'day = ' + day.strftime('%Y-%m-%d')

        # Read UKV forecast and automatically run the FO

        # extract MURK aerosol and calculate RH for each of the sites in the ceil metadata
        # reads all london model data, extracts site data, stores in single dictionary
        mod_data = FO.mod_site_extract_calc(day, ceil_metadata, modDatadir, model_type, res, 910, version=0.2, allvars=True)

        # read in CSAT Q_H obs
        datapath = datadir + 'L1/CSAT3_ECpack_KSSW_2016019_30min.nc'
        CSAT3 = eu.netCDF_read(datapath, vars='')

        time_obs = CSAT3['time']
        var_obs = CSAT3['Q_H']

        # idx heights to pull z
        zidx = 3

        # LINE PLOT
        fig, ax = plt.subplots(1, 1, figsize=(6, 2.5))

        # actual heights of UKV data
        for site in ceil_metadata.iterkeys():
            z_mod = mod_data[site]['level_height'][zidx]
            # z_kss45w = mod_data['KSS45W']['level_height'][zidx]

            # store these for correlating later
            var_mod = mod_data[site][variable][:, zidx]
            # var_mod_kss45w = mod_data['KSS45W'][variable][:, zidx]

            # store their times - same for all model areas so just use MR
            # time_mod = mod_data['MR']['time']
            time_mod = mod_data[site]['time']

            plt.plot_date(time_mod, var_mod, linestyle='-', marker='o', markersize=4, label=r'$Q_{H,'+site+'}$')
        # plt.plot_date(time_mod, var_mod_kss45w, linestyle='-', marker='o', markersize=4, label='$Q_{H,KSS45W}$')
        plt.plot_date(time_obs, var_obs, linestyle='-', marker='o', markersize=4, label=r'$Q_{H,o}$')
        plt.axhline(0, linestyle='--', alpha=0.5)

        ax.set_xlabel(r'$Time \/\/[HH:MM]$')
        ax.set_ylabel(r'$Q_{H} \/\/[W\/m^{-2}]$')
        ax.set_ylim([-10.0, 100.0])
        ax.set_xlim([day, day + dt.timedelta(days=1)])
        ax.xaxis.set_major_formatter(DateFormatter('%H:%M'))
        #plt.suptitle('Q_H z_mr=' + str(z_mr) + 'm.png')
        plt.legend(fontsize=10)
        plt.tight_layout()

        # fname = variable+'_MR' + '_' + str(z_mr) + 'm_'+ 'm_lineplot.png'
        fname = variable + '_5sites_lineplot.png'
        plt.savefig(savedir + 'point_diff/picking_heights/' + fname)

    plt.close('all')


    return

if __name__ == '__main__':
    main()























print 'END PROGRAM'
