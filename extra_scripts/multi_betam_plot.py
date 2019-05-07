"""
Quick plots of beta_m at several locations

Created by Elliott Thurs 1st September 2017
"""

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
from matplotlib.dates import DateFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable

import numpy as np
import datetime as dt

import ellUtils as eu
import ceilUtils as ceil
from forward_operator import FOUtils as FO
from forward_operator import FOconstants as FOcon

def dateList_to_datetime(dayList):

    """ Convert list of string dates into datetimes """

    datetimeDays = []

    for d in dayList:

        datetimeDays += [dt.datetime(int(d[0:4]), int(d[4:6]), int(d[6:8]))]

    return datetimeDays

if __name__ == '__main__':

    # ==============================================================================
    # Setup
    # ==============================================================================

    # which modelled data to read in
    model_type = 'UKV'
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

    # true list
    daystrList = ['20150414', '20150415', '20150421', '20150611', '20160504', '20160823', '20160911', '20161125',
                  '20161129', '20161130', '20161204']

    daystrList = ['20160119']

    days_iterate = dateList_to_datetime(daystrList)

    Z='21'

    # ==============================================================================
    # Read data
    # ==============================================================================

    # Read Ceilometer metadata

    # site list to use
    ceilsitefile = 'UKV_quickcompare.csv'
    ceil_metadata = FO.read_ceil_metadata(datadir, ceilsitefile)

    for day in days_iterate:

        print 'day = ' + day.strftime('%Y-%m-%d')

        # Read UKV forecast and automatically run the FO

        # extract MURK aerosol and calculate RH for each of the sites in the ceil metadata
        # reads all london model data, extracts site data, stores in single dictionary
        mod_data = FO.mod_site_extract_calc(day, ceil_metadata, modDatadir, model_type, res, 910, version=0.2)




        # plot

        fig, ax = plt.subplots(6, 1, figsize=(10, 6.5))

        # beta_o
        # mesh1 = ax1.pcolormesh(bsc_obs[ceil_id_full]['time'], bsc_obs[ceil_id_full]['height'], np.transpose(bsc_obs[ceil_id_full]['backscatter']),
        #                                   norm=LogNorm(vmin=1e-7, vmax=1e-5), cmap=cm.get_cmap('jet'))

        mesh = []

        for site, ax_i in zip(mod_data.iterkeys(), ax):

            data = mod_data[site]

            # beta_m
            mesh_i = ax_i.pcolormesh(data['time'], data['level_height'],
                                   np.transpose(data['backscatter']),
                                   norm=LogNorm(vmin=1e-7, vmax=1e-5), cmap=cm.get_cmap('jet'))

            # # RH
            # mesh_i = ax_i.pcolormesh(data['time'], data['level_height'],
            #                        np.transpose(data['RH']) * 100,
            #                        cmap=cm.get_cmap('Blues'), vmin=40.0, vmax=100.0)

            # m
            mesh_i = ax_i.pcolormesh(data['time'], data['level_height'],
                                   np.transpose(data['aerosol_for_visibility']),
                                   vmin=0, vmax=100, cmap=cm.get_cmap('OrRd'))  # vmin=0, vmax=100


            mesh.append(mesh_i)

            ax_i.xaxis.set_major_formatter(DateFormatter('%H:%M'))
            ax_i.yaxis.label.set_size(10)
            ax_i.xaxis.label.set_size(10)
            ax_i.set_xlim([day, day + dt.timedelta(days=1)])
            ax_i.set_ylim([0, 1000.0])
            eu.add_at(ax_i, r'$'+site+'$', loc=2)

            # remove time ticks for all but the bottom plot
            if ax_i != ax[-1]:
                ax_i.get_xaxis().set_ticks([])

        plt.subplots_adjust(right=0.8)
        # plt.suptitle(str(Z) + 'Z: m =' + str(int(m_layer_coeff*100))+' %')

        for mesh_i, ax_i in zip(mesh, ax):

            divider = make_axes_locatable(ax_i)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            plt.colorbar(mesh_i, cax=cax)



        fig.suptitle('m_MURK')
        ax0 = eu.fig_majorAxis(fig)
        ax0.set_xlabel('Time [HH:MM]', fontsize=10, labelpad=2)
        ax0.set_ylabel('Height [m]', fontsize=10, labelpad=10)

        plt.tight_layout(h_pad=0.1)

        plt.savefig(savedir + model_type + '-6sites-m_MURK_' + day.strftime('%Y%m%d') +
                    '_' + str(Z) + 'Z.png')  # filename

        plt.close(fig)



















