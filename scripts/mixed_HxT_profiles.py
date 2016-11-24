"""
Script to create sets of height x time profile plots. E.g. All 4 ceilometers for a day, all 4 \beta_m for a day etc.

Created by Fri Elliott 28th Oct 2016
"""

# ToDo Need to sort out the ceil_plot script. Trim the data a bit first to reduce the plotting time.

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
    dayStart = dt.datetime(2015, 04, 20)
    dayEnd = dt.datetime(2015, 04, 21)


    # ==============================================================================
    # Read data
    # ==============================================================================

    # 1. Read Ceilometer metadata
    # ----------------------------
    ceil_metadata = FO.read_ceil_metadata(datadir)

    # datetime range to iterate over
    days_iterate = eu.date_range(dayStart, dayEnd, 1, 'days')

    for day in days_iterate:


        # 1. Read UKV forecast in
        # -----------------------

        # extract MURK aerosol and calculate RH for each of the sites in the ceil metadata
        # (can be different locations to sites_bsc)
        # reads all london model data, extracts site data, stores in single dictionary
        mod_data = FO.mod_site_extract_calc(day, ceil_metadata, modDatadir, model_type, res)

        # 2. Read ceilometer
        # -----------------------

        # read in ALL ceilometer data, without subsampling times to match the model
        bsc_obs = FO.read_ceil_obs_all(day, site_bsc, ceilDatadir)

        # ==============================================================================
        # Plotting
        # ==============================================================================

        # plot \beta_o and \beta_m for the same site in coupled plots.

        for site in bsc_obs.keys():

            # short site id that matches the model id
            site_id = site.split('_')[-1]

            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 4))

            # obs
            mesh1, ax1 = ceil.ceil_plot_to_ax(bsc_obs[site], ax1,
                            hmax=2000, vmin=1e-7, vmax=5e-6,
                            tmin=day, tmax=day+dt.timedelta(days=1))

            # model
            mod_data[site_id]['height'] = mod_data[site_id]['level_height']
            mesh2, ax2 = ceil.ceil_plot_to_ax(mod_data[site_id], ax2,
                            hmax=2000, vmin=1e-7, vmax=5e-6,
                            tmin=day, tmax=day + dt.timedelta(days=1))

            plt.tight_layout()
            plt.subplots_adjust(left=0.1, right=1)

            # figure colourbar
            cax, kw = mpl.colorbar.make_axes([ax1, ax2])
            plt.colorbar(mesh1, cax=cax, **kw)

            ax0 = eu.fig_majorAxis(fig)
            ax0.set_title(site_id + ' - ' + day.strftime("%d/%m/%Y"))
            ax0.set_ylabel('Height [m]')
            ax0.set_xlabel('Time [HH:MM]')

            # save fig
            plt.savefig(savedir + 'dailyPlots/' + 'obs-' + site_id + '_' + model_type + '-' + site_id
                        + '_' + day.strftime("%Y%m%d") + '.png')






    print 'END PROGRAM'

    return

if __name__ == '__main__':
    main()