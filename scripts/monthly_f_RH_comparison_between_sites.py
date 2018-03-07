"""
Compare the monthly climatology of f(RH) between sites

Created by Elliott Mon 05/03/2018
Based on monthly_f_RH_creation.py
"""

import numpy as np
import ellUtils as eu
import matplotlib.pyplot as plt
import pickle
from netCDF4 import Dataset
import datetime as dt

# Reading


if __name__ == '__main__':

    # ------------------------------------------
    # Setup
    # ------------------------------------------
    # site information
    # site_ins = {'site_short': 'NK', 'site_long': 'North Kensington',
    #             'ceil_lambda': 0.905e-06, 'land-type': 'urban'}
    site_ins = {'site_short':'Ch', 'site_long': 'Chilbolton',
                'ceil_lambda': 0.905e-06, 'land-type': 'rural'}
    # site_ins = {'site_short':'Ha', 'site_long': 'Harwell',
    #             'ceil_lambda': 0.905e-06, 'land-type': 'rural'}

    ceil_lambda_nm_str = str(site_ins['ceil_lambda'] * 1e9) + 'nm'

    # User set args
    # band that read_spec_bands() uses to find the correct band
    #! Manually set
    band = 1

    # saveF(RH)?
    saveFRH = True

    # -------------------------

    # directories
    savedir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/figures/Mie/monthly_f(RH)/'
    fRHdir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/data/Mie/monthly_f(RH)/'
    specdir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/data/Mie/' \
              'monthly_f(RH)/sp_885-925_r_files/'
    pickleloaddir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/data/Mie/pickle/'


    # variables to take from file (as listed within the file) with index from BLOCK = 0
    # NOTE: data MUST be in ascending index order
    aer_index = {'Ammonium Sulphate': 1, 'Generic NaCl': 2, 'Biogenic': 3, 'Aged fossil-fuel OC': 4,
                 'Ammonium nitrate': 5}
    aer_order = ['Ammonium Sulphate', 'Generic NaCl', 'Aged fossil-fuel OC', 'Ammonium nitrate']

    aer_particles_chem = {'Ammonium Sulphate': '(NH4)2SO4', 'Generic NaCl': 'NaCl', 'Aged fossil-fuel OC': 'CORG',
                          'Ammonium nitrate': 'NH4NO3', 'Soot': 'CBLK'}

    aer_particles = ['(NH4)2SO4', 'NH4NO3', 'NaCl', 'CORG', 'CBLK']

    month_colours = {'blue': [1, 2, 12], 'green': [3, 4, 5], 'red': [6, 7, 8], 'orange': [9, 10, 11]}
    month_ls = {'-': [1, 4, 7, 10], '--': [2, 5, 8, 11], '-.': [12, 3, 6, 9]}
    # Q type to use in calculating f(RH)
    Q_type = 'extinction'
    print 'Q_type = ' + Q_type

    # ---------------------------------------------------
    # Read, Process and save f(RH)
    # ---------------------------------------------------

    f_RH_data = {}

    # read in f(RH) for each site
    for site in ['NK', 'Ch', 'Ha']:

        path = fRHdir + 'monthly_f(RH)_' + site + '_905.0nm.nc'

        # eu.netCDF_info(path)

        rawdata = eu.netCDF_read(path)

        f_RH_data[site] = rawdata['f(RH) MURK']

        if site == 'NK':
            months = rawdata['months']
            RH = rawdata['Relative Humidity']
            radii_range_nm = rawdata['radii_range_nm']

    radii_range_micron = radii_range_nm*1e-3

    # ---------------------------------------------------
    # Plotting
    # ---------------------------------------------------



    # 3. plot f(RH) for MURK, with respect to RH, for large particle sizes (> 0.4 microns)
    # wrt_radii_radii_ratio
    for month_idx in range(12):

        date_i = dt.datetime(1900, month_idx+1, 1).strftime('%b')

        # find 0.11 ahead of time
        # rad_0p11_idx = np.where(radii_range_micron == 0.11)[0][0]

        # loop through an arbitrary list of RH values to plot
        for rad_val in [0.7, 0.11, 3.0]:

            fig = plt.figure(figsize=(6, 4))

            # find where rh_val is
            _, rad_idx, _ = eu.nearest(radii_range_micron, rad_val)
            # print rad_idx

            for site in f_RH_data.iterkeys():


                # data to plot
                f_RH_plot_data = np.squeeze(f_RH_data['Ch'][month_idx, rad_idx, :]) \
                                 / f_RH_data[site][month_idx, rad_idx, :]

                plt.plot(RH*100.0, f_RH_plot_data, label=site,
                         linestyle='-')
            # plt.axvline(0.11, color='grey', alpha=0.3, linestyle='--')

            plt.xlabel('RH [%]')
            plt.ylabel('f(RH)')
            plt.ylim([0.0, 2.0])
            plt.legend(loc='top left')
            plt.suptitle(date_i + ' radii = ' + str(rad_val))
            plt.savefig(savedir + 'compare/' + 'wrt_radii_ratio' + site + '_f(RH)_MURK_'+str(month_idx+1) + \
                        '_radii_'+str(rad_val)+'.png')
            plt.close(fig)

    # 4. RATIO plot f(RH) for MURK, with respect to RH, for large particle sizes (> 0.4 microns)
    # wrt_radii_monthly_ratio
    for rad_val in [0.07, 0.11, 0.14, 0.4, 1.0, 3.0]:
        fig = plt.figure(figsize=(6, 4))

        # loop through an arbitrary list of RH values to plot

        # find where rh_val is
        rad_idx = np.where(radii_range_micron == rad_val)[0][0]
        print rad_idx

        for month_idx, month_i in enumerate(range(1, 13)):

            # get colour and linestyle for month
            for key, data in month_colours.iteritems():
                if month_i in data:
                    colour_i = key

            for key, data in month_ls.iteritems():
                if month_i in data:
                    ls = key

            # month as 3 letter str
            date_i = dt.datetime(1900, month_i, 1).strftime('%b')


            # ratio of f(RH) for this month / average across all months
            f_RH_plot_data = np.squeeze(f_RH['MURK'][month_idx, rad_idx, :]) \
                             / np.mean(f_RH['MURK'][:, rad_idx, :], axis=0)

            plt.plot(RH_int*100.0, f_RH_plot_data, label=date_i,
                     linestyle=ls, color=colour_i)
            # plt.axvline(0.11, color='grey', alpha=0.3, linestyle='--')

        plt.xlabel('RH [%]')
        plt.ylabel('f(RH)')
        plt.ylim([0.8, 1.5])
        plt.legend(loc='top left', ncol=2)
        plt.suptitle('radii: '+ str(rad_val) + ' [microns]')
        plt.savefig(savedir + 'wrt_radii_monthly_ratio/' + 'wrt_radii_ratio_' + site_ins['site_short'] + '_f(RH)_MURK_'+'%.2f' % rad_val+'.png')
        plt.close(fig)

    # # plot the data at the end so it is all neat and together
    # fig = plt.figure(figsize=(6, 4))
    #
    # for key, value in data.iteritems():
    #
    #     plt.plot(RH*100, f_RH[key], label=key, linestyle='-')
    #
    # # plt.plot(value[:, 0], f_RH['average with Aitken Sulphate'], label='average with Aitken Sulphate', linestyle='-')
    #
    # # plot the MURK one
    # # plt.plot(value[:, 0], f_RH['average'], label='average without Aitken Sulphate', color='black')
    # plt.plot(RH*100, f_RH['MURK'], label='MURK', color='black')
    #
    # # plot soot as a constant until I get f(RH) for it
    # plt.plot([0.0, 100.0], [1.0, 1.0], label='Soot')
    #
    # plt.legend(fontsize=10, loc='best')
    # plt.tick_params(axis='both', labelsize=11)
    # plt.xlabel('RH [%]', labelpad=0, fontsize=11)
    # # plt.ylabel(Q_type + ' f(RH)')
    # plt.ylabel(r'$f_{ext,rh}$', fontsize=11)
    # plt.ylim([0.0, 8.0])
    # plt.xlim([0.0, 100.0])
    # # plt.title(file_name + ': ' + band_lam_range + ' band')
    #
    # plt.tight_layout() # moved tight_layout above... was originally after the save (06/04/17)
    # plt.savefig(savedir + file_name + '_' + Q_type[0:3] + '_f_RH_all_' + band_lam_range + '_salt8.0.png')

    print 'END PROGRRAM'