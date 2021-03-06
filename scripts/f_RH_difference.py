"""
Create 4 panel plot with f(RH) differences for each of the aerosol species and the average f(RH)

! WARNING - currently fixed to do 910, including save file names. Needs updating

Created by Elliott 06/04/17
"""

import sys
sys.path.append('C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/scripts')

import numpy as np
import matplotlib.pyplot as plt
from ellUtils import fig_majorAxis
from f_RH_creation import read_aer_data, read_spec_bands, calc_f_RH
from matplotlib.ticker import FormatStrFormatter

def create_f_RH(file_paths, band, aer_index, aer_order, Q_type):

    """
    Read in data and create the f_RH data for each aerosol, as well as an average.

    :param f_RH_file:
    :param band:
    :param aer_index:
    :param aer_order:
    :return:

    Designed to only do this for one file at a time, in order to split the main file from the others.
    """

    # type saftey stuff
    if type(file_paths) == str:
        file_paths = [file_paths]
    elif type(file_paths) != list:
        raise TypeError('file_paths not given as str or list!')

    f_RH = {}

    for f_RH_file in file_paths:

        # read in the spectral band information
        spec_bands = read_spec_bands(f_RH_file)

        # wavelength range in current band
        band_idx = np.where(spec_bands['band'] == band)[0][0]
        band_lam_range = '%.0f' % (spec_bands['lower_limit'][band_idx] * 1.0e9) + '-' + \
                         '%.0f' % (spec_bands['upper_limit'][band_idx] * 1.0e9) + 'nm'

        # read the aerosol data
        data = read_aer_data(f_RH_file, aer_index, aer_order, band=band)

        # Extract RH for plotting (RH is the same across all aerosol types)
        # convert from [frac] to [%]
        RH = np.array(data[aer_order[0]][:, 0]) * 100.0

        # calculate f(RH)
        # define wavelength key using band_lam_range
        Q, f_RH[band_lam_range] = calc_f_RH(data, aer_order, Q_type=Q_type)

        # # create an average f(RH)
        # f_RH[band_lam_range]['MURK'] = np.mean(f_RH[band_lam_range].values(), axis=0)

        # create f(RH) for MURK
        f_RH[band_lam_range]['MURK'] = np.sum([0.295 * np.array(f_RH[band_lam_range]['Accum. Sulphate']),
                                              0.38 * np.array(f_RH[band_lam_range]['Aged fossil-fuel OC']),
                                              0.325 * np.array(f_RH[band_lam_range]['Ammonium nitrate'])], axis=0)

        # f_RH['average'] = np.mean((f_RH['Accum. Sulphate'], f_RH['Aged fossil-fuel OC'], f_RH['Ammonium nitrate']), axis=0)

    return f_RH, band_lam_range, RH

def calc_r_RH_difference(f_RH, f_RH_main, main_band_range):

    """
    Calculate the difference in f(RH) between the various wavelengths (f_RH) and the reference (f_RH_main)

    :param f_RH:
    :param f_RH_main:
    :param main_band_range:
    :return:
    """

    f_RH_diff = {}

    for lam_range, all_aerosol in f_RH.iteritems():

        f_RH_diff[lam_range] = {}

        for aerosol_key, aerosol_data in all_aerosol.iteritems():

            f_RH_diff[lam_range][aerosol_key] = np.array(aerosol_data) - np.array(f_RH_main[main_band_range][aerosol_key])

    return f_RH_diff

def calc_r_RH_ratio(f_RH, f_RH_main, main_band_range):
    """
    Calculate the ratio of f(RH) between the various wavelengths (f_RH) and the reference (f_RH_main)

    :param f_RH:
    :param f_RH_main:
    :param main_band_range:
    :return:
    """

    f_RH_ratio = {}

    for lam_range, all_aerosol in f_RH.iteritems():

        f_RH_ratio[lam_range] = {}

        for aerosol_key, aerosol_data in all_aerosol.iteritems():

            f_RH_ratio[lam_range][aerosol_key] = np.array(aerosol_data) / np.array(
                f_RH_main[main_band_range][aerosol_key])

    return f_RH_ratio

# plotting

def plot_absolute(f_RH_main, f_RH, RH, savedir, aer_order, aer_names, lam_colours, Q_type):

    main_lam = f_RH_main.keys()[0]

    # plot the data on a 4 panel plot
    fig, ax = plt.subplots(2, 2, figsize=(8, 6))

    # append MURK to the aer_order to create a plotting order and aer_names so it's name can be plotted
    aer_plot_order = aer_order + ['MURK']
    aer_names['MURK'] = 'MURK'

    # plot on each subplot fully, before doing the next one
    for ax_i, aerosol in zip(ax.flatten(), aer_plot_order):

        # plot all curves for a single aerosol, at a time
        for lam_range in f_RH.iterkeys():

            # ax_i.plot(RH, f_RH[lam_range][aerosol], color=lam_colours[lam_range], linestyle='--', label=lam_range)
            ax_i.plot(RH, f_RH[lam_range][aerosol], linestyle='--')

        # plot 910 nm
        # ax_i.plot(RH, f_RH_main[main_lam][aerosol], color='black', linestyle='--', label=main_lam)



        # prettify subplot
        ax_i.set_title(aer_names[aerosol])
        # ax_i.set_ylim([-0.2, 0.2])
        # ax_i.set_ylim([-1.0, 1.0])
        ax_i.set_xlim([0.0, 100.0])
        ax_i.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

    # prettify figure
    ax_main = fig_majorAxis(fig)
    ax_main.set_xlabel('RH [%]')
    plt.tight_layout()
    plt.subplots_adjust(top=0.9, right=0.8)
    ax.flatten()[1].legend(fontsize=8, bbox_to_anchor=(1.07, 1), loc=2, borderaxespad=0.0)
    fig.suptitle('Absolute value (lam_i)')
    plt.savefig(savedir + 'absolute/multi-panel-ceils_absolute_905_' + Q_type[0:3] + '_f_RH.png')
    plt.close()

    return fig

def plot_difference(f_RH_diff, RH, savedir, aer_order, aer_names, lam_colours, Q_type):

    """
    Plot the ratio of the f(RH) references compared to the main f(RH) curve

    :param f_RH_diff:
    :param RH:
    :param savedir:
    :param aer_order:
    :param lam_colours:
    :param Q_type:
    :return:
    """

    # plot the data on a 4 panel plot
    fig, ax = plt.subplots(2, 2, figsize=(8, 6))

    # append MURK to the aer_order to create a plotting order and aer_names so it's name can be plotted
    aer_plot_order = aer_order + ['MURK']
    aer_names['MURK'] = 'MURK'

    # plot on each subplot fully, before doing the next one
    for ax_i, aerosol in zip(ax.flatten(), aer_plot_order):

        # plot all curves for a single aerosol, at a time
        for lam_range in f_RH_diff.iterkeys():

            ax_i.plot(RH, f_RH_diff[lam_range][aerosol])

        # # plot all curves for a single aerosol, at a time
        # for lam_range, lam_i_colour in zip(band_order, lam_colours):
        #
        #     ax_i.plot(RH, f_RH_diff[lam_range][aerosol], color=lam_i_colour, label=lam_range)


        # prettify subplot
        ax_i.set_title(aer_names[aerosol])
        # ax_i.legend()
        # ax_i.set_ylim([-0.2, 0.2])
        ax_i.set_xlim([0.0, 100.0])

    # prettify figure
    ax_main = fig_majorAxis(fig)
    ax_main.set_xlabel('RH [%]')
    plt.tight_layout()
    plt.subplots_adjust(top=0.9, right=0.8)
    ax.flatten()[1].legend(fontsize=8, bbox_to_anchor=(1.07, 1), loc=2, borderaxespad=0.0)
    # ax_main.set_ylabel('difference \n lam_i - 910nm')

    fig.suptitle('Absolute difference (lam_i - 910nm)')
    plt.savefig(savedir + 'difference/multi-panel-ceils_diff_905_' + Q_type[0:3] + '_f_RH.png')
    plt.close()

    return fig

def plot_ratio(f_RH_ratio, RH, savedir, aer_order, aer_names, lam_colours, Q_type):

    """
    Plot the ratio of the f(RH) references compared to the main f(RH) curve

    :param f_RH_diff:
    :param RH:
    :param savedir:
    :param aer_order:
    :param lam_colours:
    :param Q_type:
    :return:
    """

    # plot the data on a 4 panel plot
    fig, ax = plt.subplots(1, 1, figsize=(5.5, 3.5))

    # append MURK to the aer_order to create a plotting order and aer_names so it's name can be plotted
    aer_plot_order = aer_order + ['MURK']
    aer_names['MURK'] = 'MURK'

    # plot on each subplot fully, before doing the next one
    for aerosol in aer_plot_order:

        # plot between ceilometers
        # plot all curves for a single aerosol, at a time
        for lam_range in f_RH_ratio.iterkeys():

            ax.plot(RH, f_RH_ratio[lam_range][aerosol], label=aerosol)

    # prettify subplot
    # ax.set_title(aer_names[aerosol])
    ax.set_ylim([0.9, 1.3])
    ax.set_xlim([0.0, 100.0])
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

    # prettify figure
    ax_main = fig_majorAxis(fig)
    ax_main.set_xlabel('RH [%]')
    # ax_main.set_ylabel(r'$f_{ext},\/\lambda\/=\/i$' + '\n' + r'$f_{ext},\/\lambda\/=\/910nm$')
    ax_main.set_ylabel(r'$\frac{f_{ext,rh}(\lambda=1064\/nm)}{f_{ext,rh}(\lambda=905\/nm)}$', labelpad=15, fontsize=15)
    plt.tight_layout()
    # plt.subplots_adjust(top=0.9, right=0.8)
    ax.legend(fontsize=9)
    # ax.flatten()[1].legend(fontsize=8, bbox_to_anchor=(1.07, 1), loc=2, borderaxespad=0.0)
    # fig.suptitle('Ratio (lam_i / 910nm)')
    plt.savefig(savedir + 'ratio/ceil_905vs1064_' + Q_type[0:3] + '_f_RH.png')
    plt.close()

    return fig


if __name__ == '__main__':

    # User set args
    # band that read_spec_bands() uses to find the correct band
    #! Manually set
    band = 1

    # -------------------------

    # directories
    savedir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/figures/Mie/f(RH)/'
    specdir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/data/Mie/spectral/'
    f_RHdir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/data/Mie/'

    # file_name = 'spec3a_sw_hadgem1_7lean_so' # original file given to me by Claire Ryder 25/01/17
    # file_name = 'sp_sw_ga7' # current UM file
    # file_name = 'sp_ew_910' # my own made file with 1 band at 910 nm

    main_file = 'sp_ew_ceil_905'
    main_path = specdir + main_file
    main_file_band = 1 # 905 - 906

    file_names = ['sp_ew_ceil_1063.95-1064.05']

    # file_names = ['sp_ew_ceil_guass_903-907',
    #               'sp_ew_ceil_guass_1062-1066']

    file_paths = [specdir + f for f in file_names]

    # index: variables to take from file (as listed within the file) with index from BLOCK = 0
    #   NOTE: data MUST be in ascending index order
    # order: plotting order
    # names: titles to apply to each subplot
    # aer_index = {'Accum. Sulphate': 1, 'Aged fossil-fuel OC': 2, 'Ammonium nitrate': 3}
    # aer_order = ['Accum. Sulphate', 'Aged fossil-fuel OC', 'Ammonium nitrate']
    # aer_names = {'Accum. Sulphate': 'Ammonium sulphate', 'Aged fossil-fuel OC': 'Organic carbon', 'Ammonium nitrate': 'Ammonium nitrate'}

    aer_index = {'Ammonium Sulphate': 1, 'Aged fossil-fuel OC': 4, 'Ammonium nitrate': 5}
    aer_order = ['Ammonium Sulphate', 'Aged fossil-fuel OC', 'Ammonium nitrate']
    aer_names = {'Ammonium Sulphate': 'Ammonium sulphate', 'Aged fossil-fuel OC': 'Organic carbon', 'Ammonium nitrate': 'Ammonium nitrate'}

    # super rigid, needs replacing
    lam_colours = {'1064nm': 'r', '905nm': 'g'}


    # Q type to use in calculating f(RH)
    Q_type = 'extinction'
    print 'Q_type = ' + Q_type

    # ---------------------------------------------------
    # Read and Process
    # ---------------------------------------------------

    # create f(RH) for the main file
    f_RH_main, main_band_range, RH = create_f_RH(main_path, main_file_band, aer_index, aer_order, Q_type)

    # create f(RH) for all none-main files
    f_RH, _, _ = create_f_RH(file_paths, band, aer_index, aer_order, Q_type)

    # take differences [reference - main] -> +ve = reference is higher; -ve = reference is lower
    # keep separate to plotting, in order to improve code modularity
    f_RH_diff = calc_r_RH_difference(f_RH, f_RH_main, main_band_range)

    # calculate the ratio [reference / main] -> shows impact on extinction coefficient.
    # if ref = 2x high, then ext coeff will be 2x higher.
    f_RH_ratio = calc_r_RH_ratio(f_RH, f_RH_main, main_band_range)

    # ---------------------------------------------------
    # Plotting
    # ---------------------------------------------------

    # make band_order = [], its a rough way to plot it but it works... this way the plotting funcitons still work
    # with the multi-band versions

    fig = plot_absolute(f_RH_main, f_RH, RH, savedir, aer_order, aer_names, lam_colours, Q_type)

    # plot the difference in f(RH): lam_i - lam_reference
    fig = plot_difference(f_RH_diff, RH, savedir, aer_order, aer_names,lam_colours, Q_type)

    # plot the ratio in f(RH): lam_i / lam_reference
    fig = plot_ratio(f_RH_ratio, RH, savedir, aer_order, aer_names, lam_colours, Q_type)

    print 'END PROGRAM'

