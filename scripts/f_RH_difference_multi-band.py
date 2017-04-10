"""
Create 4 panel plot with f(RH) differences for each of the aerosol species and the average f(RH)
Designed to deal with multiple bands

! WARNING - currently fixed to do 910, including save file names. Needs updating

Created by Elliott 06/04/17
"""

import numpy as np
import matplotlib.pyplot as plt
from ellUtils import fig_majorAxis
from f_RH_creation import read_aer_data, read_spec_bands, calc_f_RH


def create_f_RH_multi_band(file_path, bands, aer_index, aer_order, Q_type):

    """
    Read in data and create the f_RH data for each aerosol, as well as an average.

    :param f_RH_file:
    :param band:
    :param aer_index:
    :param aer_order:
    :return:

    Designed to only do this for one file at a time, in order to split the main file from the others.
    """

    f_RH = {}

    for band_i in bands:

        # read in the spectral band information
        spec_bands = read_spec_bands(file_path)

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

        # create an average f(RH)
        f_RH[band_lam_range]['average'] = np.mean(f_RH[band_lam_range].values(), axis=0)
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


def plot_difference(f_RH_diff, RH, savedir, aer_order, lam_colours, Q_type):

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
    fig, ax = plt.subplots(2, 2, figsize=(8, 8))

    # append the average to the aer_order to create a plotting order
    aer_plot_order = aer_order + ['average']

    # plot on each subplot fully, before doing the next one
    for ax_i, aerosol in zip(ax.flatten(), aer_plot_order):

        # plot all curves for a single aerosol, at a time
        for lam_range in f_RH_diff.iterkeys():

            ax_i.plot(RH, f_RH_diff[lam_range][aerosol], color=lam_colours[lam_range], label=lam_range)


        # prettify subplot
        ax_i.set_title(aerosol)
        ax_i.legend()
        # ax_i.set_ylim([-0.2, 0.2])
        ax_i.set_ylim([-1.0, 1.0])
        ax_i.set_xlim([0.0, 100.0])

    # prettify figure
    ax_main = fig_majorAxis(fig)
    ax_main.set_xlabel('RH [%]')
    # ax_main.set_ylabel('difference \n lam_i - 910nm')
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)
    fig.suptitle('Absolute difference (lam_i - 910nm)')
    plt.savefig(savedir + 'difference/multi-panel_difference_vs910'+ '_' + Q_type[0:3] + '_f_RH_1.png')
    plt.close()

    return fig


def plot_ratio(f_RH_ratio, RH, savedir, aer_order, lam_colours, Q_type):

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
    fig, ax = plt.subplots(2, 2, figsize=(8, 8))

    # append the average to the aer_order to create a plotting order
    aer_plot_order = aer_order + ['average']

    # plot on each subplot fully, before doing the next one
    for ax_i, aerosol in zip(ax.flatten(), aer_plot_order):

        # plot all curves for a single aerosol, at a time
        for lam_range in f_RH_ratio.iterkeys():

            ax_i.plot(RH, f_RH_ratio[lam_range][aerosol], color=lam_colours[lam_range], label=lam_range)

        # prettify subplot
        ax_i.set_title(aerosol)
        ax_i.legend()
        ax_i.set_ylim([0.9, 1.3])
        ax_i.set_xlim([0.0, 100.0])

    # prettify figure
    ax_main = fig_majorAxis(fig)
    ax_main.set_xlabel('RH [%]')
    # ax_main.set_ylabel('difference \n lam_i - 910nm')
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)
    fig.suptitle('Ratio (lam_i / 910nm)')
    plt.savefig(savedir + 'ratio/multi-panel_ratio_vs910'+ '_' + Q_type[0:3] + '_f_RH.png')
    plt.close()

    return fig


def main():

    # User set args
    # band that read_spec_bands() uses to find the correct band
    #! Manually set
    bands = np.arange(1, 21)

    # -------------------------

    # directories
    savedir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/figures/Mie/f(RH)/'
    specdir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/data/Mie/spectral/'
    f_RHdir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/data/Mie/'

    # file_name = 'spec3a_sw_hadgem1_7lean_so' # original file given to me by Claire Ryder 25/01/17
    # file_name = 'sp_sw_ga7' # current UM file
    # file_name = 'sp_ew_910' # my own made file with 1 band at 910 nm

    main_file = 'sp_ew_910'
    main_path = specdir + main_file

    file_names = ['sp_ew_ceil_guass_903-907',
                  'sp_ew_ceil_guass_908-912',
                  'sp_ew_ceil_guass_1062-1066']
    file_paths = [specdir + f for f in file_names]

    # variables to take from file (as listed within the file) with index from BLOCK = 0
    # NOTE: data MUST be in ascending index order
    aer_index = {'Accum. Sulphate': 1, 'Aged fossil-fuel OC': 2, 'Ammonium nitrate': 3}
    aer_order = ['Accum. Sulphate', 'Aged fossil-fuel OC', 'Ammonium nitrate']

    # super rigid, needs replacing
    lam_colours = {'1062-1066nm': 'r', '908-912nm': 'b', '903-907nm': 'g'}


    # Q type to use in calculating f(RH)
    Q_type = 'extinction'
    print 'Q_type = ' + Q_type

    # ---------------------------------------------------
    # Read and Process
    # ---------------------------------------------------

    # create f(RH) for the main file
    f_RH_main, main_band_range, RH = create_f_RH(main_path, band, aer_index, aer_order, Q_type)

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

    # plot the difference in f(RH): lam_i - lam_reference
    fig = plot_difference(f_RH_diff, RH, savedir, aer_order, lam_colours, Q_type)

    # plot the ratio in f(RH): lam_i / lam_reference
    fig = plot_ratio(f_RH_ratio, RH, savedir, aer_order, lam_colours, Q_type)

    print 'END PROGRAM'

if __name__ == '__main__':
    main()
