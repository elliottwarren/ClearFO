"""
Create 4 panel plot with f(RH) differences for each of the aerosol species and the average f(RH)
Designed to deal with multiple bands

! WARNING - currently fixed to do 910, including save file names. Needs updating

Created by Elliott 06/04/17
"""

import sys
sys.path.append('C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/scripts')

import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from ellUtils import fig_majorAxis
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib import cm

from f_RH_difference import create_f_RH, calc_f_RH, calc_r_RH_difference, calc_r_RH_ratio
from f_RH_creation import read_spec_bands, read_aer_data


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
    band_order = []

    for band_i in bands:

        # read in the spectral band information
        spec_bands = read_spec_bands(file_path)

        # wavelength range in current band
        band_idx = np.where(spec_bands['band'] == band_i)[0][0]
        band_lam_range = '%.0f' % (spec_bands['lower_limit'][band_idx] * 1.0e9) + '-' + \
                         '%.0f' % (spec_bands['upper_limit'][band_idx] * 1.0e9) + 'nm'

        # This will be used for plotting in order, later on
        band_order += [band_lam_range]

        # read the aerosol data
        data = read_aer_data(file_path, aer_index, aer_order, band=band_i)

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

    return f_RH, band_order, RH

def create_colours(intervals):

    r = 0.2
    g = 0.5

    step = 1.0 / intervals

    b_colours = np.arange(0.0, 1.0 + step, step)
    g_colours = np.arange(1.0, 0.0 - step, -step)

    colours = [[r, g_colours[i], b_colours[i]] for i in range(len(b_colours))]

    return colours

def plot_absolute(f_RH_main, f_RH, RH, savedir, aer_order, aer_names, band_order, lam_colours, Q_type):

    # get wavelength of f_RH_main
    main_lam = f_RH_main.keys()[0]

    # plot the data on a 4 panel plot
    fig, ax = plt.subplots(2, 2, figsize=(8, 6))

    # append MURK to the aer_order to create a plotting order and aer_names so it's name can be plotted
    aer_plot_order = aer_order + ['MURK']
    aer_names['MURK'] = 'MURK'

    # plot on each subplot fully, before doing the next one
    for ax_i, aerosol in zip(ax.flatten(), aer_plot_order):

        for lam_range, lam_i_colour in zip(band_order, lam_colours):

            ax_i.plot(RH, f_RH[lam_range][aerosol], color=lam_i_colour, linestyle='--', label=lam_range)

        ax_i.plot(RH, f_RH_main[main_lam][aerosol], color='black', linestyle='--', label=main_lam)


        # prettify subplot
        ax_i.set_title(aer_names[aerosol])
        # ax_i.set_ylim([-0.2, 0.2])
        # ax_i.set_ylim([-1.0, 1.0])
        ax_i.set_xlim([0.0, 100.0])

    # prettify figure
    ax_main = fig_majorAxis(fig)
    ax_main.set_xlabel('RH [%]')
    plt.tight_layout()
    plt.subplots_adjust(top=0.9, right=0.8)
    ax.flatten()[1].legend(fontsize=8, bbox_to_anchor=(1.07, 1), loc=2, borderaxespad=0.0)
    fig.suptitle('Absolute value (lam_i)')
    plt.savefig(savedir + 'absolute/multi-panel-band_absolute_vs905'+ '_' + Q_type[0:3] + '_f_RH.png')
    plt.close()

    return fig

def plot_difference(f_RH_diff, RH, savedir, aer_order, aer_names, band_order, lam_colours, Q_type):

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
        for lam_range, lam_i_colour in zip(band_order, lam_colours):

            ax_i.plot(RH, f_RH_diff[lam_range][aerosol], color=lam_i_colour, label=lam_range)


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

    fig.suptitle('Absolute difference (lam_i - 905nm)')
    plt.savefig(savedir + 'difference/multi-panel-band_difference_vs905'+ '_' + Q_type[0:3] + '_f_RH.png')
    plt.close()

    return fig

def plot_ratio(f_RH_ratio, RH, savedir, aer_order, aer_names, band_order, Q_type):

    """
    Plot the ratio of the f(RH) references compared to the main f(RH) curve

    :param f_RH_diff:
    :param RH:
    :param savedir:
    :param aer_order:
    :param Q_type:
    :return: fig
    """

    # colour map to use in plotting
    cmap_range = cm.coolwarm(np.linspace(0, 1, len(band_order)))

    # plot the data on a 4 panel plot
    fig, ax = plt.subplots(2, 2, figsize=(8, 6))

    # append MURK to the aer_order to create a plotting order and aer_names so it's name can be plotted
    aer_plot_order = aer_order + ['MURK']
    aer_names['MURK'] = 'MURK'

    # plot on each subplot fully, before doing the next one
    for ax_i, aerosol in zip(ax.flatten(), aer_plot_order):

        # plot all curves for a single aerosol, at a time
        for lam_range, cmap_range_i in zip(band_order, cmap_range):

            ax_i.plot(RH, f_RH_ratio[lam_range][aerosol], color=cmap_range_i, label=lam_range)

        # prettify subplot
        ax_i.set_title(aer_names[aerosol])
        ax_i.set_ylim([0.975, 1.025])
        ax_i.set_xlim([0.0, 100.0])
        ax_i.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))

    # prettify figure
    ax_main = fig_majorAxis(fig)
    ax_main.set_xlabel('RH [%]')
    # ax_main.set_ylabel(r'$f_{ext},\/\lambda\/=\/i$' + '\n' + r'$f_{ext},\/\lambda\/=\/910nm$')
    ax_main.set_ylabel(r'$\frac{f_{ext,rh}(\lambda=i)}{f_{ext,rh}(\lambda=905\/nm)}$', labelpad=15, fontsize=18)
    plt.tight_layout()
    plt.subplots_adjust(top=0.9, right=0.8)
    ax.flatten()[1].legend(fontsize=8, bbox_to_anchor=(1.07, 1), loc=2, borderaxespad=0.0)
    # fig.suptitle('Ratio (lam_i / 910nm)')
    plt.savefig(savedir + 'ratio/multi-panel-band_ratio_vs905'+ '_' + Q_type[0:3] + '_f_RH.png')
    plt.close()

    return fig

if __name__ == '__main__':

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

    # main_file = 'sp_ew_910'
    main_file = 'sp_ew_ceil_multi_895-915'
    main_path = specdir + main_file
    main_file_band = 11 # 905 - 906

    band_file = 'sp_ew_ceil_multi_895-915'
    band_path = specdir + band_file

    # index: variables to take from file (as listed within the file) with index from BLOCK = 0
    #   NOTE: data MUST be in ascending index order
    # order: plotting order
    # names: titles to apply to each subplot
    aer_index = {'Accum. Sulphate': 1, 'Aged fossil-fuel OC': 2, 'Ammonium nitrate': 3}
    aer_order = ['Accum. Sulphate', 'Aged fossil-fuel OC', 'Ammonium nitrate']
    aer_names = {'Accum. Sulphate': 'Ammonium sulphate', 'Aged fossil-fuel OC': 'Organic carbon', 'Ammonium nitrate': 'Ammonium nitrate'}

    # super rigid, needs replacing
    # lam_colours = {'1062-1066nm': 'r', '908-912nm': 'b', '903-907nm': 'g'}


    # Q type to use in calculating f(RH)
    Q_type = 'extinction'
    print 'Q_type = ' + Q_type

    # ---------------------------------------------------
    # Read and Process
    # ---------------------------------------------------

    # create f(RH) for the main file
    f_RH_main, main_band_range, RH = create_f_RH(main_path, main_file_band, aer_index, aer_order, Q_type)

    # create f(RH) for all bands
    f_RH, band_order, RH = create_f_RH_multi_band(band_path, bands, aer_index, aer_order, Q_type)

    # take differences [reference - main] -> +ve = reference is higher; -ve = reference is lower
    # keep separate to plotting, in order to improve code modularity
    f_RH_diff = calc_r_RH_difference(f_RH, f_RH_main, main_band_range)

    # calculate the ratio [reference / main] -> shows impact on extinction coefficient.
    # if ref = 2x high, then ext coeff will be 2x higher.
    f_RH_ratio = calc_r_RH_ratio(f_RH, f_RH_main, main_band_range)

    # ---------------------------------------------------
    # Plotting
    # ---------------------------------------------------

    # get colour gradation based on number of lamda ranges in band_order
    lam_colours = create_colours(len(band_order) - 1)

    # plot the absolute value of f(RH): lam_i
    fig = plot_absolute(f_RH_main, f_RH, RH, savedir, aer_order, aer_names, band_order, lam_colours, Q_type)

    # plot the difference in f(RH): lam_i - lam_reference
    fig = plot_difference(f_RH_diff, RH, savedir, aer_order, aer_names, band_order, lam_colours, Q_type)

    # plot the ratio in f(RH): lam_i / lam_reference
    fig = plot_ratio(f_RH_ratio, RH, savedir, aer_order, aer_names, band_order, lam_colours, Q_type)

    print 'END PROGRAM'
