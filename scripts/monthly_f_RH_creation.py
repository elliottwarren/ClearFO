"""
Create a monthly climatology of f(RH) for use in calculating the extinction efficiency (Q)

Created by Elliott 01/03/2018
Based on f_RH_creation.py
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy import interpolate
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
import pickle
from netCDF4 import Dataset
import datetime as dt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import ellUtils as eu

# Reading

def read_spec_bands(file_path):

    """
    Reads in the spectral bands from UM spectral file

    :param file_path:
    :return: spec_bands:
    """

    spec_bands = {'band': [],
                  'lower_limit': [],
                  'upper_limit': []}

    # Read data

    file = open(file_path, "r")
    line = file.readline()
    line = line.rstrip('\n\r')

    # skip until spectral line is reached
    while line != 'Limits of spectral intervals (wavelengths in m.)':
        line = file.readline()
        line = line.rstrip('\n\r')

    line = file.readline() # pass over the header

    line = file.readline() # read the first actual line
    line = line.rstrip('\n\r')


    while line != '*END':

        line = ' '.join(line.split()) # remove duplicate, trailing and leading spaces
        line_split = line.split(' ')

        spec_bands['band'] += [int(line_split[0])]
        spec_bands['lower_limit'] += [float(line_split[1])]
        spec_bands['upper_limit'] += [float(line_split[2])]

        line = file.readline()
        line = line.rstrip('\n\r')

    # convert to numpy array
    for key, value in spec_bands.iteritems():
        spec_bands[key] = np.array(value)

    file.close()

    return spec_bands

def read_aer_data(file_path, aer_index, aer_order, band=1):

    """
    Read in the aerosol data, calculate the Qs

    :param file_path:
    :param aer_index:
    :param aer_order:
    :param band
    :return:
    """

    # print "Reading from:" + file_path
    # print "Band: " + str(band)

    # define arrays
    data = {}

    for aer in aer_order:


        # print 'Getting ...' + aer

        data[aer] = []

        file = open(file_path, "r")

        line = file.readline()  # read the first actual line
        line = line.rstrip('\n\r')

        while line != '*BLOCK: TYPE =   11: SUBTYPE =    1: VERSION =    2':
            line = file.readline()
            line = line.rstrip('\n\r')

        # while line != 'Index of species = 6 Accum. Sulphate':
        while line != 'Index of species = ' + str(aer_index[aer]) + ' ' + aer:
            line = file.readline()
            line = line.rstrip('\n\r')
            line = ' '.join(line.split())

        while line != 'Band = ' + str(band):
            line = file.readline()
            line = line.rstrip('\n\r')
            line = ' '.join(line.split())


        # skip two lines
        line = file.readline()
        line = file.readline()

        # read first line of data
        line = file.readline()
        line = line.rstrip('\n\r')

        # when reached the data line...
        while (line != 'Band = ' + str(band + 1)) & (line != '*END'):

            line = ' '.join(line.split()) # remove duplicate, trailing and leading spaces
            line_split = line.split(' ')

            # extract data
            data[aer] += [line_split]

            line = file.readline()
            line = line.rstrip('\n\r')
            line = ' '.join(line.split())

        file.close()

        # convert to numpy array, list already alligned correctly for np.array()
        data[aer] = np.array(data[aer], dtype=float)


    return data

# Processing

def calc_f_RH(data, aer_order, Q_type=''):

    """
    Calculate f(RH) for the Q type
    :param data:
    :param aer_order:
    :param Q_type:
    :return:
    """

    # define arrays
    Q = {}
    f_RH = {}

    for aer in aer_order:

        if Q_type == 'extinction':
            # calculate total extinction (abs + scatt)
            Q[aer] = np.sum((data[aer][:, 1], data[aer][:, 2]), axis=0)

        elif Q_type == 'scattering':
            # extract scattering
            Q[aer] = data[aer][:, 2]

        elif Q_type == 'absorption':
            # extract absorption
            Q[aer] = data[aer][:, 1]
        else:
            raise ValueError("Incorrect Q_type value. Needs to be either 'extinction', 'absorption', or 'scattering'")

        # calculate f(RH) = Q(RH>=0)/Q(RH=0)
        f_RH[aer] = np.array([Q_RH / Q[aer][0] for Q_RH in Q[aer]])

    return Q, f_RH

# Saving

def save_fRH_netCDF(fRHdir, f_RH, radii_range_nm, RH_int, site_ins, ceil_lambda_nm_str):

    """
    Create a netCDF file for the f(RH) information to be stored in.
    :param fRHdir: dircetory to save netCDF file to (should ideally be the monthly f(RH) dir)
    :param f_RH: f(RH) curves: Dictionary variable, where keys are species and each element
                is a [3D array] with .shape(month, radii, RH)
    :param radii_range_nm: needs to match the units if changed! [currently nm]
    :param RH_int: interpolated RH values [fraction]
    :param site_ins: site information
    :param ceil_lambda_nm_str: ceilometer wavelength as str with nm on the end e.g. '905.0nm'
    :return:

    Store the MURK [month, radius, RH] AND species f(RH) curves [radius, RH]
    EW 03/2018
    """

    # create netCDF file
    filepath = fRHdir + 'monthly_f(RH)_' + site_ins['site_short'] + '_' + ceil_lambda_nm_str + '.nc'
    ncfile = Dataset(filepath, 'w')

    # create Dimensions
    ncfile.createDimension('month', len(np.arange(1, 13)))
    ncfile.createDimension('RH', len(RH_int))
    ncfile.createDimension('radii_range', len(radii_range_nm))

    # create, fill and set units for co-ordinate variables
    nc_month = ncfile.createVariable('month', np.float64, ('month',))
    nc_RH = ncfile.createVariable('RH', np.float64, ('RH',))
    nc_radii_range_nm = ncfile.createVariable('radii_range', np.float64, ('radii_range',))

    nc_month[:] = np.arange(1, 13)
    nc_RH[:] = RH_int
    nc_radii_range_nm[:] = radii_range_nm

    nc_month.units = 'month number'
    nc_RH.units = 'fraction'
    nc_radii_range_nm.units = 'nm'

    # create and fill variables
    for species_i in f_RH.keys():
        var_name = 'f(RH) ' + species_i

        if var_name == 'f(RH) MURK':
            nc_f_RH_MURK = ncfile.createVariable('f(RH) MURK', np.float64, ('month', 'radii_range', 'RH'))
            nc_f_RH_MURK[:] = f_RH['MURK']
        else:
            nc_f_RH_species_i = ncfile.createVariable(var_name, np.float64, ('radii_range', 'RH'))
            nc_f_RH_species_i[:] = f_RH[species_i]

    # extra attributes
    ncfile.history = 'Created ' + dt.datetime.now().strftime('%Y-%m-%d %H:%M') + ' GMT'

    # close file
    ncfile.close()

    print filepath + ' saved!'

    return

# plotting

def create_colours(intervals):

    """
    create a simple range of colours (RGB)
    :param intervals: number of intervals colour range will have
    :return: colours: the colour range
    """

    r = 0.2
    g = 0.5

    step = 1.0 / intervals

    b_colours = np.arange(0.0, 1.0, step)
    g_colours = np.arange(1.0, 0.0, -step)

    colours = [[r, g_colours[i], b_colours[i]] for i in range(len(b_colours))]

    return colours

if __name__ == '__main__':

    # ------------------------------------------
    # Setup
    # ------------------------------------------
    # site information
    site_ins = {'site_short': 'NK', 'site_long': 'North Kensington',
                'ceil_lambda': 0.905e-06, 'land-type': 'urban', 'SOCRATES_file_subdir': 'stdev_1.70_NK',
                'geo_stdev': '1.70'}
    # site_ins = {'site_short':'Ch', 'site_long': 'Chilbolton',
    #             'ceil_lambda': 0.905e-06, 'land-type': 'rural'}
    # site_ins = {'site_short':'Ha', 'site_long': 'Harwell',
    #             'ceil_lambda': 0.905e-06, 'land-type': 'rural'}

    ceil_lambda_nm_str = str(int(site_ins['ceil_lambda'] * 1e9)) + 'nm'

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
              'SOCRATES/'+site_ins['SOCRATES_file_subdir']+'/'
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

    # range of colours (0.1 % resolution)
    # res = 1000
    cmap = cm.coolwarm
    # colours = cm.jet(np.linspace(0, 1, res))
    # norm = mpl.colors.Normalize(vmin=0, vmax=100)

    # ---------------------------------------------------
    # Read, Process and save f(RH)
    # ---------------------------------------------------

    # read in relative volume of aerosol species
    #   created in VMachine "create_Qext_for_MURK_.py" and saved as a pickle
    # # read in any pickled S data from before
    filename = pickleloaddir + site_ins['site_short'] + '_aerosol_relative_volume.pickle'
    with open(filename, 'rb') as handle:
        pickle_load_in = pickle.load(handle)
    pm10_rel_vol = pickle_load_in['pm10_rel_vol']

    # range of radii to iterate over
    # radii_range_m = np.arange(0.005e-06, 3.685e-6 + 0.005e-06, 0.005e-06)
    radii_range_m = np.arange(0.005e-06, 7.720e-06 + 0.005e-06, 0.005e-06)
    radii_range_micron = np.arange(0.005, 7.720 + 0.005, 0.005)
    radii_range_nm = np.arange(5, 7720 + 5, 5)

    # RH to interpolate to
    RH_int = np.arange(0, 1.01, 0.01)

    # create f(RH) dictionaries to fill (including one for murk)
    #   MURK data will be in a 3D array [month, size, RH]
    #   all other speices will be 2D [size, RH] as they wont vary each month
    f_RH = {}
    f_RH['MURK'] = np.empty((12, len(radii_range_nm), 101))
    f_RH['MURK'][:] = np.nan
    for species_i in aer_particles:
        f_RH[species_i] = np.empty((len(radii_range_nm), 101))
        f_RH[species_i][:] = np.nan


    for radius_idx, radius_nm_i in enumerate(radii_range_nm):

        print radius_nm_i

        # format of radius used in the filename
        #   trying to use m or microns leads to rounding errors when making the string...
        radius_filestr = '0.%09d' % radius_nm_i

        # create filename
        filename = 'sp_885-925_r'+radius_filestr+'_stdev'+site_ins['geo_stdev']+'_num4.461e9'
        file_path = specdir + filename

        # read in the spectral band information
        spec_bands = read_spec_bands(file_path)

        # wavelength range in current band
        band_idx = np.where(spec_bands['band'] == band)[0][0]
        band_lam_range = '%.0f' % (spec_bands['lower_limit'][band_idx] * 1.0e9) + '-' + \
                         '%.0f' % (spec_bands['upper_limit'][band_idx] * 1.0e9) + 'nm'

        # read the aerosol data
        data = read_aer_data(file_path, aer_index, aer_order, band=band)

        # Extract RH (RH is the same across all aerosol types)
        RH = np.array(data[aer_order[0]][:, 0])

        # calculate f(RH) for each species
        Q, f_RH_i = calc_f_RH(data, aer_order, Q_type=Q_type)

        # add a soot f_RH
        f_RH_i['Soot'] = np.array([1.0 for i in f_RH_i['Generic NaCl']])

        # linearly interpolate f(RH) to increase resolution from 0.05 to 0.01 [fraction]
        interp_f_RH_i = {}

        # soot is always the same (fixed at 1)
        interp_f_RH_i['Soot'] = np.repeat(1.0, 101)
        f_RH['CBLK'][radius_idx, :] = interp_f_RH_i['Soot']

        # all species excluding soot
        for species_i, chem_i in aer_particles_chem.iteritems():
            f = interp1d(RH, f_RH_i[species_i], kind='linear')
            interp_f_RH_i[species_i] = f(RH_int)
            f_RH[chem_i][radius_idx, :] = interp_f_RH_i[species_i]

        # make f(RH) for murk from the interpolated f(RH)
        for month_idx in range(12):
            f_RH['MURK'][month_idx, radius_idx, :] = \
                (interp_f_RH_i['Ammonium Sulphate'] * pm10_rel_vol['(NH4)2SO4'][month_idx]) + \
                (interp_f_RH_i['Ammonium nitrate'] * pm10_rel_vol['NH4NO3'][month_idx]) + \
                (interp_f_RH_i['Aged fossil-fuel OC'] * pm10_rel_vol['CORG'][month_idx]) + \
                (interp_f_RH_i['Soot'] * pm10_rel_vol['CBLK'][month_idx]) + \
                (interp_f_RH_i['Generic NaCl'] * pm10_rel_vol['NaCl'][month_idx])


    # save f(RH) once all radii have been looped through
    if saveFRH == True:

        save_fRH_netCDF(fRHdir, f_RH, radii_range_nm, RH_int, site_ins, ceil_lambda_nm_str)


    # ---------------------------------------------------
    # Plotting
    # ---------------------------------------------------

    # 1. quick pcolor plot - seems like f(RH) has a low sensitivity to radius, above ~0.3 microns
    for month_idx in range(12):
        fig = plt.figure(figsize=(6, 4))
        plt.pcolor(radii_range_micron, RH_int * 100.0, np.transpose(f_RH['MURK'][month_idx, :, :]), vmin=1, vmax=13)
        plt.colorbar()
        ax = plt.gca()
        ax.set_xscale('log')
        plt.xlabel('radius [microns]')
        plt.ylabel('RH [%]')
        monthstr = dt.datetime(1900, month_idx+1, 1).strftime('%B')
        plt.suptitle(site_ins['site_short'] + ': f(RH) MURK' + ' - ' + monthstr)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(savedir + site_ins['site_short'] + '_f(RH)_MURK_'+str(month_idx+1))
        plt.close(fig)


    # 3. plot f(RH) for MURK, with respect to RH, for large particle sizes (> 0.4 microns)
    # wrt_radii_radii_ratio
    for month_idx in range(12):

        date_i = dt.datetime(1900, month_idx+1, 1).strftime('%b')

        fig = plt.figure(figsize=(6, 4))

        # find 0.11 ahead of time
        rad_0p11_idx = np.where(radii_range_micron == 0.11)[0][0]

        # loop through an arbitrary list of RH values to plot
        for rad_val in [0.07, 0.11, 0.14, 0.4, 1.0, 3.0]:

            # find where rh_val is
            rad_idx = np.where(radii_range_micron == rad_val)[0][0]
            # print rad_idx

            # ratio of f(RH) for this month / average across all months
            f_RH_plot_data = np.squeeze(f_RH['MURK'][month_idx, rad_idx, :]) \
                             / f_RH['MURK'][month_idx, rad_0p11_idx, :]

            plt.plot(RH_int*100.0, f_RH_plot_data, label=str(radii_range_micron[rad_idx]),
                     linestyle='-')
            # plt.axvline(0.11, color='grey', alpha=0.3, linestyle='--')

        plt.xlabel('RH [%]')
        plt.ylabel('f(RH)')
        plt.ylim([0.0, 2.0])
        plt.legend(loc='top left')
        plt.suptitle(date_i)
        plt.savefig(savedir + 'wrt_radii_radii_ratio/' + 'wrt_radii_ratio' + site_ins['site_short'] + '_f(RH)_MURK_'+str(month_idx+1))
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

    # 5. RATIO plot f(RH) for MURK, with respect to RH, for large particle sizes (> 0.4 microns)
    # wrt_radii_monthly_ratio
    # for rad_val in [0.07, 0.11, 0.14, 0.4, 1.0, 3.0]:
    for rad_val in [0.11]:

        for aer_i, aer_data in pm10_rel_vol.iteritems():

            # aer_i = 'CBLK'
            # aer_data = pm10_rel_vol['CBLK']

            # set up how to colour the lines for this aerosol
            #   what will the range be, given how variable the aerosol is across the months
            max_aer_val = aer_data.max()
            min_aer_val = aer_data.min()
            res = 1000
            # colours = cm.jet(np.linspace(0, 1, res))
            # c_range = np.linspace(0, max_aer_val, res)
            colours = cm.coolwarm(np.linspace(0, 1, res))
            c_range = np.linspace(min_aer_val, max_aer_val, res)
            norm = mpl.colors.Normalize(vmin=min_aer_val, vmax=max_aer_val)

            fig = plt.figure(figsize=(6, 4))
            ax = plt.gca()

            # loop through an arbitrary list of RH values to plot

            # find where rh_val is
            rad_idx = np.where(radii_range_micron == rad_val)[0][0]
            print rad_idx

            for month_idx, month_i in enumerate(range(1, 13)):

                # colour the month - find where along the scale, the current % species should lie
                #   res is the resolution of the colour range, so (*res) finds what colour along that it should have
                #   NOTE: 0.5 is rounded down using np.round()
                # colour_idx = int(np.round(pm10_rel_vol[aer_i][month_idx] * res))
                _, colour_idx, _ = eu.nearest(c_range, pm10_rel_vol[aer_i][month_idx])
                print colour_idx

                # # get colour and linestyle for month
                # for key, data in month_colours.iteritems():
                #     if month_i in data:
                #         colour_i = key
                #
                # for key, data in month_ls.iteritems():
                #     if month_i in data:
                #         ls = key

                # month as 3 letter str
                date_i = dt.datetime(1900, month_i, 1).strftime('%b')


                # ratio of f(RH) for this month / average across all months
                f_RH_plot_data = np.squeeze(f_RH['MURK'][month_idx, rad_idx, :]) \
                                 / np.mean(f_RH['MURK'][:, rad_idx, :], axis=0)

                handle = plt.plot(RH_int*100.0, f_RH_plot_data, label=date_i,
                         linestyle='-', color=colours[colour_idx, :])
                # plt.axvline(0.11, color='grey', alpha=0.3, linestyle='--')

            plt.xlabel('RH [%]')
            plt.ylabel('f(RH, month_i)/f(RH, average(months)')
            plt.ylim([0.8, 1.5])
            plt.legend(loc='top left', ncol=2)

            # plot colourbar
            divider = make_axes_locatable(ax)
            # cax = divider.append_axes("right", size="5%", pad=0.05)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cbar = mpl.colorbar.ColorbarBase(cax, cmap=cm.coolwarm, norm=norm)
            cbar.set_label(r'$species\/[\%]$', labelpad=-38, y=1.075, rotation=0)

            plt.suptitle('radii: '+ str(rad_val) + ' [microns]')
            plt.savefig(savedir + 'species_colour/' + 'ratio_' + site_ins['site_short'] + '_f(RH)_MURK_colour_by_'+aer_i+'.png')
            plt.close(fig)


    print 'END PROGRRAM'