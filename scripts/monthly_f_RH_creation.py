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
import pickle
from netCDF4 import Dataset
import datetime as dt

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

    print "Reading from:" + file_path
    print "Band: " + str(band)

    # define arrays
    data = {}

    for aer in aer_order:


        print 'Getting ...' + aer

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
    :param f_RH: f(RH) curves
    :param radii_range_nm: needs to match the units if changed! [currently nm]
    :param RH_int: interpolated RH values [fraction]
    :param site_ins: site information
    :param ceil_lambda_nm_str: ceilometer wavelength as str with nm on the end e.g. '905.0nm'
    :return:

    Store the MURK [month, radius, RH] AND species f(RH) curves [radius, RH]
    """

    # create netCDF file
    ncfile = Dataset(fRHdir + 'monthly_f(RH)_' + site_ins['site_short'] + '_' + ceil_lambda_nm_str + '.nc', 'w')

    # create Dimensions
    ncfile.createDimension('month', len(np.arange(1, 13)))
    ncfile.createDimension('RH', len(RH_int))
    ncfile.createDimension('radii_range', len(radii_range_nm))

    # create, fill and set units for co-ordinate variables
    nc_month = ncfile.createVariable('months', np.float64, ('month',))
    nc_RH = ncfile.createVariable('Relative Humidity', np.float64, ('RH',))
    nc_radii_range_nm = ncfile.createVariable('radii_range_nm', np.float64, ('radii_range',))

    nc_month[:] = np.arange(1, 13)
    nc_RH[:] = RH_int
    nc_radii_range_nm[:] = radii_range_nm

    nc_month.units = 'month number'
    nc_RH.units = 'fraction'
    nc_radii_range_nm.units = 'nm'

    # create and fill variables
    nc_f_RH_MURK = ncfile.createVariable('f(RH) MURK', np.float64, ('month', 'radii_range', 'RH'))
    nc_f_RH_MURK[:] = f_RH['MURK']
    for species_i in f_RH.keys():
        var_name = 'f(RH) ' + species_i
        nc_f_RH_species_i = ncfile.createVariable(var_name, np.float64, ('radii_range', 'RH'))
        nc_f_RH_species_i[:] = f_RH[species_i]

    # extra attributes
    ncfile.history = 'Created ' + dt.datetime.now().strftime('%Y-%m-%d %H:%M') + ' GMT'

    return

if __name__ == '__main__':

    # ------------------------------------------
    # Setup
    # ------------------------------------------
    # site information
    site_ins = {'site_short': 'NK', 'site_long': 'North Kensington',
                'ceil_lambda': 0.905e-06, 'land-type': 'urban'}
    # site_ins = {'site_short':'Ch', 'site_long': 'Chilbolton',
    #             'ceil_lambda': 0.905e-06, 'land-type': 'rural'}
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
    savedir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/figures/Mie/f(RH)/'
    fRHdir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/data/Mie/monthly_f(RH)/'
    specdir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/data/Mie/' \
              'monthly_f(RH)/sp_885-925_r_files/'
    pickleloaddir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/data/Mie/pickle/'


    # variables to take from file (as listed within the file) with index from BLOCK = 0
    # NOTE: data MUST be in ascending index order
    aer_index = {'Ammonium Sulphate': 1, 'Generic NaCl': 2, 'Biogenic': 3, 'Aged fossil-fuel OC': 4,
                 'Ammonium nitrate': 5}
    aer_order = ['Ammonium Sulphate', 'Generic NaCl', 'Biogenic', 'Aged fossil-fuel OC', 'Ammonium nitrate']

    aer_particles_chem = {'Ammonium Sulphate': '(NH4)2SO4', 'Generic NaCl': 'NaCl', 'Aged fossil-fuel OC': 'CORG',
                          'Ammonium nitrate': 'NH4NO3', 'Black carbon': 'CBLK'}

    aer_particles = ['(NH4)2SO4', 'NH4NO3', 'NaCl', 'CORG', 'CBLK']

    # Q type to use in calculating f(RH)
    Q_type = 'extinction'
    print 'Q_type = ' + Q_type

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
    radii_range_m = np.arange(0.005e-06, 2.885e-6, 0.005e-06)
    radii_range_micron = np.arange(0.005, 2.885, 0.005)
    radii_range_nm = np.arange(5, 2885 + 5, 5)

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

        print radius_idx

        # format of radius used in the filename
        #   trying to use m or microns leads to rounding errors when making the string...
        radius_filestr = '0.%09d' % radius_nm_i

        # create filename
        filename = 'sp_885-925_r'+radius_filestr+'_stdev1.6_num4.461e9'
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
        # all species but not including soot
        for species_i in aer_index.iterkeys():
            f = interp1d(RH, f_RH_i[species_i], kind='linear')
            interp_f_RH_i[species_i] = f(RH_int)
            f_RH[species_i][radius_idx, :] = interp_f_RH_i[species_i]


            # store the interpolated f(RH) for each species


        # make f(RH) for murk from the interpolated f(RH)
        for month_idx in range(12):
            f_RH['MURK'][month_idx, radius_idx, :] = \
                interp_f_RH_i['Ammonium Sulphate'] * pm10_rel_vol['(NH4)2SO4'][month_idx] + \
                interp_f_RH_i['Ammonium nitrate'] * pm10_rel_vol['NH4NO3'][month_idx] + \
                interp_f_RH_i['Aged fossil-fuel OC'] * pm10_rel_vol['CORG'][month_idx] + \
                interp_f_RH_i['Soot'] * pm10_rel_vol['CBLK'][month_idx] + \
                interp_f_RH_i['Generic NaCl'] * pm10_rel_vol['NaCl'][month_idx]


    # save f(RH) once all radii have been looped through
    if saveFRH == True:

        save_fRH_netCDF(fRHdir, f_RH, radii_range_nm, RH_int, site_ins, ceil_lambda_nm_str)












    # min and max possible ranges in f(RH)
    # https://books.google.co.uk/books?id=OnmyzXfS6ggC&pg=PA59&lpg=PA59&dq=mass+closure+aerosol+london&source=bl&ots=Nd9L-GBPp_&sig=0KzD_
    # USvtga8-c4fv674hPKnA1E&hl=en&sa=X&ved=0ahUKEwiDwqPeqLHUAhXIKcAKHSmsBW8Q6AEIVzAH#v=onepage&q=mass%20closure%20aerosol%20london&f=false
    # U is capitalised...

    # ---------------------------------------------------
    # Plotting
    # ---------------------------------------------------

    # plot the data at the end so it is all neat and together
    fig = plt.figure(figsize=(6, 4))

    for key, value in data.iteritems():

        plt.plot(RH*100, f_RH[key], label=key, linestyle='-')

    # plt.plot(value[:, 0], f_RH['average with Aitken Sulphate'], label='average with Aitken Sulphate', linestyle='-')

    # plot the MURK one
    # plt.plot(value[:, 0], f_RH['average'], label='average without Aitken Sulphate', color='black')
    plt.plot(RH*100, f_RH['MURK'], label='MURK', color='black')

    # plot soot as a constant until I get f(RH) for it
    plt.plot([0.0, 100.0], [1.0, 1.0], label='Soot')

    plt.legend(fontsize=10, loc='best')
    plt.tick_params(axis='both', labelsize=11)
    plt.xlabel('RH [%]', labelpad=0, fontsize=11)
    # plt.ylabel(Q_type + ' f(RH)')
    plt.ylabel(r'$f_{ext,rh}$', fontsize=11)
    plt.ylim([0.0, 8.0])
    plt.xlim([0.0, 100.0])
    # plt.title(file_name + ': ' + band_lam_range + ' band')

    plt.tight_layout() # moved tight_layout above... was originally after the save (06/04/17)
    plt.savefig(savedir + file_name + '_' + Q_type[0:3] + '_f_RH_all_' + band_lam_range + '_salt8.0.png')

    print 'END PROGRRAM'