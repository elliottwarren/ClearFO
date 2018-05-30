"""
Create the f(RH) graph and look up table for use in calculating the extinction efficiency (Q)

Created by Elliott 06/02/2017
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy import interpolate
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib

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

def read_aer_data(file_path, aer_index, aer_order, band=4):

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
        f_RH[aer] = [Q_RH / Q[aer][0] for Q_RH in Q[aer]]

    return Q, f_RH

if __name__ == '__main__':

    # User set args
    # band that read_spec_bands() uses to find the correct band
    #! Manually set
    band = 1

    # saveF(RH)?
    saveFRH = False

    # -------------------------

    # directories
    savedir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/figures/Mie/f(RH)/'
    specdir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/data/Mie/spectral/'
    f_RHdir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/data/Mie/'

    # file_name = 'spec3a_sw_hadgem1_7lean_so' # original file given to me by Claire Ryder 25/01/17
    # file_name = 'sp_sw_ga7' # current UM file
    # file_name = 'sp_ew_ceil_905'
    file_name = 'sp_895-915_r1.1e-7_stdev1.6_num8.0e9_salt8.0e-6'
    # file_name = 'sp_908-912_r1.1e-7_stdev1.6_num8.0e9' # 6 different aerosols inc. salt, biogenic and soot
    file_path = specdir + file_name

    # variables to take from file (as listed within the file) with index from BLOCK = 0
    # NOTE: data MUST be in ascending index order
    if file_name == 'spec3a_sw_hadgem1_7lean_so':
        aer_index = {'Accum. Sulphate': 6, 'Aitken Sulphate': 7, 'Aged fossil-fuel OC': 22}
        aer_order = ['Accum. Sulphate', 'Aitken Sulphate', 'Aged fossil-fuel OC']
    elif file_name == 'sp_ew_ceil_guass_903-907':
        aer_index = {'Accum. Sulphate': 1, 'Aged fossil-fuel OC': 2, 'Ammonium nitrate': 3}
        aer_order = ['Accum. Sulphate', 'Aged fossil-fuel OC', 'Ammonium nitrate']
    elif file_name == 'sp_895-915_r1.1e-7_stdev1.6_num8.0e9_salt8.0e-6':
        aer_index = {'Ammonium Sulphate': 1, 'Generic NaCl': 2, 'Biogenic': 3, 'Aged fossil-fuel OC': 4, 'Ammonium nitrate': 5}
        aer_order = ['Ammonium Sulphate', 'Generic NaCl', 'Biogenic', 'Aged fossil-fuel OC', 'Ammonium nitrate']
    elif file_name == 'sp_ew_ceil_905':
        aer_index = {'Ammonium Sulphate': 1, 'Generic NaCl': 2, 'Biogenic': 3, 'Aged fossil-fuel OC': 4,
                     'Ammonium nitrate': 5}
        aer_order = ['Ammonium Sulphate', 'Generic NaCl', 'Biogenic', 'Aged fossil-fuel OC', 'Ammonium nitrate']

    else:
        aer_index = {'Ammonium Sulphate': 2, 'Generic NaCl': 3, 'Biogenic': 4, 'Aged fossil-fuel OC': 5, 'Ammonium nitrate': 6}
        aer_order = ['Ammonium Sulphate', 'Generic NaCl', 'Biogenic', 'Aged fossil-fuel OC', 'Ammonium nitrate']
        # aer_index = {'Accum. Sulphate': 14, 'Aitken Sulphate': 15, 'Aged fossil-fuel OC': 24, 'Ammonium nitrate': 26}
        # aer_order = ['Accum. Sulphate', 'Aitken Sulphate', 'Aged fossil-fuel OC', 'Ammonium nitrate']

    # Q type to use in calculating f(RH)
    Q_type = 'extinction'
    print 'Q_type = ' + Q_type

    # ---------------------------------------------------
    # Read, Process and save f(RH)
    # ---------------------------------------------------

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

    # calculate f(RH)
    Q, f_RH = calc_f_RH(data, aer_order, Q_type=Q_type)



    # # interpolate f(RH) to increase resolution from 5% to 0.1%
    # # fit binomial
    # # m = np.polyfit(RH*100.0, f_RH_5['Generic NaCl'], 3)
    # tck = interpolate.splrep(RH*100.0, f_RH_5['Generic NaCl'], s=1)
    # ynew = interpolate.splev(np.arange(101), tck, der=0)
    #
    # (A, B), covvariances = curve_fit(lambda t, a, b: a * np.exp(b * t), RH[6:]*100.0, f_RH_5['Generic NaCl'][6:], p0=(2, 0.2))
    # f = interp1d(RH*100.0, f_RH_5['Generic NaCl'], kind='cubic') # output is a function
    # f2 = interp1d(RH * 100.0, f_RH_5['Generic NaCl'], kind='quadratic')

    # plot binomial
    #plt.plot(np.arange(101), ynew, label='cubic spline')
    #plt.plot(np.arange(101), (m[0] * (np.arange(101)**3)) + (m[1] * (np.arange(101)**2)) + (m[2]*(np.arange(101)**1)) + (m[3]), label='2deg')
    # plt.plot(np.arange(101), (m[0] * (np.arange(101)**5)) + (m[1]*(np.arange(101)**4)) + (m[2]*(np.arange(101)**3)) +
    #                          (m[3] * (np.arange(101) ** 2)) + (m[4] * (np.arange(101))) + m[5], label='5deg')
    #plt.plot(np.arange(101), A*np.exp(np.arange(101)*B), label='exp')
    #plt.plot(np.arange(101), f(np.arange(101)), label='interp1d cubic')
    #plt.plot(np.arange(101), f2(np.arange(101)), label='quadratic')
    # plt.plot(RH*100, f_RH['Generic NaCl'], label='orig')
    # plt.legend()

    # add a soot f_RH
    f_RH['Soot'] = [1.0 for i in f_RH['Generic NaCl']]

    # create an average f(RH)
    # f_RH['average with Aitken Sulphate'] = np.mean(f_RH.values(), axis=0)
    # f_RH['average'] = np.mean((f_RH['Ammonium Sulphate'], f_RH['Aged fossil-fuel OC'], f_RH['Ammonium nitrate']), axis=0)
    f_RH['MURK'] = (np.array(f_RH['Ammonium Sulphate']) * 0.295) + (0.38 * np.array(f_RH['Aged fossil-fuel OC'])) + \
    (0.325 * np.array(f_RH['Ammonium nitrate']))

    # save f(RH)
    if saveFRH == True:
        np.savetxt(f_RHdir +  file_name + '_ext_f(RH)_' + band_lam_range + '.csv',
                   np.transpose(np.vstack((RH, f_RH['MURK']))), delimiter=',', header='RH,f_RH')

    # min and max possible ranges in f(RH)
    # https://books.google.co.uk/books?id=OnmyzXfS6ggC&pg=PA59&lpg=PA59&dq=mass+closure+aerosol+london&source=bl&ots=Nd9L-GBPp_&sig=0KzD_
    # USvtga8-c4fv674hPKnA1E&hl=en&sa=X&ved=0ahUKEwiDwqPeqLHUAhXIKcAKHSmsBW8Q6AEIVzAH#v=onepage&q=mass%20closure%20aerosol%20london&f=false
    # U is capitalised...

    # Need something to represent Iron-rich Dusts (mineral matter?)
    #f_RH_Birmingham = (0.093 * np.array(f_RH['Generic NaCl'])) + (0.08 * np.array(f_RH['Soot'])) + \
    #(0.2757 * np.array(f_RH['Ammonium Sulphate'])) + (0.2757 * np.array(f_RH['Ammonium nitrate'])) + (0.2757 * np.array(f_RH['Aged fossil-fuel OC']))

    # How much bigger is Birmingham to MURK?
    # f_RH['Birmingham'] / np.array(f_RH['MURK'])

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
    plt.plot([0.0, 100.0], [1.0, 1.0], label='Black carbon')

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