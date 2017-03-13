
"""
Create the f(RH) graph and look up table for use in calculating the extinction efficiency (Q)

Created by Elliott 06/02/2017
"""

import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy

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
    :param Q_type:
    :return:
    """

    print "Reading from:" + file_path
    print "Band: " + str(band)

    # define arrays
    data = {}

    for aer in aer_order:

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
        while line != 'Band = ' + str(band + 1):


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

def main():

    # setup
    # find wavelength band
    # read file
    # find blocks
    # format data
    # create f(RH) curve using \sigma_ext(RH=0) from file as \sigma_ext(dry)

    # directories
    savedir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/figures/Mie/f(RH)/'
    datadir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/data/Mie/'
    # file_name = 'spec3a_sw_hadgem1_7lean_so' # original file given to me by Claire Ryder 25/01/17
    file_name = 'sp_sw_ga7'
    file_path = datadir + file_name

    # variables to take from file (as listed within the file) with index from BLOCK = 0
    # NOTE: data MUST be in ascending index order
    if file_name == 'spec3a_sw_hadgem1_7lean_so':
        aer_index = {'Accum. Sulphate': 6, 'Aitken Sulphate': 7, 'Aged fossil-fuel OC': 22}
        aer_order = ['Accum. Sulphate', 'Aitken Sulphate', 'Aged fossil-fuel OC']
    else:
        aer_index = {'Accum. Sulphate': 14, 'Aitken Sulphate': 15, 'Aged fossil-fuel OC': 24, 'Ammonium nitrate': 26}
        aer_order = ['Accum. Sulphate', 'Aitken Sulphate', 'Aged fossil-fuel OC', 'Ammonium nitrate']

    # Q type to use in calculating f(RH)
    Q_type = 'extinction'
    print 'Q_type = ' + Q_type

    # band range to use
    band = 3

    # ---------------------------------------------------
    # Read, Process and save f(RH)
    # ---------------------------------------------------

    # read in the spectral band information
    spec_bands = read_spec_bands(file_path)

    # wavelength range in current band
    band_idx = np.where(spec_bands['band'] == band)[0][0]
    band_lam_range = str(spec_bands['lower_limit'][band_idx] * 1.0e9)[:3] + '-' + \
                     str(spec_bands['upper_limit'][band_idx] * 1.0e9)[:3] + 'nm'



    # read the aerosol data
    data = read_aer_data(file_path, aer_index, aer_order, band=band)

    # Extract RH (RH is the same across all aerosol types)
    RH = np.array(data[aer_order[0]][:, 0])

    # calculate f(RH)
    Q, f_RH = calc_f_RH(data, aer_order, Q_type=Q_type)

    # create an average f(RH)
    f_RH['average with Aitken Sulphate'] = np.mean(f_RH.values(), axis=0)
    f_RH['average'] = np.mean((f_RH['Accum. Sulphate'], f_RH['Aged fossil-fuel OC'], f_RH['Ammonium nitrate']), axis=0)

    # save f(RH)
    # np.savetxt(datadir +  'calculated_ext_f(RH)_'+str(ceil_lam)+'nm.csv', np.transpose(np.vstack((RH, f_RH['average']))), delimiter=',', header='RH,f_RH')
    np.savetxt(datadir +  'calculated_ext_f(RH)_320-690nm.csv', np.transpose(np.vstack((RH, f_RH['average']))), delimiter=',', header='RH,f_RH')
    # ---------------------------------------------------
    # Plotting
    # ---------------------------------------------------

    # plot the data at the end so it is all neat and together
    fig = plt.figure(figsize=(6, 4))
    for key, value in data.iteritems():

        plt.plot(value[:, 0], f_RH[key], label=key, linestyle='--')

    plt.plot(value[:, 0], f_RH['average with Aitken Sulphate'], label='average with Aitken Sulphate', linestyle='-')

    # plot the average one (use RH from the last data.iteritems()
    plt.plot(value[:, 0], f_RH['average'], label='average without Aitken Sulphate', color='black')


    plt.legend(fontsize=9, loc='best')
    plt.xlabel('RH')
    plt.ylabel(Q_type + ' f(RH)')
    plt.ylim([0, 8.0])
    plt.title(file_name + ': 320-690nm band')
    plt.savefig(savedir + file_name + '_' + Q_type[0:3] + '_f_RH_320-690nm.png')
    # plt.title(file_name + ': 690-1190nm band')
    # plt.savefig(savedir + file_name + '_' + Q_type[0:3] + '_f_RH_690-1190nm.png')


if __name__ == '__main__':
    main()