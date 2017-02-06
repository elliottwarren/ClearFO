
"""
Create the f(RH) graph and look up table for use in calculating the extinction efficiency (Q)

Created by Elliott 06/02/2017
"""

import numpy as np
# import matplotlib.pyplot as plt

def main():

    # setup
    # find wavelength band
    # read file
    # find blocks
    # format data
    # create f(RH) curve using \sigma_ext(RH=0) from file as \sigma_ext(dry)

    # directories
    datadir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/data/Mie/'
    file_name = 'spec3a_sw_hadgem1_7lean_so' # original file given to me by Claire Ryder 25/01/17
    file_path = datadir + file_name

    ceil_lam = 0.91e-06

    # variables to take from file (as listed within the file) with index from BLOCK = 0
    # NOTE: data MUST be in ascending index order
    aerosols = {'Accum. Sulphate': 6, 'Aitken Sulphate': 7, 'Aged fossil-fuel OC': 22}
    aer_order = ['Accum. Sulphate', 'Aitken Sulphate', 'Aged fossil-fuel OC']


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

    # find the band for current ceilometer wavelength

    file.close()




    data = {}


    for aer in aer_order:

        data[aer] = []

        file = open(file_path, "r")

        line = file.readline()  # read the first actual line
        line = line.rstrip('\n\r')

        while line != '*BLOCK: TYPE =   11: SUBTYPE =    1: VERSION =    2':
            line = file.readline()
            line = line.rstrip('\n\r')
            line = ' '.join(line.split())

        while line != 'Index of species = 6 Accum. Sulphate':
            line = file.readline()
            line = line.rstrip('\n\r')

        while line != 'Band =    4':
            line = file.readline()
            line = line.rstrip('\n\r')

        # skip two lines
        line = file.readline()
        line = file.readline()

        # read first line of data
        line = file.readline()
        line = ' '.join(line.split()) # remove duplicate, trailing and leading spaces
        line_split = line.split(' ')

        # when reached the data line...
        while line != 'Band =    5':

            line = ' '.join(line.split()) # remove duplicate, trailing and leading spaces
            line_split = line.split(' ')

            # extract data
            data[aer] += [line_split]

            line = file.readline()
            line = line.rstrip('\n\r')

        file.close()

        # convert to numpy array, list already alligned correctly for np.array()
        data[aer] = np.array(data[aer])









    # line = line.rstrip('\n\r')
    #
    # # skip until correct block is reached
    # while line != '*BLOCK: TYPE =   11: SUBTYPE =    1: VERSION =    2':
    #     pass
    #
    # for aer in aerosols.iterkeys():
    #
    #     while line != 'Index of species =     ' + block + '  Accum. Sulphate':
    #         line = file.readline()
    #         line = line.rstrip('\n\r')













if __name__ == '__main__':
    main()