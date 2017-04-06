"""
Create 4 panel plot with f(RH) differences for each of the aerosol species and the average f(RH)

Created by Elliott 06/04/17
"""

import numpy as np
import matplotlib.pyplot as plt
import f(RH)_creation




def main():

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
    file_name = 'sp_ew_ceil_guass_903-907'
    file_path = specdir + file_name

    # variables to take from file (as listed within the file) with index from BLOCK = 0
    # NOTE: data MUST be in ascending index order
    if file_name == 'spec3a_sw_hadgem1_7lean_so':
        aer_index = {'Accum. Sulphate': 6, 'Aitken Sulphate': 7, 'Aged fossil-fuel OC': 22}
        aer_order = ['Accum. Sulphate', 'Aitken Sulphate', 'Aged fossil-fuel OC']
    else:
        aer_index = {'Accum. Sulphate': 1, 'Aged fossil-fuel OC': 2, 'Ammonium nitrate': 3}
        aer_order = ['Accum. Sulphate', 'Aged fossil-fuel OC', 'Ammonium nitrate']

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

    # create an average f(RH)
    # f_RH['average with Aitken Sulphate'] = np.mean(f_RH.values(), axis=0)
    f_RH['average'] = np.mean((f_RH['Accum. Sulphate'], f_RH['Aged fossil-fuel OC'], f_RH['Ammonium nitrate']), axis=0)
