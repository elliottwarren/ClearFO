
import os
import glob
import sys
import string

import numpy as np
import math
import ellUtils as eu

from netCDF4 import Dataset
# from read_nc_file import *
# from write_nc_file import *

import matplotlib.pyplot as plt

def load_blks(blks):

    """
    Load in all the blocks from blk files in a list, then convert the lists to single arrays for each variable

    :return: 'All the vars': p, t, t_degC, mf, ext, mass_frac, mass_frac_gkg, freq, lam
    """

    # declare variables
    abs = []
    freq = []
    lam = []

    for blk in blks:

        print 'water_vapour', blk
        fname = datadir + 'gas_00000001_blk_' + str(blk).zfill(8) + '.nc'
        f = Dataset(fname, 'r')
        p_dim = f.dimensions['pressure']
        freq_dim = f.dimensions['frequency']
        t_dim = f.dimensions['temperature']
        mf_dim = f.dimensions['mass_fraction'] # water vapour mass fraction

        # blk dimensions (quick)
        print blk, len(p_dim), len(freq_dim), len(t_dim), len(mf_dim)

        # variables
        # vars that are the same between files
        p = f.variables['pressure'][:]  # / 100.0 to conv to hPa
        t = f.variables['temperature'][:]  #
        t_degC = t - 273.15
        mass_frac = f.variables['mass_fraction'][:]
        mass_frac_gkg = f.variables['mass_fraction'][:] * 1000.0


        # vars that are different between files
        freq_Hz = f.variables['frequency'][:] # just for calculating lam
        freq.append(f.variables['frequency'][:])  # [Hz]
        abs.append(f.variables['mass_ext_coeff'][:, :, :, :])
        lam.append(c_cms / freq_Hz)

    # convert appended lists to numpy arrays
    freq = np.hstack(freq) # hstack
    abs = np.vstack(abs) # vstack (first dimension)
    lam = np.hstack(lam)

    return p, t, t_degC, abs, mass_frac, mass_frac_gkg, freq, lam

def normpdf(x, mean, sd):

    """
    Gives the probability of y, given a value of x, assuming a gaussian distribution
    :param x:
    :param mean:
    :param sd:
    :return:

    unit of p(y) is a fraction
    """

    var = float(sd) ** 2
    pi = 3.1415926
    denom = (2 * pi * var) ** .5
    num = math.exp(-(float(x) - float(mean)) ** 2 / (2 * var))

    return num / denom

def gaussian(x, mean, var):

    return (np.exp((-0.5*(np.asarray(x)-mean)**2)/var) /
            math.sqrt(2*math.pi*var))

# c_cms=2.99792e10
c_cms=299792458 # speed of light

# def water_vapour():

maindir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/'
datadir = maindir + 'data/water_vapour/'
savedir = maindir + 'figures/water_vapour/'

# load in variables from blocks
blks = range(467, 477) # blks to load # 905 nm in the middle of blk 471
p, t, t_degC, abs, mass_frac, mass_frac_gkg, freq, lam = load_blks(blks)
lam_nm = lam*1e09 # nm version

# abs shape = (frequency, mass fraction, temperature, pressure)
# tau = abs mass frac * water vapour mixing ratio * air density * depth

# idx positions for variable values used in FO
# q_ind=13 # 3%  mass fraction water vapour / air
# t_ind=13 # 280 K
q_ind=12 # 1%  mass fraction water vapour / air
t_ind=14 # 290 K / 16.85 deg C
p_ind=49 # 1100 hPa (closest to 1013 hPa that is in the file, but there is little sensitivity to p anyway)

# guassian weight
#mean = 905e-09 # nm
#sigma = 4.0
#var = 16.0


# Calculate water vapour mass absorption for a mean wavelength of 895 - 915 nm with FWHM of 4 nm.

# FWHM = Full Width Half Maximum - CL31 have FWHM of 4 (+/-2 around 905 nm)
FWHM = 4.0e-09
# sigma backcalculated from FWHM: https://en.wikipedia.org/wiki/Full_width_at_half_maximum - based on wolfram alpha
sigma = FWHM /(2 * np.sqrt(2 * np.log(2)))

# central wavelength the CL31 could be at
lam_cent = np.arange(895e-09, 915e-09, 1e-09)

# store water vapour mass absorption for the different wavelengths [np.array]
# store gaussian weights for absra plot [list]
wv_gaussians = []
wv_weighted_abs = np.empty(len(lam_cent))
wv_weighted_abs[:] = np.nan

for lam_idx, lam_i in zip(range(len(lam_cent)), lam_cent):

    # calculate gaussian across the whole set of blks for lam_i
    weights = np.array([normpdf(i, lam_i, sigma) for i in lam])

    # modify weights so they add up to 1 - append value to gaussian list
    # convex combination - https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Convex_combination_example
    wv_gaussians_i = weights/np.sum(weights)
    wv_gaussians += [wv_gaussians_i]
    # wv_gaussians += [weights] # alternative to above to quickly store the unormalised weights


    # what each value is worth (mean will have the highest weight but not 1, as all the weights need to add up to 1 in the end)
    # store value in numpy array of predefined size
    wv_weighted_abs[lam_idx] = np.sum(wv_gaussians_i * abs[:, q_ind, t_ind, p_ind])


# Plotting

n_plots=1
fig=plt.figure(1,figsize=(6,4))
ax=plt.subplot(1,1,1)
# line1=ax.plot(1.0e7/(freq[:]/c_cms),abs[:,q_ind,t_ind,p_ind], label = str(p[p_ind]))

for i in range(len(mass_frac)):

    # ax.plot(c_cms / freq[:] * 1e9, abs[:, i, t_ind, p_ind], label=str(mass_frac_gkg[i]))
    #
    # ax.plot([c_cms /freq[0]* 1e9, c_cms /freq[-1]* 1e9], [np.mean(abs[:, i, t_ind, p_ind]),np.mean(abs[:, i, t_ind, p_ind])])
    #         #,label='avg ' + str(mass_frac_gkg[i]))

    ax.plot(c_cms / freq[:] * 1e9, abs[:, i, t_ind, p_ind], label=str(mass_frac_gkg[i]))


ax.set_xlabel(r"Wavelength [nm]",fontsize=14, labelpad=0.1)
ax.set_ylabel(r"Mass absorption [m2kg-1]",fontsize=14)
ax.tick_params(labelsize=14)
# plt.ylim([0.0, 0.3])
plt.suptitle('r [g kg-1]; 1g avg = ' + str(np.mean(abs[:, -5, t_ind, p_ind])) + '; 100 g avg = ' +
             str(np.mean(abs[:, -1, t_ind, p_ind])))
plt.legend(fontsize=8)
plt.grid()

if type(blks) == list:
    if len(blks) == 1:
        plt.savefig(savedir + 'r_variable_blk' + str(blks[0]) + '.png') # single blks
    elif len(blks) > 1:
        plt.savefig(savedir + 'r_variable_blks' + str(blks[0]) + '-' + str(blks[-1]) + '.png')

plt.close(fig)

# ------------------

# plot change in water vapour mass absorption wrt central wavelength
fig=plt.figure(1,figsize=(6,4))
ax=plt.subplot(1,1,1)

ax.plot(lam_cent * 1e9, wv_weighted_abs)
ax.set_ylabel('wv absorption [m2 kg-1]')
ax.set_xlabel('wavelength [nm]')
plt.savefig(savedir + 'wv_mass_abs_wrt_centralLam.png')

plt.close(fig)
# ------------------------------------

# plot the guassians
fig=plt.figure(1,figsize=(6,4))
ax=plt.subplot(1,1,1)

for gaus_i, lam_cent_i in zip(wv_gaussians, lam_cent):

    ax.plot(lam * 1e9, gaus_i, label=str(lam_cent_i * 1e9))


ax.set_ylabel('weighting')
ax.set_xlabel('wavelength [nm]')
ax.legend()
plt.savefig(savedir + 'wv_mass_abs_gaussians_wrt_centralLam.png')

plt.close(fig)