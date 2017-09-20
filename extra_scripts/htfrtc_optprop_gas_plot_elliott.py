
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
    ext = []
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
        ext.append(f.variables['mass_ext_coeff'][:, :, :, :])
        lam.append(c_cms / freq_Hz)

    # convert appended lists to numpy arrays
    freq = np.hstack(freq) # hstack
    ext = np.vstack(ext) # vstack (first dimension)
    lam = np.hstack(lam)

    return p, t, t_degC, ext, mass_frac, mass_frac_gkg, freq, lam


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
blks = [470, 471, 472] # blks to load
p, t, t_degC, ext, mass_frac, mass_frac_gkg, freq, lam = load_blks(blks)
lam_nm = lam*1e09 # nm version


# ext shape = (frequency, mass fraction, temperature, pressure)
# tau = ext mass frac * water vapour mixing ratio * air density * depth

# idx positions for variable values used in FO
# q_ind=13 # 3%  mass fraction water vapour / air
# t_ind=13 # 280 K
q_ind=12 # 1%  mass fraction water vapour / air
t_ind=14 # 290 K / 16.85 deg C
p_ind=49 # 1100 hPa (closest to 1013 hPa that is in the file, but there is little sensitivity to p anyway)

# guassian weight
mean = 905e-09 # nm
sigma = 4.0
var = 16.0

# --------------- FWHM = 4 nm -------------------#

# FWHM = Full Width Half Maximum - CL31 have FWHM of 4 (+/-2 around 905 nm)
FWHM = 4.0e-09
# sigma backcalculated from FWHM: https://en.wikipedia.org/wiki/Full_width_at_half_maximum - based on wolfram alpha
sigma = FWHM /(2 * np.sqrt(2 * np.log(2)))

weights = np.array([normpdf(i, 905e-09, sigma) for i in lam])
# weights = weights/weights.max()
plt.plot(lam, weights)

# slow function, could wait to multiply by weights once idx extraction for plotting is done.
weighted_ext = [weights[i] * ext[i, :, :, :] for i in range(ext.shape[0])]

# create avergae using the gaussian weighted values
gaus_weight_avg = np.mean(weighted_ext,axis=0)
gaus_weight_sum = np.sum(weighted_ext,axis=0)

# gaussian weighted average
gaus_weight_avg[q_ind, t_ind, p_ind]
gaus_weight_sum[q_ind, t_ind, p_ind]

# --------------- FWHM = 8 nm -------------------#

# FWHM = Full Width Half Maximum - CL31 have FWHM of 4 (+/-2 around 905 nm)
FWHM = 8.0e-09
# sigma backcalculated from FWHM: https://en.wikipedia.org/wiki/Full_width_at_half_maximum - based on wolfram alpha
sigma = FWHM /(2 * np.sqrt(2 * np.log(2)))

weights = np.array([normpdf(i, 905e-09, sigma) for i in lam])
weights = weights/weights.max()
plt.plot(lam, weights)

# slow function, could wait to multiply by weights once idx extraction for plotting is done.
weighted_ext = [weights[i] * ext[i, :, :, :] for i in range(ext.shape[0])]

# create avergae using the gaussian weighted values
gaus_weight_avg = np.mean(weighted_ext,axis=0)

# gaussian weighted average
gaus_weight_avg[q_ind, t_ind, p_ind]

# --------------- FWHM = 8 nm -------------------#

# FWHM = Full Width Half Maximum - CL31 have FWHM of 4 (+/-2 around 905 nm)
FWHM = 4.0e-09
# sigma backcalculated from FWHM: https://en.wikipedia.org/wiki/Full_width_at_half_maximum - based on wolfram alpha
sigma = FWHM /(2 * np.sqrt(2 * np.log(2)))

weights = np.array([normpdf(i, 905e-09, sigma) for i in lam])

#all weights add up to 1
true_weights = weights/np.sum(weights)
# weights = np.ones(len(lam)) # use with weighted_ext = np.sum... to check that weighted_ext then equals the average.
# weights = weights/weights.max()

# what each value is worth (mean will have the highest weight but not 1, as all the weights need to add up to 1 in the end)
weighted_ext_true = true_weights * ext[:, q_ind, t_ind, p_ind]

weighted_average = np.sum(weighted_ext_true)

plt.plot(lam, weights)

# slow function, could wait to multiply by weights once idx extraction for plotting is done.
# weighted_ext = weights * ext[:, q_ind, t_ind, p_ind]
# from Weigner et al 2015 (eq 15) - short way to do the weighted average
weighted_ext = np.sum(weights * ext[:, q_ind, t_ind, p_ind])/np.sum(weights)

# create avergae using the gaussian weighted values
# gaus_weight_avg = np.mean(weighted_ext)

# gaussian weighted average
# gaus_weight_avg[q_ind, t_ind, p_ind]


# Plotting

n_plots=1
fig1=plt.figure(1,figsize=(6,4))
ax=plt.subplot(1,1,1)
# line1=ax.plot(1.0e7/(freq[:]/c_cms),ext[:,q_ind,t_ind,p_ind], label = str(p[p_ind]))

for i in range(len(mass_frac)):

    # ax.plot(c_cms / freq[:] * 1e9, ext[:, i, t_ind, p_ind], label=str(mass_frac_gkg[i]))
    #
    # ax.plot([c_cms /freq[0]* 1e9, c_cms /freq[-1]* 1e9], [np.mean(ext[:, i, t_ind, p_ind]),np.mean(ext[:, i, t_ind, p_ind])])
    #         #,label='avg ' + str(mass_frac_gkg[i]))

    ax.plot(c_cms / freq[:] * 1e9, ext[:, i, t_ind, p_ind], label=str(mass_frac_gkg[i]))


# ax.plot(c_cms/freq[:] * 1e9,ext[:,q_ind,t_ind,p_ind], label = str(t_degC[t_ind]))
# ax.plot(c_cms/freq[:] * 1e9,ext[:,q_ind,t_ind,p_ind], label = str(t_degC[t_ind]))
# ax.plot(c_cms/freq[:] * 1e9,ext[:,q_ind,t_ind,p_ind], label = str(t_degC[t_ind]))
# ax.plot(c_cms/freq[:] * 1e9,ext[:,q_ind,t_ind,p_ind], label = str(t_degC[t_ind]))
# ax.plot(c_cms/freq[:] * 1e9,ext[:,q_ind,t_ind,p_ind], label = str(t_degC[t_ind]))
#
# ax.axhline(np.mean(ext[:,q_ind,t_ind-1,p_ind]), label = 'avg ' + str(t_degC[t_ind]), color='red')
# ax.axhline(np.mean(ext[:,q_ind,t_ind,p_ind]), label = 'avg ' + str(t_degC[t_ind]), color = 'black')
# ax.axhline(np.mean(ext[:,q_ind,t_ind+1,p_ind]), label = 'avg ' + str(t_degC[t_ind+1]), color='cyan')
# ax.axhline(np.mean(ext[:,q_ind,t_ind+2,p_ind]), label = 'avg ' + str(t_degC[t_ind+2]), color='cyan')
# ax.axhline(np.mean(ext[:,q_ind,t_ind+3,p_ind]), label = 'avg ' + str(t_degC[t_ind+3]), color='cyan')
# diff = ext[:,q_ind,t_ind,p_ind] - ext[:,q_ind,t_ind,p_ind-1]
# line2=ax.plot(c_cms/freq[:] * 1e9,diff, label = str(p[p_ind-1]))
ax.set_xlabel(r"Wavelength [nm]",fontsize=14, labelpad=0.1)
ax.set_ylabel(r"Mass absorption [m2kg-1]",fontsize=14)
ax.tick_params(labelsize=14)
# plt.ylim([0.0, 0.3])
plt.suptitle('r [g kg-1]; 1g avg = ' + str(np.mean(ext[:, -5, t_ind, p_ind])) + '; 100 g avg = ' +
             str(np.mean(ext[:, -1, t_ind, p_ind])))
plt.legend(fontsize=8)
plt.grid()

if type(blks) == list:
    if len(blks) == 1:
        plt.savefig(savedir + 'r_variable_blk' + str(blks[0]) + '.png') # single blks
    elif len(blks) > 1:
        plt.savefig(savedir + 'r_variable_blks' + str(blks[0]) + '-' + str(blks[-1]) + '.png')
