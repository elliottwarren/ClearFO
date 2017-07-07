
import os
import glob
import sys
import string

import numpy as np

from netCDF4 import Dataset
# from read_nc_file import *
# from write_nc_file import *

import matplotlib.pyplot as plt

# c_cms=2.99792e10
c_cms=299792458

# def water_vapour():

maindir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/'
datadir = maindir + 'data/water_vapour/'
savedir = maindir + 'figures/water_vapour/'

blk=471
print 'water_vapour',blk
fname= maindir + 'gas_00000001_blk_'+str(blk).zfill(8)+'.nc'
f=Dataset(fname,'r')
p_dim=f.dimensions['pressure']
r_dim=f.dimensions['frequency']
t_dim=f.dimensions['temperature']
mf_dim=f.dimensions['mass_fraction']
print blk,len(p_dim),len(r_dim),len(t_dim),len(mf_dim)

# variables
freq=f.variables['frequency'][:] # [Hz]
p=f.variables['pressure'][:] # / 100.0 to conv to hPa
t=f.variables['temperature'][:] #
t_degC = t - 273.15
mf=f.variables['mass_fraction'][:] # [kg kg-1]

# freq, mass_frac (15,1e-8-1e-1), t (20,150-340), p (50,1.5-110000.)
ext=f.variables['mass_ext_coeff'][:,:,:,:]
mass_frac=f.variables['mass_fraction'][:]
mass_frac_gkg=f.variables['mass_fraction'][:] * 1000.0
lam = c_cms/freq
# ext shape = (frequency, mass fraction, temperature, pressure)

#tau = ext mass frac * water vapour mixing ratio * air density * depth

# return freq, ext

# freq, ext = water_vapour()

q_ind=13 # 3%  mass fraction water vapour / air
t_ind=13 # 280 K
p_ind=49 # 1100 hPa (slightly on the higher side)

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
ax.set_xlabel(r"Wavelength [nm])",fontsize=14, labelpad=0.1)
ax.set_ylabel(r"Mass absorption [m2kg-1])",fontsize=14)
ax.tick_params(labelsize=14)
# plt.ylim([0.0, 0.3])
plt.suptitle('r [g kg-1]; 1g avg = ' + str(np.mean(ext[:, -5, t_ind, p_ind])) + '; 100 g avg = ' +
             str(np.mean(ext[:, -1, t_ind, p_ind])))
plt.legend(fontsize=8)
plt.grid()
plt.savefig(savedir + 'r_variable_blk' + str(blk) + '.png')
