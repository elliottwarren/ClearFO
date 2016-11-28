"""
Script to make sensitivity plots of the different FO parameters

"""

import matplotlib.pyplot as plt
import numpy as np
import datetime as dt

from forward_operator import FOconstants as FOcon
from forward_operator import FOUtils as FO

# -----------------------------------------------------------
# Constants
# -----------------------------------------------------------

savedir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/figures/sensitivity/'

# density of aerosol (for ammonium sulphate)
rho_a = FOcon.rho_a

r0 = FOcon.r0_haywood  # radius of a 'standard' aerosol particle [m]

# (p) power used to represent the variation in aerosol particle size with mixing ratio
p = FOcon.p_aer

# B = cutils.B_activation_clark  # (B) Activation parameter (for ammonium sulphate)
B = FOcon.B_activation_haywood

# Lidar ratio [sr-1]
S = FOcon.LidarRatio['Aerosol']

# Update from Haywood et.al. (2008) using flight data from around UK
# this is the set of parameters used in the UM v8.1 (see UMDP 26)
# to compute the aerosol number N_aer in code below
N0 = FOcon.N0_aer  # (N0) standard number density of aerosol
m0 = FOcon.m0_aer  # (m0) standard mass mixing ratio of the aerosol

# For use in  eqns. 17-18 in Clark et.al. (2008):
eta = FOcon.eta
Q = FOcon.Q_ext_aer

# --------------------------------------------------

# range of murk aerosol [micrograms kg-1]
q_aer = np.arange(0, 301)

# convert micrograms kg-1 to kg/kg
q_aer_kg_kg = q_aer * 1.0e-9  # convert micrograms kg-1 to kg/kg

# relative humidity
RH = np.arange(0,101)

# ---------------------------------------------------

# Variables for model data
day = dt.datetime(2016, 05, 04)
model_type = 'UKV'
res = FOcon.model_resolution[model_type]
modDatadir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/data/' + model_type + '/'

# get model data
mod_data = FO.mod_site_extract_calc(day, {'IMU': [-0.106086, 51.52615]}, modDatadir, model_type, res)

# Take heights
z_mod = mod_data['IMU']['level_height']

# -----------------------------------------------------------
# 1. r_md
# -----------------------------------------------------------

# r_md is the dry mean volume radius. Eqn. 2 in Clark et.al. (2008)
r_md = r0 * np.power((q_aer_kg_kg / m0), p)
# r_md2 = r0 * np.power((q_aer_kg_kg / m0), 0.5)

# set up plot
fig = plt.figure(figsize=(6, 4))
ax = plt.subplot2grid((1, 1), (0, 0))

ax.plot(q_aer, r_md, label = 'p = 1/6')
# ax.plot(q_aer, r_md2, label = 'p = 1/2')

ax.set_ylabel('r_md [kg kg-1]')
ax.set_xlabel('murk aerosol [microgram kg-1]')
# plt.legend()
plt.tight_layout()

plt.savefig(savedir + 'murk_rmd.png')  # filename
plt.close(fig)

# -----------------------------------------------------------
# 2. N_0
# -----------------------------------------------------------

# Compute the aerosol number density N_aer. Eqn. 3 in Clark et.al. (2008) and Eqn 39 in UMDP 26 Large-scale precip.
N_aer = N0 * np.power((q_aer_kg_kg / m0), 1-(3*p))

# set up plot
fig = plt.figure(figsize=(6, 4))
ax = plt.subplot2grid((1, 1), (0, 0))

ax.plot(q_aer, N_aer)
ax.set_ylabel('N_0 [kg-1]')
ax.set_xlabel('murk aerosol [microgram kg-1]')
plt.tight_layout()

plt.savefig(savedir + 'murk_N0.png')  # filename
plt.close(fig)




# # -----------------------------------------------------------
# # 1. Aerosol extinction scattering
# # -----------------------------------------------------------
#
# N_aer = N0 * np.power((q_aer_kg_kg / m0), 1 - (3 * p))
#
# # r_md is the dry mean volume radius. Eqn. 2 in Clark et.al. (2008)
# r_md = r0 * np.power((q_aer_kg_kg / m0), p)
#
# # rm is the mean volume radius. Eqn. 12 in Clark et.al. (2008)
# # "When no activated particles are present an analytic solution for rm" is
# RH_crit = FOcon.RH_crit
# #    print "RH_crit = ", RH_crit
# #    print RH[1:70, 1:4]
# # mask over values less than critical
# RH_ge_RHcrit = np.ma.masked_less(RH, RH_crit)
#
# RH_factor = 0.01  # Relative Humidity in 0.38 not 38%
#
# # eq 12 - calc rm for RH greater than critical
# rm = np.ma.ones(RH.shape) - (B / np.ma.log(RH_factor * RH_ge_RHcrit))
# rm2 = np.ma.power(rm, 1. / 3.)
# rm = np.ma.array(r_md) * rm2
#
# # set rm as 0 where RH is less than crit
# rm = np.ma.MaskedArray.filled(rm, [0.0])
# where_lt_crit = np.where(np.logical_or(RH.data < RH_crit, rm == 0.0))
# rm[where_lt_crit] = r_md[where_lt_crit]
#
# # Close to activation one must solve the full equation (Kohler curve), not done in this code.
# # Assumptions made here include:
# # 1. Q_ext = scattering efficiency is independent of particle size and is assumed to be on average = 2.0
# # 2. Only time RH is taken into account is in equation 12 above in the
# # calculation of a mean radius. No attempt is made to model the variation
# # of particle growth with RH AS A FUNCTION OF SIZE of particle.
#
# # Calculate extinction coefficient
# # eqns. 17-18 in Clark et.al. (2008)
# alpha_a = (eta * Q) * np.pi * N_aer * np.power(rm, 2)
#
# # Calculate backscatter using a constant lidar ratio
# beta_unAtt = alpha_a / S
#
# # -----------------------------------------------------------
# # 2. Transmission
# # -----------------------------------------------------------
#
# # create alpha and beta coefficients for aerosol
#
# dz = np.zeros_like(z_mod)
# dz[0] = z_mod[0]
# dz[1:len(z_mod)] = z_mod[1:len(z_mod)] - z_mod[0:len(z_mod) - 1]
#
# ss = np.shape(alpha_a)
# # print "1. compute_transmission: size = ", ss[0],ss[1]
#
# integral_alpha = np.empty_like(alpha_a)
# transmission = np.empty_like(alpha_a)
# # print np.shape(integral_alpha), np.shape(transmission)
#
# # print  "2. compute_transmission: ",
# # np.shape(alpha_extinction),np.shape(dz)
#
# integral_alpha = np.cumsum(alpha_a * dz[0:len(dz)], axis=0)
#
# optical_depth = integral_alpha
# # print  "3. compute_transmission: ",alpha_extinction*dz[0:len(dz)],integral_alpha
# # transmission =  integral_alpha
# transmission = np.exp(-2.0 * integral_alpha)
#
# # -----------------------------------------------------------
# # 3. Attenuated backscatter
# # -----------------------------------------------------------
#
# # int_mod_alpha, mod_transm = compute_transmission(mod_alpha[:, 0:NN], dz)
#
# # derive modelled attenuated backscatter
# # bsc_mod = np.log10(mod_transm * mod_bscUnnAtt)
# beta_Att = transmission * beta_unAtt








print 'END PROGRAM'