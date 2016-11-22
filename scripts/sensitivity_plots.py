"""
Script to make sensitivity plots of the different FO parameters

"""

import matplotlib.pyplot as plt

from forward_operator import FOconstants as FOcon
from forward_operator import FOUtils as FO

# -----------------------------------------------------------
# Constants
# -----------------------------------------------------------
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

# Compute the aerosol number density N_aer. Eqn. 3 in Clark et.al. (2008) and Eqn 39 in UMDP 26 Large-scale precip.
q_aer_kg_kg = q_aer * 1.0e-9  # convert micrograms kg-1 to kg/kg

# -----------------------------------------------------------
# 1. r_md
# -----------------------------------------------------------


N_aer = N0 * np.power((q_aer_kg_kg / m0), 0.5)

# r_md is the dry mean volume radius. Eqn. 2 in Clark et.al. (2008)
r_md = r0 * np.power((q_aer_kg_kg / m0), p)

# rm is the mean volume radius. Eqn. 12 in Clark et.al. (2008)
# "When no activated particles are present an analytic solution for rm" is
RH_crit = FOcon.RH_crit
#    print "RH_crit = ", RH_crit
#    print RH[1:70, 1:4]
# mask over values less than critical
RH_ge_RHcrit = np.ma.masked_less(RH, RH_crit)

RH_factor = 0.01  # Relative Humidity in 0.38 not 38%

# eq 12 - calc rm for RH greater than critical
rm = np.ma.ones(RH.shape) - (B / np.ma.log(RH_factor * RH_ge_RHcrit))
rm2 = np.ma.power(rm, 1. / 3.)
rm = np.ma.array(r_md) * rm2

# set rm as 0 where RH is less than crit
rm = np.ma.MaskedArray.filled(rm, [0.0])
where_lt_crit = np.where(np.logical_or(RH.data < RH_crit, rm == 0.0))
rm[where_lt_crit] = r_md[where_lt_crit]

# Close to activation one must solve the full equation (Kohler curve), not done in this code.
# Assumptions made here include:
# 1. Q_ext = scattering efficiency is independent of particle size and is assumed to be on average = 2.0
# 2. Only time RH is taken into account is in equation 12 above in the
# calculation of a mean radius. No attempt is made to model the variation
# of particle growth with RH AS A FUNCTION OF SIZE of particle.

# Calculate extinction coefficient
# eqns. 17-18 in Clark et.al. (2008)
alpha_a = (eta * Q) * np.pi * N_aer * np.power(rm, 2)

# Calculate backscatter using a constant lidar ratio
beta_a = alpha_a / S
