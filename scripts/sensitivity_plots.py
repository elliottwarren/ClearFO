"""
Script to make sensitivity plots of the different FO parameters

"""

import matplotlib.pyplot as plt
import numpy as np

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

# range of murk aerosol [micrograms kg-1]
q_aer = np.arange(0, 301)

# convert micrograms kg-1 to kg/kg
q_aer_kg_kg = q_aer * 1.0e-9  # convert micrograms kg-1 to kg/kg

# -----------------------------------------------------------
# 1. r_md
# -----------------------------------------------------------

# r_md is the dry mean volume radius. Eqn. 2 in Clark et.al. (2008)
r_md = r0 * np.power((q_aer_kg_kg / m0), p)

# set up plot
fig = plt.figure(figsize=(6, 4))
ax = plt.subplot2grid((1, 1), (0, 0))

ax.plot(q_aer, r_md)
ax.set_ylabel('r_md [kg kg-1]')
ax.set_xlabel('murk aerosol [microgram kg-1]')
plt.tight_layout()

plt.savefig(savedir + 'murk_rmd.png')  # filename
plt.close(fig)

# -----------------------------------------------------------
# 2. N_0
# -----------------------------------------------------------

# Compute the aerosol number density N_aer. Eqn. 3 in Clark et.al. (2008) and Eqn 39 in UMDP 26 Large-scale precip.
N_aer = N0 * np.power((q_aer_kg_kg / m0), 0.5)

# set up plot
fig = plt.figure(figsize=(6, 4))
ax = plt.subplot2grid((1, 1), (0, 0))

ax.plot(q_aer, N_aer)
ax.set_ylabel('N_0 [kg-1]')
ax.set_xlabel('murk aerosol [microgram kg-1]')
plt.tight_layout()

plt.savefig(savedir + 'murk_N0.png')  # filename
plt.close(fig)



