"""
Script to show sensitivity of beta to m and RH

Created by Elliott 26/05/17

"""

import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
from copy import deepcopy

from forward_operator import FOconstants as FOcon
from forward_operator import FOUtils as FO

def create_range(value):

    step = (value[1] - value[0]) / 100.0
    var_range = np.arange(value[0], value[1] + step, step)

    return var_range

# calculate vars

def calc_Q_ext_wet(ceil_lam, r_md, RH):
    """
    Calculate Q_ext_wet using Q_ext_dry and f(RH) for current wavelength
    EW 23/02/17
    :param ceil_lam:
    :param r_md:
    :param RH:
    :return:
    """

    from ellUtils import nearest

    def read_f_RH(ceil_lam):
        """
        Read in the f_RH data from csv
        EW 21/02/17

        :param filename:
        :return: data = {RH:... f_RH:...}

        filename must be in the form of 'calculated_ext_f(RH)_[ceil_lambda]nm.csv'
        """

        # temp file name
        miedir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/data/Mie/'
        # filename = 'calculated_ext_f(RH)_' + str(ceil_lam) + 'nm.csv'
        filename = 'sp_ew_ceil_guass_908-912_ext_f(RH)_908-912nm.csv'

        # read data
        raw = np.loadtxt(miedir + filename, delimiter=',')

        f_RH = {'RH': raw[:, 0],
                'f_RH': raw[:, 1]}

        return f_RH

    def read_Q_dry_ext(ceil_lam):
        """
        Read in the Q_ext for dry murk.
        EW 21/02/17

        :param filename:
        :param lam:
        :return: Q_ext_dry = {radius:... Q_ext_dry:...}

        Requres the wavelength to be passed, just so in the future, the 910 nm file is not incorrectly used by mistake when
        it should use the file for another wavelength.
        """

        miedir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/data/Mie/'
        filename = 'calculated_Q_ext_' + str(ceil_lam) + 'nm.csv'

        raw = np.loadtxt(miedir + filename, delimiter=',')

        Q_ext_dry = {'radius': raw[:, 0],
                     'Q_ext': raw[:, 1]}

        return Q_ext_dry

    RH_factor = 0.01  # Relative Humidity in 0.38 not 38%

    # calculate Q_ext_wet
    f_RH = read_f_RH(ceil_lam)
    Q_ext_dry = read_Q_dry_ext(ceil_lam)

    # create matric of Q_ext_dry based on r_md
    Q_ext_dry_matrix = np.empty(r_md.shape)
    f_RH_matrix = np.empty(RH.shape)

    # find Q_ext dry, given the dry radius matrix
    if r_md.size != 1:
        for i in range(r_md.shape[0]):
            idx = nearest(Q_ext_dry['radius'], r_md[i])[1]
            Q_ext_dry_matrix[i] = Q_ext_dry['Q_ext'][idx]

    else:
        idx = nearest(Q_ext_dry['radius'], r_md)[1]
        Q_ext_dry_matrix = Q_ext_dry['Q_ext'][idx]

    # find f(RH), given the RH matrix
    # need RH factor as f_RH['RH'] in units of frac not percentage
    if RH.size != 1:
        for i in range(RH.shape[0]):
            idx = nearest(f_RH['RH'], RH_factor * RH[i])[1]
            f_RH_matrix[i] = f_RH['f_RH'][idx]
    else:
        idx = nearest(f_RH['RH'], RH_factor * RH)[1]
        f_RH_matrix = f_RH['f_RH'][idx]

    # calculate Q_ext_wet
    Q = Q_ext_dry_matrix * f_RH_matrix
    # print np.mean(Q_ext_dry_matrix[:,:20])

    return Q, Q_ext_dry_matrix, f_RH_matrix

def aer_ext_scat(q_aer, RH, r0 = FOcon.r0_haywood, p = FOcon.p_aer, S = FOcon.LidarRatio['Aerosol'],
                 N0=FOcon.N0_aer, m0 = FOcon.m0_aer, eta = FOcon.eta, r_md = []):
                 #, lnRH = [], BolnRH = []):

    """
    Compute aerosol extinction coefficient

    :param q_aer: aerosol mass mizing ratio [micrograms kg-1]
    :param RH: relative humidity [ratio/dimensionless]
    :param r0:
    :param B:
    :return: alpha_a: aerosol extinction coefficient
    :return: beta_a: UNattenuated backscatter
    """

    # Compute the aerosol number density N_aer. Eqn. 3 in Clark et.al. (2008) and Eqn 39 in UMDP 26 Large-scale precip.
    q_aer_kg_kg = q_aer * 1.0e-9  # convert micrograms kg-1 to kg/kg
    N_aer = N0 * np.power((q_aer_kg_kg / m0), 1-(3*p))

    # r_md is the dry mean volume radius. Eqn. 2 in Clark et.al. (2008)
    if type(r_md) == list:
        r_md = r0 * np.power((q_aer_kg_kg / m0), p)

    # rm is the mean volume radius. Eqn. 12 in Clark et.al. (2008)
    # "When no activated particles are present an analytic solution for rm" is
    RH_crit = FOcon.RH_crit

    RH_factor = 0.01  # Relative Humidity in 0.38 not 38%

    # # eq 12 - calc rm for RH greater than critical
    # # Old - original approach. No longer used with f(RH)
    # if RH.shape[0] == 1:
    #     if RH >= RH_crit:
    #         rm = np.ones(RH.shape) - (B / np.log(RH_factor * RH))
    #         rm2 = np.power(rm, 1. / 3.)
    #         rm = np.array(r_md) * rm2
    #     else:
    #         rm = r_md
    #
    # else: # RH shape > 1
    #     # set up rm array
    #     rm = np.empty(RH.shape)
    #     rm[:] = np.nan
    #
    #     # find instances above and below critical threshold
    #     gt_idx = np.where(RH >= RH_crit)[0]
    #     ls_idx = np.where(RH < RH_crit)[0]
    #
    #     # swell those above the threshold, keep rm as r_md where below.
    #     if gt_idx.size: # if not empty
    #         rm[gt_idx] = np.ones(gt_idx.shape) - (B / np.log(RH_factor * RH[gt_idx]))
    #         rm2 = np.power(rm[gt_idx], 1. / 3.)
    #         rm[gt_idx] = np.array(r_md) * rm2
    #     rm[ls_idx] = r_md

    # Close to activation one must solve the full equation (Kohler curve), not done in this code.
    # Assumptions made here include:
    # 1. Q_ext = scattering efficiency is independent of particle size and is assumed to be on average = 2.0
    # 2. Only time RH is taken into account is in equation 12 above in the
    # calculation of a mean radius. No attempt is made to model the variation
    # of particle growth with RH AS A FUNCTION OF SIZE of particle.

    # Calculate Q
    Q, Q_ext_dry_matrix, f_RH_matrix = calc_Q_ext_wet(910, r_md, RH)

    # Calculate extinction coefficient
    # eqns. 17-18 in Clark et.al. (2008)
    alpha_a = (eta * Q) * np.pi * N_aer * np.power(r_md, 2)

    # Calculate backscatter using a constant lidar ratio
    beta_a = alpha_a / S

    # put all the variables into a dinctionary, to be returned at the end of the function
    out = {'beta': beta_a,
           'alpha': alpha_a,
           'r_m': rm,
           'r_md': r_md,
           'N': N_aer,
           'Q': Q,
           'Q_ext_dry': Q_ext_dry_matrix,
           'f_RH': f_RH_matrix}

    return out


# -----------------------------------------------------------
# Constants
# -----------------------------------------------------------

savedir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/figures/sensitivity/'

# Variables for model data
day = dt.datetime(2016, 05, 04)
model_type = 'UKV'
res = FOcon.model_resolution[model_type]
modDatadir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/data/' + model_type + '/'

# ceilometer wavelength [nm]
ceil_lam = 910

# get zmod
mod_data = FO.mod_site_extract_calc(day, {'IMU': [-0.106086, 51.52615]}, modDatadir, model_type, res, ceil_lam)
z_mod = mod_data['IMU']['level_height']

# Conversion factors
micro_g_2_kg = 1.0e-9
kg_2_micro_g = 1.0e9
m_2_microns = 1.0e6
m3_2_cm3 = 1.0e-06

# fixed murk mg kg-1
q_aer_fix = np.array([10])

# create range of murk and RH
m_range = np.arange(0.0, 81.0, 1.0)
RH_range = np.array([20.0, 40.0, 60.0, 80.0, 90.0, 95.0])

# -----------------------------------------------------------
# 2. Calculating
# -----------------------------------------------------------

# Beta and Alpha
# ----------------
# do first to get the var_range out for later calculations (r_md)
# alpha, beta, alpha_c, beta_c, var_range = calc_aer_ext_vars(var_order, variables, q_aer_fix, RH_fix)

beta = {}

for RH_i in RH_range:

    # create a numpy array of the length of m_range full of RH_i to use as input to aer_ext_scat
    RH_i_array = np.empty(len(m_range))
    RH_i_array[:] = RH_i

    # m in [kg kg-1], RH in [frac]
    all_vars = aer_ext_scat(m_range, RH_i_array)
    beta[str(RH_i)] = all_vars['beta']


# Beta and alpha with different N_0
# -------------------------------------
# Fix RH and extract out specific m values, in order to compare across the different N_0s.

# Different N_0 [cm-3] (variable name style: N0 per cm3)
N0_pcm3 = [3769, 4461, 5426, 6824, 4471] # as in table 4 of paper 1
N0_pm3 = [i*1.0e06 for i in N0_pcm3]

alpha_diff_N0 = []
beta_diff_N0 = []


# for N0_pm3_i, N0_pcm3_i in zip(N0_pm3, N0_pcm3):
for N0_pm3_i in N0_pm3:

    for RH_i in [60.0]:

        # create a numpy array of the length of m_range full of RH_i to use as input to aer_ext_scat
        RH_i_array = np.empty(len(m_range))
        RH_i_array[:] = RH_i

        # m in [kg kg-1], RH in [frac]
        all_vars = aer_ext_scat(m_range, RH_i_array, N0=N0_pm3_i)

        # extract out beta at a set m
        # m_range[18] = 18.0 micrograms kg-1
        alpha_diff_N0 += [all_vars['alpha'][18]]
        beta_diff_N0 += [all_vars['beta'][18]]

# -----------------------------------------------------------
# 3. Plot each with respect to beta
# -----------------------------------------------------------

# plot beta for the range of m at different RH
fig = plt.figure(figsize=(6, 3.5))

for RH_i in RH_range:

    RH_i_str = str(RH_i)
    beta_i = beta[RH_i_str]

    if RH_i > 80.0:
        ls = '--'
    else:
        ls = '-'

    plt.plot(m_range, beta_i, label='RH = ' + RH_i_str, linestyle=ls)


plt.legend(loc='best', fancybox=True, framealpha=0.5)

# plt.xlabel('aerosol [micrograms kg^-1]')
plt.xlabel(r'$m \/\mathrm{[\mu g\/ kg^{-1}]}$', labelpad=2)
plt.ylabel(r'$\beta_{m,\/unatt} \/\mathrm{[m^{-1} sr^{-1}]}$', labelpad=2)
plt.xlim([0.0, np.max(m_range)])
plt.ylim([0.0, 6.0e-6])
ax = plt.gca()
# ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax.set_yticklabels(["{:.1e}".format(t) for t in ax.get_yticks()])

# u'${10^{-6}}$'


plt.grid()
plt.tight_layout()

fn2 = 'm_and_RH_vs_beta_v0.2.png'
#plt.savefig(savedir + fn2)

print 'END PROGRAM'