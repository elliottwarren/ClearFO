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
        filename = 'sp_ew_910_ext_f(RH)_910-910nm.csv'

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

def aer_ext_scat(q_aer, RH, r0 = FOcon.r0_haywood, p = FOcon.p_aer,
                 B=FOcon.B_activation_haywood, S = FOcon.LidarRatio['Aerosol'],
                 N0=FOcon.N0_aer, m0 = FOcon.m0_aer, eta = FOcon.eta, Q = FOcon.Q_ext_aer,
                 r_md = []):
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

    # eq 12 - calc rm for RH greater than critical
    if RH.shape[0] == 1:
        if RH >= RH_crit:
            rm = np.ones(RH.shape) - (B / np.log(RH_factor * RH))
            rm2 = np.power(rm, 1. / 3.)
            rm = np.array(r_md) * rm2
        else:
            rm = r_md

    else: # RH shape > 1
        # set up rm array
        rm = np.empty(RH.shape)
        rm[:] = np.nan

        # find instances above and below critical threshold
        gt_idx = np.where(RH >= RH_crit)[0]
        ls_idx = np.where(RH < RH_crit)[0]

        # swell those above the threshold, keep rm as r_md where below.
        if gt_idx.size: # if not empty
            rm[gt_idx] = np.ones(gt_idx.shape) - (B / np.log(RH_factor * RH[gt_idx]))
            rm2 = np.power(rm[gt_idx], 1. / 3.)
            rm[gt_idx] = np.array(r_md) * rm2
            # ToDo only works if RH shape > 1 BUT q_aer shape is also == 1!
        rm[ls_idx] = r_md

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
    alpha_a = (eta * Q) * np.pi * N_aer * np.power(rm, 2)

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

def calc_aer_ext_vars(var_order, variables, q_aer_fix, RH_fix):

    """
    calculate the extinction coefficient for each of the input

    :param variables:
    :return: alpha, beta for all vars in var order
    :return: alpha_c, beta_c for the control run with a fixed m and RH
    :return: var order, so the variables will be plotted in order afterward.
    """

    # define dictionaries to store alpha and beta values in
    alpha = {}
    beta = {}
    var_range = {}

    for key, value in variables.iteritems():

        # create ranges that have 100 elements in them
        var_range[key] = create_range(value)

    for key in var_order:

        if key == 'RH':
            out = aer_ext_scat(q_aer_fix, var_range[key])

        elif key == 'm':
            out = aer_ext_scat(var_range[key], RH_fix)
        else:
            # ** {key: var_range} = keyword argument unpacking
            out = aer_ext_scat(q_aer_fix, RH_fix, **{key: var_range[key]})

        # extract out the alpha and beta
        alpha[key] = out['alpha']
        beta[key] = out['beta']

    # control
    out_c = aer_ext_scat(q_aer_fix, RH_fix)
    alpha_c = out_c['alpha']
    beta_c = out_c['beta']

    return alpha, beta, alpha_c, beta_c, var_range

# plotting

def plot_beta(var_order, beta, beta_c, variables, var_range, clearFo_dict, savestr):

    """
    plot the beta plot

    :param var_order:
    :param beta:
    :param beta_c:
    :param variables:
    :param var_range:
    :param clearFo_dict:
    :return:
    """

    # set up plot
    fig = plt.figure(figsize=(10, 6))
    ax = plt.subplot2grid((27, 1), (0, 0), rowspan=26 - len(var_order))

    for param in var_order:

        value = beta[param]

        if (param == 'RH') | (param == 'm'):
            line, = ax.plot(np.arange(0, 101), np.log10(value[0:101]), label=variables[param][2],
                            color=variables[param][3], linewidth=2)
            line.set_dashes([10, 5, 100, 5])
        else:
            # 0:101 because computer rounding makes r0 have 102 enteries
            ax.plot(np.arange(0, 101), np.log10(value[0:101]), label=variables[param][2], color=variables[param][3])

    ax.plot([0, 100], [np.log10(beta_c), np.log10(beta_c)], label=r'$control$', color='black', linewidth=2, ls='--')

    # adjust axis so the legends can be placed on the side
    plt.tight_layout()
    fig.subplots_adjust(top=0.95, right=0.80, left=0.1, bottom=0.1)
    ax.legend(fontsize=14, bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.0)

    ax.set_ylabel(r'$log_{10}(\beta) \/\/\mathrm{[m^{-1} \/sr^{-1}]}$', fontsize=16)
    ax.yaxis.labelpad = 0
    ax.set_ylim([-6.6, -4.4])
    ax.axes.get_xaxis().set_ticks([])

    # add all the extra axes which are shown in var_order
    plot_and_add_x_axes(clearFo_dict, var_order, var_range, variables)

    plt.savefig(savedir + savestr + '_v_beta.png')  # filename

    plt.close(fig)

    return fig

def plot_and_add_x_axes(clearFo_dict, var_order, var_range, variables):

    """
    Add all the extra axes for beta or alpha plot
    :return: ax
    """

    # scaling factors
    # ['Q', 'p', 'S', 'B', 'r0', 'm0', 'N0', 'm', 'RH']

    scaling = {'N0': 1.0/1.0e6, # cm-3
               'r0': 1.0e6, # micro meters
               'm0': 1.0e9} # microgram kg-1

    # add the extra x axis
    for i, key in zip(np.arange(26 - len(var_order), 26), var_order):
        ax_i = plt.subplot2grid((27, 1), (i, 0))

        if key in scaling:
            ax_i.plot(var_range[key] * scaling[key], var_range[key] * scaling[key], label=variables[key][2], color=variables[key][3])
            ax_i.scatter(clearFo_dict[key] * scaling[key], 0.0, color=variables[key][3], edgecolors='black')
            ax_i.ticklabel_format(useOffset=False, style='plain')
            ax_i.set_xlim([var_range[key][0] * scaling[key], var_range[key][-1] * scaling[key]])
        else:
            ax_i.plot(var_range[key], var_range[key], label=variables[key][2], color=variables[key][3])
            ax_i.scatter(clearFo_dict[key], 0, color=variables[key][3], edgecolors='black')
            ax_i.set_xlim([var_range[key][0], var_range[key][-1]])

        # clean up axis
        ax_i.spines['top'].set_visible(False)
        ax_i.spines['right'].set_visible(False)
        ax_i.spines['left'].set_visible(False)
        ax_i.yaxis.set_visible(False)
        # ax_i.set_xlim([np.min(var_range[key]), np.max(var_range[key])])

        ax_i.axes.get_yaxis().set_ticks([])
        # ax_i.get_xaxis().get_major_formatter().set_scientific(False)

        # colours and remove the line
        ax_i.tick_params(axis='x', colors=variables[key][3])
        ax_i.xaxis.label.set_color(variables[key][3])
        ax_i.spines['bottom'].set_color(variables[key][3])
        ax_i.lines.pop(0)

    return


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

# fixed murk at 40 mg kg-1
q_aer_fix = np.array([10])

# create range of murk and RH
m_range = np.arange(0.0, 60.05, 0.05)
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


# -----------------------------------------------------------
# 3. Plot each with respect to beta
# -----------------------------------------------------------

# plot beta for the range of m at different RH
fig = plt.figure(figsize=(6, 3.5))

for RH_i in RH_range:

    RH_i_str = str(RH_i)
    beta_i = beta[RH_i_str]

    plt.plot(m_range, beta_i, label=RH_i_str)


plt.legend(loc='best', fancybox=True, framealpha=0.5)

# plt.xlabel('aerosol [micrograms kg^-1]')
plt.xlabel(r'$aerosol \/\mathrm{[\mu g\/ kg^{-1}]}$', labelpad=2)
plt.ylabel('beta')
#plt.xlim([0.0, 100.0])
#plt.ylim([0.0, 0.15])
plt.grid()
plt.tight_layout()

fn2 = 'm_and_RH_vs_beta.png'
plt.savefig(savedir + fn2)



print 'END PROGRAM'