"""
Script to make sensitivity plots of the different FO parameters

1. no height stuff. Just fixing all parameters/variables and varying one of them at a time
2. Constant RH and murk aerosol profiles in height to understand impact on transmission

"""

import matplotlib.pyplot as plt
import numpy as np
import datetime as dt

from forward_operator import FOconstants as FOcon
from forward_operator import FOUtils as FO

def aer_ext_scat(q_aer, RH, r0 = FOcon.r0_haywood, p = FOcon.p_aer,
                 B=FOcon.B_activation_haywood, S = FOcon.LidarRatio['Aerosol'],
                 N0=FOcon.N0_aer, m0 = FOcon.m0_aer, eta = FOcon.eta, Q = FOcon.Q_ext_aer):

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

    # Calculate extinction coefficient
    # eqns. 17-18 in Clark et.al. (2008)
    alpha_a = (eta * Q) * np.pi * N_aer * np.power(rm, 2)

    # Calculate backscatter using a constant lidar ratio
    beta_a = alpha_a / S

    return alpha_a, beta_a

def calc_aer_ext_vars(variables, q_aer_fix, RH_fix):

    """
    calculate the extinction coefficient for each of the input

    :param variables:
    :return: alpha, beta for all vars in var order
    :return: alpha_c, beta_c for the control run with a fixed m and RH
    :return: var order, so the variables will be plotted in order afterward.
    """

    var_order = ['Q', 'p', 'S', 'B', 'r0', 'm0', 'N0', 'm', 'RH']

    # define dictionaries to store alpha and beta values in
    alpha = {}
    beta = {}
    var_range = {}

    for key, value in variables.iteritems():

        step = (value[1] - value[0]) / 100.0
        var_range[key] = np.arange(value[0], value[1] + step, step)

        if key == 'RH':
            alpha[key], beta[key] = aer_ext_scat(q_aer_fix, var_range[key])
        elif key == 'm':
            alpha[key], beta[key] = aer_ext_scat(var_range[key], RH_fix)
        else:
            # ** {key: var_range} = keyword argument unpacking
            alpha[key], beta[key] = aer_ext_scat(q_aer_fix, RH_fix, **{key: var_range[key]})

    # control
    alpha_c, beta_c = aer_ext_scat(q_aer_fix, RH_fix)

    return alpha, beta, alpha_c, beta_c, var_range, var_order

# -----------------------------------------------------------
# Constants
# -----------------------------------------------------------

savedir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/figures/sensitivity/'

# Variables for model data
day = dt.datetime(2016, 05, 04)
model_type = 'UKV'
res = FOcon.model_resolution[model_type]
modDatadir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/data/' + model_type + '/'

# get model data for the height
mod_data = FO.mod_site_extract_calc(day, {'IMU': [-0.106086, 51.52615]}, modDatadir, model_type, res)
z_mod = mod_data['IMU']['level_height']

# fixed murk at 40 mg kg-1
q_aer_fix = np.array([40])

# fix RH at 60 %
RH_fix = np.array([60])

# -----------------------------------------------------------
# 1. calculate backscatters
# -----------------------------------------------------------

# [min, max]
# ToDo N0, m0 needs changing
variables = {'r0': [0.1e-6, 0.18e-6, r'$r_0 \/\mathrm{[\mu m]}$ ', ('#D822D8')], #
             'p': [0.0, 1.0/3, r'$p$','black'], #
             'B': [0.1, 1.0, r'$B$', (0.2,0.6,0.6)], #
             'S': [20.0, 70.0, r'$S \/\mathrm{[sr^{-1}]}$', ('#FFAB00')], #
             'N0': [5.0e8, 2.0e10, r'$N_0 \/\mathrm{[cm^{-3}]}$', (0.5,0.3,0.3)], #
             'm0': [1.0e-8, 3.0e-8, r'$m_0 \/\mathrm{[\mu g \/\/kg^{-1}]}$','c'], #
             'Q': [0.1, 4.0, r'$Q$','y'], #
             'RH': [0.0, 100.0, r'$RH \/\mathrm{[\%]}$','b'], #
             'm': [5.0, 150.0, r'$m \/\mathrm{[\mu g \/\/kg^{-1}]}$','r']} #

alpha, beta, alpha_c, beta_c, var_range, var_order = calc_aer_ext_vars(variables, q_aer_fix, RH_fix)

# log 10 all the beta values
for key, value in beta.iteritems():
    beta[key] = np.log10(value)


# remove alpha S as it only has one value
del alpha['S']

# # -----------------------------------------------------------
# # 2. Plot each with respect to beta
# # -----------------------------------------------------------

# -----------------------------------------------------------
# 2. Plot each with respect to beta
# -----------------------------------------------------------

for plot_variable in [beta, alpha]:

    # set up plot
    fig = plt.figure(figsize=(10, 6))
    ax = plt.subplot2grid((27, 1), (0, 0), rowspan=17)

    for param in var_order:

        value = plot_variable[param]

        if (param == 'RH') | (param == 'm'):
            line, = ax.plot(np.arange(0, 101), np.log10(value[0:101]), label=variables[param][2], color=variables[param][3], linewidth=2)
            line.set_dashes([10, 5, 100, 5])
        else:
            # 0:101 because computer rounding makes r0 have 102 enteries
            ax.plot(np.arange(0,101), np.log10(value[0:101]), label=variables[param][2], color=variables[param][3])

    ax.plot([0,100], [np.log10(beta_c), np.log10(beta_c)], label=r'$control$', color='black', linewidth=2, ls='--')

    # adjust axis so the ledgends can be placed on the side
    plt.tight_layout()
    fig.subplots_adjust(top=0.95, right=0.80, left=0.1, bottom=0.1)
    ax.legend(fontsize=14, bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.0)

    ax.set_ylabel(r'$log_{10}(\beta) \/\/\mathrm{[m^{-1} \/sr^{-1}]}$', fontsize=16)
    ax.yaxis.labelpad = 0
    ax.set_ylim([-6.6, -4.4])
    ax.axes.get_xaxis().set_ticks([])

    # add the extra x axis
    for i, key in zip(np.arange(17, 26), var_order):
        ax_i = plt.subplot2grid((27, 1), (i, 0))

        if key == 'N0': # cm-3
            ax_i.plot(var_range[key]/1.0e6, var_range[key]/1.0e6, label=variables[key][2], color=variables[key][3])
            ax_i.ticklabel_format(useOffset=False, style='plain')
            ax_i.set_xlim([var_range[key][0]/1.0e6, var_range[key][-1]/1.0e6])
        elif key == 'r0': # micro meters
            ax_i.plot(var_range[key] *1.0e6, var_range[key]*1.0e6, label=variables[key][2], color=variables[key][3])
            ax_i.ticklabel_format(useOffset=False, style='plain')
            ax_i.set_xlim([var_range[key][0]*1.0e6, var_range[key][-1]*1.0e6])
        elif key == 'm0':  # micro meters
            ax_i.plot(var_range[key] * 1.0e9, var_range[key] * 1.0e9, label=variables[key][2], color=variables[key][3])
            ax_i.ticklabel_format(useOffset=False, style='plain')
            ax_i.set_xlim([var_range[key][0] * 1.0e9, var_range[key][-1] * 1.0e9])
        else:
            ax_i.plot(var_range[key], var_range[key], label=variables[key][2], color=variables[key][3])
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

    plt.savefig(savedir + 'all_params_v_beta.png')  # filename




# plt.close(fig)


# # set up plot
# fig = plt.figure(figsize=(10, 6))
# ax = plt.subplot2grid((27, 1), (0, 0), rowspan=17)
#
# for param in var_order:
#
#     value = beta[param]
#
#     if (param == 'RH') | (param == 'm'):
#         line, = ax.plot(np.arange(0, 101), np.log10(value[0:101]), label=variables[param][2], color=variables[param][3], linewidth=2)
#         line.set_dashes([10, 5, 100, 5])
#     else:
#         # 0:101 because computer rounding makes r0 have 102 enteries
#         ax.plot(np.arange(0,101), np.log10(value[0:101]), label=variables[param][2], color=variables[param][3])
#
# ax.plot([0,100], [np.log10(beta_c), np.log10(beta_c)], label=r'$control$', color='black', linewidth=2, ls='--')
#
# # adjust axis so the ledgends can be placed on the side
# plt.tight_layout()
# fig.subplots_adjust(top=0.95, right=0.80, left=0.1, bottom=0.1)
# ax.legend(fontsize=14, bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.0)
#
# ax.set_ylabel(r'$log_{10}(\beta) \/\/\mathrm{[m^{-1} \/sr^{-1}]}$', fontsize=16)
# ax.yaxis.labelpad = 0
# ax.set_ylim([-6.6, -4.4])
# ax.axes.get_xaxis().set_ticks([])
#
# # add the extra x axis
# for i, key in zip(np.arange(17, 26), var_order):
#     ax_i = plt.subplot2grid((27, 1), (i, 0))
#
#     if key == 'N0': # cm-3
#         ax_i.plot(var_range[key]/1.0e6, var_range[key]/1.0e6, label=variables[key][2], color=variables[key][3])
#         ax_i.ticklabel_format(useOffset=False, style='plain')
#         ax_i.set_xlim([var_range[key][0]/1.0e6, var_range[key][-1]/1.0e6])
#     elif key == 'r0': # micro meters
#         ax_i.plot(var_range[key] *1.0e6, var_range[key]*1.0e6, label=variables[key][2], color=variables[key][3])
#         ax_i.ticklabel_format(useOffset=False, style='plain')
#         ax_i.set_xlim([var_range[key][0]*1.0e6, var_range[key][-1]*1.0e6])
#     elif key == 'm0':  # micro meters
#         ax_i.plot(var_range[key] * 1.0e9, var_range[key] * 1.0e9, label=variables[key][2], color=variables[key][3])
#         ax_i.ticklabel_format(useOffset=False, style='plain')
#         ax_i.set_xlim([var_range[key][0] * 1.0e9, var_range[key][-1] * 1.0e9])
#     else:
#         ax_i.plot(var_range[key], var_range[key], label=variables[key][2], color=variables[key][3])
#         ax_i.set_xlim([var_range[key][0], var_range[key][-1]])
#
#     # clean up axis
#     ax_i.spines['top'].set_visible(False)
#     ax_i.spines['right'].set_visible(False)
#     ax_i.spines['left'].set_visible(False)
#     ax_i.yaxis.set_visible(False)
#     # ax_i.set_xlim([np.min(var_range[key]), np.max(var_range[key])])
#
#     ax_i.axes.get_yaxis().set_ticks([])
#     # ax_i.get_xaxis().get_major_formatter().set_scientific(False)
#
#     # colours and remove the line
#     ax_i.tick_params(axis='x', colors=variables[key][3])
#     ax_i.xaxis.label.set_color(variables[key][3])
#     ax_i.spines['bottom'].set_color(variables[key][3])
#     ax_i.lines.pop(0)
#
# plt.savefig(savedir + 'all_params_v_beta.png')  # filename
# # plt.close(fig)


# -----------------------------------------------------------
# 3. Plot each with respect to extinction coefficient
# -----------------------------------------------------------

# set up plot
fig = plt.figure(figsize=(8, 5))
ax = plt.subplot2grid((1, 1), (0, 0))

for param, value in alpha.iteritems():

    if (param == 'RH') | (param == 'm'):
        line, = ax.plot(np.arange(0, 101), 1e3*value[0:101], label=variables[param][2], color=variables[param][3], linewidth=2)
        line.set_dashes([10, 5, 100, 5])
    else:
        # 0:101 because computer rounding makes r0 have 102 enteries
        ax.plot(np.arange(0,101), 1e3*value[0:101], label=variables[param][2], color=variables[param][3])


ax.plot([0,100], [1e3*alpha_c, 1e3*alpha_c], label=r'$control$', color='black', linewidth=2, ls='--')

# adjust axis so the ledgends can be placed on the side
plt.tight_layout()
fig.subplots_adjust(top=0.95, right=0.80, left=0.1, bottom=0.1)
ax.legend(fontsize=14, bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.0)

ax.set_ylabel(r'$\sigma_{ext} \/\/\mathrm{[km^{-1}]}$', fontsize=16)
ax.yaxis.labelpad = 0
ax.set_xlabel('step increase (%)')
ax.set_ylim([-0.0, 3.0])

plt.savefig(savedir + 'all_params_v_alpha.png')  # filename
# plt.close(fig)


print 'END PROGRAM'