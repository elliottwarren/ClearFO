"""
Script to make sensitivity plots of the different FO parameters

1. no height stuff. Just fixing all parameters/variables and varying one of them at a time
2. Constant RH and murk aerosol profiles in height to understand impact on transmission

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

# fixed murk at 40 mg kg-1
q_aer_fix = np.array([10])

# fix RH at 60 %
RH_fix = np.array([60])
RH_var = np.arange(0, 101)

# variables to plot
# all potential variables -> = ['Q', 'p', 'S', 'B', 'r0', 'm0', 'N0', 'm', 'RH']
# var_order = ['m', 'RH']
var_order = ['Q', 'p', 'S', 'B', 'r0', 'm0', 'N0', 'm', 'RH']

if len(var_order) == 9:
    savestr = 'all_params'
else:
    savestr = 'some_params'

# all FO vars in a single dict
clearFo_dict =  {'r0': FOcon.r0_haywood,
                 'p': FOcon.p_aer,
                 'B': FOcon.B_activation_haywood,
                 'S': FOcon.LidarRatio['Aerosol'],
                 'N0': FOcon.N0_aer,
                 'm0': FOcon.m0_aer,
                 'eta': FOcon.eta,
                 'Q': FOcon.Q_ext_aer,
                 'm': q_aer_fix,
                 'RH': RH_fix,
                 'ln(RH)': np.log(RH_fix),
                 'B/ln(RH)': FOcon.B_activation_haywood / np.log(RH_fix),
                 'r_md': 0.12e-6} # value when m = 40 microgram kg-1, RH = 60 %

# scaling the dependent variables
y_extras = {'r_md': [1.0e6, r'$r_{md} \/\/\mathrm{[\mu m]}$'], # micro meters
            'N': [1.0/1.0e6, r'$N \/\mathrm{[cm^{-3}]}$'], # cm-3
            'r_m': [1.0e6, r'$r_{m} \/\/\mathrm{[\mu m]}$']}  # micro meters

# -----------------------------------------------------------
# 1. calculate backscatters
# -----------------------------------------------------------

# [min, max]
# ToDo this should not be looped over, only called from within.
variables = {'r0': [0.1e-6, 0.2e-6, r'$r_0 \/\mathrm{[\mu m]}$ ', ('#D822D8')], #
             'p': [0.0, 1.0/3, r'$p$','black'], #
             'B': [0.1, 1.0, r'$B$', (0.2,0.6,0.6)], #
             'S': [20.0, 70.0, r'$S \/\mathrm{[sr^{-1}]}$', ('#FFAB00')], #
             'N0': [5.0e8, 2.0e10, r'$N_0 \/\mathrm{[cm^{-3}]}$', (0.5,0.3,0.3)], #
             'm0': [1.0e-8, 3.0e-8, r'$m_0 \/\mathrm{[\mu g \/\/kg^{-1}]}$','c'], #
             'Q': [0.1, 4.0, r'$Q$','y'], #
             'RH': [0.0, 100.0, r'$RH \/\mathrm{[\%]}$','b'], #
             'ln(RH)': [np.nan, np.nan, r'$ln(RH) \/\mathrm{[\%]}$',(0.3, 1.0, 1.0)],
             'B/ln(RH)': [np.nan, np.nan, r'$B/ln(RH) \/\mathrm{[\%]}$',(1.0, 0.5, 0.9)],
             'm': [10.0, 60.0, r'$m \/\mathrm{[\mu g \/\/kg^{-1}]}$','r'],
             'r_md': [0.8e-7, 0.26e-6, r'$r_{md} \/\mathrm{[\mu m]}$',(0.7, 0.7, 0.1)]} #

# -----------------------------------------------------------
# 2. Calculating
# -----------------------------------------------------------

# Beta and Alpha
# ----------------
# do first to get the var_range out for later calculations (r_md)
alpha, beta, alpha_c, beta_c, var_range = calc_aer_ext_vars(var_order, variables, q_aer_fix, RH_fix)

# remove alpha S as it only has one value
if len(var_order) == 9:
    del alpha['S']

# store all the variables that will have sensitivity tests wrt.
y_values = {}

# r_md
# ----------------
# just be careful with the keys here... each line of the dictionary has a parameter in it 3 times.
y_values['r_md'] = {'m': aer_ext_scat(var_range['m'], RH_fix)['r_md'],
                    'r0': aer_ext_scat(q_aer_fix, RH_fix, r0=var_range['r0'])['r_md'],
                    'm0': aer_ext_scat(q_aer_fix, RH_fix, m0=var_range['m0'])['r_md'],
                    'p': aer_ext_scat(q_aer_fix, RH_fix, p=var_range['p'])['r_md']}

# N
# ----------------
# just be careful with the keys here... each line of the dictionary has a parameter in it 3 times.
y_values['N'] = {'m': aer_ext_scat(var_range['m'], RH_fix)['N'],
                    'N0': aer_ext_scat(q_aer_fix, RH_fix, N0=var_range['N0'])['N'],
                    'm0': aer_ext_scat(q_aer_fix, RH_fix, m0=var_range['m0'])['N'],
                    'p': aer_ext_scat(q_aer_fix, RH_fix, p=var_range['p'])['N']}

#----

# r_m
# ----------------

r_m_values       = {'r_md': aer_ext_scat(q_aer_fix, RH_fix, r_md=var_range['r_md'])['r_m'],
                   'B': aer_ext_scat(q_aer_fix, RH_fix, B=var_range['B'])['r_m'],
                   'ln(RH)': aer_ext_scat(q_aer_fix, var_range['RH'])['r_m'],
                   'B/ln(RH)': aer_ext_scat(q_aer_fix, var_range['RH'])['r_m']}

# r0 vs ext
# ----------------
test = {'haywood': aer_ext_scat(var_range['m'], RH_fix, r0=FOcon.r0_haywood)['alpha'],
          'clark': aer_ext_scat(var_range['m'], RH_fix, r0=FOcon.r0_clark)['alpha']}

# r0 vs ext
# ----------------
messB = {'haywood': aer_ext_scat(var_range['m'], RH_fix, B=FOcon.B_activation_haywood)['alpha'],
          'clark': aer_ext_scat(var_range['m'], RH_fix, B=FOcon.B_activation_clark)['alpha']}

# add some extra var ranges for the r_m values
var_range['ln(RH)'] = np.log(var_range['RH'])
var_range['B/ln(RH)'] = FOcon.B_activation_haywood / np.log(var_range['RH'])

# -----------------------------------------------------------
# 3. Plot each with respect to beta
# -----------------------------------------------------------

# Beta
# -----------
fig = plot_beta(var_order, beta, beta_c, variables, var_range, clearFo_dict, savestr)

# Alpha
# -----------
# set up plot
fig = plt.figure(figsize=(10, 6))
ax = plt.subplot2grid((27, 1), (0, 0), rowspan=26 - len(var_order))

# copy array and remove S as it isn't in the alpha calculation
var_order_alpha = deepcopy(var_order)
if len(var_order) == 9:
    var_order_alpha.remove('S')

for param in var_order_alpha:

    value = alpha[param]


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
# ax.set_xlabel('step increase (%)')
ax.set_ylim([-0.0, 0.4])
ax.axes.get_xaxis().set_ticks([])

plot_and_add_x_axes(clearFo_dict, var_order_alpha, var_range, variables)

plt.savefig(savedir + savestr + '_v_alpha.png')  # filename
plt.close(fig)

# dependent variables
# ---------------------

# inc. r_md, N

# loop through each main dependent variable
for y_key, y_dict in y_values.iteritems():

    # set up figure
    fig = plt.figure(figsize=(10, 6))
    ax = plt.subplot2grid((27, 1), (0, 0), rowspan=26 - len(y_dict.keys()))

    # loop through all the independent variables for this dependent variable
    for param, value in y_dict.iteritems():

        if (param == 'RH') | (param == 'm'):
            line, = ax.plot(np.arange(0, 101), y_extras[y_key][0]*value[0:101], label=variables[param][2], color=variables[param][3], linewidth=2)
            line.set_dashes([10, 5, 100, 5])
        else:
            # 0:101 because computer rounding makes r0 have 102 enteries
            ax.plot(np.arange(0, 101), y_extras[y_key][0]*value[0:101], label=variables[param][2], color=variables[param][3])

    # control turned off for now...
    # ax.plot([0,100], [1e3*alpha_c, 1e3*alpha_c], label=r'$control$', color='black', linewidth=2, ls='--')

    # adjust axis so the ledgends can be placed on the side
    plt.tight_layout()
    fig.subplots_adjust(top=0.95, right=0.80, left=0.1, bottom=0.1)
    ax.legend(fontsize=14, bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.0)

    ax.set_ylabel(y_extras[y_key][1], fontsize=16)
    ax.yaxis.labelpad = 0
    # ax.set_xlabel('step increase (%)')
    # ax.set_ylim([-0.0, 1.5])
    ax.axes.get_xaxis().set_ticks([])

    if y_key == 'r_md':
        ax.set_ylim([0.08, 0.24])

    plot_and_add_x_axes(clearFo_dict, y_dict.keys(), var_range, variables)

    plt.savefig(savedir + savestr + '_v_' + y_key + '.png')  # filename
    plt.close(fig)

# r_m twin x!
# ---------------------
#
# # set up figure
# fig, ax = plt.subplots()
#
# axes = [ax, ax.twinx(), ax.twinx(), ax.twinx()]
#
# # move last few spines lower
# axes[-2].spines['bottom'].set_position(('axes', 1.1))
# axes[-1].spines['bottom'].set_position(('axes', 1.2))
#
# axes[-1].set_frame_on(True)
# axes[-1].patch.set_visible(False)
#
# axes[-2].set_frame_on(True)
# axes[-2].patch.set_visible(False)
#
# # loop through all the independent variables for this dependent variable
# for i, param in zip(np.arange(4), r_m_values.iterkeys()):
#
#     value = r_m_values[param]
#
#     # 0:101 because computer rounding makes r0 have 102 enteries
#     axes[i].plot(var_range[param], 1.0e6*value[0:101], label=variables[param][2], color=variables[param][3])
#
# # adjust axis so the ledgends can be placed on the side
# # plt.tight_layout()
# # fig.subplots_adjust(top=0.95, right=0.80, left=0.1, bottom=0.1)
# axes[0].legend(fontsize=14, bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.0)
#
# axes[0].set_ylabel(r'$r_{m} \/\/\mathrm{[\mu m]}$', fontsize=16)
# axes[0].yaxis.labelpad = 0
# # ax.set_xlabel('step increase (%)')
# # ax.set_ylim([-0.0, 1.5])
# axes[0].axes.get_xaxis().set_ticks([])
#
# # plot_and_add_x_axes(clearFo_dict, y_dict.keys(), var_range, variables)
#
# plt.savefig(savedir + savestr + '_v_r_m' + '.png')  # filename
# plt.close(fig)



for i, param in zip(np.arange(4),r_m_values.iterkeys()):

    value = r_m_values[param]

    fig = plt.figure(figsize=(10, 6))
    plt.plot(var_range[param], 1.0e6*value[0:101], label=variables[param][2], color=variables[param][3])
    fig.subplots_adjust(top=0.95, right=0.80, left=0.1, bottom=0.1)
    plt.legend(fontsize=14, bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.0)
    plt.savefig(savedir + str(i) + '_v_r_m' + '.png')  # filename


# r0 vs ext coeff
fig = plt.figure(figsize=(6, 3.5))
plt.plot(var_range['m'], test['clark']*1e3, label='r0 = 0.16')
plt.plot(var_range['m'], test['haywood']*1e3, label='r0 = 0.11')
plt.xlabel(r'$m \/\mathrm{[\mu g \/\/kg^{-1}]}$')
plt.ylabel(r'$\sigma_{ext} \/\/\mathrm{[km^{-1}]}$')
plt.grid()
plt.legend(loc='best', fancybox=True, framealpha=0.5)
plt.tight_layout()
plt.savefig(savedir + 'r0_v_extCoeff' + '.png')


# B vs ext coeff
fig = plt.figure(figsize=(6, 3.5))
plt.plot(var_range['m'], messB['clark']*1e3, label='B = 0.5')
plt.plot(var_range['m'], messB['haywood']*1e3, label='B = 0.14')
plt.xlabel(r'$m \/\mathrm{[\mu g \/\/kg^{-1}]}$')
plt.ylabel(r'$\sigma_{ext} \/\/\mathrm{[km^{-1}]}$')
plt.grid()
plt.legend(loc='best', fancybox=True, framealpha=0.5)
plt.tight_layout()
plt.savefig(savedir + 'B_v_extCoeff' + '.png')

print 'END PROGRAM'