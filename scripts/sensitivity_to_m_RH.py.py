"""
Plot different FO variables with respect to m and/or RH.
Display sensitivity of each to FO input variables.

Created by Cristina Charlton-Perez - 07/12/2016
Edited by Elliott Warren - 28/04/2017
"""




import numpy as np
import matplotlib.pyplot as plt  # plotting library (not Met Office)
from decimal import Decimal
from os.path import join

from forward_operator import FOUtils as FO
from forward_operator import FOconstants as FOcon


def compute_N(N0, m, m0, p=1.0 / 6.0):

    N = N0 * np.power(m / m0, (1.0 - 3.0 * p))

    return N


def compute_dry(m, m0, r0, p=1.0 / 6.0):

    dum = (m / m0)
    rdry = r0 * np.power(dum, p)

    return rdry


def compute_wet(r_md, RH, B=0.5, p=1.0 / 3.0, RH_crit=0.38):

    # Compute a wet radius where RH exceeds critical value

    RH_ge_RHcrit = np.ma.masked_less(RH, RH_crit)
    rm1 = np.ma.ones(RH.shape) - (B / np.ma.log(RH_ge_RHcrit))
    rm2 = np.ma.power(rm1, p)
    rm = np.ma.array(r_md) * rm2

    rm = np.ma.MaskedArray.filled(rm, [0.0])
    #where_lt_crit = np.where(np.logical_or(RH < RH_crit, rm == 0.0))
#     print 'where_lt_crit = ', where_lt_crit
#     print type(where_lt_crit)
# #     logical_and(my_array > 3, my_array < 7)
#     if np.where(rm == 0.0):
#         print "RM EQUALS ZERO!!!!!"
#         print rm

    rwet = np.where(RH < RH_crit, r_md, rm)

    return rwet


def plot_r_md(m, m0_aer, r0, p_aer, N0_aer, kg_2_micro_g, m_2_microns, savedir):

    """ plot m vs r_md (dry mean radius) """

    fig = plt.figure(figsize=(6, 3.5))

    r_md = compute_dry(m, m0_aer, r0, p=p_aer)
    plt.plot(m * kg_2_micro_g, m_2_microns * r_md, 'k', label='p = 1/6')

    # add line to show where m0 is.
    plt.plot([m0_aer * kg_2_micro_g, m0_aer * kg_2_micro_g], [np.min(m * kg_2_micro_g), np.max(m * kg_2_micro_g)],
             color='red', ls='--')

    ax2 = plt.gca()
    # plt.legend()
    # ax2.legend(loc='best', fancybox=True, framealpha=0.5)

    # plt.xlabel('aerosol [micrograms kg^-1]')
    ax2.set_xlabel(r'$aerosol \/\mathrm{[\mu g\/ kg^{-1}]}$', labelpad=2)
    plt.ylabel(r'$r_{md} \/\/\mathrm{[\mu m]}$')
    plt.xlim([0.0, 100.0])
    plt.ylim([0.0, 0.15])
    plt.grid()
    plt.tight_layout()

    fn2 = 'Vary_p_wrt_r0_' + str(N0_aer) + '_m0_' + str(m0_aer) + '_r0_' + str(r0) + '.png'
    plt.savefig(join(savedir, fn2))

    return fig


def plot_N(m, m0_aer, r0, p_aer, N0_aer, kg_2_micro_g, m3_2_cm3, savedir):

    """ plot N as a function of m"""

    fig = plt.figure(figsize=(6, 3.5))

    N_aerosol = compute_N(N0_aer, m, m0_aer, p=p_aer)
    plt.plot(m * kg_2_micro_g, N_aerosol * m3_2_cm3, 'g', linewidth=2.0, label='p = 1/6')

    plt.plot([m0_aer * kg_2_micro_g, m0_aer * kg_2_micro_g], [0.0, 20000.0],
             color='red', ls='--')

    ax2 = plt.gca()

    ax2.set_xlabel(r'$aerosol \/\mathrm{[\mu g\/ kg^{-1}]}$', labelpad=2)
    plt.ylabel(r'$N \/\/\mathrm{[cm^{-3}]}$')
    plt.grid()
    plt.xlim([np.min(m) * kg_2_micro_g, np.max(m) * kg_2_micro_g])
    plt.ylim([0.0, 20000.0])
    plt.tight_layout()

    fn = 'Vary_p_' + str(N0_aer) + '_m0_' + str(m0_aer) + '_r0_' + str(r0) + '.png'
    plt.savefig(join(savedir, fn))



    return


def main():

    # Path to save figures
    # savepath = '//data/dynamic/frpz/Aerosol_model_sensitivity'
    savedir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/figures/extraSensitivity/'

    # Conversion factors
    micro_g_2_kg = 1.0e-9
    kg_2_micro_g = 1.0e9
    m_2_microns = 1.0e6
    m3_2_cm3 = 1.0e-06

    # read in the FO constants

    RH_crit = 38.0 # critical value of relative humidity to start swelling the aerosol
    rho_a = FOcon.rho_a  # density of aerosol (for ammonium sulphate) [kg m-3]
    r0 = FOcon.r0_haywood   # radius of a 'standard' aerosol particle [m]
    p_aer = FOcon.p_aer
    B = FOcon.B_activation_haywood
    N0_aer = FOcon.N0_aer # (N0) standard number density of aerosol [m-3]
    m0_aer = FOcon.m0_aer # (m0) standard mass mixing ratio of the aerosol [kg kg-1]
    eta = FOcon.eta

    # FO set variable ranges
    m = np.arange(0, 101, 1) * micro_g_2_kg
    RH = np.arange(0, 101, 1) * 0.01


    # Plotting functions

    fig = plot_r_md(m, m0_aer, r0, p_aer, N0_aer, kg_2_micro_g, m_2_microns, savedir)

    fig = plot_N(m, m0_aer, r0, p_aer, N0_aer, kg_2_micro_g, m3_2_cm3, savedir)






####

    fig = plt.figure(figsize=(6, 3.5))
    r_md = compute_dry(m, m0_aer, r0, p=p_aer)
    print len(m), len(RH)
    # list_RH = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]
    # m_matrix = np.ones((len(m), len(list_RH)))
    #
    # colors = ["r-", "m-",  "g-", "c-", "b-", "k-",
    #           "r--", "m--", "g--", "c--", "b--", "k--"]

    list_RH = [0.4, 0.6, 0.8, 0.9, 0.95, 0.99]
    m_matrix = np.ones((len(m), len(list_RH)))

    colors = ["k-", "m-", "g-", "b-",
              "k--", "m--", "g--", "b--",]

    ncolor = 0

    plt.plot([m0_aer * kg_2_micro_g, m0_aer * kg_2_micro_g], [np.min(m * kg_2_micro_g), np.max(m * kg_2_micro_g)],
             color='red', ls='--')

    for iRH in list_RH:

        print "iRH = ", iRH
        r_wet = compute_wet(r_md, iRH * np.ones(r_md.shape), B=B, p=p_aer, RH_crit=0.38)
        m_matrix[:, iRH] = r_wet[:]

        print colors[ncolor]
        plt.plot(m * kg_2_micro_g, m_2_microns * r_wet,
                 colors[ncolor], label="RH=" + str(iRH))
        #plt.legend()
        ncolor += 1

    plt.ylabel(r'$r_{m} \/\/\mathrm{[\mu m]}$')
    plt.xlabel(r'$aerosol \/\mathrm{[\mu g\/ kg^{-1}]}$', labelpad=2)
    # Make transparent legend
    #ax = plt.gca()
    plt.legend(loc='best', fancybox=True, framealpha=0.8, fontsize=9)

    plt.xlim([0.0, 100.0])
    plt.ylim([0.0, 0.25])
    plt.grid()
    plt.tight_layout()

    fn = 'Vary_RH_N0_' + str(N0_aer) + '_m0_' + str(m0_aer) + '_r0_' + str(r0) + '.png'
    plt.savefig(join(savedir, fn))

#####
    # plt.show()

    plt.close('all')

if __name__ == '__main__':
    main()
