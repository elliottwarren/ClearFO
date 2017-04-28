"""
Plot different FO variables with respect to m and/or RH.
Display sensitivity of each to FO input variables.

Created by Cristina Charlton-Perez - 07/12/2016
Edited by Elliott Warren - 28/04/2017
"""




import iris
import numpy as np
import matplotlib.pyplot as plt  # plotting library (not Met Office)
import math
from decimal import Decimal
from os.path import join


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


def main():

    # Path to save figures
    # savepath = '//data/dynamic/frpz/Aerosol_model_sensitivity'
    savedir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/figures/extraSensitivity/'

    # Conversion factors
    micro_g_2_kg = 1.0e-9
    kg_2_micro_g = 1.0e9
    m_2_microns = 1.0e6

    # critical value of relative humidity to start swelling the aerosol
    RH_crit = 38.0

    # Parameters from Pete Clarks (2008) visibility paper:
    rho_a = 1.7e3  # density of aerosol (for ammonium sulphate) [kg m-3]
    r0 = 1.1e-7   # radius of a 'standard' aerosol particle [m]

    p_aer = 1. / 6  # (p) power used to represent the variation in aerosol
    # particle size with mixing ratio
    # (p) power used to represent the variation in aerosol particle size
    # with mixing ratio UMDP 26 page 26 equation 39

    # Clark activation value
    B_activation = 0.5  # (B) Activation parameter (for ammonium sulphate)
    # Haywood activation value
    #B_activation = 0.14

    # Update from Haywood et.al. (2008) using flight data from around UK
    # this is the set of parameters used in the UM v8.1 (see UMDP 26)
    # to compute the aerosol number N_aer in code below
    N0_aer = 2.0e9  # (N0) standard number density of aerosol [m-3]
    # (m0) standard mass mixing ratio of the aerosol [kg kg-1]
    m0_aer = 1.8958e-8

    # N0_aer = 5.0e8 # (N0) standard number density of aerosol [m-3]
    # m0_aer = 1.6e-8 # (m0) standard mass mixing ratio of the aerosol [kg
    # kg-1]

    # For use in  eqns. 17-18 in Clark et.al. (2008):
    eta = 0.75
    Q_ext_aer = 2.0

    # m = np.arange(100) * 0.60 * micro_g_2_kg
    m = np.arange(100) * 0.60 * micro_g_2_kg
    RH = np.arange(100) * 0.01

    fig, ax1 = plt.subplots()

    ax1.plot(m * kg_2_micro_g, 'b-')
    ax1.plot(m / m0_aer, 'g-')
    ax1.set_xlabel('element')
    # Make the y-axis label and tick labels match the line color.
    ax1.set_ylabel('aerosol [micrograms kg-1]', color='b')
    for tl in ax1.get_yticklabels():
        tl.set_color('b')

    ax2 = ax1.twinx()
    ax2.plot(RH * 100.0, 'r.')
    ax2.set_ylabel('RH [%]', color='r')
    for tl in ax2.get_yticklabels():
        tl.set_color('r')
##########################
    fnum = 10
    # plt.figure(fnum)
    #
    # r_md = compute_dry(m, m0_aer, r0, p=1.0 / 6.0)
    # plt.plot(m, m_2_microns * r_md, 'k', label='p = 1/6')
    #
    # r_md = compute_dry(m, m0_aer, r0, p=1.0 / 4.0)
    # plt.plot(m, m_2_microns * r_md, 'c', label='p = 1/4')
    #
    # r_md = compute_dry(m, m0_aer, r0, p=1.0 / 3.0)
    # plt.plot(m, m_2_microns * r_md, 'b', label='p = 1/3')
    #
    # r_md = compute_dry(m, m0_aer, r0, p=1.0 / 2.0)
    # plt.plot(m, m_2_microns * r_md, 'm', label='p = 1/2')
    #
    # ax2 = plt.gca()
    # # plt.legend()
    # ax2.legend(loc='best', fancybox=True, framealpha=0.5)
    #
    # plt.xlabel("aerosol [kg kg^-1]")
    # plt.ylabel('r_m dry [microns]')
    # plt.grid()
    #
    # #titulo = 'Dry radius as a power function of aerosol with variable exponent p'
    # titulo = "r0 = " + str(r0) + " m0 = " + str(m0_aer)
    # print titulo
    # plt.title(titulo)

    fnum += 1
    plt.figure(fnum, figsize=(6, 3.5))

    r_md = compute_dry(m, m0_aer, r0, p=1.0 / 6.0)
    plt.plot(m * kg_2_micro_g, m_2_microns * r_md, 'k', label='p = 1/6')

    r_md = compute_dry(m, m0_aer, r0, p=1.0 / 4.0)
    plt.plot(m * kg_2_micro_g, m_2_microns * r_md, 'c', label='p = 1/4')

    r_md = compute_dry(m, m0_aer, r0, p=1.0 / 3.0)
    plt.plot(m * kg_2_micro_g, m_2_microns * r_md, 'b', label='p = 1/3')

    r_md = compute_dry(m, m0_aer, r0, p=1.0 / 2.0)
    plt.plot(m * kg_2_micro_g, m_2_microns * r_md, 'm', label='p = 1/2')

    ax2 = plt.gca()
    # plt.legend()
    ax2.legend(loc='best', fancybox=True, framealpha=0.5)

    # plt.xlabel('aerosol [micrograms kg^-1]')
    ax2.set_xlabel('aerosol [micrograms kg^-1]', labelpad=2)
    plt.ylabel('r_m dry [microns]')
    plt.grid()
    plt.tight_layout()

    #titulo = 'Dry radius as a power function of aerosol with variable exponent p'
    titulo = "r0 = " + str(r0) + " m0 = " + str(m0_aer)
    print titulo
    plt.title(titulo)

    fn2 = 'Vary_p_wrt_r0_' + str(N0_aer) + '_m0_' + str(m0_aer) + '_r0_' + str(r0) + '.png'
    plt.savefig(join(savedir, fn2))

    fnum += 1
    plt.figure(fnum)
    lw1 = 2.0

    N_aerosol = compute_N(N0_aer, m, m0_aer, p=1.0 / 9.0)
    plt.plot(m * kg_2_micro_g, N_aerosol, 'g',
             linewidth=lw1, label='p = 1/9')

    N_aerosol = compute_N(N0_aer, m, m0_aer, p=1.0 / 8.0)
    plt.plot(m * kg_2_micro_g, N_aerosol, 'y',
             linewidth=lw1, label='p = 1/8')

    N_aerosol = compute_N(N0_aer, m, m0_aer, p=1.0 / 6.0)
    plt.plot(m * kg_2_micro_g, N_aerosol, 'k',
             linewidth=lw1, label='p = 1/6')

    N_aerosol = compute_N(N0_aer, m, m0_aer, p=1.0 / 4.0)
    plt.plot(m * kg_2_micro_g, N_aerosol, 'c',
             linewidth=lw1, label='p = 1/4')

    N_aerosol = compute_N(N0_aer, m, m0_aer, p=1.0 / 3.0)
    plt.plot(m * kg_2_micro_g, N_aerosol, 'b',
             linewidth=lw1, label='p = 1/3')

    # N_aerosol = compute_N(N0_aer, m, m0_aer, p=1.0 / 2.0)
    # plt.plot(m * kg_2_micro_g, N_aerosol, 'r',
    #          linewidth=lw1, label='p = 1/2')

    ax2 = plt.gca()
    # plt.legend()
    ax2.legend(loc='best', fancybox=True, framealpha=0.5)

    plt.xlabel("aerosol [micrograms kg^-1]")
    plt.ylabel('N [m-3]')
    plt.grid()

    #titulo = 'Dry radius as a power function of aerosol with variable exponent p'
    N0_str = '%.2E' % Decimal(N0_aer)
    titulo = "N0 = " + N0_str + " m0 = " + str(m0_aer)
    print titulo
    plt.title(titulo)
    fn = 'Vary_p_' + str(N0_aer) + '_m0_' + str(m0_aer) + '_r0_' + str(r0) + '.png'
    plt.savefig(join(savedir, fn))

####
    r_md = compute_dry(m, m0_aer, r0, p=1.0 / 6.0)
    print len(m), len(RH)
    list_RH = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]
    m_matrix = np.ones((len(m), len(list_RH)))

    colors = ["r-", "m-",  "g-", "c-", "b-", "k-",
              "r--", "m--", "g--", "c--", "b--", "k--"]
    ncolor = 0
    fnum += 1
    plt.figure(fnum)
    for iRH in list_RH:
        print "iRH = ", iRH
        r_wet = compute_wet(
            r_md, iRH * np.ones(r_md.shape), B=0.5, p=1.0 / 3.0, RH_crit=0.38)
        m_matrix[:, iRH] = r_wet[:]

        print colors[ncolor]
        plt.plot(m * kg_2_micro_g, m_2_microns * r_wet,
                 colors[ncolor], label="RH=" + str(iRH))
        plt.ylabel('r_m wet radius [microns]')
        plt.xlabel('aerosol [micrograms per kg]')
        ncolor += 1

    # Make transparent legend
    ax = plt.gca()
    ax.legend(loc='best', fancybox=True, framealpha=0.5)

    plt.grid()
    titulo = "r0 = " + str(r0) + " m0 = " + str(m0_aer)
    print titulo
    plt.title(titulo)
    fn = 'Vary_RH_N0_' + str(N0_aer) + '_m0_' + str(m0_aer) + '_r0_' + str(r0) + '.png'
    # plt.savefig(savedir)
    plt.savefig(join(savedir, fn))

#####
    # plt.show()


if __name__ == '__main__':
    main()
