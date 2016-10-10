

"""
Compute the attenuated modelled backscatter by calling the process_modelled_data() function

Created on Wed 14/09/16
"""

import cristinaCeilUtils as cutils
import numpy as np

# compute modelled backscatter

def aer_ext_scat(q_aer, RH, r0 = cutils.r0_haywood, B = cutils.B_activation_haywood):

    rho_a = 1.7e3  # density of aerosol (for ammonium sulphate) [kg m-3]

    r0 = cutils.r0_haywood  # radius of a 'standard' aerosol particle [m]
    # r0 = cutils.r0_clark

    # (p) power used to represent the variation in aerosol particle size with mixing ratio
    p = cutils.p_aer

    # Haywood (updated)
    # B = cutils.B_activation_clark  # (B) Activation parameter (for ammonium sulphate)
    B = cutils.B_activation_haywood

    # lidar ratio
    LR_d = cutils.LidarRatio_Dictionary()
    # lidar ratio in rain, to convert from extinction to backscatter
    S = LR_d['Aerosol'] # aerosol (60 sr)

    # Update from Haywood et.al. (2008) using flight data from around UK
    # this is the set of parameters used in the UM v8.1 (see UMDP 26)
    # to compute the aerosol number N_aer in code below
    N0 = cutils.N0_aer  # (N0) standard number density of aerosol
    m0 = cutils.m0_aer  # (m0) standard mass mixing ratio of the aerosol

    # For use in  eqns. 17-18 in Clark et.al. (2008):
    eta = cutils.eta
    Q = cutils.Q_ext_aer

    # Compute the aerosol number density N_aer. Eqn. 3 in Clark et.al. (2008) and Eqn 39 in UMDP 26 Large-scale precip.
    q_aer_kg_kg = q_aer * 1.0e-9  # convert micrograms kg-1 to kg/kg
    N_aer = N0 * np.power((q_aer_kg_kg / m0), 0.5)

    # r_md is the dry mean volume radius. Eqn. 2 in Clark et.al. (2008)
    r_md = r0 * np.power((q_aer_kg_kg / m0), p)

    # rm is the mean volume radius. Eqn. 12 in Clark et.al. (2008)
    # "When no activated particles are present an analytic solution for rm" is
    RH_crit = cutils.RH_crit
    #    print "RH_crit = ", RH_crit
    #    print RH[1:70, 1:4]
    # mask over values less than critical
    RH_ge_RHcrit = np.ma.masked_less(RH, RH_crit)

    RH_factor = 0.01  # Relative Humidity in 0.38 not 38%
    # RH_factor=1.0

    # eq 12 - calc rm for RH greater than critical
    rm = np.ma.ones(RH.shape) - (B / np.ma.log(RH_factor * RH_ge_RHcrit))
    rm2 = np.ma.power(rm, 1. / 3.)
    rm = np.ma.array(r_md) * rm2

    # set rm as 0 where RH is less than crit
    rm = np.ma.MaskedArray.filled(rm, [0.0])
    where_lt_crit = np.where(np.logical_or(RH.data < RH_crit, rm == 0.0))
    #    print 'where_lt_crit = ', where_lt_crit
    #    print type(where_lt_crit)
    ##     logical_and(my_array > 3, my_array < 7)
    #    if np.where(rm == 0.0):
    #        print "RM EQUALS ZERO!!!!!"
    #        print rm

    #    print 'Has it worked?', type(rm)
    rm[where_lt_crit] = r_md[where_lt_crit]

    temp = np.transpose(rm)

    # Close to activation one must solve the full equation (Kohler curve), not done in this code.
    # Assumptions made here include:
    # 1. Q_ext = scattering efficiency is independent of particle size and is assumed to be on average = 2.0
    # 2. Only time RH is taken into account is in equation 12 above in the
    # calculation of a mean radius. No attempt is made to model the variation
    # of particle growth with RH AS A FUNCTION OF SIZE of particle.

    # eqns. 17-18 in Clark et.al. (2008)
    alpha_a = (eta * Q) * np.pi * N_aer * np.power(rm, 2)
    print 'type of aerosol ext array; ', type(alpha_a)
    print alpha_a[0, :]

    beta_a = alpha_a / S
    # beta_a = alpha_a / 20.0 # lidar ratio for marine aerosol
    return alpha_a, beta_a


def compute_transmission(alpha_extinction, dz):

    ss = np.shape(alpha_extinction)
    # print "1. compute_transmission: size = ", ss[0],ss[1]

    integral_alpha = np.empty_like(alpha_extinction)
    transmission = np.empty_like(alpha_extinction)
    # print np.shape(integral_alpha), np.shape(transmission)

    # print  "2. compute_transmission: ",
    # np.shape(alpha_extinction),np.shape(dz)

    integral_alpha = np.cumsum(alpha_extinction * dz[0:len(dz)], axis=0)
    # print  "3. compute_transmission: ",alpha_extinction*dz[0:len(dz)],integral_alpha
    # transmission =  integral_alpha
    transmission = np.exp(-2.0 * integral_alpha)

    return integral_alpha, transmission


def process_modelled_data(aer_mod, rh_mod, z_mod):
    """
    Compute the attenuated modelled backscatter using beta and the transmissivity

    :param aer_mod:
    :param rh_mod:
    :param z_mod:
    :return: bsc_mod
    """

    # create alpha and beta coefficients for aerosol
    mod_alpha, mod_bscUnnAtt = aer_ext_scat(aer_mod, 100.0 * rh_mod)

    # TODO figure out how to make this better...
    NN = 28800 # number of enteries in time
    dz = np.zeros_like(z_mod)
    dz[0] = z_mod[0]
    dz[1:len(z_mod)] = z_mod[1:len(z_mod)] - z_mod[0:len(z_mod)-1]

    # integrated alpha and transmission for each height
    int_mod_alpha, mod_transm = compute_transmission(mod_alpha, dz)
    # int_mod_alpha, mod_transm = compute_transmission(mod_alpha[:, 0:NN], dz)

    # derive modelled attenuated backscatter
    bsc_mod = np.log10(mod_transm * mod_bscUnnAtt)
    return bsc_mod