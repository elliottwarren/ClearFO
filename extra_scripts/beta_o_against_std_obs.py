"""
Fast script to create a multi-plot of 1) observed backscatter with 2) RH and 3) wind/rain

Just showing that the ceilometer observations respond to changing meteorological conditions

beta_o smoothed
incoming solar radiation (for cloud spotting)
rainfall (rain spotting)
rh (aerosol swelling/shrinking)
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
from matplotlib.dates import DateFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from ellUtils import ellUtils as eu
import datetime as dt
import numpy as np
from ceilUtils import ceilUtils as ceil
from forward_operator import FOconstants as FOcon

def read_met_obs(metdatadir, daystrList):
    """
    Read in the URAO met obs file (might need to prep file before hand to remove 2400 hr instance and turn into
    0000
    :return: met_obs (dict)
    """

    # filename
    met_obs_raw = eu.csv_read(metdatadir + daystrList[0] + '.csv')

    # sort dates
    YYYYMMDD = [met_obs_raw[i][0] for i in range(2, len(met_obs_raw))]
    HHMM = ['{:04d}'.format(int(met_obs_raw[i][1])) for i in range(2, len(met_obs_raw))]
    comb_date = [YYYYMMDD[i] + HHMM[i] for i in range(len(HHMM))]
    met_obs_time = np.array([dt.datetime(int(d[0:4]), int(d[4:6]), int(d[6:8]), int(d[8:10]), int(d[10:12]))
                             for d in comb_date])

    # sort data into a dictionary
    met_obs = {}
    for i, ob in enumerate(met_obs_raw[0]):
        met_obs[ob.replace(' ', '')] = [met_obs_raw[j][i] for j in range(2, len(met_obs_raw))]
    met_obs['time'] = met_obs_time

    return met_obs

def read_aer_obs(aerdatadir, day):
    """
    Read aerosol observations
    :param aerdatadir:
    :param day:
    :return:
    """

    filename = aerdatadir + '/' + day.strftime('%Y%m%d') + '_PM10_HollowayRoad.csv'
    aer_obs_raw = eu.csv_read(filename)

    dates_raw = [aer_obs_raw[i][2] for i in range(1, len(aer_obs_raw))]
    dates = [dt.datetime.strptime(i, '%d/%m/%Y %H:%M') for i in dates_raw]

    pm10_raw = [float(aer_obs_raw[i][3]) for i in range(1, len(aer_obs_raw))]
    aer_obs = {'pm10': pm10_raw, 'time': dates}

    return aer_obs

if __name__ == '__main__':

    # ---------------------------
    # Setup
    # ---------------------------

    # directories
    maindir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/'
    datadir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/data/'
    ceildatadir = datadir +'L1/'
    metdatadir = datadir +'URAO/'
    aerdatadir = datadir + 'ERG/'
    savedir = maindir + 'figures/dailyplots/'

    #daystrList = ['20170709'] # URAO
    # daystrList = ['20170808'] # IMU
    daystrList = ['20160913'] # IMU

    day = eu.dateList_to_datetime(daystrList)[0]

    # site = 'URAO'
    # ceil_id = 'CL31-F'
    # ceil_id_full = ceil_id + '_' + site
    # #site_bsc = {ceil_id_full: FOcon.site_bsc[ceil_id_full]}
    # site_bsc = {ceil_id_full: 1.5}

    site = 'IMU'
    ceil_id = 'CL31-A'
    ceil_id_full = ceil_id + '_' + site
    site_bsc = {ceil_id_full: FOcon.site_bsc[ceil_id_full]}

    ftype='BSC'

    # --------------------------------
    # Read
    # --------------------------------

    # 1. observed attenuated backscatter

    # times to match to, so the time between days will line up
    start = dt.datetime(day.year, day.month, day.day, 0, 0, 0)
    end = start + dt.timedelta(days=1)
    time_match = eu.date_range(start, end, 15, 'seconds')

    bsc_obs = ceil.read_all_ceils_BSC(day, site_bsc, ceildatadir, calib=False, timeMatch=time_match, var_type='beta_tR')

    #2. met obs
    if site == 'URAO':
        met_obs = read_met_obs(metdatadir, daystrList) # URAO
    elif site == 'IMU':
        met_obs = eu.netCDF_read(datadir+'/L1/Davis_IMU_'+day.strftime('%Y%j')+'_1min.nc') #IMU
        rad_obs = eu.netCDF_read(datadir+'/L1/PSP_IMU_'+day.strftime('%Y%j')+'_15min.nc') #IMU
        met_obs['RH'] = np.squeeze(met_obs['RH']) # get rid of unecessary dimensions
        rad_obs['Kdn'] = np.squeeze(rad_obs['Kdn'])

    # 3. aerosol obs
    # PM10 (if IMU site) - Holloway Road (nearby)
    if site == 'IMU':
        aer_obs = read_aer_obs(aerdatadir, day)

    # ----------------------
    # Process
    # ----------------------

    # only keep up to just above 3000 m
    _, idx, _ = eu.nearest(bsc_obs[ceil_id_full]['height'], 3000)
    bsc_obs[ceil_id_full]['height'] = bsc_obs[ceil_id_full]['height'][:idx+1]
    bsc_obs[ceil_id_full]['backscatter'] = bsc_obs[ceil_id_full]['backscatter'][:,:idx + 1]


    # -------------------------------
    # Plotting
    # -------------------------------

    fig, ax = plt.subplots(3, 1, figsize=(8, 6))

    # beta_o
    mesh1 = ax[0].pcolormesh(bsc_obs[ceil_id_full]['time'], bsc_obs[ceil_id_full]['height'],
                           np.transpose(bsc_obs[ceil_id_full]['backscatter']),
                           norm=LogNorm(vmin=1e-7, vmax=1e-5), cmap=cm.get_cmap('jet'))
    eu.add_at(ax[0], 'a)', loc=2, frameon=True)

    # URAO
    # ax[1].plot_date(met_obs['time'], met_obs['Sdw'], color='orange', linewidth=2, fmt='-', ls='-')
    # eu.add_at(ax[1], 'b)', loc=2, frameon=True)
    # ax[2].plot_date(met_obs['time'], met_obs['RH'], color='blue', linewidth=2, fmt='-', ls='-')
    # eu.add_at(ax[2], 'c)', loc=2, frameon=True)

    # ax[1].plot_date(rad_obs['time'], rad_obs['Kdn'], color='orange', linewidth=2, fmt='-', ls='-')
    # eu.add_at(ax[1], 'b)', loc=2, frameon=True)
    ax[1].plot_date(aer_obs['time'], aer_obs['pm10'], color='red', linewidth=2, fmt='-', ls='-')
    eu.add_at(ax[1], 'b)', loc=2, frameon=True)
    ax[2].plot_date(met_obs['time'], met_obs['RH'], color='blue', linewidth=2, fmt='-', ls='-')
    eu.add_at(ax[2], 'c)', loc=2, frameon=True)

    # prettify
    ax[0].autoscale(enable=True, axis='both', tight=True)

    # add colorbar
    divider = make_axes_locatable(ax[0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(mesh1, cax=cax)

    # labels
    for ax_i in ax:
        ax_i.xaxis.set_major_formatter(DateFormatter('%H:%M'))

    ax[0].set_ylabel('Height [m]', fontsize=13)
    # ax[1].set_ylabel('incoming\nshortwave '+r'$\mathrm{[W\/m_{-2}]}$', fontsize=13)
    ax[1].set_ylabel('PM' + r'$_{10}$' +' '+ r'$\mathrm{[\mu g\/\/ m^{-3}]}$', fontsize=13)
    ax[2].set_ylabel('RH [%]', fontsize=13)

    ax[0].xaxis.set_major_formatter(DateFormatter('%H:%M'))
    ax[-1].set_xlabel('Time [HH:MM]')

    # ensure radiation window spans 24 hours (like RH)
    ax[1].set_xlim([bsc_obs[ceil_id_full]['time'][0], bsc_obs[ceil_id_full]['time'][-1]])
    ax[2].set_xlim([bsc_obs[ceil_id_full]['time'][0], bsc_obs[ceil_id_full]['time'][-1]])

    plt.suptitle(site+' CL31')

    plt.tight_layout() # must come before axis resizing

    # adjust other subplots to be the same width as plot 1
    # [left, bottom, width, height]
    box1 = ax[0].get_position()
    for ax_i in ax[1:]:
        box2 = ax_i.get_position()
        ax_i.set_position([box2.x0, box2.y0,
                          box1.width, box1.height])

    plt.savefig(savedir+'/'+site+'_'+day.strftime('%Y%m%d')+'.png')