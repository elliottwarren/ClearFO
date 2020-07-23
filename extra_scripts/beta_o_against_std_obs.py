"""
Fast script to create a multi-plot of 1) observed backscatter with 2) RH and 3) wind/rain

Just showing that the ceilometer observations respond to changing meteorological conditions

beta_o smoothed
incoming solar radiation (for cloud spotting)
rainfall (rain spotting)
rh (aerosol swelling/shrinking)
"""

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from scipy.stats import spearmanr
from scipy.stats import pearsonr
from ellUtils import ellUtils as eu
import datetime as dt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
from matplotlib.dates import  DateFormatter

if __name__ == '__main__':

    # ---------------------------
    # Setup
    # ---------------------------

    # directories
    maindir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/'
    datadir = 'C:/Users/Elliott/Documents/PhD Reading/PhD Research/Aerosol Backscatter/clearFO/data/'
    savedir = maindir + 'figures/dailyplots/'

    daystrList = ['20170709']

    days_iterate = eu.dateList_to_datetime(daystrList)

    # --------------------------------
    # Read
    # --------------------------------

    # 1. observed attenuated backscatter


