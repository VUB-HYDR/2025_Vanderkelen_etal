## main processing script

import os
import xarray as xr
import pickle as pk
import time
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl

from settings import *
from load_manip import *
from exposure import * 
from utils import *
from main_function import *

# Script to plot countries
from plotting import *

# set own plotting parameters
set_plot_param()


# extreme event
global flags

flags = {}
flags['extr'] = 'waterdeficit'   # 0: all
                                  # 2: cropfailedarea
                                  # 3: waterscarcity
                                  # 4: falkenmark 
                                  # 5: water deficit
flags['mode'] = 'exposure'   # whether to calculate exposure, duration, intensity or all. 

flags['runs'] = 0           # 0: do not process ISIMIP runs (i.e. load runs pickle)
                            # 1: process ISIMIP runs (i.e. produce and save runs as pickle)
flags['mask'] = 0           # 0: do not process country data (i.e. load masks pickle)
                            # 1: process country data (i.e. produce and save masks as pickle)
flags['exposure'] = 1     # 0: do not process ISIMIP runs to compute exposure (i.e. load exposure pickle)
                            # 1: process ISIMIP runs to compute exposure (i.e. produce and save exposure as pickle)
flags['exposure_pic'] = 0   # 0: do not process ISIMIP runs to compute picontrol exposure (i.e. load exposure pickle)
                            # 1: process ISIMIP runs to compute picontrol exposure (i.e. produce and save exposure as pickle)
# TODO: add rest of flags

flags['extr'] = 'waterdeficit'   # 0: all

ds_waterdeficit, ds_waterdeficit_perregion, ds_landfrac = do_lifetime_analysis(flags)

flags['extr'] = 'withdrawal'   # 0: all

ds_withdrawal, ds_withdrawal_perregion, ds_landfrac = do_lifetime_analysis(flags)