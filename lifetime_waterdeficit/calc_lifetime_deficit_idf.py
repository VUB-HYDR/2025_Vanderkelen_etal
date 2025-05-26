# Calculating script for liftetime exposure


import os
import xarray as xr
import pickle as pk
import time
import matplotlib as mpl
import matplotlib.colors as colors
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from settings import *
from load_manip import *
from exposure import * 
from utils import *
from main_function import *

# Script to plot countries
from plotting import *

# extreme event
global flags

flags = {}

flags['mode'] = 'frequency'   # whether to calculate exposure, duration, intensity or all. 

flags['runs'] = 0          # 0: do not process ISIMIP runs (i.e. load runs pickle)
                            # 1: process ISIMIP runs (i.e. produce and save runs as pickle)
flags['mask'] = 1           # 0: do not process country data (i.e. load masks pickle)
                            # 1: process country data (i.e. produce and save masks as pickle)
flags['exposure'] = 1     # 0: do not process ISIMIP runs to compute exposure (i.e. load exposure pickle)
                            # 1: process ISIMIP runs to compute exposure (i.e. produce and save exposure as pickle)
flags['exposure_pic'] = 0   # 0: do not process ISIMIP runs to compute picontrol exposure (i.e. load exposure pickle)
                            # 1: process ISIMIP runs to compute picontrol exposure (i.e. produce and save exposure as pickle)
flags['expectancy'] = 1     # 0: Load varying life expectancies
                            # 1: Constant life expectancies at 1960 levels

    
    
# select mask
#flags['mode'] = 'duration'   # whether to calculate exposure, duration, intensity or all. 
flags['mode'] = 'frequency'   # whether to calculate exposure, duration, intensity or all. 

# load lifetime waterdeficit duration
flags['extr'] = 'waterdeficitduration' 
#flags['extr'] = 'waterdeficitintensity' 

ds_duration, ds_duration_perregion, d_duration_allruns = do_lifetime_analysis(flags)

