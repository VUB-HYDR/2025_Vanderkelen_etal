# ---------------------------------------------------------------
# Utils functions for the water scarcity add-on
# ----------------------------------------------------------------

import numpy as np
import xarray as xr
import pandas as pd
import geopandas as gpd
import pickle as pk
from scipy import interpolate
import regionmask
import glob

from settings import *
from exposure import *
init()


# Calculate statistics for landfraction dictonaries
def calc_landfrac_stats(d_landfrac_peryear_perregion, var_tag):

    # reshape dictionary of data arrays into data arrays with runs as dimension
    da_landfrac_peryear_perregion = xr.concat([
        xr.concat([d_landfrac_peryear_perregion[x][v] for v in d_landfrac_peryear_perregion[x].keys()], dim='region') for x in  d_landfrac_peryear_perregion.keys()], dim='runs')

    # rename and order dimensions
    da_landfrac_peryear_perregion = da_landfrac_peryear_perregion.assign_coords({'region':list(d_landfrac_peryear_perregion[1].keys()), 'runs':list(d_landfrac_peryear_perregion.keys())}).rename({'time':'birth_year'}).transpose('runs', 'birth_year', 'region')
    
    # calculate stats in dataset
    ds_landfrac_stats = calc_exposure_mmm_xr(
    da_landfrac_peryear_perregion,
    'region',
    var_tag,
    flag_df2da = False
    )
    

    # take stats
    da_exposure_mmm_global = da_landfrac_peryear_perregion.mean(dim=('runs','region'))
    da_exposure_std_global = da_landfrac_peryear_perregion.std(dim=('runs','region'))
    da_exposure_lqntl_global = da_landfrac_peryear_perregion.quantile(
        q=0.25,
        dim=('runs','region'),
        method='inverted_cdf'
    )
    da_exposure_uqntl_global = da_landfrac_peryear_perregion.quantile(
        q=0.75,
        dim=('runs','region'),
        method='inverted_cdf'
    )
    
    # assemble into dataset
    ds_global_stats = xr.Dataset(
        data_vars={
            'mmm_global_{}'.format(var_tag): (['birth_year'],da_exposure_mmm_global.data),
            'std_global_{}'.format(var_tag): (['birth_year'],da_exposure_std_global.data),
            'lqntl_global_{}'.format(var_tag): (['birth_year'],da_exposure_lqntl_global.data),
            'uqntl_global_{}'.format(var_tag): (['birth_year'],da_exposure_uqntl_global.data),
        },
        coords={
            'birth_year': ('birth_year',da_exposure_uqntl_global.birth_year.data) }
    )
    
    return xr.merge([ds_landfrac_stats, ds_global_stats]).rename({'birth_year':'time'})

  