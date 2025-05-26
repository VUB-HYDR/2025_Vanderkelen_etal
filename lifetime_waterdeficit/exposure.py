# ---------------------------------------------------------------
# Functions to compute exposure
# ----------------------------------------------------------------


from operator import index
import numpy as np
import xarray as xr
import pandas as pd
import geopandas as gpd
import pickle as pk
from scipy import interpolate
import regionmask
import glob
import time
import matplotlib.pyplot as plt
from settings import *

#%% ----------------------------------------------------------------
# bootstrapping function 

def resample(
    da, 
    resample_dim,
    life_extent,
):
    """Resample with replacement in dimension ``resample_dim``. https://climpred.readthedocs.io/en/stable/_modules/climpred/bootstrap.html

    Args:
        initialized (xr.Dataset): input xr.Dataset to be resampled.
        resample_dim (str): dimension to resample along.
        life_extent (int): number of years per lifetime
        
    Returns:
        xr.Dataset: resampled along ``resample_dim``.

    """
    to_be_resampled = da[resample_dim].values
    smp = np.random.choice(to_be_resampled, life_extent)
    smp_da = da.sel({resample_dim: smp})
    smp_da[resample_dim] = np.arange(1960,1960+life_extent)
    return smp_da

#%% ----------------------------------------------------------------
# *improved function to compute extreme event exposure across a person's lifetime
def calc_life_exposure(
    df_exposure,
    df_life_expectancy,
    col,
):

    # initialisation variables
    phase, ages, age_young, age_ref, age_range, year_ref, year_start, birth_years, year_end, year_range, RCP2GMT_maxdiff_threshold, year_start_GMT_ref, year_end_GMT_ref, model_names = init()

    
    # initialise birth years 
    exposure_birthyears_percountry = np.empty(len(df_life_expectancy))

    for i, birth_year in enumerate(df_life_expectancy.index):

        life_expectancy = df_life_expectancy.loc[birth_year,col] 

        # define death year based on life expectancy
        death_year = birth_year + np.floor(life_expectancy)

        # integrate exposure over full years lived
        exposure_birthyears_percountry[i] = df_exposure.loc[birth_year:death_year,col].sum()

        # add exposure during last (partial) year
        exposure_birthyears_percountry[i] = exposure_birthyears_percountry[i] + \
            df_exposure.loc[death_year+1,col].sum() * \
                (life_expectancy - np.floor(life_expectancy))

    # a series for each column to somehow group into a dataframe
    exposure_birthyears_percountry = pd.Series(
        exposure_birthyears_percountry,
        index=df_life_expectancy.index,
        name=col,
    )

    return exposure_birthyears_percountry

#%% ----------------------------------------------------------------
# *improved function to compute extreme event average intensity across a person's lifetime
def calc_life_intensity_average(
    df_exposure,
    df_life_expectancy,
    col,
):

    # initialisation variables
    phase, ages, age_young, age_ref, age_range, year_ref, year_start, birth_years, year_end, year_range, RCP2GMT_maxdiff_threshold, year_start_GMT_ref, year_end_GMT_ref, model_names = init()

    
    # initialise birth years 
    exposure_birthyears_percountry = np.empty(len(df_life_expectancy))

    for i, birth_year in enumerate(df_life_expectancy.index):

        life_expectancy = df_life_expectancy.loc[birth_year,col] 

        # define death year based on life expectancy
        death_year = birth_year + np.floor(life_expectancy)

        # integrate exposure over full years lived
        exposure_birthyears_percountry[i] = df_exposure.loc[birth_year:death_year,col].mean()

    # a series for each column to somehow group into a dataframe
    exposure_birthyears_percountry = pd.Series(
        exposure_birthyears_percountry,
        index=df_life_expectancy.index,
        name=col,
    )

    return exposure_birthyears_percountry

#%% ----------------------------------------------------------------
# Function to load the dryland mask based on the simulation
def load_drylandmask_old(
    year_start,
    year_end,
    model,
    forcing,
    rcpscenario):
    
    
    # initialisation variables
    phase, ages, age_young, age_ref, age_range, year_ref, year_start, birth_years, year_end, year_range, RCP2GMT_maxdiff_threshold, year_start_GMT_ref, year_end_GMT_ref, model_names = init()

    # load 2D model constants
    da_drylandmask_hist = xr.open_dataset('./data/isimip/2b/dryland_mask/'+model.lower()+'_'+forcing+'_historical_drylandmask_global_annual_1861_2005.nc4')['mask']
    
    da_drylandmask_fut = xr.open_dataset('./data/isimip/2b/dryland_mask/'+model.lower()+'_'+forcing+'_'+rcpscenario+'_drylandmask_global_annual_2006_2099.nc4')['mask']

    # manually adjust time dimension in both data arrays (because original times could not be decoded)
    #da_population_histsoc['time'] = np.arange(1861,2006)
    #da_population_ssp2soc['time'] = np.arange(2006,2101)
    # concatenate historical and future data
    da_drylandmask = xr.concat([da_drylandmask_hist, da_drylandmask_fut], dim='year') 


    # if needed, repeat last year until entire period of interest is covered
    if np.nanmax(da_drylandmask.year) < year_end:
        population_10y_mean = da_drylandmask.loc[-10:,:,:].mean(dim='year').expand_dims(dim='year',axis=0) # repeat average of last 10 years (i.e. end-9 to end ==> 2090:2099)
        population_10y = population_10y_mean > 0.5
        for year in range(np.nanmax(da_drylandmask.year)+1,year_end+1): 
            da_drylandmask = xr.concat([da_drylandmask,population_10y.assign_coords(year = [year])], dim='year')

    # retain only period of interest
    da_drylandmask = da_drylandmask.sel(year=slice(year_start,year_end))

    return da_drylandmask


#%% ----------------------------------------------------------------
# calculated weighted fieldmean per country mask
def calc_weighted_fldmean(
    da, 
    weights, 
    countries_mask, 
    ind_country, 
    flag_region,
):

    # one country provided, easy masking
    if not flag_region : 
        da_masked = da.where(countries_mask == ind_country)
    
    # if more countries are provided, combine the different masks 
    else: 
        
        if len(ind_country) > 1:
            
            mask = xr.DataArray(
                np.in1d(countries_mask,ind_country).reshape(countries_mask.shape),
                dims=countries_mask.dims,
                coords=countries_mask.coords,
            )
            da_masked = da.where(mask)
    
    da_weighted_fldmean = da_masked.weighted(weights).mean(dim=("lat", "lon"))
    
    return da_weighted_fldmean


#%% ----------------------------------------------------------------
# get member countries per region
def get_countries_of_region(
    region, 
    df_countries,
): 

    # Get list of member countries from region
    member_countries = df_countries.loc[df_countries['region']==region]['name'].values

    # not region but income group
    if len(member_countries) == 0: 
        member_countries = df_countries.loc[df_countries['incomegroup']==region]['name'].values

    # get all countries for the world
    if region == 'World':
        member_countries = df_countries['name'].values

    return member_countries    

#%% ----------------------------------------------------------------
# function to compute multi-model mean across ISIMIP simulations based on mf_exposure_mmm.m
def calc_exposure_mmm_xr(
    d_exposure,
    dim_1_name,
    var_tag,
    flag_df2da = True
):
    
    if flag_df2da:    
        # change hist+RCP exposure dictionary of runs to data array
        da_exposure = xr.concat(
            [xr.DataArray(v).rename({'dim_0':'birth_year','dim_1':dim_1_name}) for v in d_exposure.values()],
            dim='runs',
        ).assign_coords({'runs':list(d_exposure.keys())})
    else: 
        da_exposure = d_exposure

    # take stats
    da_exposure_mmm = da_exposure.mean(dim='runs')
    da_exposure_std = da_exposure.std(dim='runs')
    da_exposure_lqntl = da_exposure.quantile(
        q=0.25,
        dim='runs',
        method='inverted_cdf'
    )
    da_exposure_uqntl = da_exposure.quantile(
        q=0.75,
        dim='runs',
        method='inverted_cdf'
    )

    # get EMF of stats (divide by 60 yr old / 1960 cohort)
    #da_exposure_mmm_EMF = da_exposure_mmm / da_exposure_mmm.sel(birth_year=1960)
    #da_exposure_lqntl_EMF = da_exposure_lqntl / da_exposure_mmm.sel(birth_year=1960)
    #da_exposure_uqntl_EMF = da_exposure_uqntl / da_exposure_mmm.sel(birth_year=1960)

    # get delta's first and then calc statistics
    da_exposure_mmm_delta = np.repeat(np.expand_dims((da_exposure.sel(birth_year=2020) - da_exposure.sel(birth_year=1960)).mean(dim='runs'), axis=0),len(da_exposure.birth_year) ,axis=0)
    da_exposure_std_delta= np.repeat(np.expand_dims((da_exposure.sel(birth_year=2020) - da_exposure.sel(birth_year=1960)).std(dim='runs'), axis=0),len(da_exposure.birth_year) ,axis=0)
    da_exposure_lqntl_delta = np.repeat(np.expand_dims((da_exposure.sel(birth_year=2020) - da_exposure.sel(birth_year=1960)).quantile(
        q=0.25,
        dim='runs',
        method='inverted_cdf'
    ), axis=0),len(da_exposure.birth_year) ,axis=0)
    da_exposure_uqntl_delta = np.repeat(np.expand_dims((da_exposure.sel(birth_year=2020) - da_exposure.sel(birth_year=1960)).quantile(
        q=0.75,
        dim='runs',
        method='inverted_cdf'
    ), axis=0),len(da_exposure.birth_year) ,axis=0)



    # get delta of stats (subtract the 60 yr old / 1960 cohort)
    #da_exposure_mmm_delta = da_exposure_mmm - da_exposure_mmm.sel(birth_year=1960)
    #da_exposure_lqntl_delta = da_exposure_lqntl - da_exposure_mmm.sel(birth_year=1960)
    #da_exposure_uqntl_delta = da_exposure_uqntl - da_exposure_mmm.sel(birth_year=1960)    
    #da_exposure_std_delta = da_exposure_std - da_exposure_mmm.sel(birth_year=1960)    

    # do the same for gobal variables
    da_exposure_mmm_global = da_exposure.mean(dim=('runs',dim_1_name))
    da_exposure_std_global = da_exposure.std(dim=('runs',dim_1_name))
    da_exposure_lqntl_global = da_exposure.quantile(
        q=0.25,
        dim=('runs',dim_1_name),
        method='inverted_cdf'
    )
    da_exposure_uqntl_global = da_exposure.quantile(
        q=0.75,
        dim=('runs',dim_1_name),
        method='inverted_cdf'
    )

    # get EMF of stats (divide by 60 yr old / 1960 cohort)
    #da_exposure_mmm_EMF_global = da_exposure_mmm_global / da_exposure_mmm_global.sel(birth_year=1960)
    #da_exposure_lqntl_EMF_global = da_exposure_lqntl_global / da_exposure_mmm_global.sel(birth_year=1960)
    #da_exposure_uqntl_EMF_global = da_exposure_uqntl_global / da_exposure_mmm_global.sel(birth_year=1960)

    # get delta of stats (subtract the 60 yr old / 1960 cohort)
    da_exposure_mmm_delta_global = da_exposure_mmm_global - da_exposure_mmm_global.sel(birth_year=1960)
    da_exposure_lqntl_delta_global = da_exposure_lqntl_global - da_exposure_mmm_global.sel(birth_year=1960)
    da_exposure_uqntl_delta_global = da_exposure_uqntl_global - da_exposure_mmm_global.sel(birth_year=1960)    
    da_exposure_std_delta_global = da_exposure_std_global - da_exposure_mmm_global.sel(birth_year=1960)    


    # assemble into dataset
    ds_exposure_stats = xr.Dataset(
        data_vars={
            'mmm_{}'.format(var_tag): (['birth_year',dim_1_name],da_exposure_mmm.data),
            'std_{}'.format(var_tag): (['birth_year',dim_1_name],da_exposure_std.data),
            'lqntl_{}'.format(var_tag): (['birth_year',dim_1_name],da_exposure_lqntl.data),
            'uqntl_{}'.format(var_tag): (['birth_year',dim_1_name],da_exposure_uqntl.data),
            'mmm_delta_{}'.format(var_tag): (['birth_year',dim_1_name],da_exposure_mmm_delta.data),
            'lqntl_delta_{}'.format(var_tag): (['birth_year',dim_1_name],da_exposure_lqntl_delta.data),
            'uqntl_delta_{}'.format(var_tag): (['birth_year',dim_1_name],da_exposure_uqntl_delta.data),
            'std_delta_{}'.format(var_tag): (['birth_year',dim_1_name],da_exposure_std_delta.data),            
            'mmm_global_{}'.format(var_tag): (['birth_year'],da_exposure_mmm_global.data),
            'std_global_{}'.format(var_tag): (['birth_year'],da_exposure_std_global.data),
            'lqntl_global_{}'.format(var_tag): (['birth_year'],da_exposure_lqntl_global.data),
            'uqntl_global_{}'.format(var_tag): (['birth_year'],da_exposure_uqntl_global.data),
            'mmm_delta_global_{}'.format(var_tag): (['birth_year'],da_exposure_mmm_delta_global.data),
            'lqntl_delta_global_{}'.format(var_tag): (['birth_year'],da_exposure_lqntl_delta_global.data),
            'uqntl_delta_global_{}'.format(var_tag): (['birth_year'],da_exposure_uqntl_delta_global.data),            
            'std_delta_global_{}'.format(var_tag): (['birth_year'],da_exposure_std_delta_global.data),            
        },
        coords={
            'birth_year': ('birth_year',da_exposure.birth_year.data),
            dim_1_name: (dim_1_name,da_exposure[dim_1_name].data)}
    )

    return ds_exposure_stats

def calc_exposure_mmm_mitigation_xr(d_exposure_rcp26runs, d_exposure_rcp60runs, 
    dim_1_name,
    var_tag,
    flag_df2da = True, 
                                    
):

    if flag_df2da:    
        # change hist+RCP exposure dictionary of runs to data array
        da_exposure_26 = xr.concat(
            [xr.DataArray(v).rename({'dim_0':'birth_year','dim_1':dim_1_name}) for v in d_exposure_rcp26runs.values()],
            dim='runs').assign_coords({'runs':list(d_exposure_rcp26runs.keys())})
        da_exposure_60 = xr.concat(
            [xr.DataArray(v).rename({'dim_0':'birth_year','dim_1':dim_1_name}) for v in d_exposure_rcp60runs.values()],
            dim='runs').assign_coords({'runs':list(d_exposure_rcp60runs.keys())})    

    else: 
        da_exposure_26 = d_exposure_rcp26runs
        da_exposure_60 = d_exposure_rcp60runs

    # take stats
    da_mitigation = da_exposure_60
    da_mitigation.values = da_exposure_60.values - da_exposure_26.values

    da_exposure_mmm = da_mitigation.mean(dim='runs')
    da_exposure_std = da_mitigation.std(dim='runs')
    da_exposure_lqntl = da_mitigation.quantile(
        q=0.25,
        dim='runs',
        method='inverted_cdf'
    )
    da_exposure_uqntl = da_mitigation.quantile(
        q=0.75,
        dim='runs',
        method='inverted_cdf'
    )


    # assemble into dataset
    ds_exposure_stats_mitigation = xr.Dataset(
        data_vars={
            'mmm_mitigation_{}'.format(var_tag): (['birth_year',dim_1_name],da_exposure_mmm.data),
            'std_mitigation_{}'.format(var_tag): (['birth_year',dim_1_name],da_exposure_std.data),
            'lqntl_mitigation_{}'.format(var_tag): (['birth_year',dim_1_name],da_exposure_lqntl.data),
            'uqntl_mitigation_{}'.format(var_tag): (['birth_year',dim_1_name],da_exposure_uqntl.data)},
        coords={
            'birth_year': ('birth_year',da_exposure_60.birth_year.data),
            dim_1_name: (dim_1_name,da_exposure_60[dim_1_name].data)}
    )
    
    return ds_exposure_stats_mitigation
#%% ----------------------------------------------------------------
# function to compute multi-model mean across ISIMIP simulations based on mf_exposure_mmm.m
def calc_exposure_mmm_pic_xr(
    d_exposure_pic,
    dim_1_name,
    var_tag,
):        
        
    # concat pic data array from dict of separate arrays
    da_exposure_pic = xr.concat(
        [v for v in d_exposure_pic.values()],
        dim='runs',    
    ).assign_coords({'runs':list(d_exposure_pic.keys())})

    # runs and lifetimes (lifetimes from boostrapping) redundant, so compile together
    da_exposure_pic = da_exposure_pic.stack(
        pic_lifetimes=['runs','lifetimes'],
    )

    # pic exposure stats for EMF
    da_exposure_pic_mmm = da_exposure_pic.mean(dim='pic_lifetimes')
    da_exposure_pic_std = da_exposure_pic.std(dim='pic_lifetimes')
    da_exposure_pic_lqntl = da_exposure_pic.quantile(
        q=0.25,
        dim='pic_lifetimes',
        method='inverted_cdf',
    )
    da_exposure_pic_uqntl = da_exposure_pic.quantile(
        q=0.75,
        dim='pic_lifetimes',
        method='inverted_cdf',
    )    

    # pic quantile for birth cohort exposure emergence
    da_exposure_pic_ext = da_exposure_pic.quantile(
        q=0.99,
        dim='pic_lifetimes',
        # method='inverted_cdf',
    )
    
    # assemble into dataset
    ds_exposure_pic_stats = xr.Dataset(
        data_vars={
            'mmm_{}'.format(var_tag): ([dim_1_name],da_exposure_pic_mmm.data),
            'std_{}'.format(var_tag): ([dim_1_name],da_exposure_pic_std.data),
            'lqntl_{}'.format(var_tag): ([dim_1_name],da_exposure_pic_lqntl.data),
            'uqntl_{}'.format(var_tag): ([dim_1_name],da_exposure_pic_uqntl.data),
            'ext': ([dim_1_name],da_exposure_pic_ext.data),
        },
        coords={
            dim_1_name: (dim_1_name,da_exposure_pic[dim_1_name].data),
        }
    )

    return ds_exposure_pic_stats


#%% ----------------------------------------------------------------
# convert Area Fraction Affected (AFA) to 
# per-country number of extremes affecting one individual across life span
def calc_exposure(
    grid_area,
    d_regions,
    d_isimip_meta, 
    df_birthyears_regions, 
    df_countries, 
    countries_regions, 
    countries_mask, 
    da_population, 
    df_life_expectancy_5,
    d_all_cohorts,
    modes,

):
        phase, ages, age_young, age_ref, age_range, year_ref, year_start, birth_years, year_end, year_range, RCP2GMT_maxdiff_threshold, year_start_GMT_ref, year_end_GMT_ref, model_names = init()

        d_exposure_perrun_RCP     = {}
        d_exposure_perrun_15      = {}
        d_exposure_perrun_20      = {}
        d_exposure_perrun_NDC     = {}
        d_exposure_perrun_R26eval = {} 
        d_cohort_exposure = {}
        
        d_landfrac_peryear_perregion = {}
        d_exposure_perregion_perrun_RCP = {}
        
        # unpack region information
        df_birthyears_regions = d_regions['birth_years']
        d_cohort_weights_regions = d_regions['cohort_size']    
        da_cohort_size = xr.DataArray(
            # np.asarray(list(d_cohort_size.values())),
            np.asarray([v for k,v in d_all_cohorts.items() if k in list(df_countries['name'])]),
            coords={
                'country': ('country', list(df_countries['name'])),
                'time': ('time', year_range),
                'ages': ('ages', np.arange(104,-1,-1)),
            },
            dims=[
                'country',
                'time',
                'ages',
            ]
        ) 
        
        
        # --------------------------------------------------------------------
        # load mask and apply it on population dataset
        
        if modes[0] in ['dryland', 'arid', 'hyperarid', 'semiarid', 'subhumid', 'popdensity1000', 'popdensity500', 'temperatedrysummer', 'temperatedrywinter', 'aridsteppe']:
            mask = xr.open_dataset('./data/masks/'+modes[0]+'.nc')['mask']
            da_population = da_population.where(mask).fillna(0)
        elif modes[0] in ['urban', 'rural']:
            popshare = xr.open_dataset('./data/masks/'+modes[0]+'share_2020.nc', decode_times=False )['mask']
            da_population = (da_population * popshare).fillna(0)
            # to be removed
            da_population.to_dataset(name='mask').to_netcdf('./data/masks/'+modes[0]+'_popshare_tocheck.nc')
            

        
        
        # loop over simulations
        for i in list(d_isimip_meta.keys()): 

            print('simulation '+str(i)+ ' of '+str(len(d_isimip_meta)))
            print(d_isimip_meta[i]['extreme'])

            # load AFA data of that run
            with open('./data/pickles/isimip_AFA_{}_{}_{}.pkl'.format(d_isimip_meta[i]['extreme'],d_isimip_meta[i]['mode'],str(i)), 'rb') as f:
                da_AFA = pd.read_pickle(f)
          
            # --------------------------------------------------------------------
            # apply dryland mask --- OLD with dryland mask based on ISIMIP
            
            #if modes=='dryland': 
                # load dryland mask related to simulation
            #     dryland_mask = load_drylandmask(year_start, year_end, d_isimip_meta[i]['model'], d_isimip_meta[i]['gcm'], d_isimip_meta[i]['rcp'])
            #     da_population = da_population.where(dryland_mask).fillna(0)
                
            # --------------------------------------------------------------------
            # per country 

            # initialise dicts
            d_exposure_peryear_percountry = {}

            # get spatial average
            for j, country in enumerate(df_countries['name']):

                print('processing country '+str(j+1)+' of '+str(len(df_countries)), end='\r')
                
                # calculate mean per country weighted by population
                ind_country = countries_regions.map_keys(country)

                # historical + RCP simulations
                d_exposure_peryear_percountry[country] = calc_weighted_fldmean( 
                    da_AFA,
                    da_population, 
                    countries_mask, 
                    ind_country, 
                    flag_region= False,
                )
                
            da_exposure_peryear_percountry = xr.DataArray(
                list(d_exposure_peryear_percountry.values()),
                coords={
                    'country': ('country', list(d_exposure_peryear_percountry.keys())),
                    'time': ('time', da_AFA.time.values),
                },
                dims=[
                    'country',
                    'time',
                ],
            )

            d_cohort_exposure[i]= da_exposure_peryear_percountry * da_cohort_size
                
            # --------------------------------------------------------------------
            # convert dict to dataframe for vectorizing and integrate exposures then map to GMTs  
            
            frame = {k:v.values for k,v in d_exposure_peryear_percountry.items()}
            df_exposure = pd.DataFrame(frame,index=np.arange(1960,2114))           

            # apply calc life exposure to columns (countries) of df_exposure
            d_exposure_perrun_RCP[i] = df_exposure.apply(
                lambda col: calc_life_exposure(
                    df_exposure,
                    df_life_expectancy_5,
                    col.name,
                ),
                axis=0,
            )
            
            if modes[0] == 'intensity': # not exposure (sum) but intensity (average)
                # apply calc life exposure to columns (countries) of df_exposure
                d_exposure_perrun_RCP[i] = df_exposure.apply(
                    lambda col: calc_life_intensity_average(
                        df_exposure,
                        df_life_expectancy_5,
                        col.name,
                    ),
                    axis=0,
                )
            
            
            # if max threshold criteria met, run gmt mapping
            if d_isimip_meta[i]['GMT_15_valid']:
                
                d_exposure_perrun_15[i] = df_exposure.apply(
                    lambda col: calc_life_exposure(
                        df_exposure.reindex(df_exposure.index[d_isimip_meta[i]['ind_RCP2GMT_15']]).set_index(df_exposure.index),
                        df_life_expectancy_5,
                        col.name,
                    ),
                    axis=0,
                )
                
            if d_isimip_meta[i]['GMT_20_valid']:
                
                d_exposure_perrun_20[i] = df_exposure.apply(
                    lambda col: calc_life_exposure(
                        df_exposure.reindex(df_exposure.index[d_isimip_meta[i]['ind_RCP2GMT_20']]).set_index(df_exposure.index),
                        df_life_expectancy_5,
                        col.name,
                    ),
                    axis=0,
                )
                
            if d_isimip_meta[i]['GMT_NDC_valid']:
                
                d_exposure_perrun_NDC[i] = df_exposure.apply(
                    lambda col: calc_life_exposure(
                        df_exposure.reindex(df_exposure.index[d_isimip_meta[i]['ind_RCP2GMT_NDC']]).set_index(df_exposure.index),
                        df_life_expectancy_5,
                        col.name,
                    ),
                    axis=0,
                )
            
            if d_isimip_meta[i]['GMT_R26eval_valid']:
                
                d_exposure_perrun_R26eval[i] = df_exposure.apply(
                    lambda col: calc_life_exposure(
                        df_exposure.reindex(df_exposure.index[d_isimip_meta[i]['ind_RCP2GMT_R26eval']]).set_index(df_exposure.index),
                        df_life_expectancy_5,
                        col.name,
                    ),
                    axis=0,
                )
            
            # --------------------------------------------------------------------
            # per region
            #  

            print('')

            # initialise dictionaries
            d_landfrac_peryear_perregion[i] = {}
            d_exposure_perregion_RCP = {}

            # loop over regions
            for k, region in enumerate(df_birthyears_regions.columns): 
                
                print('processing region '+str(k+1)+' of '+str(len(df_birthyears_regions.columns)), end='\r')

                # Get list of member countries from region - with seperate treatment for world (luke: now inside get_countries_of_regions func)
                member_countries = get_countries_of_region(region, df_countries)
        
                # get spatial average of landfraction: historical + RCP simulations
                ind_countries = countries_regions.map_keys(member_countries)

                #print('calculating landfrac')
                d_landfrac_peryear_perregion[i][region] = calc_weighted_fldmean(
                    da_AFA, 
                    grid_area, 
                    countries_mask, 
                    ind_countries, 
                    flag_region=True,
                )

                #print('calculating cohort weights')
                # filter cohort weights to only keep countries within mask 
                d_cohort_weights_regions[region] = d_cohort_weights_regions[region].loc[:,d_cohort_weights_regions[region].columns.isin(df_countries.index)]
                
                # get weighted spatial average for all member countries per region (Luke: but this is not necessarily spatial; d_cohorts_weights_regions has 2020 cohort sizes per country)
                d_exposure_perregion_RCP[region] = (d_exposure_perrun_RCP[i].loc[:,member_countries] * d_cohort_weights_regions[region].values).sum(axis=1) /\
                    np.nansum(d_cohort_weights_regions[region].values, axis=1)
 
            print('')
            
            # save exposures for every run
            d_exposure_perregion_perrun_RCP[i]  = pd.DataFrame(d_exposure_perregion_RCP)
            
        da_exposure_cohort = xr.concat(
            [v for v in d_cohort_exposure.values()],
            dim='runs',
        ).assign_coords({'runs':list(d_cohort_exposure.keys())})

        # --------------------------------------------------------------------
        # save workspave in pickles
        #  

        # save pickles
        print()
        print('Saving processed exposures')

        # pack exposure information
        d_exposure = {
            'exposure_perrun_RCP' : d_exposure_perrun_RCP, 
            'exposure_perrun_15' : d_exposure_perrun_15,
            'exposure_perrun_20' : d_exposure_perrun_20,
            'exposure_perrun_NDC' : d_exposure_perrun_NDC,
            'exposure_perregion_perrun_RCP' : d_exposure_perregion_perrun_RCP, 
            'landfrac_peryear_perregion' : d_landfrac_peryear_perregion,
            'exposure_per_cohort': da_exposure_cohort,
        }

        with open('./data/pickles/exposure_{}_{}.pkl'.format(d_isimip_meta[i]['extreme'],d_isimip_meta[i]['mode']), 'wb') as f:
            pk.dump(d_exposure,f)

        return d_exposure_perrun_RCP, d_exposure_perregion_perrun_RCP, d_exposure_perrun_15, d_exposure_perrun_20, d_exposure_perrun_NDC, da_exposure_cohort
        
#%% ----------------------------------------------------------------
# convert PIC Area Fraction Affected (AFA) to 
# per-country number of extremes affecting one individual across life span
def calc_exposure_pic(
    grid_area,
    d_regions,
    d_pic_meta, 
    df_birthyears_regions, 
    df_countries, 
    countries_regions, 
    countries_mask, 
    da_population, 
    df_life_expectancy_5, 
):

        d_exposure_perrun_pic = {}       
        d_landfrac_peryear_perregion_pic = {}
        d_exposure_perregion_perrun_pic = {}
        
        # unpack region information
        df_birthyears_regions = d_regions['birth_years']
        d_cohort_weights_regions = d_regions['cohort_size']                
        
        # loop over simulations
        for n,i in enumerate(list(d_pic_meta.keys())):

            print('simulation '+str(n+1)+ ' of '+str(len(d_pic_meta)))

            # load AFA data of that run
            with open('./data/pickles/isimip_AFA_pic_{}_{}.pkl'.format(d_pic_meta[i]['extreme'],str(i)), 'rb') as f:
                da_AFA_pic = pd.read_pickle(f)
            
            # get 1960 life expectancy
            life_expectancy_1960 = xr.DataArray(
                df_life_expectancy_5.loc[1960].values,
                coords={
                    'country': ('country', df_life_expectancy_5.columns)
                }
            )            
            
            # --------------------------------------------------------------------
            # per country 
            # start_time = time.time()
            d_exposure_peryear_percountry_pic = {}
            
            # get spatial average
            for j, country in enumerate(df_countries['name']): # with other stuff running, this loop took 91 minutes
                # therefore consider first doing the weighted mean and then boot strapping? does that make sense?

                print('processing country '+str(j+1)+' of '+str(len(df_countries)), end='\r')
                # calculate mean per country weighted by population
                ind_country = countries_regions.map_keys(country)

                # corresponding picontrol - assume constant 1960 population density (this line takes about 16h by itself)
                d_exposure_peryear_percountry_pic[country] = calc_weighted_fldmean(
                    da_AFA_pic, 
                    da_population[0,:,:], # earliest year used for weights
                    countries_mask, 
                    ind_country, 
                    flag_region= False,
                )
                
            da_exposure_pic = xr.DataArray(
                list(d_exposure_peryear_percountry_pic.values()),
                coords={
                    'country': ('country', list(d_exposure_peryear_percountry_pic.keys())),
                    'time': ('time', da_AFA_pic.time.values),
                },
                dims=[
                    'country',
                    'time',
                ],
            )

            # bootstrap native pic exposed area data
            life_extent=82 # max 1960 life expectancy is 81, therefore bootstrap lifetimes of 82 years
            nboots=100 # number of lifetimes to be bootstrapped should go to settings
            resample_dim='time'
            da_exposure_pic = xr.concat([resample(da_exposure_pic,resample_dim,life_extent) for i in range(nboots)],dim='lifetimes')
            
            # save pic2 data for checking
            with open('./data/pickles/da_exposure_pic.pkl', 'wb') as f: # note; 'with' handles file stream closing
                pk.dump(da_exposure_pic,f)            
            
            # --------------------------------------------------------------------
            # substitute calc_life_exposure because we are only doing the 1960 cohort
            d_exposure_perrun_pic[i] = da_exposure_pic.where(da_exposure_pic.time < 1960 + np.floor(life_expectancy_1960)).sum(dim='time') + \
                da_exposure_pic.where(da_exposure_pic.time == 1960 + np.floor(life_expectancy_1960)).sum(dim='time') * \
                    (life_expectancy_1960 - np.floor(life_expectancy_1960))

            # --------------------------------------------------------------------
            # per region
            #  

            print('')

            # initialise dictionaries
            d_landfrac_peryear_perregion_pic[i] = {}
            d_exposure_perregion_pic = {}

            # loop over regions
            for k, region in enumerate(df_birthyears_regions.columns): 
                
                print('processing region '+str(k+1)+' of '+str(len(df_birthyears_regions.columns)), end='\r')

                # Get list of member countries from region - with seperate treatment for world (luke: now inside get_countries_of_regions func)
                member_countries = get_countries_of_region(region, df_countries)
        
                # get spatial average of landfraction: historical + RCP simulations
                ind_countries = countries_regions.map_keys(member_countries)

                #print('calculating landfrac') # don't need to bootstrap more samples of lifetimes from here because this is just PIC landfrac affected
                d_landfrac_peryear_perregion_pic[i][region] = calc_weighted_fldmean(
                    da_AFA_pic, 
                    grid_area, 
                    countries_mask, 
                    ind_countries, 
                    flag_region=True,
                )

                #print('calculating cohort weights')
                # filter cohort weights to only keep countries within mask 
                d_cohort_weights_regions[region] = d_cohort_weights_regions[region].loc[:,d_cohort_weights_regions[region].columns.isin(df_countries.index)]
                da_cwr_1960 = xr.DataArray(
                    d_cohort_weights_regions[region].loc[60],
                    coords={
                        'country': ('country', member_countries),
                    }
                )
                
                # get weighted spatial average for all member countries per region
                d_exposure_perregion_pic[region] = (d_exposure_perrun_pic[i].sel(country=member_countries) * da_cwr_1960).sum(axis=1) /\
                    da_cwr_1960.sum(dim='country')

            # save exposures for every run
            df_exposure_perregion_pic = pd.DataFrame(d_exposure_perregion_pic)
            da_exposure_perregion_pic = xr.DataArray(df_exposure_perregion_pic).rename({'dim_0':'lifetimes','dim_1':'region'})
            d_exposure_perregion_perrun_pic[i] = da_exposure_perregion_pic
            print('')
            
        # --------------------------------------------------------------------
        # save workspave in pickles
        #  

        # save pickles
        print()
        print('Saving processed exposures')

        # pack region information
        d_exposure = {
            'exposure_perrun' : d_exposure_perrun_pic, 
            'exposure_perregion_perrun' : d_exposure_perregion_perrun_pic, 
            'landfrac_peryear_perregion' : d_landfrac_peryear_perregion_pic 
        }

        with open('./data/pickles/exposure_pic_{}.pkl'.format(d_pic_meta[i]['extreme']), 'wb') as f:
            pk.dump(d_exposure,f)

        return d_exposure_perrun_pic, d_exposure_perregion_perrun_pic

    # calculate percentage of lifetime demand not met on all runs and calculate statistics afterwards. 
def calc_pctdeficit_allruns(d_deficit,d_withdrawal):
    return {k: d_deficit[k]/d_withdrawal[k] * 100 for k in d_withdrawal.keys() & d_deficit}

    
# calculate pct water deficit per region (weighted average based on cohort size) based on waterdeficit/withdrawal *100 for every country
# (This is alternative to what is done in the main function: calculating the waterdeficit and withdrawal per region based on countries and only then taking the ratio)
def calc_pctwaterdeficit_perregion_perrun(d_pct_waterdeficit_perrun):

    d_regions = pd.read_pickle(open('./data/pickles/region_info.pkl', 'rb'))
    d_cohort_weights_regions = d_regions['cohort_size']
    
    d_countries = pd.read_pickle(open('./data/pickles/country_info.pkl', 'rb'))
    df_countries = d_countries['info_pop']
    
    d_exposure_perregion_perrun = {}

    for i in d_pct_waterdeficit_perrun.keys():

        d_exposure_perregion = {}


        for region in d_cohort_weights_regions.keys():

            member_countries = get_countries_of_region(region, df_countries)  

            exposure_perrun =  d_pct_waterdeficit_perrun[i]

            # filter cohort weights to only keep countries within mask 
            d_cohort_weights_regions[region] = d_cohort_weights_regions[region].loc[:,d_cohort_weights_regions[region].columns.isin(df_countries.index)]

            # get weighted spatial average for all member countries per region (Luke: but this is not necessarily spatial; d_cohorts_weights_regions has 2020 cohort sizes per country)
            d_exposure_perregion[region] = (exposure_perrun.loc[:,member_countries] * d_cohort_weights_regions[region].values).sum(axis=1) /\
                np.nansum(d_cohort_weights_regions[region].values, axis=1)

        d_exposure_perregion_perrun[i]  = pd.DataFrame(d_exposure_perregion)

    return d_exposure_perregion_perrun

#%% ----------------------------------------------------------------
# convert nmonths and nevents to duration

def calc_duration(
    grid_area,
    d_regions,
    d_isimip_meta, 
    df_birthyears_regions, 
    df_countries, 
    countries_regions, 
    countries_mask, 
    da_population, 
    df_life_expectancy_5,
    d_all_cohorts,
    modes,

):

# duration function

    phase, ages, age_young, age_ref, age_range, year_ref, year_start, birth_years, year_end, year_range, RCP2GMT_maxdiff_threshold, year_start_GMT_ref, year_end_GMT_ref, model_names = init()

    d_duration_perrun_RCP     = {}
    d_duration_perregion_perrun_RCP = {}
    d_duration_percountry = {}
    d_landfrac_peryear_perregion = {}
    da_exposure_cohort = {}


    # unpack region information
    df_birthyears_regions = d_regions['birth_years']
    d_cohort_weights_regions = d_regions['cohort_size']    
    da_cohort_size = xr.DataArray(
        # np.asarray(list(d_cohort_size.values())),
        np.asarray([v for k,v in d_all_cohorts.items() if k in list(df_countries['name'])]),
        coords={
            'country': ('country', list(df_countries['name'])),
            'time': ('time', year_range),
            'ages': ('ages', np.arange(104,-1,-1)),
        },
        dims=[
            'country',
            'time',
            'ages',
        ]
    ) 



    # loop over simulations
    for i in list(d_isimip_meta.keys()): 

        print('simulation '+str(i)+ ' of '+str(len(d_isimip_meta)))

        # load AFA data of that run
        with open('./data/pickles/isimip_AFA_{}_{}_{}.pkl'.format(d_isimip_meta[i]['extreme'],d_isimip_meta[i]['mode'],str(i)), 'rb') as f:
            da_AFA = pd.read_pickle(f)

        # --------------------------------------------------------------------
        # apply dryland mask --- OLD with dryland mask based on ISIMIP

        #if modes=='dryland': 
            # load dryland mask related to simulation
        #     dryland_mask = load_drylandmask(year_start, year_end, d_isimip_meta[i]['model'], d_isimip_meta[i]['gcm'], d_isimip_meta[i]['rcp'])
        #     da_population = da_population.where(dryland_mask).fillna(0)

        # --------------------------------------------------------------------
        # per country 

        # initialise dicts
        d_exposure_peryear_percountry = {}

        # get spatial average
        for j, country in enumerate(df_countries['name']):

            print('processing country '+str(j+1)+' of '+str(len(df_countries)), end='\r')

            # calculate mean per country weighted by population
            ind_country = countries_regions.map_keys(country)

            # calculate per country and birthyear the mean duration of events. 
            da_nevents_country = da_AFA['nevents'].where(countries_mask == ind_country, drop=True)
            da_nmonths_country = da_AFA['nmonths'].where(countries_mask == ind_country, drop=True)


            # initialise birth years 
            list_da_duration_birthyears_percountry = []
            list_da_nmonths_birthyears_percountry = []

            d_duration_perbirthyear = {} 
            for birth_year in df_life_expectancy_5.index:

                df_life_expectancy = df_life_expectancy_5.loc[birth_year,country] 

                # define death year based on life expectancy
                death_year = birth_year + np.floor(df_life_expectancy)

                # integrate exposure over full years lived
                nevents_birthyears_percountry = da_nevents_country.sel(time=slice(birth_year,death_year)).sum('time')
                nmonths_birthyears_percountry = da_nmonths_country.sel(time=slice(birth_year,death_year)).sum('time')
                duration_birthyears_percountry = nmonths_birthyears_percountry / nevents_birthyears_percountry
                # correct for gridcells where there is no event (nevents == 0, which lead to infinite durations 
                duration_birthyears_percountry = duration_birthyears_percountry.where(duration_birthyears_percountry!=np.inf,  np.nan)
                list_da_duration_birthyears_percountry.append(duration_birthyears_percountry)

                # historical + RCP simulations
                da_duration = calc_weighted_fldmean( 
                    duration_birthyears_percountry,
                    da_population.where(duration_birthyears_percountry, drop=True).sel(time=2020).squeeze().drop('time'), 
                    countries_mask, 
                    ind_country, 
                    flag_region= False)

                d_duration_perbirthyear[birth_year] = float(da_duration.values)
            d_duration_percountry[country] = d_duration_perbirthyear

        duration_birthyears_percountry = pd.DataFrame.from_dict(d_duration_percountry,orient='index')

        d_duration_perrun_RCP[i] = duration_birthyears_percountry.transpose()



        # --------------------------------------------------------------------
        # per region
        #  

        print('')

        # initialise dictionaries
        d_landfrac_peryear_perregion[i] = {}
        d_duration_perregion_RCP = {}

        # loop over regions
        for k, region in enumerate(df_birthyears_regions.columns): 

            print('processing region '+str(k+1)+' of '+str(len(df_birthyears_regions.columns)), end='\r')

            # Get list of member countries from region - with seperate treatment for world (luke: now inside get_countries_of_regions func)
            member_countries = get_countries_of_region(region, df_countries)

            # get spatial average of landfraction: historical + RCP simulations
            ind_countries = countries_regions.map_keys(member_countries)

            #print('calculating cohort weights')
            # filter cohort weights to only keep countries within mask 
            d_cohort_weights_regions[region] = d_cohort_weights_regions[region].loc[:,d_cohort_weights_regions[region].columns.isin(df_countries.index)]

            # get weighted spatial average for all member countries per region (Luke: but this is not necessarily spatial; d_cohorts_weights_regions has 2020 cohort sizes per country)
            d_duration_perregion_RCP[region] = (d_duration_perrun_RCP[i].loc[:,member_countries] * d_cohort_weights_regions[region].values).sum(axis=1) /\
                np.nansum(d_cohort_weights_regions[region].values, axis=1)

        print('')

        # save exposures for every run
        d_duration_perregion_perrun_RCP[i]  = pd.DataFrame(d_duration_perregion_RCP)



    # --------------------------------------------------------------------
    # save workspave in pickles
    #  

    # save pickles
    print()
    print('Saving processed durations')

    # pack exposure information
    d_exposure = {
        'exposure_perrun_RCP' : d_duration_perrun_RCP, 
        'exposure_perregion_perrun_RCP' : d_duration_perregion_perrun_RCP, 
        'landfrac_peryear_perregion' : d_landfrac_peryear_perregion,
        'exposure_per_cohort': da_exposure_cohort,
        'exposure_perrun_15' :d_duration_perrun_RCP,
        'exposure_perrun_20' : d_duration_perrun_RCP,
        'exposure_perrun_NDC' : d_duration_perrun_RCP,    
    }


    with open('./data/pickles/exposure_{}_{}.pkl'.format(d_isimip_meta[i]['extreme'],d_isimip_meta[i]['mode']), 'wb') as f:
        pk.dump(d_exposure,f)


    return d_duration_perrun_RCP, d_duration_perregion_perrun_RCP
