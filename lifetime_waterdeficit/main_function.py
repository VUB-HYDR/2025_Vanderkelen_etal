#!/usr/bin/env python3
# ---------------------------------------------------------------
# Main FUNCTION to postprocess and visualise lifetime exposure data
# Exactly the same as the main script, but written as a function. 

# Python translation of the MATLAB scripts of Thiery et al. (2021)
# https://github.com/VUB-HYDR/2021_Thiery_etal_Science
# ----------------------------------------------------------------


import os
import xarray as xr
import pickle as pk
import time
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import mapclassify as mc
from copy import deepcopy as cp

from settings import *
from load_manip import *
from exposure import * 
from utils import *


def do_lifetime_analysis(flags):

    # 1. Intialise
    
    scriptsdir = os.getcwd()


    # set global variables
    phase, ages, age_young, age_ref, age_range, year_ref, year_start, birth_years, year_end, year_range, RCP2GMT_maxdiff_threshold, year_start_GMT_ref, year_end_GMT_ref, model_names = init()

    # set extremes based on flag (this needs to happen here as it uses the flags dict defined above)
    extremes = set_extremes(flags)

    # set modes based on flag 
    modes = set_modes(flags)
    
    
    # 2. Load and manipulate data

    # Load global mean temperature projections
    global df_GMT_15, df_GMT_20, df_GMT_NDC

    df_GMT_15, df_GMT_20, df_GMT_NDC = load_GMT(
        year_start,
        year_end,
    ) 


    #  Load and manipulate life expectancy, cohort and mortality data

    if flags['mask']: # load data and do calculations

        print('Processing country info')

        # load worldbank and unwpp data
        meta, worldbank, unwpp = load_worldbank_unwpp_data()

        # unpack values
        df_countries, df_regions = meta
        df_worldbank_country, df_worldbank_region = worldbank
        df_unwpp_country, df_unwpp_region = unwpp

        # manipulate worldbank and unwpp data to get birth year and life expectancy values
        df_birthyears, df_life_expectancy_5 = get_life_expectancies(
            df_worldbank_country, 
            df_unwpp_country,
        )

        # load population size per age cohort data
        wcde = load_wcde_data() 

        # interpolate population size per age cohort data to our ages (0-60)
        d_cohort_size = get_cohortsize_countries(
            wcde, 
            df_countries, 
            df_GMT_15,
        )

        d_all_cohorts = get_all_cohorts(
            wcde, 
            df_countries, 
            df_GMT_15,
        )

        # do the same for the regions; get life expectancy, birth years and cohort weights per region, as well as countries per region
        d_region_countries, df_birthyears_regions, df_life_expectancy_5_regions, d_cohort_weights_regions = get_regions_data(
            df_countries, 
            df_regions, 
            df_worldbank_region, 
            df_unwpp_region, 
            d_cohort_size,
        )
        
        # apply option of constant life expectancy
        if flags['expectancy']:
            df_life_expectancy_5.loc[year_start+1:year_end] = df_life_expectancy_5.loc[year_start].values 
            df_life_expectancy_5_regions.loc[year_start+1:year_end] = df_life_expectancy_5_regions.loc[year_start].values            
            
        # --------------------------------------------------------------------
        # Load population and country masks, and mask population per country

        # Load SSP population totals 
        da_population = load_population(
            year_start,
            year_end,
        )

        gdf_country_borders = gpd.read_file('./data/natural_earth/Cultural_10m/Countries/ne_10m_admin_0_countries.shp'); 

        # mask population totals per country  and save country regions object and countries mask
        df_countries, countries_regions, countries_mask, gdf_country_borders = get_mask_population(
            da_population, 
            gdf_country_borders, 
            df_countries,
        ) 

        # pack country information
        d_countries = {
            'info_pop' : df_countries, 
            'borders' : gdf_country_borders,
            'population_map' : da_population,
            'birth_years' : df_birthyears,
            'life_expectancy_5': df_life_expectancy_5, 
            'cohort_size' : d_cohort_size, 
            'all_cohorts' : d_all_cohorts,
            'mask' : (countries_regions,countries_mask),
        }

        # pack region information
        d_regions = {
            'birth_years' : df_birthyears_regions,
            'life_expectancy_5': df_life_expectancy_5_regions, 
            'cohort_size' : d_cohort_weights_regions,
        }

        # save metadata dictionary as a pickle
        print('Saving country and region data')

        if not os.path.isdir('./data/pickles'):
            os.mkdir('./data/pickles')
        with open('./data/pickles/country_info.pkl', 'wb') as f: # note; 'with' handles file stream closing
            pk.dump(d_countries,f)
        with open('./data/pickles/region_info.pkl', 'wb') as f:
            pk.dump(d_regions,f)

    else: # load processed country data

        print('Loading processed country and region data')

        # load country pickle
        d_countries = pd.read_pickle(open('./data/pickles/country_info.pkl', 'rb'))

        # unpack country information
        df_countries = d_countries['info_pop']
        gdf_country_borders = d_countries['borders']
        da_population = d_countries['population_map']
        df_birthyears = d_countries['birth_years']
        df_life_expectancy_5 = d_countries['life_expectancy_5']
        d_cohort_size = d_countries['cohort_size']
        d_all_cohorts = d_countries['all_cohorts']
        countries_regions, countries_mask = d_countries['mask']

        # load regions pickle
        d_regions = pd.read_pickle(open('./data/pickles/region_info.pkl', 'rb'))

        # unpack region information
        df_birthyears_regions = d_regions['birth_years']
        df_life_expectancy_5_regions = d_regions['life_expectancy_5']
        d_cohort_weights_regions = d_regions['cohort_size']
        
        # constant life expectancy
        if flags['expectancy']:
            df_life_expectancy_5.loc[year_start+1:year_end] = df_life_expectancy_5.loc[year_start].values 
            df_life_expectancy_5_regions.loc[year_start+1:year_end] = df_life_expectancy_5_regions.loc[year_start].values            
            


    # 4.3  Load ISIMIP model data
    global grid_area
    grid_area = xr.open_dataarray('./data/isimip/clm45_area.nc4')


    d_isimip_meta, d_pic_meta = load_isimip(
    phase,
    flags['runs'], 
    flags['exposure_pic'],
    extremes, 
    modes,
    model_names,
    df_GMT_15,
    df_GMT_20,
    df_GMT_NDC,
    year_end_GMT_ref, 
    year_start_GMT_ref, 
    year_end, 
    year_start, 
    RCP2GMT_maxdiff_threshold)
    


    # 5. Compute exposure per lifetime

    if flags['exposure']: 

        start_time = time.time()

        if modes[0] != 'duration': 

            # calculate exposure per country and per region and save data (takes 23 mins)
            d_exposure_perrun_RCP, d_exposure_perregion_perrun_RCP, d_exposure_perrun_15, d_exposure_perrun_20, d_exposure_perrun_NDC, da_exposure_cohort = calc_exposure(
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
            )

            print("--- {} minutes ---".format(
                np.floor((time.time() - start_time) / 60),
                )
                  )

            # pack country information
            d_exposure = {}
            d_exposure['exposure_perrun_RCP'] =  d_exposure_perrun_RCP
            d_exposure['exposure_perrun_15'] = d_exposure_perrun_15 
            d_exposure['exposure_perrun_20'] = d_exposure_perrun_20
            d_exposure['exposure_perrun_NDC'] = d_exposure_perrun_NDC
            d_exposure['exposure_per_cohort'] = da_exposure_cohort

        elif modes[0] == 'duration':


            d_duration_perrun_RCP, d_duration_perregion_perrun_RCP = calc_duration(
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
                modes)        

            print("--- {} minutes ---".format(
                np.floor((time.time() - start_time) / 60),
                )
                  )

            # pack country information
            d_exposure = {}
            d_exposure['exposure_perrun_RCP'] =  d_duration_perrun_RCP

            # these are placeholders, as not used anymore
            d_exposure['exposure_perrun_15'] = d_duration_perrun_RCP 
            d_exposure['exposure_perrun_20'] = d_duration_perrun_RCP
            d_exposure['exposure_perrun_NDC'] = d_duration_perrun_RCP
            d_exposure['exposure_per_cohort'] = d_duration_perrun_RCP


    else: # load processed country data

        print('Loading processed exposures')

        # load country pickle
        with open('./data/pickles/exposure_{}_{}.pkl'.format(d_isimip_meta[list(d_isimip_meta.keys())[0]]['extreme'],d_isimip_meta[list(d_isimip_meta.keys())[0]]['mode']), 'rb') as f:
            d_exposure = pd.read_pickle(f)

        # unpack country information
        d_exposure_perrun_RCP = d_exposure['exposure_perrun_RCP']
        d_exposure_perrun_15 = d_exposure['exposure_perrun_15']
        d_exposure_perrun_20 = d_exposure['exposure_perrun_20']
        d_exposure_perrun_NDC = d_exposure['exposure_perrun_NDC']
        da_exposure_cohort = d_exposure['exposure_per_cohort']

        # unpack region information
        d_exposure_perregion_perrun_RCP = d_exposure['exposure_perregion_perrun_RCP']
        d_landfrac_peryear_perregion = d_exposure['landfrac_peryear_perregion']





    # 5.3 Compile hist+RCP and pic for EMF  analysis
    # call function to compute mmm, std, qntl for exposure (also 99.99 % of pic as "ext")
    ds_exposure_RCP = calc_exposure_mmm_xr(
        d_exposure_perrun_RCP,
        'country',
        'RCP',
    )
    
        
    ds_exposure_15 = calc_exposure_mmm_xr(
        d_exposure_perrun_15,
        'country',
        '15',
    )
    ds_exposure_20 = calc_exposure_mmm_xr(
        d_exposure_perrun_20,
        'country',
        '20',
    )
    ds_exposure_NDC = calc_exposure_mmm_xr(
        d_exposure_perrun_NDC,
        'country',
        'NDC',
    )
    ds_exposure_perregion = calc_exposure_mmm_xr(
        d_exposure_perregion_perrun_RCP,
        'region',
        'RCP',
    )

    if flags['exposure_pic']:

        ds_exposure_pic = calc_exposure_mmm_pic_xr(
            d_exposure_perrun_pic,
            'country',
            'pic',
        )
        ds_exposure_pic_perregion = calc_exposure_mmm_pic_xr(
            d_exposure_perregion_perrun_pic,
            'region',
            'pic',
        )

    # pool all datasets for different trajectories
    ds_exposure = xr.merge([
        ds_exposure_RCP,
        ds_exposure_15,
        ds_exposure_20,
        ds_exposure_NDC,
    ])


    # Do the same for landfrac
    #ds_landfrac_stats_RCP = calc_landfrac_stats(d_landfrac_peryear_perregion, 'RCP')

    #ds_landfrac_stats_15 = calc_landfrac_stats(d_landfrac_peryear_perregion, '15')
    #ds_landfrac_stats_20 = calc_landfrac_stats(d_landfrac_peryear_perregion, '20')
    #ds_landfrac_stats_NDC = calc_landfrac_stats(d_landfrac_peryear_perregion, 'NDC')

    # pool all datasets for different trajectories
    #ds_landfrac = xr.merge([
    #    ds_exposure_RCP,
    #    ds_landfrac_stats_15,
    #    ds_landfrac_stats_20,
    #    ds_landfrac_stats_NDC,
    #])

    # Info on return variables: 
    #  ds_exposure and ds_exposure_perregion are datasets with statistics on exposure (mean, quantiles, std)
    #  d_exposure is the dictionary with the exposures of all runs together, to calculate own stats. 

    return ds_exposure, ds_exposure_perregion, d_exposure


# Main function to calculate water deficit (%) based on RCP scenarios
def calc_lifetime_waterdeficit(flags): 
    
    # load lifetime waterdeficit
    flags['extr'] = 'waterdeficit'   # 0: all
    ds_waterdeficit, ds_waterdeficit_perregion, d_waterdeficit_allruns = do_lifetime_analysis(flags)
    
    
    # load lifetime waterwithdrawal
    flags['extr'] = 'withdrawal'   # 0: all

    ds_withdrawal, ds_withdrawal_perregion, d_withdrawal_allruns = do_lifetime_analysis(flags)
    
    
    ## find indices of RCP scenarios per simulation

    # unpack info from isimip simulations - to know RCP scenario
    with open('./data/pickles/isimip_metadata_{}_{}.pkl'.format(flags['extr'],flags['mode']), 'rb') as f:
        d_isimip_meta = pd.read_pickle(f)

    ind_rcp26 = []
    ind_rcp60 = []
    for i in d_isimip_meta.keys(): 
        if d_isimip_meta[i]['rcp'] == 'rcp26':
            ind_rcp26.append(i)
        elif d_isimip_meta[i]['rcp'] == 'rcp60':
            ind_rcp60.append(i)

    # filter simulations based on RCP scenario
    d_waterdeficit_rcp26runs = {i: d_waterdeficit_allruns['exposure_perrun_RCP'][i] for i in ind_rcp26}
    d_waterdeficit_rcp60runs = {i: d_waterdeficit_allruns['exposure_perrun_RCP'][i] for i in ind_rcp60}
    d_withdrawal_rcp26runs = {i: d_withdrawal_allruns['exposure_perrun_RCP'][i] for i in ind_rcp26}
    d_withdrawal_rcp60runs = {i: d_withdrawal_allruns['exposure_perrun_RCP'][i] for i in ind_rcp60}

    
    # per country
    
    d_pct_waterdeficit_perrun_RCP26 = calc_pctdeficit_allruns(d_waterdeficit_rcp26runs,d_withdrawal_rcp26runs)
    d_pct_waterdeficit_perrun_RCP60 = calc_pctdeficit_allruns(d_waterdeficit_rcp60runs,d_withdrawal_rcp60runs)
    d_pct_waterdeficit_perrun_RCP   = calc_pctdeficit_allruns(d_waterdeficit_allruns['exposure_perrun_RCP'],d_withdrawal_allruns['exposure_perrun_RCP'])

    # average of both RCPs

    ds_pct_waterdeficit_RCP26 = calc_exposure_mmm_xr(d_pct_waterdeficit_perrun_RCP26, 'country', 'RCP26' )
    ds_pct_waterdeficit_RCP60 = calc_exposure_mmm_xr(d_pct_waterdeficit_perrun_RCP60, 'country', 'RCP60' )
    ds_pct_waterdeficit_RCP   = calc_exposure_mmm_xr(d_pct_waterdeficit_perrun_RCP, 'country', 'RCP' )

    ds_pct_waterdeficit_mitigation = calc_exposure_mmm_mitigation_xr(d_pct_waterdeficit_perrun_RCP26, d_pct_waterdeficit_perrun_RCP60, 'country','' )

    # pool all datasets for different trajectories
    ds_pct_waterdeficit = xr.merge([  ds_pct_waterdeficit_RCP26,  ds_pct_waterdeficit_RCP60, ds_pct_waterdeficit_RCP, ds_pct_waterdeficit_mitigation])

    # per region (calculating based on country withdrawal and deficit vs averaging country % water deficits directly)
    # (This is alternative to what is done in the main function: calculating the waterdeficit and withdrawal per region based on countries and only then taking the ratio)

    # calculate pct water deficit per region (weighted average based on cohort size) based on waterdeficit/withdrawal *100 for every country
    # (This is alternative to what is done in the main function: calculating the waterdeficit and withdrawal per region based on countries and only then taking the ratio)
    d_pct_waterdeficit_perregion_perrun_RCP26 = calc_pctwaterdeficit_perregion_perrun(d_pct_waterdeficit_perrun_RCP26)
    d_pct_waterdeficit_perregion_perrun_RCP60 = calc_pctwaterdeficit_perregion_perrun(d_pct_waterdeficit_perrun_RCP60)
    d_pct_waterdeficit_perregion_perrun_RCP   = calc_pctwaterdeficit_perregion_perrun(d_pct_waterdeficit_perrun_RCP)

    # calculate multi-model statistics
    ds_pct_waterdeficit_perregion_RCP26 = calc_exposure_mmm_xr(d_pct_waterdeficit_perregion_perrun_RCP26, 'region', 'RCP26' )
    ds_pct_waterdeficit_perregion_RCP60 = calc_exposure_mmm_xr(d_pct_waterdeficit_perregion_perrun_RCP60, 'region', 'RCP60' )
    ds_pct_waterdeficit_perregion_RCP   = calc_exposure_mmm_xr(d_pct_waterdeficit_perregion_perrun_RCP  , 'region', 'RCP' )


        
    ds_pct_waterdeficit_perregion = xr.merge([ ds_pct_waterdeficit_perregion_RCP26,  ds_pct_waterdeficit_perregion_RCP60, ds_pct_waterdeficit_perregion_RCP])
                                    
    return ds_pct_waterdeficit, ds_pct_waterdeficit_perregion, d_pct_waterdeficit_perrun_RCP26, d_pct_waterdeficit_perrun_RCP60


# Main function to calculate water deficit based on different GMT levels (old)
def calc_lifetime_waterdeficit_GMT(flags): 
    # load deficit
    flags['extr'] = 'waterdeficit'   # 0: all
    ds_waterdeficit, ds_waterdeficit_perregion, d_waterdeficit_allruns = do_lifetime_analysis(flags)

    # load withdrawal
    flags['extr'] = 'withdrawal'   # 0: all
    ds_withdrawal, ds_withdrawal_perregion, d_withdrawal_allruns = do_lifetime_analysis(flags)


    # per country
    d_pct_waterdeficit_perrun_RCP = calc_pctdeficit_allruns(d_waterdeficit_allruns['exposure_perrun_RCP'],d_withdrawal_allruns['exposure_perrun_RCP'])

    d_pct_waterdeficit_perrun_15  = calc_pctdeficit_allruns(d_waterdeficit_allruns['exposure_perrun_15'], d_withdrawal_allruns['exposure_perrun_15'])
    d_pct_waterdeficit_perrun_20  = calc_pctdeficit_allruns(d_waterdeficit_allruns['exposure_perrun_20'], d_withdrawal_allruns['exposure_perrun_20'])
    d_pct_waterdeficit_perrun_NDC = calc_pctdeficit_allruns(d_waterdeficit_allruns['exposure_perrun_NDC'],d_withdrawal_allruns['exposure_perrun_NDC'])

    ds_pct_waterdeficit_RCP = calc_exposure_mmm_xr(d_pct_waterdeficit_perrun_RCP, 'country', 'RCP' )
    ds_pct_waterdeficit_20 = calc_exposure_mmm_xr(d_pct_waterdeficit_perrun_20, 'country', '20' )
    ds_pct_waterdeficit_15 = calc_exposure_mmm_xr(d_pct_waterdeficit_perrun_15, 'country', '15' )
    ds_pct_waterdeficit_NDC = calc_exposure_mmm_xr(d_pct_waterdeficit_perrun_NDC, 'country', 'NDC' )

    # pool all datasets for different trajectories
    ds_pct_waterdeficit = xr.merge([ ds_pct_waterdeficit_RCP, ds_pct_waterdeficit_15, ds_pct_waterdeficit_20, ds_pct_waterdeficit_NDC])



    # per region (calculating based on country withdrawal and deficit vs averaging country % water deficits directly)
    # (This is alternative to what is done in the main function: calculating the waterdeficit and withdrawal per region based on countries and only then taking the ratio)

    # calculate pct water deficit per region (weighted average based on cohort size) based on waterdeficit/withdrawal *100 for every country
    # (This is alternative to what is done in the main function: calculating the waterdeficit and withdrawal per region based on countries and only then taking the ratio)
    d_pct_waterdeficit_perregion_perrun_RCP = calc_pctwaterdeficit_perregion_perrun(d_pct_waterdeficit_perrun_RCP)
    d_pct_waterdeficit_perregion_perrun_15  = calc_pctwaterdeficit_perregion_perrun(d_pct_waterdeficit_perrun_15)
    d_pct_waterdeficit_perregion_perrun_20  = calc_pctwaterdeficit_perregion_perrun(d_pct_waterdeficit_perrun_20)
    d_pct_waterdeficit_perregion_perrun_NDC = calc_pctwaterdeficit_perregion_perrun(d_pct_waterdeficit_perrun_NDC)


    # calculate multi-model statistics
    ds_pct_waterdeficit_perregion_RCP = calc_exposure_mmm_xr(d_pct_waterdeficit_perregion_perrun_RCP, 'region', 'RCP' )
    ds_pct_waterdeficit_perregion_15  = calc_exposure_mmm_xr(d_pct_waterdeficit_perregion_perrun_15, 'region', '15' )
    ds_pct_waterdeficit_perregion_20  = calc_exposure_mmm_xr(d_pct_waterdeficit_perregion_perrun_20, 'region', '20' )
    ds_pct_waterdeficit_perregion_NDC = calc_exposure_mmm_xr(d_pct_waterdeficit_perregion_perrun_NDC, 'region', 'NDC' )

    ds_pct_waterdeficit_perregion = xr.merge([ ds_pct_waterdeficit_perregion_RCP, ds_pct_waterdeficit_perregion_15, ds_pct_waterdeficit_perregion_20, ds_pct_waterdeficit_perregion_NDC])

    return ds_pct_waterdeficit, ds_pct_waterdeficit_perregion


# Main function to calculate absolute water deficit (m3/lifetime) based on RCP scenarios
def calc_lifetime_waterdeficit_absolute(flags): 
    
    # load lifetime waterdeficit
    flags['extr'] = 'waterdeficit'   # 0: all
    ds_waterdeficit, ds_waterdeficit_perregion, d_waterdeficit_allruns = do_lifetime_analysis(flags) 
    
    ## find indices of RCP scenarios per simulation

    # unpack info from isimip simulations - to know RCP scenario
    with open('./data/pickles/isimip_metadata_{}_{}.pkl'.format(flags['extr'],'exposure'), 'rb') as f:
        d_isimip_meta = pd.read_pickle(f)

    ind_rcp26 = []
    ind_rcp60 = []
    for i in d_isimip_meta.keys(): 
        if d_isimip_meta[i]['rcp'] == 'rcp26':
            ind_rcp26.append(i)
        elif d_isimip_meta[i]['rcp'] == 'rcp60':
            ind_rcp60.append(i)

    # filter simulations based on RCP scenario
    d_waterdeficit_rcp26runs = {i: d_waterdeficit_allruns['exposure_perrun_RCP'][i] for i in ind_rcp26}
    d_waterdeficit_rcp60runs = {i: d_waterdeficit_allruns['exposure_perrun_RCP'][i] for i in ind_rcp60}

    # average of both RCPs

    ds_waterdeficit_RCP26 = calc_exposure_mmm_xr(d_waterdeficit_rcp26runs, 'country', 'RCP26' )
    ds_waterdeficit_RCP60 = calc_exposure_mmm_xr(d_waterdeficit_rcp60runs, 'country', 'RCP60' )

    # pool all datasets for different trajectories
    ds_waterdeficit = xr.merge([  ds_waterdeficit_RCP26,  ds_waterdeficit_RCP60])
                                    
    return ds_waterdeficit


# Main function to calculate mean duration of deficit 
def calc_lifetime_deficit_duration(flags): 
    
    # load lifetime waterdeficit
    ds_duration, ds_duration_perregion, d_duration_allruns = do_lifetime_analysis(flags)
    
    ## find indices of RCP scenarios per simulation

    # unpack info from isimip simulations - to know RCP scenario
    with open('./data/pickles/isimip_metadata_{}_{}.pkl'.format(flags['extr'],flags['mode']), 'rb') as f:
        d_isimip_meta = pd.read_pickle(f)

    ind_rcp26 = []
    ind_rcp60 = []
    for i in d_isimip_meta.keys(): 
        if d_isimip_meta[i]['rcp'] == 'rcp26':
            ind_rcp26.append(i)
        elif d_isimip_meta[i]['rcp'] == 'rcp60':
            ind_rcp60.append(i)

    # filter simulations based on RCP scenario
    d_duration_rcp26runs = {i: d_duration_allruns['exposure_perrun_RCP'][i] for i in ind_rcp26}
    d_duration_rcp60runs = {i: d_duration_allruns['exposure_perrun_RCP'][i] for i in ind_rcp60}


    # average of both RCPs

    ds_duration_RCP26 = calc_exposure_mmm_xr(d_duration_rcp26runs, 'country', 'RCP26' )
    ds_duration_RCP60 = calc_exposure_mmm_xr(d_duration_rcp60runs, 'country', 'RCP60' )

    # pool all datasets for different trajectories
    ds_duration_both = xr.merge([  ds_duration_RCP26,  ds_duration_RCP60])

    # per region (calculating based on country withdrawal and deficit vs averaging country % water deficits directly)
    # (This is alternative to what is done in the main function: calculating the waterdeficit and withdrawal per region based on countries and only then taking the ratio)

    # calculate pct water deficit per region (weighted average based on cohort size) based on waterdeficit/withdrawal *100 for every country
    # (This is alternative to what is done in the main function: calculating the waterdeficit and withdrawal per region based on countries and only then taking the ratio)
    d_duration_perregion_perrun_RCP26 = calc_pctwaterdeficit_perregion_perrun(d_duration_rcp26runs)
    d_duration_perregion_perrun_RCP60 = calc_pctwaterdeficit_perregion_perrun(d_duration_rcp60runs)

    # calculate multi-model statistics
    ds_duration_perregion_RCP26 = calc_exposure_mmm_xr(d_duration_perregion_perrun_RCP26, 'region', 'RCP26' )
    ds_duration_perregion_RCP60 = calc_exposure_mmm_xr(d_duration_perregion_perrun_RCP60, 'region', 'RCP60' )

    ds_duration_perregion_both = xr.merge([ ds_duration_perregion_RCP26,  ds_duration_perregion_RCP60])
                                    
    return ds_duration_both, ds_duration_perregion_both, d_duration_rcp26runs, d_duration_rcp60runs
