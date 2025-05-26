# Python script with global variables to be used across functions and scripts

import os
import pandas as pd

# Initialise based on ISIMIP phase
# find way to communicate this from main script? 
phase = '2b'

# for isimip 2b
soc_scenario= '2005soc' # 2005soc alternative
soc_scenario= 'histsoc' # 2005soc alternative

# for isimip 3b - tried with 2015soc
histsoc_scenario = 'histsoc' #'2015soc' # default is histsoc



if phase == '3b':

    # 1. PATHS

    # on hydra: 
    try: 

        if os.getlogin() =='vsc10055': 
            datadir = '/data/brussel/vo/000/bvo00012/data/dataset/ISIMIP/ISIMIP3b/'

        # on laptop: 
        if os.getlogin() == 'ivand': 
            datadir = './data/ISIMIP3b/'

    # on a compute node this doesn't work always, just hard code path in that case. 
    except: 
            datadir = '/data/brussel/vo/000/bvo00012/data/dataset/ISIMIP/ISIMIP3b/'


    # scripts directory
    scriptsdir = os.getcwd()

    # directory to store output
    outdir = scriptsdir+'/data/'

    # directory with river network info
    routingdir = datadir + 'InputData/geo_conditions/river_routing/'

    # directory where lifetime framework reads from
    lifetimedir = '/scratch/brussel/vo/000/bvo00012/vsc10055/waterscarcity/lifetime_exposure_isimip/data/isimip/'

    # directory where lifetime framework reads from
    lifetimerootdir = '/scratch/brussel/vo/000/bvo00012/vsc10055/waterscarcity/lifetime_exposure_isimip/'
    
    # 2. MODEL SETTINGS

    # define scenarios
    scenarios = ['historical', 'ssp126', 'ssp370', 'ssp585']

    # dictionary defining in which scenario folder scenario is located
    d_scenario_folder = {'historical' : 'historical',
                        'ssp126'     : 'future'    , 
                        'ssp370'     : 'future'    ,
                        'ssp585'     : 'future'    } 


    # dictionary defining which socio-economic scenario to pick based on scenario
    d_soc_scenario = {'historical' : 'histsoc',
                    'ssp126'     : '2015soc', 
                    'ssp370'     : '2015soc',
                    'ssp585'     : '2015soc'} 
    

    # dictionary defining period based on scenario
    d_senario_period = {'historical' : '1850_2014',
                    'ssp126'     : '2015_2100', 
                    'ssp370'     : '2015_2100',
                    'ssp585'     : '2015_2100'} 

    # initialise reference period for computing GMT anomalies
    year_start_GMT_ref = 1850
    year_end_GMT_ref = 1900

    timestep = 'monthly'


    # variable names for plotting purposes
    variable_name = {'waterscarcity':'Water Scarcity Index', 
                     'falkenmark': 'Falkenmark Index'}
    
    
    # manually define the time periods for isimip2b files, as these times cannot be read in
    d_time_ncfile = {'historical': pd.date_range(start='1850-01-01', end='2015-01-01', freq='M'),
                     'rcp26': pd.date_range(start='2015-01-01', end='2101-01-01', freq='M'),
                     'rcp60':pd.date_range(start='2015-01-01', end='2101-01-01', freq='M') }
    
    
    # test with 2015 soc scenario for the historical period: overwrite settings made above
    if histsoc_scenario == '2015soc':
        d_soc_scenario = {'historical' : '2015soc',
                          'ssp126'     : '2015soc', 
                          'ssp370'     : '2015soc',
                          'ssp585'     : '2015soc'} 
        
        lifetimedir = '/scratch/brussel/vo/000/bvo00012/vsc10055/waterscarcity/lifetime_exposure_isimip/data/isimip/2015soc/'

        
    
#ISIMIP2b    
elif phase == '2b': 

    # 1. PATHS

    # on hydra: 
    try: 

        if os.getlogin() =='vsc10055': 
            datadir = '/data/brussel/vo/000/bvo00012/data/dataset/ISIMIP/ISIMIP2b/'

        # on laptop: 
        if os.getlogin() == 'ivand': 
            datadir = './data/ISIMIP2b/'

    # on a compute node this doesn't work always, just hard code path in that case. 
    except: 
        datadir = '/data/brussel/vo/000/bvo00012/data/dataset/ISIMIP/ISIMIP2b/'


    # scripts directory
    scriptsdir = os.getcwd()

    # directory to store output
    outdir = scriptsdir+'/data/2b/'

    # directory with river network info
    routingdir = '/data/brussel/vo/000/bvo00012/data/dataset/ISIMIP/ISIMIP3b/InputData/geo_conditions/river_routing/'
    
    # directory where lifetime framework reads from
    lifetimedir = '/scratch/brussel/vo/000/bvo00012/vsc10055/waterscarcity/lifetime_exposure_isimip/data/isimip/2b/'

    # directory where lifetime framework reads from
    lifetimerootdir = '/scratch/brussel/vo/000/bvo00012/vsc10055/waterscarcity/lifetime_exposure_isimip/'
    
    # 2. MODEL SETTINGS

    # define scenarios
    scenarios = ['historical', 'rcp26', 'rcp60']

    # dictionary defining in which scenario folder scenario is located
    d_scenario_folder = {'historical' : 'historical',
                        'rcp26'     : 'future'    , 
                        'rcp60'     : 'future'    } 

    if soc_scenario == '2005soc': 
        
        # dictionary defining which socio-economic scenario to pick based on scenario
        d_soc_scenario = {'historical' : '2005soc',
                        'rcp26'     : '2005soc', 
                        'rcp60'     : '2005soc'} 

        lifetimedir = '/scratch/brussel/vo/000/bvo00012/vsc10055/waterscarcity/lifetime_exposure_isimip/data/isimip/2b/2005soc/'
    else: 
        # dictionary defining which socio-economic scenario to pick based on scenario
        d_soc_scenario = {'historical' : 'histsoc',
                        'rcp26'     : 'rcp26soc', 
                        'rcp60'     : 'rcp60soc'} 
        
        
    # dictionary defining period based on scenario
    d_senario_period = {'historical' : '1861_2005',
                    'rcp26'     : '2006_2099', 
                    'rcp60'     : '2006_2099',
                    'ssp585'     : '2006_2099'} 

    # initialise reference period for computing GMT anomalies
    year_start_GMT_ref = 1850
    year_end_GMT_ref = 1900

    timestep = 'monthly'


    # variable names for plotting purposes
    variable_name = {'waterscarcity':'Water Scarcity Index', 
                     'falkenmark': 'Falkenmark Index'}    
    
    
    # manually define the time periods for isimip2b files, as these times cannot be read in
    d_time_ncfile = {'historical': pd.date_range(start='1861-01-01', end='2006-01-01', freq='M'),
                     'rcp26': pd.date_range(start='2006-01-01', end='2100-01-01', freq='M'),
                     'rcp60':pd.date_range(start='2006-01-01', end='2100-01-01', freq='M') }
    # manually define the time periods for isimip2b files, as these times cannot be read in
    d_time_ncfile_annual = {'historical': pd.date_range(start='1861-01-01', end='2006-01-01', freq='Y'),
                     'rcp26': pd.date_range(start='2006-01-01', end='2100-01-01', freq='Y'),
                     'rcp60':pd.date_range(start='2006-01-01', end='2100-01-01', freq='Y') }
    sec_per_year = 365*60*60