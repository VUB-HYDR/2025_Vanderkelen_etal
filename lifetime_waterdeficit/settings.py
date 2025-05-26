# ----------------------------------------------------------------
# Settings
# These are global variables to be used throughout the whole project
# ----------------------------------------------------------------

import numpy as np





#%% ----------------------------------------------------------------
def init(): 

    # Initialise based on ISIMIP phase
    global phase
    
    phase = '2b'
    
    # initialise age and associated time period of interest
    global ages, age_young, age_ref, age_range, year_ref, year_start, birth_years, year_end, year_range
    ages        = np.arange(60,-1,-1)
    age_young   = 0
    age_ref     = np.nanmax(ages)
    age_range   = np.arange(0,105)
    year_ref    = 2020
    year_start  = year_ref - age_ref
    birth_years = np.arange(year_start,year_ref+1)     
    year_end    = 2113                            # based on maximum life expectancy reported in UN WPP
    year_range  = np.arange(year_start,year_end+1)


    # initialise age groups
    # (https://www.carbonbrief.org/analysis-why-children-must-emit-eight-times-less-co2-than-their-grandparents)
    # (https://www.pewresearch.org/fact-tank/2019/01/17/where-millennials-end-and-generation-z-begins/)
    global agegroups
    agegroups = {
        'Boomers' : (1950, 1965),
        'Gen X' : (1965, 1981),
        'Millenials' : (1981, 1997),
        'Gen Z' : (1997, 2020)
    }


    # initialise reference period for computing GMT anomalies
    global year_start_GMT_ref, year_end_GMT_ref
    year_start_GMT_ref = 1850
    year_end_GMT_ref = 1900

    global extremes_legend_dict
    extremes_legend_dict = {
        'waterscarcity' : 'Water scarcity index',
        'falkenmark': 'Falkenmark index',
        'waterdeficit': 'Water deficit',
        'withdrawal': 'Water withdrawal',
        'all' : 'All'               
    }
    
    
    # initialise model names
    global model_names

    
    if phase == '3b':
        model_names = {
          'waterscarcity'  : ['CWatM', 'WaterGAP2-2e','H08' ], # ,  'WaterGAP2-2e'],# 'H08']
          'falkenmark' : ['CWatM', 'WaterGAP2-2e', 'H08'],#
          'waterdeficit' : ['CWatM', 'WaterGAP2-2e', 'H08'] ,#
          'withdrawal'  : ['CWatM', 'WaterGAP2-2e', 'H08'],
        }
    elif phase == '2b':
        model_names = {
          'waterscarcity'  : ['LPJmL','MATSIRO', 'H08' , 'CWatM' ], # ,  'WaterGAP2-2e'],# 'H08']
          'falkenmark' : ['LPJmL','MATSIRO', 'H08' , 'CWatM'],#
          'waterdeficit' : ['LPJmL','MATSIRO', 'H08' , 'CWatM'] ,#
          'withdrawal'  : ['LPJmL','MATSIRO', 'H08' , 'CWatM'], 
           'waterdeficitcstpop'  :  ['LPJmL','MATSIRO', 'H08' , 'CWatM'],
           'waterdeficitduration'  :  ['LPJmL','MATSIRO', 'H08' , 'CWatM'],
           'waterdeficitintensity'  :  ['LPJmL','MATSIRO', 'H08' , 'CWatM'],
           'withdrawalcstpop'  :  ['LPJmL','MATSIRO', 'H08' , 'CWatM'],
           'indww'  :  [ 'H08'],
           'domww'  :  ['H08'],
           'irrww'  :  ['LPJmL','MATSIRO', 'H08' , 'CWatM'],
           'irrwwwaterdeficit'  :  ['LPJmL','MATSIRO', 'H08' , 'CWatM'],
           'indwwwaterdeficit'  :  ['LPJmL','MATSIRO', 'H08' , 'CWatM'],
           'domwwwaterdeficit'  :  ['LPJmL','MATSIRO', 'H08' , 'CWatM'],
           'h08waterdeficitcstpop':  [ 'H08' ],
           'h08withdrawalcstpop':  [ 'H08' ],
           'test':  [ 'H08' ],
       }

    # Set threshold maximum T difference between RCP and GMT trajectories
    # i.e. any run with T difference exceeding this threshold is excluded
    # year-to-year jumps in GMT larger than 0.1, so using a 0.1 maxdiff threshold erronously removes runs
    # used to be 0.5, but then rcp2.6 is used for high-warming levels
    # Anything between 0.1 and 0.2 removes RCP2.6 in NDC scenarios (see histograms of maxdiff_NDC)
    # take 0.2 to have more data in BE scenarios and hence smooth EMF curves in BE plot
    global RCP2GMT_maxdiff_threshold
    RCP2GMT_maxdiff_threshold = 0.2 # [K]


    # set kernel x-values
    global kernel_x
    kernel_x = np.arange(1,50.5,0.5)
    
    return phase, ages, age_young, age_ref, age_range, year_ref, year_start, birth_years, year_end, year_range, RCP2GMT_maxdiff_threshold, year_start_GMT_ref, year_end_GMT_ref, model_names


#%% ----------------------------------------------------------------
# set extremes based on flag (this needs to happen here as it uses the flags dict defined above)
def set_extremes(flags):

    global extremes
    
    if not flags['extr'] == 'all': # processing for single extreme
        extremes = [flags['extr']]
    else: 
        extremes = [
            'waterscarcity',
            'falkenmark',
        ]
    return extremes


#%% ----------------------------------------------------------------
# set modes (exposure, intensity, duration) based on flag (this needs to happen here as it uses the flags dict defined above)
def set_modes(flags):

    global modes
    
    if not flags['mode'] == 'all': # processing for single extreme
        modes = [flags['mode']]
    else: 
        extremes = [
            'exposure',
            'intensity',
            'duration'
        ]
    return modes

