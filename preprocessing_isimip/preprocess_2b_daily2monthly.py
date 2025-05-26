# Calculate monthly means from daily for H08 and LPJmL

import os
from init import *
from functions import *

# settings

variables = ['dis', 'qtot']
models = ['H08']#, 'LPJmL',]
forcings  = ['gfdl-esm2m', 'ipsl-cm5a-lr','hadgem2-es', 'miroc5']

d_scenario_years = {'historical' : '1861_2005',
                    'rcp26'     : '2006_2099'    , 
                    'rcp60'     : '2006_2099'    } 


scenarios = ['historical', 'rcp26', 'rcp60']
scenarios = [ 'historical']


for model in models: 
    for variable in variables:    

        for forcing in forcings: 

            for scenario in scenarios:

                filedir, filenames = get_isimip_simulation_name_dir(variable, datadir, model, forcing, scenario, return_single=False, timestep='daily')
                filenames.sort()

                print('Number of files: '+str(len(filenames)))
                filename_monthly_merged = model.lower()+'_'+forcing+'_ewembi_'+scenario+'_'+d_soc_scenario[scenario]+'_co2_'+variable+'_global_monthly_'+d_scenario_years[scenario]+'.nc4'
                filenames_monthly = [filedir+filename.replace('daily','monthly') for filename in filenames]
                
                if not os.path.exists(filedir+filename_monthly_merged): 
                    print('calculating monthly means for '+model + ' ' + forcing + ' '+scenario)
                    
                    for i,filename in enumerate(filenames): 
                        
                        if not os.path.exists(filedir+filename.replace('daily','monthly')): 
                            
                            os.system('cdo monmean '+filedir+filename+' '+filedir+filename.replace('daily','monthly'))
                    
                    # merge all monthly files
                    print('merging '+model + ' ' + forcing + ' '+scenario)
                    os.system('cdo mergetime '+' '.join(filenames_monthly) +' '+filedir+filename_monthly_merged)
                
                os.system( 'rm '+' '.join(filenames_monthly))
      
