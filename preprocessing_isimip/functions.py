# ----------------------------------------------------------------
# Functions for water scarcity
#
# ----------------------------------------------------------------

import os 
import xarray as xr
import numpy as np
    
from init import *
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import glob
from netCDF4 import Dataset
import regionmask


# Disable printing
def block_print():
    sys.stdout = open(os.devnull, 'w')

# Restore printing
def enable_print():
    sys.stdout = sys.__stdout__calc

# get filename and location of simulation for one model + forcing combination for all defined variables
def get_isimip_simulation_name_dir(variable, datadir, model, forcing, scenario, return_single=True, timestep='monthly'):

    # variable define end
    variable = variable +'_'
    
    # open historical data
    modeldir   = datadir + 'OutputData/water_global/'+model+'/'+forcing+'/' 


    filedir = modeldir + d_scenario_folder[scenario] +'/'


    # retrieve all filenames in directory
    filenames_all = os.listdir(modeldir + d_scenario_folder[scenario])


    # keep only filename for scenario, socioeconomic scenario, variable and timestep
    filename_filtered = [s for s in filenames_all if ((scenario in s) and (d_soc_scenario[scenario] in s) and (variable in s) and (timestep in s))] 


    # if filename exists in filtered version
    if len(filename_filtered) > 0: 
        
        if return_single == True:
            filename = filename_filtered[0]

            # check if file exists
            if os.path.exists((filedir + filename).encode('unicode_escape')): 
                return filedir, filename

            else: 
                print(filedir + filename + ' does not exist.')

                return ('', 'false')
        else: 
            return filedir, filename_filtered

    else: 
        print('No '+ model+' '+forcing+' ' + scenario+ ' simulations for '+ variable)

        return ('', 'false')

# get filename and location of simulation for output simulations
def get_simulation_name_dir(variable, datadir, model, forcing, scenario, timestep='monthly', socscen = True):

    # variable define end
    variable = variable +'_'
    
    timestep = timestep+'_'
    
    filedir = datadir + '/'+model.lower()+'/'

    # retrieve all filenames in directory
    filenames_all = os.listdir(filedir)

    # keep only filename for scenario, socioeconomic scenario, variable and timestep
    filename_filtered = [s for s in filenames_all if ((forcing in s) and (scenario in s)  and (variable in s) and (timestep in s))] 

    # also filter based on socioeconomic scenario (this is done in intermediate, but not in lifetime dir)
    if socscen:
        filename_filtered = [s for s in filename_filtered if (d_soc_scenario[scenario] in s)] 
        
    # if filename exists in filtered version
    if len(filename_filtered) > 0: 
        filename = filename_filtered[0]

        # check if file exists
        if os.path.exists((filedir + filename).encode('unicode_escape')): 

            return filedir, filename

        else: 
            print(filedir + filename + ' does not exist.')

            return (filedir, 'false')

    else: 
        print('No '+ model+' '+forcing+' ' + scenario+ ' simulations for '+ variable)

        return (filedir, 'false')


# Load simulations for one model + forcing combination for all defined variables
def load_simulation(variable, datadir, model, forcing, scenario, flag_postprocessed=False):

    # for isimip 2b, decode times needs to be false, the phase variable itself is set in the init file
    if phase == '3b': 
        decode_times = True 
    else: 
        decode_times = False
    
    # load raw isimip simulation
    if not flag_postprocessed:


        # open historical data
        modeldir   = datadir + 'OutputData/water_global/'+model+'/'+forcing+'/'            
        filedir = modeldir + d_scenario_folder[scenario] +'/'
        

    # load postprocessed simulation (from other directory structure)
    else: 
        
        filedir = './data/intermediate/'+model.lower()+'/'
        if phase == '2b': 
            filedir = './data/2b/intermediate/'+model.lower()+'/'
        
    # retrieve all filenames in directory
    filenames_all = os.listdir(filedir)   
    
    # keep only filename for scenario, socioeconomic scenario, variable and timestep
    filename_filtered = [s for s in filenames_all if ((scenario in s) and (d_soc_scenario[scenario] in s) and (variable+'_' in s) and (timestep in s) and (forcing+'_' in s))] 


    
    # problem with dis_upstream and dis_allupstream as the variable in the filename and the variable in the ncfile do not match -- to solve by renaming the file, solve now by manually changing. 
    if variable == 'dis_directupstream': 
        variable = 'dis_upstream'
        
        
    # if filename exists in filtered version
    if len(filename_filtered) > 0: 
        filename = filename_filtered[0]

        # check if file exists
        if os.path.exists((filedir + filename).encode('unicode_escape')): 

            ds = xr.open_dataset(filedir + filename, decode_times=decode_times)
            #print(filedir+filename)
            da = ds[variable]

            # for isimip 2b, manually assign time periods, as these can not be read in 
            if phase=='2b': 
                
                da['time'] = d_time_ncfile[scenario]
            return da
            
        else: 
            print(filedir + filename + ' does not exist.')

            return []

    else: 
        print('No '+ model+' '+forcing+' ' + scenario+ ' simulations for '+ variable)

        return []

    
# load simulations in dict pper modelforcing and scenario
def load_simulations_permodelforcing(variable, models, forcings, scenarios, datadir, flag_postprocessed=False):
    
    d_da_permodelforcing =  {}
    
    # loop over models
    for model in models:
    
        # loop over GCM forcings
        for forcing in forcings: 
            
            d_da_perscenario = {}
            
            # loop over scenarios
            for scenario in scenarios:

                # load simulation per variable
                da = load_simulation(variable, datadir, model, forcing, scenario, flag_postprocessed)

                d_da_perscenario[scenario] = da

            d_da_permodelforcing[model+'_'+forcing] = d_da_perscenario
            
    return d_da_permodelforcing


# function to create a mask based on lon, lat and ar6 regions
# returns dataarray with 3D mask with dimensions lon, lat and region
def create_regionsmask(da):

    regions_mask = regionmask.defined_regions.ar6.land.mask_3D(da.lon, da.lat)
    
    return regions_mask

# Preprocess potential total water withdrawal and save them as intermediary datafiles
def process_ptotww(models, forcings, scenarios): 


    for model in models:

        for forcing in forcings: 

            for scenario in scenarios:

                intermediatedir = outdir+'intermediate/'+model.lower()+'/'
                
                if model == 'CWatM' or model == 'WaterGAP2-2e': 

                    filedir, filename = get_isimip_simulation_name_dir('ptotww', datadir, model, forcing, scenario)

                    # only process file if not existing yet
                    if not os.path.exists(intermediatedir+filename): 

                        print('processing ptotww for '+model + ' ' + forcing + ' '+scenario)

                        # copy over ptotww to intermediate use
                        os.system('cp '+filedir+filename+' '+intermediatedir+filename)


                if model == 'H08': 
                    
                    filedir, filename = get_isimip_simulation_name_dir('amanuse', datadir, model, forcing, scenario)

                    if not os.path.exists(intermediatedir+filename.replace('amanuse', 'ptotww')): 
                        print('processing ptotww for '+model + ' ' + forcing + ' '+scenario)

                    # load water use of seperate sectors
                        amanuse = load_simulation('amanuse', datadir, model, forcing, scenario)
                        adomuse = load_simulation('adomuse', datadir, model, forcing, scenario)
                        airruse = load_simulation('airruse', datadir, model, forcing, scenario)

                        # calculate total potential water withdrawal
                        ptotww = amanuse + adomuse + airruse
                        ptotww.attrs = {'long_name': 'Total Potential Water Withdrawal (all sectors)',
                                        'standard_name': 'ptotww',
                                        'units': 'kg m-2 s-1', 
                                        'comment': 'Sum of actual water withdrawal of sectors: manufacturing, domestic and irrigation'}

                        # Save the variable into intermediate netCDF file
                        ptotww.to_dataset(name='ptotww').to_netcdf(intermediatedir+filename.replace('amanuse', 'ptotww'))

                        
# Calculate total potential water withdrawal for ISIMIP2b simulations
def calc_ptotww_2b(models, forcings, scenarios, da_cellarea): 

    # calculate potential total water withdrawal for isimip2b simulations, based on defenitions in main_2b noteboook
    for model in models:

        for forcing in forcings: 

            for scenario in scenarios:

                intermediatedir = outdir+'intermediate/'+model.lower()+'/'

                if model == 'H08': 

                    filedir, filename = get_isimip_simulation_name_dir('amanww', datadir, model, forcing, scenario)

                    if not os.path.exists(intermediatedir+filename.replace('amanww', 'ptotww')): 
                        print('processing ptotww for '+model + ' ' + forcing + ' '+scenario)

                    # load water use of seperate sectors
                        amanww = load_simulation('amanww', datadir, model, forcing, scenario)
                        adomww = load_simulation('adomww', datadir, model, forcing, scenario)
                        pirrww = load_simulation('pirrww', datadir, model, forcing, scenario)

                        # calculate total potential water withdrawal
                        ptotww = amanww + adomww + pirrww
                        ptotww.attrs = {'long_name': 'Total Potential Water Withdrawal (all sectors)',
                                        'standard_name': 'ptotww',
                                        'units': 'kg m-2 s-1', 
                                        'comment': 'Sum of actual water withdrawal of sectors: amanww + adomww + pirrww'}

                        # Save the variable into intermediate netCDF file
                        ptotww.to_dataset(name='ptotww').to_netcdf(intermediatedir+filename.replace('amanww', 'ptotww'))

                
                elif model == 'CWatM': 

                    filedir, filename = get_isimip_simulation_name_dir('adomuse', datadir, model, forcing, scenario)

                    if not os.path.exists(intermediatedir+filename.replace('adomuse', 'ptotww')): 
                        print('processing ptotww for '+model + ' ' + forcing + ' '+scenario)

                    # load water use of seperate sectors
                        adomuse = load_simulation('adomuse', datadir, model, forcing, scenario)
                        ainduse = load_simulation('ainduse', datadir, model, forcing, scenario,)
                        pirrww = load_simulation('pirrww', datadir, model, forcing, scenario)

                        # calculate total potential water withdrawal
                        ptotww = adomuse + ainduse + pirrww
                        ptotww.attrs = {'long_name': 'Total Potential Water Withdrawal (all sectors)',
                                        'standard_name': 'ptotww',
                                        'units': 'kg m-2 s-1', 
                                        'comment': 'Sum of actual water withdrawal of sectors: ainduse + adomuse + pirrww'}

                        # Save the variable into intermediate netCDF file
                        ptotww.to_dataset(name='ptotww').to_netcdf(intermediatedir+filename.replace('adomuse', 'ptotww'))
                
                
                # TO DO: check individual values!!!!
                elif model == 'LPJmL' or model == 'MATSIRO': 
                    filedir, filename = get_isimip_simulation_name_dir('pirrww', datadir, model, forcing, scenario)

                    if not os.path.exists(intermediatedir+filename.replace('pirrww', 'ptotww')): 
                        print('processing ptotww for '+model + ' ' + forcing + ' '+scenario)

                        # load potential irrigation water withdrawal (available for both)
                        pirrww = load_simulation('pirrww', datadir, model, forcing, scenario)

                        if scenario == 'historical': 
                            inputdir = datadir+'/InputData/water_abstraction/histsoc/'
                            domww_m3yr = xr.open_dataset(inputdir+'domww_histsoc_annual_1901-2005.nc', decode_times=False)['domww']
                            domww_m3yr = domww_m3yr.assign_coords(time = pd.date_range(start='1901', end='2006', freq='Y'))
                            indww_m3yr = xr.open_dataset(inputdir+'indww_histsoc_annual_1901-2005.nc', decode_times=False)['indww']
                            indww_m3yr = indww_m3yr.assign_coords(time = pd.date_range(start='1901', end='2006', freq='Y'))

                            # extend historical period based on pirr period
                            domww_m3yr = match_time_period(domww_m3yr,pirrww)
                            indww_m3yr = match_time_period(indww_m3yr,pirrww)

                        else: 
                            inputdir= datadir+'/InputData/water_abstraction/'+scenario+'soc/'

                            domww_m3yr = xr.open_dataset(inputdir+'domww_rcp26soc_annual_2006-2099.nc', decode_times=False)['domww']
                            domww_m3yr = domww_m3yr.assign_coords(time = pd.date_range(start='2006', end='2100', freq='Y'))
                            indww_m3yr = xr.open_dataset(inputdir+'indww_rcp26soc_annual_2006-2099.nc', decode_times=False)['indww']
                            indww_m3yr = indww_m3yr.assign_coords(time = pd.date_range(start='2006', end='2100', freq='Y'))


                        # transform the units from m3/year to mm/m2s
                        domww = (domww_m3yr * 1000 /(da_cellarea * sec_per_year))
                        domww = domww.where(domww<1e34,np.nan)
                        indww = (indww_m3yr * 1000 /(da_cellarea * sec_per_year))
                        indww = indww.where(indww<1e34,np.nan)
                        
                        # Convert input water aithdrawal to to monthly. 
                        domww = conv_annual2monthly(domww)
                        indww = conv_annual2monthly(indww)
                      
                        # calculate total potential water withdrawal
                        ptotww = domww + indww + pirrww
                        ptotww.attrs = {'long_name': 'Total Potential Water Withdrawal (all sectors)',
                                        'standard_name': 'ptotww',
                                        'units': 'kg m-2 s-1', 
                                        'comment': 'Sum of actual water withdrawal of sectors: domww + indww + pirrww (domww and indww from ISIMIP input data)'}

                        # Save the variable into intermediate netCDF file
                        ptotww.to_dataset(name='ptotww').to_netcdf(intermediatedir+filename.replace('pirrww', 'ptotww'))
                        
# calculate total water withdrawal, using the corrected raw simulations for H08 (new future simulations, and corrected historical simulations)
# the other models use dom and ind demands from H08. 
def calc_ptotww_2b_corrected_h08(models, forcings, scenarios, da_cellarea):

    # calculate potential total water withdrawal for isimip2b simulations, based on defenitions in main_2b noteboook
    for model in models:

        for forcing in forcings: 

            for scenario in scenarios:

                intermediatedir = outdir+'intermediate/'+model.lower()+'/'

                if model == 'H08': 

                    filedir, filename = get_isimip_simulation_name_dir('pirrww', datadir, model, forcing, scenario)
                    
                    if scenario != 'historical': 
                        flag_postprocessed_pirrww_h08 = True
                    else: 
                        flag_postprocessed_pirrww_h08 = False
                        
                    
                    
                    if not os.path.exists(intermediatedir+filename.replace('pirrww', 'ptotww').replace(forcing, forcing+'_ewembi')): 
                        print('processing ptotww for '+model + ' ' + forcing + ' '+scenario)


                        # load water use of separate sectors
                        amanww = load_simulation('amanww', datadir, model, forcing, scenario, flag_postprocessed=True)
                        adomww = load_simulation('adomww', datadir, model, forcing, scenario, flag_postprocessed=True)
                        pirrww = load_simulation('pirrww', datadir, model, forcing, scenario, flag_postprocessed=flag_postprocessed_pirrww_h08)

                        # calculate total potential water withdrawal
                        ptotww = amanww + adomww + pirrww
                        ptotww.attrs = {'long_name': 'Total Potential Water Withdrawal (all sectors)',
                                        'standard_name': 'ptotww',
                                        'units': 'kg m-2 s-1', 
                                        'comment': 'Sum of actual water withdrawal of sectors: amanww + adomww + pirrww'}

                        # Save the variable into intermediate netCDF file
                        ptotww.to_dataset(name='ptotww').to_netcdf(intermediatedir+filename.replace('pirrww', 'ptotww').replace(forcing, forcing+'_ewembi'))


                
                elif model == 'CWatM': 

                    filedir, filename = get_isimip_simulation_name_dir('pirrww', datadir, model, forcing, scenario)

                    if not os.path.exists(intermediatedir+filename.replace('pirrww', 'ptotww')): 
                        print('processing ptotww for '+model + ' ' + forcing + ' '+scenario)

                    # load water use of seperate sectors

                        amanww = load_simulation('amanww', datadir, 'H08', forcing, scenario, flag_postprocessed=True)
                        adomww = load_simulation('adomww', datadir, 'H08', forcing, scenario, flag_postprocessed=True)                        
                        pirrww = load_simulation('pirrww', datadir, model, forcing, scenario)

                        # calculate total potential water withdrawal
                        ptotww = amanww + adomww + pirrww
                        ptotww.attrs = {'long_name': 'Total Potential Water Withdrawal (all sectors)',
                                        'standard_name': 'ptotww',
                                        'units': 'kg m-2 s-1', 
                                        'comment': 'Sum of actual water withdrawal of sectors: amanww + adomww + pirrww,  amanww and adomww from H08 simulation'}

                        # Save the variable into intermediate netCDF file
                        ptotww.to_dataset(name='ptotww').to_netcdf(intermediatedir+filename.replace('pirrww', 'ptotww'))
                
                
                # TO DO: check individual values!!!!
                elif model == 'LPJmL' or model == 'MATSIRO': 
                    filedir, filename = get_isimip_simulation_name_dir('pirrww', datadir, model, forcing, scenario)

                    if not os.path.exists(intermediatedir+filename.replace('pirrww', 'ptotww')): 
                        print('processing ptotww for '+model + ' ' + forcing + ' '+scenario)

                        # load potential irrigation water withdrawal (available for both)
                        pirrww = load_simulation('pirrww', datadir, model, forcing, scenario)
                        amanww = load_simulation('amanww', datadir, 'H08', forcing, scenario, flag_postprocessed=True)
                        adomww = load_simulation('adomww', datadir, 'H08', forcing, scenario, flag_postprocessed=True)                        

                        # calculate total potential water withdrawal
                        ptotww = adomww + amanww + pirrww
                        ptotww.attrs = {'long_name': 'Total Potential Water Withdrawal (all sectors)',
                                        'standard_name': 'ptotww',
                                        'units': 'kg m-2 s-1', 
                                        'comment': 'Sum of actual water withdrawal of sectors: domww + indww + pirrww (domww from adomww and indww from amanww from H08 simulation'}

                        # Save the variable into intermediate netCDF file
                        ptotww.to_dataset(name='ptotww').to_netcdf(intermediatedir+filename.replace('pirrww', 'ptotww'))      
                        
# calculate water scarcity index according to Veldkamp et al., 2017 NComm
# this method deviates in the calculation of water availability
def calc_waterscarcity_index(model, forcing, scenario, landmask, da_cellarea, flag_save=True):  
     
    # load simulations

    # load grid cell runoff
    qtot = load_simulation('qtot', datadir, model, forcing, scenario)

    # load discharge into gricell from directly upstream grid cells
    dis_upstream = load_simulation('dis_directupstream', datadir, model, forcing, scenario, flag_postprocessed=True)

    # calculations

    # calculate water availability (q_avail)
    q_avail = dis_upstream * 1000/da_cellarea + qtot # UNITS! dis is in m³/s, convert to mm/m²s

    del qtot, dis_upstream
    # determine environmental flow requirement (EFR) as defined in Veldkamp et al., 2017 Supplementary info
    # transforming units from m³/s to kg/m²s

    # THIS CAUSES THE KERNEL TO DIE
    # load discharge from grid cell to determine environmental flow requirements -- currently not used!!
    # dis  = load_simulation('dis' , datadir, model, forcing, scenario)
    
   
    efr = get_efr(q_avail) # *1000/da_cellarea

    # del dis
    
    
    # load water withdrawal data (demand)
    ptotww = load_simulation('ptotww', datadir, model, forcing, scenario, flag_postprocessed=True)
    
    # calculate water scarcity index and mask for the land values
    waterscarcity_index = (ptotww / (q_avail- efr))
    
    del ptotww
    waterscarcity_index = waterscarcity_index.where(landmask)
    
    if flag_save: 
        # save water scarcity index
        filedir, filename = get_isimip_simulation_name_dir('qtot', datadir, model, forcing, scenario)
        intermediatedir = outdir+'/intermediate/'+model.lower()+'/'

        waterscarcity_index.to_dataset(name='wsindex').to_netcdf(intermediatedir+filename.replace('qtot', 'waterscarcity_index'))
        q_avail.to_dataset(name='q_avail').to_netcdf(intermediatedir+filename.replace('qtot', 'q_avail'))

    return waterscarcity_index



## CHECK THIS FUCTION AGAIN, KERNEL ALWAYS DIES
def get_efr(qtot):

    # compute yearly means
    qtot_ymean = qtot.groupby('time.year').mean()

    # intialise lowflow, mediumflow and highflow 
    mask_lowflow    = xr.DataArray(data = np.zeros_like(qtot.values),coords=qtot.coords, dims=qtot.dims)
    mask_mediumflow = xr.DataArray(data = np.zeros_like(qtot.values),coords=qtot.coords, dims=qtot.dims)
    mask_highflow   = xr.DataArray(data = np.zeros_like(qtot.values),coords=qtot.coords, dims=qtot.dims)


    # for every year, calculate masks with corresponding annual means
    
    for group_name, group_da in qtot.groupby('time.year'):
        mask_lowflow[mask_lowflow.time.dt.year == group_name] = group_da <= 0.4*qtot_ymean.where(qtot_ymean.year==group_name, drop=True).squeeze()
        mask_mediumflow[mask_mediumflow.time.dt.year == group_name] = group_da > 0.8*qtot_ymean.where(qtot_ymean.year==group_name, drop=True).squeeze()

    mask_highflow = mask_lowflow + mask_mediumflow == 0 # mediumflow is in between high flow and low flow


    # determine environmental flow conditions
    efr =   qtot

    # apply inverse masks, as where function assigns values to places where mask is not true
    efr = efr.where(mask_lowflow==0   ,  0.6 * qtot)
    efr = efr.where(mask_mediumflow==0,  0.3 * qtot)
    efr = efr.where(mask_highflow==0  ,  0.45 * qtot)
    
    return efr


# calculate water index exposure - for monthly data and save to output folder

def calc_scarcity_exposure(waterscarcity_index, threshold, landmask, timestep, outdir, model, forcing, scenario, flag_save=True):

    # water scarcity occurs when water scarcity index is larger than 1. 
    waterscarcity_exposure = (waterscarcity_index > threshold).astype('float')

    #waterscarcity_exposure =  waterscarcity_exposure.where(landmask)
    
    if flag_save: 
        # convert data array to dataset
        ds_waterscarcity_exposure = waterscarcity_exposure.to_dataset(name='exposure')


        # if annual, rename time dimensions from year to time for consistency
        if timestep == 'annual': 
            ds_waterscarcity_exposure = ds_waterscarcity_exposure.rename_dims({'year':'time'})

        # define filename and output directory of model
        filename_out = model.lower()+'_'+forcing+'_'+scenario+'_waterscarcity_global_'+timestep+'_landarea_'+d_senario_period[scenario]+'.nc4'
        model_outdir = outdir+'/waterscarcity/'+model.lower()+'/'


        # print to screen
        print('saving '+filename_out)


        # Create a new directory because it does not exist 
        if not os.path.exists(model_outdir):
            os.makedirs(model_outdir)


        # save to output location
        ds_waterscarcity_exposure.to_netcdf(model_outdir+filename_out)

    else:

        return waterscarcity_exposure
    
# calculate water scarcity duration    
# calculate water scarcity duration    
def calc_scarcity_duration(waterscarcity_index, threshold, landmask, outdir, model, forcing, scenario, flag_save=True):

    # water scarcity occurs when water scarcity index is larger than 1. 
    waterscarcity_exposure = waterscarcity_index > threshold

    # calculate the different variables related to duration

    # number of water scarce months
    waterscarcity_nmonths = waterscarcity_exposure.groupby('time.year').sum()

    # if first month is water scarce
    waterscarcity_first = (waterscarcity_exposure.where(waterscarcity_exposure.time.dt.month== 1) ==1).groupby('time.year').sum()

    # if last month is water scarce 
    waterscarcity_last = (waterscarcity_exposure.where(waterscarcity_exposure.time.dt.month== 12) ==1).groupby('time.year').sum()


    # calculate number of water scarce events
    nevents_list = []
    for year in np.unique(waterscarcity_exposure.time.dt.year.values): 

        waterscarcity_year = waterscarcity_exposure.where(waterscarcity_exposure.time.dt.year== year, drop=True)

        # calculate number of events (row of booleans, if difference is taken, -1 (1-0) represents end of event. 
        nevents_year = ((waterscarcity_year*1).diff(dim='time') ==-1).sum("time")

        # this does not include the end of the year, so add that as separate event if last month is 1
        nevents_year = nevents_year.where(waterscarcity_year[-1,:,:]==0, nevents_year+1)


        # correct for events that are multiyear spanning. 

        nevents_list.append(nevents_year)

    waterscarcity_nevents = xr.concat(nevents_list, dim="year").drop('time').assign_coords({"year":waterscarcity_nmonths.year.values})


    # loactions with multi-year events. Subtract 1 here!!! for nevents!!!!

    boolean_last = waterscarcity_last.shift(year=1)==1
    boolean_first = waterscarcity_first==1
    waterscarcity_nevents.where( np.invert(boolean_last*boolean_first), waterscarcity_nevents-1)


    ds_waterscarcity_duration = waterscarcity_nmonths.to_dataset(name='nmonths')
    ds_waterscarcity_duration['nevents'] = waterscarcity_nevents



    if flag_save: 

        # rename time dimensions from year to time for consistency

        ds_waterscarcity_duration = ds_waterscarcity_duration.rename_dims({'year':'time'}).where(landmask).drop('year')

        # define filename and output directory of model
        filename_out = model.lower()+'_'+forcing+'_'+scenario+'_waterscarcityduration_global_annual_landarea_'+d_senario_period[scenario]+'.nc4'

        model_outdir = outdir+'/waterscarcity/'+model.lower()+'/'


        # print to screen
        print('saving '+filename_out)


        # Create a new directory because it does not exist 
        if not os.path.exists(model_outdir):
            os.makedirs(model_outdir)


        # save to output location
        ds_waterscarcity_duration.to_netcdf(model_outdir+filename_out)

    else:

        return ds_waterscarcity_duration        

# convert monthly water scarcity exposure to annual water scarcity
def waterscarcity_monthly_to_annual(models, forcings, scenarios, aggregation_mode = 'binary' ): 
    # loop over models
    for model in models:

        # loop over GCM forcings
        for forcing in forcings: 

            for scenario in scenarios:

                # convert monthly scale to annual time scale using cdo (for performance)
                # open monthly file
                filedir, filename = get_simulation_name_dir('waterscarcity', lifetimedir+'waterscarcity', model, forcing, scenario, 'monthly', socscen = False)
                
                if os.path.exists(filedir+filename): 
                    
                    if os.path.exists(filedir+filename.replace('monthly','annual')): 
                        os.system('rm '+filedir+filename.replace('monthly','annual'))
                        
                    if aggregation_mode=='binary': 
                        
                        os.system('cdo yearmax '+filedir+filename+' '+filedir+filename.replace('monthly','annual'))
                    elif aggregation_mode == 'sum': 
                        os.system('cdo yearsum '+filedir+filename+' '+filedir+filename.replace('monthly','annualsum'))
                       
                    

# save new variable to lifetimedirectory
def save_to_lifetimedir(da, variable, model, forcing, scenario, annual=False):
    # save water scarcity index
    filedir, filename = get_isimip_simulation_name_dir('qtot', datadir, model, forcing, scenario)
    save_dir =  lifetimedir+variable+'/'+model.lower()+'/'
    
    # Create a new directory because it does not exist 
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    filename = save_dir+filename.replace('qtot', variable)
    
    if annual: 
        filename = filename.replace('monthly','annual')
        
    if os.path.exists(filename): 
        os.system('rm '+filename)
        
        
    da.to_dataset(name=variable).to_netcdf(filename)


# convert monthly to annual for intermediate variable
def intermediate_monthly_to_annual(variable, models, forcings, scenarios, binary=True): 
    # loop over models
    for model in models:

        # loop over GCM forcings
        for forcing in forcings: 

            for scenario in scenarios:

                # convert monthly scale to annual time scale using cdo (for performance)
                # open monthly file
                filedir, filename = get_simulation_name_dir(variable, outdir+'intermediate', model, forcing, scenario, 'monthly')
                
                if os.path.exists(filedir+filename): 
                    
                    if os.path.exists(filedir+filename.replace('monthly','annual')): 
                        os.system('rm '+filedir+filename.replace('monthly','annual'))
                        
                    os.system('cdo yearmax '+filedir+filename+' '+filedir+filename.replace('monthly','annual'))
                    
# open hist and future population datasets and merge them to get one series
def merge_pop(popdir, fn_pop_hist, fn_pop_fut ):

    pop_hist  = xr.open_dataset(fn_pop_hist, decode_times=False)
    pop_fut  = xr.open_dataset(fn_pop_fut, decode_times=False)

    pop = pop_hist.merge(pop_fut)['number_of_people']
    pop = pop_hist.merge(pop_fut)['number_of_people']

    pop = xr.concat([pop,pop[-1,:,:]], dim='time')
    pop = pop.assign_coords(time = pd.date_range(start='1861', end='2101', freq='Y'))
    
    pop.to_dataset(name='pop').to_netcdf(outdir+'/intermediate/'+'population_histsoc_rcp26soc_0p5deg_annual_1861-2100.nc4')

# extend population data from yearly to monthly by using the same constant value
def pop_yearly_to_monthly(popdir, fn_pop_hist, fn_pop_fut ): 
    
    pop_hist  = xr.open_dataset(fn_pop_hist, decode_times=False)
    pop_fut  = xr.open_dataset(fn_pop_fut, decode_times=False)

    pop = pop_hist.merge(pop_fut)['number_of_people']
    pop = pop_hist.merge(pop_fut)['number_of_people']

    pop = xr.concat([pop,pop[-1,:,:]], dim='time')
    pop = pop.assign_coords(time = pd.date_range(start='1861', end='2100', freq='Y'))
    pop_values = pop.values

    for y in range(pop_values.shape[0]): 

        pop_values_monthly = np.empty((12,pop_values.shape[1],pop_values.shape[2]))

        pop_values_monthly = np.tile(pop_values[y,:,:],(12,1,1))

        if y ==0: 

            pop_values_new = pop_values_monthly
        else: 

            pop_values_new =  np.concatenate((pop_values_new, pop_values_monthly),axis=0)

    pop_monthly = xr.DataArray(data=pop_values_new, dims=pop.dims, coords= {'lat':pop.lat, 'lon':pop.lon, 'time':pd.date_range(start='1861', end='2100', freq='M') })
    pop_monthly.to_dataset(name='pop').to_netcdf(outdir+'/intermediate/'+'population_histsoc_rcp26soc_0p5deg_monthly_1861-2100.nc4')
    pop.to_dataset(name='pop').to_netcdf(outdir+'/intermediate/'+'population_histsoc_rcp26soc_0p5deg_annual_1861-2100.nc4')


def load_pop(fn_pop): 
    
    if phase == '2b': 
        popdir = scriptsdir+'/data/intermediate/'
    else: 
        popdir = outdir+'/intermediate/'
    return xr.open_dataset(popdir+fn_pop)['pop']   


# convert annual dataset to monthly, and all monthly values per year add up to yearly total

def conv_annual2monthly(da): 

    values = da.values

    # time is third dimension -> move to first
    if np.shape(values)[2]<np.shape(values)[0]:
        values = np.swapaxes(values,0,2)
        values = np.swapaxes(values,1,2) # make sure lat lon are in right order too


    for y in range(values.shape[0]): 

        values_monthly = np.empty((12,values.shape[1],values.shape[2]))

        values_monthly = np.tile(values[y,:,:]/12,(12,1,1)) # divide annual value by 12 to get monthly value

        if y ==0: 

            values_new = values_monthly
        else: 

            values_new =  np.concatenate((values_new, values_monthly),axis=0)

    da_monthly = xr.DataArray(data=values_new, dims=('time','lat', 'lon'), coords= {'lat':da.lat, 'lon':da.lon, 'time':pd.date_range(start=str(da.time.dt.year.min().values), end=str(da.time.dt.year.max().values+1), freq='M') })
    return da_monthly


# extend the da's time period back into time (based on da_target first year) by concatenating the first year
# returns yearly values
def match_time_period(da,da_target):

    # get beginning and end of file
    start_year = str(da.time.dt.year.min().values)
    end_year = str(da.time.dt.year.max().values)

    # what the real period needs to be
    target_start_year = str(da_target.time.dt.year.min().values)
    target_end_year = str(da_target.time.dt.year.max().values)

    if target_start_year < start_year: 

        nyears = int(start_year) - int(target_start_year)

        for n in range(1,nyears+1): 

            # workaround: time dimesion of da switches places after concatenating
            if n == 1: 
                da =  xr.concat([da[0,:,:],da], dim='time')

            else: 
                da =  xr.concat([da[:,:,0],da], dim='time')
        # assign right time dimension (therefore reduce from monthly, to annual scale)
        da['time'] = pd.date_range(start=str(target_start_year), end=str(int(target_end_year)+1), freq='Y')
        
        return da

# calculate the falkenmark index: water availability per capita (m³/cap)
# annual timescale
def calc_falkenmark_index(model, forcing, scenario, pop, da_cellarea, landmask, flag_save=True):

    # load grid cell runoff
    qtot = load_simulation('qtot', datadir, model, forcing, scenario)

    # load discharge into gricell from directly upstream grid cells
    dis_upstream = load_simulation('dis_directupstream', datadir, model, forcing, scenario, flag_postprocessed=True)

    # load discharge from grid cell to determine environmental flow requirements
    dis  = load_simulation('dis' , datadir, model, forcing, scenario)


    # calculations

    # calculate water availability (q_avail)
    q_avail_withoutefr = dis_upstream + qtot * da_cellarea  / 1000 # UNITS! dis is in m³/s, convert to mm/m²
    
    # subtract envrironmental flow requirement of water availability
    q_avail =  q_avail_withoutefr - get_efr(q_avail_withoutefr)
    q_avail = q_avail.assign_coords(time = d_time_ncfile[scenario])


    days_in_month = []
    for q_avail_month in q_avail.time.values:
        days_in_month.append(pd.to_datetime(q_avail_month).days_in_month)

    secs_in_month = np.expand_dims(np.array(days_in_month) * 24*60*60, axis=(1,2))

    q_avail_annual = (q_avail * secs_in_month).groupby('time.year').sum()

    # fix the time dimension. 
    q_avail_annual = q_avail_annual.rename({'year':'time'}).assign_coords(time = pd.date_range(start=str(q_avail_annual.year.min().values), end=str(q_avail_annual.year.max().values+1), freq='Y'))
    #
    pop_selected = pop.where(landmask)
    falkenmark_index = q_avail_annual.where(pop_selected) / pop_selected
    
    #falkenmark_index = falkenmark_index.where(pop>0, q_avail_annual).where(landmask)

    if flag_save: 
        # save water scarcity index
        filedir, filename = get_isimip_simulation_name_dir('qtot', datadir, model, forcing, scenario)
        intermediatedir = outdir+'/intermediate/'+model.lower()+'/'

        start_year = str(falkenmark_index.time.min().values.astype('datetime64[Y]').astype(int) + 1970)
        end_year = str(falkenmark_index.time.max().values.astype('datetime64[Y]').astype(int) + 1970)
        filename_new = (filename[:-12:] + start_year + '_'+end_year+'.nc').replace('qtot', 'falkenmark_index')

        falkenmark_index.to_dataset(name='falkenmark_index').to_netcdf(intermediatedir+filename_new)


        return falkenmark_index
    
    
# calculate exposure if water scarce annual water availabily (m³/cap) is less than the threshold
def calc_falkenmark_exposure(falkenmark_index, threshold, landmask,  outdir, model, forcing, scenario, flag_save=True):
    
        # water scarcity occurs when falkenmark index is smaller than threshold
    falkenmark_exposure = falkenmark_index.groupby('time.year').sum() < threshold

    falkenmark_exposure =  falkenmark_exposure.where(falkenmark_index.max(dim='time')).where(landmask)
    
    if flag_save: 
        # convert data array to dataset
        ds_falkenmark_exposure = falkenmark_exposure.to_dataset(name='exposure')


        ds_falkenmark_exposure = ds_falkenmark_exposure.rename_dims({'year':'time'})
        ds_falkenmark_exposure = ds_falkenmark_exposure.rename({'year':'time'})

        start_year = str(falkenmark_index.time.min().values.astype('datetime64[Y]').astype(int) + 1970)
        end_year = str(falkenmark_index.time.max().values.astype('datetime64[Y]').astype(int) + 1970)
        
        # define filename and output directory of model
        filename_out = model.lower()+'_'+forcing+'_'+scenario+'_falkenmark_global_annual_landarea_'+start_year+'_'+end_year+'.nc4'
        model_outdir = outdir+'/falkenmark/'+model.lower()+'/'


        # print to screen
        print('saving '+filename_out)


        # Create a new directory because it does not exist 
        if not os.path.exists(model_outdir):
            os.makedirs(model_outdir)


        # save to output location
        ds_falkenmark_exposure.to_netcdf(model_outdir+filename_out)


    return falkenmark_exposure


# load waterscarcity exposure and calculate the total area exposed to water scarcity per forcing, scenario and model
def calc_totalarea_exposed(variable, landmask, da_cellarea, models, scenarios, forcings, timestep='annual'):

    # load total land area exposed
    total_landarea = da_cellarea.where(landmask).sum()

    # loop over models
    model_list = []
    for model in models:    
        scenario_list = []
        for scenario in scenarios:
            ds_forcings = xr.Dataset()
            # loop over GCM forcings

            area_exposed_list = []
            for forcing in forcings: 
                print(model +' '+scenario+' '+forcing)

                # define filename and output directory of model

                filedir, filename = get_simulation_name_dir(variable, lifetimedir+variable, model, forcing, scenario, timestep, socscen=False)

                if filename!='false': 
                    # save to output location
                    ds_waterscarcity_exposure = xr.open_dataset(filedir+filename, engine='netcdf4')
                    da_waterscarcity_exposure = ds_waterscarcity_exposure['exposure'].where(landmask)
                    area_exposed = (da_waterscarcity_exposure * da_cellarea).sum(dim=('lat','lon'))/total_landarea *100

                else: 
                    if scenario == 'historical': 
                        empty_array = np.empty_like(area_exposed_list[0].values)

                    else: 
                        empty_array = np.empty_like(area_exposed_list[1].values)
                    
                    empty_array[:] = np.nan    
                    area_exposed = xr.DataArray(empty_array, dims=('time'))
                area_exposed_list.append(area_exposed)


            ds_scenario = xr.concat(area_exposed_list, dim='forcing').assign_coords({'forcing':forcings})
            scenario_list.append(ds_scenario)                                         
        ds_model = xr.concat(scenario_list, dim='scenario').assign_coords({'scenario':scenarios})
        model_list.append(ds_model)

    ds_area_exposed_full = xr.concat(model_list,dim='model').assign_coords({'model':models})
    return ds_area_exposed_full.to_dataset(name='area_exposed')

# load waterscarcity exposure and calculate the total and relative population exposed to water scarcity per forcing, scenario and model
def calc_totalpop_exposed(variable, landmask, pop, models, scenarios, forcings):

    # load total land area exposed
    total_pop = pop.where(landmask).sum(dim=('lat','lon'))

    # loop over models
    model_list = []
    model_list_rel = []
    for model in models:    
        scenario_list = []
        scenario_list_rel = []
        for scenario in scenarios:
            ds_forcings = xr.Dataset()
            # loop over GCM forcings

            pop_exposed_list = []
            pop_exposed_rel_list = []

            for forcing in forcings: 
                print(model +' '+scenario+' '+forcing)

                # define filename and output directory of model

                filedir, filename = get_simulation_name_dir(variable, lifetimedir+variable, model, forcing, scenario, timestep, socscen=False)

                if filename!='false': 
                    # save to output location
                    ds_waterscarcity_exposure = xr.open_dataset(filedir+filename, engine='netcdf4')
                    da_waterscarcity_exposure = ds_waterscarcity_exposure['exposure'].where(landmask)
                    # as pop is annual and exposure likely monthlx, convert time axis to datetime to make them comparable to each other
                    da_waterscarcity_exposure['time'] = da_waterscarcity_exposure.indexes['time'].to_datetimeindex()
                    
                    pop_exposed = (da_waterscarcity_exposure * pop).sum(dim=('lat','lon'))
                    pop_exposed_rel = (pop_exposed / total_pop) * 100

                else: 
                    if scenario == 'historical': 
                        empty_array = np.empty_like(pop_exposed_list[0].values)
                        empty_array_rel = np.empty_like(pop_exposed_rellist[0].values)

                    else: 
                        empty_array = np.empty_like(pop_exposed_list[1].values)
                        empty_array_rel = np.empty_like(pop_exposed_list_rel[1].values)

                    empty_array[:] = np.nan
                    empty_array_rel[:] = np.nan    

                    pop_exposed = xr.DataArray(empty_array, dims=('time'))
                    pop_exposed_rel = xr.DataArray(empty_array_rel, dims=('time'))

                pop_exposed_list.append(pop_exposed)
                pop_exposed_rel_list.append(pop_exposed_rel)

            ds_scenario = xr.concat(pop_exposed_list, dim='forcing').assign_coords({'forcing':forcings})
            ds_scenario_rel = xr.concat(pop_exposed_rel_list, dim='forcing').assign_coords({'forcing':forcings})

            scenario_list.append(ds_scenario)      
            scenario_list_rel.append(ds_scenario_rel)                                         

        ds_model = xr.concat(scenario_list, dim='scenario').assign_coords({'scenario':scenarios})
        model_list.append(ds_model)

        ds_model_rel = xr.concat(scenario_list_rel, dim='scenario').assign_coords({'scenario':scenarios})
        model_list_rel.append(ds_model_rel)

    ds_pop_exposed_full = xr.concat(model_list,dim='model').assign_coords({'model':models})
    ds_pop_exposed_rel_full = xr.concat(model_list_rel,dim='model').assign_coords({'model':models})

    return ds_pop_exposed_full.to_dataset(name='pop_exposed'), ds_pop_exposed_rel_full.to_dataset(name='pop_exposed')


def plot_global_land_exposure_mmm(variable, landmask, da_cellarea, models, scenarios, forcings): 

    # calculate land area exposed
    ds_area_exposed_full = calc_totalarea_exposed(variable, landmask, da_cellarea, models, scenarios, forcings)

    # per model and scenario
    fig, ax = plt.subplots(figsize=(12,6))

    # loop over models
    for model in models:

        for scenario in scenarios:

            ds_area_exposed_full['area_exposed'].sel({'model': model, 'scenario':scenario}).mean(dim='forcing').plot(ax=ax, label = model+' ' +scenario)
    ax.legend()
    ax.set_ylabel('% of land area')
    ax.set_title('')
    ax.set_title('Fraction of land area annually exposed to water scarcity '+variable_name[variable], loc='right');
    ax.grid(color='lightgray', linestyle='-', linewidth=1, alpha=0.5)

def plot_global_land_exposure_permodelscenario(variable, landmask, da_cellarea, models, scenarios, forcings): 
   
    # calculate land area exposed
    ds_area_exposed_full = calc_totalarea_exposed(variable, landmask, da_cellarea, models, scenarios, forcings)

    # per model and scenario
    fig, ax = plt.subplots(figsize=(12,6))


    for scenario in scenarios:

        mean = ds_area_exposed_full['area_exposed'].sel({'scenario':scenario}).mean(dim=('forcing','model'))
        std = ds_area_exposed_full['area_exposed'].sel({'scenario':scenario}).std(dim=('forcing','model'))
        mean.plot(ax=ax, label =scenario)
        ax.fill_between(mean.time.values, mean-std,mean+std, alpha=0.2)

    ax.legend()
    ax.set_ylabel('% of land area')
    ax.set_title('')
    ax.set_title('Fraction of land area annually exposed to water scarcity '+variable_name[variable], loc='right');
    ax.grid(color='lightgray', linestyle='-', linewidth=1, alpha=0.5)

    
# load waterscarcity exposure in dataset format per forcing, scenario and model
def load_exposure_ds(variable, models, scenarios, forcings):


    # loop over models
    model_list = []
    for model in models:    
        scenario_list = []
        for scenario in scenarios:
            ds_forcings = xr.Dataset()
            # loop over GCM forcings

            da_list = []
            for forcing in forcings: 

                # define filename and output directory of model
                
                filedir, filename = get_simulation_name_dir(variable, lifetimedir+variable, model, forcing, scenario, 'annual', socscen=False)

                # deal with files not present - fill with nans
                if filename!='false': 
                    # save to output location
                    ds_waterscarcity_exposure = xr.open_dataset(filedir+filename, engine='netcdf4')
                    da_waterscarcity_exposure = ds_waterscarcity_exposure['exposure']

                else: 
                    if scenario == 'historical': 
                        empty_array = np.empty_like(da_list[0].values)

                    else: 
                        empty_array = np.empty_like(da_list[1].values)
                    
                    empty_array[:] = np.nan    
                    da_waterscarcity_exposure = xr.DataArray(empty_array, dims=('time'))
                    
                da_list.append(da_waterscarcity_exposure)


            ds_scenario = xr.concat(da_list, dim='forcing').assign_coords({'forcing':forcings})
            scenario_list.append(ds_scenario)                                         
        ds_model = xr.concat(scenario_list, dim='scenario').assign_coords({'scenario':scenarios})
        model_list.append(ds_model)

    ds_exposure = xr.concat(model_list,dim='model').assign_coords({'model':models})
    return ds_exposure


# calculate the value of reference GMT to calculate anomalies
def get_GMT_ref(file_names_gmt, year_start_GMT_ref, year_end_GMT_ref): 

    if phase =='2b':
        
        picontrol_name = '_piControl_'
    else: 
        picontrol_name = '_picontrol_'
    file_name_gmt = [s for s in file_names_gmt if '_historical_' in s]
    file_name_gmt_pic = [s for s in file_names_gmt if picontrol_name in s]

    GMT_his = pd.read_csv(
        file_name_gmt[0],
        delim_whitespace=True,
        skiprows=1,
        header=None).rename(columns={0:'year',1:'tas'}).set_index('year')

    GMT_pic = pd.read_csv(
        file_name_gmt_pic[0],
        delim_whitespace=True, 
        skiprows=1, 
        header=None).rename(columns={0:'year',1:'tas'}).set_index('year')

    GMT_ref = pd.concat([GMT_pic.loc[year_start_GMT_ref:np.min(GMT_his.index)-1,:], GMT_his.loc[:year_end_GMT_ref,:]]).mean().values[0]
    
    return GMT_ref


# get the GMT anomalies based on the forcing and scenario
def get_GMT_anomalies(forcing, scenario):
    
    if phase =='2b':
        
        if forcing == 'hadgem2-es': # .upper() method doesn't work for HadGEM2-ES on linux server (only Windows works here)
            forcing_dirname ='HadGEM2-ES' # ignore running mean files
        else:
            forcing_dirname = forcing.upper()

        file_names_gmt = glob.glob(datadir+'/DerivedInputData/globalmeans/GCM_atmosphere_biascorrected_tas/'+forcing_dirname+'/*tas*fldmean.yearmean.txt') # ignore running mean files
    
    else: # 3b
        # load GMT for rcp and historical period - note that these data are in different files

        if forcing == 'hadgem2-es': # .upper() method doesn't work for HadGEM2-ES on linux server (only Windows works here)
            forcing_dirname ='HadGEM2-ES' # ignore running mean files
        elif forcing == 'canesm5': # .upper() method doesn't work for HadGEM2-ES on linux server (only Windows works here)
            forcing_dirname = 'CanESM5'
        elif forcing == 'ec-earth3': # .upper() method doesn't work for HadGEM2-ES on linux server (only Windows works here)
            forcing_dirname = 'EC-Earth3' # ignore running mean files
        else:
            forcing_dirname = forcing.upper()


        file_names_gmt = glob.glob(lifetimedir+'/DerivedInputData/globalmeans/tas/'+forcing_dirname+'/*tas_global_annual*.txt') # ignore running mean files
    
    file_name_gmt = [s for s in file_names_gmt if scenario in s]

    GMT = pd.read_csv(
        file_name_gmt[0],
        delim_whitespace=True,
        skiprows=1,
        header=None).rename(columns={0:'year',1:'tas'}).set_index('year')


    # concatenate historical and future data
    #df_GMT = pd.concat([GMT_his,GMT_fut])

    # calculate reference GMT - use data from pic until 1861 and from his from then onwards
    GMT_ref = get_GMT_ref(file_names_gmt, year_start_GMT_ref, year_end_GMT_ref)

    # convert GMT from absolute values to anomalies
    df_GMT = GMT - GMT_ref

    return df_GMT

# get a dataset of GMT anomalies per model, forcing and scenario
def get_GMT_xr(models, scenarios, forcings):
    # GMT xr
    model_list = []
    for model in models:
        scenario_list = []
        for scenario in scenarios: 
            forcing_list = []
            for forcing in forcings: 

                df_GMT = get_GMT_anomalies(forcing, scenario)
                forcing_list.append(df_GMT.to_xarray())

            ds_scenario = xr.concat(forcing_list, dim='forcing').assign_coords({'forcing':forcings})
            scenario_list.append(ds_scenario)
        ds_model = xr.concat(scenario_list, dim='scenario').assign_coords({'scenario':scenarios})
        model_list.append(ds_model)
    ds_GMT = xr.concat(model_list, dim='model').assign_coords({'model':models}).rename({'year':'time'})
    
    return ds_GMT

# get a dataset of GMT anomalies per model, forcing and scenario, only for future scenarios (hist concatenated)
def get_GMT_concat_xr(models, scenarios, forcings):
    # GMT xr
    model_list = []
    for model in models:

        forcing_list = []
        for forcing in forcings: 

            scenario_list = []
            for scenario in scenarios[1:]: 

                df_GMT_fut = get_GMT_anomalies(forcing, scenario)
                df_GMT_hist = get_GMT_anomalies(forcing, 'historical')
                df_GMT = pd.concat([df_GMT_hist, df_GMT_fut])

                scenario_list.append(df_GMT.to_xarray())

            ds_scenario = xr.concat(scenario_list, dim='scenario').assign_coords({'scenario':scenarios[1:]})
            forcing_list.append(ds_scenario)
        ds_model = xr.concat(forcing_list, dim='forcing').assign_coords({'forcing':forcings})
        model_list.append(ds_model)
    ds_GMT = xr.concat(model_list, dim='model').assign_coords({'model':models}).rename({'year':'time'})
    
    return ds_GMT


def plot_exposure_GMT_scaling(variable, landmask, da_cellarea, models, scenarios, forcings, timestep='monthly'):
    # load area exposed dataset and GMT
    ds_area_exposed = calc_totalarea_exposed(variable, landmask, da_cellarea, models, scenarios, forcings)
    ds_GMT = get_GMT_xr(models, scenarios, forcings)

    if variable == 'falkenmark': 
        # account for differing period
        ds_GMT = ds_GMT.where(ds_GMT['time'] >= ds_area_exposed.time[0], drop=True)

    # add GMT to area exposed dataset
    ds_GMT.coords['time'] = ds_area_exposed.time
    ds_area_exposed['GMT'] = ds_GMT['tas']

    ds_area_exposed_mmm = ds_area_exposed.mean(dim=('forcing','model'))
    ds_area_exposed_min = ds_area_exposed.min(dim=('forcing','model'))
    ds_area_exposed_max = ds_area_exposed.max(dim=('forcing','model'))

    # do plotting
    fig, ax = plt.subplots(figsize=(12,6))

    for scenario in scenarios:
        ds_sel = ds_area_exposed_mmm.sel({'scenario':scenario})
        ds_sel_min = ds_area_exposed_min.sel({'scenario':scenario})
        ds_sel_max = ds_area_exposed_max.sel({'scenario':scenario})

        ax.plot(ds_sel['GMT'],ds_sel['area_exposed'], label=scenario)
        ax.fill_between(ds_sel['GMT'], ds_sel_min['area_exposed'],ds_sel_max['area_exposed'], alpha=0.2)
        ax.legend(frameon=False)
        ax.set_title('GMT scaling of '+variable_name[variable], loc='right')
        ax.set_ylabel('Global land area exposed to water scarcity (%)')
        ax.set_xlabel('$\Delta$ GMT (K)')
        ax.grid(color='lightgray', linestyle='-', linewidth=1, alpha=0.5)
        
# load waterscarcity exposure and calculate the total area exposed to water scarcity per forcing, scenario and model
# load waterscarcity exposure and calculate the total area exposed to water scarcity per forcing, scenario and model
# load waterscarcity exposure and calculate the total area exposed to water scarcity per forcing, scenario and model
# load waterscarcity exposure and calculate the total area exposed to water scarcity per forcing, scenario and model
def get_exposure_map_GMT_level(variable, GMT_level, GMT_offset, landmask, models, scenarios, forcings, flag_intermediate):


    # load GMTs
    ds_GMT = get_GMT_concat_xr(models, scenarios, forcings)
    
    # look up time of GMT level reached
    time_mask = ((ds_GMT>GMT_level-GMT_offset) & (ds_GMT<GMT_level+GMT_offset))['tas']
    ds_lookup = ds_GMT['tas'].where(time_mask)>-5
    
    # loop over models
    model_list = []
    for model in models:    
        forcing_list = []

        for forcing in forcings: 


            map_list = []
            scenarios_used = []

            for scenario in scenarios[1:]:
                if scenario in ds_lookup['scenario'].values: 

                    scenarios_used.append(scenario)

                    # define filename and output directory of model
                    if flag_intermediate: 
                        var_dir = outdir+'intermediate/'
                        da_var = variable
                        socscen = True
                    else: 
                        var_dir = lifetimedir+variable
                        if phase == '2b':
                            var_dir = lifetimedir+ '2b/'+variable                    
                        da_var = 'exposure'
                        socscen=False

                    filedir, filename = get_simulation_name_dir(variable, var_dir, model, forcing, scenario, 'annual', socscen=False)
                    filedir_hist, filename_hist = get_simulation_name_dir(variable, var_dir, model, forcing, 'historical', 'annual', socscen=False)

                    if filename!='false': 
                        # save to output location
                        ds_waterscarcity_exposure_fut = xr.open_dataset(filedir+filename, engine='netcdf4')
                        ds_waterscarcity_exposure_hist = xr.open_dataset(filedir_hist+filename_hist, engine='netcdf4')
                        ds_waterscarcity_exposure = xr.concat([ds_waterscarcity_exposure_fut,ds_waterscarcity_exposure_hist], dim='time')
                        # select years to calculate average from corresponding to GMT average

                        selected_years = ds_lookup.time.where(ds_lookup.sel({'forcing':forcing,'scenario':scenario,'model':model}), drop=True).values
                        
                        if variable == 'falkenmark': 
                            years_mask = ds_waterscarcity_exposure['time'].isin(selected_years) 
                        else: 
                            years_mask = ds_waterscarcity_exposure['time'].dt.year.isin(selected_years) 
                            
                        da_map = ds_waterscarcity_exposure[da_var].where(years_mask).mean('time')*100
                        
                    else: 
                        if scenario == 'historical': 
                            empty_array = np.empty_like(map_list[0].values)

                        else: 
                            empty_array = np.empty_like(map_list[1].values)

                        empty_array[:] = np.nan    
                        da_map = xr.DataArray(empty_array, dims=('time'))

                    map_list.append(da_map)

                            
            ds_scenario = xr.concat(map_list, dim='scenario').assign_coords({'scenario':scenarios_used})
            forcing_list.append(ds_scenario)    
                                
        ds_model = xr.concat(forcing_list, dim='forcing') .assign_coords({'forcing':forcings})
            
        model_list.append(ds_model)

    da_map_full = xr.concat(model_list,dim='model').assign_coords({'model':models}).where(landmask)
    
    return da_map_full

def calc_mean_perbasin(da, landmask):
    print('calculating means per basin')
    #da_basins = xr.open_dataset(outdir+'intermediate/basins.nc')['basin']
    da_basins = xr.open_dataset(outdir+'upstream_calc/basins.nc')['basin']
    
    # only keep basins with more than 10 grid cells to reduce time
    basins, basin_counts = np.unique(da_basins, return_counts=True)
    basin_numbers = basins[basin_counts>10]
    da_basin_means = []

    for number in basin_numbers: 
        mean_basin = da.where(da_basins==number).mean(dim=('lat','lon','forcing','scenario','model'))
        da_basin_means.append(((da_basins==number) * mean_basin.values.item()))

    da_basins_concat = xr.concat(da_basin_means, dim='basins')
    da_perbasin = da_basins_concat.sum(dim='basins').where(landmask)

    return da_perbasin

def calc_exposure_map(variable, GMT_offset, GMT_level, landmask, models, scenarios, forcings, flag_intermediate=False, flag_perbasin=False, flag_relative = False): 
    # calculate exposure map for 0° GMT warming
    print('calculating ref')
    da_map_ref = get_exposure_map_GMT_level(variable, 0, GMT_offset, landmask, models, scenarios, forcings, flag_intermediate) #['area_exposed']

    
    # calculate exposure map for 0° GMT warming
    print('calculating fut')
    da_map_fut = get_exposure_map_GMT_level(variable, GMT_level, GMT_offset, landmask, models, scenarios, forcings, flag_intermediate) #['area_exposed']

    # calculate difference and multi-model mean
    delta_map = (da_map_fut - da_map_ref)
    
    if flag_relative: 
        delta_map = delta_map/da_map_ref *100

    if flag_perbasin: 
        delta_mmm = calc_mean_perbasin(delta_map, landmask)
    else: 
        delta_mmm = delta_map.mean(dim=('forcing','model','scenario'))

    return delta_mmm


# function to plot exposure map at predefined GMT levels
def plot_exposure_map(variable, GMT_offset, GMT_level, landmask, models, scenarios, forcings, flag_intermediate=False, flag_perbasin=False): 
    
    delta_mmm = calc_exposure_map(variable, GMT_offset, GMT_level, landmask, models, scenarios, forcings, flag_intermediate=True, flag_perbasin=False)

        
    # figure plotting
    fig = plt.figure(figsize=(12,5))
    ax = plt.subplot(111)


    if not flag_intermediate: 
        delta_mmm.plot(ax=ax, cbar_kwargs={'label': '$\Delta$ land area exposed (%)', 'fraction': 0.02, 'pad': 0.04})
        ax.set_ylabel('')
        ax.set_xlabel('')
        ax.set_title('$\Delta$ land area exposed at +'+str(GMT_level)+'° GMT for the '+variable_name[variable], loc='right')

    else: 
        delta_mmm.plot(ax=ax, cbar_kwargs={'label': '$\Delta$ water availability mm/s', 'fraction': 0.02, 'pad': 0.04}, cmap='RdBu') # vmax=0.05, vmin = -0.05,
        ax.set_ylabel('')
        ax.set_xlabel('')
        ax.set_title('$\Delta$ water availability at +'+str(GMT_level)+'° GMT ', loc='right')


        
        
# load intermediate variable per forcing, scenario and model
def load_intermediate_variable_ds(variable, mask,  models, scenarios, forcings):
   

    if phase == '3b': 
        decode_times = True 
    else: 
        decode_times = False
        
        
    # loop over models
    model_list = []
    for model in models:
        print('Loading model '+model)
        scenario_list = []
        for scenario in scenarios:
            ds_forcings = xr.Dataset()
            # loop over GCM forcings

            area_exposed_list = []
            for forcing in forcings: 

                # define filename and output directory of model

                filedir, filename = get_simulation_name_dir(variable, outdir+'intermediate/', model, forcing, scenario, 'monthly')
                if filename!='false': 
                    # save to output location
                    ds = xr.open_dataset(filedir+filename, engine='netcdf4', decode_times=decode_times)[variable]
                    
                    # for isimip 2b, manually assign time periods, as these can not be read in 
                    if phase=='2b': 
                        ds['time'] = d_time_ncfile[scenario]
                        
                    ds_var = ds.groupby('time.year').mean()
                    
                    ds_var_tosave = ds_var.where(mask).mean(dim=('lat','lon'))
                    
                    #select one grid cell 
                    #    ds_var_tosave = ds_var.sel({'lat':mode[0], 'lon':mode[1]})
                    del ds_var
                else: 
                    if scenario == 'historical': 
                        empty_array = np.empty_like(area_exposed_list[0].values)

                    else: 
                        empty_array = np.empty_like(area_exposed_list[1].values)
                    
                    empty_array[:] = np.nan    
                    ds_var_tosave = xr.DataArray(empty_array, dims=('time'))
                
                area_exposed_list.append(ds_var_tosave)


            ds_scenario = xr.concat(area_exposed_list, dim='forcing').assign_coords({'forcing':forcings})
            scenario_list.append(ds_scenario)                                         
        ds_model = xr.concat(scenario_list, dim='scenario').assign_coords({'scenario':scenarios})
        model_list.append(ds_model)

    ds_var_full = xr.concat(model_list,dim='model').assign_coords({'model':models})
    return ds_var_full


### plot the intermediate variable in time series per model and scenario

def plot_intermediate_variable(variable, mask, models, scenarios, forcings): 
    
    ds = load_intermediate_variable_ds(variable, mask,  models, scenarios, forcings)
    

    title_scale = 'Mean '+str(mask['names'].values)+' '


    if variable == 'q_avail': 
        title_var = 'water availability '
    elif variable == 'ptotww':
        title_var = 'water withdrawal '
    elif variable == 'waterscarcity_index':
        title_var = 'water scarcity index '
        
    var = 'ws_index'
    fig, ax = plt.subplots(figsize=(12,6))

    # loop over models
    for i, scenario in enumerate(scenarios):
            da_sel = ds.sel({'scenario':scenario})

            da_sel.mean(dim=('forcing','model')).plot(ax=ax, label =scenario)
            ax.fill_between(da_sel.year.values,da_sel.min(dim=('forcing','model')),da_sel.max(dim=('forcing','model')), alpha=0.3)

    ax.legend()
    ax.set_ylabel('mm/m²')
    ax.set_title('')
    ax.set_title(title_scale+title_var +' (multi-model mean)', loc='right');


    fig, ax = plt.subplots(figsize=(12,6))

    for model in models: 
        # loop over models
            da_sel = ds.sel({'model':model})

            da_sel.mean(dim=('forcing','scenario')).plot(ax=ax, label = model)
            ax.fill_between(da_sel.year.values,da_sel.min(dim=('forcing','scenario')),da_sel.max(dim=('forcing','scenario')), alpha=0.3)

    ax.legend()
    ax.set_ylabel('mm/m²')
    ax.set_title('')
    ax.set_title(title_scale+title_var +' (per model)', loc='right');
    

    
    
# save basins file per numbers
def save_basinfile(): 
    # -----------------

    basinfile = routingdir+"ddm30_basins_cru_neva.nc"

    # ---------------------------------------
    # load datasets

    # Basins from DDM30 - ISIMIP2
    basins = np.array(np.zeros((360,720)), dtype=np.int64)
    nf1 = Dataset(basinfile, 'r')
    prefix = list(nf1.variables.items())[-1][0]
    basin=nf1.variables[prefix][:,:].data
    basin = basin.astype(np.int64)
    basin[basin>10e15] = 0
    basin[basin<0] = 0
    basin = np.flipud(basin)
    #basin = basin.astype(np.int64)
    basins[12:292,:]= basin
    nf1.close()

    # open land mask
    landmask  = xr.open_dataarray(routingdir +'ddm30_basins_cru_neva.nc') >0

    da_basins = xr.DataArray(np.flipud(basin), dims=( 'lat','lon'), coords = { 'lat':landmask.lat, 'lon':landmask.lon})
    da_basins = da_basins.where(da_basins>0).where(da_basins<70000)

    da_basins = (da_basins/1000).round()
    da_basins.to_dataset(name='basin').to_netcdf(outdir+'intermediate/basins.nc')
    

    

# ---------------------------------------------------------------------
# 3. Functions to plot
# ---------------------------------------------------------------------

def set_plot_param():
    """Set my own customized plotting parameters"""
    
    import matplotlib as mpl
    mpl.rc('axes',edgecolor='grey')
    mpl.rc('axes',labelcolor='dimgrey')
    mpl.rc('xtick',color='dimgrey')
    mpl.rc('xtick',labelsize=12)
    mpl.rc('ytick',color='dimgrey')
    mpl.rc('ytick',labelsize=12)
    mpl.rc('axes',titlesize=14)
    mpl.rc('axes',labelsize=12)
    mpl.rc('legend',fontsize='large')
    mpl.rc('text',color='dimgrey')


def plot_ts(da):
    """ plot timeseries (of spatial mean) """
    da_ts = da.mean(dim=('lon','lat'))
    da_ts.plot()
    plt.title(da.long_name + '('+da.units+')')
    plt.xlim(da_ts.time[0].values,da_ts.time[-1].values)


def plot_ymean(da):
    """calculate and plot annual mean timeseries """
    # check if input da is already ymean, otherwise do calculation 
    if len(da) < 500: 
        da_ymean = da.mean(dim=('lon','lat'))
    else: 
        da_ymean = da.mean(dim=('lon','lat')).groupby('time.year').mean('time')
    
    xlims = (da_ymean.year[0].values,da_ymean.year[-1].values)
    da_tseries = da_ymean.plot(xlim=xlims)
    plt.title(da.long_name, pad=5)
    plt.ylabel(da.name+' [' + da.units + ']')
    #plt.plot([da_ymean.year[0],da_ymean.year[-1]], [0,0], linewidth=1, color='gray')

# calculate and plot annual sum timeseries 
def plot_ysum(da):
    
    da_ymean = da.sum(dim=('lon','lat')).groupby('time.year').mean('time')
    xlims = (da_ymean.year[0].values,da_ymean.year[-1].values)
    da_tseries = da_ymean.plot(xlim=xlims)
    plt.title(da.long_name+' [' + da.units + ']' )


# plot timmean per selected region
def plot_yts_sel_regions(da_to_mask, selected_regions):
    mask = regionmask.defined_regions.srex.mask(da_to_mask)
    
    # annual means are already calculated
    if len(da_to_mask) < 50:
        da_mask_ts = da_to_mask.groupby(mask).mean('stacked_lat_lon').sel(region=selected_regions)
    else: 
        da_mask_ts = da_to_mask.groupby(mask).mean('stacked_lat_lon').groupby('time.year').mean().sel(region=selected_regions)

    # add abbreviations and names
    abbrevs = regionmask.defined_regions.srex[da_mask_ts.region.values].abbrevs
    names = regionmask.defined_regions.srex[da_mask_ts.region.values].names
    da_mask_ts.coords['abbrevs'] = ('region', abbrevs)
    da_mask_ts.coords['names'] = ('region', names)
    
    f, axes = plt.subplots(3, 2, figsize=(8,5))
    f.suptitle(da_to_mask.name, fontsize=14)

    low = da_mask_ts.min()
    high = da_mask_ts.max()
    for i in range(len(selected_regions)):
        (nx,ny) = axes.shape
        if i < nx : ax = axes[i,0]
        else      : ax = axes[i-nx,1]
            
        ts_region = da_mask_ts.isel(region=i)
        ts_region.plot(ax=ax)
        ax.set_title(da_mask_ts.isel(region=i).names.values)
        ax.set_ylim(low,high)
        ax.set_ylabel('('+da_to_mask.units+')')
        ax.set_xlim(da_mask_ts.year[0],da_mask_ts.year[-1])
        ax.plot([da_mask_ts.year[0],da_mask_ts.year[-1]], [0,0], linewidth=1, color='gray')
    
    plt.setp(axes, xlabel="")

    f.tight_layout()


# plot global map of difference 
def plot_delta_map(da_delta, plot_regions=False, vlims=False, calcsum=False, cmap='BrBG'):
    
    # calculate annual sum instead of mean (precip)
    if calcsum: 
        da_delta_ysum = da_delta.groupby('time.year').sum()
        da_delta_mean = da_delta_ysum.mean('year')
        da_delta_mean.attrs['units'] = 'mm/year'
    # only one value
    elif len(da_delta.dims) < 3: 
        da_delta_mean = da_delta
    # annual means already taken
    elif len(da_delta) < 50:
        da_delta_mean = da_delta.mean('year')
    else:
        da_delta_mean = da_delta.mean('time')
    
    fig = plt.figure(figsize=(12,5))
    ax = plt.subplot(111)
        
    # limiting values for plotting are given    
    if vlims==False: 
        da_delta_mean.plot(ax=ax, cmap=cmap, cbar_kwargs={'label': da_delta.name+' ('+da_delta.units+')', 'fraction': 0.02, 'pad': 0.04})
    else: 
        da_delta_mean.plot(ax=ax, cmap=cmap, vmin=vlims[0], vmax=vlims[1], extend='both',  cbar_kwargs={'label': da_delta.name+' ('+da_delta.units+')', 'fraction': 0.02, 'pad': 0.04}, add_labels=False)
        
    ax.set_title(da_delta.long_name, loc='right')
    ax.coastlines(color='dimgray', linewidth=0.5)
    # exclude Antactica from plot
    #ax.set_extent((-180,180,-63,90), crs=proj) 

    if plot_regions: regionmask.defined_regions.srex.plot(ax=ax,add_ocean=False, coastlines=False, add_label=False) #label='abbrev'
    return fig, ax



f, ax = plt.subplots()

# plot global map of difference 
def plot_delta_map_noax(ax, da_delta, plot_regions=False, vlims=False, calcsum=False, cmap='BrBG'):
    """plot difference maps without creating a figure within function"""
    # calculate annual sum instead of mean (precip)
    if calcsum: 
        da_delta_ysum = da_delta.groupby('time.year').sum()
        da_delta_mean = da_delta_ysum.mean('year')
        da_delta_mean.attrs['units'] = 'mm/year'
    # only one value
    elif len(da_delta.dims) < 3: 
        da_delta_mean = da_delta
    # annual means already taken
    elif len(da_delta) < 50:
        da_delta_mean = da_delta.mean('year')
    else:
        da_delta_mean = da_delta.mean('time')
    
    # limiting values for plotting are given    
    if vlims==False: 
        da_delta_mean.plot(ax=ax, cmap=cmap, cbar_kwargs={'label': da_delta.name+' ('+da_delta.units+')'})
    else: 
        da_delta_mean.plot(ax=ax, cmap=cmap, vmin=vlims[0], vmax=vlims[1], extend='both',  cbar_kwargs={'label': da_delta.name+' ('+da_delta.units+')'}, add_labels=False)
        
    ax.set_title(da_delta.long_name, loc='right')
    ax.coastlines(color='dimgray', linewidth=0.5)
    # exclude Antactica from plot
    ax.set_extent((-180,180,-63,90)) 

    if plot_regions: regionmask.defined_regions.srex.plot(ax=ax,add_ocean=False, coastlines=False, add_label=False) #label='abbrev'
    return ax


# save new variable to lifetimedirectory in right format
def save_to_lifetimedir(da, variable, model, forcing, scenario, annual=False):
    # save water scarcity index
    filedir, filename = get_simulation_name_dir('waterscarcity', lifetimedir+'waterscarcity', model, forcing, scenario, timestep, socscen=False)    
    save_dir =  lifetimedir+variable+'/'+model.lower()+'/'
    
    # Create a new directory because it does not exist 
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    filename = save_dir+filename.replace('waterscarcity', variable).replace('1850','1861')
    
    if annual: 
        filename = filename.replace('monthly','annual')
        
    if os.path.exists(filename): 
        os.system('rm '+filename)
        
        
    da.to_dataset(name=variable).to_netcdf(filename)


# calculate and save water deficit, and water withdrawal
def calc_water_deficit(model, forcing, scenario, da_cellarea, landmask, pop):

    # load grid cell runoff
    qtot = load_simulation('qtot', datadir, model, forcing, scenario)

    # load discharge into gricell from directly upstream grid cells
    dis_upstream = load_simulation('dis_directupstream', datadir, model, forcing, scenario, flag_postprocessed=True)

    # calculate water availability (q_avail)
    q_avail_withoutefr = dis_upstream * 1000/da_cellarea + qtot # UNITS! dis is in m³/s, convert to mm/m²s

    # subtract envrironmental flow requirement of water availability
    q_avail =  q_avail_withoutefr - get_efr(q_avail_withoutefr)
    del  dis_upstream
    
    # put hard condition of setting negative water availabilities to 0 (this is a problem for the simulations of H08)
    q_avail = q_avail.where(q_avail>0,0)

    # load water withdrawal data (demand)
    ptotww = load_simulation('ptotww', datadir, model, forcing, scenario, flag_postprocessed=True)

    # solve this more ellegantly

    # fix the time dimension. 
    qtot = qtot.assign_coords(time = d_time_ncfile[scenario])
 
    days_in_month = []
    for q_avail_month in qtot.time.values:
        days_in_month.append(pd.to_datetime(q_avail_month).days_in_month)

    del qtot
    
    secs_in_month = np.expand_dims(np.array(days_in_month) * 24*60*60, axis=(1,2))
    
    pop = pop.where(landmask) 

    # hardcoded bit to get time dimension right
    if phase == '3b': 
        if scenario == 'historical': 
            time_ds = xr.open_dataset(lifetimedir+'falkenmark/cwatm/cwatm_gfdl-esm4_historical_falkenmark_global_annual_landarea_1861_2014.nc4')['exposure'].time
        else: 
            time_ds = xr.open_dataset(lifetimedir+'falkenmark/cwatm/cwatm_gfdl-esm4_ssp126_falkenmark_global_annual_landarea_2015_2100.nc4')['exposure'].time
    elif phase == '2b': 
        if scenario == 'historical': 
            time_ds = xr.open_dataset(lifetimedir+'falkenmark/lpjml/lpjml_gfdl-esm2m_historical_falkenmark_global_annual_landarea_1861_2005.nc4')['exposure'].time
        else: 
            time_ds = xr.open_dataset(lifetimedir+'falkenmark/lpjml/lpjml_gfdl-esm2m_rcp26_falkenmark_global_annual_landarea_2006_2099.nc4')['exposure'].time

    # calculate water deficit
    water_deficit = ((ptotww - q_avail) * da_cellarea/1000  * secs_in_month)  #m³
    
    del q_avail    
    withdrawal = (ptotww*da_cellarea/1000 * secs_in_month).groupby('time.year').sum().rename({'year':'time'})

    water_deficit_annual = water_deficit.where(water_deficit>0).groupby('time.year').sum().rename({'year':'time'})
    
    if phase == '2b': 
        # fix the time dimension. 
        if scenario == 'historical': 
            water_deficit_annual = water_deficit_annual.assign_coords(time = pd.date_range(start='1861', end='2006', freq='Y'))
            withdrawal  = withdrawal.assign_coords(time = pd.date_range(start='1861', end='2006', freq='Y'))
        else: 
            water_deficit_annual = water_deficit_annual.assign_coords(time = pd.date_range(start='2006', end='2100', freq='Y'))
            withdrawal  = withdrawal.assign_coords(time = pd.date_range(start='2006', end='2100', freq='Y'))
    elif phase == '3b': 
        # fix the time dimension. 
        if scenario == 'historical': 
            water_deficit_annual = water_deficit_annual.assign_coords(time = pd.date_range(start='1850', end='2015', freq='Y'))
            withdrawal  = withdrawal.assign_coords(time = pd.date_range(start='1850', end='2015', freq='Y'))
        else: 
            water_deficit_annual = water_deficit_annual.assign_coords(time = pd.date_range(start='2015', end='2101', freq='Y'))
            withdrawal  = withdrawal.assign_coords(time = pd.date_range(start='2015', end='2101', freq='Y'))
            
            
    water_deficit_perperson = (water_deficit_annual.where(pop)/pop).assign_coords(time= time_ds.values) 

    withdrawal_perperson = (withdrawal.where(pop) / pop).assign_coords(time= time_ds.values)


    # save into new variables
    save_to_lifetimedir(water_deficit_perperson, 'waterdeficit', model, forcing, scenario, annual=True)
    save_to_lifetimedir(withdrawal_perperson, 'withdrawal', model, forcing, scenario, annual = True)    
    
    
    
# total water withdrawal for evaluation purposes
def calc_total_water_withdrawal(model, forcing, scenario, da_cellarea, landmask, pop):

    # load water withdrawal data (demand)
    ptotww = load_simulation('ptotww', datadir, model, forcing, scenario, flag_postprocessed=True)
    
    days_in_month = []
    for q_avail_month in ptotww.time.values:
        days_in_month.append(pd.to_datetime(q_avail_month).days_in_month)

    
    secs_in_month = np.expand_dims(np.array(days_in_month) * 24*60*60, axis=(1,2))
    
    # hardcoded bit to get time dimension right
    if phase == '3b': 
        if scenario == 'historical': 
            time_ds = xr.open_dataset(lifetimedir+'falkenmark/cwatm/cwatm_gfdl-esm4_historical_falkenmark_global_annual_landarea_1861_2014.nc4')['exposure'].time
        else: 
            time_ds = xr.open_dataset(lifetimedir+'falkenmark/cwatm/cwatm_gfdl-esm4_ssp126_falkenmark_global_annual_landarea_2015_2100.nc4')['exposure'].time
    elif phase == '2b': 
        if scenario == 'historical': 
            time_ds = xr.open_dataset(lifetimedir+'falkenmark/lpjml/lpjml_gfdl-esm2m_historical_falkenmark_global_annual_landarea_1861_2005.nc4')['exposure'].time
        else: 
            time_ds = xr.open_dataset(lifetimedir+'falkenmark/lpjml/lpjml_gfdl-esm2m_rcp26_falkenmark_global_annual_landarea_2006_2099.nc4')['exposure'].time

    withdrawal = (ptotww*da_cellarea/1000 * secs_in_month).groupby('time.year').sum().rename({'year':'time'})

    
    if phase == '2b': 
        # fix the time dimension. 
        if scenario == 'historical': 
            withdrawal  = withdrawal.assign_coords(time = pd.date_range(start='1861', end='2006', freq='Y'))
        else: 
            withdrawal  = withdrawal.assign_coords(time = pd.date_range(start='2006', end='2100', freq='Y'))

    withdrawal_total = withdrawal.assign_coords(time= time_ds.values)


    # save into new variables
    save_to_lifetimedir(withdrawal_total, 'totalwithdrawal', model, forcing, scenario, annual = True)
    
    
# Calculate water withdrawal per sector for ISIMIP2b simulations and save then in intermediate director
def calc_ptotww_2b_persector(models, forcings, scenarios, da_cellarea): 

    # calculate potential total water withdrawal for isimip2b simulations, besed on defenitions in main_2b noteboook
    for model in models:
        for forcing in forcings: 
            for scenario in scenarios:

                intermediatedir = outdir+'intermediate/'+model.lower()+'/'

                if model == 'H08': 

                    filedir, filename = get_isimip_simulation_name_dir('amanww', datadir, model, forcing, scenario)
                    
                    if not os.path.exists(intermediatedir+filename.replace('amanww', 'indww')): 
                        print('processing indww for '+model + ' ' + forcing + ' '+scenario)  
                        indww = load_simulation('amanww', datadir, model, forcing, scenario)
                        indww.to_dataset(name='indww').to_netcdf(intermediatedir+filename.replace('amanww', 'indww'))
                    
                    if not os.path.exists(intermediatedir+filename.replace('amanww', 'domww')): 
                        print('processing domww for '+model + ' ' + forcing + ' '+scenario)
                        domww = load_simulation('adomww', datadir, model, forcing, scenario)
                        domww.to_dataset(name='domww').to_netcdf(intermediatedir+filename.replace('amanww', 'domww'))
 
                    if not os.path.exists(intermediatedir+filename.replace('amanww', 'irrww')): 
                        print('processing irrww for '+model + ' ' + forcing + ' '+scenario)
                        irrww = load_simulation('pirrww', datadir, model, forcing, scenario)
                        irrww.to_dataset(name='irrww').to_netcdf(intermediatedir+filename.replace('amanww', 'irrww'))

                
                elif model == 'CWatM': 
                    filedir, filename = get_isimip_simulation_name_dir('adomuse', datadir, model, forcing, scenario)

                    if not os.path.exists(intermediatedir+filename.replace('adomuse', 'domww')): 
                        print('processing indww for '+model + ' ' + forcing + ' '+scenario)                      
                        domww = load_simulation('adomuse', datadir, model, forcing, scenario)
                        domww.to_dataset(name='domww').to_netcdf(intermediatedir+filename.replace('adomuse', 'domww'))
                        
                    if not os.path.exists(intermediatedir+filename.replace('adomuse', 'indww')):
                        print('processing indww for '+model + ' ' + forcing + ' '+scenario)                      
                        indww = load_simulation('ainduse', datadir, model, forcing, scenario,)
                        indww.to_dataset(name='indww').to_netcdf(intermediatedir+filename.replace('adomuse', 'indww'))
                        
                    if not os.path.exists(intermediatedir+filename.replace('adomuse', 'irrww')): 
                        print('processing indww for '+model + ' ' + forcing + ' '+scenario)  
                        irrww = load_simulation('pirrww', datadir, model, forcing, scenario)                    
                        irrww.to_dataset(name='irrww').to_netcdf(intermediatedir+filename.replace('adomuse', 'irrww'))


                elif model == 'LPJmL' or model == 'MATSIRO': 
                    
                    filedir, filename = get_isimip_simulation_name_dir('pirrww', datadir, model, forcing, scenario)
                    
                    if not os.path.exists(intermediatedir+filename.replace('pirrww', 'irrww')): 
                        print('processing irrww for '+model + ' ' + forcing + ' '+scenario)  

                        # load potential irrigation water withdrawal (available for both)
                        irrww = load_simulation('pirrww', datadir, model, forcing, scenario)
                        irrww.to_dataset(name='irrww').to_netcdf(intermediatedir+filename.replace('pirrww', 'irrww'))

                        
                    if not (os.path.exists(intermediatedir+filename.replace('pirrww', 'indww')) or os.path.exists(intermediatedir+filename.replace('pirrww', 'domww'))): 
                        print('processing indww and domww for '+model + ' ' + forcing + ' '+scenario)  
                        
                        irrww= xr.open_dataset(intermediatedir+filename.replace('pirrww', 'irrww'))

                        
                        if scenario == 'historical': 
                            inputdir = datadir+'/InputData/water_abstraction/histsoc/'
                            domww_m3yr = xr.open_dataset(inputdir+'domww_histsoc_annual_1901-2005.nc', decode_times=False)['domww']
                            domww_m3yr = domww_m3yr.assign_coords(time = pd.date_range(start='1901', end='2006', freq='Y'))
                            indww_m3yr = xr.open_dataset(inputdir+'indww_histsoc_annual_1901-2005.nc', decode_times=False)['indww']
                            indww_m3yr = indww_m3yr.assign_coords(time = pd.date_range(start='1901', end='2006', freq='Y'))

                            # extend historical period based on pirr period
                            domww_m3yr = match_time_period(domww_m3yr, irrww)
                            indww_m3yr = match_time_period(indww_m3yr, irrww)

                        else: 
                            inputdir= datadir+'/InputData/water_abstraction/'+scenario+'soc/'

                            domww_m3yr = xr.open_dataset(inputdir+'domww_rcp26soc_annual_2006-2099.nc', decode_times=False)['domww']
                            domww_m3yr = domww_m3yr.assign_coords(time = pd.date_range(start='2006', end='2100', freq='Y'))
                            indww_m3yr = xr.open_dataset(inputdir+'indww_rcp26soc_annual_2006-2099.nc', decode_times=False)['indww']
                            indww_m3yr = indww_m3yr.assign_coords(time = pd.date_range(start='2006', end='2100', freq='Y'))


                        # transform the units from m3/year to mm/m2s
                        domww = (domww_m3yr * 1000 /(da_cellarea * sec_per_year))
                        domww = domww.where(domww<1e34,np.nan)
                        indww = (indww_m3yr * 1000 /(da_cellarea * sec_per_year))
                        indww = indww.where(indww<1e34,np.nan)

                        # Convert input water aithdrawal to monthly. 
                        domww = conv_annual2monthly(domww)
                        indww = conv_annual2monthly(indww)

                        # Save the variable into intermediate netCDF file
                        indww.to_dataset(name='indww').to_netcdf(intermediatedir+filename.replace('pirrww', 'indww'))
                        domww.to_dataset(name='domww').to_netcdf(intermediatedir+filename.replace('pirrww', 'domww'))

                        # water withdrawal per sector

def calc_water_withdrawal_persector(var, model, forcing, scenario, da_cellarea, landmask, pop):

    ww = load_simulation(var, datadir, model, forcing, scenario, flag_postprocessed=True)

    days_in_month = []
    for q_avail_month in ww.time.values:
        days_in_month.append(pd.to_datetime(q_avail_month).days_in_month)

    secs_in_month = np.expand_dims(np.array(days_in_month) * 24*60*60, axis=(1,2))

    if scenario == 'historical': 
        time_ds = xr.open_dataset(lifetimedir+'falkenmark/lpjml/lpjml_gfdl-esm2m_historical_falkenmark_global_annual_landarea_1861_2005.nc4')['exposure'].time
    else: 
        time_ds = xr.open_dataset(lifetimedir+'falkenmark/lpjml/lpjml_gfdl-esm2m_rcp26_falkenmark_global_annual_landarea_2006_2099.nc4')['exposure'].time


    withdrawal = (ww*da_cellarea/1000 * secs_in_month).groupby('time.year').sum().rename({'year':'time'})



    # fix the time dimension. 
    if scenario == 'historical': 
        withdrawal  = withdrawal.assign_coords(time = pd.date_range(start='1861', end='2006', freq='Y'))

    else: 
        withdrawal  = withdrawal.assign_coords(time = pd.date_range(start='2006', end='2100', freq='Y'))

    # save into new variables
    save_to_lifetimedir(withdrawal, var, model, forcing, scenario, annual = True)