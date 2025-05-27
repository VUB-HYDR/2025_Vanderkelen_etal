
# ----------------------------------------------------------------
# Functions to do the plotting and cohort postprocessing for plotting
# ----------------------------------------------------------------


import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import matplotlib as mpl
import pickle as pk
import xarray as xr
import numpy as np


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


    
    
# function to load dataframe with population totals per birth year and world region
def get_cohortsizes_perregion(): 
    
    # load regions pickle
    d_regions = pd.read_pickle(open('./data/pickles/region_info.pkl', 'rb'))

    # unpack region information
    d_cohort_weights_regions = d_regions['cohort_size']

    # preprocess cohorts
    d_cohort_per_region = {}
    for key in d_cohort_weights_regions:
        d_cohort_per_region[key] = d_cohort_weights_regions[key].sum(axis=1) *1e3
    
    cohort_per_region = pd.concat(d_cohort_per_region,axis=1)
        
    return cohort_per_region


# function to group the cohort classes into predefined classes (takes the sum over the classes)
def group_cohorts_into_classes(df, cohort_classes, cohort_classes_label):

    df_cohort_grouped = pd.DataFrame()

    for i,cohort_class in enumerate(cohort_classes): 
        
        if len(cohort_class)>1: 
        
            df_grouped = df.loc[df.index.isin(list(np.arange(cohort_class[0],cohort_class[1]+1)))].sum(axis=0).to_frame(cohort_classes_label[i])
        else: 
            df_grouped = df.loc[df.index.isin(cohort_class)].sum(axis=0).to_frame(cohort_classes_label[i])


        if i == 0 : 
                df_cohort_grouped = df_grouped
        else:     
                df_cohort_grouped = df_cohort_grouped.join(df_grouped)

    return df_cohort_grouped.T



















    













