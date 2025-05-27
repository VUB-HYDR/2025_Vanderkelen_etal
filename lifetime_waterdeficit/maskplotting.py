# ---------------------------------------------------------------
# Functions to make easy plots for masked results
# ----------------------------------------------------------------

import os
import xarray as xr
import pickle as pk
import time
import matplotlib as mpl
from matplotlib.patches import Rectangle

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as colors

from settings import *
from load_manip import *
from exposure import * 
from utils import *
from main_function import *

# Script to plot countries
from plotting import *

# set own plotting parameters
set_plot_param()
map_proj = ccrs.Robinson(central_longitude=0, globe=None)


labelfontsize=16
titlesize=16
labelsize=titlesize
mpl.rc('xtick',labelsize=labelfontsize)
mpl.rc('ytick',labelsize=labelfontsize)
mpl.rc('axes',titlesize=titlesize)
mpl.rc('axes',labelsize=labelsize)
mpl.rc('legend',fontsize='large')


# extreme event
global flags

flags = {}
flags['extr'] = 'waterdeficit'   # 0: all
                                  # 2: cropfailedarea
                                  # 3: waterscarcity
                                  # 4: falkenmark 
                                  # 5: water deficit

flags['runs'] = 0           # 0: do not process ISIMIP runs (i.e. load runs pickle)
                            # 1: process ISIMIP runs (i.e. produce and save runs as pickle)
flags['mask'] = 0           # 0: do not process country data (i.e. load masks pickle)
                            # 1: process country data (i.e. produce and save masks as pickle)
flags['exposure'] = 0     # 0: do not process ISIMIP runs to compute exposure (i.e. load exposure pickle)
                            # 1: process ISIMIP runs to compute exposure (i.e. produce and save exposure as pickle)
flags['exposure_pic'] = 0   # 0: do not process ISIMIP runs to compute picontrol exposure (i.e. load exposure pickle)
                            # 1: process ISIMIP runs to compute picontrol exposure (i.e. produce and save exposure as pickle)


    
# settings for plotting
    
flag_region = 'world'    
regions = {'world': ['East Asia & Pacific', 'Europe & Central Asia', 'Latin America & Caribbean', 'Middle East & North Africa', 'North America', 'South Asia', 'Sub-Saharan Africa'], 
           'income' :[  'Lower middle income','Low income', 'High income','Upper middle income',]}

regions_abbrevs = {'world': ['EASP','EUCA','LAMC','MENA','NAM','SAS','SSA'], 
                  'income' : [  'Lower middle income','Low income', 'High income','Upper middle income',]}

#colors = {'world': ['wheat', 'lightsteelblue','palegoldenrod','thistle','darkseagreen','lightpink','paleturquoise'],
#         'income' : ['lightsteelblue','darkseagreen','wheat','lightpink']}

#colors = {'world': {'Europe & Central Asia' : 'tab:blue', 'Latin America & Caribbean':'tab:orange','East Asia & Pacific': 'tab:green','South Asia': 'tab:pink','Sub-Saharan Africa': 'tab:red','North America':'tab:purple',  'Middle East & North Africa': 'tab:brown'},
#         'income' : {'Upper middle income': 'tab:blue','Lower middle income': 'tab:orange','High income' : 'tab:green','Low income':'tab:red'}}

rcps = ['RCP60','RCP26']
rcp_text = {'RCP60': 'RCP 6.0', 'RCP26' : 'RCP 2.6'  }

panellabels = ['a.','b.','c.','d.','e.','f.','g.','h.','i.','j.','k.','l.','m.','n.']




def plot_lineplot(ds_pct_waterdeficit_perregion, mode):

    
    colors_toplot = {'world': {'Europe & Central Asia' : 'tab:blue', 'Latin America & Caribbean':'tab:orange','East Asia & Pacific': 'tab:green','South Asia': 'tab:pink','Sub-Saharan Africa'
                               :'tab:red','North America':'tab:purple',  'Middle East & North Africa': 'tab:brown'},
             'income' : {'Upper middle income': 'tab:blue','Lower middle income': 'tab:orange','High income' : 'tab:green','Low income':'tab:red'}}

    flag_region = 'world'    
    regions = {'world': ['East Asia & Pacific', 'Europe & Central Asia', 'Latin America & Caribbean', 'Middle East & North Africa', 'North America', 'South Asia', 'Sub-Saharan Africa'], 
               'income' :[  'Lower middle income','Low income', 'High income','Upper middle income',]}


    # plot delta compared to 1960 birth cohort

    fig = plt.figure(figsize=(12,6))

    ax = plt.subplot2grid(shape=(1,1), loc=(0,0), colspan=1)


    flag_uncertainty = 'std' # std

    var_type = ''#'delta_'

    rcp = 'RCP60'
    var = 'mmm_'+var_type+rcp


    # sort regions based on values for 2020 birth year (for plotting uncertainty bar)
    da_tosort = ds_pct_waterdeficit_perregion.sel(birth_year = 2020)[var]
    regions_sorted = da_tosort.sortby(da_tosort, ascending=False)['region'].values
    regions_toloop = [region for region in regions_sorted if region in regions[flag_region]]

    for j, region in enumerate(regions[flag_region]):  
                    

            ds_pct_waterdeficit_perregion.sel({'region':region})[var].plot(ax=ax, label=region, linewidth=2, color = colors_toplot[flag_region][region])

            # draw uncertainty rectangles

            # get values
            end_year = ds_pct_waterdeficit_perregion.birth_year.max()
            lqntl = ds_pct_waterdeficit_perregion.sel({'region':region})['lqntl_'+var_type+rcp].sel(birth_year = end_year).values
            uqntl = ds_pct_waterdeficit_perregion.sel({'region':region})['uqntl_'+var_type+rcp].sel(birth_year = end_year).values

            std = ds_pct_waterdeficit_perregion.sel({'region':region})['std_'+var_type+rcp].sel(birth_year = end_year).values
            mean = ds_pct_waterdeficit_perregion.sel({'region':region})['mmm_'+var_type+rcp].sel(birth_year = end_year).values

            if flag_uncertainty == 'std': 
                lower = mean - std 
                upper = mean + std                
            else: 
                lower = lqntl 
                upper = uqntl

            width = upper - lower

            # define width of the mean bar
            mean_width = 0.5
            offset = 0.5

            # quantile        
            ax.add_patch(Rectangle((end_year+1+j+j*offset,lower), 1, width, facecolor=colors_toplot[flag_region][region], clip_on=False,linewidth = 0, alpha=0.5))
            ax.add_patch(Rectangle((end_year+1+j+j*offset,mean-mean_width/2), 1, mean_width, facecolor=colors_toplot[flag_region][region], clip_on=False,linewidth = 0))     



    ax.legend(bbox_to_anchor=(1.52,0.5), loc='right',frameon=False)
    ax.legend(loc='upper left',frameon=False, ncol=2, fontsize=14)

    ax.set_title(mode, loc='right')
    ax.set_title('', loc='center')
    #ax.set_title(panellabels[i], loc='left')

    ax.set_xlim((1960,2020))
    ax.set_ylim((0,100))

    ax.set_xlabel('Birth year')
    ax.set_ylabel('Lifetime water deficit (%)')
    ax.grid(color='lightgray', alpha=0.5)
    ax.spines[['right', 'top']].set_visible(False)        

    fig.tight_layout()


def plot_lineplot_globalmean(ds_pct_waterdeficit_perregion, mode):

    
    colors_toplot = {'world': {'Europe & Central Asia' : 'tab:blue', 'Latin America & Caribbean':'tab:orange','East Asia & Pacific': 'tab:green','South Asia': 'tab:pink','Sub-Saharan Africa'
                               :'tab:red','North America':'tab:purple',  'Middle East & North Africa': 'tab:brown'},
             'income' : {'Upper middle income': 'tab:blue','Lower middle income': 'tab:orange','High income' : 'tab:green','Low income':'tab:red'}}

    flag_region = 'world'    


    # plot delta compared to 1960 birth cohort

    fig = plt.figure(figsize=(8,6))

    ax = plt.subplot2grid(shape=(1,1), loc=(0,0), colspan=1)


    flag_uncertainty = 'std' # std

    var_type = ''#'delta_'

    rcp = 'RCP60'
    var = 'mmm_'+var_type+rcp


    # sort regions based on values for 2020 birth year (for plotting uncertainty bar)
    da_tosort = ds_pct_waterdeficit_perregion.sel(birth_year = 2020)[var]
    regions_sorted = da_tosort.sortby(da_tosort, ascending=False)['region'].values
    regions_toloop = [region for region in regions_sorted if region in regions[flag_region]]

    for j, region in enumerate(['World']):  
                    

            ds_pct_waterdeficit_perregion.sel({'region':region})[var].plot(ax=ax, label=region, linewidth=2)

            # draw uncertainty rectangles

            # get values
            end_year = ds_pct_waterdeficit_perregion.birth_year.max()
            lqntl = ds_pct_waterdeficit_perregion.sel({'region':region})['lqntl_'+var_type+rcp].sel(birth_year = end_year).values
            uqntl = ds_pct_waterdeficit_perregion.sel({'region':region})['uqntl_'+var_type+rcp].sel(birth_year = end_year).values

            std = ds_pct_waterdeficit_perregion.sel({'region':region})['std_'+var_type+rcp].sel(birth_year = end_year).values
            mean = ds_pct_waterdeficit_perregion.sel({'region':region})['mmm_'+var_type+rcp].sel(birth_year = end_year).values

            if flag_uncertainty == 'std': 
                lower = mean - std 
                upper = mean + std                
            else: 
                lower = lqntl 
                upper = uqntl

            width = upper - lower

            # define width of the mean bar
            mean_width = 0.5
            offset = 0.5

            # quantile        
            ax.add_patch(Rectangle((end_year+1+j+j*offset,lower), 1, width,  clip_on=False,linewidth = 0, alpha=0.5))
            ax.add_patch(Rectangle((end_year+1+j+j*offset,mean-mean_width/2), 1, mean_width,  clip_on=False,linewidth = 0))     



    ax.legend(bbox_to_anchor=(1.52,0.5), loc='right',frameon=False)
    ax.legend(loc='upper left',frameon=False, ncol=2, fontsize=14)

    ax.set_title(mode, loc='right')
    ax.set_title('', loc='center')
    #ax.set_title(panellabels[i], loc='left')

    ax.set_xlim((1960,2020))
    ax.set_ylim((0,100))

    ax.set_xlabel('Birth year')
    ax.set_ylabel('Lifetime water deficit (%)')
    ax.grid(color='lightgray', alpha=0.5)
    ax.spines[['right', 'top']].set_visible(False)        

    fig.tight_layout()    
def plot_lollipop(ds_pct_waterdeficit_perregion, mode):

    fig, ax = plt.subplots(1,1, figsize=(6,4))
    birth_year = 2020


    df_rcp26 = ds_pct_waterdeficit_perregion[ 'mmm_RCP26'].sel(birth_year=birth_year).to_dataframe().drop(['birth_year'],axis=1)
    values_rcp60 = ds_pct_waterdeficit_perregion[ 'mmm_RCP60'].sel(birth_year=birth_year).values

    values = np.stack((df_rcp26.index.values,np.squeeze(df_rcp26.values),values_rcp60))
    df = pd.DataFrame(values.T, columns=['regions','RCP26','RCP60'])

    df = df[df["regions"].isin(regions[flag_region])]

    regions_ordered = df.sort_values(by='RCP26',ascending=True)["regions"].values

    df['regions_cat'] = pd.Categorical(df['regions'], categories=regions_ordered, ordered=True)
    df = df.sort_values(by='regions_cat')

    my_range=range(1,len(df.index)+1)


    ax.hlines(y=my_range, xmin=df['RCP26'], xmax=df['RCP60'], color='grey', alpha=0.4)
    ax.scatter(df['RCP26'], my_range, color='green', alpha=0.7, s=50, label='RCP 2.6')
    ax.scatter(df['RCP60'], my_range, color='tab:orange', alpha=0.7 , s=50,label='RCP 6.0')
    ax.legend(loc='lower right');

    ax.set_yticks(my_range)
    ax.set_yticklabels(df['regions']); 

    ax.spines[['right', 'left', 'top']].set_visible(False)    

    ax.set_xlabel('Lifetime water deficit (%)')
    ax.set_title(mode, loc='right', fontsize=18);
    ax.grid(color='lightgray', alpha=0.5)
    ax.set_xlim((15,85));

    fig.tight_layout()



def plot_generation_map(ds_pct_waterdeficit, mode): 
    # load country borders (from pickles)
    d_countries = pk.load(open('./data/pickles/country_info.pkl', 'rb'))

    # unpack country information
    gdf_country_borders = d_countries['borders']


    variable = 'mmm_RCP60'
    da = ds_pct_waterdeficit[variable]
    df_2 = da.sel({'birth_year' :2020}).to_dataframe()
    df = da.sel({'birth_year' : 1960}).to_dataframe()

    d_waterdeficit = (df_2-df)

    legend_label = '$\Delta$ relative water deficit (%)'
    gdf_emf = gdf_country_borders.merge(d_waterdeficit, left_index=True, right_on='country')
    bounds = np.array([-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30])

    fig, ax = plt.subplots(figsize=(18,10), subplot_kw={'projection':map_proj})
    ax.axis('off')
    ax.coastlines(color='lightgray',linewidth=0.5)
    ax.add_feature(cfeature.BORDERS,color='lightgray',linewidth=0.5)

    cax = ax.inset_axes((1.02, 0.1, 0.03, 0.9)); #make a color bar axis

    gdf_emf.plot(ax=ax, column=variable, legend = True,  cmap='BrBG_r',  cax=cax, legend_kwds={'label': legend_label}, vmin=-30, vmax=30, transform=ccrs.PlateCarree(), norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256, extend='both'))
    ax.set_title('$\Delta$ Relative lifetime water deficit (2020-1960) '+mode,  loc='right');
    
    
def plot_mitigation_map(ds_pct_waterdeficit, mode):
    
    
    # load country borders (from pickles)
    d_countries = pk.load(open('./data/pickles/country_info.pkl', 'rb'))

    # unpack country information
    gdf_country_borders = d_countries['borders']
    
    variable = 'delta_60_26'

    d_waterdeficit_60= ds_pct_waterdeficit['mmm_RCP60'].sel({'birth_year' : 2020}).to_dataframe()
    d_waterdeficit_26= ds_pct_waterdeficit['mmm_RCP26'].sel({'birth_year' : 2020}).to_dataframe()

    d_waterdeficit_60_26 = (d_waterdeficit_60['mmm_RCP60'] - d_waterdeficit_26['mmm_RCP26']).rename(variable)


    legend_label = '$\Delta$ relative water deficit (%)'
    gdf_emf = gdf_country_borders.merge(d_waterdeficit_60_26, left_index=True, right_on='country')
    bounds = np.array([-15,-12.5,-10,-7.5,-5,-2.5,0,2.5,5,7.5,10,12.5,15])
    #bounds = np.array([-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30])


    fig, ax = plt.subplots(figsize=(18,10), subplot_kw={'projection':map_proj})
    ax.axis('off')
    ax.coastlines(color='lightgray',linewidth=0.5)
    ax.add_feature(cfeature.BORDERS,color='lightgray',linewidth=0.5)

    cax = ax.inset_axes((1.02, 0.1, 0.03, 0.9)); #make a color bar axis

    gdf_emf.plot(ax=ax, column=variable, legend = True,  cmap='RdBu_r',  cax=cax, legend_kwds={'label': legend_label}, vmin=-30, vmax=30, transform=ccrs.PlateCarree(), norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256, extend='both'))
    ax.set_title('$\Delta$ Relative lifetime water deficit for person born in 2020, RCP 6.0 - RCP 2.6 '+mode,  loc='right');
    
    
    
    
def plot_barplots(dryland_mode, ds_pct_waterdeficit):

    BE_colors = ['gold','darkorange','crimson', 'maroon', 'purple']
    flag_region = 'world'    
    regions = ['East Asia & Pacific', 'Europe & Central Asia', 'Latin America & Caribbean', 'Middle East & North Africa', 'North America', 'South Asia', 'Sub-Saharan Africa']

    # load country borders (from pickles)

    d_countries = pk.load(open('./data/pickles/country_info.pkl', 'rb'))
    df_countries = d_countries['info_pop']

    # unpack country information
    gdf_country_borders = d_countries['borders']

    # load cohort info
    d_regions = pk.load(open('./data/pickles/region_info.pkl', 'rb'))




    # calculated weighted fieldmean per country mask
    def calc_dryland_pop(
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

        da_weighted_fldmean = da_masked.weighted(weights).sum(dim=("lat", "lon"))

        return da_weighted_fldmean



    # calculated weighted fieldmean per country mask
    def calc_total_pop(
        da, 
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

        da_weighted_fldmean = da_masked.sum(dim=("lat", "lon"))

        return da_weighted_fldmean


    # calculate population fraction in drylands

    d_countries = pk.load(open('./data/pickles/country_info.pkl', 'rb'))
    df_countries = d_countries['info_pop']
    # load cohort info
    d_regions = pk.load(open('./data/pickles/region_info.pkl', 'rb'))

    # Load necessary data
    countries_regions, countries_mask = d_countries['mask']

    # load griddded population
    da_pop =  load_population(year_start,year_end)

    # load multi-model mean dryland mask
    if dryland_mode == 'urban':
        dryland_mask = xr.open_dataset('./data/masks/urbanshare_2020.nc', decode_times=False)['mask'].fillna(0)
        
    
    elif dryland_mode == 'exposure' :
        dryland_mask = (da_pop > 0).fillna(0)

    else:     
        dryland_mask = xr.open_dataset('./data/masks/'+dryland_mode+'.nc')['mask'].fillna(0)

    # per country
    frac_popindryland = {}

    for j, country in enumerate(df_countries['name']):

        print('processing country '+str(j+1)+' of '+str(len(df_countries)), end='\r')

        ind_country = countries_regions.map_keys(country)

        flag_region = False

        totpop_country = calc_total_pop(da_pop,countries_mask, ind_country, flag_region)
        drylandpop_country = calc_dryland_pop(da_pop, dryland_mask, countries_mask, ind_country, flag_region)

        frac_popindryland[country] = drylandpop_country/totpop_country



    # loop over regions
    for k, region in enumerate(d_regions['birth_years'].columns): 

        print('processing region '+str(k+1)+' of '+str(len(d_regions['birth_years'].columns)), end='\r')

        # Get list of member countries from region - with seperate treatment for world (luke: now inside get_countries_of_regions func)
        member_countries = get_countries_of_region(region, df_countries)

        # get spatial average of landfraction: historical + RCP simulations
        ind_countries = countries_regions.map_keys(member_countries)


        flag_region = True

        totpop_region = calc_total_pop(da_pop, countries_mask, ind_countries, flag_region)
        drylandpop_region = calc_dryland_pop(da_pop, dryland_mask, countries_mask, ind_countries, flag_region)

        frac_popindryland[region] = drylandpop_region/totpop_region


    # get cohort size per country
    df_cohortsize_percountry_2020 = pd.DataFrame() 

    for country in d_countries['cohort_size']:

        if country in frac_popindryland.keys():
            coef = frac_popindryland[country].sel(time=2020).values
        else: 
            coef = 0

        df_cohortsize_percountry_2020[country] =  d_countries['cohort_size'][country].loc[2020] * 1000 * coef


    # Define cohort classes and labels
    cohort_classes       = [[0,9], [10,19], [20,29], [30,39], [40,49], [50,59]]
    cohort_classes_label = ['0-9', '10-19', '20-29', '30-39', '40-49', '50-59']

    # get the cohortsizes per birth year, per region
    df_pop_perregion = get_cohortsizes_perregion()

    # group the cohorts into classes (calculate sum of # people within cohorts)
    df_cohort_grouped = group_cohorts_into_classes(df_pop_perregion, cohort_classes, cohort_classes_label).iloc[::-1]

    # calculate share of every cohort class for every region
    cohort_grouped_share = df_cohort_grouped.div(df_cohort_grouped.sum(axis=0), axis=1)

    df_countries 
    # calculate country weights based on total region population

    d_country_weights= {region: df_countries[df_countries['region']==region]['population'].div(float(df_countries.groupby('region').sum().loc[region].values)) for region in regions}
    d_country_weights['World'] = df_countries['population'] / df_countries['population'].sum()


    #### RCP 60

    # convert data array into dataframe for plotting
    var = 'mmm_RCP60'
    df_pct_waterdeficit = pd.DataFrame(ds_pct_waterdeficit[var].values, columns = ds_pct_waterdeficit[var]['country'], index = df_cohortsize_percountry_2020.index)

    var = 'lqntl_RCP60'
    df_pct_waterdeficit_lqntl = pd.DataFrame(ds_pct_waterdeficit[var].values, columns = ds_pct_waterdeficit[var]['country'], index = df_cohortsize_percountry_2020.index)

    var = 'uqntl_RCP60'
    df_pct_waterdeficit_uqntl = pd.DataFrame(ds_pct_waterdeficit[var].values, columns = ds_pct_waterdeficit[var]['country'], index = df_cohortsize_percountry_2020.index)



    # per birth year and country calculate number of people above certain lifetime water deficit threshold and then calculate weighted average per region and globally

    thresholds = [50,75,90] # %
    thresholds = [50,60,70,80,90] # %

    d_npeople_ot = {}

    for threshold in thresholds:

        df_cohortsize_percountry_overthreshold = df_cohortsize_percountry_2020.where(df_pct_waterdeficit > threshold)

        df_npeople_ot_perbirthyear = pd.DataFrame()

        # group the countries by region
        for region in regions: 

            countries_in_region = list(df_countries[df_countries['region'] == region].index)

            # calculate the number of people with lifetime water deficit over the threshold per region (based on country average values)
            # THIS CALCULATION SHOULD BE THE WEIGHTED SUM
            # df_npeople_ot_perbirthyear[region] = (df_cohortsize_percountry_overthreshold[countries_in_region] * d_country_weights[region]).sum(axis =1)

            # this should not be weighted but just be summed? 
            df_npeople_ot_perbirthyear[region] = (df_cohortsize_percountry_overthreshold[countries_in_region]).sum(axis =1)


        # for the world, weighted based on all countries (not average per region!!!!)
        df_npeople_ot_perbirthyear['World'] =  (df_cohortsize_percountry_overthreshold * d_country_weights['World']).sum(axis =1)
        df_npeople_ot_perbirthyear['World'] =  (df_cohortsize_percountry_overthreshold).sum(axis =1)


        # group number of people per cohort 
        d_npeople_ot[threshold] = group_cohorts_into_classes(df_npeople_ot_perbirthyear, cohort_classes, cohort_classes_label).iloc[::-1]

    # turn dictionary into dataset  
    ds_npeople = xr.concat(
        [xr.DataArray(v).rename({'dim_0':'cohort_class','dim_1':'region'}) for v in d_npeople_ot.values()],
        dim='threshold').assign_coords({'threshold':list(d_npeople_ot.keys())})


    # million people
    thresholds_to_show = thresholds #[:-1]

    region = 'World'
    # define burning ember colors

    fig, axes = plt.subplots(1,2, figsize=(16,5))
    axes = axes.flatten()
    ax = axes[0]
    ax.set_ylabel('million people')
    ax.set_xlabel('Age cohorts');
    #ax.set_title('Number of people per age cohort worldwide',loc='right')
    ax.spines[['right', 'top']].set_visible(False)    
    ax.set_title(dryland_mode+' RCP 6.0', loc='right', fontsize=18);


    # get people per threshold for region

    df_toplot = df_npeople_thresholds_perregion = pd.DataFrame(ds_npeople.sel(region=region).values.T, index=  ds_npeople.cohort_class.values ,columns =ds_npeople.threshold.values) *1e-6

    df_toplot.iloc[:,0].plot.bar(ax=ax, color=BE_colors[0], label=str(df_toplot.keys()[0])+' % water deficit', rot=0)

    for n in range(1,len(thresholds_to_show)):
        df_toplot.iloc[:,n].plot.bar(ax=ax, color=BE_colors[n],  label=str(df_toplot.keys()[n])+' % water deficit', alpha=0.5, rot=0)

    leg = ax.legend(bbox_to_anchor=(1,1), loc="upper left",frameon=False)    
    for lh in leg.legendHandles: lh.set_alpha(1)


    # get people multiplication factor
    pmf = df_toplot/df_toplot.iloc[0,:]

    n=0
    for j,p in enumerate(ax.patches):
        if j < 6: 
            ax.annotate('x '+str(np.round(pmf.iloc[:,0].values,1)[j]), (p.get_x() + 0.1, p.get_height() + 15))
        elif j > 6 and j < 12: 

            ax.annotate('x '+str(np.round(pmf.iloc[:,1].values,1)[j-6]), (p.get_x() + 0.1, p.get_height() + 15))
        elif j > 12 and j < 18: 

            ax.annotate('x '+str(np.round(pmf.iloc[:,1].values,1)[j-12]), (p.get_x() + 0.1, p.get_height() + 10))

            
    #### RCP 26

    # convert data array into dataframe for plotting
    var = 'mmm_RCP26'
    df_pct_waterdeficit = pd.DataFrame(ds_pct_waterdeficit[var].values, columns = ds_pct_waterdeficit[var]['country'], index = df_cohortsize_percountry_2020.index)

    var = 'lqntl_RCP26'
    df_pct_waterdeficit_lqntl = pd.DataFrame(ds_pct_waterdeficit[var].values, columns = ds_pct_waterdeficit[var]['country'], index = df_cohortsize_percountry_2020.index)

    var = 'uqntl_RCP26'
    df_pct_waterdeficit_uqntl = pd.DataFrame(ds_pct_waterdeficit[var].values, columns = ds_pct_waterdeficit[var]['country'], index = df_cohortsize_percountry_2020.index)



    # per birth year and country calculate number of people above certain lifetime water deficit threshold and then calculate weighted average per region and globally

    thresholds = [50,75,90] # %
    thresholds = [50,60,70,80,90] # %

    d_npeople_ot = {}

    for threshold in thresholds:

        df_cohortsize_percountry_overthreshold = df_cohortsize_percountry_2020.where(df_pct_waterdeficit > threshold)

        df_npeople_ot_perbirthyear = pd.DataFrame()

        # group the countries by region
        for region in regions: 

            countries_in_region = list(df_countries[df_countries['region'] == region].index)

            # calculate the number of people with lifetime water deficit over the threshold per region (based on country average values)
            # THIS CALCULATION SHOULD BE THE WEIGHTED SUM
            # df_npeople_ot_perbirthyear[region] = (df_cohortsize_percountry_overthreshold[countries_in_region] * d_country_weights[region]).sum(axis =1)

            # this should not be weighted but just be summed? 
            df_npeople_ot_perbirthyear[region] = (df_cohortsize_percountry_overthreshold[countries_in_region]).sum(axis =1)


        # for the world, weighted based on all countries (not average per region!!!!)
        df_npeople_ot_perbirthyear['World'] =  (df_cohortsize_percountry_overthreshold * d_country_weights['World']).sum(axis =1)
        df_npeople_ot_perbirthyear['World'] =  (df_cohortsize_percountry_overthreshold).sum(axis =1)


        # group number of people per cohort 
        d_npeople_ot[threshold] = group_cohorts_into_classes(df_npeople_ot_perbirthyear, cohort_classes, cohort_classes_label).iloc[::-1]

    # turn dictionary into dataset  
    ds_npeople = xr.concat(
        [xr.DataArray(v).rename({'dim_0':'cohort_class','dim_1':'region'}) for v in d_npeople_ot.values()],
        dim='threshold').assign_coords({'threshold':list(d_npeople_ot.keys())})


    # million people
    thresholds_to_show = thresholds #[:-1]

    region = 'World'
    # define burning ember colors

    ax = axes[1]
    ax.set_ylabel('million people')
    ax.set_xlabel('Age cohorts');
    #ax.set_title('Number of people per age cohort worldwide',loc='right')
    ax.spines[['right', 'top']].set_visible(False)    
    ax.set_title(dryland_mode+' RCP 2.6', loc='right', fontsize=18);


    # get people per threshold for region

    df_toplot = df_npeople_thresholds_perregion = pd.DataFrame(ds_npeople.sel(region=region).values.T, index=  ds_npeople.cohort_class.values ,columns =ds_npeople.threshold.values) *1e-6

    df_toplot.iloc[:,0].plot.bar(ax=ax, color=BE_colors[0], label=str(df_toplot.keys()[0])+' % water deficit', rot=0)

    for n in range(1,len(thresholds_to_show)):
        df_toplot.iloc[:,n].plot.bar(ax=ax, color=BE_colors[n],  label=str(df_toplot.keys()[n])+' % water deficit', alpha=0.5, rot=0)

    leg = ax.legend(bbox_to_anchor=(1,1), loc="upper left",frameon=False)    
    for lh in leg.legendHandles: lh.set_alpha(1)


    # get people multiplication factor
    pmf = df_toplot/df_toplot.iloc[0,:]

    n=0
    for j,p in enumerate(ax.patches):
        if j < 6: 
            ax.annotate('x '+str(np.round(pmf.iloc[:,0].values,1)[j]), (p.get_x() + 0.1, p.get_height() + 15))
        elif j > 6 and j < 12: 

            ax.annotate('x '+str(np.round(pmf.iloc[:,1].values,1)[j-6]), (p.get_x() + 0.1, p.get_height() + 15))
        elif j > 12 and j < 18: 

            ax.annotate('x '+str(np.round(pmf.iloc[:,1].values,1)[j-12]), (p.get_x() + 0.1, p.get_height() + 10))
    fig.tight_layout()
    
    
def get_globalnumbers(ds_pct_waterdeficit_perregion, mode):


    df_2020 = ds_pct_waterdeficit_perregion[ 'mmm_RCP60'].sel(birth_year=2020).to_dataframe().drop(['birth_year'],axis=1)
    values_1960 = ds_pct_waterdeficit_perregion[ 'mmm_RCP60'].sel(birth_year=1960).values

    values = np.stack((df_2020.index.values,np.squeeze(df_2020.values),values_1960))
    df = pd.DataFrame(values.T, columns=['regions','2020','1960'])

    df = df[df["regions"]=='World']

    print('Global mean lifetime water deficit (%)')
    print(' ')
    print('birthyear 1960: '+str(np.round(df['1960'].values[0],2))+ ' %')
    print('birthyear 2020: '+str(np.round(df['2020'].values[0],2))+ ' %')
    print(' ')


    birth_year = 2020
    df_rcp26 = ds_pct_waterdeficit_perregion[ 'mmm_RCP26'].sel(birth_year=birth_year).to_dataframe().drop(['birth_year'],axis=1)
    values_rcp60 = ds_pct_waterdeficit_perregion[ 'mmm_RCP60'].sel(birth_year=birth_year).values

    values = np.stack((df_rcp26.index.values,np.squeeze(df_rcp26.values),values_rcp60))
    df = pd.DataFrame(values.T, columns=['regions','RCP26','RCP60'])

    df = df[df["regions"]=='World']
    print('RCP 2.6: '+str(np.round(df['RCP26'].values[0],2))+ ' %')
    print('RCP 6.0: '+str(np.round(df['RCP60'].values[0],2))+ ' %')
