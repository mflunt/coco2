#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 15:18:49 2022

Read in different data files

For each site create a dataframe with CO total and CO2 total in 

Also CO2 ff 

Also worth reading in obs for CO2 and CO. 

@author: mlunt
"""
import xarray
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
import json
import glob
#import inv_mod

def open_ds(fname):
    with xarray.open_dataset(fname) as ds:
        ds.load()
    return ds

def read_csv_data(fname_in, dict_si, species, start, end):
    """
    Read ICOS csv data files
    s
    Return a time indexed dataframe with obs and std dev.
    """
    
    df_temp = pd.read_csv(fname_in, skiprows=dict_si["header_lines"]-1,  sep=';')
    df = (df_temp.where((df_temp["Flag"] == "O") | (df_temp["Flag"] == "U"))).dropna(axis=0, how="all")
    
    a_time=pd.to_datetime(df[['Year', 'Month', 'Day', 'Hour']])
    d_data = {"obs":df[species.lower()].values,
                "std":df["Stdev"].values}
    
    df_site = pd.DataFrame(data = d_data,
              index=a_time)
    df_filt = df_site[slice(start, end)]
    
    ds_filt = xarray.Dataset.from_dataframe(df_filt)
    ds2 = ds_filt.rename({"index": "time"})
    
    return ds2  #, df_site

def read_nc_data(fname_in, species, start, end):
    """
    Read DECC network ntecdf data
    Return a time indexed dataframe with obs and std dev.
    """
    
    ds_temp = open_ds(fname_in)
    ds = ds_temp.sel(time=slice(start, end))
    
    #df_filt = ds[species.lower()]#.to_dataframe(name="obs")
    #df_filt["std"] = ds[species.lower()+"_variability"].values
    
    return ds

def get_site_output(species, sites, run_dir, run_date, start_date, end_date, keys, 
                    read_obs=False, perturb_chem=False, edgar_ratio=None):
    """
    Combine the enesmble data for individual sites into one array
    
    """
    
    y_mod_dict={}
    y_obs_dict={}
    
    for site in sites:
      
        run_str_site = run_dir + "GC_output_" + species + "_" + site + "_" + run_date + ".csv" 
        
        fname_ens = glob.glob(run_str_site)

        if len(fname_ens)<1:
            print ("Exiting loop, no data for " + site + " " + species)
            continue     #skip loopif site file doesn't exist
    
        df_site1_temp = pd.read_csv(fname_ens[0], index_col="time", parse_dates=True)
        df_site1 = df_site1_temp[start_date:end_date]
     
        
        if perturb_chem == True and species == "CO":
            df_site1["COCH4"] = df_site1["COCH4"]*2.
            df_site1["CONMVOC"] = df_site1["CONMVOC"]*2.
        
        site_mod = df_site1[keys].sum(axis=1)
        #site_obs = df_site1["obs"]
        #site_std = df_site1["std"]
        #site_slope = df_site1["slope"]
        
        # For EDGAR need to multiply each column of df_site1 by ratio
        if edgar_ratio is not None:
            bf_keys = []
            for xi in range(109):
                bf_keys.append(species + "_E" + str(xi+1))
            #print(bf_keys)
            
            if species == "CO":
                co_keys = keys[:-109]
                #print(co_keys)
                site_mod_edgar = df_site1[co_keys].sum(axis=1)
            else:
                site_mod_edgar = df_site1[bf_keys[0]]*0.
            for xi, key in enumerate(bf_keys):
                site_mod_edgar = site_mod_edgar + df_site1[key]*edgar_ratio[xi]
                
            y_mod_dict[site + "_edgar"] = site_mod_edgar
            
        y_mod_dict[site] = site_mod
        y_obs_dict[site] = df_site1[["obs", "std", "slope"]]
        #y_obs_dict[site + "_std"] = site_std
        #y_obs_dict[site+ "_slope"] = site_slope
        
    
    if read_obs ==True:
        return y_mod_dict, y_obs_dict
    else:
        return y_mod_dict


#%%
start_date = "20170101"
#end_date = "20170131"
end_date = "20191231"

#run_dates = ["20171201", "20180101"]

pd_run_dates = pd.date_range(start_date, end_date, freq="MS")
run_dates = pd_run_dates.strftime('%Y%m%d')  

pd_end_dates = pd.date_range(start_date, end_date, freq="M")
end_dates = pd_end_dates.strftime('%Y%m%d')  

#run_dates = ["20180101", "20180201"]

species_list=["CO2", "CO"]


#sites = ["GAT", "HTM", "IPR", 
#        "KRE",  "NOR", "OPE",
#        "SAC" ,  "SVB", "TRN", 
#        "TAC",  "BSD", "HFD"]

#sites=["TAC"] # , "HFD"]

sites = ["GAT", "HTM",  
        "KRE",  "NOR", "OPE",
        "SAC" ,  "SVB", "TRN", 
        "TAC",  "BSD", "HFD"]

run_groups = ["gee", "reco" ,"ff"]

nbasis=109
spc_names=[]
spc_names_coff=[]
spc_names_co0=["COIC", "CO_BC1", "CO_BC2", "CO_BC3", "CO_BC4", "COCH4", "CONMVOC"  ]
spc_names_co=["CO_BC1", "CO_BC2", "CO_BC3", "CO_BC4", "COCH4", "CONMVOC"  ]
for xi in range(nbasis):
    spc_names.append("CO2" + "_E" + str(xi+1))
    spc_names_coff.append("CO" + "_E" + str(xi+1))
    spc_names_co0.append("CO" + "_E" + str(xi+1))
    spc_names_co.append("CO" + "_E" + str(xi+1))

data_root =  "/home/mlunt/datastore/DARE/gc_output/inv_runs/"
out_dir = "/home/mlunt/datastore/DARE/ml_data/model_data/"

fname_edgar =  "/home/mlunt/datastore/EDGAR/edgar_tno_bf_ratios.csv"
df_ratio = pd.read_csv(fname_edgar)
edgar_tno_ratio_co2 = df_ratio["co2_ratio"].values
edgar_tno_ratio_co = df_ratio["co_ratio"].values
#%%


bc_site_dir = data_root + "bc_CO2_201701" + "/"  
# Read in CO2 BC separately
y_bc_co2, y_obs_co2 = get_site_output("CO2", sites, bc_site_dir, "20170101", start_date, end_date,
                             ["CO2_BC1"], read_obs=True)

y_mod_sum_list=[]
y_mod_sum_list2=[]
y_mod_ff_list=[]
y_mod_coff_list=[]
y_mod_coff_list2=[]
y_obs_co_list=[]
y_mod_coff_only_list=[]
for run_date in run_dates:
    
    run_name_ff =  "ff_" + "CO2" + "_" + run_date   
    co2ff_out_dir = data_root + run_name_ff + "/"  
    
    # Read in fixed Basis_function names
    y_mod_ff_temp = get_site_output("CO2", sites, co2ff_out_dir, run_date, start_date, end_date,
                                                   spc_names, read_obs=False, edgar_ratio=edgar_tno_ratio_co2)
    
    run_name_gee =  "gee_" + "CO2" + "_" + run_date   
    co2gee_out_dir = data_root + run_name_gee + "/" 
    
    y_mod_gee_temp = get_site_output("CO2", sites, co2gee_out_dir, run_date, start_date, end_date,
                                                   spc_names, read_obs=False)
    
    run_name_reco =  "reco_" + "CO2" + "_" + run_date   
    co2reco_out_dir = data_root + run_name_reco + "/" 
    
    y_mod_reco_temp = get_site_output("CO2", sites, co2reco_out_dir, run_date, start_date, end_date,
                                                   spc_names, read_obs=False)
    
    
    run_name_coff =  "ff_" + "CO" + "_" + run_date   
    coff_out_dir = data_root + run_name_coff + "/"  
    
    # Read in CO BC and ff all at once (and IC)
    
    
    if run_date == run_dates[0]:
        y_mod_coff_temp, y_obs_co_temp = get_site_output("CO", sites, coff_out_dir, run_date, start_date, end_date,
                                                       spc_names_co0, read_obs=True, 
                                                       edgar_ratio = edgar_tno_ratio_co)
        
        y_mod_co_pert_temp = get_site_output("CO", sites, coff_out_dir, run_date, start_date, end_date,
                                                       spc_names_co0, read_obs=False, perturb_chem=True,
                                                       edgar_ratio = edgar_tno_ratio_co)
        
        y_mod_coff_only = get_site_output("CO", sites, coff_out_dir, run_date, start_date, end_date,
                                                       spc_names_coff, read_obs=False, 
                                                       edgar_ratio = edgar_tno_ratio_co)
        
    else:
        y_mod_coff_temp, y_obs_co_temp = get_site_output("CO", sites, coff_out_dir, run_date, start_date, end_date,
                                                       spc_names_co, read_obs=True, perturb_chem=False,
                                                       edgar_ratio = edgar_tno_ratio_co)
    
    
        y_mod_co_pert_temp = get_site_output("CO", sites, coff_out_dir, run_date, start_date, end_date,
                                                       spc_names_co, read_obs=False, perturb_chem=True,
                                                       edgar_ratio = edgar_tno_ratio_co)
        
        y_mod_coff_only= get_site_output("CO", sites, coff_out_dir, run_date, start_date, end_date,
                                                       spc_names_coff, read_obs=False, perturb_chem=False,
                                                       edgar_ratio = edgar_tno_ratio_co)
    

    y_mod_sum_rdate={}
    y_mod_sum_rdate2={}
    for site in sites:
        if site in y_mod_ff_temp.keys():
            y_mod_sum_rdate[site] = (y_mod_ff_temp[site] + 
                           y_mod_gee_temp[site] + 
                           y_mod_reco_temp[site] )
            
            # Create another variable with different bio emissions to see how ml model performs
            y_mod_sum_rdate2[site] = (y_mod_ff_temp[site] + 
                           y_mod_gee_temp[site]*1.05 + 
                           y_mod_reco_temp[site]*1.1 )
        
                       
     # Add model sum and ff to a list   
     
    y_mod_sum_list.append(y_mod_sum_rdate)
    y_mod_sum_list2.append(y_mod_sum_rdate2)
    y_mod_ff_list.append(y_mod_ff_temp)
    
    y_mod_coff_list.append(y_mod_coff_temp)
    y_mod_coff_list2.append(y_mod_co_pert_temp)
    y_mod_coff_only_list.append(y_mod_coff_only)
    y_obs_co_list.append(y_obs_co_temp)
    
    # Also need t oadd CO - bc, IC, chem and bfs 


#%%

y_mod_sum_site={}
y_mod_sum_site2={}
y_mod_ff_site={}
y_mod_edgar_site={}
y_mod_co_site={}
y_mod_co_site2={}
y_mod_co_edgar_site={}
y_obs_co_site={}
y_std_co_site={}
y_slope_co_site={}
y_mod_coff_site={}
for site in sites:
    
    dum_list_sum=[]
    dum_list_sum2=[]
    dum_list_ff=[]
    dum_list_edgar=[]
    dum_list_coff=[]
    dum_list_coff2=[]
    dum_list_co_obs=[]
    dum_list_co_std=[]
    dum_list_co_slope=[]
    dum_list_coff_ed=[]
    dum_list_coff_only=[]
    for ri, run_date in enumerate(run_dates):
        if site in y_mod_sum_list[ri].keys():
            dum_list_sum.append(y_mod_sum_list[ri][site].rename(run_date))
            dum_list_sum2.append(y_mod_sum_list2[ri][site].rename(run_date))
            dum_list_ff.append(y_mod_ff_list[ri][site].rename(run_date))
            dum_list_edgar.append(y_mod_ff_list[ri][site + "_edgar"].rename(run_date))
        if site in y_mod_coff_list[ri].keys():
            dum_list_coff.append(y_mod_coff_list[ri][site].rename(run_date))
            dum_list_coff2.append(y_mod_coff_list2[ri][site].rename(run_date))
            dum_list_coff_ed.append(y_mod_coff_list[ri][site+ "_edgar"].rename(run_date))
        
            dum_list_coff_only.append(y_mod_coff_only_list[ri][site].rename(run_date))
        
            co_obs_dum = y_obs_co_list[ri][site]["obs"][run_date:end_dates[ri]]
            co_std_dum = y_obs_co_list[ri][site]["std"][run_date:end_dates[ri]]
            co_slope_dum = y_obs_co_list[ri][site]["slope"][run_date:end_dates[ri]]
            dum_list_co_obs.append(co_obs_dum)
            dum_list_co_std.append(co_std_dum)
            dum_list_co_slope.append(co_slope_dum)
            
    y_sum_temp =pd.concat(dum_list_sum, axis=1)
    y_mod_sum_site[site] = y_sum_temp.fillna(0).sum(axis=1) 
    
    y_sum_temp2 =pd.concat(dum_list_sum2, axis=1)
    y_mod_sum_site2[site] = y_sum_temp2.fillna(0).sum(axis=1) 
    
    y_ff_temp =pd.concat(dum_list_ff, axis=1)
    y_mod_ff_site[site] = y_ff_temp.fillna(0).sum(axis=1) 
    
    y_edgar_temp =pd.concat(dum_list_edgar, axis=1)
    y_mod_edgar_site[site] = y_edgar_temp.fillna(0).sum(axis=1) 
    
    y_coff_temp =pd.concat(dum_list_coff, axis=1)
    y_mod_co_site[site] = y_coff_temp.fillna(0).sum(axis=1) 
    
    y_coff_temp2 =pd.concat(dum_list_coff2, axis=1)
    y_mod_co_site2[site] = y_coff_temp2.fillna(0).sum(axis=1) 
    
    y_obs_co_site[site] =pd.concat(dum_list_co_obs, axis=0)
    y_std_co_site[site] =pd.concat(dum_list_co_std, axis=0)
    y_slope_co_site[site] =pd.concat(dum_list_co_slope, axis=0)
    
    y_coff_temp_ed =pd.concat(dum_list_coff_ed, axis=1)
    y_mod_co_edgar_site[site] = y_coff_temp_ed.fillna(0).sum(axis=1) 
  
    y_coff_only_temp =pd.concat(dum_list_coff_only, axis=1)
    y_mod_coff_site[site] = y_coff_only_temp.fillna(0).sum(axis=1) 
    
#%%

y_mod_co2_site={}
for site in sites:
    y_mod_co2_site = y_mod_sum_site[site] + y_bc_co2[site]
    y_mod_co2_site2 = y_mod_sum_site2[site] + y_bc_co2[site]

    df_site_out = pd.DataFrame()
    
    #df_site_out["co2_mod"] = y_mod_co2_site[site]

    
    df_site_out["co2ff_mod"] = y_mod_ff_site[site]
    df_site_out["co2ff_mod_edgar"] = y_mod_edgar_site[site]
    df_site_out["co2_obs"] = y_obs_co2[site]["obs"]
    df_site_out["co2_obs_std"] = y_obs_co2[site]["std"]
    df_site_out["co2_obs_slope"] = y_obs_co2[site]["slope"]
    
    #df_site_out["co_mod"] = y_mod_co_site[site]
    df_site_out["co_obs"] = y_obs_co_site[site]
    df_site_out["co_obs_std"] = y_std_co_site[site]
    df_site_out["co_obs_slope"] = y_slope_co_site[site]
    
    df_site_out["coff_mod" ] = y_mod_coff_site[site]
    
    
    # Correct model to fit observed monthly mean
    co2_mod = y_mod_co2_site
    co_mod = y_mod_co_site[site]
    co_mod_edgar = y_mod_co_edgar_site[site]
    
    co2_mod2 = y_mod_co2_site2
    co_mod2 = y_mod_co_site2[site]
    
    co2_mod_mn = co2_mod.groupby(pd.Grouper(freq='M')).transform('mean')
    co2_obs_mn = df_site_out["co2_obs"].groupby(pd.Grouper(freq='M')).transform('mean')
    diff_co2 = co2_obs_mn - co2_mod_mn
    
    co2_mod_mn2 = co2_mod2.groupby(pd.Grouper(freq='M')).transform('mean')
    diff2_co2 = co2_obs_mn - co2_mod_mn2
    
    co_mod_mn = co_mod.groupby(pd.Grouper(freq='M')).transform('mean')
    co_obs_mn = df_site_out["co_obs"].groupby(pd.Grouper(freq='M')).transform('mean')
    diff_co = co_obs_mn - co_mod_mn
    
    co_mod_mn2 = co_mod2.groupby(pd.Grouper(freq='M')).transform('mean')
    diff2_co = co_obs_mn - co_mod_mn2
    
    
    # Now add corrected model to dataset
    df_site_out["co2_mod"] = co2_mod + diff_co2
    df_site_out["co_mod"] = co_mod + diff_co
    
    df_site_out["co2_mod_perturb"] = co2_mod2 + diff2_co2
    df_site_out["co_mod_perturb"] = co_mod2 + diff2_co
    
    df_site_out["co_mod_edgar"] = co_mod_edgar + diff_co
    
    # Correct CO2 and CO model each month to fit observed background variation
    
    fname_out = out_dir + site + "_GC_output.csv"
    df_site_out.to_csv(fname_out)
