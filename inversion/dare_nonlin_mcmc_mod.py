#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 24 09:48:42 2021

MOdule file to set up MCMC code

Basically to keep line size down from actual run script

@author: mlunt
"""
import numpy as np
import xarray
import matplotlib.pyplot as plt
import pandas as pd
import json
import time as run_time
import test_nonlin as mcmc_mod
import datetime as dt
import glob

def get_nsigma_y(y_time, y_site, y_spc,start_date, end_date, nmeasure, 
                 bl_period,sites, species, nsites_spc):
      
    
    d0=pd.to_datetime(start_date)
    d1=pd.to_datetime(end_date)
    delta = d1 - d0
    ndays = delta.days + 1
    
    #bl_period = ndays/1.
    y_bl=np.zeros((nmeasure))
    
    nsigma=0
    nsigma_max = np.int(np.ceil(ndays/np.float(bl_period)))
    ntime_stn=np.zeros((nsites_spc))
    
    ydim1=0
    
    for spc in species:
    
        for si,site in enumerate(sites):
            
            wh_site_spc = np.where( (y_site==site) & (y_spc == spc))[0]
            if len(wh_site_spc >0):
                
                #mf_time_temp = y_time[y_site==site]
                mf_time_temp = y_time[wh_site_spc]
                
                #fp_data_H3 = fp_data_H[sites[si]].dropna("time", how="all")        
                nsigma_stn=0
                bl_start=d0
                #mf_time_temp=fp_data_H3.time.values
                #mf_time_temp=fp_data_H[sites[si]].time.values
                mf_time_temp2=pd.to_datetime(mf_time_temp)
                
                ntime_stn[si]=len(mf_time_temp)
                
                for ti in range(nsigma_max):
                    bl_end=bl_start+dt.timedelta(days=bl_period)
                    
                    wh=np.where(np.logical_and(mf_time_temp2>=bl_start,
                                               mf_time_temp2<bl_end))
                                                  
                    if len(wh[0]) > 0:
                        # This line is a problem - it assumes everything is sequential
                        #y_bl[wh+np.sum(ntime_stn[:si],dtype=np.uint16)]=nsigma_stn+nsigma
                        #nsigma_stn+=1
                        
                        y_bl[wh_site_spc[wh]] = nsigma_stn+nsigma
                        nsigma_stn+=1
                        
                    bl_start=bl_start+dt.timedelta(days=bl_period)
                    
                    n_obs = len(wh[0])
                    if n_obs > ydim1:
                        ydim1 = n_obs*1
            
                nsigma+=nsigma_stn
    
    # INDEX R
    #print(ydim1, nsigma)
    #return y_bl, ydim1,nsigma
    #return y_bl, ydim1, nsigma
    R_indices = np.zeros((ydim1,nsigma), dtype=np.uint16)
    for ii in range(nsigma):      
        wh_bl=np.where(y_bl == ii)
        nwh=len(wh_bl[0])
        R_indices[:nwh,ii]=wh_bl[0]+1
        if nwh < ydim1:
            R_indices[nwh:,ii]=np.max(wh_bl)+1
    
    ydim2=nsigma*1
    
    return R_indices, ydim1, ydim2 

#%%
def run_mcmc(species_list, start_date, end_date, nIC, nBC, ngroups, nstate0, nruns,
             bl_period, sites, nIt, nsub, burn_in, nobs_gas1, k_ef, ef_indices,
             beta, ad_ap, ef_ap, H, y_obs, y_std,y_time,y_site,y_mod_ap, y_spc,
             ad_pdf_all, ad_param1, ad_param2, stepsize_ad, 
             ef_pdf_all, ef_param1, ef_param2, stepsize_ef, 
             sigma_model_pdf, sigma_model_ap, sigma_model_hparam1_0, sigma_model_hparam2_0,
             stepsize_sig_y_0, 
             para_temp=False, 
             ): 


    """
    Function to define all arrays needed for MCMC from the input scalars
    
    Also calls mcmc_mod to do MCMC
    
    Just need to figure out how to account for different run_dates 
    and different species groups CO, Co2gee, reco ff etc. 
    
    Use a y_spc variable to separate out where species are the same.
    
    Sites needs to be a list of both sites[CO2] and sites[CO]
    
    To include R corrrelations between CO adn CO2 need to somehow redefine deltatime.
    Essentially define as 2x tau for where y_time_co  = y_time_co2 and y_site =
    
    """
    ########################################################
    # Define full x hyperparam vecotrs for size k
    
    nobs= len(y_obs)
    nstate = len(ad_ap)
    nbeta = len(beta)
    #nsites = len(sites)
    n0_ap = y_mod_ap - y_obs # This has to be this way round to be consistent with fortran script

    nobs_site_list = []
    nobs_list2 = []
    #for ss, si in enumerate(sites):
    #    wh_site = np.ravel(np.where(y_site == si))
    #    nobs_site[ss]=len(wh_site)
        
    for spc in species_list:
        for ss, si in enumerate(sites):
            wh_site = np.where((y_site == si) & (y_spc == spc))[0]
            
            if len(wh_site) >0:
                nobs_site_list.append(len(wh_site))
                    
    nobs_site = np.asarray(nobs_site_list)
    nsites_spc = len(nobs_site) # This should be legth of sites[CO2] + sites[CO]
      
    # Define R indices
    R_indices, ydim1, ydim2 = get_nsigma_y(y_time,y_site, y_spc, start_date, end_date, 
                                       nobs, bl_period,sites,species_list, nsites_spc)

    #out_dict={"deltatime": deltatime}
    #          "ydim1": ydim1,
    #          "ydim2": ydim2}
    #return out_dict

    ################################################
    # Define sigma model params
    sigma_model0=np.zeros((ydim2))
    sigma_model0[:]=sigma_model_ap  
    
    sigma_measure = y_std.copy()
    error_structure = sigma_measure*0.+1.
    
    sigma_model_hparam1 = np.zeros((ydim2))+sigma_model_hparam1_0
    sigma_model_hparam2 = np.zeros((ydim2))+sigma_model_hparam2_0
    stepsize_sigma_y = np.zeros((ydim2))  + stepsize_sig_y_0
    
    ########################################################
    # Set-up all variables for parallel tempering:
    ad=np.zeros((nstate,nbeta))
    ef=np.zeros((k_ef,nbeta))
    sigma_model = np.zeros((ydim2,nbeta))
    n0=np.zeros((nobs,nbeta))  
    for ib in range(nbeta):  
        ad[:,ib]=ad_ap.copy()  
        ef[:,ib]=ef_ap.copy()  
        n0[:,ib] = n0_ap.copy()
        sigma_model[:,ib]=sigma_model0.copy()

    sigma_model_pdf_all = np.zeros((ydim2), dtype=np.int8) + sigma_model_pdf
    
    # Call mcmc script
    nit_sub=nIt/nsub
    if para_temp == True:
        para_temp_in = 1
    else:
        para_temp_in = 0
        
        
#    out_dict={"R_indices": R_indices,
#              "nsites_spc": nsites_spc,
#              "x_ap": x_ap,
#              "pdf_param1_ap": pdf_param1_ap,
#              "pdf_param2_ap": pdf_param2_ap}
#    return out_dict
        
   
    #%%
    # MCMC version
    print ("Starting MCMC...")
    
    startt = run_time.time()
    
    # y_corr options
    ad_it, ef_it, y_post_it, sigma_model_it, sigma_y_it, \
    n0T_out, accept_ad, reject_ad, accept_ef, reject_ef, \
    accept_sigma_y, reject_sigma_y, accept_swap, reject_swap, \
    tot_acc_ad, tot_acc_ef, tot_acc_sigma_y, \
    accept_ad_all, reject_ad_all, accept_ef_all, reject_ef_all, \
    accept_sigma_y_all, reject_sigma_y_all = mcmc_mod.mcmc_uncorr.hbtdmcmc(
    beta, ad, ef, H,y_obs,n0, 
    ad_param1, ad_param2, ef_param1, ef_param2,  sigma_model, sigma_measure, 
    R_indices, sigma_model_hparam1, sigma_model_hparam2, stepsize_sigma_y, sigma_model_pdf, 
    para_temp_in, ad_pdf_all, ef_pdf_all, burn_in, 
    stepsize_ad, stepsize_ef, ef_indices, nobs_gas1, 
    nIt, nsub, nit_sub, 
    nbeta, nstate, k_ef, nobs,  ydim1, ydim2)    
        
    
    endt = run_time.time()
    print ("MCMC complete in ", (endt-startt), " seconds")
    
    ratio_ad = 1.*accept_ad/(accept_ad+reject_ad)
    ratio_ef = 1.*accept_ef/(accept_ef+reject_ef)
    ratio_sigma_y = accept_sigma_y/(accept_sigma_y+reject_sigma_y)
    #ratio_all = accept_all/(accept_all+reject_all)
    
    #ds_out = xarray.Dataset(#"x_it": (["nstate", "nIt"], x_it),
    #                        "y_it": (["nobs","nIt"], y_post_it)
                                     
    out_dict = {"ad_it": ad_it,
                "ef_it": ef_it, 
                "y_it": y_post_it,
                "sigma_model_it": sigma_model_it,
                "sigma_y_it": sigma_y_it,
                "ratio_ad": ratio_ad,
                "ratio_ef": ratio_ef,
                "ratio_sigma_y": ratio_sigma_y,
                "tot_acc_ad":tot_acc_ad,
                "tot_acc_ef":tot_acc_ef,

                }                                 
    
    return out_dict

def combine_site_output(species, sites, run_dir, run_date, 
                        start_date, end_date, keys, daytime=False, read_obs=True):
    """
    Combine the enesmble data for individual sites into one array
    
    """
    
    y_mod_all_dict={}
    for key in keys:
        y_mod_all_dict[key]=[]
    
    if read_obs == True:
        y_obs_list=[]
        y_std_list=[]
        y_time_list=[]
        y_site_list=[]
        y_lat_list=[]
        y_lon_list=[]
        y_slope_list=[]
    
    for site in sites:
      
        run_str_site = run_dir + "GC_output_" + species + "_" + site + "_" + run_date + ".csv" 
        
        fname_ens = glob.glob(run_str_site)

        if len(fname_ens)<1:
            print ("Exiting loop, no data for " + site + " " + species)
            continue     #skip loopif site file doesn't exist
    
        df_site1_temp = pd.read_csv(fname_ens[0], index_col="time", parse_dates=True)
        df_site1 = df_site1_temp[start_date:end_date]
        
        if daytime == True:
            dummy1  = df_site1[(df_site1.index.hour>9) & (df_site1.index.hour < 18)]
            df_site1 = dummy1
            
        for key in keys:
            
            if key in list(df_site1.keys()):
                y_mod_all_dict[key].append(df_site1[key].values)

        if read_obs == True:
            
            y_obs_list.append(df_site1["obs"].values)
            y_std_list.append(df_site1["std"].values)
            
            y_time_list.append(df_site1.index)
            
            nobs_site = len(df_site1["obs"].values)
            site_str = []
            for ti in range(nobs_site):
                site_str.append(site)
                
            y_site_list.append(np.asarray(site_str))
            
            y_lat_list.append(df_site1["lat"].values)
            y_lon_list.append(df_site1["lon"].values)
            
            y_slope_list.append(df_site1["slope"].values)
            

    y_mod_all_sites=[]
    for key in keys:
        if len(y_mod_all_dict[key])>0:
            y_mod_all_sites.append(np.hstack(y_mod_all_dict[key]))
        
    y_mod_all_out = np.stack(y_mod_all_sites)
    
    if read_obs == True:
        
        y_dict_out={"obs": np.hstack(y_obs_list),
                    "std": np.hstack(y_std_list),
                    "time": np.hstack(y_time_list),
                    "site": np.hstack(y_site_list),
                    "lat": np.hstack(y_lat_list),
                    "lon": np.hstack(y_lon_list),
                    "slope": np.hstack(y_slope_list)}
        
        return y_mod_all_out, y_dict_out
    
    return y_mod_all_out