#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 11:08:38 2020

Code for CO2-only inversion

Read in site files and combine 

Try firstly solving for bulk CO2 (CO2 + CO2 ff together)

Then try solving separately, but make sure it works first.- i.e with pseudo data

@author: mlunt
"""
import numpy as np
import xarray
import matplotlib.pyplot as plt
import glob
#import re
import pandas as pd
import json
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colorbar as cbar
#from areagrid import areagrid
from dateutil.relativedelta import relativedelta
import ml_inv_mod as inv_mod
import subprocess
import argparse

def open_ds(fname):
    with xarray.open_dataset(fname) as ds:
        ds.load()
    return ds

def combine_site_obs(sites, obs_dir, 
                        start_date, end_date, daytime=False):
    """
    Combine the enesmble data for individual sites into one array
    
    """
    y_obs_list=[]
    y_time_list=[]
    y_site_list=[]
    site_obs_out={}
    for site in sites:
      
        run_str_site = obs_dir + site + "_co2ff_xgb.csv"  
        
        fname_ens = glob.glob(run_str_site)

        if len(fname_ens)<1:
            print ("Exiting loop, no data for " + site + " " + "co2ff")
            continue     #skip loopif site file doesn't exist
    
        df_site1_temp = pd.read_csv(fname_ens[0], index_col="time", parse_dates=True)
        df_site1 = df_site1_temp[start_date:end_date]
        
        if daytime == True:
            dummy1  = df_site1[(df_site1.index.hour>9) & (df_site1.index.hour < 18)]
            df_site1 = dummy1
            
        if len(df_site1)<1:
            print ("Exiting loop, no data for " + site + " " + "co2ff")
            continue    
            
        y_obs_list.append(df_site1["co2ff_obs"].values)
        y_time_list.append(df_site1.index)
        
        nobs_site = len(df_site1["co2ff_obs"].values)
        site_str = []
        for ti in range(nobs_site):
            site_str.append(site)
            
        y_site_list.append(np.asarray(site_str))
        
        site_obs_out[site]= df_site1
    
    y_dict_out={"obs": np.hstack(y_obs_list), 
                "time": np.hstack(y_time_list),
                "site": np.hstack(y_site_list)}
        
    return y_dict_out,site_obs_out

def combine_site_output(species, sites, run_dir, run_date, site_obs, 
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
        #df_site1 = df_site1_temp[start_date:end_date]
        
        #df_site1 = df_site1_temp.loc[site_obs[site].index]
        df_site1 = pd.merge(site_obs[site],df_site1_temp, how="inner", left_index=True, right_index=True)
        
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

#%%
#parser = argparse.ArgumentParser(description='This is a demo script by Mark.')              
#parser.add_argument("start", help="Start date string yyyymmdd") 
#parser.add_argument("end", help="End date string yyyymmdd") 
#parser.add_argument("version", help="Run version name")
#args = parser.parse_args()
#
#start_date = args.start
#end_date = args.end
#version = args.version


start_date = "20180701"
end_date = "20180731"
version = "test_mcmc_co_only"


data_root = "/home/mlunt/datastore/DARE/gc_output/inv_runs/"
run_groups = ["ff"]
#run_name1 = "test_CO2" 
#run_name2 = "test_CO2ff"

nbasis = 109
species= "CO2"
nBC=1


daytime=False
well_mixed=True
mixing_thresh= 0.2   # Permitted std. dev of model concs in lowest 5 model layers. 


scale_prior =False
monthly_ratios = np.asarray([1.199219, 1.140589, 1.073386, 1.012469, 0.918819, 0.817042, 0.768065,
       0.832761, 0.928363, 1.030775, 1.109421, 1.169091]) # Monthly ratios from UK monthly prior / UK prior annual mean 
month = int(start_date[4:6])
prior_scaling = 1./monthly_ratios[month-1]  #*1.5   # Make flat prior 50% larger

sites_in = ["GAT", "HTM",  
        "KRE",  "NOR", "OPE",
        "SAC" ,  "SVB", "TRN", 
        "TAC",  "BSD", "HFD"]

#%%
# Setup MCMC parameters here
#################
###############################
# Define MCMC parameters
para_temp=False   # 0 = off 1 = on

nIt = 5000
burn_in = 5000
nsub = 10

# Paramters for mcmc inversion
#################################
# Define nbeta

nbeta = 2   # Number of parallel chains - needs to be defined even if no Parallel tempering
series = np.linspace(0.,1.,num=nbeta)
beta_1 = np.exp(0 + (np.log(250)-0)*series)
beta= 1./beta_1

####################################################
# Define stepsizes
stepsize_sig_y_0 = 0.2  # units of ppb
stepsize_0 = 0.5
stepsize_pdf_p1_0 = 0.1
stepsize_pdf_p2_0 = 0.1
#####################################################
# Set up hyper-parameters

x_pdf = 3   # 3  = lognormal
pdf_param1_pdf = 2
pdf_param2_pdf=2

pdf_param1_0 = 1.
pdf_param2_0 = 0.5

pdf_p1_hparam1_0 = pdf_param1_0*1.
pdf_p1_hparam2_0 = 0.1
pdf_p2_hparam1_0 = pdf_param2_0*1.
pdf_p2_hparam2_0 = 0.2

#########################################################
# Sey up sigma_y hyper-parameters

#sigma_model_ap = 25.    # Was 5.
sigma_model_ap = 10.    # Was 5.

sigma_y_period = 31   # Length of each estimation period fro sigma_y in days. 
                     # The longer it is the quicker inversion will run
sigma_model_pdf = 1   # Unifrom
#sigma_model_pdf = 3   # lognormal

#sigma_model_hparam1_0 = sigma_model_ap*1.   # Units == ppb
#sigma_model_hparam2_0 = sigma_model_ap/5.   # units = ppb 

sigma_model_hparam1_0 = 0.5   # Units == ppb
sigma_model_hparam2_0 = 20.   # units = ppb   # Was 30

####################################################

#%%

#sites= ["TAC", "RGL", "BSD", "HFD"]
#sites=["NOR"]
inv_out_dir = "/home/mlunt/datastore/DARE/ml_Data/inv_outputs/"  + version + "/"

BC_names = ["bias"]
 

spc_names=[]
for xi in range(nbasis):
    #spc_names.append("CO2R" + str(xi+1).zfill(2) )
    spc_names.append(species + "_E" + str(xi+1))
 

spc_mod_dict1 ={}
spc_mod_dict2 ={}
for spc in spc_names:
    spc_mod_dict1[spc]=[]
    spc_mod_dict2[spc]=[]
 
# Set up run dates
#This assumes each assimilation window is 1 month
# Set up which dates to read in. Max 1 assim window and 1 lag windows

if start_date == "20170101":
    run_dates= ["20170101"]
else:    
    pd_run_date0  = pd.to_datetime(start_date) - relativedelta(months=1)
    run_date0 = pd_run_date0.strftime("%Y%m%d")
    run_dates   = [run_date0, start_date]
        
nruns = len(run_dates)

    
ngroups = len(run_groups)

#%%
# Read in data from ml_Data directory
obs_dir = "/home/mlunt/datastore/DARE/ml_data/xgb_output/"

y_data_dict, site_obs = combine_site_obs(sites_in, obs_dir, 
                        start_date, end_date, daytime=False)
sites_data = list(site_obs.keys())
# Need to use these data as where to inddex CO2 model somehow
# Migth be easier to save a monthly dataframe for each site and then pass that in to combine_site_output
#%%
H_basis_rdate={}
H_bc_rdate = {}
x_group_rdate={}
for run_date in run_dates:
   # H_basis_rdate[run_date] = []
    
    x_group_rdate[run_date]=[]
   
    
    
    for xi in range(nBC):
        x_group_rdate[run_date].append("bias")
        
        
    for run_str in run_groups:
        run_name=run_str + "_" + species + "_" + run_date   
        site_out_dir = data_root + run_name + "/"  
        
     
        # Read in fixed Basis_function names
        y_basis_stack, y_data_dum = combine_site_output(species, sites_data, site_out_dir, run_date, site_obs,
                        start_date, end_date, spc_names, 
                        daytime=daytime, read_obs=True)
        
        H_basis_rdate[run_date]=y_basis_stack
        
       
        for xi in range(nbasis):
            x_group_rdate[run_date].append(run_str)
        
    


#%%

# Aggregate model outputs into one Jacobian matrix 
nobs_unfilt = len(y_data_dict["obs"])
nstate = nbasis*nruns + nBC
nstate0 = nbasis

H_unfilt = np.zeros((nstate,nobs_unfilt))

if nBC == 1:
    H_unfilt[0,:] = 1.
    
for ri, run_date in enumerate(run_dates):
   
        H_unfilt[nBC + ri*(nstate0): nBC+ ri*(nstate0)+nbasis,:] = H_basis_rdate[run_date]
    
x_group_list=[]
for run_date in run_dates:
    x_group_list.append(x_group_rdate[run_date])
   
x_group = np.hstack(x_group_list)


#%%

# If filtering data by where model is well-mixed set to true
if well_mixed == True:
    y_slope = y_data_dum["slope"]  # slope is the std dev of total CO or CO2 species in the bottom 5 levels of the model
    wh_mixed = np.where(np.abs(y_slope) < mixing_thresh)[0]
    
    print(species, len(y_data_dict["obs"]))
    print("Well mixed:")
    print(species, len(wh_mixed))
    
    y_obs = y_data_dict["obs"][wh_mixed]
    y_std = np.abs(y_obs)*0.1
    y_time = y_data_dict["time"][wh_mixed]
    y_site = y_data_dict["site"][wh_mixed]
   
          
    H = np.transpose(H_unfilt[:,wh_mixed])
    
else:       
    y_obs = y_data_dict["obs"]
    y_std = np.abs(y_obs)*0.1
    y_time = y_data_dict["time"]
    y_site = y_data_dict["site"]
 
    H= np.transpose(H_unfilt)
    
nobs = len(y_obs)  

y_spc_list=[]
for ti in range(nobs):
    y_spc_list.append(species)
y_spc = np.asarray(y_spc_list)    

sites = np.unique(y_site)

# End of obs read section    
    


#%%
sigma_model = 0.5
sigma_y = np.sqrt(y_std**2 + sigma_model**2)

x_ap = np.zeros((nstate))+1.
x_ap[0] = 0.
x_uncert = x_ap.copy()*0.2
x_uncert[0] = 0.05


# Reset x_prior if scale prior is true
if scale_prior == True:
    wh_ff = np.where(x_group =="ff")
    x_ap[wh_ff] = prior_scaling
    

#x_uncert[nspc+1:]=0.05
y_mod_ap = np.dot(H, x_ap)

R_inv = np.diag(1./sigma_y**2)
P_inv = np.diag(1./x_uncert**2)
#%%
# Inversion code

HRH = np.dot(np.dot(np.transpose(H),R_inv),H) 
LHS = np.linalg.inv(HRH + P_inv)
HRy = np.dot(np.dot(np.transpose(H),R_inv),y_obs) 
Px = np.dot(P_inv,x_ap)
RHS = HRy + Px

x_post= np.dot(LHS,RHS)
y_post = np.dot(H,x_post)

x_post_uncert = np.sqrt(np.diag(LHS)) 

# Write an f2py mcmc routine to do this. 
# Include uncertainty in sigma_y. 
# Also include error correlations in distance and time
# Use code I wrote for CH4 inversions. That should work I think.


#%%
"""
###################################
# Call MCMC routine here
pdf_param1_ap = x_ap * pdf_param1_0
pdf_param2_ap = x_ap*0. + pdf_param2_0
pdf_p1_hparam1 = x_ap * pdf_p1_hparam1_0
pdf_p1_hparam2 = x_ap*0 + pdf_p1_hparam2_0
pdf_p2_hparam1 = pdf_param2_ap * 1.
pdf_p2_hparam2 = pdf_param2_ap * pdf_p2_hparam2_0

stepsize = (x_ap *0. +1)* stepsize_0
stepsize_pdf_p1 = pdf_param1_ap * stepsize_pdf_p1_0
stepsize_pdf_p2 = pdf_param2_ap * stepsize_pdf_p2_0

x_pdf_all = np.zeros((nstate), dtype=np.int8) + x_pdf

mcmc_out_dict = dare_hbmcmc_mod.run_mcmc(["CO"], start_date, end_date, nIC, nBC, ngroups, nstate0, nruns,
             sigma_y_period, sites, nIt, nsub, burn_in,
             beta, x_ap, H, y_obs, y_std,y_time,y_site,y_mod_ap, y_spc,
             x_pdf_all, pdf_param1_pdf, pdf_param2_pdf, pdf_param1_ap, pdf_param2_ap, 
             pdf_p1_hparam1, pdf_p1_hparam2,
             pdf_p2_hparam1, pdf_p2_hparam2,
             stepsize, stepsize_pdf_p1, stepsize_pdf_p2,
             sigma_model_pdf, sigma_model_ap, sigma_model_hparam1_0, sigma_model_hparam2_0,
             stepsize_sig_y_0, tau,
             para_temp=False
             )

#%%
x_it= mcmc_out_dict["x_it"]
y_post_it = mcmc_out_dict["y_it"]
sigma_y_it = mcmc_out_dict["sigma_y_it"]
x_post_mcmc = np.mean(x_it,axis=1)
y_post_mcmc = np.mean(y_post_it,axis=1)
pdf_param2_mean = np.mean(mcmc_out_dict["pdf_param2_it"], axis=1)
ratio_out = mcmc_out_dict["ratio"]

ratio_sig_y_out = mcmc_out_dict["ratio_sigma_y"]
"""

fig,ax = plt.subplots()

ax.plot(y_mod_ap, label="Model prior")
ax.plot(y_obs, 'o', markersize=2, label="ICOS and DECC obs")
ax.plot(y_post, label = "Model posterior")

#ax.plot(y_post_mcmc,  label = "MCMC posterior")
#ax.plot(y_obs, 'o', markersize=2, label="ICOS and DECC obs")

ax.set_ylabel("CO2 (ppm)")
ax.set_xlabel("Observation count")

leg = ax.legend()



#%%
"""
wh_cobc = np.where((x_group == "BC") | (x_group =="IC"))[0]
wh_cochem = np.where((x_group == "CONMVOC") | (x_group =="COCH4"))[0]
y_ap_cobc = np.dot(H[:,wh_cobc],x_ap[wh_cobc])
y_post_cobc = np.dot(H[:,wh_cobc],x_post_mcmc[wh_cobc])
y_ap_cochem = np.dot(H[:,wh_cochem],x_ap[wh_cochem])
y_post_cochem = np.dot(H[:,wh_cochem],x_post_mcmc[wh_cochem])
# Save outputs
# x_post
# x_uncert
# y_post

# Split into different groups

ds_out = xarray.Dataset({#"x_mean":(["regions"], x_post),
                         "x_it":(["regions", "nIt"], x_it),
                         
                         "y_post_it":(["nobs", "nIt"], y_post_it),
                         "y_mod_ap":(["nobs"], y_mod_ap),
                         "y_site": (["nobs"], y_site),                       
                         "y_obs":(["nobs"], y_obs),
                         #"sigma_y":(["nobs"], sigma_y),
                         "sigma_y_it":(["nobs", "nIt"], sigma_y_it),
                         "y_time": (["nobs"], y_time),   
                         "y_obs_std": (["nobs"], y_std),
                         #"x_uncert":(["regions"], x_post_uncert),
                         "x_ap_uncert":(["regions"], pdf_param2_mean),
                         "groups":(["ngroups"], run_groups),
                         "run_dates":(["nruns"], run_dates),
                         "species":(["nspecies"], ["CO"]),
                         "acc_ratio": (["regions"], ratio_out),
                         "sigma_y_ratio":  ratio_sig_y_out,
                         "sites": (["nsites"], sites),
                         "x_group": (["regions"], x_group),
                         "y_cobc_ap":(["nobs"], y_ap_cobc),
                         "y_cobc_post":(["nobs"], y_post_cobc),
                         "y_cochem_ap":(["nobs"], y_ap_cochem),
                         "y_cochem_post":(["nobs"], y_post_cochem),
                         "tau": tau,
                         "nIC":nIC,
                         "nBC": nBC,
                         "burn_in": burn_in,
                         "nsub": nsub,
                         "nchem": nchem
                         },
                         )
#ds_out.attrs["comment"] = "Output from 31/08/21 with lognormal sigma_model and prior 20ppb to eb consistent with NAME rjmcmc"

for ri, run_str in enumerate(run_groups):
    
    #x_group = x_post[(nBC + nbasis*ngroups)*(nruns-1) + nIC +nBC + nbasis*ri: 
    #    (nBC + nbasis*ngroups)*(nruns-1) + nIC +nBC + nbasis*(ri+1)]
    #x_uncert_group = x_post_uncert[(nBC + nbasis*ngroups)*(nruns-1) + nIC +nBC + nbasis*ri: 
    #    (nBC + nbasis*ngroups)*(nruns-1) + nIC +nBC + nbasis*(ri+1)]

    x_ap_group = x_ap[(nBC + nbasis*ngroups)*(nruns-1) + nIC +nBC +nchem+ nbasis*ri: 
        (nBC + nbasis*ngroups)*(nruns-1) + nIC +nBC + nchem + nbasis*(ri+1)]
    x_it_group = x_it[(nBC + nbasis*ngroups)*(nruns-1) + nIC +nBC + nchem + nbasis*ri: 
        (nBC + nbasis*ngroups)*(nruns-1) + nIC +nBC +nchem + nbasis*(ri+1),:]
        
    pdf_p2_group = pdf_param2_mean[(nBC + nbasis*ngroups)*(nruns-1) + nIC +nBC + nchem + nbasis*ri: 
        (nBC + nbasis*ngroups)*(nruns-1) + nIC +nBC +nchem + nbasis*(ri+1)]
        
    ds_out["x_it_" + species.lower() + run_str] = (("basis", "nIt"), x_it_group)
    ds_out["x_ap_" + species.lower() + run_str] = (("basis"), x_ap_group)
    
    ds_out["x_ap_uncert_" + species.lower() + run_str] = (("basis"), pdf_p2_group)
    
    ds_out["x_it_" + species.lower() + run_str].encoding['least_significant_digit'] = 3 
    ds_out["x_ap_" + species.lower() + run_str].encoding['least_significant_digit'] = 2 
    ds_out["x_ap_uncert_" + species.lower() + run_str].encoding['least_significant_digit'] = 2 
    #ds_out["x_post_" + run_str] = (("basis"), x_group)
    #ds_out["x_uncert_" + run_str] = (("basis"), x_uncert_group)

    # Also save out y_mod_ap_group to file
    # And y_post_group  - want to see how much each term is contributing to measurements at different sites


fname_out = inv_out_dir + "inv_out_"  + version + "_" +  start_date + ".nc"

if len(glob.glob(inv_out_dir)) ==0:
    subprocess.call(["mkdir", '-p', inv_out_dir])

for key in list(ds_out.keys()):    
    ds_out[key].encoding['zlib'] = True 

for key in ["y_post_it", "x_it", "y_mod_ap", "y_obs",
            "sigma_y_it", "x_ap_uncert", "acc_ratio", "sigma_y_ratio"]:
    ds_out[key].encoding['least_significant_digit'] = 3 

ds_out.to_netcdf(path=fname_out) 

# Write a separate script to calculate grid cell and country emissions and uncertainties
# Then plot up national emissions from that. and compare to EnKF results

#%%
# Assign x values to a basis function map and plot
"""

# Read in basis function file
fname_mask = "/home/mlunt/ceph/verify/model_settings/DARE_runs/inv_runs/masks/DARE_basis_functions_0.25x0.3125.nc"
ds_mask = open_ds(fname_mask)

mask_keys = []
for xi in range(nbasis):
    mask_keys.append("bf_"+ str(xi+1).zfill(2))

scale_map={}
for ri, run_str in enumerate(run_groups):
    scale_map[run_str] = ds_mask.bf_01[0,:,:]*0.

    for xi, reg in enumerate(mask_keys):
        scale_map[run_str] = scale_map[run_str] + ds_mask[reg][0,:,:]*(x_post[nBC + nbasis*ngroups*(nruns-1) 
                                                                        + nbasis*ri + xi])
    
#%%
json_dir = "/home/mlunt/programs/Python/DARE/inv_code/"
with open(json_dir + "dare_site_info.json") as f:
    site_info=json.load(f)

#%%

map_data = scale_map["ff"]
lon1 = ds_mask.lon.values
lat1 = ds_mask.lat.values
dlon1 = lon1[1]-lon1[0]
dlat1 = lat1[1]-lat1[0]

proj = ccrs.PlateCarree()
fig2,ax2=plt.subplots(subplot_kw=dict(projection=proj),figsize=(8,8))
h2 = ax2.pcolormesh(lon1-dlon1/2., lat1-dlat1/2., map_data, transform=ccrs.PlateCarree(), cmap='RdBu_r', vmin=0.5, vmax=1.5)
ax2.coastlines()
ax2.add_feature(cfeature.BORDERS)

for site in sites:
    
    dict_si = site_info[site]
    
    lat_si = dict_si["latitude"]
    lon_si = dict_si["longitude"]
    
    ax2.scatter(lon_si,lat_si, s=12, color='black')


cax2,kw2 = cbar.make_axes(ax2,location='bottom',pad=0.05,shrink=0.7)
out=fig2.colorbar(h2,cax=cax2,**kw2, label="CO2 flux (kg/m2/s)")

#area=areagrid(lat1,lon1)

#total_post = (emis_map*area).sum()*60.*60.*24.*30./1.e12
#total_ap = (emis_map_ap*area).sum()*60.*60.*24.*30./1.e12

