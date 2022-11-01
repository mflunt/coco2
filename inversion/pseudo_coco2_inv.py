#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 11:08:38 2020

Code for CO2-only inversion

Read in site files and combine 

Try firstly solving for bulk CO2 (CO2 + CO2 ff together)

Then try solving separately, but make sure it works first.- i.e with pseudo data

Need to add MCMC.

Also need to add switches to either solve CO2 and Co as one. 
OR add offset CO term to inversion

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
#import inv_mod
import subprocess
import argparse
import dare_nonlin_mcmc_mod as dare_mod

def open_ds(fname):
    with xarray.open_dataset(fname) as ds:
        ds.load()
    return ds

#%%
parser = argparse.ArgumentParser(description='This is a demo script by Mark.')              
parser.add_argument("start", help="Start date string yyyymmdd") 
parser.add_argument("end", help="End date string yyyymmdd") 
parser.add_argument("version", help="Run version name")
args = parser.parse_args()

start_date = args.start
end_date = args.end
version = args.version


#start_date = "20180301"
#end_date = "20180331"
#version = "pseudo_coco2_corr1_sigma2"

coco2_corr = "weak"   # "perfect", "strong", "edgar" or "weak"

noise_co2 = 2.  # CO2 data noise in ppm
noise_co = 10. # CO data noise in ppb  2,5,or 10

data_root = "/home/mlunt/datastore/DARE/gc_output/inv_runs/"
#species = ["CO2","CO2","CO2","CO" ]
run_groups = ["gee", "reco", "ff"]
#run_name1 = "test_CO2" 
#run_name2 = "test_CO2ff"

species = "CO2"
nbasis = 109
nBC=4
nchem = 2
start0 = start_date

sigma_mod_spc={"CO2":2,
               "CO": 4}

co2_only=False
fixed_coco2_ratio=False
joint_bio=True

daytime=False
well_mixed=True
mixing_thresh= {"CO2":0.2,
                "CO": 1} # Permitted std. dev of model concs in lowest 5 model layers. 

scale_prior_bio = False
scale_prior = False
monthly_ratios = np.asarray([1.199219, 1.140589, 1.073386, 1.012469, 0.918819, 0.817042, 0.768065,
       0.832761, 0.928363, 1.030775, 1.109421, 1.169091]) # Monthly ratios from UK monthly prior / UK prior annual mean 
month = int(start_date[4:6])
#prior_scaling = 1./monthly_ratios[month-1] *1.2
prior_scaling = 1.2
prior_bio_scaling = 0.5

#sites= ["CMN", "GAT", "HPB", "HTM", "IPR", "JFJ",
#        "KIT", "KRE", "LIN", "NOR", "OPE", "PUY",
#        "SAC" , "SMR", "SVB", "TOH", "TRN", "UTO",
#        "TAC", "RGL", "BSD", "HFD"]

sites_in={}
sites_in["CO2"]= ["GAT",  "HTM", 
        "KIT", "KRE", "LIN", "NOR", "OPE", 
        "SAC" , "SVB", "TOH", "TRN", 
        "TAC", "RGL", "BSD", "HFD",
        "CBW", "MHD",   #"CRP", "MLH"
        "HUN", "BRM", "OHP"]

sites_in["CO"] = ["GAT", "HTM",  
        "KRE",  "NOR", "OPE",
        "SAC" ,  "SVB", "TRN", 
        "TAC",  "BSD", "HFD"]


# Setup MCMC parameters here
#################
###############################
# Define MCMC parameters
para_temp=False   # 0 = off 1 = on

nIt = 25000
burn_in = 25000
nsub = 50

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
stepsize_0 = 0.2

###################################################
# Set up hyper-parameters
ad_pdf = 3   # 3  = lognormal
ad_param1_0 = 1.
ad_param2_0 = 0.4

ef_pdf = 3   # 3  = lognormal
ef_param1_0 = 1.
ef_param2_0 = 0.2

########################################################
# Sey up sigma_y hyper-parameters

#sigma_model_ap = 25.    # Was 5.
sigma_model_ap = 2.    # Was 5.

sigma_y_period = 31   # Length of each estimation period fro sigma_y in days. 
                     # The longer it is the quicker inversion will run
sigma_model_pdf = 1   # Unifrom
#sigma_model_pdf = 3   # lognormal

#sigma_model_hparam1_0 = sigma_model_ap*1.   # Units == ppb
#sigma_model_hparam2_0 = sigma_model_ap/5.   # units = ppb 

sigma_model_hparam1_0 = 0.5   # Units == ppb
sigma_model_hparam2_0 = 20.   # units = ppb   # Was 30

##################################################
#%%

inv_out_dir = "/home/mlunt/datastore/DARE/inv_outputs/inv_runs/"  + version + "/"

BC_names = []
BC_names_co=[]
for xi in range(nBC):
    BC_names.append(species + "_BC" + str(xi+1)) 
    BC_names_co.append("CO" + "_BC" + str(xi+1)) 

BC_names_co.append("CONMVOC")
BC_names_co.append("COCH4")


spc_names=[]
spc_names_co = []
for xi in range(nbasis):
    #spc_names.append("CO2R" + str(xi+1).zfill(2) )
    spc_names.append(species + "_E" + str(xi+1))
    spc_names_co.append("CO" + "_E" + str(xi+1))
 

spc_mod_dict1 ={}
spc_mod_dict2 ={}
for spc in spc_names:
    spc_mod_dict1[spc]=[]
    spc_mod_dict2[spc]=[]
 
# Set up run dates
#This assumes each assimilation window is 1 month
# Set up which dates to read in. Max 1 assim window and 1 lag windows

if start_date == start0:
    run_dates= [start_date]
else:    
    pd_run_date0  = pd.to_datetime(start_date) - relativedelta(months=1)
    run_date0 = pd_run_date0.strftime("%Y%m%d")
    run_dates   = [run_date0, start_date]
        
nruns = len(run_dates)

if start0 in run_dates:
    nIC = 2
else:
    nIC=0
    
ngroups = len(run_groups)

#%%

H_basis_co2_rdate={}
H_bc_co2_rdate = {}
x_group_rdate={}
H_basis_co_rdate={}
H_bc_co_rdate = {}
for run_date in run_dates:
   
    H_basis_co2_rdate[run_date] = []
    #H_basis_co_rdate[run_date] = []
    
    if run_date == start0:
        x_group_rdate[run_date]=["IC_CO2", "IC_CO"]
    else:
        x_group_rdate[run_date]=[]
    #H_bc_run_date[run_date] = []
    
    
    for xi in range(nBC):
        x_group_rdate[run_date].append("BC_"+species)
        
    for xi in range(nBC):
        x_group_rdate[run_date].append("BC_"+"CO")
        
    x_group_rdate[run_date].append("CONMVOC")
    x_group_rdate[run_date].append("COCH4")
        
    for run_str in run_groups:
        run_name=run_str + "_" + species + "_" + run_date   
        site_out_dir = data_root + run_name + "/"  
        
        if run_str == "ff":
            if run_date == start_date:
                y_bc_stack, y_data_dict = dare_mod.combine_site_output(species, sites_in[species], site_out_dir, run_date, 
                            start_date, end_date, BC_names, 
                            daytime=daytime, read_obs=True)
            else:
                y_bc_stack = dare_mod.combine_site_output(species, sites_in[species], site_out_dir, run_date, 
                            start_date, end_date, BC_names, 
                            daytime=daytime, read_obs=False)
                
            H_bc_co2_rdate[run_date] = y_bc_stack
            
        # Read in fixed Basis_function names
        y_basis_stack = dare_mod.combine_site_output(species, sites_in[species], site_out_dir, run_date, 
                        start_date, end_date, spc_names, 
                        daytime=daytime, read_obs=False)
        
        H_basis_co2_rdate[run_date].append(y_basis_stack)
        
        
        if joint_bio == True:
            if run_str =="gee":
                for xi in range(nbasis):
                    x_group_rdate[run_date].append("bio")
            elif run_str == "ff":
                for xi in range(nbasis):
                    x_group_rdate[run_date].append("ff")
        else:
            for xi in range(nbasis):
                x_group_rdate[run_date].append(run_str)
        
#        for xi in range(nbasis):
#            x_group_rdate[run_date].append(run_str)
        
    # Read in IC term if needed - outside of run_str loop
    run_name_co="ff_CO_" + run_date   
    site_out_dir_co = data_root + run_name_co + "/" 
    if run_date == start0:
        y_ic_co2_stack = dare_mod.combine_site_output(species, sites_in[species], site_out_dir, run_date, 
                    start_date, end_date, [species + "IC"], 
                    daytime=daytime, read_obs=False)
        
        y_IC_co2 = y_ic_co2_stack.squeeze()
        
        y_ic_co_stack = dare_mod.combine_site_output("CO", sites_in["CO"], site_out_dir_co, run_date, 
                    start_date, end_date, ["COIC"], 
                    daytime=daytime, read_obs=False)
        
        y_IC_co = y_ic_co_stack.squeeze()
        
    # Also read in CO files
    if run_date == start_date:
        y_bc_co_stack, y_co_data_dict = dare_mod.combine_site_output("CO", sites_in["CO"], site_out_dir_co, run_date, 
                            start_date, end_date, BC_names_co, 
                            daytime=daytime, read_obs=True)
    else:
        y_bc_co_stack = dare_mod.combine_site_output("CO", sites_in["CO"], site_out_dir_co, run_date, 
                            start_date, end_date, BC_names_co, 
                            daytime=daytime, read_obs=False)
                
    H_bc_co_rdate[run_date] = y_bc_co_stack
    
    # Read in fixed Basis_function names
    y_basis_co_stack = dare_mod.combine_site_output("CO", sites_in["CO"], site_out_dir_co, run_date, 
                    start_date, end_date, spc_names_co, 
                    daytime=daytime, read_obs=False)
    
    H_basis_co_rdate[run_date] = y_basis_co_stack
    
    #if fixed_coco2_ratio==False:
    #    for xi in range(nbasis):
    #        x_group_rdate[run_date].append("coff")
    
    #for xi in range(nbasis):
    #    x_group_rdate[run_date].append("COff")

#%%    

if joint_bio ==True:
    run_groups = ["bio", "ff"]
    ngroups = len(run_groups)

# Aggregate model outputs into one Jacobian matrix 
        
nobs_unfilt_co2 = len(y_data_dict["obs"]) 
nobs_unfilt_co = len(y_co_data_dict["obs"]) 
nobs_unfilt = nobs_unfilt_co2 + nobs_unfilt_co


nstate = (nBC*2 + nchem + nbasis*(ngroups))*nruns + nIC
nstate0 = nBC*2+nchem +nbasis*(ngroups)

H_unfilt = np.zeros((nstate,nobs_unfilt))


if nIC > 0:
    H_unfilt[0,:nobs_unfilt_co2] = y_IC_co2
    H_unfilt[1,nobs_unfilt_co2:] = y_IC_co
    
for ri, run_date in enumerate(run_dates):
    H_unfilt[nIC + ri*(nstate0): nIC + ri*(nstate0)+nBC,:nobs_unfilt_co2] = H_bc_co2_rdate[run_date]
    H_unfilt[nIC + ri*(nstate0)+nBC: nIC + ri*(nstate0)+nBC*2 +nchem, nobs_unfilt_co2:] = H_bc_co_rdate[run_date]
    
    for xi in range(ngroups):
        
        if joint_bio==True:
            if xi ==0:
                H_unfilt[nIC + ri*(nstate0)+nBC*2 + nchem + nbasis*xi:
                nIC+ ri*(nstate0)+nBC*2 + nchem +nbasis*(xi+1),:nobs_unfilt_co2] = (
                    H_basis_co2_rdate[run_date][xi] + H_basis_co2_rdate[run_date][xi+1])
            elif xi ==1:
                H_unfilt[nIC + ri*(nstate0)+nBC*2 + nchem + nbasis*xi:
                nIC+ ri*(nstate0)+nBC*2 + nchem+nbasis*(xi+1), :nobs_unfilt_co2] = H_basis_co2_rdate[run_date][xi+1]
        else:    
            H_unfilt[nIC + ri*(nstate0)+nBC*2 +nchem + nbasis*xi: nIC+ ri*(nstate0)+nBC*2 + nchem +nbasis*(xi+1),
                 :nobs_unfilt_co2] = H_basis_co2_rdate[run_date][xi]
        
        
        #H_unfilt[nIC + ri*(nstate0)+nBC*2 +nchem + nbasis*xi: nIC+ ri*(nstate0)+nBC*2 + nchem +nbasis*(xi+1),
        #         :nobs_unfilt_co2] = H_basis_co2_rdate[run_date][xi]
    if co2_only == False:
        H_unfilt[nIC + ri*(nstate0)+nBC*2 +nchem + nbasis*(ngroups-1): nIC+ ri*(nstate0)+nBC*2 + nchem +nbasis*ngroups,
             nobs_unfilt_co2:] = H_basis_co_rdate[run_date]
    
#    H_unfilt[nIC + ri*(nstate0)+nBC*2 +nchem + nbasis*(ngroups-1): nIC+ ri*(nstate0)+nBC*2 + nchem +nbasis*ngroups,
#             nobs_unfilt_co2:] = H_basis_co_rdate[run_date]
    
    #if fixed_coco2_ratio ==False:
    #    H_unfilt[nIC + ri*(nstate0)+nBC*2 + nchem + nbasis*ngroups: nIC+ ri*(nstate0)+nBC*2 + nchem +nbasis*(ngroups+1),
    #         nobs_unfilt_co2:] = H_basis_co_rdate[run_date]
    
x_group_list=[]
for run_date in run_dates:
    x_group_list.append(x_group_rdate[run_date])
   
x_group = np.hstack(x_group_list)
#%%

# Concat co2 and co dataset

y_obs_temp = np.hstack([y_data_dict["obs"],y_co_data_dict["obs"]]) 
y_std_temp = np.hstack([y_data_dict["std"],y_co_data_dict["std"]]) 
y_time_temp = np.hstack([y_data_dict["time"],y_co_data_dict["time"]]) 
y_site_temp = np.hstack([y_data_dict["site"],y_co_data_dict["site"]]) 

# If filtering data by where model is well-mixed set to true
if well_mixed == True:
    y_slope_co2 = y_data_dict["slope"]  # slope is the std dev of total CO or CO2 species in the bottom 5 levels of the model
    wh_mixed_co2 = np.where(np.abs(y_slope_co2) < mixing_thresh["CO2"])[0]
    
    y_slope_co = y_co_data_dict["slope"]  # slope is the std dev of total CO or CO2 species in the bottom 5 levels of the model
    wh_mixed_co = np.where(np.abs(y_slope_co) < mixing_thresh["CO"])[0]
    
    #print(species, len(y_data_dict["obs"]))
    #print("Well mixed:")
    #print(species, len(wh_mixed))
    
    wh_mixed = np.hstack([wh_mixed_co2, wh_mixed_co+nobs_unfilt_co2])
    
    y_obs = y_obs_temp[wh_mixed]
    y_std = y_std_temp[wh_mixed]
    y_time = y_time_temp[wh_mixed]
    y_site = y_site_temp[wh_mixed]
    #y_lat_co2 = y_data_dict["lat"][wh_mixed_co2]
    #y_lon_co2 = y_data_dict["lon"][wh_mixed_co2]
    
    if start_date == start0:
        y_IC = y_ic_co2_stack.squeeze()[wh_mixed_co2]
      
    H = np.transpose(H_unfilt[:,wh_mixed])
    
    nobs_co2 = len(wh_mixed_co2)
    nobs_co = len(wh_mixed_co)
    
else:       
    y_obs = y_obs_temp.copy()
    y_std = y_std_temp.copy()
    y_time = y_time_temp.copy()
    y_site = y_site_temp.copy()
    #y_lat = y_data_dict["lat"]
    #y_lon = y_data_dict["lon"]
    
    H= np.transpose(H_unfilt)
    
    nobs_co2 = nobs_unfilt_co2
    nobs_co = nobs_unfilt_co
    
nobs = len(y_obs)  

y_spc_list=[]
for ti in range(nobs_co2):
    y_spc_list.append(species)
for ti in range(nobs_co):
    y_spc_list.append("CO")
    
y_spc = np.asarray(y_spc_list)  

sites = np.unique(y_site)

# End of obs read section    

#%%  

#H_group={}
#y_mod_group={}
#for xi, run_str in enumerate(run_prefixes): 
#    
#    y_mod_group[run_str] = np.zeros((nobs))
#    for run_date in run_dates:
#        y_mod_group[run_str] = y_mod_group[run_str] + H_basis_rdate[run_date][xi][:,wh_mixed].sum(axis=0)

    #y_mod_group[run_str] = np.sum(H_group[run_str],axis=1)   
#y_mod_ff = np.sum(H_ff,axis=1)    
#y_mod_ap = np.sum(H,axis=1) 

#%%
sigma_model = y_std.copy()*0.
sigma_model[:nobs_co2] = sigma_mod_spc["CO2"]
sigma_model[nobs_co2:] = sigma_mod_spc["CO"]
#sigma_y = np.sqrt(y_std**2 + sigma_model**2)

#wh_coff = np.where(x_group == "coff")[0]
wh_ff = np.where(x_group =="ff")[0]
wh_co2gee = np.where(x_group =="gee")[0]
wh_co2reco = np.where(x_group =="reco")[0]


ad_ap = np.zeros((nstate))+1.

####################
# Define EF dimensions and prior values
k_ef = len(wh_ff)   # Define ef in terms of this (since only applies to ff terms)
ef_indices = wh_ff + 1  # Add 1 for Fortran
ef_ap = np.zeros((k_ef))+1.
####################

# Reset x_prior if scale prior is true
if scale_prior == True:
    #wh_ff = np.where(x_group =="ff")
    ad_ap[wh_ff] = prior_scaling


if scale_prior_bio == True:
    wh_bio = np.where((x_group =="reco") | (x_group =="gee") | (x_group == "bio"))
    ad_ap[wh_bio] = prior_bio_scaling

y_mod_ap = np.dot(H,ad_ap)


#%%
# Set up AD and EF pdf parameters and stepsizes
ad_pdf_param1 = ad_ap * ad_param1_0
ad_pdf_param2 = ad_ap*0. + ad_param2_0
stepsize_ad = (ad_ap *0. +1)* stepsize_0
ad_pdf_all = np.zeros((nstate), dtype=np.int8) + ad_pdf

ef_pdf_param1 = ef_ap * ef_param1_0
ef_pdf_param2 = ef_ap*0. + ef_param2_0
stepsize_ef = (ef_ap *0. +1)* stepsize_0
ef_pdf_all = np.zeros((k_ef), dtype=np.int8) + ef_pdf

if fixed_coco2_ratio == True:
    stepsize_ef = stepsize_ef*0.
    ef_pdf_param2 = ef_pdf_param2*0.001

# Reduce uncertainty on COff -i.e. force more of a positive correlation between CO2ff and COffs
#pdf_param2_ap[wh_coff] =  0.1
#pdf_p2_hparam1[wh_coff] = 0.1

# Increase uncertainty on COff -i.e. force more of no correlation between CO2ff and COffs
#pdf_param2_ap[wh_coff] =  2.
#pdf_p2_hparam1[wh_coff] = 2.

if joint_bio == True:
    # Increase uncertainty on joint bio term to make it a fair comparison to gee/reco approach
    wh_bio = np.where(x_group == "bio")[0]
    ad_pdf_param2[wh_bio] = ad_param2_0*4.
    ad_pdf_all[wh_bio] = 2

#%%
# Set up y_true and y_obs
# Set up pseudo obs

fname_ratio  = "/home/mlunt/datastore/EDGAR/edgar_tno_bf_ratios.csv"
df_ratio = pd.read_csv(fname_ratio)
x_ratio_ff = df_ratio["co2_ratio"]
x_ratio_coff = df_ratio["co_ratio"]

seed = 10
np.random.seed(seed)  # Fix the seed for repeatability

x_true = np.zeros((nstate)) +1.
#wh_ff = np.where(x_group =="ff")[0]
#x_true[wh_ff] = 1.2
#rand_ff = np.random.normal(loc=1., scale = ad_param2_0, size = nbasis)
rand_gee = np.random.normal(loc=1., scale = ad_param2_0, size = nbasis)
rand_reco = np.random.normal(loc=1., scale = ad_param2_0, size = nbasis)

rand_ff = x_ratio_ff.copy()

if joint_bio == False:
    wh_reco = np.where(x_group =="reco")[0]
    wh_gee = np.where(x_group =="gee")[0]
    x_true[wh_reco] = rand_reco
    x_true[wh_gee] = rand_gee
else:
    wh_bio = np.where(x_group =="bio")[0]
    #x_true[wh_bio] = rand_reco
    
    H_bio_sum = np.sum(H[:nobs_co2,wh_bio], axis=1)
    if H_bio_sum.mean()>=0:
        x_true[wh_bio] = rand_reco + 1.1
    else:
        x_true[wh_bio] = (rand_reco-1)*-1 + 0.9
    
#x_true[wh_ff] = x_ratio_ff.values
x_true[wh_ff] = rand_ff
x_true[x_true<0] = 0.2

# Use realistic noise values - get from inversion outputs
# co2 rmse = 3 ppm maybe use 2 ppm
# co rmse = 8 ppb  # maybe use 5 ppb

y_true = np.dot(H, x_true)

y_co2ff = np.dot(H[:,wh_ff], rand_ff)
y_co2bio = np.dot(H[:,wh_bio], x_true[wh_bio])

y_co2ff_noise = np.random.normal(loc=0., scale=np.abs(y_co2ff)*0.2, size = nobs)
y_co2bio_noise = np.random.normal(loc=0., scale=np.abs(y_co2bio)*0.5, size = nobs)
y_noise = y_co2ff_noise + y_co2bio_noise

#y_noise = np.random.normal(loc=0., scale=noise_co2, size = nobs)
y_obs = y_true + y_noise

# This allows me to separate co treu from CO2ff true
H_co_temp = np.zeros((nobs_unfilt_co, (nBC+ nchem + nbasis)*nruns + nIC//2))
if nIC >0:
    H_co_temp[:,0] = y_IC_co

H_co_temp[:, nIC//2:nBC + nIC//2+nchem] = np.transpose(H_bc_co_rdate[run_date])
H_co_temp[:, nIC//2+nBC+nchem:] = np.transpose(H_basis_co_rdate[run_date])
H_co = H_co_temp[wh_mixed_co,:]

x_true_co = np.zeros((nbasis+nchem+nBC+nIC//2)) +1.
#x_true2_co = np.zeros((nbasis+nchem+nBC+nIC//2)) +1.
#x_true_co[-nbasis:] = 1.1  + np.random.normal(loc=0, scale=0.05, size=nbasis)  # Was 1.2

# This code for strong CO:CO2 correlations
if coco2_corr == "strong":
    x_true_co[-nbasis:] = (rand_ff-1)*0.8 +1. + np.random.normal(loc=0, scale=0.2, size=nbasis)

elif coco2_corr == "perfect":
    # This code for EDDGAR CO ratios
    x_true_co[-nbasis:] = rand_ff.copy()

elif coco2_corr == "edgar":
    # This code for EDDGAR CO ratios
    x_true_co[-nbasis:] = x_ratio_coff.copy()

else:  
    # Use this code ofr poorly correlated CO:CO2 errors
    x_true_co[-nbasis:] = (rand_ff-1.)*0.25 + 1. + np.random.normal(loc=0, scale=ef_param2_0, size=nbasis)
    #x_true_co[-nbasis:] = (x_ratio_ff-1.)*0.25 + 1. + np.random.normal(loc=0, scale=ef_param2_0, size=nbasis)

y_true_co = np.dot(H_co,x_true_co)

y_noise_co = np.random.normal(loc=0., scale=noise_co, size = nobs_co)  # low_noise = 0.2; normal noise = 2.

y_obs_co = y_true_co + y_noise_co

y_obs[nobs_co2:] = y_obs_co

corr_ff = np.corrcoef(rand_ff, x_true_co[-nbasis:])
#corr_ff2 = np.corrcoef(rand_ff, x_true2_co[-nbasis:])


#%%
# Need to make x[co2ff] have gaussian pdf
# Also check the uncertaintites are correctly defined
mcmc_out_dict = dare_mod.run_mcmc(["CO2", "CO"], start_date, end_date, nIC, nBC, ngroups, nstate0, nruns,
             sigma_y_period, sites, nIt, nsub, burn_in, nobs_co2,
             k_ef, ef_indices,
             beta, ad_ap, ef_ap, H, y_obs, y_std,y_time,y_site,y_mod_ap, y_spc,
             ad_pdf_all, ad_pdf_param1, ad_pdf_param2, stepsize_ad,
             ef_pdf_all, ef_pdf_param1, ef_pdf_param2, stepsize_ef,
             sigma_model_pdf, sigma_model_ap, sigma_model_hparam1_0, sigma_model_hparam2_0,
             stepsize_sig_y_0, 
             para_temp=False
             )

#%%
ad_it= mcmc_out_dict["ad_it"]
ef_it = mcmc_out_dict["ef_it"]
y_post_it = mcmc_out_dict["y_it"]
sigma_y_it = mcmc_out_dict["sigma_y_it"]
#x_post_mcmc = np.mean(x_it,axis=1)
y_post_mcmc = np.mean(y_post_it,axis=1)
ratio_ad = mcmc_out_dict["ratio_ad"]
ratio_ef = mcmc_out_dict["ratio_ef"]
ratio_sig_y_out = mcmc_out_dict["ratio_sigma_y"]

ad_post_mean = ad_it.mean(axis=1)
#%%

#fig,ax = plt.subplots()
#
#ax.plot(y_mod_ap, label="Model prior")
#ax.plot(y_obs, 'o', markersize=2, label="ICOS and DECC obs")
###ax.plot(y_post, label = "Model posterior")
##ax.plot(y_obs, 'o', markersize=2, label="ICOS and DECC obs")
#ax.plot(y_post_mcmc,  label = "MCMC posterior")
#ax.set_ylabel("CO2 (ppm)")
#ax.set_xlabel("Observation count")

#leg = ax.legend()

#%%
# calculate gains
#gain_all = 1-np.sum(np.abs(x_post_mcmc - x_true))/np.sum(np.abs(x_ap-x_true))

wh_ff03 = np.where((ad_it.std(axis=1)<=0.4) & (x_group == "ff"))
#wh_reco03 = np.where((ad_it.std(axis=1)<=0.4) & (x_group == "reco"))
#wh_gee03 = np.where((ad_it.std(axis=1)<=0.4) & (x_group == "gee"))

gain_ff = 1-np.sum(np.abs(ad_post_mean[wh_ff] - x_true[wh_ff]))/np.sum(np.abs(ad_ap[wh_ff]-x_true[wh_ff]))
#gain_reco = 1-np.sum(np.abs(ad_post_mean[wh_reco03] - x_true[wh_reco03]))/np.sum(np.abs(ad_ap[wh_reco03]-x_true[wh_reco03]))
#gain_gee = 1-np.sum(np.abs(ad_post_mean[wh_gee03] - x_true[wh_gee03]))/np.sum(np.abs(ad_ap[wh_gee03]-x_true[wh_gee03]))

# FInd CO_chem indices and outptu prior and post y_cochem
wh_cochem = np.where((x_group == "CONMVOC") | (x_group =="COCH4"))[0]
wh_cobc = np.where((x_group == "BC_CO") | (x_group =="IC_CO"))[0]
wh_co2bc = np.where((x_group == "BC_CO2") | (x_group =="IC_CO2"))[0]

y_ap_cobc = np.dot(H[:,wh_cobc],ad_ap[wh_cobc])
y_post_cobc = np.dot(H[:,wh_cobc],ad_post_mean[wh_cobc])
y_ap_cochem = np.dot(H[:,wh_cochem],ad_ap[wh_cochem])
y_post_cochem = np.dot(H[:,wh_cochem],ad_post_mean[wh_cochem])

y_ap_co2bc = np.dot(H[:,wh_co2bc],ad_ap[wh_co2bc])
y_post_co2bc = np.dot(H[:,wh_co2bc],ad_post_mean[wh_co2bc])
#%%

# Save outputs
# x_post
# x_uncert
# y_post

# Split into different groups

ds_out = xarray.Dataset({#"x_mean":(["regions"], x_post),
                         "x_it":(["regions", "nIt"], ad_it),
                         
                         "y_post_it":(["nobs", "nIt"], y_post_it),
                         "y_mod_ap":(["nobs"], y_mod_ap),
                         "y_site": (["nobs"], y_site),                       
                         "y_obs":(["nobs"], y_obs),
                         "y_true":(["nobs"], y_true),
                         #"sigma_y":(["nobs"], sigma_y),
                         "sigma_y_it":(["nobs", "nIt"], sigma_y_it),
                         "y_time": (["nobs"], y_time),   
                         "y_obs_std": (["nobs"], y_std),
                         #"x_uncert":(["regions"], x_post_uncert),
                         "groups":(["ngroups"], run_groups),
                         "run_dates":(["nruns"], run_dates),
                         "species":(["nspecies"], ["CO2"]),
                         "acc_ratio_ad": (["regions"], ratio_ad),
                         "sigma_y_ratio":  ratio_sig_y_out,
                         "sites": (["nsites"], sites),
                         "x_group": (["regions"], x_group),
                         "y_co2bc_ap":(["nobs"], y_ap_co2bc),
                         "y_co2bc_post":(["nobs"], y_post_co2bc),
                         "y_cobc_ap":(["nobs"], y_ap_cobc),
                         "y_cobc_post":(["nobs"], y_post_cobc),
                         "y_cochem_ap":(["nobs"], y_ap_cochem),
                         "y_cochem_post":(["nobs"], y_post_cochem),
                         
                         "nIC":nIC,
                         "nBC": nBC,
                         "nchem": nchem,
                         "burn_in": burn_in,
                         "nsub": nsub,
                         "nobs_co2": nobs_co2,
                         "coco2_corr": corr_ff[0,1],

                         "noise_co": noise_co,
                         "noise_co2": noise_co2
                         },
                         )
#ds_out.attrs["comment"] = "Output from 31/08/21 with lognormal sigma_model and prior 20ppb to eb consistent with NAME rjmcmc"

x_true_coff = x_true_co[-nbasis:]
ds_out["x_true_" + "coff"] = (("basis"), x_true_coff)
for ri, run_str in enumerate(run_groups):
    
    #x_group = x_post[(nBC + nbasis*ngroups)*(nruns-1) + nIC +nBC + nbasis*ri: 
    #    (nBC + nbasis*ngroups)*(nruns-1) + nIC +nBC + nbasis*(ri+1)]
    #x_uncert_group = x_post_uncert[(nBC + nbasis*ngroups)*(nruns-1) + nIC +nBC + nbasis*ri: 
    #    (nBC + nbasis*ngroups)*(nruns-1) + nIC +nBC + nbasis*(ri+1)]

    x_ap_group = ad_ap[nstate0*(nruns-1) + nIC +nBC*2 + nchem + nbasis*ri: 
        nstate0*(nruns-1) + nIC +nBC*2 + nchem + nbasis*(ri+1)]

    x_true_group = x_true[nstate0*(nruns-1) + nIC +nBC*2 + nchem + nbasis*ri: 
        nstate0*(nruns-1) + nIC +nBC*2 + nchem + nbasis*(ri+1)]
    x_it_group = ad_it[nstate0*(nruns-1) + nIC +nBC*2 + nchem + nbasis*ri: 
        nstate0*(nruns-1) + nIC +nBC*2 + nchem + nbasis*(ri+1),:]
        
    pdf_p2_group = ad_pdf_param2[nstate0*(nruns-1) + nIC +nBC*2 + nchem + nbasis*ri: 
        nstate0*(nruns-1) + nIC +nBC*2 + nchem + nbasis*(ri+1)]
        
    ds_out["x_it_" + species.lower() + run_str] = (("basis", "nIt"), x_it_group)
    ds_out["x_ap_" + species.lower() + run_str] = (("basis"), x_ap_group)
    ds_out["x_true_" + species.lower() + run_str] = (("basis"), x_true_group)
    
    ds_out["x_ap_uncert_" + species.lower() + run_str] = (("basis"), pdf_p2_group)
    
    ds_out["x_it_" + species.lower() + run_str].encoding['least_significant_digit'] = 3 
    ds_out["x_ap_" + species.lower() + run_str].encoding['least_significant_digit'] = 2 
    ds_out["x_true_" + species.lower() + run_str].encoding['least_significant_digit'] = 2 
    ds_out["x_ap_uncert_" + species.lower() + run_str].encoding['least_significant_digit'] = 2 
    #ds_out["x_post_" + run_str] = (("basis"), x_group)
    #ds_out["x_uncert_" + run_str] = (("basis"), x_uncert_group)


    if run_str == "ff":
        #if fixed_coco2_ratio==True:
        #    ds_out["x_it_" + "coff"] =  (("basis", "nIt"), x_it_group)
        #    ds_out["x_ap_" + "coff"] = (("basis"), x_ap_group)
        #    ds_out["x_ap_uncert_" + "coff"] = (("basis"), pdf_p2_group)
        #else:
            # This assumes that coff comes after ff
            #x_ap_coff = ad_ap[nstate0*(nruns-1) + nIC +nBC*2 +nchem + nbasis*(ri+1):]
            ef_ap_coff = ef_ap[nbasis*(nruns-1):]
            # (nBC*2 + nchem + nbasis*ngroups)*(nruns-1) + nIC +nBC*2 + nchem + nbasis*(ri+2)]
            #ad_it_coff = ad_it[nstate0*(nruns-1) + nIC +nBC*2 +nchem + nbasis*(ri+1):,:]
            ef_it_coff = ef_it[nbasis*(nruns-1):,:]
            #ad_p2_coff = ad_pdf_param2[nstate0*(nruns-1) + nIC +nBC*2 +nchem + nbasis*(ri+1):]
            ef_p2_coff = ef_pdf_param2[nbasis*(nruns-1):]

            coff_p2_ap = (ef_p2_coff/ef_ap_coff + pdf_p2_group/x_ap_group)*x_ap_group*ef_ap_coff
       
            ds_out["x_it_" + "coff"] =  (("basis", "nIt"), x_it_group*ef_it_coff)
            ds_out["x_ap_" + "coff"] = (("basis"), x_ap_group*ef_ap_coff)
            ds_out["x_ap_uncert_" + "coff"] = (("basis"), coff_p2_ap)
            ds_out["ef_it_coff"] = (("basis","nIt"), ef_it_coff)
            

#corr_out=np.zeros((nbasis))
#for xi in range(nbasis):
#    corr_out[xi] = np.corrcoef(ds_out.x_it_co2ff[xi,:], ds_out.x_it_coff[0,:])[0,1]


fname_out = inv_out_dir + "inv_out_"  + version + "_" +  start_date + ".nc"

if len(glob.glob(inv_out_dir)) ==0:
    subprocess.call(["mkdir", '-p', inv_out_dir])

for key in list(ds_out.keys()):    
    ds_out[key].encoding['zlib'] = True 

for key in ["y_post_it", "x_it", "y_mod_ap", "y_obs",
            "sigma_y_it", "acc_ratio_ad", "sigma_y_ratio"]:
    ds_out[key].encoding['least_significant_digit'] = 3 

ds_out.to_netcdf(path=fname_out) 

