#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 16:49:15 2021

Create posterior scaling maps, emission maps and country totals from 
DARE inv outputs

Do for each group member (e.g. gee, reco, ff)

Plot outputs and compare to Enkf outputs in separate script

Need to calculate prior uncertainties as well based on std

Add total formation from COchem (CH4 and VOCs) into outputs

Is this in HEMCO?

@author: mlunt
"""

import xarray
import numpy as np
import matplotlib.pyplot as plt
import regionmask
import areagrid
import pandas as pd

def open_ds(fname):
    with xarray.open_dataset(fname) as ds:
        ds.load()
    return ds

#version  = "test_mcmc_coco2bio2"
#version  = "test_mcmc_coco2_50pc_uncert"
#version  = "test_mcmc_co2_50pc_uncert"

#version  = "test_mcmc_co_only"

#version  = "test_mcmc_coco2_seq"
#version  = "test_mcmc_coco2_rcorr"
#version  = "test_mcmc_coco2_v2"
#version  = "test_mcmc_co_only_mix05"

#version  = "test_mcmc_coco2_coff01"

#version  = "nonlin_coco2"
#version  = "nonlin_coco2_fixed"
#version  = "nonlin_co2_only"
#version  = "lin_coco2"
#version  = "nonlin_co2_only_flat_ap"
#version  = "nonlin_coco2_flat_ap"
#version  = "lin_co2_only"
#version  = "lin_co2_only_5pc_uncert"

#version = "nonlin_coco2_v2"  # v2 ignores TAC
#version = "nonlin_coco2_v3"  # v3 with new DECC data
#version = "nonlin_coco2_fixed_v3"  # v3 with new DECC data
#version = "lin_co2_only_v3"  # v3 with new DECC data
#version = "nonlin_coco2_v4"  # v4 has stricter mixing threshold requirements
#version = "nonlin_coco2_daytime"
#version = "nonlin_coco2_bio_sigma"
#version = "nonlin_coco2_fix_gee"

#version = "co2_only_cams_BC"
#version = "co2_only_cams_BC_ff10"
#version = "co2_only_cams_BC_nonuk"
#version = "co2_only_cams_BC_daytime"

#version = "co2_only_inland3"
#version = "co2_only_new_mix_noIRE"
#version = "paper_co2_only"
#version = "paper_coco2"

#version = "co2_only_ocean"
version = "co2_only_newBC"

#version = "coco2_20pcap"
#version = "coco2_inland"


#species=["CO2"]
species=["CO2", "CO"]
#species=["CO"]
#run_dates= [ "20180301", "20180401"]

#run_dates = ["20180101", "20180201"]#,
       # "20180601", "20180701", ]#, "20180201", "20180301"]

run_dates = ["20180101", "20180201", "20180301",
             "20180401", "20180501", "20180601",
             "20180701", "20180801", "20180901", 
             "20181001", "20181101", "20181201"]  

#run_dates = ["20170101", "20170201", "20170301",
#             "20170401", "20170501", "20170601",
#             "20170701", "20170801", "20170901", 
#             "20171001", "20171101", "20171201",    
#             "20180101", "20180201", "20180301",
#             "20180401", "20180501", "20180601",
#             "20180701", "20180801", "20180901", 
#             "20181001", "20181101", "20181201",
#             "20190101", "20190201", "20190301",
#             "20190401", "20190501", "20190601",
#             "20190701", "20190801", "20190901", 
#             "20191001", "20191101", "20191201"] 

#spc_groups = {"CO2":["gee", "reco", "ff"],
#              "CO": ["ff"]}

spc_groups = {"CO2":["bio", "ff"],
              "CO": ["ff"]}

run_groups=[]
for spc in species: 
    group_list = spc_groups[spc]
    for group in group_list:
        run_groups.append(spc.lower() + group)
    
hemco_root = "/home/mlunt/ceph/verify/model_settings/DARE_runs/inv_runs/"
post_root = "/home/mlunt/datastore/DARE/inv_outputs/inv_runs/"
inv_out_dir = post_root + version + "/"


countries = ["United Kingdom", "Ireland", "France", "Germany",
             "Belgium", "Netherlands", "Switzerland", "Italy",
             "Denmark", "Norway", "Sweden", "Poland"]

pd_rdates = pd.to_datetime(run_dates)

#%%

# 1. Read in basis function file
fname_mask = "/home/mlunt/ceph/verify/model_settings/DARE_runs/inv_runs/masks/DARE_basis_functions_0.25x0.3125.nc"
ds_mask = open_ds(fname_mask)

nbasis=109
mask_keys = []
for xi in range(nbasis):
    mask_keys.append("bf_"+ str(xi+1).zfill(2))

emis_keys={}
for spc in species:
    emis_keys[spc]=[]

for spc in species:
    for xi in range(nbasis):
        emis_keys[spc].append(spc + "_E"+ str(xi+1))

#%%
# 2. 
# Loop through each run date. Retrieve the prior emissions field for each group
emis_ap_list={}
emis_ap_uncert_list={}
emis_post_it_list={}
#emis_uncert_list={}
scale_post_it_list={}
bf_post_it_list={}
#scale_uncert_list={}

bf_corr_coco2_list=[]

for group in run_groups:
    emis_ap_list[group]=[]
    emis_ap_uncert_list[group]=[]
    emis_post_it_list[group]=[]
    #emis_uncert_list[group]=[]
    scale_post_it_list[group]=[]
    bf_post_it_list[group]=[]
    #scale_uncert_list[group]=[]

coch4_it_list=[]
covoc_it_list=[]
coch4_ap_list=[]
covoc_ap_list=[]

for run_date in run_dates:
    
    fname_inv_out = inv_out_dir + "inv_out_" + version + "_" + run_date + ".nc" 
    # Read in output file
    ds_inv_out  = open_ds(fname_inv_out)
    
    emis_ap={}
    emis_ap_uncert={}
    emis_post_it={}
    scale_factor_it={}
    bf_post_it={}
    
    #for xi, reg in enumerate(emis_keys[spc]):
    #    basis_corr_co_co2[xi] = np.corrcoef(ds_inv_out["x_it_coff"][xi,:], ds_inv_out["x_it_co2ff"][xi,:])
    
    for spc in species:
    
        for group in spc_groups[spc]:
    
            if group == "bio":
                run_name1 =  "gee" + "_" + spc + "_" + run_date
                hemco_dir1 = hemco_root + run_name1 + "/OutputDir/window_output/"
                fname_hco_diag1 = hemco_dir1 + "HEMCO_diagnostics." + run_date + "0000.nc"
                ds_hemco = open_ds(fname_hco_diag1)
                
                run_name2 =  "reco" + "_" + spc + "_" + run_date
                hemco_dir2 = hemco_root + run_name2 + "/OutputDir/window_output/"
                fname_hco_diag2 = hemco_dir2 + "HEMCO_diagnostics." + run_date + "0000.nc"
                ds_hemco2 = open_ds(fname_hco_diag2)
            
            else:
                run_name =  group + "_" + spc + "_" + run_date
                hemco_dir = hemco_root + run_name + "/OutputDir/window_output/"
                fname_hco_diag = hemco_dir + "HEMCO_diagnostics." + run_date + "0000.nc"
                ds_hemco = open_ds(fname_hco_diag)
    
            emis_ap[group] = ds_hemco[emis_keys[spc][0]]*0.
            emis_ap_uncert[group] = ds_hemco[emis_keys[spc][0]]*0.
            emis_post_it[group] = ds_hemco[emis_keys[spc][0]]*0.
            #emis_uncert[group] = ds_hemco[emis_keys[0]]*0.
            #scale_factor_it[group] = ds_mask[mask_keys[0]]*0.
    
            bf_post_it[group] = ds_inv_out["x_it_" + spc.lower() + group]*0.
            # Now loop through each basis function
            # Map out hte prior and posterior emissions
            for xi, reg in enumerate(emis_keys[spc]):
                
                if group == "bio":
                    emis_bf = ds_hemco[reg] + ds_hemco2[reg]
                else:
                    emis_bf = ds_hemco[reg]
                #emis_bf = ds_hemco[reg]
                mask_bf = ds_mask[mask_keys[xi]].squeeze()
                emis_ap[group] = emis_ap[group] + emis_bf*ds_inv_out["x_ap_" + spc.lower() +group][xi]
                emis_ap_uncert[group] = emis_ap_uncert[group] + emis_bf*ds_inv_out["x_ap_uncert_" + spc.lower() +group][xi]
                emis_post_it[group] = emis_post_it[group] + emis_bf*ds_inv_out["x_it_" + spc.lower() + group][xi,:]
                #emis_uncert[group] = emis_post[group] + emis_bf*ds_inv_out["x_uncert_" + group][xi]
                #scale_factor_it[group] = scale_factor_it[group] + mask_bf*ds_inv_out["x_it_" + spc.lower() + group][xi,:]
                
                if xi ==0:
                    area_bf = areagrid.areagrid(emis_bf.lat.values, emis_bf.lon.values)
                emis_bf_area = emis_bf*area_bf 
                bf_post_it[group][xi,:] = emis_bf_area.sum()*ds_inv_out["x_it_" + spc.lower() + group][xi,:]
                
            scale_factor_it[group] = ds_inv_out["x_it_" + spc.lower() + group]
                
        # need to combine each grouped run_date map into a list I guess
    
            emis_ap_list[spc.lower() + group].append(emis_ap[group])   
            emis_ap_uncert_list[spc.lower() + group].append(emis_ap_uncert[group])   
            emis_post_it_list[spc.lower() +group].append(emis_post_it[group])  
            #emis_uncert_list[group].append(emis_uncert[group])  
            scale_post_it_list[spc.lower() + group].append(scale_factor_it[group])
            
            bf_post_it_list[spc.lower() +group].append(bf_post_it[group])  
           
        if spc == "CO":
            
            x_group = ds_inv_out.x_group
            wh_coch4 = np.where(x_group == "COCH4")[0][-1]
            wh_covoc = np.where(x_group == "CONMVOC")[0][-1]
            
            coch4_post_it = ds_hemco["COCH4"]*ds_inv_out["x_it"][wh_coch4,:]
            covoc_post_it = ds_hemco["CONMVOC"]*ds_inv_out["x_it"][wh_covoc,:]
            
            coch4_it_list.append(coch4_post_it)
            covoc_it_list.append(covoc_post_it)
            coch4_ap_list.append(ds_hemco["COCH4"])
            covoc_ap_list.append(ds_hemco["CONMVOC"])
    
    # For COco2 inversion calculate posterior correlations between bfs 
    # Caluclate btoh CO:CO2ff corrs as wel as CO2ff:CO2bio
    #But for bio will need to sum emissions for each bf to calculate it.
    #bf_corr_coco2 = np.zeros((nbasis))
    #for xi in range(nbasis):
    #    bf_corr_coco2[xi] = np.corrcoef(ds_inv_out["x_it_co2ff"][xi,:], ds_inv_out["x_it_coff"][xi,:])[0,1]
    
    #bf_corr_coco2_list.append(bf_corr_coco2)
# Create xarray dataarrays for each group
    
    
da_emis_ap={}
da_emis_ap_uncert={}
da_emis_post_it={}
da_scale_post_it={}
da_bf_post_it={}
for group in run_groups:
    da_emis_ap[group] = xarray.concat(emis_ap_list[group], dim="time")
    da_emis_ap_uncert[group] = xarray.concat(emis_ap_uncert_list[group], dim="time")
    da_emis_post_it[group] = xarray.concat(emis_post_it_list[group], dim="time")
    da_scale_post_it[group] = xarray.concat(scale_post_it_list[group], dim="time")
    da_bf_post_it[group] = xarray.concat(bf_post_it_list[group], dim="time")

da_coch4_ap = xarray.concat(coch4_ap_list, dim="time")
da_covoc_ap = xarray.concat(covoc_ap_list, dim="time")
da_coch4_it = xarray.concat(coch4_it_list, dim="time")
da_covoc_it = xarray.concat(covoc_it_list, dim="time")


bf_corr_coco2 = np.zeros((len(run_dates), nbasis))
bf_corr_bioff = np.zeros((len(run_dates), nbasis))
for ti in range(len(run_dates)):
    for xi in range(nbasis):
        if "CO" in species and "CO2" in species:
            bf_corr_coco2[ti,xi] = np.corrcoef(da_bf_post_it["co2ff"][ti,xi,:],da_bf_post_it["coff"][ti,xi,:] )[0,1]
        
        if "CO2" in species:
            if "co2gee" in run_groups:
                dum_bio = da_bf_post_it["co2reco"][ti,xi,:] + da_bf_post_it["co2gee"][ti,xi,:]
            else:
                dum_bio = da_bf_post_it["co2bio"][ti,xi,:]
            bf_corr_bioff[ti,xi] = np.corrcoef(da_bf_post_it["co2ff"][ti,xi,:],dum_bio )[0,1]

bf_std={}
bf_mean={}
bf_emis_mean={}
for group in run_groups:
    bf_std[group] = da_scale_post_it[group].std(dim="nIt")
    bf_mean[group] = da_scale_post_it[group].mean(dim="nIt")    
    bf_emis_mean[group] = da_bf_post_it[group].mean(dim="nIt")*60*60*24*365/1.e9    # kg/s to Tg/yr
#bf_corr_coco2 = xarray.concat(bf_corr_coco2_list, dim="time")
#%%

t_convert = 60*60*24*365.

lon_state = ds_hemco.lon.values
lat_state = ds_hemco.lat.values

# 4. With prior and posterior maps, calculate the national prior and posterior emissions and uncertainties    
mask = regionmask.defined_regions.natural_earth.countries_110.mask(lon_state, lat_state, xarray=True)
names = regionmask.defined_regions.natural_earth.countries_110.names

names2 = np.asarray(names)

wh_country={}
for country in countries:
    c_idx = np.where(names2 == country)
    basis_cnt = mask.where(mask == c_idx)
    basis_cnt_v = np.ravel(basis_cnt)
    wh_cnt = np.where(np.isfinite(basis_cnt_v))
    wh_country[country] = wh_cnt[0]
    

area = areagrid.areagrid(lat_state,lon_state)
# Turn area into an xarray dataarray for multiplication
da_area = xarray.DataArray(area, coords=[lat_state,lon_state], dims=["lat", "lon"])

#for ti,run_date in enumerate(run_dates):
    
emis_cnt_ap={}
uncert_cnt_ap={}
emis_cnt_it={}
#scale_cnt_it={}

nIt = len(ds_inv_out.nIt)

ncountries = len(countries)
ntime = len(run_dates)
for group in run_groups:
    emis_cnt_ap[group] = xarray.DataArray(np.zeros((ncountries,ntime)), 
               coords=[countries, pd_rdates], dims=["country", "time"])
    uncert_cnt_ap[group] = xarray.DataArray(np.zeros((ncountries,ntime)), 
               coords=[countries, pd_rdates], dims=["country", "time"])
    
    
    emis_cnt_it[group] = xarray.DataArray(np.zeros((ncountries,ntime, nIt)),
               coords=[countries, pd_rdates, np.arange(nIt)], dims=["country", "time", "nIt"])
    #scale_cnt_it[group] = xarray.DataArray(np.zeros((ncountries,ntime, nIt)),
    #            coords=[countries, pd_rdates, np.arange(nIt)], dims=["country", "time", "nIt"])
 
    
emis_cnt_ap_coch4 = xarray.DataArray(np.zeros((ncountries,ntime)), 
               coords=[countries, pd_rdates], dims=["country", "time"])
emis_cnt_mean_coch4 = xarray.DataArray(np.zeros((ncountries,ntime)), 
               coords=[countries, pd_rdates], dims=["country", "time"])

emis_cnt_ap_covoc = xarray.DataArray(np.zeros((ncountries,ntime)), 
               coords=[countries, pd_rdates], dims=["country", "time"])
emis_cnt_mean_covoc = xarray.DataArray(np.zeros((ncountries,ntime)), 
               coords=[countries, pd_rdates], dims=["country", "time"])
    
for ci,country in enumerate(countries):
        
    #wh = wh_country[country]
    c_idx = np.where(names2 == country)
    
    for group in run_groups:
        emis_cnt_dum1 = da_emis_ap[group].where(mask==c_idx)
        emis_cnt_dum2 = da_emis_post_it[group].where(mask==c_idx)
        uncert_cnt_dum1 = da_emis_ap_uncert[group].where(mask==c_idx)
    
        emis_cnt_ap[group][ci,:] = (emis_cnt_dum1*da_area).sum(dim=["lat","lon"])*t_convert/1.e9 # Tg/yr
        uncert_cnt_ap[group][ci,:] = (uncert_cnt_dum1*da_area).sum(dim=["lat","lon"])*t_convert/1.e9 # Tg/yr
        emis_cnt_it[group][ci,:,:] = (emis_cnt_dum2*da_area).sum(dim=["lat","lon"])*t_convert/1.e9 # Tg/yr
        
        emis_cnt_coch4_dum1 = da_coch4_ap.where(mask==c_idx)
        emis_cnt_coch4_dum2 = (da_coch4_it.mean(dim="nIt")).where(mask==c_idx)
        
        emis_cnt_covoc_dum1 = da_covoc_ap.where(mask==c_idx)
        emis_cnt_covoc_dum2 = (da_covoc_it.mean(dim="nIt")).where(mask==c_idx)
        
        emis_cnt_ap_coch4[ci,:] = (emis_cnt_coch4_dum1*da_area).sum(dim=["lat","lon"])*t_convert/1.e9 
        emis_cnt_mean_coch4[ci,:] = (emis_cnt_coch4_dum2*da_area).sum(dim=["lat","lon"])*t_convert/1.e9 
        
        emis_cnt_ap_covoc[ci,:] = (emis_cnt_covoc_dum1*da_area).sum(dim=["lat","lon"])*t_convert/1.e9 
        emis_cnt_mean_covoc[ci,:] = (emis_cnt_covoc_dum2*da_area).sum(dim=["lat","lon"])*t_convert/1.e9 
        #scale_cnt_it[group][ci,:,:] = np.sqrt(((scale_cnt_dum2*da_area)**2).sum(dim=["lat","lon"]))
    
# Calculate correlation coefficients between groups for each country
        # i.e. want to see how independent they are.
        
corr_gee_reco = np.zeros((ncountries, ntime)) 
corr_gee_ff = np.zeros((ncountries, ntime)) 
corr_reco_ff = np.zeros((ncountries, ntime))    
corr_bio_ff = np.zeros((ncountries, ntime))  
corr_coff_ff = np.zeros((ncountries, ntime))        

for ci,country in enumerate(countries):
    
    for ti in range(ntime):
        
        if "co2gee" in run_groups and "co2reco" in run_groups:
            corr_gee_reco[ci,ti] = np.corrcoef(emis_cnt_it["co2gee"][ci,ti,:], emis_cnt_it["co2reco"][ci,ti,:])[0,1]
            corr_gee_ff[ci,ti] = np.corrcoef(emis_cnt_it["co2gee"][ci,ti,:], emis_cnt_it["co2ff"][ci,ti,:])[0,1]
            corr_reco_ff[ci,ti] = np.corrcoef(emis_cnt_it["co2reco"][ci,ti,:], emis_cnt_it["co2ff"][ci,ti,:])[0,1]
            
            bio_cnt_it = emis_cnt_it["co2gee"] + emis_cnt_it["co2reco"]
            corr_bio_ff[ci,ti] = np.corrcoef(bio_cnt_it[ci,ti,:], emis_cnt_it["co2ff"][ci,ti,:])[0,1]
        elif "co2bio" in run_groups:
            corr_bio_ff[ci,ti] = np.corrcoef(emis_cnt_it["co2bio"][ci,ti,:], emis_cnt_it["co2ff"][ci,ti,:])[0,1]
            
        if "CO" in species and "CO2" in species:
            corr_coff_ff[ci,ti] = np.corrcoef(emis_cnt_it["coff"][ci,ti,:], emis_cnt_it["co2ff"][ci,ti,:])[0,1]
            
            
#%%
# Calculate COchem totals for the domain

coch4_ap =  (da_coch4_ap*da_area).sum(dim=["lat","lon"])*t_convert/1.e9  
covoc_ap =  (da_covoc_ap*da_area).sum(dim=["lat","lon"])*t_convert/1.e9  

coch4_it =  (da_coch4_it*da_area).sum(dim=["lat","lon"])*t_convert/1.e9  
covoc_it =  (da_covoc_it*da_area).sum(dim=["lat","lon"])*t_convert/1.e9         
            
coch4_mean = coch4_it.mean(dim="nIt")
covoc_mean = covoc_it.mean(dim="nIt")

coch4_an_it = coch4_it.groupby('time.year').mean(dim="time")
covoc_an_it = covoc_it.groupby('time.year').mean(dim="time")

coch4_pc_list=[]
covoc_pc_list=[]
coch4_pc_an_list=[]
covoc_pc_an_list=[]
for pc in [0.025, 0.16, 0.84,0.975]:
    coch4_pc_dum = coch4_it.quantile(pc,dim="nIt")
    coch4_pc_an_dum  = coch4_an_it.quantile(pc,dim="nIt")
    
    covoc_pc_dum = covoc_it.quantile(pc,dim="nIt")
    covoc_pc_an_dum  = covoc_an_it.quantile(pc,dim="nIt")
    
    coch4_pc_list.append(coch4_pc_dum)
    coch4_pc_an_list.append(coch4_pc_an_dum)
    
    covoc_pc_list.append(covoc_pc_dum)
    covoc_pc_an_list.append(covoc_pc_an_dum)
        
coch4_pc = xarray.concat(coch4_pc_list, dim="quantile")
coch4_pc_an = xarray.concat(coch4_pc_an_list, dim="quantile")

covoc_pc = xarray.concat(coch4_pc_list, dim="quantile")
covoc_pc_an = xarray.concat(coch4_pc_an_list, dim="quantile")

for ci,country in enumerate(countries):
        
    #wh = wh_country[country]
    c_idx = np.where(names2 == country)

#%%

# Calculate percentiles to save output space. i.e just save mean and 2.5,16,84 and 97.5 percentiles
emis_cnt_mean={}
scale_cnt_mean={}
emis_cnt_pc={}
#scale_cnt_pc={}

emis_map_mean={}
scale_map_mean={}

emis_map_pc={}
#scale_map_pc={}
for group in run_groups:
    
    emis_cnt_mean[group] = emis_cnt_it[group].mean(dim="nIt") 
    #scale_cnt_mean[group] = scale_cnt_it[group].mean(dim="nIt") 
    emis_map_mean[group]  = da_emis_post_it[group].mean(dim="nIt")
    scale_map_mean[group]  = da_scale_post_it[group].mean(dim="nIt")
    
    emis_cnt_pc_list=[]
    scale_cnt_pc_list=[]
    for pc in [0.025, 0.16, 0.84,0.975]:
    #for pc in [2.5, 16, 84,97.5]:        
        
        # Need to work out how to combiine into one array
        #emis_cnt_pc_dum  = np.percentile(emis_cnt_it[group],pc,axis=2)
        emis_cnt_pc_dum  = emis_cnt_it[group].quantile(pc,dim="nIt")
        #scale_cnt_pc_dum  = scale_cnt_it[group].quantile(pc,dim="nIt")
    
        emis_cnt_pc_list.append(emis_cnt_pc_dum)
        #scale_cnt_pc_list.append(scale_cnt_pc_dum)
        
        #emis_cnt_pc[group]  = emis_cnt_it[group].quantile(pc,dim="nIt") 
        #scale_cnt_pc[group]  = scale_cnt_it[group].quantile(pc,dim="nIt") 
        
        #emis_map_pc[group] = da_emis_post_it[group].quantile(pc, dim="nIt")

    emis_cnt_pc[group] = xarray.concat(emis_cnt_pc_list, dim="quantile")
    #scale_cnt_pc[group] = xarray.concat(scale_cnt_pc_list, dim="percentile")


#%%
# calculate Annual means and annual percentiles based on emis_cnt_it
emis_cnt_pc_an={}
for group in run_groups:   
    emis_cnt_it_an = emis_cnt_it[group].groupby('time.year').mean(dim="time")
    emis_cnt_pc_an_list=[]
    for pc in [0.025, 0.16, 0.84,0.975]:
        emis_cnt_pc_an_dum  = emis_cnt_it_an.quantile(pc,dim="nIt")
        emis_cnt_pc_an_list.append(emis_cnt_pc_an_dum)

    emis_cnt_pc_an[group] = xarray.concat(emis_cnt_pc_an_list, dim="quantile")
    
#%%
# Also should add reco and gee together in this script to calculate bio uncertainties
# Should do the same for CO2sum as well so I have correct uncertainties for these values both annual and monthly
co2sum_cnt_pc_list=[]
co2sum_cnt_pc_an_list=[]    
if "co2gee" in run_groups:
    emis_cnt_it_co2bio = emis_cnt_it["co2gee"] + emis_cnt_it["co2reco"]
    emis_cnt_it_co2sum = emis_cnt_it["co2gee"] + emis_cnt_it["co2reco"] + emis_cnt_it["co2ff"]
    
    emis_cnt_mean_co2bio = emis_cnt_it_co2bio.mean(dim="nIt")
    emis_cnt_mean_co2sum = emis_cnt_it_co2sum.mean(dim="nIt")
    
    emis_cnt_it_an_co2bio = emis_cnt_it_co2bio.groupby('time.year').mean(dim="time")
    emis_cnt_it_an_co2sum = emis_cnt_it_co2sum.groupby('time.year').mean(dim="time")
    
    co2bio_cnt_pc_list=[]
    co2bio_cnt_pc_an_list=[]
    for pc in [0.025, 0.16, 0.84,0.975]:
        co2bio_cnt_pc_dum = emis_cnt_it_co2bio.quantile(pc,dim="nIt")
        co2sum_cnt_pc_dum = emis_cnt_it_co2sum.quantile(pc,dim="nIt")
        
        co2bio_cnt_pc_an_dum  = emis_cnt_it_an_co2bio.quantile(pc,dim="nIt")
        co2sum_cnt_pc_an_dum  = emis_cnt_it_an_co2sum.quantile(pc,dim="nIt")
        
        co2bio_cnt_pc_list.append(co2bio_cnt_pc_dum)
        co2sum_cnt_pc_list.append(co2sum_cnt_pc_dum)
        
        co2bio_cnt_pc_an_list.append(co2bio_cnt_pc_an_dum)
        co2sum_cnt_pc_an_list.append(co2sum_cnt_pc_an_dum)
        
    emis_cnt_pc_co2bio = xarray.concat(co2bio_cnt_pc_list, dim="quantile")
    emis_cnt_pc_co2sum = xarray.concat(co2sum_cnt_pc_list, dim="quantile")
    
    emis_cnt_pc_an_co2bio = xarray.concat(co2bio_cnt_pc_an_list, dim="quantile")
    emis_cnt_pc_an_co2sum = xarray.concat(co2sum_cnt_pc_an_list, dim="quantile")
    
elif "CO2" in species:
    emis_cnt_it_co2sum = emis_cnt_it["co2bio"] + emis_cnt_it["co2ff"]
    emis_cnt_mean_co2sum = emis_cnt_it_co2sum.mean(dim="nIt")
    emis_cnt_it_an_co2sum = emis_cnt_it_co2sum.groupby('time.year').mean(dim="time")
    
    for pc in [0.025, 0.16, 0.84,0.975]:
        co2sum_cnt_pc_dum = emis_cnt_it_co2sum.quantile(pc,dim="nIt")
        co2sum_cnt_pc_an_dum  = emis_cnt_it_an_co2sum.quantile(pc,dim="nIt")
        
        co2sum_cnt_pc_list.append(co2sum_cnt_pc_dum)
        co2sum_cnt_pc_an_list.append(co2sum_cnt_pc_an_dum)
        
    emis_cnt_pc_co2sum = xarray.concat(co2sum_cnt_pc_list, dim="quantile")
    emis_cnt_pc_an_co2sum = xarray.concat(co2sum_cnt_pc_an_list, dim="quantile")
    

#%%

# 5. Save maps and countries to an output file
ds_out = xarray.Dataset()

for group in run_groups:
    
    ds_out["country_ap_"  + group] =  emis_cnt_ap[group]
    ds_out["country_ap_std_"  + group] =  uncert_cnt_ap[group]
    ds_out["country_mean_"  + group] =  emis_cnt_mean[group]
    ds_out["country_pc_"  + group] =  emis_cnt_pc[group].transpose("country", "time", "quantile")
    ds_out["emis_mean_"  + group] = emis_map_mean[group]
    ds_out["emis_ap_" + group] = da_emis_ap[group]
    ds_out["country_pc_annual_"  + group] =  emis_cnt_pc_an[group].transpose("country", "year", "quantile")

    ds_out["country_mean_"  + group].attrs["units"] = "Tg/yr"
    ds_out["country_pc_"  + group].attrs["units"] = "Tg/yr"
    ds_out["country_ap_"  + group].attrs["units"] = "Tg/yr"
    
    ds_out["emis_mean_"  + group].attrs["units"] = "kg/m2/s"
    ds_out["emis_ap_"  + group].attrs["units"] = "kg/m2/s"

    ds_out["bf_mean_" + group] = bf_mean[group].transpose("basis", "time")
    ds_out["bf_std_" + group] = bf_std[group].transpose("basis", "time")
    
    ds_out["bf_emis_mean_" + group] = bf_emis_mean[group].transpose("basis", "time")
    ds_out["bf_emis_mean_" + group].attrs["units"] = "Tg/yr"
    

if "CO2" in species:
    ds_out["country_mean_co2sum"] = emis_cnt_mean_co2sum
    ds_out["country_pc_co2sum"] = emis_cnt_pc_co2sum
    ds_out["country_pc_annual_co2sum"] = emis_cnt_pc_an_co2sum

if "co2gee" in run_groups:
    ds_out["country_mean_co2bio"] = emis_cnt_mean_co2bio
    ds_out["country_pc_co2bio"] = emis_cnt_pc_co2bio.transpose("country", "time", "quantile")
    ds_out["country_pc_annual_co2bio"] = emis_cnt_pc_an_co2bio.transpose("country", "year", "quantile")


#pd_rdates = pd.to_datetime(run_dates)
#ds_out.coords["time"] = pd_rdates
#ds_out.coords["country"] = countries
if "co2gee" in run_groups and "co2reco" in run_groups:
    ds_out["corr_gee_reco"] = (("country", "time"), corr_gee_reco)
    ds_out["corr_gee_ff"] = (("country", "time"), corr_gee_ff)
    ds_out["corr_reco_ff"] = (("country", "time"), corr_reco_ff)

if "CO2" in species:
    ds_out["corr_bio_ff"] = (("country", "time"), corr_bio_ff)
    ds_out["bf_corr_bio_ff"] = (("basis", "time"), np.transpose(bf_corr_bioff) )

if "CO" in species and "CO2" in species:
    ds_out["corr_coff_ff"] = (("country", "time"), corr_coff_ff)
    ds_out["bf_corr_coco2"] = (("basis", "time"), np.transpose(bf_corr_coco2) )

if "CO" in species:
    ds_out["coch4_mean"] =  coch4_mean
    ds_out["conmvoc_mean"] = covoc_mean
    
    ds_out["coch4_ap"] = coch4_ap
    ds_out["conmvoc_ap"] =  covoc_ap
    
    ds_out["coch4_pc"] = coch4_pc.transpose("time", "quantile")
    ds_out["conmvoc_pc"] = covoc_pc.transpose("time", "quantile")
    
    ds_out["coch4_pc_annual"] = coch4_pc_an.transpose("year", "quantile")
    ds_out["conmvoc_pc_annual"] = covoc_pc_an.transpose("year", "quantile")
    
    ds_out["country_mean_coch4"] = emis_cnt_mean_coch4
    ds_out["country_mean_conmvoc"] = emis_cnt_mean_covoc
    
    ds_out["country_ap_coch4"] = emis_cnt_ap_coch4
    ds_out["country_ap_conmvoc"] = emis_cnt_ap_covoc
    

fname_out = post_root + "emissions/" + "post_emissions_" + version +  ".nc"
for key in list(ds_out.keys()):
    ds_out[key].encoding['zlib'] = True 
ds_out.to_netcdf(path=fname_out)


