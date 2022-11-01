#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 11:31:42 2020

Extract u, v from met files.

Want the area average wind direction (u and v) for 12 each day. 
Take this from bottom 5 levels or so maybe?  

Having looked at winds maybe go bottom 6 layers.Goes to about 725 m up. 
Assume a mid-day PBL of at leat 1000 m and this is probably fine.

Pick time index =4

Save as an annual file - so need to work out hpw to loop through months and days

@author: mlunt
"""

import numpy as np
import xarray
import matplotlib.pyplot as plt
import pandas as pd
import json
from datetime import timedelta

def open_ds(fname):
    with xarray.open_dataset(fname) as ds:
        ds.load()
    return ds

data_root = "/geos/d21/GC_DATA/ExtData/GEOS_0.25x0.3125_EU/GEOS_FP/"

year = "2017"

start_date = year + "0101"
end_date = year + "1231"

#dates= [ "20190321"]

pd_dates = pd.date_range(start_date,end_date, freq='D')

json_dir = "/home/mlunt/programs/Python/DARE/inv_code/"
# Read in json file
with open(json_dir + "dare_site_info.json") as f:
    site_info=json.load(f)
    
# Use sites which have both CO and CO2. Just focus on TAC for now.
#sites = ["TAC"]

sites = ["GAT", "HTM",  
        "KRE",  "NOR", "OPE",
        "SAC" ,  "SVB", "TRN", 
          "BSD", "HFD"]   # TAC - no need already done

lat_site={}
lon_site={}
levi = {}
for site in sites:
    dict_si = site_info[site]
    network_si = dict_si["network"]
    
    lat_site[site] = dict_si["latitude"]
    lon_site[site] = dict_si["longitude"]
    alt_si = dict_si["height_magl"][0]
   
    if alt_si < 50: # 1st level goes to ~ 120 m 
        levi[site] =0
    elif  50 <= alt_si < 200: # 2nd leve lgos to ~ 250 m
        levi[site]=1
    elif 200 <= alt_si <400: # lEvel 3 top around 400 m 
        levi[site]=2
    else:
        levi[site]=3

t_offset = timedelta(hours=1, minutes=30)

U_list={}
V_list={}
pblh_list={}

for site in sites:
    U_list[site]=[]
    V_list[site]=[]
    pblh_list[site]=[]

#psurf_list=[]
for pd_date in pd_dates:

    date_str = pd_date.strftime("%Y%m%d")
    month = date_str[4:6]
    
    data_dir = data_root + year + "/" + month + "/"
    
    #/geos/d21/GC_DATA/ExtData/GEOS_0.25x0.3125_EU/GEOS_FP/2017/01/
    #GEOSFP.20170131.A3dyn.025x03125.EU.nc
    
    file_str1 = "GEOSFP." + date_str + ".A3dyn." + "025x03125.EU.nc"
    file_str2 = "GEOSFP." + date_str + ".A1." + "025x03125.EU.nc"
    #file_str3 = "GEOSFP." + date_str + ".I3." + "025x03125.EU.nc"
    
    
    fname1 = data_dir + file_str1
    fname2 = data_dir + file_str2
    #fname3 = data_dir + file_str3
    ds_uv = open_ds(fname1)
    ds_pbl= open_ds(fname2)
    #ds_ps = open_ds(fname3)
    
    for site in sites:
        U = ds_uv.U.sel(lon=lon_site[site], lat=lat_site[site], method="nearest")
        V = ds_uv.V.sel(lon=lon_site[site], lat=lat_site[site], method="nearest")
    
        pblh = ds_pbl.PBLH.sel(lon=lon_site[site], lat=lat_site[site], method="nearest")
    #psurf = ds_ps.PS.sel(lon=lon_site, lat=lat_site, method="nearest")
    
        U_site = U[:,levi[site]]
        V_site = V[:,levi[site]]

        pblh_3hr = pblh.resample(time='3H', loffset = t_offset).mean()
    
        U_list[site].append(U_site)
        V_list[site].append(V_site)
        
        pblh_list[site].append(pblh_3hr)
        #psurf_list.append(psurf)
    
   
    
#%%
for site in sites:
    U_all = xarray.concat(U_list[site], dim="time")  
    V_all = xarray.concat(V_list[site], dim="time")  
    pblh_all = xarray.concat(pblh_list[site], dim="time")  
    #psurf_all = xarray.concat(psurf_list, dim="time")    
        
    ds_out  = xarray.Dataset()
    
    ds_out["U"] = U_all
    ds_out["V"] = V_all
    
    ds_out["pblh"] = pblh_all
    #ds_out["psurf"] = psurf_all
        
    # Convert to pandas dataframe
    df_out = ds_out.to_dataframe()
        
    fname_out = "/home/mlunt/datastore/DARE/ml_data/met/GEOSFP_025x03125_" + site + "_" + year + ".csv"
        
    df_out.to_csv(fname_out)

