#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 16:35:37 2022

Calculate methane Dunbar landfill emissions from Ferry data

@author: mlunt
"""

import numpy as np
import xarray
import pandas as pd
import matplotlib.pyplot as plt
from acrg_grid import haversine
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io.img_tiles import GoogleTiles as OSM

def open_ds(fname, group=None):
    with xarray.open_dataset(fname, group) as ds:
        ds.load()
    return ds

def calc_angle(loc1, loc2, positive = False):
    """Calculate angle between two lat,lon points
    Args:
        loc1: Corrdinate pair for first point [lat,lon]
        loc2: Coordinate pair for second point [lat,lon]
    Returns:
        angle_rad: Angle in radians
        angle_deg: Angle in degrees
    """
    
    lat1_rad = np.radians(loc1[0])
    lon1_rad = np.radians(loc1[1])
    
    lat2_rad = np.radians(loc2[0])
    lon2_rad = np.radians(loc2[1])
    
    X = np.cos(lat2_rad) * np.sin(lon2_rad - lon1_rad)
    Y = np.cos(lat1_rad) * np.sin(lat2_rad)  - np.sin(lat1_rad) * np.cos(lat2_rad) * np.cos(lon2_rad - lon1_rad)
    
    angle_rad = np.arctan2(X,Y)
    angle_deg = np.degrees(angle_rad)
    
    if positive == True:
        
        wh_neg = np.where(angle_deg < 0)
        angle_deg[wh_neg] = angle_deg[wh_neg] + 360
        angle_rad[wh_neg] = angle_rad[wh_neg] + 2*np.pi
    
    return angle_rad, angle_deg

data_dir = "/home/mlunt/datastore/GAUGE/ferry/"
fname_in = data_dir + "ferry_dunbar_dwind_2015_2017.nc"


lat_dun = 55.98
lon_dun = -2.47

ds = open_ds(fname_in)

day_str_all = pd.to_datetime(ds.time.values).strftime('%Y%m%d')
date_arr = np.unique(day_str_all)
pd_date_arr =pd.to_datetime(date_arr)
#date_arr = np.unique(ds.time)

# Loop through each day
# Find where wind direction and Edinburgh angle intersect (as before)
# Need to define a background level from either side of where Edinburgh is (and ignore any other spikes)
# Use difference between them to define the enhancement.

#%%
diff=[]
diff_ch4=[]
diff_dates=[]
ed_indices=[]
pres_edi=[]
temp_edi=[]
ws_edi=[]

lon_min_list=[]
lon_max_list=[]
lat_min_list=[]
lat_max_list=[]

time_list=[]
cos_theta_list=[]
for date in date_arr[:-3]:
    ds_day = ds.sel(time=date)
    
    if len(ds_day.time) >0:
        
#        wh_good = np.where( (ds_day.methane > 1850) & 
#                   (ds_day.carbon_dioxide > 390) & 
#                   (ds_day.latitude > 40) & 
#                    (ds_day.longitude > -10) &
#                    (ds_day.wind_speed >= 0) & 
#                    (ds_day.wind_direction >= 0))[0]
        
        co2_day = ds_day.carbon_dioxide
        ch4_day = ds_day.methane
        lat_day = ds_day.latitude
        lon_day = ds_day.longitude
        ws_day = ds_day.wind_speed
        wd_day = ds_day.wind_direction
        course_day = ds_day.ship_course
        
        pres_day =  ds_day.air_pressure.where(ds_day.air_pressure >0)
        temp_day =  ds_day.air_temperature.where(ds_day.air_temperature >-100)

        rad_dun,deg_dun = calc_angle([lat_day,lon_day], [lat_dun, lon_dun], positive=True)
        
        # Might need to account for ship course when calculating wind direction. 
        
    
        wh_dun = np.where((np.abs(wd_day - deg_dun )< 20) & (ws_day > 2)  &
                             (deg_dun < 270) & (np.abs(wd_day -225) < 45) &
                             (lon_day >= -3.35) )[0]
        
       
        #wh_bg = np.where((np.abs(wd_day - deg_edi > 10)) & (ws_day > 2) )[0]
        # Maybe just define measurements to east as background
        
        #co2_edi = co2_day[wh_dun]
        #co2_bg = (co2_day[wh_dun[0]-5:wh_dun[0]].mean() + co2_day[wh_dun[-1]+1:wh_dun[-1]+5].mean())/2
        
        
        # Need to llok at mole fractions either side of wh_dun to get enhancement above non-Edi values 

                
        if len(wh_dun) >1:
            ch4_dun = ch4_day[wh_dun]
            ch4_bg = (ch4_day[wh_dun[0]-3:wh_dun[0]].mean() + ch4_day[wh_dun[-1]+1:wh_dun[-1]+4].mean())/2
        
        
            pres_edi.append(pres_day[wh_dun].values.mean() )
            temp_edi.append(temp_day[wh_dun].values.mean() )
            ws_edi.append(ws_day[wh_dun].values.mean() )
            
            diff.append(np.mean(ch4_dun.values) - np.mean(ch4_bg.values))
            #diff_ch4.append(np.mean(ch4_edi.values) - np.mean(ch4_bg.values))
            diff_dates.append(date)
            ed_indices.append(wh_dun)
            
            lon_min_list.append(lon_day.values[wh_dun[0]])
            lon_max_list.append(lon_day.values[wh_dun[-1]])
            
            lat_min_list.append(lat_day.values[wh_dun[0]])
            lat_max_list.append(lat_day.values[wh_dun[-1]])
            
            time_list.append(ds_day.time.values[wh_dun[0]])
            
            theta = np.mean(course_day[wh_dun] - (ws_day[wh_dun]-90) )
            cos_theta_list.append( np.cos(np.radians(theta)) ) 
            
            
diff_arr = np.asarray(diff)
diff_date_arr = np.asarray(diff_dates)

diff_plus = diff_arr[diff_arr >2]
date_plus = diff_date_arr[diff_arr >2]

ws_edi_arr = np.asarray(ws_edi)
pres_edi_arr = np.asarray(pres_edi)
temp_edi_arr = np.asarray(temp_edi)
time_edi_arr = np.asarray(time_list)

ws_plus = ws_edi_arr[diff_arr>2]
pres_plus = pres_edi_arr[diff_arr>2]
temp_plus = temp_edi_arr[diff_arr>2]

time_plus = time_edi_arr[diff_arr>2]
pd_time_plus = pd.to_datetime(time_plus)
hour_plus = pd_time_plus.hour

lon_min_arr = np.asarray(lon_min_list)
lon_max_arr = np.asarray(lon_max_list)

lat_min_arr = np.asarray(lat_min_list)
lat_max_arr = np.asarray(lat_max_list)

lon_min_plus = lon_min_arr[diff_arr>2]
lon_max_plus = lon_max_arr[diff_arr>2]

lat_min_plus = lat_min_arr[diff_arr>2]
lat_max_plus = lat_max_arr[diff_arr>2]

cos_theta_arr = np.asarray(cos_theta_list)
cos_theta_plus = cos_theta_arr[diff_arr>2]
#diff_arr_ch4 = np.asarray(diff_ch4)

#diff_plus_ch4 = diff_arr_ch4[diff_arr_ch4 >10]
#date_plus_ch4 = diff_date_arr[diff_arr_ch4 >10]

#%%
# Convert ppm enhamcenet to g/m3

# air density
# Calculate from Temp and pressure
R_spec = 287.05   # Specific gas constant for air

air_dense = pres_plus*1.e2/(R_spec*(temp_plus+273))  # kg/m3

air_dense = np.nan_to_num(air_dense)
air_dense[air_dense ==0] = 1.225

C_enh = diff_plus/1.e9 * 16/28.97 * air_dense*1000   # kg/m3
# Now calculate flux from this - how? Do I just tak ave enhancemnet and multiply by distance and height?
# Try it

# Calculate mean eind speed in perp direction
#ws_perp = ws_ex[ind0:ind1].mean()


# Calculate width of enhancement plume
nplus = len(diff_plus)
dist = np.zeros((nplus))
for ti in range(nplus):
    dist[ti] = haversine.distance([lat_min_plus[ti], lon_min_plus[ti]], 
        [lat_max_plus[ti], lon_max_plus[ti]])

delta_x = np.abs(cos_theta_plus)*dist

#%%
# Need to get boundary layer height

# Read in from MERRA2 or GEOS-FP files - only have 0.25x0.3125_EU for 2015 and 2017. 
# Or can maybe get from NAME runs? - Don't go far enough north
plon  = (lon_min_plus+lon_max_plus)/2
plat  = (lat_min_plus+lat_max_plus)/2
dz= 200  # Boundary layer height - make it up

met_root = "/geos/d21/GC_DATA/ExtData/GEOS_0.25x0.3125_EU/GEOS_FP/"
pblh_list=[]
pblh_list2=[]
for ti,date in enumerate(date_plus):
    
    met_dir = met_root + date[:4] + "/" + date[4:6] + "/"
    fname_met = met_dir + "GEOSFP." + date + ".A1.025x03125.EU.nc"
    
    ds_met = open_ds(fname_met)
    
    
    
    pblh_list.append(ds_met.PBLH.sel(lon = lon_dun, lat = lat_dun, time = time_plus[ti], method="nearest").values)
    pblh_list2.append(ds_met.PBLH.sel(lon = plon[ti], lat = plat[ti], time = time_plus[ti], method="nearest").values)
    
pblh_plus = np.asarray(pblh_list)
pblh_plus2 = np.asarray(pblh_list2)
#%%

F = C_enh * ws_plus * dist*1000 *dz   # Should be in g/s
F_annual = F*60*60*24*365/1.e9   # In Gg/yr = kt

F2 = C_enh * ws_plus * delta_x*1000 *pblh_plus  #pblh_plus   # Should be in g/s
F_annual2 = F2*60*60*24*365/1.e9  

# For Dunbar median emission are 2.6 kt /yr if  I assume PBLH of 200 m. 
# Total UK waste emissions should be around 35% of 2000 kt (Gg) so 700 kt

# If emissions are 300 tonnes what would enhancement be under certain conditions
F_theory = 0.3*1e9/60/60/24/365
delta_x_theory = 2.5*1000/2.5  # km of transect
ws_theory = 10
pblh_theory = np.arange(100,1100, 100)

enh_theory = F_theory /ws_theory/delta_x_theory / pblh_theory
diff_theory = enh_theory *28.97/16.04/1.225/1000.*1e9

#%%
cline = cfeature.NaturalEarthFeature('physical', 'coastline', '10m',
                                     facecolor="None")

# 21
for pdate in date_plus[29:30]:
    #pdate = "20160624"
    ds_day2 = ds.sel(time=pdate)
    
    plon = ds_day2.longitude
    plat = ds_day2.latitude
    pco2 = ds_day2.carbon_dioxide
    pch4 = ds_day2.methane
    pwd = ds_day2.wind_direction
    pws = ds_day2.wind_speed

    u = pws * np.cos(np.radians(270-pwd)) # Have to switch to math coords relative to x axis
    v = pws * np.sin(np.radians(270-pwd))
    
    cmin = np.percentile(pch4,5)
    cmax = np.percentile(pch4,95)
    
    cline = cfeature.NaturalEarthFeature('physical', 'coastline', '10m',
                                         facecolor="None")
    
    proj = ccrs.PlateCarree()
    fig,ax = plt.subplots(subplot_kw=dict(projection=proj),figsize=(8,8))
    
    #ax.set_extent((-6,6,50,57))
    ax.set_extent((-3.5,-0.5,54.7,56.2))
    ax.scatter(plon,plat, c=pch4, cmap="RdYlBu_r", vmin = cmin, vmax=cmax)
    
    ax.barbs(plon,plat, u,v, length=5)
    #ax.coastlines()
    ax.add_feature(cline, edgecolor='black')
    
    ax.scatter(lon_dun,lat_dun, s=8)
    
    
#%%
# Plot on google maps background
# 20170818 is a good example # index 2019
# 
    
tiler = OSM(style='satellite')
mercator = tiler.crs

idx = 29

ds_day2 = ds.sel(time=date_plus[idx])
plon = ds_day2.longitude
plat = ds_day2.latitude
pco2 = ds_day2.carbon_dioxide
pch4 = ds_day2.methane
pwd = ds_day2.wind_direction
pws = ds_day2.wind_speed

u = pws * np.cos(np.radians(270-pwd)) # Have to switch to math coords relative to x axis
v = pws * np.sin(np.radians(270-pwd))

fig,ax = plt.subplots(subplot_kw=dict(projection=proj),figsize=(8,8))
ax.set_extent((-2.6, -2.2, 55.85, 56.1))
zoom = 12  #6 or 9 depending on zoom
ax.add_image(tiler, zoom )

h3 = ax.scatter(plon,plat, c=pch4, cmap="RdYlBu_r", vmin = 1960, vmax=2000)
ax.barbs(plon,plat, u,v, length=8)
ax.scatter(lon_dun,lat_dun, s=20, color='red')

cbaxes = fig.add_axes([0.20, 0.2, 0.6, 0.02]) 
##[left, bottom, width, height],
cb = plt.colorbar(h3, cax = cbaxes, orientation='horizontal', extend='both', label = 'CH4 mole fraction (ppb)') 
ax.set_title(date_plus[idx])