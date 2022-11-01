#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 13:50:53 2022

Detect point source/small area plumes in ferry data 

@author: mlunt
"""

import numpy as np
import xarray
import areagrid
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colorbar as cbar
import glob
from acrg_grid import haversine

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

fname = data_dir + "picarro-vaisala_finlandia-seaways_20170101.nc"

ds_temp = open_ds(fname)

ds = ds_temp.where((ds_temp.latitude > 52.4) & (ds_temp.longitude > -3.) &
                   (ds_temp.methane > 1850) &  (ds_temp.carbon_dioxide > 380) & 
                   (ds_temp.carbon_dioxide < 450) & 
                   (ds_temp.wind_speed >= 4) & (ds_temp.wind_direction >= 0),
                   drop=True)

#ds_west = ds.where((ds.wind_direction >= 180) & (ds.wind_direction < 360) &
#                   (ds_temp.wind_speed >= 2), drop=True)
#%%
# Also filter out times when ship course is opposite to wind direction (or the same)

angle_diff  = np.abs(ds.ship_course - ds.wind_direction)


ds_good = ds.where(np.abs(angle_diff - 180) > 30, drop=True)

#%%
# Groupby day to get enhancement
co2 = ds_good.carbon_dioxide
co2_day = co2.groupby('time.dayofyear').min()
co2_plus = co2.groupby('time.dayofyear') - co2_day

ch4 = ds_good.methane
ch4_day = ch4.groupby('time.dayofyear').min()
ch4_plus = ch4.groupby('time.dayofyear') - ch4_day

lon = ds_good.longitude
lat= ds_good.latitude

ws = ds_good.wind_speed
wd = ds_good.wind_direction

sort_idx= np.argsort(co2_plus).values
co2_plus_sort = co2_plus[sort_idx]
lat_sort = lat[sort_idx]
lon_sort = lon[sort_idx]

sort_idx2= np.argsort(ch4_plus).values
ch4_plus_sort = ch4_plus[sort_idx2]
lat_sort2 = lat[sort_idx2]
lon_sort2 = lon[sort_idx2]

#%%
# Use a z-score based detection algorithm
# Limit data to between Dunbar and 

lag = 5
threshold =2.5

nobs = len(ch4)

window_mean=np.zeros((nobs))
window_std=np.zeros((nobs))
z=np.zeros((nobs))

window_mean_co2=np.zeros((nobs))
window_std_co2=np.zeros((nobs))
z_co2=np.zeros((nobs))
for ti in range(lag,nobs-lag):
    
    window_mean[ti] = ch4[ti-lag:ti+lag].mean()
    window_std[ti] = ch4[ti-lag:ti+lag].std()
    
    z[ti] = (ch4[ti]-window_mean[ti])/window_std[ti]
    
    window_mean_co2[ti] = co2[ti-lag:ti+lag].mean()
    window_std_co2[ti] = co2[ti-lag:ti+lag].std()
    
    z_co2[ti] = (co2[ti]-window_mean_co2[ti])/window_std_co2[ti]


z_sort = z[sort_idx2]
z_binary = z.copy()*0.
z_binary[z >threshold] = 1


z_binary_co2 = z_co2.copy()*0.
z_binary_co2[z_co2 >threshold] = 1

lat_plus_co2 = lat[z_binary_co2==1]
lon_plus_co2 = lon[z_binary_co2==1]

ws_plus_co2 = ws[z_binary_co2==1]
wd_plus_co2 = wd[z_binary_co2==1]

lat_plus = lat[z_binary==1]
lon_plus = lon[z_binary==1]

ws_plus = ws[z_binary==1]
wd_plus = wd[z_binary==1]

u = ws_plus * np.cos(np.radians(270-wd_plus)) # Have to switch to math coords relative to x axis
v = ws_plus * np.sin(np.radians(270-wd_plus))

u_co2 = ws_plus_co2 * np.cos(np.radians(270-wd_plus_co2)) # Have to switch to math coords relative to x axis
v_co2 = ws_plus_co2 * np.sin(np.radians(270-wd_plus_co2))
#%%
cline = cfeature.NaturalEarthFeature('physical', 'coastline', '10m',
                                     facecolor="None")

cmin = 0
cmax = 40

proj = ccrs.PlateCarree()
fig,ax = plt.subplots(subplot_kw=dict(projection=proj),figsize=(8,8))

#ax.set_extent((-6,6,50,57))
#ax.set_extent((-3.5,-0.5,54.7,56.2))
#ax.scatter(lon_sort,lat_sort, c=co2_plus_sort, cmap="RdYlBu_r", vmin = cmin, vmax=cmax)
ax.scatter(lon_plus_co2,lat_plus_co2)


ax.barbs(lon_plus_co2,lat_plus_co2, u_co2,v_co2, length=5)
#ax.coastlines()
ax.add_feature(cline, edgecolor='black')

#%%
cmin2=0
cmax2 = 2
proj = ccrs.PlateCarree()
fig,ax = plt.subplots(subplot_kw=dict(projection=proj),figsize=(8,8))

#ax.set_extent((-6,6,50,57))
#ax.set_extent((-3.5,-0.5,54.7,56.2))
#ax.scatter(lon_sort2,lat_sort2, c=z_sort, cmap="RdYlBu_r", vmin = cmin2, vmax=cmax2)

ax.scatter(lon_plus,lat_plus)

ax.barbs(lon_plus,lat_plus, u,v, length=5)

ax.add_feature(cline, edgecolor='black')
