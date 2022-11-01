#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 12:32:13 2022

Read in met data and model data for a given site

Concat together on same time indices - make sure everything alligned.
Put into one dataframe to use for ML

Feature engineering - U + v into ws and direction;
                    - also might be worth using CO:CO2 ratio
                    - convert time index into hour, day, month, year
                    
Clean up data - remove any missing indices. 

Model - try random forest and xgboost

@author: mlunt
"""

import glob
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.pipeline import Pipeline
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error
from sklearn.feature_selection import mutual_info_regression
from sklearn.model_selection import GridSearchCV
from xgboost import XGBRegressor

def make_mi_scores(X, y):
    mi_scores = mutual_info_regression(X, y)
    mi_scores = pd.Series(mi_scores, name="MI Scores", index=X.columns)
    mi_scores = mi_scores.sort_values(ascending=False)
    return mi_scores

#%%
met_dir = "/home/mlunt/datastore/DARE/ml_data/met/"
model_dir = "/home/mlunt/datastore/DARE/ml_data/model_data/"
out_dir = "/home/mlunt/datastore/DARE/ml_data/xgb_output/"
sites = ["GAT", "HTM",  
        "KRE",  "NOR", "OPE",
        "SAC" ,  "SVB", "TRN",
          "BSD", "HFD", "TAC"]  
#sites = ["HFD"] #, "BSD"]
#site = "TAC"
for site in sites:

    # Read in met data
    fnames_met = glob.glob(met_dir + "GEOSFP_025x03125_" + site + "*.csv")
    
    fnames_met.sort()
    df_met_list=[]
    for fname in fnames_met:
        
        df_dum = pd.read_csv(fname, index_col = "time", parse_dates=True)
        df_met_list.append(df_dum)
    
    df_met = pd.concat(df_met_list)
    
    # Read in model output data
    
    fname_model =  model_dir + site + "_GC_output.csv"
    
    df_mod = pd.read_csv(fname_model, index_col="time", parse_dates=True)
    
    #%%
    
    df_join = pd.merge(df_mod,df_met, how="inner", left_index=True, right_index=True)
    df_join2 = df_join.dropna(axis=0, how="any")
    
    # Calculate ws and wd
    ws = np.sqrt((df_join2["U"]**2 + df_join2["V"]**2))
    wd = (np.arctan2(df_join2["V"], df_join2["U"]))  # np.rad2deg(np.arctan2(df_join2["V"], df_join2["U"]))
    
    
    #df_add = pd.DataFrame( "wd": wd.rename("wd"))
    
    df_join2["wd"]  = wd.rename("wd")
    df_join2["ws"] = ws.rename("ws")
    # Extract year, day of year and hour value from time
    df_join2["year"] = df_join2.index.year
    df_join2["doy"] = df_join2.index.dayofyear
    df_join2["hour"] = df_join2.index.hour
    
    #df_train  = df_join2["20170101":"20191130"]
    #df_test = df_join2["20191201":]
    
    x_cols = ["co2_mod", "co_mod", "U", "V", "pblh", "ws", "wd", "year", "doy" ,"hour"]
    #y = df_train["co2ff_mod"]
    #X = df_train[x_cols]
    
    #y_test = df_test["co2ff_mod"]
    #X_test = df_test[x_cols]
    
    train_X, val_X, train_y, val_y = train_test_split(df_join2[x_cols], df_join2["co2ff_mod"], random_state = 0)
    
    #%%
    model = RandomForestRegressor(n_estimators=100, random_state=0)
    my_pipeline = Pipeline(steps=[('model', model) ])
    
    #scores = -1 * cross_val_score(my_pipeline, train_X, train_y,
    #                              cv=5,
    #                              scoring='neg_mean_absolute_error')
    
    #print("MAE scores:\n", scores)
    
    my_pipeline.fit(train_X,train_y)
    y_predict = my_pipeline.predict(val_X)
    
    mi_scores = make_mi_scores(train_X,train_y)
    
    mae = mean_absolute_error(val_y, y_predict)
    print(site + " MAE = ", mae)
    
    
    #%%
    # Use gridsearchCV to optimize random forest parameters
    param_grid = { 
    'n_estimators': [50,100,200, 400],
    'max_features': ['auto', 'sqrt', 'log2'],
    'max_depth' : [None, 2, 4,8]}
    
#    rfc_model =  RandomForestRegressor(n_estimators=100, random_state=0)
#    CV_rfc = GridSearchCV(estimator=rfc_model, param_grid=param_grid, cv= 5)
#    CV_rfc.fit(train_X, train_y)
#    print(CV_rfc.best_params_)
    
    #best params = {'max_depth': 8, 'max_features': 'auto', 'n_estimators': 400}

    model_opt = RandomForestRegressor(n_estimators=400, random_state=0, max_features='auto',
                                      max_depth=None)
    
    model_opt.fit(train_X,train_y)
    y_predict_opt = model_opt.predict(val_X)
    
    
    mae2 = mean_absolute_error(val_y, y_predict_opt)
    print(site + " optimized RF MAE = ", mae2)
    
    
    #%%
    # Test xgboost performance
    my_model = XGBRegressor(n_estimators=500, learning_rate=0.05, n_jobs=32, random_state=0)
    my_model.fit(train_X, train_y, 
             early_stopping_rounds=5, 
             eval_set=[(val_X, val_y)],
             verbose=False)
    
    y_pred3 = my_model.predict(val_X)
    mae3 = mean_absolute_error(val_y, y_pred3)
    print(site + " XGBoost MAE = ", mae3)
    
    #%%
    
    # Replace Co2_mod and co_mod with observed values instead. How well does Co2ff compare to model fitted on model data. 
    #val_x2 = 
    x_cols2 = ["co2_obs", "co_obs", "U", "V", "pblh", "ws", "wd", "year", "doy" ,"hour"]
    X2 = df_join2[x_cols2]
    X2 = X2.rename(columns = {"co2_obs": "co2_mod",
               "co_obs": "co_mod"})
        
#    x_cols2 = ["co2_mod_perturb", "co_mod_perturb", "U", "V", "pblh", "ws", "wd", "year", "doy" ,"hour"]
#    X2 = df_join2[x_cols2]
#    X2 = X2.rename(columns = {"co2_mod_perturb": "co2_mod",
#               "co_mod_perturb": "co_mod"})
#    X2["co2_smod"] = X2["co2_mod"] - df_join2["co2ff_mod"].values*0.1
    #X2["co_mod"] = X2["co_mod"] - df_join2["co_mod"].values*0.25
        
    co2ff_pred = my_pipeline.predict(X2)
    co2ff_pred3 = my_model.predict(X2)
    co2ff_true = df_join2["co2ff_mod"]
    
    
    # Performance worsens when I have new "true" value and different co2bio and cochem contributions. 
    # But maybe just see how inversion performs using "real" data for now. 
    
    df_out = pd.DataFrame(data = {"co2ff_obs": pd.Series(co2ff_pred3, index = co2ff_true.index)})
    
    fname_out = out_dir + site + "_co2ff_xgb.csv"
    
    df_out.to_csv(fname_out)
    