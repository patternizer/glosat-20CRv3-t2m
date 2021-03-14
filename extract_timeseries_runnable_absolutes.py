#!/usr/bin/env python

#-----------------------------------------------------------------------
# PROGRAM: extract_timeseries_runnable_absolutes.py
#-----------------------------------------------------------------------
# Version 0.1
# 12 March, 2021
# Dr Michael Taylor
# https://patternizer.github.io
# patternizer AT gmail DOT com
#-----------------------------------------------------------------------

# Data libraries
import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime

# Plotting libraries:
import matplotlib
matplotlib.use('agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker

# Maths libraries:
from scipy.interpolate import griddata
from scipy import spatial
from math import radians, cos, sin, asin, sqrt

# OS libraries:
import os, sys
from  optparse import OptionParser
#import argparse

# Silence library version notifications
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

def locate_timeseries(target_lat,target_lon):

    #-----------------------------------------------------------------------------
    # SETTINGS
    #-----------------------------------------------------------------------------

    fontsize = 16
    use_anomalies = False
    use_yearly = False
    make_plot = True

    if target_lon < 0:
        target_lon += 360.0
    targetstr = 'lon_' + str(round(target_lon,2)) + '_' + 'lat_' + str(round(target_lat,2))

    if use_anomalies == True:
        if use_yearly == True:
            dirstr = 'ANOMALIES/yearly/'
            typestr = 'anomaly_yearly'
        else:
            dirstr = 'ANOMALIES'        
            typestr = 'anomaly'    
    else:
        if use_yearly == True:
            dirstr = 'ABSOLUTES/yearly/'
            typestr = 'absolute_yearly'
        else:
            dirstr = 'ABSOLUTES'    
            typestr = 'absolute'

    outdir = 'EXTRACT_RUN'
    ensemble_mean = dirstr+'/'+'ensemble_mean'+'_'+typestr+'.nc'
    ensemble_sd = dirstr+'/ensemble_sd'+'_'+typestr+'.nc'
    ensemble_min = dirstr+'/ensemble_min'+'_'+typestr+'.nc'
    ensemble_max = dirstr+'/ensemble_max'+'_'+typestr+'.nc'
    ensemble_pctl_05 = dirstr+'/ensemble_pctl_05'+'_'+typestr+'.nc'
    ensemble_pctl_10 = dirstr+'/ensemble_pctl_10'+'_'+typestr+'.nc'
    ensemble_pctl_25 = dirstr+'/ensemble_pctl_25'+'_'+typestr+'.nc'
    ensemble_pctl_50 = dirstr+'/ensemble_pctl_50'+'_'+typestr+'.nc'
    ensemble_pctl_75 = dirstr+'/ensemble_pctl_75'+'_'+typestr+'.nc'
    ensemble_pctl_90 = dirstr+'/ensemble_pctl_90'+'_'+typestr+'.nc'
    ensemble_pctl_95 = dirstr+'/ensemble_pctl_95'+'_'+typestr+'.nc'

    #-----------------------------------------------------------------------------
    # METHODS
    #-----------------------------------------------------------------------------

    def haversine(lat1, lon1, lat2, lon2):

        # convert decimal degrees to radians 
        lat1, lon1, lat2, lon2 = map(radians, [lat1, lon1, lat2, lon2])

        # haversine formula 
        dlon = lon2 - lon1 
        dlat = lat2 - lat1 
        a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
        c = 2 * asin(sqrt(a)) 
        # Radius of earth is 6371 km
        km = 6371* c 
        return km

    def find_nearest(lat, lon, df):

        # Find nearest cell in array
        distances = df.apply(lambda row: haversine(lat, lon, row['lat'], row['lon']), axis=1)
        return df.loc[distances.idxmin(),:]

    #-----------------------------------------------------------------------------
    # LOAD: ENSEMBLE MEAN TIMESERIES COORDS + UNPACK
    #-----------------------------------------------------------------------------

    ds_mean = xr.open_dataset(ensemble_mean, decode_cf=True)
    lat = ds_mean.lat
    lon = ds_mean.lon
    X,Y = np.meshgrid(lon,lat)
    N = len(lat)*len(lon)
    x = X.reshape(N)
    y = Y.reshape(N)
    dlatlon = pd.DataFrame({'lon':x, 'lat':y}, index=range(N))

    #-----------------------------------------------------------------------------
    # FIND CLOSEST GRIDCELL (using Haversine distance) + CALC DISTANCE TO CENTRE
    #-----------------------------------------------------------------------------

    pt = [target_lat,target_lon]  
    query = find_nearest(pt[0],pt[1],dlatlon)
    nearest_lat = query.lat
    nearest_lon = query.lon
    nearest_cell = query.name
    nearest_lat_idx = np.where(ds_mean.lat==nearest_lat)[0][0] # --> 20CRv3 grid lat idx
    nearest_lon_idx = np.where(ds_mean.lon==nearest_lon)[0][0] # --> 20CRv3 grid lon idx
    delta_lon = (x[1]-x[0])/2
    delta_lat = (y[1]-y[0])/2
    distance_to_gridcell_centre = haversine(target_lat, target_lon, nearest_lat+delta_lat, nearest_lon+delta_lon)

    #-----------------------------------------------------------------------------
    # EXTRACT GRIDCELL TIMESERIES + CONVERT TO DEGREES CENTIGRADE (not SD) + SAVE
    #-----------------------------------------------------------------------------

    if use_anomalies == True:
        ts_mean = ds_mean.TMP2m[:,0,nearest_lat_idx,nearest_lon_idx]
    else:
        ts_mean = ds_mean.TMP2m[:,0,nearest_lat_idx,nearest_lon_idx] - 273.15
    t = ds_mean.TMP2m[:,0,nearest_lat_idx,nearest_lon_idx].time
    t_mean = pd.date_range(start=str(t[0].values)[0:4], periods=len(t), freq='M')   
    df_mean = pd.DataFrame({'t':t_mean,'ts':ts_mean})
    df_mean.to_csv(outdir+'/'+'ensemble_mean_' + targetstr + '.csv')

    ds_sd = xr.open_dataset(ensemble_sd, decode_cf=True); 
    ts_sd = ds_sd.TMP2m[:,0,nearest_lat_idx,nearest_lon_idx]
    t = ds_sd.TMP2m[:,0,nearest_lat_idx,nearest_lon_idx].time
    t_sd = pd.date_range(start=str(t[0].values)[0:4], periods=len(t), freq='M')   
    df_sd = pd.DataFrame({'t':t_sd,'ts':ts_sd})
    df_sd.to_csv(outdir+'/'+'ensemble_sd_' + targetstr + '.csv')

    if os.path.exists(ensemble_pctl_05): 
        use_pctl_05 = True
        ds_pctl_05 = xr.open_dataset(ensemble_pctl_05, decode_cf=True); 
        if use_anomalies == True:
            ts_pctl_05 = ds_pctl_05.TMP2m[:,0,nearest_lat_idx,nearest_lon_idx]
        else:
            ts_pctl_05 = ds_pctl_05.TMP2m[:,0,nearest_lat_idx,nearest_lon_idx] - 273.15
        t = ds_pctl_05.TMP2m[:,0,nearest_lat_idx,nearest_lon_idx].time
        t_pctl_05 = pd.date_range(start=str(t[0].values)[0:4], periods=len(t), freq='M')   
        df_pctl_05 = pd.DataFrame({'t':t_pctl_05,'ts':ts_pctl_05})
        df_pctl_05.to_csv(outdir+'/'+'ensemble_pctl_05_' + targetstr + '.csv')
    else: 
        use_pctl_05 = False
    if os.path.exists(ensemble_pctl_10): 
        use_pctl_10 = True
        ds_pctl_10 = xr.open_dataset(ensemble_pctl_10, decode_cf=True); 
        if use_anomalies == True:
            ts_pctl_10 = ds_pctl_10.TMP2m[:,0,nearest_lat_idx,nearest_lon_idx]
        else:
            ts_pctl_10 = ds_pctl_10.TMP2m[:,0,nearest_lat_idx,nearest_lon_idx] - 273.15
        t = ds_pctl_10.TMP2m[:,0,nearest_lat_idx,nearest_lon_idx].time
        t_pctl_10 = pd.date_range(start=str(t[0].values)[0:4], periods=len(t), freq='M')   
        df_pctl_10 = pd.DataFrame({'t':t_pctl_10,'ts':ts_pctl_10})
        df_pctl_10.to_csv(outdir+'/'+'ensemble_pctl_10_' + targetstr + '.csv')
    else: 
        use_pctl_10 = False
    if os.path.exists(ensemble_pctl_25): 
        use_pctl_25 = True
        ds_pctl_25 = xr.open_dataset(ensemble_pctl_25, decode_cf=True); 
        if use_anomalies == True:
            ts_pctl_25 = ds_pctl_25.TMP2m[:,0,nearest_lat_idx,nearest_lon_idx]
        else:
            ts_pctl_25 = ds_pctl_25.TMP2m[:,0,nearest_lat_idx,nearest_lon_idx] - 273.15
        t = ds_pctl_25.TMP2m[:,0,nearest_lat_idx,nearest_lon_idx].time
        t_pctl_25 = pd.date_range(start=str(t[0].values)[0:4], periods=len(t), freq='M')   
        df_pctl_25 = pd.DataFrame({'t':t_pctl_25,'ts':ts_pctl_25})
        df_pctl_25.to_csv(outdir+'/'+'ensemble_pctl_25_' + targetstr + '.csv')
    else: 
        use_pctl_25 = False
    if os.path.exists(ensemble_pctl_50): 
        use_pctl_50 = True
        ds_pctl_50 = xr.open_dataset(ensemble_pctl_50, decode_cf=True); 
        if use_anomalies == True:
            ts_pctl_50 = ds_pctl_50.TMP2m[:,0,nearest_lat_idx,nearest_lon_idx]
        else:
            ts_pctl_50 = ds_pctl_50.TMP2m[:,0,nearest_lat_idx,nearest_lon_idx] - 273.15
        t = ds_pctl_50.TMP2m[:,0,nearest_lat_idx,nearest_lon_idx].time
        t_pctl_50 = pd.date_range(start=str(t[0].values)[0:4], periods=len(t), freq='M')   
        df_pctl_50 = pd.DataFrame({'t':t_pctl_50,'ts':ts_pctl_50})
        df_pctl_50.to_csv(outdir+'/'+'ensemble_pctl_50_' + targetstr + '.csv')
    else: 
        use_pctl_50 = False
    if os.path.exists(ensemble_pctl_75): 
        use_pctl_75 = True
        ds_pctl_75 = xr.open_dataset(ensemble_pctl_75, decode_cf=True); 
        if use_anomalies == True:
            ts_pctl_75 = ds_pctl_75.TMP2m[:,0,nearest_lat_idx,nearest_lon_idx]
        else:
            ts_pctl_75 = ds_pctl_75.TMP2m[:,0,nearest_lat_idx,nearest_lon_idx] - 273.15
        t = ds_pctl_75.TMP2m[:,0,nearest_lat_idx,nearest_lon_idx].time
        t_pctl_75 = pd.date_range(start=str(t[0].values)[0:4], periods=len(t), freq='M')   
        df_pctl_75 = pd.DataFrame({'t':t_pctl_75,'ts':ts_pctl_75})
        df_pctl_75.to_csv(outdir+'/'+'ensemble_pctl_75_' + targetstr + '.csv')
    else: 
        use_pctl_75 = False
    if os.path.exists(ensemble_pctl_90): 
        use_pctl_90 = True
        ds_pctl_90 = xr.open_dataset(ensemble_pctl_90, decode_cf=True); 
        if use_anomalies == True:
            ts_pctl_90 = ds_pctl_90.TMP2m[:,0,nearest_lat_idx,nearest_lon_idx]
        else:
            ts_pctl_90 = ds_pctl_90.TMP2m[:,0,nearest_lat_idx,nearest_lon_idx] - 273.15
        t = ds_pctl_90.TMP2m[:,0,nearest_lat_idx,nearest_lon_idx].time
        t_pctl_90 = pd.date_range(start=str(t[0].values)[0:4], periods=len(t), freq='M')   
        df_pctl_90 = pd.DataFrame({'t':t_pctl_90,'ts':ts_pctl_90})
        df_pctl_90.to_csv(outdir+'/'+'ensemble_pctl_90_' + targetstr + '.csv')
    else: 
        use_pctl_90 = False
    if os.path.exists(ensemble_pctl_95): 
        use_pctl_95 = True
        ds_pctl_95 = xr.open_dataset(ensemble_pctl_95, decode_cf=True); 
        if use_anomalies == True:
            ts_pctl_95 = ds_pctl_95.TMP2m[:,0,nearest_lat_idx,nearest_lon_idx]
        else:
            ts_pctl_95 = ds_pctl_95.TMP2m[:,0,nearest_lat_idx,nearest_lon_idx] - 273.15
        t = ds_pctl_95.TMP2m[:,0,nearest_lat_idx,nearest_lon_idx].time
        t_pctl_95 = pd.date_range(start=str(t[0].values)[0:4], periods=len(t), freq='M')   
        df_pctl_95 = pd.DataFrame({'t':t_pctl_95,'ts':ts_pctl_95})
        df_pctl_95.to_csv(outdir+'/'+'ensemble_pctl_95_' + targetstr + '.csv')
    else: 
        use_pctl_95 = False

    # TO DO: CONVERT TO CRUTEM FORMAT + SAVE

    if make_plot == True:

        #-----------------------------------------------------------------------------
        # PLOT: timeseries ensemble median + credible intervals
        #-----------------------------------------------------------------------------

        figstr = outdir+'/' + targetstr + '.png'
        titlestr = '('+str(round(target_lon,3))+','+str(round(target_lat,3))+'): 20CRv3 timeseries extract at gridcell ('+str(round(nearest_lon,3))+','+str(round(nearest_lat,3))+') separation='+str(round(distance_to_gridcell_centre,3))+'km'

        fig,ax = plt.subplots(figsize=(15,10))
        if use_pctl_05 & use_pctl_95 == True:
            plt.fill_between(df_pctl_05.t, df_pctl_05.ts.rolling(360).mean(), df_pctl_95.ts.rolling(360).mean(), alpha=0.25, color='red', lw=1, label='5%-95% credible interval') # smooth on 30-year climate scale
        if use_pctl_10 & use_pctl_90 == True:
            plt.fill_between(df_pctl_10.t, df_pctl_10.ts.rolling(360).mean(), df_pctl_90.ts.rolling(360).mean(), alpha=0.5, color='red', lw=1, label='10%-90% credible interval') # smooth on 30-year climate scale
        if use_pctl_25 & use_pctl_75 == True:
            plt.fill_between(df_pctl_25.t, df_pctl_25.ts.rolling(360).mean(), df_pctl_75.ts.rolling(360).mean(), alpha=0.75, color='red', lw=1, label='25%-75% credible interval') # smooth on 30-year climate scale
        if use_pctl_50 == True:
            plt.step(df_pctl_50.t,df_pctl_50.ts.rolling(360).mean(), color='black', lw=3, label='median') # smooth on 30-year climate scale
        plt.step(df_mean.t,df_mean.ts.rolling(360).mean(), color='blue', lw=1, linestyle='dashed', label='mean') # smooth on 30-year climate scale            
        plt.fill_between(df_mean.t, df_mean.ts.rolling(360).mean()+df_sd.ts.rolling(360).mean(), df_mean.ts.rolling(360).mean()-df_sd.ts.rolling(360).mean(), alpha=0.1, color='blue', lw=1, label='$\pm$ 1SD') # smooth on 30-year climate scale
        plt.step(df_mean.t,df_mean.ts.rolling(360).mean()+df_sd.ts.rolling(360).mean(), color='blue', lw=1, linestyle='dashdot') # smooth on 30-year climate scale
        plt.step(df_mean.t,df_mean.ts.rolling(360).mean()-df_sd.ts.rolling(360).mean(), color='blue', lw=1, linestyle='dashdot') # smooth on 30-year climate scale
        plt.tick_params(labelsize=fontsize)
        ax.xaxis.grid(True, which='major')        
        ax.yaxis.grid(True, which='major')        
        plt.legend(loc=4, ncol=1, fontsize=fontsize)
        plt.xlabel("Year", fontsize=fontsize)
        if use_anomalies == True:
            plt.ylabel("Temperature anomaly, $\mathrm{\degree}C$", fontsize=fontsize)
        else:
            plt.ylabel("Absolute temperature, $\mathrm{\degree}C$", fontsize=fontsize)
        plt.title(titlestr, fontsize=fontsize)
        plt.savefig(figstr)
        plt.close('all')

    return distance_to_gridcell_centre 

# -----------------------------------------------------------------------------
if __name__ == "__main__":
    
    parser = OptionParser("usage: %prog [options] lat lon")
    (options, args) = parser.parse_args()
    if len(args) != 2:
        parser.error("incorrect number of arguments: please enter lat and lon separated by spaces")

    target_lat = float(args[0])
    target_lon = float(args[1])
    print('(lat,lon)=',str(target_lat),str(target_lon))
    distance_to_gridcell_centre = locate_timeseries(target_lat,target_lon)
    print('Distance to gridcell centre = ',str(round(distance_to_gridcell_centre,2)),' km')

    # ------------------------
    print('** END')

