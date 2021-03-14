#!/usr/bin/env python

#-----------------------------------------------------------------------
# PROGRAM: extract_timeseries_and_compare_absolutes.py
#-----------------------------------------------------------------------
# Version 0.2
# 12 March, 2021
# Dr Michael Taylor
# https://patternizer.github.io
# patternizer AT gmail DOT com
#-----------------------------------------------------------------------

# Dataframe libraries:
import numpy as np
import pandas as pd
import xarray as xr

# Datetime libraries:
from datetime import datetime
import nc_time_axis
import cftime

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

def locate_timeseries(station_code, station_name):

    #-----------------------------------------------------------------------------
    # SETTINGS
    #-----------------------------------------------------------------------------

    fontsize = 16
    use_anomalies = False
    use_yearly = True
    make_plot = True

    if use_anomalies == True:
        df_in = pd.read_pickle('df_anom.pkl', compression='bz2')        
        if use_yearly == True:
            dirstr = 'ANOMALIES/yearly/'
            typestr = 'anomaly_yearly'
        else:
            dirstr = 'ANOMALIES'        
            typestr = 'anomaly'    
    else:
        df_in = pd.read_pickle('df_temp.pkl', compression='bz2')
        if use_yearly == True:
            dirstr = 'ABSOLUTES/yearly/'
            typestr = 'absolute_yearly'
        else:
            dirstr = 'ABSOLUTES'    
            typestr = 'absolute'

    if not station_code:
        targets = df_in[df_in['stationname'].str.contains(station_name, case = False)].reset_index(drop=True)
        if len(targets) == 0:
            print('station name not in archive, returning ...')
            return
        maxlen = 0; i=0
        for j in range(len(targets['stationcode'].unique())):
            tslen = len(targets[targets['stationcode']==targets['stationcode'].unique()[j]])
            if tslen > maxlen:
                maxlen = tslen
                i=j
        target = targets[targets['stationcode']==targets['stationcode'].unique()[i]]
        target_lat = target['stationlat'].unique()[0]
        target_lon = target['stationlon'].unique()[0]
        targetstr = target['stationcode'].unique()[0] + '_' + (target['stationname'].unique()[0].lower()).replace(" ", "_") + '_' + typestr        
    else:
        target = df_in[df_in['stationcode']==station_code].reset_index(drop=True)
        target_lat = target['stationlat'].unique()[0]
        target_lon = target['stationlon'].unique()[0]
        targetstr = station_code + '_' + (target['stationname'].unique()[0].lower()).replace(" ", "_") + '_' + typestr 
    if target_lon < 0: target_lon += 360.0

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

    ds = xr.open_dataset(ensemble_mean, decode_cf=True)
    lat = ds.lat
    lon = ds.lon
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
    nearest_lat_idx = np.where(ds.lat==nearest_lat)[0][0] # --> 20CRv3 grid lat idx
    nearest_lon_idx = np.where(ds.lon==nearest_lon)[0][0] # --> 20CRv3 grid lon idx
    delta_lon = (x[1]-x[0])/2
    delta_lat = (y[1]-y[0])/2
    distance_to_gridcell_centre = haversine(target_lat, target_lon, nearest_lat+delta_lat, nearest_lon+delta_lon)

    #-----------------------------------------------------------------------------
    # EXTRACT GRIDCELL TIMESERIES + CONVERT TO DEGREES CENTIGRADE
    #-----------------------------------------------------------------------------

    if use_anomalies == True:
        ts_ensemble = ds.TMP2m[:,0,nearest_lat_idx,nearest_lon_idx]
    else:
        ts_ensemble = ds.TMP2m[:,0,nearest_lat_idx,nearest_lon_idx] - 273.15
    t = ds.TMP2m[:,0,nearest_lat_idx,nearest_lon_idx].time
    startyear = '1806'
    if use_yearly == True:
        t_ensemble = pd.date_range(start=startyear, periods=len(t), freq='A')   
        ts_station = np.mean(np.array(target.groupby('year').mean().iloc[:,0:11]),axis=1) 
        # FIX: pandas bug <1678 and >2262 calendar limit for long rescued station timeseries
        if target['year'].iloc[0] > 1678:
            t_station = pd.date_range(start=str(target['year'].iloc[0]), periods=len(ts_station), freq='A')   
            ts_station = pd.Series(ts_station, index=t_station)
        else:
            t_station_xr = xr.cftime_range(start=str(target['year'].iloc[0]), periods=len(ts_station), freq='A', calendar="gregorian")     
            t_station = [float(t_station_xr[i].year) for i in range(len(t_station_xr))]
            ts_station = pd.Series(ts_station, index=t_station)
    else:
        t_ensemble = pd.date_range(start=startyear, periods=len(t), freq='M')   
        ts_station = np.array(target.groupby('year').mean().iloc[:,0:11]).ravel() 
        # FIX: pandas <1678 and >2262 calendar limit for long rescued station timeseries
        if target['year'].iloc[0] > 1678:
            t_station = pd.date_range(start=str(target['year'].iloc[0]), periods=len(ts_station), freq='M')          
        else:
            t_station_xr = xr.cftime_range(start=str(target['year'].iloc[0]), periods=len(ts_station), freq='M', calendar='gregorian')     
            year = [t_station_xr[i].year for i in range(len(t_station_xr))]
            year_frac = []
            for i in range(len(t_station_xr)):
                if i%12 == 0:
                    istart = i
                    iend = istart+11   
                    frac = np.cumsum([t_station_xr[istart+j].day for j in range(12)])
                    year_frac += list(frac/frac[-1])
                else:
                    i += 1
            year_decimal = [float(year[i])+year_frac[i] for i in range(len(year))]    
            t_station = year_decimal

    #-----------------------------------------------------------------------------
    # CONVERT TO CRUTEM FORMAT + SAVE
    #-----------------------------------------------------------------------------

    if make_plot == True:

        #-----------------------------------------------------------------------------
        # PLOT: extracted 20CRv3 timeseries versus GloSAT.p03
        #-----------------------------------------------------------------------------

        figstr = outdir + '/' + targetstr + '.png'
        titlestr = '('+str(round(target_lon,3))+','+str(round(target_lat,3))+'): 20CRv3 timeseries extract at gridcell ('+str(round(nearest_lon,3))+','+str(round(nearest_lat,3))+') separation='+str(round(distance_to_gridcell_centre,3))+'km'

        fig,ax = plt.subplots(figsize=(15,10))
        plt.step(t_ensemble,ts_ensemble, color='red', label='20CRv3')
        plt.step(t_station, ts_station, color='blue', label='GloSAT')
        plt.tick_params(labelsize=fontsize)
        ax.xaxis.grid(True, which='major')        
        ax.yaxis.grid(True, which='major')        
        plt.legend(loc=2, ncol=1, fontsize=fontsize)
        plt.xlabel("Year", fontsize=fontsize)
        if use_anomalies == True:
            plt.ylabel("Temperature anomaly, $\mathrm{\degree}C$", fontsize=fontsize)
        else:
            plt.ylabel("Absolute temperature, $\mathrm{\degree}C$", fontsize=fontsize)
        plt.title(titlestr, fontsize=fontsize)
        plt.savefig(figstr)
        plt.close('all')

# -----------------------------------------------------------------------------
if __name__ == "__main__":
    
    parser = OptionParser("usage: %prog [options] station_code station_name")
    (options, args) = parser.parse_args()
    if len(args) == 0:
        parser.error("incorrect number of arguments: please enter at least a station_code and/or a station_name separated by spaces")
        station_code = []
        station_name = []
    elif len(args) == 1:
        if args[0].isdigit():
            station_code = args[0]
            station_name = []
        else:
            station_code = []
            station_name = args[0]
    elif len(args) == 2:
        if args[0].isdigit():
            station_code = args[0]
            station_name = args[1]
        else:
            station_code = args[1]
            station_name = args[0]
    else:
        parser.error("incorrect number of arguments: you entered more than a station_code and/or a station_name separated by spaces please check and re-do")
        station_code = []
        station_name = []

    locate_timeseries(station_code, station_name)

    # ------------------------
    print('** END')

