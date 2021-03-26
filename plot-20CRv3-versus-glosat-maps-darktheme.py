#! /usr/bin python

#------------------------------------------------------------------------------
# PROGRAM: plot-prelim-maps.py
#------------------------------------------------------------------------------
# Version 0.2
# 25 March, 2021
# Michael Taylor
# https://patternizer.github.io
# patternizer AT gmail DOT com
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# IMPORT PYTHON LIBRARIES
#------------------------------------------------------------------------------
import numpy as np
import numpy.ma as ma
import pandas as pd
import xarray as xr
from datetime import datetime
import calendar as cal
# Plotting libraries:
import matplotlib
import matplotlib.pyplot as plt; plt.close('all')
import matplotlib.cm as cm
from matplotlib import rcParams
from matplotlib import colors as mcol
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
import matplotlib.ticker as mticker
# Mapping libraries:
import cartopy
import cartopy.crs as ccrs
from cartopy.io import shapereader
import cartopy.feature as cf
from cartopy.util import add_cyclic_point
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

# Silence library version notifications
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

# Calculate current time for polar plots

now = datetime.now()
currentmn = str(now.month)
if now.day == 1:
    currentdy = str(cal.monthrange(now.year,now.month-1)[1])
    currentmn = str(now.month-1)
else:
    currentdy = str(now.day-1)
if int(currentdy) < 10:
    currentdy = '0' + currentdy    
currentyr = str(now.year)
if int(currentmn) < 10:
    currentmn = '0' + currentmn
titletime = str(currentdy) + '/' + currentmn + '/' + currentyr

#------------------------------------------------------------------------------
# SETTINGS: 
#------------------------------------------------------------------------------

use_obs = 'GloSAT'
#use_obs = 'CRUTEM'
#use_obs = 'HadCRUT'
use_minmax = False
use_horizontal_colorbar = False
use_mask = False
use_gridliner = True
cmap = 'coolwarm'
fontsize = 14
top_pad = 0.9

projection = 'robinson'

if projection == 'platecarree': p = ccrs.PlateCarree(central_longitude=0); threshold = 0
if projection == 'mollweide': p = ccrs.Mollweide(central_longitude=0); threshold = 1e6
if projection == 'robinson': p = ccrs.Robinson(central_longitude=0); threshold = 0
if projection == 'equalearth': p = ccrs.EqualEarth(central_longitude=0); threshold = 0
if projection == 'geostationary': p = ccrs.Geostationary(central_longitude=0); threshold = 0
if projection == 'goodehomolosine': p = ccrs.InterruptedGoodeHomolosine(central_longitude=0); threshold = 0
if projection == 'europp': p = ccrs.EuroPP(); threshold = 0
if projection == 'northpolarstereo': p = ccrs.NorthPolarStereo(); threshold = 0
if projection == 'southpolarstereo': p = ccrs.SouthPolarStereo(); threshold = 0
if projection == 'lambertconformal': p = ccrs.LambertConformal(central_longitude=0); threshold = 0

#----------------------------------------------------------------------------
# DARK BACKGROUND THEME
#----------------------------------------------------------------------------

matplotlib.rcParams['text.usetex'] = False
rcParams['font.family'] = ['DejaVu Sans']
rcParams['font.sans-serif'] = ['Avant Garde']
plt.rc('text',color='white')
plt.rc('lines',color='white')
plt.rc('patch',edgecolor='white')
plt.rc('grid',color='lightgray')
plt.rc('xtick',color='white')
plt.rc('ytick',color='white')
plt.rc('axes',edgecolor='lightgray')
plt.rc('axes',facecolor='black')
plt.rc('axes',labelcolor='white')
plt.rc('figure',facecolor='black')
plt.rc('figure',edgecolor='black')
plt.rc('savefig',edgecolor='black')
plt.rc('savefig',facecolor='black')

#----------------------------------------------------------------------------
# CREDITS
#----------------------------------------------------------------------------

if use_obs == 'GloSAT':
    titlestr = 'GloSAT.p03 (non-infilled) versus 20CRv3'        
    sourcestr1 = r'$\bf{GloSAT}$' + ' ' + r'$\bf{p03}$' + ': https://crudata.uea.ac.uk/cru/data/temperature/'        
elif use_obs == 'CRUTEM':
    titlestr = 'CRUTEM5 (non-infilled) versus 20CRv3'        
    sourcestr1 = r'$\bf{CRUTEM}$' + ' ' + r'$\bf{5.0.1.0}$' + ': https://crudata.uea.ac.uk/cru/data/temperature/'        
else:
    titlestr = 'HadCRUT5 (non-infilled ensemble mean) versus 20CRv3'        
    sourcestr1 = r'$\bf{HadCRUT}$' + ' ' + r'$\bf{5.0.1.0}$' + ': https://www.metoffice.gov.uk/hadobs/hadcrut5/'        
sourcestr2 = r'$\bf{20CRv3}$' + ' ' + r'$\bf{TMP2m}$' + ': https://psl.noaa.gov/data/20thC$_{-}$Rean/'        
baselinestr = r'$\bf{Baseline}$' + ': 1961-1990'        
authorstr = r'$\bf{Graphic}$' + ': Michael Taylor, CRU/UEA' + ' -- ' + titletime

#------------------------------------------------------------------------------
# I/O: GloSAT.prelim01_reqSD_alternativegrid-178101-201912.nc
#------------------------------------------------------------------------------

# load .nc file (netcdf4) into xarray

if use_obs == 'GloSAT':
    filename_GloSAT = 'DATA/GloSAT.prelim03_reqSD_standardgrid-178101-202012.nc'
elif use_obs == 'CRUTEM':
    filename_GloSAT = 'DATA/CRUTEM.5.0.1.0.anomalies.nc'
else:
    filename_GloSAT = 'DATA/HadCRUT.5.0.1.0.anomalies.ensemble_mean.nc'
filename_20CRv3 = 'ensemble_anomalies_5x5_mean.nc'

ds_GloSAT = xr.open_dataset(filename_GloSAT, decode_cf=True) 
lat_GloSAT = np.array(ds_GloSAT.latitude)
lon_GloSAT = np.array(ds_GloSAT.longitude)
time_GloSAT = np.array(ds_GloSAT.time)
t_GloSAT = [ str(time_GloSAT[i])[0:7] for i in range(len(time_GloSAT)) ]

ds_20CRv3 = xr.open_dataset(filename_20CRv3, decode_cf=True) 
lat_20CRv3 = np.array(ds_20CRv3.lat)
lon_20CRv3 = np.array(ds_20CRv3.lon) - 180
time_20CRv3 = pd.date_range(start='1806', periods=len(np.array(ds_20CRv3.time)), freq='M') 
t_20CRv3 = [ str(time_20CRv3[i])[0:7] for i in range(len(time_20CRv3)) ]

if use_obs == 'GloSAT':
    timeshift = -300 # months between 1806-01 (start of 20CRv3) and 1781-01 (start of GloSAT.p03)
    nstart = 2519 # --> 1806-01 (i.e. (165x12)+11 months before 2015-12)
elif use_obs == 'CRUTEM':
    timeshift = 528 # months between 1806-01 (start of 20CRv3) and 1850-01 (start of CRUTEM5)
    nstart = 1991 # --> 1850-01 (i.e. (165x12)+11 months before 2015-12)
else:
    timeshift = 528 # months between 1806-01 (start of 20CRv3) and 1850-01 (start of CRUTEM5)
    nstart = 1991 # --> 1850-01 (i.e. (165x12)+11 months before 2015-12)

#nstart = 2519 # --> 1806-01 
#nstart = 1991 # --> 1850-01 (i.e. (165x12)+11 months before 2015-12)
#nstart = 1631 # --> 1880-01 
#nstart = 1391 # --> 1900-01
#nstart = 791  # --> 1950-01
nstart = 191  # --> 2000-01

nGloSAT = 2880
timeshift_CRUTEM = 828
timeshift_HadCRUT = 828
timeshift_20CRv3 = 300

x, y = np.meshgrid(lon_GloSAT, lat_GloSAT)    

for i in range(len(time_20CRv3)-timeshift-nstart-1,len(time_20CRv3)-timeshift-nstart): 

    filestr = use_obs.lower() + '_' + "temperature_anomaly" + "_" + t_GloSAT[i] + ".png"
    if use_obs == 'GloSAT':
        titlestr1 = 'GloSAT.p03: ' + r'$\bf{{{a}}}$'.format(a=t_GloSAT[i])
    elif use_obs == 'CRUTEM':
        titlestr1 = 'CRUTEM5.0.1.0: ' + r'$\bf{{{a}}}$'.format(a=t_GloSAT[i])
    else:
        titlestr1 = 'HadCRUT5.0.1.0: ' + r'$\bf{{{a}}}$'.format(a=t_GloSAT[i])
    titlestr2 = '20CRv3: ' + r'$\bf{{{a}}}$'.format(a=t_GloSAT[i])
    colorbarstr = r'Temperature anomaly (from 1961-1990) [$^{\circ}$C]'

    fig  = plt.figure(figsize=(15,10))

    ax1 = fig.add_subplot(211, projection=p)

    if use_obs == 'GloSAT':
        v = ds_GloSAT.temperature_anomaly[i,:,:]
        make_plot = True
    elif use_obs == 'CRUTEM':    
        if i >= timeshift_CRUTEM:
            v = ds_GloSAT.tas[i-timeshift_CRUTEM,:,:]
            make_plot = True
        else:
            make_plot = False
    else:
        if i >= timeshift_HadCRUT:
            v = ds_GloSAT.tas_mean[i-timeshift_HadCRUT,:,:]
            make_plot = True
        else:
            make_plot = False
    if use_minmax == True: 
        vmin = np.min(v); vmax = np.max(v)
    else: 
        vmin = -5.0; vmax = +5.0    
    if make_plot == True:
        g = v.plot( ax = ax1, transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax, cmap=cmap, cbar_kwargs={'orientation':'vertical','extend':'both','shrink':0.8, 'pad':0.1})         
        cb = g.colorbar; cb.ax.tick_params(labelsize=fontsize); cb.set_label(label = r'T2m Anomaly [$^{\circ}$C]', size=fontsize)
    else:
        v = ds_20CRv3.TMP2m[0,0,:,:] * np.nan
        g = v.plot( ax = ax1, transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax, cmap=cmap, cbar_kwargs={'orientation':'vertical','extend':'both','shrink':0.8, 'pad':0.1})         
        cb = g.colorbar; cb.ax.tick_params(labelsize=fontsize); cb.set_label(label = r'T2m Anomaly [$^{\circ}$C]', size=fontsize)
    ax1.set_global()        
    ax1.coastlines(color='grey')
    ax1.set_title(titlestr1, fontsize=fontsize)    

    if use_gridliner == True:
        
        gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='white', alpha=0.2, linestyle='-')
        gl.xlabels_top = False; gl.xlabels_bottom = False; gl.ylabels_left = False; gl.ylabels_right = False
        gl.xlines = True; gl.ylines = True
        gl.xlocator = mticker.FixedLocator(np.linspace(-180,180,73)) # every 5 degrees
        gl.ylocator = mticker.FixedLocator(np.linspace(-90,90,37))   # every 5 degrees
        gl.xformatter = LONGITUDE_FORMATTER; gl.yformatter = LATITUDE_FORMATTER

    if use_mask == True:

        g = ccrs.Geodetic()
        trans = ax1.projection.transform_points(g, x, y)
        x0 = trans[:,:,0]
        x1 = trans[:,:,1]    
        for mask in (x0>threshold,x0<=threshold):
            im1 = ax1.pcolor(ma.masked_where(mask, x0), ma.masked_where(mask, x1), ma.masked_where(mask, v), vmin=vmin, vmax=vmax, transform=ax1.projection, cmap=cmap) 
        im1.set_clim(vmin,vmax)    
        
    ax2 = fig.add_subplot(212, projection=p)

    if i >= timeshift_20CRv3:
        v = ds_20CRv3.TMP2m[i-timeshift_20CRv3,0,:,:]        
        make_plot = True
    else:            
        make_plot = False
    if use_minmax == True:    
        vmin = np.min(v); vmax = np.max(v)
    else:
        vmin = -5.0; vmax = +5.0
    if make_plot == True:
        g = v.plot( ax = ax2, transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax, cmap=cmap, cbar_kwargs={'orientation':'vertical','extend':'both','shrink':0.8, 'pad':0.1})         
        cb = g.colorbar; cb.ax.tick_params(labelsize=fontsize); cb.set_label(label = r'T2m Anomaly [$^{\circ}$C]', size=fontsize)
    else:
        v = ds_20CRv3.TMP2m[0,0,:,:] * np.nan
        g = v.plot( ax = ax2, transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax, cmap=cmap, cbar_kwargs={'orientation':'vertical','extend':'both','shrink':0.8, 'pad':0.1})         
        cb = g.colorbar; cb.ax.tick_params(labelsize=fontsize); cb.set_label(label = r'T2m Anomaly [$^{\circ}$C]', size=fontsize)
    ax2.set_global()
    ax2.coastlines(color='grey')
    ax2.set_title(titlestr2, fontsize=fontsize)    

    if use_gridliner == True:
        gl = ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='black', alpha=0.2, linestyle='-')
        gl.xlabels_top = False; gl.xlabels_bottom = False; gl.ylabels_left = False; gl.ylabels_right = False
        gl.xlines = True; gl.ylines = True
        gl.xlocator = mticker.FixedLocator(np.linspace(-180,180,73)) # every 5 degrees
        gl.ylocator = mticker.FixedLocator(np.linspace(-90,90,37))   # every 5 degrees
        gl.xformatter = LONGITUDE_FORMATTER; gl.yformatter = LATITUDE_FORMATTER

    if use_mask == True:
        g = ccrs.Geodetic()
        trans = ax2.projection.transform_points(g, x, y)
        x0 = trans[:,:,0]
        x1 = trans[:,:,1]    
        for mask in (x0>threshold,x0<=threshold):
            im2 = ax2.pcolor(ma.masked_where(mask, x0), ma.masked_where(mask, x1), ma.masked_where(mask, v), vmin=vmin, vmax=vmax, transform=ax2.projection, cmap=cmap) 
        im2.set_clim(vmin,vmax)    

    fig.suptitle(titlestr, fontsize=24, color='white', fontweight='bold')        
    plt.annotate(sourcestr1, xy=(200,110), xycoords='figure pixels', color='white', fontsize=fontsize) 
    plt.annotate(sourcestr2, xy=(200,80), xycoords='figure pixels', color='white', fontsize=fontsize) 
    plt.annotate(baselinestr, xy=(200,50), xycoords='figure pixels', color='white', fontsize=fontsize)   
    plt.annotate(authorstr, xy=(200,20), xycoords='figure pixels', color='white', fontsize=fontsize, bbox=dict(boxstyle="square, pad=0.3", fc='black', edgecolor='white', linewidth=0.2))     
    fig.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=None, hspace=None)
    plt.savefig(filestr)
    plt.close('all')
    
#------------------------------------------------------------------------------
print('** END')
