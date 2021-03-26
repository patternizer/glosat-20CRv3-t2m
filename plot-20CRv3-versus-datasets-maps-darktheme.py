#! /usr/bin python

#------------------------------------------------------------------------------
# PROGRAM: plot-prelim-maps.py
#------------------------------------------------------------------------------
# Version 0.3
# 26 March, 2021
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

use_minmax = False
use_horizontal_colorbar = False
cmap = 'coolwarm'
fontsize = 14
cbstr = r'2m temperature anomaly [$^{\circ}$C]'

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

sourcestr_GloSAT = r'$\bf{GloSAT}$' + ' ' + r'$\bf{p03}$' + ': https://crudata.uea.ac.uk/cru/data/temperature/'        
sourcestr_CRUTEM = r'$\bf{CRUTEM}$' + ' ' + r'$\bf{5.0.1.0}$' + ': https://crudata.uea.ac.uk/cru/data/temperature/'        
sourcestr_HadCRUT = r'$\bf{HadCRUT}$' + ' ' + r'$\bf{5.0.1.0}$' + ': https://www.metoffice.gov.uk/hadobs/hadcrut5/'        
sourcestr_20CRv3 = r'$\bf{20CRv3}$' + ' ' + r'$\bf{TMP2m}$' + ': https://psl.noaa.gov/data/20thC$_{-}$Rean/'        
baselinestr = r'$\bf{Baseline}$' + ': 1961-1990'        
authorstr = r'$\bf{Graphic}$' + ': Michael Taylor, CRU/UEA' + ' -- ' + titletime

#------------------------------------------------------------------------------
# I/O: GloSAT.prelim01_reqSD_alternativegrid-178101-201912.nc
#------------------------------------------------------------------------------

# load .nc file (netcdf4) into xarray

filename_GloSAT = 'DATA/GloSAT.prelim03_reqSD_standardgrid-178101-202012.nc'
filename_CRUTEM = 'DATA/CRUTEM.5.0.1.0.anomalies.nc'
filename_HadCRUT = 'DATA/HadCRUT.5.0.1.0.anomalies.ensemble_mean.nc'
filename_20CRv3 = 'ensemble_anomalies_5x5_mean.nc'

ds_GloSAT = xr.open_dataset(filename_GloSAT, decode_cf=True) 
lat_GloSAT = np.array(ds_GloSAT.latitude)
lon_GloSAT = np.array(ds_GloSAT.longitude)
time_GloSAT = np.array(ds_GloSAT.time)
t_GloSAT = [ str(time_GloSAT[i])[0:7] for i in range(len(time_GloSAT)) ]

ds_CRUTEM = xr.open_dataset(filename_CRUTEM, decode_cf=True) 
lat_CRUTEM = np.array(ds_CRUTEM.latitude)
lon_CRUTEM = np.array(ds_CRUTEM.longitude)
time_CRUTEM = np.array(ds_CRUTEM.time)
t_CRUTEM = [ str(time_CRUTEM[i])[0:7] for i in range(len(time_CRUTEM)) ]

ds_HadCRUT = xr.open_dataset(filename_HadCRUT, decode_cf=True) 
lat_HadCRUT = np.array(ds_HadCRUT.latitude)
lon_HadCRUT = np.array(ds_HadCRUT.longitude)
time_HadCRUT = np.array(ds_HadCRUT.time)
t_HadCRUT = [ str(time_HadCRUT[i])[0:7] for i in range(len(time_HadCRUT)) ]

ds_20CRv3 = xr.open_dataset(filename_20CRv3, decode_cf=True) 
lat_20CRv3 = np.array(ds_20CRv3.lat)
lon_20CRv3 = np.array(ds_20CRv3.lon) - 180
time_20CRv3 = pd.date_range(start='1806', periods=len(np.array(ds_20CRv3.time)), freq='M') 
t_20CRv3 = [ str(time_20CRv3[i])[0:7] for i in range(len(time_20CRv3)) ]

n_GloSAT = len(t_GloSAT)    # 2880 months: 1781-01 - 2020-12
n_CRUTEM = len(t_CRUTEM)    # 2052 months: 1850-01 - 2020-12
n_HadCRUT = len(t_HadCRUT)  # 2052 months: 1850-01 - 2020-12
n_20CRv3 = len(t_20CRv3)    # 2520 months: 1806-01 - 2015-12
timeshift_CRUTEM = 828      #  828 months: 1781-01 - 1850-01
timeshift_HadCRUT = 828     #  828 months: 1781-01 - 1850-01
timeshift_20CRv3 = 300      #  300 months: 1781-01 - 1806-01

def make_plot(axi,v,vmin,vmax,cbstr,titlestr,cmap,fontsize):

    g = v.plot(ax=axi, transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax, cmap=cmap, cbar_kwargs={'orientation':'vertical','extend':'both','shrink':0.7, 'pad':0.1})         
    cb = g.colorbar; cb.ax.tick_params(labelsize=fontsize); cb.set_label(label=cbstr, size=fontsize); cb.remove()
    axi.set_global()        
    axi.coastlines(color='grey')
    axi.set_title(titlestr, fontsize=fontsize)    
    gl = axi.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='white', alpha=0.2, linestyle='-')
    gl.xlabels_top = False; gl.xlabels_bottom = False; gl.ylabels_left = False; gl.ylabels_right = False
    gl.xlines = True; gl.ylines = True
    gl.xlocator = mticker.FixedLocator(np.linspace(-180,180,73)) # every 5 degrees
    gl.ylocator = mticker.FixedLocator(np.linspace(-90,90,37))   # every 5 degrees
    gl.xformatter = LONGITUDE_FORMATTER; gl.yformatter = LATITUDE_FORMATTER

    return g

for i in range(2820,n_GloSAT):

    filestr = "temperature_anomaly" + "_" + t_GloSAT[i] + ".png"
    titlestr = 'Gridded (5x5) observations versus reanalysis: ' + r'$\bf{{{a}}}$'.format(a=t_GloSAT[i])    

    fig, axs = plt.subplots(2,2, figsize=(15,10), subplot_kw=dict(projection=p))
    # PLOT: GloSAT.p03
    titlestr_GloSAT = 'GloSAT.p03 (non-infilled)' 
    v = ds_GloSAT.temperature_anomaly[i,:,:]
    if use_minmax == True: vmin = np.min(v); vmax = np.max(v)
    else: vmin = -5.0; vmax = +5.0    
    g = make_plot(axs[0,0],v,vmin,vmax,cbstr,titlestr_GloSAT,cmap,fontsize)
    # PLOT: 20CRv3        
    titlestr_20CRv3 = '20CRv3 (regridded)'
    if ((i >= timeshift_20CRv3) & (i <= n_20CRv3)): v = ds_20CRv3.TMP2m[i-timeshift_20CRv3,0,:,:]        
    else: v = ds_20CRv3.TMP2m[0,0,:,:]*np.nan
    if use_minmax == True: vmin = np.min(v); vmax = np.max(v)
    else: vmin = -5.0; vmax = +5.0
    g = make_plot(axs[1,0],v,vmin,vmax,cbstr,titlestr_20CRv3,cmap,fontsize)
    # PLOT: CRUTEM5        
    titlestr_CRUTEM = 'CRUTEM5.0.1.0 (non-infilled)'     
    if i >= timeshift_CRUTEM: v = ds_CRUTEM.tas[i-timeshift_CRUTEM,:,:]        
    else: v = ds_20CRv3.TMP2m[0,0,:,:] * np.nan
    if use_minmax == True: vmin = np.min(v); vmax = np.max(v)
    else: vmin = -5.0; vmax = +5.0    
    g = make_plot(axs[0,1],v,vmin,vmax,cbstr,titlestr_CRUTEM,cmap,fontsize)
    # PLOT: HadCRUT5        
    titlestr_HadCRUT = 'HadCRUT5.0.1.0 (non-infilled, ensmean)'       
    if i >= timeshift_HadCRUT: v = ds_HadCRUT.tas_mean[i-timeshift_HadCRUT,:,:]        
    else: v = ds_20CRv3.TMP2m[0,0,:,:]*np.nan
    if use_minmax == True: vmin = np.min(v); vmax = np.max(v)
    else: vmin = -5.0; vmax = +5.0    
    g = make_plot(axs[1,1],v,vmin,vmax,cbstr,titlestr_HadCRUT,cmap,fontsize)

    fig.suptitle(titlestr, fontsize=30, color='white', fontweight='bold')        
    plt.annotate(sourcestr_GloSAT, xy=(200,150), xycoords='figure pixels', color='white', fontsize=fontsize) 
    plt.annotate(sourcestr_CRUTEM, xy=(200,125), xycoords='figure pixels', color='white', fontsize=fontsize) 
    plt.annotate(sourcestr_HadCRUT, xy=(200,100), xycoords='figure pixels', color='white', fontsize=fontsize) 
    plt.annotate(sourcestr_20CRv3, xy=(200,75), xycoords='figure pixels', color='white', fontsize=fontsize) 
    plt.annotate(baselinestr, xy=(200,50), xycoords='figure pixels', color='white', fontsize=fontsize)   
    plt.annotate(authorstr, xy=(900,50), xycoords='figure pixels', color='white', fontsize=fontsize, bbox=dict(boxstyle="square, pad=0.3", fc='black', edgecolor='white', linewidth=0.2))     
    fig.subplots_adjust(left=None, bottom=0.2, right=0.95, top=None, wspace=None, hspace=None)

    cb = fig.colorbar(g, ax=axs.ravel().tolist(), shrink=0.6, extend='both')
    cb.set_label(cbstr, rotation=90, labelpad=20, fontsize=fontsize)
    cb.ax.tick_params(labelsize=fontsize)
    cb.set_ticks(np.linspace(vmin,vmax,11))

    plt.savefig(filestr)
    plt.close('all')
    
#------------------------------------------------------------------------------
print('** END')
