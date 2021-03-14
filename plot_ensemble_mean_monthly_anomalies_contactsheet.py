#! /usr/bin python

#------------------------------------------------------------------------------
# PROGRAM: plot_ensemble_mean_monthly_anomalies_contactsheet.py
#------------------------------------------------------------------------------
# Version 0.2
# 7 March, 2021
# Michael Taylor
# michael DOT a DOT taylor AT uea DOT ac DOT uk 
#------------------------------------------------------------------------------

import os, glob
import imageio
import numpy as np
import numpy.ma as ma
import pandas as pd
import xarray as xr
from datetime import datetime
import calendar as cal
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import colors as mcol
from matplotlib.cm import ScalarMappable
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.io import shapereader
import cartopy.feature as cf
from cartopy.util import add_cyclic_point
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import seaborn

import warnings
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

#----------------------------------------------------------------------------
# SETTINGS
#----------------------------------------------------------------------------

fontsize = 14
cmap = 'RdBu_r'
# cmap = 'gist_earth', # green-brown-white
# cmap = 'gist_yarg',  # grey-black (high contrast)
# cmap = 'gist_ncar', # lime-orange-white (high contrast)
# cmap = 'nipy_spectral', # teal-orange-lightgrey (high contrast)
#projection = 'Robinson'
projection = 'Orthographic'
make_gif = False

#----------------------------------------------------------------------------
# DARK BACKGROUND THEME
#----------------------------------------------------------------------------

matplotlib.rcParams['text.usetex'] = False
#rcParams['font.family'] = 'sans-serif'
rcParams['font.family'] = ['DejaVu Sans']
rcParams['font.sans-serif'] = ['Avant Garde']
plt.rc('text',color='black')
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
# PLOT: Ensemble mean temperature anomalies
#----------------------------------------------------------------------------

ds = xr.open_dataset('ABSOLUTES/ensemble_mean_absolute.nc', decode_cf=True)
par = ds.TMP2m 
lon = ds.lon 
lat = ds.lat 
time = ds.time

dn = xr.open_dataset('NORMALS/ensemble_mean_normals.nc', decode_cf=True)
normals = dn.TMP2m 

N = np.floor(par.shape[0]/12).astype(int)
for i in range(N):
        
    if projection == 'Orthographic':            
        p = ccrs.Orthographic(0, 0); threshold=0
    elif projection == 'Robinson':        
        p = ccrs.Robinson(central_longitude=0); threshold=0
        
    titlestr = str(time[i*12])[35:39]        
    figstr = 't2m_' + str(time[i*12])[35:39] +'.png'        
    datastr = r'$\bf{Data}$' + ': 20CRv3 T2m'        
    sourcestr = r'$\bf{Source}$' + ': https://psl.noaa.gov/data/20thC$_{-}$Rean/'        
    baselinestr = r'$\bf{Baseline}$' + ': 1961-1990'        
    authorstr = r'$\bf{Graphic}$' + ': Michael Taylor, CRU/UEA (@MichaelTaylorEO)' + ' -- ' + titletime

    fig, axs = plt.subplots(3,4, figsize=(15,10), subplot_kw=dict(projection=p))                
    for j in range(12):

        if j == 0: r=0; c=0
        elif j == 1: r=0; c=1
        elif j == 2: r=0; c=2
        elif j == 3: r=0; c=3
        elif j == 4: r=1; c=0
        elif j == 5: r=1; c=1
        elif j == 6: r=1; c=2
        elif j == 7: r=1; c=3
        elif j == 8: r=2; c=0
        elif j == 9: r=2; c=1
        elif j == 10: r=2; c=2
        elif j == 11: r=2; c=3

        v = par[i*12+j,0,:,:] - normals[j,0,:,:]
        vmin = -8.0
        vmax = 8.0
        x, y = np.meshgrid(lon,lat)        
        g = v.plot( ax = axs[r,c], transform=ccrs.PlateCarree(), vmin = vmin, vmax = vmax, cmap=cmap, cbar_kwargs={'orientation':'vertical','extend':'both','shrink':0.8, 'pad':0.1, 'label':r'T2m Anomaly [$^{\circ}$C]'}) 
        cb = g.colorbar
#       if (j != 3) & (j != 7) & (j != 11):
#           cb.remove()   
#       cb.set_label(label='T2m Anomaly [°C]', fontsize=fontsize)
        cb.ax.tick_params(labelsize=fontsize)

        if projection == 'Orthographic':            
            top_pad = 0.85       
        elif projection == 'Robinson':        
#           axs[r,c].set_extent([-180,180,-90,90], ccrs.PlateCarree())
            top_pad = 1.0
               
        parallels = np.arange(-90,90,30)
        meridians = np.arange(-180,180,30)
        gl = axs[r,c].gridlines(crs=ccrs.PlateCarree(), xlocs=meridians, ylocs=parallels, linestyle="dotted", linewidth=0.7, color='black', alpha=0.5)
        axs[r,c].add_feature(cf.LAND, facecolor='linen')
        axs[r,c].add_feature(cf.COASTLINE, edgecolor="black", linewidth=0.7)
        axs[r,c].add_feature(cf.BORDERS, edgecolor="black", linewidth=0.5)        
#       axs[r,c].set_title(str(time[i*12+j])[35:39]+'-'+str(j+1).zfill(2), fontsize=fontsize, color='white', fontweight='bold') 
        axs[r,c].set_title(str(j+1).zfill(2), fontsize=fontsize, color='white', y=1.0, fontweight='bold') 
        
    fig.suptitle(titlestr, fontsize=36, color='white', fontweight='bold')        
    plt.annotate(datastr, xy=(150,90), xycoords='figure pixels', color='white', fontsize=fontsize) 
    plt.annotate(sourcestr, xy=(150,60), xycoords='figure pixels', color='white', fontsize=fontsize) 
    plt.annotate(baselinestr, xy=(150,30), xycoords='figure pixels', color='white', fontsize=fontsize)   
    plt.annotate(authorstr, xy=(750,30), xycoords='figure pixels', color='white', fontsize=fontsize, bbox=dict(boxstyle="square, pad=0.3", fc='black', edgecolor='white', linewidth=0.2))     
    fig.subplots_adjust(left=None, bottom=1.0-top_pad, right=None, top=top_pad, wspace=None, hspace=None)
    plt.savefig(figstr)
    plt.close()

#if make_gif == True:

#    images = sorted(glob.glob('t2m_*.png'))
#    var = [imageio.imread(file) for file in images]
#    imageio.mimsave('t2m.gif', var, fps = 10)

# COVERT GIF to MP4
# ffmpeg -i t2m.gif -movflags faststart -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" t2m.mp4

# -----------------------------------------------------------------------------
print('** END')
