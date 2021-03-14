#! /usr/bin python

#------------------------------------------------------------------------------
# PROGRAM: plot_ensemble_mean_monthly_anomalies_months1-12.py
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

fontsize = 16
cmap = 'RdBu_r'
# cmap = 'gist_earth', # green-brown-white
# cmap = 'gist_yarg',  # grey-black (high contrast)
# cmap = 'gist_ncar', # lime-orange-white (high contrast)
# cmap = 'nipy_spectral', # teal-orange-lightgrey (high contrast)
projection = 'Robinson'
#projection = 'Orthographic'
make_gif = False

#----------------------------------------------------------------------------
# DARK BACKGROUND THEME
#----------------------------------------------------------------------------

matplotlib.rcParams['text.usetex'] = False
rcParams['font.family'] = 'sans-serif'
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
N=1
for i in range(N):

    for j in range(12):

        titlestr = str(time[i*12+j])[35:39]+'-'+str(j+1).zfill(2)

        fig = plt.figure(figsize=(15,10))  

        if projection == 'Orthographic':
            p = ccrs.Orthographic(0, 0); threshold=0
            axs = plt.axes(projection=p)
            top_pad = 0.8
        elif projection == 'Robinson':
            p = ccrs.Robinson(central_longitude=0); threshold=0
            axs = plt.axes(projection=p)
            axs.set_extent([-180,180,-90,90], ccrs.PlateCarree())
            top_pad = 1.0

        parallels = np.arange(-90,90,30)
        meridians = np.arange(-180,180,30)
        gl = axs.gridlines(crs=ccrs.PlateCarree(), xlocs=meridians, ylocs=parallels, linestyle="dotted", linewidth=0.7, color='black', alpha=0.5)   

        x, y = np.meshgrid(lon,lat)        
        v = par[i*12+j,0,:,:] - normals[j,0,:,:]
        vmin = -8.0; vmax = 8.0
        g = v.plot( ax = axs, transform=ccrs.PlateCarree(), vmin = vmin, vmax = vmax, cmap=cmap, cbar_kwargs={'orientation':'vertical','extend':'both','shrink':0.6, 'pad':0.1})                   
        cb = g.colorbar   
        cb.set_label(label=r'2m Temperature Anomaly [$^{\circ}$C]', fontsize=fontsize)
        cb.ax.tick_params(labelsize=fontsize)
        axs.add_feature(cf.LAND, facecolor='linen')
        axs.add_feature(cf.COASTLINE, edgecolor="black", linewidth=0.7)
        axs.add_feature(cf.BORDERS, edgecolor="black", linewidth=0.5)        
        axs.set_title(titlestr, fontsize=36, color='white', y=1.08, fontweight='bold')
#       fig.suptitle(titlestr, fontsize=36, color='white', fontweight='bold')
        plt.annotate(r'$\bf{Data}$' + ': 20CRv3', xy=(150,150), xycoords='figure pixels', color='white', fontsize=fontsize) 
        plt.annotate(r'$\bf{Source}$' + ': https://psl.noaa.gov/data/20thC$_{-}$Rean/', xy=(150,120), xycoords='figure pixels', color='white', fontsize=fontsize) 
        plt.annotate(r'$\bf{Baseline}$' + ': 1961-1990', xy=(150,90), xycoords='figure pixels', color='white', fontsize=fontsize) 
        plt.annotate(r'$\bf{Graphic}$' + ': Michael Taylor, CRU/UEA (@MichaelTaylorEO)' + ' -- ' + titletime, xy=(150,60), xycoords='figure pixels', color='white', fontsize=fontsize, bbox=dict(boxstyle="square, pad=0.3", fc='black', edgecolor='white', linewidth=0.2))     
        fig.subplots_adjust(top=top_pad, bottom=1.0-top_pad)
        plt.savefig('t2m_' + str(time[i*12+j])[35:39] + '_' + str(j+1).zfill(2) + '.png')
        plt.close()

if make_gif == True:

    for j in range(12):

        images = sorted(glob.glob('PLOTS/Separate-months-Robinson-RdBu_r/PNG/t2m_*' + str(j+1).zfill(2) + '.png'))
        var = [imageio.imread(file) for file in images]
        imageio.mimsave('t2m_' + str(j+1).zfill(2) + '.gif', var, fps = 10)
        
    images = sorted(glob.glob('PLOTS/Separate-months-Robinson-RdBu_r/PNG/t2m_*' + '.png'))   
    var = [imageio.imread(file) for file in images]    
    imageio.mimsave('t2m.gif', var, fps = 10)

# COVERT GIF to MP4
# ffmpeg -i t2m.gif -movflags faststart -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" t2m.mp4

# -----------------------------------------------------------------------------
print('** END')
