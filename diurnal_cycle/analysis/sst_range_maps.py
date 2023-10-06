import sys
import os
import re
import xarray as xr
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.dates as mdates
from matplotlib import cm, gridspec, rcParams, colors
from matplotlib.offsetbox import AnchoredText
from mpl_toolkits.axes_grid1 import AxesGrid
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.geoaxes import GeoAxes
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import sys
import utils

ddir = '/lustre/f2/dev/ncep/Jack.Reeveseyre/diurnal_cycle/'
yyyy_mm = '2005_01' # e.g., '2005_01'
ds = xr.open_dataset(ddir + 'ocn_' + yyyy_mm + '_grid_meanDiurnalCycle.nc')
sstdr = ds.sel(st_ocean=0.5).max(dim='hour') - ds.sel(st_ocean=0.5).min(dim='hour')

def switch_lon_lims(lon_list, min_lon=0.0):
    result = (lon_list - min_lon) % 360.0 + min_lon
    return result

fig = plt.figure(figsize=(8, 6))
data_crs = ccrs.PlateCarree()
plot_crs = ccrs.Robinson(central_longitude=-150)
ax = fig.add_subplot(111, projection=plot_crs)
ax.set_global()
gl = ax.gridlines(draw_labels=True,
		  xlocs=np.arange(-300, 151, 60),
		  ylocs=np.arange(-90, 91, 30),
		  linewidths=0.3)
gl.right_labels = True
gl.xformatter = LongitudeFormatter(zero_direction_label=False,
				   degree_symbol='')
gl.yformatter = LatitudeFormatter(degree_symbol='')
ax.coastlines(zorder=5)
ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=4)

cmap_b = plt.get_cmap('magma')
norm_b = colors.BoundaryNorm(np.arange(0.0, 1.01, 0.1),
			     ncolors=cmap_b.N)

p = ax.pcolormesh(switch_lon_lims(sstdr['geolon_t'].data, -180.0),
		          sstdr['geolat_t'].data,
		          sstdr['temp'].isel(time=0).data,
		          transform=data_crs,
		          cmap=cmap_b, norm=norm_b, shading='nearest')
ax.set_title('TEMPERATURE DIURNAL CYCLE RANGE | ' + yyyy_mm)

fig.colorbar(p,
	     label='SST DIURNAL RANGE (K)',
	     ax=ax, orientation='horizontal',shrink=0.5)

# Save file.
plotdir = '/ncrc/home2/Jack.Reeveseyre/cfs-analysis/diurnal_cycle/plots/'
plotfn = 'map_SST_RANGE_' + yyyy_mm + '.'
plotfileformat='png'
plt.savefig(plotdir + plotfn + plotfileformat, format=plotfileformat,
	    dpi=400)
