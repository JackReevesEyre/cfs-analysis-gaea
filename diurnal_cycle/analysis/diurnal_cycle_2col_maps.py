"""Plots multiple maps of diurnal cycle metrics on the same page.
  
All metrics. Two columns, one month per column. 

Usage:
    $ python diurnal_cycle_2col_maps.py plot_var plot_month1 plot_month2
where:
    plot_month1 and plot_month2 are individual months (1-12) 
    plot_var is one of {'TEMP', 'CUR'}
"""
#------------------------------------------------------------------------------
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
from diurnal_range_maps import plot_dets, approx_local_time, switch_lon_lims
#------------------------------------------------------------------------------

def main(plot_var, plot_month1, plot_month2):
    
    """ For reference:
    Atmo dimensions:    (time: 1, hour: 24, height: 1, lat: 190, lon: 384)
    Ocean dimensions:   (time: 1, hour: 24, st_ocean: 20, yt_ocean: 410, xt_ocean: 720,
                         yu_ocean: 410, xu_ocean: 720)
    """

    # Define some useful constants.
    ddir = '/lustre/f2/dev/ncep/Jack.Reeveseyre/diurnal_cycle/'
    atmo_chunks = {'time':48, 'hour':24, 'height':1, 'lat':190, 'lon':384}
    ocn_chunks = {'time':48, 'hour':24, 'st_ocean':20, 'yt_ocean':10,
                  'xt_ocean':20, 'yu_ocean':10, 'xu_ocean':20}
    
    # Specify depth used as reference in current shear calculation.
    curr_shear_depth = 60.0
    
    # Construct month strings.
    month_names = ['JAN','FEB','MAR','APR','MAY','JUN',
                   'JUL','AUG','SEP','OCT','NOV','DEC']
    month_str1 = month_names[int(plot_month1) - 1]
    month_str2 = month_names[int(plot_month2) - 1]
    
    # Construct file name and open it if it exists.
    dicy_fn1 = 'meanDiurnalCycle_metrics_' + plot_var + '_' + month_str1 + '.nc'
    dicy_fn2 = 'meanDiurnalCycle_metrics_' + plot_var + '_' + month_str2 + '.nc'
    ds1 = xr.open_dataset(ddir + dicy_fn1, decode_timedelta=False)
    ds2 = xr.open_dataset(ddir + dicy_fn2, decode_timedelta=False)

    # Convert all velocities to cm s-1.
    if plot_var == 'CUR':
        for v in ['RANGE', 'RANGE_25']:
            ds[v] = ds[v]*100.0
            ds[v].attrs['units'] = 'cm s-1'
    
    #----------------------------------------------------------------------------
    # Plot figure.
    
    fig_letters = 'abcdefghijklmnopqrstuvwxyz'
    
    # Define whether and where to mask fields according to diurnal range.
    range_mask = 0.0 #plot_dets(plot_var, 'mask_lim')
    
    # Set up figure and axes.
    data_crs = ccrs.PlateCarree()
    plot_crs = ccrs.Robinson(central_longitude=-150)
    axes_class = (GeoAxes,
                  dict(map_projection=plot_crs))
    fig = plt.figure(figsize=(12, 8))
    axgr = AxesGrid(fig, 111, axes_class=axes_class,
                    nrows_ncols=(4, 2),
                    direction='column',
                    axes_pad=0.25,
                    cbar_location='right',
                    cbar_mode='each',
                    cbar_pad=0.1,
                    cbar_size='2%',
                    label_mode='')  # note the empty label_mode
    
    # Calculate lat and lon meshes
    lonm = switch_lon_lims(ds1[plot_dets(plot_var, 'lonname')].data, -180.0)
    latm = ds1[plot_dets(plot_var, 'latname')].data
    
    # Create list of data arrays to plot - dictates order.
    # Note with direction='column', index goes down 1st column, then 2nd col etc.
    ds_list = [ds1, ds1, ds1, ds1,
               ds2, ds2, ds2, ds2]
    ptype_list = ['RANGE', 'PHASE', 'THRESH_DEPTH', 'THRESH_DELAY',
                  'RANGE', 'PHASE', 'THRESH_DEPTH', 'THRESH_DELAY']
    
    # Loop over axes and plot each.
    for i, ax in enumerate(axgr):
        
        # Select data for this map.
        ptype = ptype_list[i]
        ds = ds_list[i]
        if i < 4:
            mon = month_str1
        else:
            mon = month_str2
            
        # Set plot details.
        ax.set_global()
        gl = ax.gridlines(draw_labels=True,
                          xlocs=np.arange(-300, 151, 60),
                          ylocs=np.arange(-90, 91, 30),
                          linewidths=0.3)
        if i in [0,4]:
            gl.top_labels = True
            gl.bottom_labels = True
        else:
            gl.top_labels = False
            gl.bottom_labels = True
        gl.left_labels = True
        gl.xformatter = LongitudeFormatter(zero_direction_label=False,
                                           degree_symbol='')
        gl.xlabel_style = {'size':10, 'rotation':60}
        gl.ylabel_style = {'size':8, 'rotation':0}
        gl.yformatter = LatitudeFormatter(degree_symbol='')
        ax.coastlines(zorder=5)
        ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=4)
        
        # Define color map.
        cmap_b = plt.get_cmap(plot_dets(plot_var, 'cmap', ptype))
        norm_b = colors.BoundaryNorm(
            np.arange(plot_dets(plot_var, 'min', ptype),
                      plot_dets(plot_var, 'max', ptype)*1.01,
                      plot_dets(plot_var, 'step', ptype)),
            ncolors=cmap_b.N
        )
        
        # Plot the data.
        if ((ptype == 'RANGE') or (range_mask == 0.0)):
            # Unmasked version.
            p = ax.pcolormesh(lonm, latm,
	        	      ds[ptype].data,
                              transform=data_crs,
                              cmap=cmap_b, norm=norm_b, shading='nearest')
        else:
            # Masked version.
            import numpy.ma as ma
            dsm = ma.masked_where(ds['RANGE'].data < range_mask, ds[ptype].data)
            ptr = ax.pcolormesh(lonm, latm,
	        	        ds[ptype].data, alpha=0.25,
                                transform=data_crs,
                                cmap=cmap_b, norm=norm_b, shading='nearest')
            p = ax.pcolormesh(lonm, latm,
	        	      dsm,
                              transform=data_crs,
                              cmap=cmap_b, norm=norm_b, shading='nearest')
    
        # Titles.
        if ptype in ['RANGE']:
            ax.set_title(plot_dets(plot_var, 'plotname') + ', ' + mon)
        else:
            pass

        # Colorbars.
        if i < 4:
            cb = axgr.cbar_axes[i].colorbar(p)
            cb.remove()
        else:
            cb = axgr.cbar_axes[i].colorbar(p)
            cb.set_label(ptype +
                         ' (' + plot_dets(plot_var, 'units', ptype) + ')')
        
        # Add figure letters (a, b, c, etc.).
        ax.text(0.1, 1.0, fig_letters[i],
                transform=ax.transAxes,
                horizontalalignment='center', verticalalignment='center',
                bbox={'facecolor':'white', 'alpha':0.9, 'pad':0.1,
                      'boxstyle':'round'})
    
    # Save file.
    plotdir = '/ncrc/home2/Jack.Reeveseyre/cfs-analysis/diurnal_cycle/plots/'
    plotfn = 'maps_2col_' + plot_var + '_' + month_str1 + '_' + month_str2 + '.'
    plotfileformat='png'
    plt.savefig(plotdir + plotfn + plotfileformat, format=plotfileformat,
                dpi=400)
    #
    return


if __name__ == "__main__": 
    main(sys.argv[1], sys.argv[2], sys.argv[3])

