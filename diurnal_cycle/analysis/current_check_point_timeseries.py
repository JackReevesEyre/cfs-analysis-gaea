"""Plots maps of current and wind components.
Usage:
    $ python current_check_point_timeseries.py plot_month plot_var
where:
    plot_month is 'ANN' 
               or an individual month (1-12) 
               or a comma separated list of months (e.g., "1,2,12")
    plot_var is one of {u, v, tau_x, tau_y}
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
#----------------
from diurnal_range_maps import switch_lon_lims, get_filelist
#------------------------------------------------------------------------------

def main(plot_month, plot_var):
    
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
    
    # Construct month string.
    if plot_month == 'ANN':
        month_str = 'ANN'
    else:
        month_names = ['JAN','FEB','MAR','APR','MAY','JUN',
                       'JUL','AUG','SEP','OCT','NOV','DEC']
        plot_month_list = np.array(plot_month.split(','), dtype=int)
        if len(plot_month_list) == 1:
            month_str = month_names[plot_month_list[0] - 1]
        else:
            month_str = ''.join([month_names[m-1][0]
                                 for m in plot_month_list])
            if month_str == 'JFD':
                month_str = 'DJF'
    
    # Load data.
    lon_range = 5.0
    if plot_var in ['tau_x', 'tau_y']:
        file_list = get_filelist('surface', plot_month, ddir)
        ds = xr.open_mfdataset(file_list, chunks=ocn_chunks)[[plot_var]]
        ds.attrs['input_files'] = file_list
        dsr = ds.mean(dim='time', keep_attrs=True).max(dim='hour', keep_attrs=True) - \
            ds.mean(dim='time', keep_attrs=True).min(dim='hour', keep_attrs=True)
        max_point = dsr[plot_var].where(dsr[plot_var] == dsr[plot_var].max(), drop=True)
        plot_title = 'diurnal cycle of ' + plot_dets(plot_var, 'plotname') + ', ' + month_str
        ds = ds[plot_var].sel(
            xu_ocean=max_point.xu_ocean.data[0],                                                           
            yu_ocean=slice(max_point.yu_ocean.data[0] - lon_range,
                           max_point.yu_ocean.data[0] + lon_range)
            ).mean(dim='time').compute()
    elif plot_var in ['u', 'v']:
        file_list = get_filelist('ocn', plot_month, ddir)
        ds = xr.open_mfdataset(file_list, chunks=ocn_chunks)[[plot_var]]
        ds.attrs['input_files'] = file_list
        dsr = ds.isel(st_ocean=0).mean(dim='time', keep_attrs=True).max(dim='hour', keep_attrs=True) - \
            ds.isel(st_ocean=0).mean(dim='time', keep_attrs=True).min(dim='hour', keep_attrs=True)
        max_point = dsr[plot_var].where(dsr[plot_var] == dsr[plot_var].max(), drop=True)
        ds = ds.sel(
            xu_ocean=max_point.xu_ocean.data[0], 
            yu_ocean=slice(max_point.yu_ocean.data[0] - lon_range,
                           max_point.yu_ocean.data[0] + lon_range),
            st_ocean=0.5).mean(dim='time').compute()
        ds[plot_var] = ds[plot_var]*100.0
        ds[plot_var].attrs['units'] = 'cm s-1'
        ds[plot_var].attrs['long_name'] = {'u':'eastward current', 
                                           'v':'northward current'}[plot_var]
        plot_title = 'diurnal cycle of ' + \
            plot_dets(plot_var, 'plotname') + ' at ' + \
            str(ds.st_ocean.data) + ' m, ' + month_str
    else:
        sys.exit('plot_var not recognized.')
    
    # Set up figure and axes.
    fig, ax = plt.subplots(1, 1, figsize=(7, 5))
    ax.set_xlim(-0.5, 23.5)
    ax.set_xlabel('hour (UTC)')

    # Define color map.
    cmap_b = cm.get_cmap('tab20')

    # Plot data.
    for ilat, lat in enumerate(ds.yu_ocean.data):
        lnwid = 1.5 if lat == max_point.yu_ocean.data[0] else 0.75
        ax.plot(ds.hour.data,
                ds[plot_var].sel(yu_ocean=lat).data,
                color=cmap_b(ilat), linewidth=lnwid,
                label="{:+.2f}".format(ds.geolat_c.sel(yu_ocean=lat).data))
    
    # Finish up figure.
    ax.set_title(plot_title, loc='left')
    ax.set_title('longitude: ' + 
                 str(ds.geolon_c.sel(yu_ocean=max_point.yu_ocean.data[0])), 
                 loc='right')
    ax.legend(ncol=1, loc='right')
    
    # Save file.
    plotdir = '/ncrc/home2/Jack.Reeveseyre/cfs-analysis/diurnal_cycle/plots/'
    plotfn = 'current_check_point_timeseries' + \
        '_' + plot_var + '_' + month_str + '.'
    plotfileformat='png'
    plt.savefig(plotdir + plotfn + plotfileformat, format=plotfileformat,
                dpi=400)
    #
    return


def plot_dets(vn, att):
    det_dict = {
        'tau_x':{
	    'latname':'geolat_c',
	    'lonname':'geolon_c',
	    'plotname':'eastward wind stress',
	    'min':-0.3,
	    'max':0.3,
	    'step':0.05,
	    'units':'N m-2',
	    'cmap':'BrBg'
	},
        'tau_y':{
	    'latname':'geolat_c',
	    'lonname':'geolon_c',
	    'plotname':'northward wind stress',
	    'min':-0.2,
	    'max':0.2,
	    'step':0.025,
	    'units':'N m-2',
	    'cmap':'BrBg'
	},
        'u':{
	    'latname':'geolat_c',
	    'lonname':'geolon_c',
	    'plotname':'eastward current',
	    'min':-100.0,
	    'max':100.0,
	    'step':20.0,
	    'units':'cm s-1',
	    'cmap':'BrBg'
	},
        'v':{
	    'latname':'geolat_c',
	    'lonname':'geolon_c',
	    'plotname':'northward current',
	    'min':-30.0,
	    'max':30.0,
	    'step':5.0,
	    'units':'cm s-1',
	    'cmap':'BrBg'
	}
    }
    #
    result = det_dict[vn][att]
    #
    return result


if __name__ == "__main__": 
    main(sys.argv[1], sys.argv[2])
