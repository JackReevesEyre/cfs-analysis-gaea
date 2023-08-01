"""Plots maps of current and wind components.
Usage:
    $ python current_daily_range_maps.py plot_date plot_N plot_var
where:
    plot_date is formatted like yyyy_mm
    plot_N is an integer in the range [1,31]
    plot_var is one of {u, v}
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

def load_ds(plot_date, plot_N, plot_var):

    """ For reference:
    Atmo dimensions:    (time: 1, hour: 24, height: 1, lat: 190, lon: 384)
    Ocean dimensions:   (time: 1, hour: 24, st_ocean: 20, yt_ocean: 410, xt_ocean: 720,
                         yu_ocean: 410, xu_ocean: 720)
    """

    # Define some useful constants.
    data_dir = '/lustre/f2/dev/ncep/Jack.Reeveseyre/CFSm501hr/CFSm501hr1980010200/DATA/'
    fn_root = 'ocn_' + plot_date
    atmo_chunks = {'time':48, 'hour':24, 'height':1, 'lat':190, 'lon':384}
    ocn_chunks = {'time':48, 'hour':24, 'st_ocean':20, 'yt_ocean':10,
                  'xt_ocean':20, 'yu_ocean':10, 'xu_ocean':20}

    # Load data.
    N_days = np.int64(plot_N)
    if plot_var in ['u', 'v']:
        ds = xr.open_mfdataset(data_dir + fn_root + '*.nc',
                               chunks=ocn_chunks)[[plot_var]]
    else:
        sys.exit('plot_var not recognized.')

    return ds


def main(ds_in, plot_date, plot_N, plot_var):
    
    """ For reference:
    Atmo dimensions:    (time: 1, hour: 24, height: 1, lat: 190, lon: 384)
    Ocean dimensions:   (time: 1, hour: 24, st_ocean: 20, yt_ocean: 410, xt_ocean: 720,
                         yu_ocean: 410, xu_ocean: 720)
    """

    # Define some useful constants.
    data_dir = '/lustre/f2/dev/ncep/Jack.Reeveseyre/CFSm501hr/CFSm501hr1980010200/DATA/' 
    fn_root = 'ocn_' + plot_date 
    atmo_chunks = {'time':48, 'hour':24, 'height':1, 'lat':190, 'lon':384}
    ocn_chunks = {'time':48, 'hour':24, 'st_ocean':20, 'yt_ocean':10,
                  'xt_ocean':20, 'yu_ocean':10, 'xu_ocean':20}
    N_days = np.int64(plot_N)
    
    # Process data.
    ds = ds_in.copy()
    if plot_var in ['u', 'v']:
        ds = ds.isel(st_ocean=0, time=slice(0, N_days*24))\
                .groupby('time.hour').mean(dim='time', keep_attrs=True)
        ds = ds.max(dim='hour', keep_attrs=True) - \
            ds.min(dim='hour', keep_attrs=True)
        ds[plot_var] = ds[plot_var]*100.0
        ds[plot_var].attrs['units'] = 'cm s-1'
        ds[plot_var].attrs['long_name'] = {'u':'eastward current',
                                           'v':'northward current'}[plot_var]
        print(ds)
        plot_title = plot_dets(plot_var, 'plotname') + ' ' +\
            plot_N + '-day'\
            ' diurnal range at ' + \
            str(ds.st_ocean.data) + ' m, ' + plot_date
    else:
        sys.exit('plot_var not recognized.')
    
    # Set up figure and axes.
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
    
    # Define color map.
    cmap_b = plt.get_cmap(plot_dets(plot_var, 'cmap'))
    norm_b = colors.BoundaryNorm(
        np.arange(plot_dets(plot_var, 'min'),
                  plot_dets(plot_var, 'max')*1.01,
                  plot_dets(plot_var, 'step')),
        ncolors=cmap_b.N
    )
    p = ax.pcolormesh(switch_lon_lims(ds[plot_dets(plot_var, 'lonname')].data,
                                      -180.0),
                      ds[plot_dets(plot_var, 'latname')].data,
	                  ds[plot_var].data,
                      transform=data_crs,
                      cmap=cmap_b, norm=norm_b, shading='nearest')
    
    # Finish up figure.
    ax.set_title(plot_title)
    fig.colorbar(p,
                 label='diurnal range (' +
		         plot_dets(plot_var, 'units') + ')',
                 extend='max',
                 ax=ax, orientation='horizontal',shrink=0.5)
    
    # Save file.
    plotdir = '/ncrc/home2/Jack.Reeveseyre/cfs-analysis/diurnal_cycle/plots/'
    plotfn = 'current_' + plot_N + 'day_range_map_' + \
        plot_var + '_' + plot_date + '.'
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
	    'min':0.0,
	    'max':50.0,
	    'step':5.0,
	    'units':'cm s-1',
	    'cmap':'magma'
	},
        'v':{
	    'latname':'geolat_c',
	    'lonname':'geolon_c',
	    'plotname':'northward current',
	    'min':0.0,
	    'max':30.0,
	    'step':2.0,
	    'units':'cm s-1',
	    'cmap':'magma'
	}
    }
    #
    result = det_dict[vn][att]
    #
    return result


if __name__ == "__main__":
    ds = load_ds(sys.argv[1], sys.argv[2], sys.argv[3])
    main(ds, sys.argv[1], sys.argv[2], sys.argv[3])
