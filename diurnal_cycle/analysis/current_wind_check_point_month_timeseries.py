"""Plots maps of current and wind components.
Usage:
    $ python current_wind_check_point_month_timeseries.py plot_date
where:
    plot_date is formatted like yyyy_mm
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

def load_ds(plot_date):

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
    var_list = ['u', 'v', 'temp', 'tau_x', 'tau_y']
    ds = xr.open_mfdataset(data_dir + fn_root + '*.nc',
                           chunks=ocn_chunks)[var_list]\
           .sel(st_ocean=0.5, 
                yu_ocean=[-35.0, -30.0,-25.0], 
                yt_ocean=[-35.0, -30.0,-25.0], 
                xu_ocean=-20.0,
                xt_ocean=-20.0, 
                method='nearest').load()
    ds['time'] = ds.indexes['time'].to_datetimeindex()

    return ds


def main(ds_in, plot_date):
    
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
    
    # Process data.
    ds = ds_in.copy()
    for plot_var in ['u', 'v']: 
        ds[plot_var] = ds[plot_var]*100.0
        ds[plot_var].attrs['units'] = 'cm s-1'
        ds[plot_var].attrs['long_name'] = {'u':'eastward current',
                                           'v':'northward current'}[plot_var]
    ds['tau_mag'] = np.sqrt(ds['tau_x']**2 + ds['tau_y']**2)
    ds['tau_mag'].attrs['units'] = ds['tau_x'].attrs['units']
    ds['tau_mag'].attrs['long_name'] = 'wind stress magnitude' 
    
    # Set up figure and axes.
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(6, 10), sharex=True)
    ax4.set_xlabel('time')
    ax4.set_xticks(np.arange(ds.time.min().data,
                             ds.time.max().data + np.timedelta64(2,'h'),
                             np.timedelta64(10, 'D')),
                   minor=False)
    ax4.set_xticks(np.arange(ds.time.min().data,
                             ds.time.max().data + np.timedelta64(2,'h'),
                             np.timedelta64(1, 'D')),
                   minor=True)

    # Define color map.
    cmap_b = cm.get_cmap('tab10')
    
    # Define axis labels for each plot.
    ax1.set_ylabel('u current velocity (cm s-1)')
    ax2.set_ylabel('v current velocity (cm s-1)')
    ax3.set_ylabel('SST (deg C)')
    ax4.set_ylabel(r'$|\tau|$' + ' (N m-2)')

    # Plot data.
    for ilat, lat in enumerate(ds.yu_ocean.data):
        if np.all(np.isnan(ds['u'].sel(yu_ocean=lat).data)):
            continue
        else:
            ax1.plot(ds.time.data,
                     ds['u'].sel(yu_ocean=lat).data,
                     color=cmap_b(ilat), linewidth=1.0,
                     label="{:+.2f}".format(ds.geolat_c.sel(yu_ocean=lat).data))
            ax2.plot(ds.time.data,
                     ds['v'].sel(yu_ocean=lat).data,
                     color=cmap_b(ilat), linewidth=1.0,
                     label="{:+.2f}".format(ds.geolat_c.sel(yu_ocean=lat).data))
            ax3.plot(ds.time.data,
                     ds['temp'].sel(yt_ocean=lat, method='nearest').data,
                     color=cmap_b(ilat), linewidth=1.0,
                     label="{:+.2f}".format(ds.geolat_t\
                                            .sel(yt_ocean=lat, method='nearest')\
                                            .data))
            ax4.plot(ds.time.data,
                     ds['tau_mag'].sel(yu_ocean=lat).data,
                     color=cmap_b(ilat), linewidth=1.0,
                     label="{:+.2f}".format(ds.geolat_c.sel(yu_ocean=lat).data))

    # Finish up figure.
    ax1.set_title(plot_date, loc='left')
    ax1.set_title('longitude: ' + 
                  str(ds.geolon_c.mean().data), 
                  loc='right')
    ax1.legend(ncol=3, loc='lower left') #, bbox_to_anchor=(1.0, 0.5))
    
    # Save file.
    plotdir = '/ncrc/home2/Jack.Reeveseyre/cfs-analysis/diurnal_cycle/plots/'
    plotfn = 'current_month_timeseries_' + plot_date + '.'
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
    ds = load_ds(sys.argv[1])
    main(ds, sys.argv[1])

    

