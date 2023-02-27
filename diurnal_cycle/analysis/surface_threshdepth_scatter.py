"""Plots scatter plots of surface diurnal range vs. threshold depth. 
  
One plot for temperature, one for current.


Usage:
    $ python surface_threshdepth_scatter.py plot_month
where:
    plot_month is an individual months (1-12)  
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

def main(plot_month):
    
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
    month_str = month_names[int(plot_month) - 1]
    
    # Construct file name and open it if it exists.
    dicy_fn1 = 'meanDiurnalCycle_metrics_TEMP_' + month_str + '.nc'
    dicy_fn2 = 'meanDiurnalCycle_metrics_CUR_' + month_str + '.nc'
    ds1 = xr.open_dataset(ddir + dicy_fn1, decode_timedelta=False)
    ds2 = xr.open_dataset(ddir + dicy_fn2, decode_timedelta=False)

    # Convert all velocities to cm s-1.
    for v in ['RANGE', 'RANGE_25']:
        ds2[v] = ds2[v]*100.0
        ds2[v].attrs['units'] = 'cm s-1'
    
    #----------------------------------------------------------------------------
    # Plot figure.
    
    fig_letters = 'abcdefghijklmnopqrstuvwxyz'
    
    # Define whether and where to mask fields according to diurnal range.
    range_mask = 0.0 #plot_dets(plot_var, 'mask_lim')
    
    # Set up figure and axes.
    fig, axs = plt.subplots(1, 2,
                            sharey=True, sharex=False,
                            figsize=(8, 5))     
    
    # Divide up data by latitude (for colors on scatter plot).
    trop_lim = 25.0
    pol_lim = 35.0
    ds1_t = ds1.where(np.abs(ds1.geolat_t) <= trop_lim, other=np.nan)
    ds1_m = ds1.where((np.abs(ds1.geolat_t) > trop_lim) & 
                      (np.abs(ds1.geolat_t) < pol_lim), other=np.nan)
    ds1_p = ds1.where(np.abs(ds1.geolat_t) >= pol_lim, other=np.nan)
    ds2_t = ds2.where(np.abs(ds2.geolat_c) <= trop_lim, other=np.nan)
    ds2_m = ds2.where((np.abs(ds2.geolat_c) > trop_lim) &
                      (np.abs(ds2.geolat_c) < pol_lim), other=np.nan)
    ds2_p = ds2.where(np.abs(ds2.geolat_c) >= pol_lim, other=np.nan)

    # Define details by dataset.
    ds_list = [ds1_t, ds1_m, ds1_p, ds2_t, ds2_p]
    col_list = ['red', 'black', 'blue', 'red', 'blue']
    size_list = [0.2, 0.1, 0.2, 0.2, 0.2]
    marker_list = ['o', '.', 's', 'o', 's']
    label_list = ['|lat| <= ' + str(trop_lim),
                  str(trop_lim) + ' < |lat| < ' + str(pol_lim),
                  '|lat| > ' + str(pol_lim),
                  '|lat| <= ' + str(trop_lim),
                  '|lat| > ' + str(pol_lim)]
    axno_list = [0, 0, 0, 1, 1]
    bin_width_list = [0.05, 0.05, 0.05, 1.0, 1.0]

    # Plot data.
    for ids, ds in enumerate(ds_list):
        ds_xcenters, ds_ymedians, ds_y5p, ds_y95p = \
            bin_stats(ds['RANGE'], ds['THRESH_DEPTH'], 
                      bin_width=bin_width_list[ids], 
                      bin_min_count=20)
        axs[axno_list[ids]].errorbar(
            ds_xcenters, ds_ymedians,
            yerr=[ds_ymedians - ds_y5p,
                  ds_ymedians + ds_y95p],
            fmt='None', ecolor=col_list[ids],
            label=label_list[ids]
        )
    """
    axs[0].scatter(ds1_t['RANGE'], 
                   ds1_t['THRESH_DEPTH'],
                   c='red',
                   s=0.2, marker='o', alpha=0.02,
                   label='|lat| <= ' + str(trop_lim))
    axs[0].scatter(ds1_m['RANGE'],
                   ds1_m['THRESH_DEPTH'],
                   c='black',
                   s=0.1, marker='.', alpha=0.02,
                   label=str(trop_lim) + ' < |lat| < ' + str(pol_lim))
    axs[0].scatter(ds1_p['RANGE'],
                   ds1_p['THRESH_DEPTH'],
                   c='blue',
                   s=0.2, marker='s', alpha=0.02,
                   label='|lat| > ' + str(pol_lim))
    axs[1].scatter(ds2_t['RANGE'],
                   ds2_t['THRESH_DEPTH'],
                   c='red',
                   s=0.2, marker='o', alpha=0.02,
                   label='|lat| <= ' + str(trop_lim))
    axs[1].scatter(ds2_p['RANGE'],
                   ds2_p['THRESH_DEPTH'],
                   c='blue',
                   s=0.2, marker='s', alpha=0.02,
                   label='|lat| > ' + str(pol_lim))"""
    
    # Set plot details.
    axs[0].set_xlabel('SST diurnal range (K)')
    axs[0].set_ylabel('threshold depth (m)')
    axs[1].set_xlabel('along-wind current diurnal range (' + r'$cm~s^{-1}$' + ')')
    axs[0].set_title(month_str, loc='left')

    # Add legend.
    axs[0].legend(loc='upper right', frameon=False)

    # Add figure letters (a, b, c, etc.).
    axs[0].text(0.1, 0.9, fig_letters[0],
                transform=axs[0].transAxes,
                horizontalalignment='center', verticalalignment='center',
                bbox={'facecolor':'white', 'alpha':0.9, 'pad':0.1,
                      'boxstyle':'round'})
    axs[1].text(0.1, 0.9, fig_letters[1],
                transform=axs[1].transAxes,
                horizontalalignment='center', verticalalignment='center',
                bbox={'facecolor':'white', 'alpha':0.9, 'pad':0.1,
                      'boxstyle':'round'})

    # Save file.
    plotdir = '/ncrc/home2/Jack.Reeveseyre/cfs-analysis/diurnal_cycle/plots/'
    plotfn = 'surface_treshdepth_scatter_' + month_str + '.'
    plotfileformat='png'
    plt.savefig(plotdir + plotfn + plotfileformat, format=plotfileformat,
                dpi=400)
    #
    return


def bin_stats(x, y, bin_width=1.0, bin_min_count=20):
    bin_edges = np.arange(0.0, np.ceil(np.nanmax(x)) + 0.0001,
                          bin_width)
    bin_c = 0.5*(bin_edges[:-1] + bin_edges[1:])
    bin_m = np.nan*np.zeros(len(bin_c))
    bin_5 = np.nan*np.zeros(len(bin_c))
    bin_95 = np.nan*np.zeros(len(bin_c))
    for i_b in range(len(bin_c)):
        quants = np.nanquantile(
            y.where((x >= bin_edges[i_b]) &
                    (x < bin_edges[i_b + 1])),
            np.array([0.05, 0.5, 0.95])
        )
        count = np.sum(~np.isnan(
            y.where((x >= bin_edges[i_b]) &
              (x < bin_edges[i_b + 1]))
        ))
        if count >= bin_min_count:
            bin_5[i_b] = quants[0]
            bin_m[i_b] = quants[1]
            bin_95[i_b] = quants[2]
    #
    return bin_c, bin_m, bin_5, bin_95


if __name__ == "__main__": 
    main(sys.argv[1])

