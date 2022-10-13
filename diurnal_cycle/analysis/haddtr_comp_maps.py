"""Plots maps of SST diurnal range and compares with observations.
Usage:
    $ python haddtr_comp_maps.py plot_month
where:
    plot_month is 'ANN' 
               or an individual month (JAN, FEB, etc.) 
               or a season ("DJF", etc.)
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
#------------------------------------------------------------------------------

def main(plot_month):


    # Load data:
    ddir = '/lustre/f2/dev/ncep/Jack.Reeveseyre/diurnal_cycle/'
    odir = '/ncrc/home2/Jack.Reeveseyre/data/haddtr/'
    # interpolated model
    dsm = xr.open_dataset(ddir + 'regrid/' +
                          'meanDiurnalCycle_metrics_TEMP_' +
                          plot_month + '_5x5_bilin.nc')\
          .rename(name_dict={'yt_ocean':'lat', 'xt_ocean':'lon'})
    # raw model
    dsr = xr.open_dataset(ddir +
                          'meanDiurnalCycle_metrics_TEMP_' +
                          plot_month + '.nc')
    # observations
    if plot_month == 'ANN':
        obs_time_period = 'annual'
    elif plot_month in ['DJF', 'MAM', 'JJA', 'SON']:
        obs_time_period = 'seasonal'
    else:
        obs_time_period = 'monthly'
    dso = xr.open_dataset(odir + obs_time_period + '_climatology.nc')

    # Extract and calculate required fields.
    time_name_dict = {
        'annual':{'ANN':'annual'},
        'seasonal':{'DJF':'Winter (DJF)', 'MAM':'Spring (MAM)',
                    'JJA':'Summer (JJA)', 'SON':'Autumn (SON)'},
        'monthly':{'JAN':'1', 'FEB':'2', 'MAR':'3',
                   'APR':'4', 'MAY':'5', 'JUN':'6',
                   'JUL':'7', 'AUG':'8', 'SEP':'9',
                   'OCT':'10', 'NOV':'11', 'DEC':'12'}
    }
    dtro = dso['dtr'].isel(time_index=(
        dso['time_name'] ==time_name_dict[obs_time_period][plot_month]
    )).squeeze()
    dtr_diff = dsm['RANGE'] - dtro.data
    plot_data = [
        {'data':dsm['RANGE'],
         'title':'CFSm501 interpolated',
         'colb':'abs',
         'grid':'had'},
        {'data':dsr['RANGE'],
         'title':'CFSm501 model grid',
         'colb':'abs',
         'grid':'mom'},
        {'data':dtro,
         'title':'HadDTR observations',
         'colb':'abs',
         'grid':'had'},
        {'data':dtr_diff,
         'title':'CFSm501 minus HadDTR',
         'colb':'diff',
         'grid':'had'}
    ]
    import pdb; pdb.set_trace()
    #-----------
    # Plot maps.
    
    # Set up figure and axes.
    data_crs = ccrs.PlateCarree()
    plot_crs = ccrs.Robinson(central_longitude=-150)
    fig = plt.figure(figsize=(9, 6))
    axes_class = (GeoAxes,
                  dict(map_projection=plot_crs))
    axgr = AxesGrid(fig, 111, axes_class=axes_class,
                    nrows_ncols=(2, 2),
                    axes_pad=(0.25,1.0),
                    cbar_location='bottom',
                    cbar_mode='edge',
                    cbar_pad=0.25,
                    cbar_size='2%',
                    label_mode='')
    
    # Define color maps.
    cmap_b = plt.get_cmap('magma')
    norm_b = colors.BoundaryNorm(np.arange(0.0, 1.01, 0.1),
                                 ncolors=cmap_b.N)
    cmap_d = plt.get_cmap('BrBG')
    norm_d = colors.BoundaryNorm(np.arange(-0.5, 0.51, 0.1),
                                 ncolors=cmap_d.N)
    
    # Loop over axes and plot each.
    plots=  {}
    for i, ax in enumerate(axgr):
        print(i)
        plot_dict = plot_data[i]
        if plot_dict['grid'] == 'had':
            latm, lonm = np.meshgrid(plot_dict['data'].coords['lat'],
                                     plot_dict['data'].coords['lon'],
                                     indexing='ij')
        else:
            latm = dsr.geolat_t.data
            lonm = switch_lon_lims(dsr.geolon_t.data, -180.0)
         
        # Add features.
        ax.set_global()
        gl = ax.gridlines(draw_labels=True,
                          ylabel_style={'rotation':45.0},
                          xlocs=np.arange(-300, 151, 60),
                          ylocs=np.arange(-90, 91, 30),
                          linewidths=0.3)
        gl.xformatter = LongitudeFormatter(zero_direction_label=False,
                                           degree_symbol='')
        gl.yformatter = LatitudeFormatter(degree_symbol='')
        ax.coastlines(zorder=5)
        ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=4)
        
        # Plot the data.
        plots[plot_dict['colb']] = ax.pcolormesh(
            lonm, latm,
	    plot_dict['data'].data,
            transform=data_crs,
            cmap=cmap_d if plot_dict['colb'] == 'diff' else cmap_b,
            norm=norm_d if plot_dict['colb'] == 'diff' else norm_b,
            shading='nearest')
        
        # Add title.
        ax.set_title(plot_dict['title'])
        
        # Sort ticklabels.
        if i in [0, 1]:
           gl.xlabels_top=True
           gl.xlabels_bottom=False
        else:
           gl.xlabels_top=False 
           gl.xlabels_bottom=True
        if i in [0, 2]:
           gl.ylabels_left=True
           gl.ylabels_right=False
        else:
           gl.ylabels_left=False
           gl.ylabels_right=True
    
    # Add colorbars.
    cba = axgr.cbar_axes[0].colorbar(plots['abs'], 
                                     pad=0.35, shrink=0.6, aspect=10,
                                     extend='max',
                                     label='SST diurnal range (K)')
    cbd = axgr.cbar_axes[1].colorbar(plots['diff'], 
                                     pad=0.35, shrink=0.6, aspect=10,
                                     extend='both',
                                     label='SST diurnal range difference (K)')
    
    # Add over all title.
    fig.suptitle('SST diurnal range, ' + plot_month)
    
    # Save file.
    plotdir = '/ncrc/home2/Jack.Reeveseyre/cfs-analysis/diurnal_cycle/plots/'
    plotfn = 'map_haddtr_comp_' + plot_month + '.'
    plotfileformat='png'
    plt.savefig(plotdir + plotfn + plotfileformat, format=plotfileformat,
                dpi=400)
    #
    return


def switch_lon_lims(lon_list, min_lon=0.0):
    result = (lon_list - min_lon) % 360.0 + min_lon
    return result


#------------------------------------------------------------------------------
if __name__ == "__main__":
    main(sys.argv[1])
