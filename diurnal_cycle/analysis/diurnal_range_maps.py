"""Plots maps of diurnal cycle metrics.

Usage:
    $ python diurnal_range_maps.py plot_month plot_var plot_type

where:
    plot_month is 'ANN' 
               or an individual month (1-12) 
               or a comma separated list of months (e.g., "1,2,12")
    plot_var is one of {'TEMP', 'CUR', 'WIND'}
    plot_type is one of {'RANGE', 'PHASE', 'HALF_DEPTH', 'PHASE_DELAY', 'ALL'}
        (Note that HALF_DEPTH and PHASE_DELAY aren't plotted for WIND.)
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

def main(plot_month, plot_var, plot_type):
    
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
    
    # Construct file name and open it if it exists.
    dicy_fn = 'meanDiurnalCycle_metrics_' + plot_var + '_' + month_str + '.nc'
    try:
        ds = xr.open_dataset(ddir + dicy_fn, decode_timedelta=False)
    except:
        # If it doesn't exist, construct it:
        
        # Load data.
        if plot_var == 'WIND':
            file_list = get_filelist('atmo', plot_month, ddir)
            ds = xr.open_mfdataset(file_list)[['UGRD', 'VGRD']].sel(height=10.0)
            # Take time average.
            ds = ds.mean(dim='time', keep_attrs=True)
            ds.attrs['input_files'] = file_list
        elif plot_var == 'CUR':
            surface_file_list = get_filelist('surface', plot_month, ddir)
            # Or maybe should be:
            #surface_file_list = get_filelist('surface', 'ANN', ddir)
            ds_stress_rec_mean = xr.open_mfdataset(surface_file_list)[['tau_x', 'tau_y']]\
            .mean(dim=['time','hour'], keep_attrs=True)
            file_list = get_filelist('ocn', plot_month, ddir)
            ds = xr.open_mfdataset(file_list, chunks=ocn_chunks)[['u', 'v']]
            # Calculate down-wind current.
            ds_stress_rec_mean = utils.wind_dir_twice(ds_stress_rec_mean, 'tau_x', 'tau_y')
            ds = utils.rotation_xy(ds, 'u', 'v',
                                   np.radians(ds_stress_rec_mean['arctan_VoverU'].data),
                                   'math', exp_dim=(0,1))
            # Take time average.
            ds = ds.mean(dim='time', keep_attrs=True)
            # Convert to shear, relative to a fixed depth.
            z = ds.st_ocean.data
            st_ocean_base = np.nanmax(z[z <= curr_shear_depth])
            xr.set_options(keep_attrs=True)
            for varname in ['u','v','U_streamwise','U_crossstream']:
                ds[varname] = ds[varname] \
                    - ds[varname].sel(st_ocean=st_ocean_base)
            xr.set_options(keep_attrs=False)
            ds.attrs['input_files'] = file_list
            ds_var = 'U_streamwise'
            thresh_depth = 0.02 # m s-1
        elif plot_var == 'TEMP':
            file_list = get_filelist('ocn', plot_month, ddir)
            ds = xr.open_mfdataset(file_list, chunks=ocn_chunks)[['temp']]
            # Take time average.
            ds = ds.mean(dim='time', keep_attrs=True)
            ds.attrs['input_files'] = file_list
            ds_var = 'temp'
            thresh_depth = 0.1 # K
        else:
            sys.exit('plot_var not recognized.')
         
        # Change NaN to an unrealistic value.
        # (Some of the functions below cannot handle all-NaN slices.)
        nan_hr_mask = xr.where(
            np.isnan(ds[ds_var].isel(st_ocean=0)).sum(dim='hour') < 24,
            1, np.nan 
        )
        unrealistic_value = -9999.0
        ds[ds_var] = ds[ds_var].where(~np.isnan(ds[ds_var]), 
                                      other=unrealistic_value)
        
        # Calculate statistics (range etc.)
        if plot_var in ['TEMP', 'CUR']:
            RANGE_3D = (ds[ds_var].max(dim='hour') 
                        - ds[ds_var].min(dim='hour')).compute()
            iPHASE = ds[ds_var].argmax(dim='hour').compute()
            local_m_utc = approx_local_time(ds)
            PHASE_3D = (ds.hour[iPHASE] + local_m_utc) % 24
            iHALF = (RANGE_3D >= 0.5*RANGE_3D.isel(st_ocean=0))\
                .cumprod(dim='st_ocean').sum(dim='st_ocean') - 1
            iHALF = iHALF.where(iHALF >= 0, other=0)
            iTHRESH = (RANGE_3D >= thresh_depth)\
                .cumprod(dim='st_ocean').sum(dim='st_ocean') - 1
            iTHRESH = iTHRESH.where(iTHRESH >= 0, other=0)
            ds['RANGE'] = RANGE_3D.isel(st_ocean=0)*nan_hr_mask
            ds['PHASE'] = PHASE_3D.isel(st_ocean=0)*nan_hr_mask
            ds['RANGE_25'] = RANGE_3D.sel(st_ocean=25.0)*nan_hr_mask
            ds['PHASE_25'] = PHASE_3D.sel(st_ocean=25.0)*nan_hr_mask
            ds['HALF_DEPTH'] = ds.st_ocean[iHALF]*nan_hr_mask
            ds['PHASE_DELAY'] = (PHASE_3D.isel(st_ocean=iHALF) 
                                 - PHASE_3D.isel(st_ocean=0))*nan_hr_mask % 24
            ds['THRESH_DEPTH'] = ds.st_ocean[iTHRESH]*nan_hr_mask
            ds['THRESH_DELAY'] = (PHASE_3D.isel(st_ocean=iTHRESH)
                                  - PHASE_3D.isel(st_ocean=0))*nan_hr_mask % 24
            ds_out = ds[['RANGE', 'PHASE', 'HALF_DEPTH', 'PHASE_DELAY',
                         'RANGE_25', 'PHASE_25', 'THRESH_DEPTH', 'THRESH_DELAY']]
        else:
            ds['RANGE'] = (ds[ds_var].max(dim='hour') 
                           - ds[ds_var].min(dim='hour')).compute()*nan_hr_mask
            local_m_utc = approx_local_time(ds)
            iPHASE = ds[ds_var].argmax(dim='hour').compute()
            ds['PHASE'] = (ds.hour[iPHASE] + local_m_utc) % 24
            ds_out = ds[['RANGE', 'PHASE']]
        
        # Save out data.
        ds_out['RANGE'].attrs = {
            'long_name':'range of mean diurnal cycle',
            'units':ds[ds_var].attrs['units'],
            'variable':ds[ds_var].attrs['long_name'],
            'cell_methods':month_str + ' mean over years & months for each hour of day; (max - min) over hours of day.',
            'time_period':month_str + ' 2002-2005' 
        }
        ds_out['PHASE'].attrs = {
            'long_name':'local time of diurnal cycle maximum',
            'units':'hours',
            'variable':ds[ds_var].attrs['long_name'],
            'cell_methods':month_str + ' mean over years & months for each hour of day; local time of max over hours of day.',
            'time_period':month_str + ' 2002-2005'
        }
        if plot_var in ['TEMP', 'CUR']:
            ds_out['HALF_DEPTH'].attrs = {
                'long_name':'depth at which range of mean diurnal cycle is half of surface range',
                'units':'m',
                'variable':ds[ds_var].attrs['long_name'],
                'cell_methods':month_str + ' mean over years & months for each hour of day; (max - min) over hours of day; deepest grid level i such that range[i] >= range[surface].',
                'time_period':month_str + ' 2002-2005'
            }
            ds_out['PHASE_DELAY'].attrs = {
                'long_name':'delay of time of diurnal cycle maximum at HALF_DEPTH after maximum at surface',
                'units':'hours',
                'variable':ds[ds_var].attrs['long_name'],
                'cell_methods':month_str + ' mean over years & months for each hour of day ; (max - min) over hours of day.',
                'time_period':month_str + ' 2002-2005'
            }
            ds_out['RANGE_25'].attrs = {
                'long_name':'range of mean diurnal cycle at 25 m depth',
                'units':ds[ds_var].attrs['units'],
                'variable':ds[ds_var].attrs['long_name'],
                'cell_methods':month_str +  ' mean over years & months for each hour of day; (max - min) over hours of day.',
                'time_period':month_str + ' 2002-2005'
            }
            ds_out['PHASE_25'].attrs = {
                'long_name':'local time of diurnal cycle maximum at 25 m depth',
                'units':'hours',
                'variable':ds[ds_var].attrs['long_name'],
                'cell_methods':month_str + ' mean over years & months for each hour of day ; (max - min) over hours of day.',
                'time_period':month_str + ' 2002-2005'
            }
            ds_out['THRESH_DEPTH'].attrs = {
                'long_name':'deepest level at which range of mean diurnal cycle is >= ' + str(thresh_depth),
                'units':'m',
                'threshold':str(thresh_depth) + ' ' + ds[ds_var].attrs['units'],
                'variable':ds[ds_var].attrs['long_name'],
                'cell_methods':month_str +  ' mean over years & months for each hour of day; (max - min) over hours of day; deepest grid level i such that range[i] >=' + str(thresh_depth) + ' ' + ds[ds_var].attrs['units'],
                'time_period':month_str + ' 2002-2005'
            }
            ds_out['THRESH_DELAY'].attrs = {
                'long_name':'delay of time of diurnal cycle maximum at THRESH_DEPTH after maximum at surface',
                'units':'hours',
                'threshold':str(thresh_depth) + ' ' + ds[ds_var].attrs['units'],
                'variable':ds[ds_var].attrs['long_name'],
                'cell_methods':month_str + ' mean over years & months for each hour of day ; (max - min) over hours of day.',
                'time_period':month_str + ' 2002-2005'
            }
        fill_val_f64 = 9.969209968386869e+36
        nc_enc = {}
        for k in ds_out.keys():
            nc_enc[k] = {'dtype':'float64', '_FillValue':fill_val_f64}
        print(ds_out)
        ds_out.to_netcdf(path=ddir + dicy_fn,
                         mode='w',format='NETCDF4',
                         encoding=nc_enc)

    # Add PHASE in UTC.
    #print('WARNING: converting PHASE to UTC.')
    #local_m_utc = approx_local_time(ds)
    #ds['PHASE'] = (ds['PHASE'] - local_m_utc) % 24
    
    # Convert all velocities to cm s-1.
    if plot_var == 'CUR':
        for v in ['RANGE', 'RANGE_25']:
            ds[v] = ds[v]*100.0
            ds[v].attrs['units'] = 'cm s-1'
    
    # Plot.
    if plot_type == 'ALL':
        plot_all_maps(ds, plot_var, month_str)
    else:
        plot_one_map(ds, plot_type, plot_var, month_str)
    return


def plot_one_map(ds, ptype, pvar, mon):
    
    # Define whether and where to mask fields according to diurnal range.
    range_mask = 0.0 #plot_dets(pvar, 'mask_lim')
    
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
    cmap_b = plt.get_cmap(plot_dets(pvar, 'cmap', ptype))
    norm_b = colors.BoundaryNorm(np.arange(plot_dets(pvar, 'min', ptype),
                                           plot_dets(pvar, 'max', ptype)*1.01,
                                           plot_dets(pvar, 'step', ptype)),
                                 ncolors=cmap_b.N)
    
    # Plot the data.
    if ((ptype == 'RANGE') or (range_mask == 0.0)):
        p = ax.pcolormesh(switch_lon_lims(ds[plot_dets(pvar, 'lonname')].data, -180.0),
                          ds[plot_dets(pvar, 'latname')].data,
	        	  ds[ptype].data,
                          transform=data_crs,
                          cmap=cmap_b, norm=norm_b, shading='nearest')
    else:
        import numpy.ma as ma
        dsm = ma.masked_where(ds['RANGE'].data < range_mask, ds[ptype].data)
        ptr = ax.pcolormesh(switch_lon_lims(ds[plot_dets(pvar, 'lonname')].data, -180.0),
                            ds[plot_dets(pvar, 'latname')].data,
	        	    ds[ptype].data, alpha=0.25,
                            transform=data_crs,
                            cmap=cmap_b, norm=norm_b, shading='nearest')
        p = ax.pcolormesh(switch_lon_lims(ds[plot_dets(pvar, 'lonname')].data, -180.0),
                          ds[plot_dets(pvar, 'latname')].data,
	        	  dsm,
                          transform=data_crs,
                          cmap=cmap_b, norm=norm_b, shading='nearest')
    
    # Finish up figure.
    if ptype in ['THRESH_DEPTH', 'THRESH_DELAY']:
        ax.set_title(plot_dets(pvar, 'plotname') + ' DIURNAL CYCLE ' + ptype 
                     + ', ' + mon + '\nthreshold: ' + ds[ptype].attrs['threshold'])
    else:
        ax.set_title(plot_dets(pvar, 'plotname') + ' DIURNAL CYCLE ' + ptype
                     + ', ' + mon)
    fig.colorbar(p,
                 label=ptype + ' (' + plot_dets(pvar, 'units', ptype) + ')',
                 ax=ax, orientation='horizontal',shrink=0.5)
    
    # Save file.
    plotdir = '/ncrc/home2/Jack.Reeveseyre/cfs-analysis/diurnal_cycle/plots/'
    plotfn = 'map_' + pvar + '_' + ptype + '_' + mon + '.'
    plotfileformat='png'
    plt.savefig(plotdir + plotfn + plotfileformat, format=plotfileformat,
                dpi=400)
    #
    return


def plot_all_maps(ds, var, mon):

    # Set up figure and axes.
    fig = plt.figure(figsize=(8, 6))
    #
    return


def switch_lon_lims(lon_list, min_lon=0.0):
    result = (lon_list - min_lon) % 360.0 + min_lon
    return result


def approx_local_time(ds):
    if 'geolon_t' in ds.keys():
        lon = ds['geolon_t']
    elif 'geolon_c' in ds.keys():
        lon = ds['geolon_c']
    else:
        print('longitude variable not found.')
        print(ds)
        lon_name = input('Please enter the name of the longitude varible:\n')
        lon = ds[lon_name]
    lt = switch_lon_lims(lon, min_lon=-180.0)/15.0
    return lt


def get_filelist(comp, plot_month, ddir):
    if comp == 'atmo':
        comp_str = 'atmo'
        end_str = ')_grid_meanDiurnalCycle.nc'
    elif comp == 'ocn':
        comp_str = 'ocn'
        end_str = ')_grid_meanDiurnalCycle.nc'
    elif comp == 'surface':
        comp_str = 'ocn'
        end_str = ')_grid_surface_meanDiurnalCycle.nc'
    else:
        sys.exit('model component not recognized, comp must be {atmo, surface, ocn}')

    # Construct regex.
    if plot_month == 'ANN':
        pattern = '^' + comp_str + '.*_grid_meanDiurnalCycle.nc'
    else:
        plot_month_list = np.array(plot_month.split(','), dtype=int)
        pattern = '^' + comp_str + '_.*_(' + \
            '|'.join("{:02d}".format(m) for m in plot_month_list) + \
            end_str

    # Get list of matching files.
    result = []
    allfiles = os.listdir(ddir)
    allfiles.sort()
    for fn in allfiles:
        if re.match(pattern, fn):
            result.append(ddir + fn)
    #
    return result


def plot_dets(vn, att, pty=None):
    det_dict = {
        'TEMP':{
            'mask_lim':0.1,
	    'latname':'geolat_t',
	    'lonname':'geolon_t',
	    'plotname':'TEMPERATURE',
	    'minRANGE':0.0,
	    'maxRANGE':1.0,
            'stepRANGE':0.1,
	    'unitsRANGE':'K',
	    'cmapRANGE':'magma',
            'extendRANGE':'max',
	    'minRANGE_25':0.0,
            'maxRANGE_25':0.5,
            'stepRANGE_25':0.05,
            'unitsRANGE_25':'K',
            'cmapRANGE_25':'magma',
            'extendRANGE_25':'max',
            'minPHASE':0.0,
	    'maxPHASE':24.0,
	    'stepPHASE':1.0,
	    'cmapPHASE':'twilight',
	    'unitsPHASE':'local time',
            'extendPHASE':'neither',
            'minPHASE_25':0.0,
            'maxPHASE_25':24.0,
            'stepPHASE_25':1.0,
            'cmapPHASE_25':'twilight',
            'unitsPHASE_25':'local time',
            'extendPHASE':'neither',
	    'minHALF_DEPTH':0.0,
	    'maxHALF_DEPTH':20.0,
	    'stepHALF_DEPTH':2.0,
	    'cmapHALF_DEPTH':'magma',
	    'unitsHALF_DEPTH':'m',
            'extendHALF_DEPTH':'max',
	    'minPHASE_DELAY':0.0,
	    'maxPHASE_DELAY':6.0,
	    'stepPHASE_DELAY':0.5,
	    'cmapPHASE_DELAY':'magma',
	    'unitsPHASE_DELAY':'hours',
            'extendPHASE_DELAY':'neither',
            'minTHRESH_DEPTH':0.0,
            'maxTHRESH_DEPTH':20.0,
            'stepTHRESH_DEPTH':2.0,
            'cmapTHRESH_DEPTH':'magma',
            'unitsTHRESH_DEPTH':'m',
            'extendTHRESH_DEPTH':'max',
            'minTHRESH_DELAY':0.0,
            'maxTHRESH_DELAY':6.0,
            'stepTHRESH_DELAY':0.5,
            'cmapTHRESH_DELAY':'magma',
            'unitsTHRESH_DELAY':'hours',
            'extendTHRESH_DELAY':'max'

	},
        'CUR':{
            'mask_lim':2.0,
            'latname':'geolat_c',
            'lonname':'geolon_c',
            'plotname':'ALONG-WIND CURRENT',
            'minRANGE':0.0,
            'maxRANGE':30.0,
            'stepRANGE':2.0,
            'unitsRANGE':'cm s-1',
            'cmapRANGE':'magma',
            'extendRANGE':'max',
            'minRANGE_25':0.0,
            'maxRANGE_25':0.2,
            'stepRANGE_25':0.01,
            'unitsRANGE_25':'cm s-1',
            'cmapRANGE_25':'magma',
            'extendRANGE_25':'max',
            'minPHASE':0.0,
            'maxPHASE':24.0,
            'stepPHASE':1.0,
            'cmapPHASE':'twilight',
            'unitsPHASE':'local time',
            'extendPHASE':'neither',
            'minPHASE_25':0.0,
            'maxPHASE_25':24.0,
            'stepPHASE_25':1.0,
            'cmapPHASE_25':'twilight',
            'unitsPHASE_25':'local time',
            'extendPHASE_25':'neither',
            'minHALF_DEPTH':0.0,
            'maxHALF_DEPTH':20.0,
            'stepHALF_DEPTH':2.0,
            'cmapHALF_DEPTH':'magma',
            'unitsHALF_DEPTH':'m',
            'extendHALF_DEPTH':'max',
            'minPHASE_DELAY':0.0,
            'maxPHASE_DELAY':6.0,
            'stepPHASE_DELAY':0.5,
            'cmapPHASE_DELAY':'magma',
            'unitsPHASE_DELAY':'hours',
            'extendPHASE_DELAY':'max',
            'minTHRESH_DEPTH':0.0,
            'maxTHRESH_DEPTH':40.0,
            'stepTHRESH_DEPTH':4.0,
            'cmapTHRESH_DEPTH':'magma',
            'unitsTHRESH_DEPTH':'m',
            'extendTHRESH_DEPTH':'max',
            'minTHRESH_DELAY':0.0,
            'maxTHRESH_DELAY':24.0,
            'stepTHRESH_DELAY':2.0,
            'cmapTHRESH_DELAY':'magma',
            'unitsTHRESH_DELAY':'hours',
            'extendTHRESH_DELAY':'neither'

        }
    }
    #
    if att in ['min', 'max', 'cmap', 'step', 'units', 'extend']:
        result = det_dict[vn][att+pty]
    else:
        result = det_dict[vn][att]
    #
    return result


if __name__ == "__main__": 
    main(sys.argv[1], sys.argv[2], sys.argv[3])
