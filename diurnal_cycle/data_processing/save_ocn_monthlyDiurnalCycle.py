import sys
import xarray as xr
import numpy as np
import pandas as pd
import datetime

def main(year, month, xy_method):

    print('Processing ' + str(year) + '-' + str(month) + '...')
    
    # Spatial details.
    z_lim = 100.0
    # Lats/lons only used if xy_method == 'points':
    interp_method = 'nearest'
    lats = np.array([-8.0, -5.0, -2.0, 0.0, 2.0, 5.0, 8.0, 9.0])
    lons = np.array([137.0, 147.0, 156.0, 165.0, 180.0, 190.0, 
                     205.0, 220.0, 235.0, 250.0, 265.0])
    if np.max(lats) > 60.0:
        print('Warning: interpolation to points not currently supported north of 60N.')
    lons = switch_lon_lims(lons, min_lon=-280.0)
    
    if xy_method == 'points':
        t1 = datetime.datetime.now()
        ds = get_mon_di_cy_points(year, month, 
                                  lats, lons, z_lim,
                                  interp_method)
        t2 = datetime.datetime.now()
        dt2 = t2 - t1
        print(f'Time taken to read with open_mfdataset: {dt2}')
    elif xy_method == 'grid':
        t1 = datetime.datetime.now()
        ds = get_mon_di_cy_grid(year, month, z_lim)
        t2 = datetime.datetime.now()
        dt2 = t2 - t1
        print(f'Time taken to read with open_mfdataset: {dt2}')
    else:
        sys.exit('xy_method should be one of {points, grid}')
    
    # Define encoding.
    fill_val_f64 = 9.969209968386869e+36
    fill_val_i32 = np.iinfo(np.int32).min + 2
    fill_val_i64 = np.iinfo(np.int64).min + 2
    t_units = 'seconds since 1970-01-01 00:00:00'
    t_cal = 'standard'
    encdg = {}
    for k in ds.data_vars.keys():
        encdg[k] = {'_FillValue':fill_val_f64, 'dtype':'float64'}
    for k in ds.coords.keys():
        encdg[k] = {'_FillValue':fill_val_f64, 'dtype':'float64'}
    encdg['time']['units'] = t_units
    encdg['time']['calendar'] = t_cal

    # Save file.
    out_dir = '/lustre/f2/dev/ncep/Jack.Reeveseyre/diurnal_cycle/'
    out_fn = 'ocn_' + "{:0>4d}".format(year) + '_' + "{:0>2d}".format(month) \
        + '_' + xy_method + '_meanDiurnalCycle.nc'
    ds.to_netcdf(path=out_dir + out_fn,
                 mode='w',format='NETCDF4',
                 encoding=encdg,unlimited_dims=['time'])
    #
    return 


def get_mon_di_cy_points(yr, mn, lat_list, lon_list, max_depth,
                         select_how, interp_how):

    var_list = ['temp', 'salt', 'u', 'v']
    data_dir = '/lustre/f2/dev/ncep/JieShun.Zhu/CFSm501hr/CFSm501hr1980010200/DATA/'
    fn_root = 'ocn_' + "{:0>4d}".format(yr) +'_' + "{:0>2d}".format(mn)
    
    # Loads all files for one month. 
    ds_all = xr.open_mfdataset(data_dir + fn_root + '*.nc')  

    # Interpolated to points and calculate mean diurnal cycle.
    ds_sub = ds_all[var_list].interp(
        xt_ocean=lon_list, yt_ocean=lat_list,
        xu_ocean=lon_list, yu_ocean=lat_list,
        method=interp_how
    )
    ds_sub = ds_sub.sel(st_ocean=slice(0,max_depth))
    di_cy = ds_sub.groupby('time.hour')\
            .mean('time').chunk({'hour':24}).load()
    di_cy = di_cy.assign_coords({'time':np.datetime64(str(yr) + '-' + "{:0>2d}".format(mn), 's')})
    di_cy = di_cy.expand_dims(dim='time', axis=0)

    # Copy/add attributes.
    for v in var_list:
        di_cy[v].attrs = ds_sub[v].attrs
        di_cy[v].attrs['cell_methods'] = \
            'hour: mean within hours, time: mean over days (comment: monthly mean diurnal cycle), [xy][ut]_ocean: ' + interp_how
    di_cy.attrs = ds_sub.attrs
    di_cy.attrs['postprocessing'] = 'monthly mean diurnal cycle, spatially interpolated to points using ' + interp_how 

    ds_all.close()
    
    return di_cy


def get_mon_di_cy_grid(yr, mn, max_depth):

    var_list = ['temp', 'salt', 'u', 'v']
    data_dir = '/lustre/f2/dev/ncep/JieShun.Zhu/CFSm501hr/CFSm501hr1980010200/DATA/'
    fn_root = 'ocn_' + "{:0>4d}".format(yr) +'_' + "{:0>2d}".format(mn)

    # Loads all files for one month. 
    ds_all = xr.open_mfdataset(data_dir + fn_root + '*.nc')

    # Calculate mean diurnal cycle.
    ds_sub = ds_all[var_list]
    ds_sub = ds_sub.sel(st_ocean=slice(0,max_depth))
    di_cy = ds_sub.groupby('time.hour')\
            .mean('time').chunk({'hour':24}).load()
    di_cy = di_cy.assign_coords({'time':np.datetime64(str(yr) + '-' + "{:0>2d}".format(mn), 's')})
    di_cy = di_cy.expand_dims(dim='time', axis=0)

    # Copy/add attributes.
    for v in var_list:
        di_cy[v].attrs = ds_sub[v].attrs
        di_cy[v].attrs['cell_methods'] = \
            'hour: mean within hours, time: mean over days (comment: monthly mean diurnal cycle)'
    di_cy.attrs = ds_sub.attrs
    di_cy.attrs['postprocessing'] = 'monthly mean diurnal cycle'

    ds_all.close()

    return di_cy


def get_mon_di_cy_manual(yr, mn, lat_list, lon_list, max_depth,
                         select_how, interp_how):
    
    # Set up.
    var_list = ['temp', 'salt', 'u', 'v']
    data_dir = '/lustre/f2/dev/ncep/JieShun.Zhu/CFSm501hr/CFSm501hr1980010200/DATA/'
    time_array = np.arange(
        np.datetime64(str(yr) + '-' + "{:0>2d}".format(mn)), 
        np.datetime64(str(yr) + '-' + "{:0>2d}".format(mn)) + np.timedelta64(1,'M'), 
        dtype='datetime64[h]'
    )
    pddt = pd.DatetimeIndex(time_array)
    
    # Load data.
    ds_sub = xr.concat(
        [xr.open_dataset(data_dir + 'ocn_' + 
                         "{:0>4d}".format(tt.year) + '_' +
                         "{:0>2d}".format(tt.month) + '_' +
                         "{:0>2d}".format(tt.day) + '_' +
                         "{:0>2d}".format(tt.hour) + '.nc')[var_list]\
             .sel(st_ocean=slice(0,max_depth))\
             .interp(xt_ocean=lon_list, yt_ocean=lat_list,
                     xu_ocean=lon_list, yu_ocean=lat_list,
                     method=interp_how) 
         for tt in pddt],
        dim='time'
    )
    di_cy = ds_sub.groupby('time.hour').mean('time').chunk({'hour':24}).load() 
    di_cy = di_cy.assign_coords({'time':np.datetime64(str(yr) + '-' + "{:0>2d}".format(mn), 's')})
    di_cy = di_cy.expand_dims(dim='time', axis=0)

    # Copy/add attributes.
    for v in var_list:
        di_cy[v].attrs = ds_sub[v].attrs
        di_cy[v].attrs['cell_methods'] = \
            'hour: mean within hours, time: mean over days (comment: monthly mean diurnal cycle), [xy][ut]_ocean: ' + interp_how
    di_cy.attrs = ds_sub.attrs
    di_cy.attrs['postprocessing'] = 'monthly mean diurnal cycle, spatially interpolated to points using ' + interp_how
 
    ds_sub.close()
    
    return di_cy
    

def switch_lon_lims(lon_list, min_lon=0.0):
    result = (lon_list - min_lon) % 360.0 + min_lon
    return result


if __name__ == '__main__':
    cla = sys.argv
    xy = cla[1]
    for yr in np.arange(2002, 2006):
        for mn in np.arange(1,13):
            main(yr, mn, xy)
