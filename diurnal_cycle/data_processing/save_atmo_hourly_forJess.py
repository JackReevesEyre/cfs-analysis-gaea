import xarray as xr
import numpy as np
import pandas as pd
import datetime
import sys
import re

def main(xy_method, h_start, h_end, indir, outdir, f_list):
    
    ref_time = np.datetime64('1980-01-02T00:00')
    t_start = ref_time + np.timedelta64(h_start, 'h')
    t_end = ref_time + np.timedelta64(h_end, 'h')
    t_array_test = np.arange(t_start, t_end + np.timedelta64(1,'h'), 
                             dtype='datetime64[h]')
    ys = pd.to_datetime(str(t_start)).strftime('%Y')
    ye = pd.to_datetime(str(t_end)).strftime('%Y')
    ms = pd.to_datetime(str(t_start)).strftime('%m')
    me = pd.to_datetime(str(t_end)).strftime('%m')
    if ((ys != ye) | (ms != me)):
        print('Specified file list contains data from more than one month!')
    else:
        year = int(ys)
        month = int(ms)
    
    print('Processing ' + str(year) + '-' + str(month) + '...')
   
    # Spatial details:
    z_lim = 100.0
    # Only needed if xy_method == 'points':
    interp_method = 'nearest' 
    lats = np.array([-8.0, -5.0, -2.0, 0.0, 2.0, 5.0, 8.0, 9.0])
    lons = np.array([137.0, 147.0, 156.0, 165.0, 180.0, 190.0, 
                     205.0, 220.0, 235.0, 250.0, 265.0])
    if np.max(lats) > 60.0:
        print('Warning: interpolation to points not currently supported north of 60N.') 
    
    if xy_method == 'points':
        t1 = datetime.datetime.now()
        ds = get_hourly_points(year, month, 
                               indir, f_list, 
                               lats, lons, z_lim,
                               interp_method,
                               t_array_test)
        t2 = datetime.datetime.now()
        dt2 = t2 - t1
        print(f'Time taken to read with open_mfdataset: {dt2}')
    elif xy_method == 'grid':        
        t1 = datetime.datetime.now()
        ds = get_mon_di_cy_grid(year, month,
                                indir, f_list,
                                z_lim, t_array_test)
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
    out_fn = 'atmo_' + "{:0>4d}".format(year) + '_' + "{:0>2d}".format(month) + \
        '_' + xy_method + '_all_vars_hourly.nc'
    ds.to_netcdf(path=outdir + out_fn,
                 mode='w',format='NETCDF4',
                 encoding=encdg,unlimited_dims=['time'])
    #
    return 


def time_from_filename(ds):
    ds_out = ds.copy()
    rs = re.compile('flxf([0-9]*).gdas2.*')
    ds_out = ds_out.assign_coords(time=np.array(
        [np.float64(rs.search(ds.encoding['source'].split('/')[-1]).group(1))]
    ))
    ds_out.time.attrs = ds.time.attrs
    return ds_out


def get_hourly_points(yr, mn, data_dir, flist, 
                         lat_list, lon_list, max_depth,
                         interp_how, test_times): 
    
    # Loads all files for one month.
    if flist == 'None':
        flist_air = data_dir + 'flxf*_air.nc'
        flist_sfc = data_dir + 'flxf*_sfc.nc'
    else:
        flist_air = [data_dir + fn + "_air.nc" 
                     for fn in flist.rstrip().split(" ")]
        flist_sfc = [data_dir + fn + "_sfc.nc"
                     for fn in flist.rstrip().split(" ")]
    ds_air = xr.open_mfdataset(
        flist_air,
        engine='netcdf4',
        preprocess=time_from_filename,
        decode_times=False
    )  
    ds_sfc = xr.open_mfdataset(
        flist_sfc,
        engine='netcdf4',
        preprocess=time_from_filename,
        decode_times=False
    )
    ds_air = ds_air.rename({'var33':'UGRD_10m', 'var34':'VGRD_10m',
                            'var11':'TMP_2m', 'var51':'SPFH_2m'})
    ds_sfc = ds_sfc.rename({'var124':'UFLX', 'var125':'VFLX',
                            'var122':'SHTFL', 'var121':'LHTFL',
                            'var204':'DSWRF', 'var205':'DLWRF',
                            'var211':'USWRF', 'var212':'ULWRF',
                            'var59':'PRATE', 'var1':'PRES',
                            'var11':'TMP_SFC'})
    ds_all = ds_air.sel(height=10.0, height_2=2.0)\
            .drop(['height','height_2'])\
            .merge(ds_sfc) 
    
    # Decode and test time array.
    ds_all = xr.decode_cf(ds_all, decode_times=True) 
    if np.all(np.equal(ds_all.time.data, test_times)):
        print('New time variable is as expected...')
    else:
        print('New time variable does not match test values !!')
        import pdb; pdb.set_trace()

    # Interpolate to points.
    ds_sub = ds_all.interp(
        lat=lat_list, lon=lon_list,
        method=interp_how
    )
    
    # Calculate mean diurnal cycle.
    di_cy = ds_sub 
    
    # Copy/add attributes.
    for v in ds_all.data_vars.keys():
        di_cy[v].attrs = ds_sub[v].attrs 
    di_cy.attrs = ds_sub.attrs
    di_cy.attrs['postprocessing'] = 'hourly data, spatially interpolated to points using ' + interp_how 
    
    # Rename variables and add attributes.
    di_cy['UGRD_10m'].attrs['standard_name'] = 'eastward_wind'
    di_cy['UGRD_10m'].attrs['long_name'] = 'eastward_component_of_wind_velocity'
    di_cy['UGRD_10m'].attrs['units'] = 'm s-1'
    di_cy['UGRD_10m'].attrs['level_type'] = 'height above surface'
    di_cy['UGRD_10m'].attrs['height_meters'] = 10.0
    di_cy['VGRD_10m'].attrs['standard_name'] = 'northward_wind'
    di_cy['VGRD_10m'].attrs['long_name'] = 'northward_component_of_wind_velocity'
    di_cy['VGRD_10m'].attrs['units'] = 'm s-1'
    di_cy['VGRD_10m'].attrs['level_type'] = 'height above surface'
    di_cy['VGRD_10m'].attrs['height_meters'] = 10.0
    di_cy['TMP_2m'].attrs['standard_name'] = 'air_temperature'
    di_cy['TMP_2m'].attrs['long_name'] = 'air_temperature'
    di_cy['TMP_2m'].attrs['units'] = 'K'
    di_cy['TMP_2m'].attrs['level_type'] = 'height above surface'
    di_cy['TMP_2m'].attrs['height_meters'] = 2.0
    di_cy['SPFH_2m'].attrs['standard_name'] = 'specific_humidity'
    di_cy['SPFH_2m'].attrs['long_name'] = 'specific_humidity'
    di_cy['SPFH_2m'].attrs['units'] = 'kg kg-1'
    di_cy['SPFH_2m'].attrs['level_type'] = 'height above surface'
    di_cy['SPFH_2m'].attrs['height_meters'] = 2.0
    di_cy['UFLX'].attrs['standard_name'] = 'surface_downward_eastward_stress'
    di_cy['UFLX'].attrs['long_name'] = 'surface_downward_eastward_stress'
    di_cy['UFLX'].attrs['units'] = 'N m-2'
    di_cy['UFLX'].attrs['level_type'] = 'surface'
    di_cy['VFLX'].attrs['long_name'] = 'surface_downward_northward_stress'
    di_cy['VFLX'].attrs['units'] = 'N m-2' 
    di_cy['VFLX'].attrs['standard_name'] = 'surface_downward_northward_stress'
    di_cy['VFLX'].attrs['level_type'] = 'surface'
    di_cy['SHTFL'].attrs['long_name'] = 'surface_upward_sensible_heat_flux'
    di_cy['SHTFL'].attrs['units'] = 'W m-2'
    di_cy['SHTFL'].attrs['standard_name'] = 'surface_upward_sensible_heat_flux'
    di_cy['SHTFL'].attrs['level_type'] = 'surface'
    di_cy['LHTFL'].attrs['long_name'] = 'surface_upward_latent_heat_flux'
    di_cy['LHTFL'].attrs['units'] = 'W m-2'
    di_cy['LHTFL'].attrs['standard_name'] = 'surface_upward_latent_heat_flux'
    di_cy['LHTFL'].attrs['level_type'] = 'surface'
    di_cy['DSWRF'].attrs['long_name'] = 'surface_downwelling_shortwave_flux'
    di_cy['DSWRF'].attrs['units'] = 'W m-2'
    di_cy['DSWRF'].attrs['standard_name'] = 'surface_downwelling_shortwave_flux'
    di_cy['DSWRF'].attrs['level_type'] = 'surface'
    di_cy['DLWRF'].attrs['long_name'] = 'surface_downwelling_longwave_flux'
    di_cy['DLWRF'].attrs['units'] = 'W m-2'
    di_cy['DLWRF'].attrs['standard_name'] = 'surface_downwelling_longwave_flux'
    di_cy['DLWRF'].attrs['level_type'] = 'surface'
    di_cy['USWRF'].attrs['long_name'] = 'surface_upwelling_shortwave_flux'
    di_cy['USWRF'].attrs['units'] = 'W m-2'
    di_cy['USWRF'].attrs['standard_name'] = 'surface_upwelling_shortwave_flux'
    di_cy['USWRF'].attrs['level_type'] = 'surface'
    di_cy['ULWRF'].attrs['long_name'] = 'surface_upwelling_longwave_flux'
    di_cy['ULWRF'].attrs['units'] = 'W m-2'
    di_cy['ULWRF'].attrs['standard_name'] = 'surface_upwelling_longwave_flux'
    di_cy['ULWRF'].attrs['level_type'] = 'surface'
    di_cy['PRATE'].attrs['long_name'] = 'precipitation_flux'
    di_cy['PRATE'].attrs['units'] = 'kg m-2 s-2'
    di_cy['PRATE'].attrs['standard_name'] = 'precipitation_flux'
    di_cy['PRATE'].attrs['level_type'] = 'surface'
    di_cy['PRES'].attrs['long_name'] = 'air_pressure'
    di_cy['PRES'].attrs['units'] = 'Pa'
    di_cy['PRES'].attrs['standard_name'] = 'air_pressure'
    di_cy['PRES'].attrs['level_type'] = 'surface'
    di_cy['TMP_SFC'].attrs['long_name'] = 'surface_temperature'
    di_cy['TMP_SFC'].attrs['units'] = 'K'
    di_cy['TMP_SFC'].attrs['standard_name'] = 'surface_temperature'
    di_cy['TMP_SFC'].attrs['level_type'] = 'surface'
    
    ds_all.close()
    
    return di_cy


def get_mon_di_cy_grid(yr, mn, data_dir, flist,
                       max_depth, test_times):

    var_list = ['var33', 'var34']

    # Loads all files for one month.
    if flist == 'None':
        flist_arg = data_dir + 'flxf*nc'
    else:
        flist_arg = [data_dir + fn + ".nc"
                     for fn in flist.rstrip().split(" ")]
    ds_all = xr.open_mfdataset(
        flist_arg,
        engine='netcdf4',
        preprocess=time_from_filename,
        decode_times=False
    )

    # Decode and test time array.
    ds_all = xr.decode_cf(ds_all, decode_times=True)
    if np.all(np.equal(ds_all.time.data, test_times)):
        print('New time variable is as expected...')
    else:
        print('New time variable does not match test values !!')
        import pdb; pdb.set_trace()

    # Calculate mean diurnal cycle.
    ds_sub = ds_all[var_list]
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

    # Rename variables and add attributes.
    di_cy = di_cy.rename({'var33':'UGRD', 'var34':'VGRD'})
    di_cy['UGRD'].attrs['standard_name'] = 'eastward_wind'
    di_cy['UGRD'].attrs['long_name'] = 'eastward_component_of_wind_velocity'
    di_cy['UGRD'].attrs['units'] = 'm s-1'
    di_cy['VGRD'].attrs['standard_name'] = 'northward_wind'
    di_cy['VGRD'].attrs['long_name'] = 'northward_component_of_wind_velocity'
    di_cy['VGRD'].attrs['units'] = 'm s-1'

    ds_all.close()

    return di_cy


def switch_lon_lims(lon_list, min_lon=0.0):
    result = (lon_list - min_lon) % 360.0 + min_lon
    return result


if __name__ == '__main__': 
    cla = sys.argv
    main(cla[1], cla[2], cla[3], cla[4], cla[5], cla[6])
