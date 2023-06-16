import numpy as np
import xarray as xr
import pathlib

dy_dir = '/lustre/f2/dev/ncep/JieShun.Zhu/CFSm501dy/CFSm501dy1980010200/DATA/'
hr_dir = '/lustre/f2/dev/ncep/JieShun.Zhu/CFSm501hr/CFSm501hr1980010200/DATA/'

# Create lat and lon mask from a single day file.
ds1 = xr.open_dataset(dy_dir + 'ocn_1982_01_01_00.nc')
pacmask = ((ds1.geolat_t >= -5.0) & (ds1.geolat_t <= 5.0) & (ds1.geolon_t >= -230.0) & (ds1.geolon_t <= -80.0))
pacmaskx = pacmask.any(dim='yt_ocean')
pacmasky = pacmask.any(dim='xt_ocean')

# Get data from multiple files for this masked area, and average to monthly.
mon_list = []
for yr in np.arange(1981, 2006):
    for mn in np.arange(1,13):
        if yr >= 2002:
            data_dir = hr_dir
        else:
            data_dir = dy_dir
        print('Opening ' + data_dir + 'ocn_' + "{:0>4d}".format(yr) + 
              '_' + "{:0>2d}".format(mn) + '*.nc')
        try:
            tall = xr.open_mfdataset(data_dir + 
                                     'ocn_' + "{:0>4d}".format(yr) + 
                                     '_' + "{:0>2d}".format(mn) + '*.nc')\
                   .temp.sel(st_ocean=0.5)[dict(xt_ocean=pacmaskx, yt_ocean=pacmasky)]\
                   .mean(dim='yt_ocean').mean(dim='time')\
                   .assign_coords({'time':np.datetime64("{:0>4d}".format(yr) + '-' + 
                                                        "{:0>2d}".format(mn) + 
                                                        '-15T00:00:00')})
            mon_list.append(tall)
        except:
            print('No data for this month.')

# Join months together.
ds_mon = xr.concat(mon_list, dim='time').to_dataset()

# Sort attributes and encoding.
ds_mon.temp.attrs = ds1.temp.attrs
ds_mon.time.attrs = {
    'long_name':'time',
    'standard_name':'time',
    'description':'time at 00 UTC on 15th of the month'
}
ds_mon.attrs = {
    'description':'CFSm501 monthly mean, meridional mean surface (top level) temperature from daily data',
    'contact':'Jack Reeves Eyre (jack.reeveseyre@noaa.gov)'
}
fill_val_f64 = 9.969209968386869e+36
fill_val_i32 = np.iinfo(np.int32).min + 2
fill_val_i64 = np.iinfo(np.int64).min + 2
t_units = 'seconds since 1970-01-01 00:00:00'
t_cal = 'standard'
encdg = {}
for k in ds_mon.data_vars.keys():
    encdg[k] = {'_FillValue':fill_val_f64, 'dtype':'float64'}
for k in ds_mon.coords.keys():
    encdg[k] = {'_FillValue':fill_val_f64, 'dtype':'float64'}
encdg['time']['units'] = t_units
encdg['time']['calendar'] = t_cal

# Save out file.
out_dir = '/lustre/f2/dev/ncep/Jack.Reeveseyre/climate_indices/'
pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)
out_fn = 'ocn_sfc_temp_pacific_5S_5N_monthly_mean.nc'
ds_mon.to_netcdf(path = out_dir + out_fn,
                 mode='w', format='NETCDF4',
                 encoding=encdg, unlimited_dims=['time'])

