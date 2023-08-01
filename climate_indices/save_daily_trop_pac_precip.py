import numpy as np
import xarray as xr
import pandas as pd
import pathlib
import sys

data_year = int(sys.argv[1])

data_dir = '/lustre/f2/dev/ncep/JieShun.Zhu/CFSm501hr/CFSm501hr1980010200/DATA/'

# Create lat and lon mask from a single day file.
ds1 = xr.open_dataset(data_dir + 'ocn_2002_01_01_00.nc')
pacmask = ((ds1.geolat_t >= -10.0) & (ds1.geolat_t <= 10.0))
pacmaskx = pacmask.any(dim='yt_ocean')
pacmasky = pacmask.any(dim='xt_ocean')

# Get data from multiple files for this masked area, and average to monthly.
dy_list = []
for dt64 in np.arange(str(data_year) + '-01-01', 
                      str(data_year + 1) + '-01-01', 
                      dtype='datetime64[D]'):
    yyyy = pd.to_datetime(str(dt64)).strftime('%Y')
    mm = pd.to_datetime(str(dt64)).strftime('%m')
    dd = pd.to_datetime(str(dt64)).strftime('%d')
    print('Opening ocn_' + yyyy + '_' + mm + '_' + dd + '*.nc')
    tall = xr.open_mfdataset(data_dir + 'ocn_' + yyyy + '_' + 
                             mm + '_' + dd + '*.nc')\
       .lprec[dict(xt_ocean=pacmaskx, yt_ocean=pacmasky)]\
       .mean(dim='yt_ocean').mean(dim='time')\
       .assign_coords({'time':dt64 + np.timedelta64(12,'h')})
    dy_list.append(tall.load())

# Join months together.
ds_dy = xr.concat(dy_list, dim='time').to_dataset()

# Sort attributes and encoding.
ds_dy.lprec.attrs = ds1.lprec.attrs
ds_dy.time.attrs = {
    'long_name':'time',
    'standard_name':'time',
    'description':'time at 12 UTC each day'
}
ds_dy.attrs = {
    'description':'CFSm501 daily mean, meridional mean precipitation from daily data',
    'contact':'Jack Reeves Eyre (jack.reeveseyre@noaa.gov)'
}
fill_val_f64 = 9.969209968386869e+36
fill_val_i32 = np.iinfo(np.int32).min + 2
fill_val_i64 = np.iinfo(np.int64).min + 2
t_units = 'seconds since 1970-01-01 00:00:00'
t_cal = 'standard'
encdg = {}
for k in ds_dy.data_vars.keys():
    encdg[k] = {'_FillValue':fill_val_f64, 'dtype':'float64'}
for k in ds_dy.coords.keys():
    encdg[k] = {'_FillValue':fill_val_f64, 'dtype':'float64'}
encdg['time']['units'] = t_units
encdg['time']['calendar'] = t_cal

# Save out file.
out_dir = '/lustre/f2/dev/ncep/Jack.Reeveseyre/climate_indices/'
pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)
out_fn = 'ocn_precip_pacific_10S_10N_' + str(data_year) + '_daily_mean.nc'
print('Saving file ' + out_dir + out_fn)
ds_dy.to_netcdf(path = out_dir + out_fn,
                mode='w', format='NETCDF4',
                encoding=encdg, unlimited_dims=['time'])
