import xarray as xr
import pandas as pd
import numpy as np

def read_haddtr(clim_type):
    """Reads HadDTR (diurnal range of SST); returns xarray dataset.

    Input:
        clim_type -- [string] One of 'monthly', 'seasonal', 'annual'
    """
    
    # -------------------- Open file. --------------------
    #ddir = '/home/jackre/Data/gridded-obs/HadDTR/'
    ddir = '/ncrc/home2/Jack.Reeveseyre/data/haddtr/'
    fn = clim_type + '_climatology.txt'
    fc = open(ddir + fn, 'r')

    # -------------------- Create empty data array. --------------------
    lon = np.arange(-177.5, 177.7, 5.0)
    lat = np.arange(67.5, -67.6, -5.0)
    nt = {'monthly':12, 'seasonal':4, 'annual':1}[clim_type]
    time = np.arange(nt)
    n_bounds = np.arange(2)
    year_bounds = np.full((nt, len(n_bounds)), np.nan, dtype=np.int64)
    lat_bounds = np.full((len(lat), len(n_bounds)), np.nan, dtype=np.float64)
    lon_bounds = np.full((len(lon), len(n_bounds)), np.nan, dtype=np.float64)
    ds = xr.Dataset(
        data_vars=dict(
            dtr=(['time_index', 'lat', 'lon'],
                 np.full((nt, len(lat), len(lon)), np.nan)),
            time_name=(['time_index'],
                       np.full((nt), 'missing', dtype='object'))
        ),
        coords=dict(
            lon=lon,
            lat=lat,
            time_index=time,
            n_bounds=n_bounds,
            climatology_bounds=(['time_index', 'n_bounds'], year_bounds),
            lat_bounds=(['lat', 'n_bounds'], lat_bounds),
            lon_bounds=(['lon', 'n_bounds'], lon_bounds)
        )
    )

    # -------------------- Read data and insert into data array. ---------------
    it = 0
    while it < nt:

        # Get header for this time step.
        lc = fc.readline().strip().split()
        while lc == []:
            lc = fc.readline().strip().split()
        header = {lcc.split('=')[0]:lcc.split('=')[1]
                  for lcc in lc if len(lcc.split('=')) == 2}
        ds['climatology_bounds'].loc[dict(time_index=it, n_bounds=0)] = \
            int(header['start_year'])
        ds['climatology_bounds'].loc[dict(time_index=it, n_bounds=1)] = \
            int(header['end_year'])
        if clim_type == 'monthly':
            ds['time_name'].loc[dict(time_index=it)] = header['month']
        elif clim_type == 'seasonal':
            ds['time_name'].loc[dict(time_index=it)] = header['season'] + \
                ' ' + ''.join([lcc for lcc in lc if len(lcc.split('=')) == 1])
        else:
           ds['time_name'].loc[dict(time_index=it)] = 'annual'
        assert int(header['rows']) == len(lat), \
            'Number of data rows inconsistent with expected lat array.'
        assert int(header['columns']) == len(lon), \
            'Number of data columns inconsistent with expected lon array.'

        # Loop over rows and put data into array.
        ir = 0
        while ir < len(lat):
            lc = fc.readline().strip().split()
            ds['dtr'].loc[dict(time_index=it, lat=ds.lat[ir])] = \
                np.array(lc, dtype=np.float64)
            #
            ir = ir + 1
        #
        it = it + 1
        
    # Convert -99.99 to nans.
    ds['dtr'] = ds.dtr.where(ds.dtr != -99.99, np.nan)
    
    # Add lat and lon bounds.
    dlat = ds['lat'].data[1] - ds['lat'].data[0] 
    dlon = ds['lon'].data[1] - ds['lon'].data[0]
    ds['lat_bounds'].loc[dict(n_bounds=0)] = ds['lat'].data - (dlat/2.0)
    ds['lat_bounds'].loc[dict(n_bounds=1)] = ds['lat'].data + (dlat/2.0)
    ds['lon_bounds'].loc[dict(n_bounds=0)] = ds['lon'].data - (dlon/2.0)
    ds['lon_bounds'].loc[dict(n_bounds=1)] = ds['lon'].data + (dlon/2.0)

    # -------------------- Add attributes. --------------------
    ds['lat'].attrs = {'units':'degrees_north',
                       'standard_name':'latitude',
                       'description':'grid cell center',
                       'bounds':'lat_bounds'}
    ds['lon'].attrs = {'units':'degrees_east',
                       'standard_name':'longitude',
                       'description':'grid cell center',
                       'bounds':'lon_bounds'}
    ds['lat_bounds'].attrs = {'units':'degrees_north',
                              'standard_name':'latitude',
                              'description':'grid cell edges'}
    ds['lon_bounds'].attrs = {'units':'degrees_east',
                              'standard_name':'longitude',
                              'description':'grid cell edges'}
    ds['time_index'].attrs = {'long_name':'which time step this data is for',
                              'instance_dimension':'time_name',
                              'bounds':'climatology_bounds'}
    ds['climatology_bounds'].attrs = {'units':'year'}
    ds['dtr'].attrs = {'units':'K',
                       'long_name':'mean_diurnal_range_of_sea_surface_temperature',
                       'cell_methods':'time: mean over years, time: range over days, area: mean.'}
    if clim_type == 'monthly':
        ds['time_name'].attrs = {'climatology_type':'monthly',
                                 'units':'1=January, 2=February, ..., 12=December'}
    elif clim_type == 'seasonal':
        ds['time_name'].attrs = {'climatology_type':'seasonal',
                                 'units':'season_description'}
    else:
        ds['time_name'].attrs = {'climatology_type':'annual',
                                 'units':''}
    ds.attrs = {
        'title':'Hadley Centre climatology of the diurnal temperature range of the sea surface based on the ICOADS',
        'institution':'Met Office Hadley Center',
        'source':'https://www.metoffice.gov.uk/hadobs/haddtr/',
        'references':'Kennedy J.J., Brohan P. and S.F.B.Tett, 2007: A global climatology of the diurnal variations in sea-surface temperature and implications for MSU temperature trends. Geophys. Res. Lett., Vol. 34, No. 5, L05712 doi:10.1029/2006GL028920',
        'history':'Date accessed: 2022-10-04.'
    }
    
    #
    return ds


def save_haddtr_netcdf(clim_type):

    ds = read_haddtr(clim_type)
    
    #ddir = '/home/jackre/Data/gridded-obs/HadDTR/'
    ddir = '/ncrc/home2/Jack.Reeveseyre/data/haddtr/'
    fn_out = clim_type + '_climatology.nc'
    
    # Define encoding.
    nc_enc = {}
    fill_val_f64 = 9.969209968386869e+36
    fill_val_i32 = np.iinfo(np.int32).min + 2
    fill_val_i64 = np.iinfo(np.int64).min + 2
    for k in ds.keys():
        if ds[k].dtype == np.float64:
            nc_enc[k] = {'_FillValue':fill_val_f64}
        elif ds[k].dtype in (np.int64, np.int32):
            nc_enc[k] = {'_FillValue':fill_val_i32, 'dtype':'int32'}
        elif ds[k].dtype == 'object':
            nc_enc[k] = {'dtype':'S1', '_Encoding':'utf_8'}
        else:
            print('manual encoding not defined for data type ' +
                  str(ds[k].dtype) + ' (' + k + ').')
    for k in ds.coords.keys():
        if ds[k].dtype == np.float64:
            nc_enc[k] = {'_FillValue':fill_val_f64}
        elif ds[k].dtype in (np.int64, np.int32):
            nc_enc[k] = {'_FillValue':fill_val_i32, 'dtype':'int32'}
        elif ds[k].dtype == 'object':
            nc_enc[k] = {'dtype':'S1', '_Encoding':'utf_8'}
        else:
            print('manual encoding not defined for data type ' +
                  str(ds[k].dtype) + ' (' + k + ').')
            
    # Write file.
    ds.to_netcdf(
        path=ddir + fn_out,
        mode='w',
        format='NETCDF3_CLASSIC',
        encoding=nc_enc,
        unlimited_dims=None
    )
    #
    return

