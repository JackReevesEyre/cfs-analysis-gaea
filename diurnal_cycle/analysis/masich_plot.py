import xarray as xr
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import sys
import utils

def main(plot_month):
    
    # Specify which locations to plot for.
    latlon_pairs = np.array([[2.0, 220.0], [0.0, 220.0]])

    # Load data.
    ddir = '/lustre/f2/dev/ncep/Jack.Reeveseyre/diurnal_cycle/'
    ds_atmo = xr.open_mfdataset(ddir + "atmo_*_points_meanDiurnalCycle.nc")\
            .chunk({'time':48}).sel(height=10.0)
    ds_atmo_rec_mean = ds_atmo.mean(dim=['time','hour'], keep_attrs=True)
    ds_ocean = xr.open_mfdataset(ddir + "ocn_*_meanDiurnalCycle.nc")\
            .chunk({'time':48})
    
    # Take time average (optionally for a subset of months).
    if plot_month == 'ALL':
        ds_atmo = ds_atmo.mean(dim='time', keep_attrs=True)
        ds_ocean = ds_ocean.mean(dim='time', keep_attrs=True)
        month_str = 'ANN'
    else:
        ds_atmo = ds_atmo.sel(time=(ds_atmo.time.dt.month == plot_month))\
            .mean(dim='time', keep_attrs=True)
        ds_ocean = ds_ocean.sel(time=(ds_ocean.time.dt.month == plot_month))\
            .mean(dim='time', keep_attrs=True)
        month_names = ['JAN','FEB','MAR','APR','MAY','JUN',
                       'JUL','AUG','SEP','OCT','NOV','DEC']
        month_str = month_names[plot_month - 1]

    # Calculate down-wind current.
    ds_atmo_rec_mean = utils.wind_dir_twice(ds_atmo_rec_mean, 'UGRD', 'VGRD')
    ds_ocean = utils.rotation_xy(ds_ocean, 'u', 'v', 
                                 np.radians(ds_atmo_rec_mean['arctan_VoverU'].data), 
                                 'math', exp_dim=(0,1))

    # Loop over locations and plot.
    for i in range(latlon_pairs.shape[0]): 
        ma_plot(ds_ocean.sel(yt_ocean=latlon_pairs[i,0],
                             yu_ocean=latlon_pairs[i,0],
                             xt_ocean=switch_lon_lims(latlon_pairs[i,1], 
                                                      min_lon=-280.0),
                             xu_ocean=switch_lon_lims(latlon_pairs[i,1], 
                                                      min_lon=-280.0)), 
                ds_atmo.sel(lat=latlon_pairs[i,0],
                            lon=latlon_pairs[i,1]),
                month_str)

    return


def switch_lon_lims(lon_list, min_lon=0.0):
    result = (lon_list - min_lon) % 360.0 + min_lon
    return result


def ma_plot(dso, dsa, month, max_depth=60.0): 
    
    # Format lat and lon strings.
    if dsa['lat'] >= 0.0:
        lat_str = "{:.1f}".format(dsa['lat'].data) + 'N'
    else:
        lat_str = "{:.1f}".format(np.abs(dsa['lat'].data)) + 'S'
    if dsa['lon'] > 180.0:
        lon_str = "{:.1f}".format(360.0 - dsa['lon'].data) + 'W'
    else:
        lon_str = "{:.1f}".format(dso['lon'].data) + 'E'
    
    # Change hour to local time.
    hours_ahead = np.around(dsa['lon'].data/15.0)
    #dsa = dsa.assign_coords(hour=((dsa.hour + hours_ahead) % 24))
    #dso = dso.assign_coords(hour=((dso.hour + hours_ahead) % 24))
    dsa = dsa.roll(hour=int(hours_ahead), roll_coords=False)
    dso = dso.roll(hour=int(hours_ahead), roll_coords=False)

    # Calculate anomalies.
    z = dso.st_ocean.data
    st_ocean_base = np.nanmax(z[z <= max_depth])
    dso_anom = dso.copy()
    for varname in ['u','v','U_streamwise','U_crossstream']:
        dso_anom[varname] = dso_anom[varname] \
            - dso_anom[varname].sel(st_ocean=st_ocean_base)
        dso_anom[varname] = 100.0*dso_anom[varname]
        dso_anom[varname].attrs['units'] = 'cm s-1'
    dso_anom = dso_anom - dso_anom.mean(dim='hour')

    # Set up main figure.
    fig = plt.figure(figsize=(6,4))
    gs = GridSpec(2, 2, height_ratios=[4,1], hspace=0.3)
    axs = fig.add_subplot(gs[0,:])
    axs.set_title(lat_str + ' ' + lon_str, loc='left')
    axs.set_title(month, loc='right')
    axs.set_xlim(-0.5, 23.5)
    ymin = 0.0
    ymax = max_depth
    axs.set_ylim(ymax, ymin)
    axs.set_xlabel('time (local)')
    axs.set_ylabel('depth (m)')
    axs.tick_params(axis='x', which='both', bottom=True, top=True)
    axs.tick_params(axis='y', which='both', left=True, right=True)
    cbl = fig.add_subplot(gs[1,0])
    cbr = fig.add_subplot(gs[1,1])
    
    # Plot temperature data.
    tcol_min = -0.15
    tcol_max = 0.15
    tbounds = np.linspace(tcol_min, tcol_max, 11)
    tnorm = mpl.colors.BoundaryNorm(boundaries=tbounds, ncolors=256, extend='both')
    x2d, y2d = np.meshgrid(dso_anom['hour'], dso_anom['st_ocean'], 
                           indexing='ij')
    pcm = axs.pcolormesh(x2d, y2d, dso_anom['temp'], 
                         shading='nearest',
                         cmap='RdBu_r', norm=tnorm)
    plt.colorbar(pcm, cax=cbl, extend='both', orientation='horizontal',
                 label='temperature anomaly (' + r'$^{\circ}$' + 'C)')

    # Plot velocity data.
    ucol_min = -3.0
    ucol_max = 3.0
    ubounds = np.linspace(ucol_min, ucol_max, 11)
    unorm = mpl.colors.BoundaryNorm(boundaries=ubounds, ncolors=256, extend='both')
    clevels = np.arange(-3.0, 3.1, 0.5)
    cltypes = np.array(clevels, dtype='object_')
    cltypes[:] = 'solid'
    cltypes[len(clevels)//2] = 'dashed'
    cont = axs.contour(x2d, y2d, dso_anom['U_streamwise'],
                       levels=clevels,
                       cmap='PuOr', norm=unorm,
                       linewidths=0.7, linestyles=cltypes)
    axs.clabel(cont, np.array([-3, -2, -1, 1, 2, 3]), fmt='%d', 
               fontsize=7.0, colors='k')
    plt.colorbar(cont, cax=cbr, extend='both', orientation='horizontal', 
                 label='along-wind velocity anomaly (' + r'$cm\ s^{-1}$' + ')')
    
    # Save file.
    plotdir = '/ncrc/home2/Jack.Reeveseyre/cfs-analysis/diurnal_cycle/plots/'
    plotfn = 'masich_plot_' + lat_str + '_' + lon_str + '_' + month + '.'
    plotfileformat='pdf'
    plt.savefig(plotdir + plotfn + plotfileformat, format=plotfileformat)
    return


if __name__ == "__main__":
    if len(sys.argv) == 1:
        plot_month = 'ALL'
    else:
        plot_month = int(sys.argv[1])
    main(plot_month)
