import xarray as xr
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def main():
    
    # Specify which locations to plot for.
    latlon_pairs = np.array([[2.0, 220.0], [0.0, 220.0]])

    # Load data.
    ddir = '/lustre/f2/dev/ncep/Jack.Reeveseyre/diurnal_cycle/'
    ds_atmo = xr.open_mfdataset(ddir + "atmo_*_points_meanDiurnalCycle.nc")\
            .chunk({'time':48}).mean(dim='time').sel(height=10.0)
    ds_ocean = xr.open_mfdataset(ddir + "ocn_*_meanDiurnalCycle.nc")\
            .chunk({'time':48}).mean(dim='time')
    
    # Loop over locations and plot.
    for i in range(latlon_pairs.shape[0]): 
        ck_plot(ds_ocean.sel(yt_ocean=latlon_pairs[i,0],
                             yu_ocean=latlon_pairs[i,0],
                             xt_ocean=switch_lon_lims(latlon_pairs[i,1], 
                                                      min_lon=-280.0),
                             xu_ocean=switch_lon_lims(latlon_pairs[i,1], 
                                                      min_lon=-280.0)), 
                ds_atmo.sel(lat=latlon_pairs[i,0],
                            lon=latlon_pairs[i,1]))

    return


def switch_lon_lims(lon_list, min_lon=0.0):
    result = (lon_list - min_lon) % 360.0 + min_lon
    return result


def ck_plot(dso, dsa): 
    
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

    # Set up figure.
    fig, axs = plt.subplots(1,1, figsize=(6,4), clip_on=False)
    axs.set_title(lat_str + ' ' + lon_str, loc='left', pad=20)
    axs.set_xlim(-0.5, 23.5)
    ymin = 0.0
    ymax = 26.0
    axs.set_ylim(ymax, ymin)
    axs.set_xlabel('time (local)')
    axs.set_ylabel('depth (m)')
    axs.tick_params(axis='x', which='both', bottom=True, top=True)
    axs.tick_params(axis='y', which='both', left=True, right=True)
    
    # Work out color scale.
    col_min = dso['temp'].sel(st_ocean=slice(ymin, ymax)).min()
    col_max = dso['temp'].sel(st_ocean=slice(ymin, ymax)).max()
    bounds = np.linspace(np.floor(col_min*10.0)/10.0, 
                         np.ceil(col_max*10.0)/10.0, 11)
    norm = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=256, extend='both')
    # Add ocean data.
    x2d, y2d = np.meshgrid(dso['hour'], dso['st_ocean'], indexing='ij')
    pcm = axs.pcolormesh(x2d, y2d, dso['temp'], shading='nearest',
                         cmap='RdBu_r', norm=norm)
    plt.colorbar(pcm, location='right')
    dso_base = dso.sel(st_ocean=25.0)
    axs.quiver(x2d, y2d, dso['u'] - dso_base['u'], 
               dso['v'] - dso_base['v'])
    
    # Add atmo data.
    axs.quiver(dsa['hour'], 0.0*xr.ones_like(dsa['hour']), 
               dsa['UGRD'], dsa['VGRD'])

    # Save file.
    plotdir = '/ncrc/home2/Jack.Reeveseyre/cfs-analysis/diurnal_cycle/plots/'
    plotfn = 'croninkessler_plot_' + lat_str + '_' + lon_str + '.'
    plotfileformat='pdf'
    plt.savefig(plotdir + plotfn + plotfileformat, format=plotfileformat)
    return


if __name__ == "__main__":
    main()
