import xarray as xr
import numpy as np
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
    
    # Set up figure.
    fig, axs = plt.subplots(1,1, figsize=(5,8))
    axs.set_title(lat_str + ' ' + lon_str, loc='left', pad=30)
    axs.set_xlim(-0.5, 23.5)
    axs.set_ylim(25.0, 0.0)
    axs.set_xlabel('time (UTC)')
    axs.set_ylabel('depth (m)')
    axs.tick_params(axis='x', which='both', bottom=True, top=True)
    axs.tick_params(axis='y', which='both', left=True, right=True)
    
    # Add ocean data.
    x2d, y2d = np.meshgrid(dso['hour'], dso['st_ocean'], indexing='ij')
    axs.pcolormesh(x2d, y2d, dso['temp'], shading='nearest')
    axs.quiver(x2d, y2d, dso['u'], dso['v'])
    
    # Add atmo data.
    axs.quiver(dsa['hour'], -5.0, dsa['UGRD'], dsa['VGRD'])
    import pdb; pdb.set_trace()

    # Save file.
    plotdir = '/ncrc/home2/Jack.Reeveseyre/plots/diurnal_cycle/'
    plotfn = 'croninkessler_plot_' + lat_str + '_' + lon_str + '.'
    plotfileformat='pdf'
    plt.savefig(plotdir + plotfn + plotfileformat, format=plotfileformat,
                dpi=500)
    return


if __name__ == "__main__":
    main()
