import numpy as np
import xarray as xr
import pathlib
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import sys

data_dir = '/lustre/f2/dev/ncep/Jack.Reeveseyre/climate_indices/'
plot_dir = '/ncrc/home2/Jack.Reeveseyre/cfs-analysis/climate_indices/plots/'
fn = 'ocn_sfc_temp_pacific_5S_5N_monthly_mean.nc'
ds = xr.open_dataset(data_dir + fn)
#ds = ds.sel(time=slice('2001-01-01', '2005-12-31'))

# Calculate anomalies.
ds_clim = ds.groupby('time.month').mean()
ds_anom = ds.groupby('time.month') - ds_clim

# Plot data.
fig, ax = plt.subplots(figsize=(8, 10))
ax.invert_yaxis()
cf = ax.pcolormesh(ds_anom.xt_ocean,
                   ds_anom.time,
                   ds_anom.temp, shading='nearest',
                   cmap='RdBu_r', vmin=-2.0, vmax=2.0)
cbar = plt.colorbar(cf, orientation='horizontal', pad=0.04, aspect=50)
cbar.set_label('K')

# Plot details.
ax.set_yticks(ds_anom.time[0::24] - np.timedelta64(14, 'D'))
ax.set_yticks(ds_anom.time[0::12] - np.timedelta64(14, 'D'), minor=True)
ax.set_xticks(np.arange(-360.0, 361.0, 30.0))
ax.xaxis.set_major_formatter(LongitudeFormatter(zero_direction_label=False,
                                                degree_symbol=''))
ax.set_xticks(np.arange(-360.0, 361.0, 15.0), minor=True)
ax.set_xlim(-230.0, -80.0)
plt.title('5S-5N SST anomaly', loc='left')

# Save figure.
plot_fn = 'hovmoller_ocn_sfc_temp_pacific_5S_5N_monthly_anomaly'
plot_format = 'png'
fig.savefig(plot_dir + plot_fn + '.' + plot_format,
            format=plot_format, dpi=300)
