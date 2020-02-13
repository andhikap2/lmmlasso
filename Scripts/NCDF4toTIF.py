#ncdf4 things

from netCDF4 import Dataset
import numpy as np
import rasterio
from rasterio.transform import from_origin


my_example_nc_file = 'tn_ens_mean_0.1deg_reg_v19.0eHOM.nc'
fh = Dataset(my_example_nc_file, mode='r')

lons = fh.variables['longitude'][:]
lats = fh.variables['latitude'][:]
tn = fh.variables['tn'][:] #load data for the 100th day
tn_units = fh.variables['tn'].units


fh.close() #close the file once you're done getting the info from it

import matplotlib.pyplot as plt 
from mpl_toolkits.basemap import Basemap

# Get some parameters for the Stereographic Projection
lon_0 = lons.mean()
lat_0 = lats.mean()

m = Basemap(width=5000000,height=3500000,
            resolution='l',projection='stere',\
            lat_ts=40,lat_0=lat_0,lon_0=lon_0)

lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)

# Plot Data
cs = m.pcolor(xi,yi,np.squeeze(tn))

# Add Grid Lines
m.drawparallels(np.arange(-80., 81., 10.), labels=[1,0,0,0], fontsize=10)
m.drawmeridians(np.arange(-180., 181., 10.), labels=[0,0,0,1], fontsize=10)

# Add Coastlines, States, and Country Boundaries
m.drawcoastlines()
m.drawstates()
m.drawcountries()

# Add Colorbar
cbar = m.colorbar(cs, location='bottom', pad="10%")
cbar.set_label(tn_units)

# Add Title
plt.title('Minimum Temperature')

plt.show()




'''loading ncdf4 and converting it into a tiff file'''

from netCDF4 import Dataset
import numpy as np
import rasterio
from rasterio.transform import from_origin
import rioxarray
import xarray

xds = xarray.open_dataset('tx_ens_mean_0.1deg_reg_v19.0eHOM.nc')
xds["tx"].rio.to_raster('tx_ens_mean_0.1deg_reg_v19.0eHOM.tif')



'''loading ncdf4 and converting it into a tiff file'''

from netCDF4 import Dataset
import numpy as np
import rasterio
from rasterio.transform import from_origin
import rioxarray
import xarray

xds = xarray.open_dataset('tn_ens_mean_0.1deg_reg_v19.0eHOM.nc.1')
xds["tn"].rio.to_raster('tn_ens_mean_0.1deg_reg_v19.0eHOM.tif')


lradin <- brick("tn_ens_mean_0.1deg_reg_v19.0eHOM.nc.1")

# Register GDAL format drivers and configuration options with a
# context manager.
with rasterio.Env():
    profile = src.profile
    profile.update(
        dtype=rasterio.uint8,
        count=25171,
        compress='lzw')
    with rasterio.open('example.tif', 'w', **profile) as dst:
        dst.write(array.astype(rasterio.uint8), (0:25171))

# At the end of the ``with rasterio.Env()`` block, context
# manager exits and all drivers are de-registered.