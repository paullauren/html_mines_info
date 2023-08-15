import pystac_client
import planetary_computer
import rasterio as rio
import rasterio.mask
from matplotlib import pyplot as plt
import numpy as np
import rich
import pandas as pd
import xarray as xr
import stackstac
from IPython.display import display
import geopandas as gpd
import rioxarray
from rasterio import features
import time
import rasterstats
from mpl_toolkits.axes_grid1 import make_axes_locatable

catalog = pystac_client.Client.open(
    "https://planetarycomputer.microsoft.com/api/stac/v1",
    modifier=planetary_computer.sign_inplace,
    )

item = catalog.get_collection('sentinel-2-l2a').get_item("S2B_MSIL2A_20230813T013659_R117_T53SNU_20230813T063926")

aoi_bounds = [135.824621, 34.955573, 136.004300, 35.084910]
bands = ['B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B11', 'B12']
FILL_VALUE = 2**16-1

ds = stackstac.stack(
                    item,
                    assets = bands,
                    resolution=10,
                    dtype="uint16",
                    fill_value=FILL_VALUE,
                    bounds_latlon=aoi_bounds,
                    )

nir = ds.sel(band="B08").astype('float')
swir = ds.sel(band="B11").astype('float')
red = ds.sel(band="B04").astype('float')
green = ds.sel(band="B03").astype('float')
rededge = ds.sel(band="B05").astype('float')
blue = ds.sel(band="B02").astype('float')

#Chlorophyll vegetation index
cvi = nir*red/green**2
#Global vegetation moisture index
gvmi = ((nir+0.1)-(swir+0.02))/((nir+0.1)+(swir+0.02))
#Blue-wide dynamic range vegetation index
bwdrvi = (0.1*nir-blue)/(0.1*nir+blue)

cvi_np=cvi.to_numpy()
gvmi_np=gvmi.to_numpy()
bwdrvi_np=bwdrvi.to_numpy()

def normalize(array, vmax=False):
    """Normalizes numpy arrays into scale 0.0 - 1.0"""
    if vmax:
        array[array>vmax] = vmax
    array_min, array_max = array.min(), array.max()
    return ((array - array_min)/(array_max - array_min))

cvi_n = normalize(cvi_np)
gvmi_n = normalize(gvmi_np)
bwdrvi_n = normalize(bwdrvi_np)

rgb = np.stack((cvi_n, gvmi_n, bwdrvi_n))
rgb_natcol = np.moveaxis(rgb.squeeze(), 0, -1)

plt.figure()

fig, axs = plt.subplots(2, 2)
axs[0, 0].imshow(cvi.squeeze(), cmap="RdYlGn")
axs[0, 0].title.set_text('CVI')
divider = make_axes_locatable(axs[0,0])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(axs[0, 0].imshow(cvi.squeeze(), cmap="RdYlGn"), cax=cax, orientation='vertical')

axs[0, 1].imshow(gvmi.squeeze(), cmap="RdYlGn")
axs[0, 1].title.set_text('GVMI')
divider = make_axes_locatable(axs[0,1])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(axs[0, 1].imshow(gvmi.squeeze(), cmap="RdYlGn"), cax=cax, orientation='vertical')

axs[1, 0].imshow(bwdrvi.squeeze(), cmap="RdYlGn")
axs[1, 0].title.set_text('BWDRVI')
divider = make_axes_locatable(axs[1,0])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(axs[1, 0].imshow(bwdrvi.squeeze(), cmap="RdYlGn"), cax=cax, orientation='vertical')

axs[1, 1].imshow(rgb_natcol)
axs[1, 1].title.set_text('Composite of CVI, GVMI and BWDRVI')
divider = make_axes_locatable(axs[1,1])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(axs[1, 1].imshow(rgb_natcol), cmap="RdYlGn", cax=cax, orientation='vertical')

plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.4,
                    hspace=0.4)

plt.show()
