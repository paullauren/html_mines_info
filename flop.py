import pystac_client
import planetary_computer
import rasterio as rio
import rasterio.mask
import mpl_scatter_density
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
import matplotlib
from math import *

catalog = pystac_client.Client.open(
    "https://planetarycomputer.microsoft.com/api/stac/v1",
    modifier=planetary_computer.sign_inplace,
    )

la_teste_de_buch = [("S2B_MSIL2A_20220617T105629_R094_T30TXQ_20220619T094751", "June 2022"), ("S2A_MSIL2A_20220801T105631_R094_T30TXQ_20220802T084538","August 2022"), ("S2B_MSIL2A_20230722T105629_R094_T30TXQ_20230722T155341","July 2023")]
landiras = [("S2B_MSIL2A_20220617T105629_R094_T30TXQ_20220619T094751","June 2022"), ("S2A_MSIL2A_20220801T105631_R094_T30TXQ_20220802T084538","August 2022"), ("S2A_MSIL2A_20230727T105621_R094_T30TXQ_20230727T201740","July 2023")]

aoi_ltdb = [-1.290052, 44.419104, -1.082548, 44.614408]
#aoi_landiras = [-0.693548, 44.466738, -0.369085, 44.596085]
aoi_landiras = [-0.693548, 44.43, -0.369085, 44.596085]

"""
##NBR
fig, axs = plt.subplots(1, 3)

for i in range(3):
    name, date = landiras[i]

    item = catalog.get_collection('sentinel-2-l2a').get_item(name)

    boi = ['B02','B03','B04','B05','B06','B07','B08','B11','B12']

    FILL_VALUE = 2**16-1

    ds = stackstac.stack(
                item,
                assets = boi,
                resolution=20,
                dtype="uint16",
                fill_value=FILL_VALUE,
                bounds_latlon=aoi_landiras,
                    )
    
    nir = ds.sel(band="B08").astype('float')
    sir = ds.sel(band="B12").astype('float')

    nbr = (nir-sir)/(nir+sir)

    axs[i].imshow(nbr.squeeze(), cmap="RdYlGn")
    axs[i].title.set_text(f"NBR: {date}")
    divider = make_axes_locatable(axs[i])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(axs[i].imshow(nbr.squeeze(), cmap="RdYlGn", vmin=-0.5, vmax=0.5), cax=cax, orientation='vertical')

plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.4,
                    hspace=0.4)

plt.show()
"""

##dNBR

item = catalog.get_collection('sentinel-2-l2a').get_item("S2A_MSIL2A_20220801T105631_R094_T30TXQ_20220802T084538")
item2 = catalog.get_collection('sentinel-2-l2a').get_item("S2A_MSIL2A_20230727T105621_R094_T30TXQ_20230727T201740")
boi = ['B02','B03','B04','B05','B06','B07','B08','B11','B12']

FILL_VALUE = 2**16-1

ds = stackstac.stack(
                item,
                assets = boi,
                resolution=20,
                dtype="uint16",
                fill_value=FILL_VALUE,
                bounds_latlon=aoi_ltdb,
                    )
ds2 = stackstac.stack(
                item2,
                assets = boi,
                resolution=20,
                dtype="uint16",
                fill_value=FILL_VALUE,
                bounds_latlon=aoi_ltdb,
                    )

nir = ds.sel(band="B08").astype('float')
sir = ds.sel(band="B12").astype('float')

nbr = (nir-sir)/(nir+sir)
nbr_n = nbr.to_numpy()

nir2 = ds2.sel(band="B08").astype('float')
sir2 = ds2.sel(band="B12").astype('float')


nbr2 = (nir2-sir2)/(nir2+sir2)
nbr2_n = nbr2.to_numpy()

dNBR = nbr_n - nbr2_n

dNBR_n = dNBR[0]

"""
plt.imshow(dNBR_n.squeeze(), cmap="RdYlGn_r")
plt.colorbar(cmap="RdYlGn_r")
plt.title("dNBR (June 2022-August 2022)")
plt.show()
"""

"""
cmap = matplotlib.colors.ListedColormap(['olivedrab','lightgreen','green','red'])
cmap.set_over('purple')
cmap.set_under('white')
bounds = [-0.5, -0.25, -0.1, 0.1, 1.3]        
norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)  

fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'xticks': [], 'yticks': []})
cax = ax.imshow(dNBR_n, cmap=cmap, norm = norm)
plt.title('Burn Recovery Map (dNBR)')
cbar = fig.colorbar(cax, ax=ax, fraction=0.035, pad=0.04, ticks=[-0.375, -0.175, 0, 0.6])
cbar.ax.set_yticklabels(['High restored', 'Low restored', 'Unburned', 'Burnt'])
plt.show()
"""

l,c = dNBR_n.shape
burnt = dNBR_n[-0.1>dNBR_n]

print(f"Recovered area: {len(burnt)*400/10000}ha")

"""
burnt = dNBR_n[dNBR_n > 0.1]
print(f"Mean: {np.mean(burnt)}")
print(f"STD: {np.std(burnt)}")
print(f"Median: {np.median(burnt)}")

plt.hist(burnt)
plt.axvline(np.median(burnt), color='k', linestyle='dotted', linewidth=1)
plt.axvline(burnt.mean(), color='k', linestyle='solid', linewidth=1)
plt.axvline(burnt.mean()-burnt.std(), color='k', linestyle='dashed', linewidth=1)
plt.axvline(burnt.mean()+burnt.std(), color='k', linestyle='dashed', linewidth=1)
plt.xlabel("NBR")
plt.ylabel("Count")
plt.title("Statistical distribution of NBR index on the burnt area")
plt.show()


print(f"Burnt area: {c*400}m²")
"""
"""
##Landcover
import planetary_computer
import pystac_client

# Open the Planetary Computer STAC API
catalog = pystac_client.Client.open(
    "https://planetarycomputer.microsoft.com/api/stac/v1/",
    modifier=planetary_computer.sign_inplace,
)
collection = catalog.get_collection("esa-cci-lc")
collection

latitude = 39.50
longitude = -98.35

search = catalog.search(collections=collection, bbox=aoi_landiras, datetime="2012-06")
items = list(search.items())

import odc.stac
ds = odc.stac.load(items, bbox=aoi_landiras)
landcovers = ds["lccs_class"].to_numpy()[0]
l,c = landcovers.shape

b_landcovers = [landcovers[x][y] for x in range(l) for y in range(c) if -0.1>dNBR_n[x][y]]        

unique, counts = np.unique(b_landcovers, return_counts=True)
print(dict(zip(unique, counts)))

percentage = np.array([element/np.sum(counts) for element in counts])
print(dict(zip(unique, percentage)))
"""