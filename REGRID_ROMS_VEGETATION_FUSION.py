"""
REGRID_ROMS_VEGETATION_FUSION

REgrid NOPP/C-CAP Fusion Product to ROMS grid. 


"""


#Import Necessary Modules
import xarray as xr
#import dask as da
#import pandas as pd
import numpy as np
#import requests
#import cartopy
import time,os
#import geoviews as gv
#import rasterio as rs
import rioxarray
#import shapefile
#import matplotlib.pyplot as plt
#import cartopy.feature as cfeature
#from cartopy.feature import ShapelyFeature
#import cartopy.io.shapereader as shpreader
#from shapely import geometry
#import datashader
#import hvplot.xarray
#import cmocean.cm as cmo
#import xcmocean
import xesmf as xe
#from rioxarray import merge
#import geopandas
#from shapely.geometry import box
#Dask Initialize
#from dask.distributed import Client
#import dask.array as da
#client = Client()


# INPUT GRID FILE OPTIONS
gfile=r'/home/hunter/roms/NOPP/grids/michael_grd3.nc'
inifile=r'/home/hunter/roms/NOPP/DEM/SANIBEL_ROMS_FORECAST_11R_JUN2021_MAX3_CUDEM.nc'
vegfile=r'/home/hunter/roms/NOPP/vegetation/FPH_NOPP_Fusion/FPH_NOPP_Fusion_SUBSET2.tif'



#Extract Grid Information
grd = xr.open_dataset(gfile,chunks={'eta_rho':900,'xi_rho':600}) 
grd=grd.set_coords(('lat_rho','lon_rho'))
grdH=grd.h.load()
grdH = grdH.rename({"lon_rho": "lon", "lat_rho": "lat"})

glon=np.concatenate((grdH.lon[0,:].values,grdH.lon[:,-1].values,np.flip(grdH.lon[-1,:].values),np.flip(grdH.lon[:,0].values)))  
glat=np.concatenate((grdH.lat[0,:].values,grdH.lat[:,-1].values,np.flip(grdH.lat[-1,:].values),np.flip(grdH.lat[:,0].values)))  
llbounds=d = np.column_stack((glon,glat))

# Get SUBSETTED CUDEM 9th arc-sec data 
print('Reading and reprojecting data')
st = time.time()
img=rioxarray.open_rasterio(vegfile,masked=True)
img=img.squeeze()
img = img.rio.reproject("EPSG:4326")

varint = xr.Dataset({"lat": (["eta_rho","xi_rho"], grd.lat_rho.values), "lon": (["eta_rho","xi_rho"], grd.lon_rho.values)})
varint=varint.set_coords(('lat','lon'))
print('Regridding')
regridder = xe.Regridder(img, varint, "nearest_s2d", locstream_out=False)
imgout=regridder(img, keep_attrs=True)
 # get the end time
et = time.time()
    # get the execution time
elapsed_time = et - st
print('Execution time:', elapsed_time, 'seconds')



#tmpH=newH['h'].values
#tmpH[tmpH<depthcut]=depthcut
#wrtiting data put. 
#nc = netCDF4.Dataset(newgfile, "r+", format="NETCDF4")
#h=nc['h']
#h[:,:]=tmpH
#nc.close()

