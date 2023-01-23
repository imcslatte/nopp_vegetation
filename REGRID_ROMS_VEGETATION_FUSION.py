
"""
REGRID_ROMS_VEGETATION_FUSION.py

Download VEGETATION FUSION data from file, maps the land use classifications to ROMS plant_* vaiables, and then regrids them to a given ROMS grid. 

usage:

1) edit relevant variables in the section labeled below.*
2) execute file using something like:
    % conda activate "RELEVANT_EVIRONEMENT"
    % python REGRID_ROMS_VEGETATION_FUSION.py

*- see README
Created by Elias Hunter, hunter@marine.rutgers.edu, 1/23/2023 
Code for mapping classifcation to plant parameters is from Chris Sherwood (csherwood@usgs.gov) and John Warner (jcwarner@usgs.gov)
 
"""

import xarray as xr
import numpy as np
import cartopy
import time
import rioxarray
import matplotlib.pyplot as plt
import cartopy 
import xesmf as xe
from osgeo import gdal
import shutil,netCDF4


print('VEGETATION REGRID STARTED')


########################################################################
#EDIT BETWEEN HERE
########################################################################
# INPUT GRID FILE 
file=r'/home/hunter/roms/NOPP/grids/michael_grd3.nc'
# INPUT Vegetation Classification FILE 
infile=r'/home/hunter/roms/NOPP/data/FPH_NOPP_Fusion/FPH_NOPP_Fusion.tif'
# INPUT Initialization FILE 
inifile=r'/home/hunter/roms/NOPP/grids/michael_grd3_ini2.nc'
# Output Initialization FILE 
newinifile=r'/home/hunter/roms/NOPP/grids/michael_grd3_ini3.nc'
# Temporary subset geotiff file  
tmptiff = r'/home/hunter/roms/NOPP/data/FPH_NOPP_Fusion/tmp.tif' 

#Padding for geotiff subset. 
dx=0.01
dy=0.01

########################################################################
#AND HERE
########################################################################


print('Getting grid')
proj = cartopy.crs.Mercator(central_longitude=-98)
pc = cartopy.crs.PlateCarree()
grd1 = xr.open_dataset(file,chunks={'eta_rho':900,'xi_rho':600}) 
grd1=grd1.set_coords(('lat_rho','lon_rho'))

glon=np.concatenate((grd1.lon_rho[0,:].values,grd1.lon_rho[:,-1].values,np.flip(grd1.lon_rho[-1,:].values),np.flip(grd1.lon_rho[:,0].values)))  
glat=np.concatenate((grd1.lat_rho[0,:].values,grd1.lat_rho[:,-1].values,np.flip(grd1.lat_rho[-1,:].values),np.flip(grd1.lat_rho[:,0].values)))  

#Subset limits
ulx = np.min(glon)-dx
uly = np.max(glat)+dy
lrx = np.max(glon)+dx
lry = np.min(glat)-dy
#Subset geotiff using gdal Translate
ds = gdal.Open(infile)
ds = gdal.Translate(tmptiff, ds, projWin = [ulx, uly, lrx, lry],projWinSRS="EPSG:4326")
ds = None


# Read in subsetted geotiff
img=rioxarray.open_rasterio(tmptiff,masked=True)
img=img.squeeze()
img = img.rio.reproject("EPSG:4326")


#regrid subsetted file using xesmf and a nearest neighbor interpolation
st = time.time()    
varint = xr.Dataset({"lat": (["eta_rho","xi_rho"], grd1.lat_rho.values), "lon": (["eta_rho","xi_rho"], grd1.lon_rho.values)})
varint=varint.set_coords(('lat','lon'))
regridder = xe.Regridder(img, varint, "nearest_s2d", locstream_out=False)
y=regridder(img, keep_attrs=True)
et = time.time()
# get the execution time
elapsed_time = et - st
print('REGRID execution time:', elapsed_time, 'seconds')



# Lookup up table to specify plant parameters from CCAP classifications. 
# lookup table for converting C-CAP landcover to a simplified classification


#CCLass: indices to tuples
iDescription = 0
iPlantHeight_m = 1
iPlantThick_m = 2
iPlantdens_m2 = 3
iPLatmDiam_m = 4

#   Id: (Name, plant_height(m), plant_thickness(m), plant_density(stems/m2), plant_diameter(m) )   
CClass = {
    0: ("NoData"    ,        0,   0,      0,   0 ), 
    1: ("open water",        0,   0,      0,   0), 
    2: ("sandy"     ,        0,   0,      0,   0),       
    3: ("grassy veg",        0.5, 0.005,  1.0, 0.03),
    4: ("woody veg" ,        2.0, 0.005, 10.0, 0.03),
    5: ("marshy veg",        0.5, 0.005,  1.0, 0.03),
    6: ("devel open",        0,   0,      0,   0),    
    7: ("devel structures", 10,   1,      1,   1)
    }

# CCAP landcover classes
# Tuple contains name and corresponding COAWST class. Those classes could be reassigned.
# Classes we would not expect (e.g., tundra) return 0 as a warning something might be amiss
CCAP = { 0: ( "Background" , 0 ),
         1: ( "Unclassified", 0 ),
         2: ( "Developed, High Intensity", 7 ),
         3: ( "Developed, Medium Intensity", 7 ),
         4: ( "Devleoped, Low Intensity", 7 ),
         5: ( "Developed, Open Space", 6 ),
         6: ( "Cultivated Crops", 3 ),
         7: ( "Pasture/Hay", 3 ),
         8: ( "Grassland/Herbaceous", 3),
         9: ( "Deciduous Forest", 4 ),
        10: ( "Evergreen Forest", 4 ),
        11: ( "Mixed Forest", 4 ),
        12: ( "Scrub/Shrub", 4 ),
        13: ( "Palustrine Forested Wetland", 5 ),
        14: ( "Palustrine Scrub/Shrub Wetland", 5 ),
        15: ( "Palustrine Emergent Wetland (Peristent)", 5 ),
        16: ( "Estuarine Forested Wetland", 4 ),
        17: ( "Estuarine Scrub/Shrub Wetland", 4, ),
        18: ( "Estuarine Emergent Wetland", 5 ),
        19: ( "Unconsolidated Shore", 2 ),
        20: ( "Barren Land", 2 ),
        21: ( "Open Water", 1 ),
        22: ( "Palustrine Aquatic Bed", 1 ),
        23: ( "Estuarine Aquatic Bed", 1 ),
        24: ( "Tundra", 0 ), 
        25: ( "Perennial Ice/Snow", 0 ),        
       }


s = np.shape(y)
print("Shape: ",s)
plantheight= np.nan*np.zeros_like(y)
plantthick= np.nan*np.zeros_like(y)
plantdens= np.nan*np.zeros_like(y)
plantdiam= np.nan*np.zeros_like(y)
for j in range (s[1]):
    for i in range (s[0]):
        plantheight[i,j] = CClass[CCAP[y.data[i,j]][1]][iPlantHeight_m]
        plantthick[i,j] = CClass[CCAP[y.data[i,j]][1]][iPlantThick_m]
        plantdens[i,j] = CClass[CCAP[y.data[i,j]][1]][iPlantdens_m2]
        plantdiam[i,j] = CClass[CCAP[y.data[i,j]][1]][iPLatmDiam_m]



#Create new initialization  file
shutil.copy(inifile,newinifile)


#wrtite plant paramters to new initilization file.
nc = netCDF4.Dataset(newinifile, "r+", format="NETCDF4")
PH=nc['plant_height']
PD=nc['plant_density']
PDI=nc['plant_diameter']
PT=nc['plant_thickness']


PH[:,:]=plantheight
PD[:,:]=plantdens
PDI[:,:]=plantdiam
PT[:,:]=plantthick


nc.close()

print('VEGETATION REGRID FINISHED')