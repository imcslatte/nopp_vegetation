
"""
REGRID_ROMS_VEGETATION_FUSION_v2.py

Download VEGETATION FUSION data from file, maps the land use classifications to ROMS plant_* vaiables, and then regrids them to a given ROMS grid. 

usage:

1) edit relevant variables in the section labeled below.*
2) execute file using something like:
    % conda activate "RELEVANT_EVIRONEMENT"
    % python REGRID_ROMS_VEGETATION_FUSION.py

*- see README
Created by Elias Hunter, hunter@marine.rutgers.edu, 2/13/2023 
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
file=r'http://icoast.rc.ufl.edu/thredds/dodsC/L2_grd/SNBL/SANIBEL_ROMS_FORECAST_11R_JUN2021_MAX3.nc'
# INPUT Vegetation Classification FILE 

infile=[r'https://cmgp-sfm-public-read-bucket.s3.us-west-2.amazonaws.com/NOPP/WPFL_NOPP_Fusion_COG.tif',
        r'https://cmgp-sfm-public-read-bucket.s3.us-west-2.amazonaws.com/NOPP/FLEV_NOPP_Fusion_COG.tif']
# INPUT Initialization FILE 
inifile=r'/home/hunter/roms/NOPP/grids/L1_L2_SANIBEL_run4_ini.nc'
# Output Initialization FILE 
newinifile=r'/home/hunter/roms/NOPP/grids/L1_L2_SANIBEL_run4_ini_VEG.nc'
# Temporary subset geotiff file  
tmptiff = r'/home/hunter/roms/NOPP/data/tmp/tmp'

#Padding for geotiff subset. 
dx=0.01
dy=0.01

#
dxy=50 #meters
dlat2m=111e3
dlat=dxy/dlat2m 

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
st = time.time()   
tmpfiles=[]
for ind,fi in enumerate(infile):
    print(fi)
    ofile=f'{tmptiff}_{ind}.tif'
    print(ofile)
    ds = gdal.Open(fi)
    ds = gdal.Translate(ofile, ds, projWin = [ulx, uly, lrx, lry],projWinSRS="EPSG:4326")
    gdal.Warp(ofile,ofile,dstSRS='EPSG:4326')
    ds = None
    tmpfiles.append(ofile)
et = time.time()
elapsed_time = et - st
print('Subset execution time:', elapsed_time, 'seconds')

# Read in subsetted geotiff
st = time.time()    
imgs=[];
for ifi in tmpfiles:
    print(ifi)
    img=rioxarray.open_rasterio(ifi)
    img=img.squeeze()
    img = img.rename({"x": "lon", "y": "lat"})
    imgs.append(img)


et = time.time()
# get the execution time
elapsed_time = et - st

#regrid subsetted file using xesmf and a nearest neighbor interpolation
st = time.time()    
dims=grd1.sizes
tlats= grd1['lat_rho'].values
tlons= grd1['lon_rho'].values
vID=np.zeros([dims['eta_rho'],dims['xi_rho']])-9.0

for img in imgs:
    for i in range(0,dims['eta_rho']-1):
        for j in range(0,dims['xi_rho']-1):
     #       st=time.time()
            dlon=dlat*np.cos(np.pi*tlats[i,j]/180.0)
            tmp=img.sel(lat=slice(tlats[i,j]+dlat,tlats[i,j]-dlat),lon=slice(tlons[i,j]-dlon,tlons[i,j]+dlon))
 
            if tmp.sizes['lon']==0 or tmp.sizes['lat']==0:
                continue
            tmp.load()
  #          tmp.where(tmp==128,0.0)   
            tmpvid=np.bincount(tmp.values.flatten()).argmax()
            if  vID[i,j]>0.0:
                continue
            
            tmpvid=np.bincount(tmp.values.flatten()).argmax()
            vID[i,j]=tmpvid
            

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
        128: ( "Unclassified", 0 ),
        -9: ( "Unclassified", 0 ),       
       }


s = np.shape(vID)
print("Shape: ",s)
plantheight= np.nan*np.zeros_like(vID)
plantthick= np.nan*np.zeros_like(vID)
plantdens= np.nan*np.zeros_like(vID)
plantdiam= np.nan*np.zeros_like(vID)

for j in range (s[1]):
    for i in range (s[0]):
        plantheight[i,j] = CClass[CCAP[vID.data[i,j]][1]][iPlantHeight_m]
        plantthick[i,j] = CClass[CCAP[vID.data[i,j]][1]][iPlantThick_m]
        plantdens[i,j] = CClass[CCAP[vID.data[i,j]][1]][iPlantdens_m2]
        plantdiam[i,j] = CClass[CCAP[vID.data[i,j]][1]][iPLatmDiam_m]
        
        
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