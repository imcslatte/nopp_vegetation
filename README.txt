README
Usage notes for REGRID_ROMS_VEGETATION_FUSION.py. This script reads and existing vegetation file (dowloaded from https://drive.google.com/drive/u/0/folders/1auf2adLJvvG7ZI57UwmfaSn5EGuJMYIy), then regrids the land use classification data to a supplied ROMS grid. The clasifcation is then mapped to the plant_density, plant_height, plant_diameter 


1) Download and uzip the relevant land use data from the google drive. There is currently no automated way to do this, as the data exists on a google drive.  

2) Edit the file REGRID_ROMS_VEGETATION_FUSION.py as needed. Specifically edit the variables in the following code block:


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


file: input grid file. This must be a local file. 
infile: INPUT Vegetation Classification FILE. Geotiff format 
inifile: ROMS initialization file containg original plant_* variables
newinifile: New ROMS initialization file with updated plant_* variables
dx and dy are padding scales used to define the subset footprint. 

3) Activate relevant conda environment. Note this must be run on either linux or maxOS as it uses xesmf, which is not available on Windows.    

  % activate conda VEGATATION_ENVIRON
  
  necessary conda packages:
	xarray
	numpy
	cartopy
	time
	rioxarray
	matplotlib
	cartopy 
	xesmf
	osgeo
	shutil
	netCDF4

 4) run  REGRID_ROMS_VEGETATION_FUSION.py
 
    % python  REGRID_ROMS_VEGETATION_FUSION.py 
    
    or
    
    % nohup python  REGRID_ROMS_VEGETATION_FUSION.py > log.dat & 
    
    
5) Thius will generate a new initilization file, newinifile  
 
 author: Elias Hunter, hunter @ marine.rutgers.edu 1/23/2023.  
 Chris Sherwood (csherwood@usgs.gov) and John Warner (jcwarner@usgs.gov) provided code for mapping the land use classifications to plant parameters. 