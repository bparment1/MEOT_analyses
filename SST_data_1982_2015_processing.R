#################################    EOT AND PCA ROTATION  #######################################
########################### SPACE-TIME VARIABILITY  #############################################
#This script examines the 1982-2016 dataset downloaded from NOAA:
#http://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.html
#
#The goal is to run EOT, MEOT and PCA on the updated SST dataset using the IDRISI and the R package "remote".
#
#AUTHOR: Benoit Parmentier                                                                       #
#DATE CREATED:07/07/2016 
#DATE MODIFIED: 07/07/2016
#
#PROJECT: MEOT/EOT climate variability extraction
#
##################################################################################################
#
###Loading r library and packages

library(raster)                             # loading the raster package
library(gtools)                             # loading ...
library(sp)                                 # spatial objects in R
library(gplots)                             # 
library(rgdal)                              # gdal driver for R
library(RColorBrewer)                       # color scheme, palettes used for plotting
library(gdata)                              # read different format (including .xlsx)
library(plotrix)                            # plot options and functions including plotCI
library(rasterVis)                          # raster visualization
library(gridExtra)
library(latticeExtra)
library(colorRamps)                         # contains matlab.like palette
library(lsr)
library(psych)
library(GPArotation)
library(zoo)
library(xts)
library(remote)                            # EOT implementation in R/cpp

#################################################
###### Functions  used in the script  ##########

create_dir_fun <- function(outDir,out_suffix){
  #if out_suffix is not null then append out_suffix string
  if(!is.null(out_suffix)){
    out_name <- paste("output_",out_suffix,sep="")
    outDir <- file.path(outDir,out_name)
  }
  #create if does not exists
  if(!file.exists(outDir)){
    dir.create(outDir)
  }
  return(outDir)
}

load_obj <- function(f){
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

#infile1_function <- file.path("/home/bparmentier/Google Drive/Papers_writing_MEOT/R_scripts/",
#                              "EOT_PCA_rotation_functions_09152015.R")
#source(infile1_function)

#############################################
######## Parameters and arguments  ########

out_suffix <- "sst_07072016"
inDir <- "/home/bparmentier/Google Drive/Papers_writing_MEOT/MEOT_paper/SST_data_update_1982_2015"
setwd(inDir)

#outDir <- "/Users/benoitparmentier/Dropbox/Data/Dissertation_paper2_04142012"
#outDir <- "/home/bparmentier/Google Drive/Papers_writing_MEOT/EOT_paper"
outDir <- inDir

create_outDir_param = TRUE

#Create output directory

if(create_outDir_param==TRUE){  
  outDir <- create_dir_fun(outDir,out_suffix)
  setwd(outDir)
}else{
  setwd(outDir) #use previoulsy defined directory
}

########################################################
##############  Start of th script  ##############

### PART 0: READ IN DATASETS RELATED TO TELECONNECTION AND PREVIOUS LOADINGS


system("gdalinfo '/home/bparmentier/Google\ Drive/Papers_writing_MEOT/MEOT_paper/SST_data_update_1982_2015/sst.mnmean.nc'")
r_filename <- ("/home/bparmentier/Google Drive/Papers_writing_MEOT/MEOT_paper/SST_data_update_1982_2015/sst.mnmean.nc")

r_sst <- brick(r_filename)

system("gdal_translate -of GTiff '/home/bparmentier/Google\ Drive/Papers_writing_MEOT/MEOT_paper/SST_data_update_1982_2015/sst.mnmean.nc' test.tif")
r_sst <- brick("test.tif")
NAvalue(r_sst) <- 
plot(r_sst,y=1)

#http://geog.uoregon.edu/bartlein/courses/geog607/Rmd/netCDF_01.htm

#ncFile <- open.ncdf(nc)
#print(ncFile) 
#vars <- names(ncFile$var)[1:12] # I'll try to use these variable names later to make a list of bricks
