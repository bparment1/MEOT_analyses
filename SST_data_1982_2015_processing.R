#################################    EOT AND PCA ROTATION  #######################################
########################### SPACE-TIME VARIABILITY  #############################################
#This script examines the 1982-2016 dataset downloaded from NOAA:
#http://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.html
#
#The goal is to run EOT, MEOT and PCA on the updated SST dataset using the IDRISI and the R package "remote".
#
#AUTHOR: Benoit Parmentier                                                                       #
#DATE CREATED:07/07/2016 
#DATE MODIFIED: 07/08/2016
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

CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84 # CONST 2
proj_str<- CRS_WGS84 #param 2
CRS_reg <- CRS_WGS84 # PARAM 3

file_format <- ".rst" #PARAM 4
NA_value <- -9999 #PARAM5
NA_value_SST <- 32767
NA_flag_val <- NA_value #PARAM6
out_suffix <-"NEST_prism_07032016" #output suffix for the files and ouptu folder #PARAM 7
create_out_dir_param=TRUE #PARAM8
num_cores <- 4 #PARAM 9

rainfall_dir <- "/home/bparmentier/Google Drive/NEST_Data" #PARAM 10
station_data_fname <- file.path("/home/bparmentier/Google Drive/NEST_Data/", "WQ_TECS_Q.txt") #PARAM 11,DMR
#station_data_fname <- file.path("/home/bparmentier/Google Drive/NEST/", "MHB_data_2006-2015.csv") #PARAM 11

#years_to_process <- 2003:2016
years_to_process <- 1982:2016
#start_date <- "2012-01-01" #PARAM 12
#end_date <- "2012-12-31" #PARAM 13 #should process by year!!!
var_name <- "sst" #PARAM 14, Name of variable of interest: bacteria measurement (DMR data)
scaling <- 1/0.0099999998

r_mask_filename <- "/home/bparmentier/Google Drive/Papers_writing_MEOT/MEOT_paper/SST_data_update_1982_2015/lsmask.nc"
  
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
#http://geog.uoregon.edu/bartlein/courses/geog607/Rmd/netCDF_01.htm


system("gdalinfo '/home/bparmentier/Google\ Drive/Papers_writing_MEOT/MEOT_paper/SST_data_update_1982_2015/sst.mnmean.nc'")
r_filename <- ("/home/bparmentier/Google Drive/Papers_writing_MEOT/MEOT_paper/SST_data_update_1982_2015/sst.mnmean.nc")

system("gdal_translate -of GTiff '/home/bparmentier/Google\ Drive/Papers_writing_MEOT/MEOT_paper/SST_data_update_1982_2015/sst.mnmean.nc' test.tif")
r_sst <- brick("test.tif")
r_sst <- subset(r_sst,1:408)

r_sst_name <- unlist(lapply(1982:2015,FUN=function(year_val){paste0(year_val,"_",1:12)}))
names(r_sst) <- r_sst_name

r_sst <- setMinMax(r_sst)

### Assign dates and names
#NAvalue(r_sst)
#NAvalue(r_sst) <- NA_value_SST

#r_mask <- raster(r_mask_filename)
system("gdal_translate -of GTiff '/home/bparmentier/Google\ Drive/Papers_writing_MEOT/MEOT_paper/SST_data_update_1982_2015/lsmask.nc' lsmask.tif")
r_mask <- raster("lsmask.tif")
NAvalue(r_mask) <- 0
plot(r_mask)

#filename(r_stt)
#raster_name <- file.path(out_dir_str,paste(raster_name_tmp,"_masked.tif",sep=""))
raster_name <- file.path(outDir,paste("sst_1982_2015","_masked.tif",sep=""))
r_sst_m <- mask(r_sst,r_mask,filename=raster_name,overwrite=TRUE)

plot(r_sst_m,y=c(1,7))
filename(r_sst_m)
inMemory(r_sst_m)

### rescale ##
r_sst_m <- r_sst_m*1/(scaling)
raster_name <- file.path(outDir,paste("sst_1982_2015","_masked_scaled.tif",sep=""))
writeRaster(r_sst_m,filename=raster_name,overwrite=T)

#writeRaster(r_sst_m,filename=)

##Quick movie

for(i in 1:24){
  plot(r_sst_m,y=i)
}

### change coordinate systems??
if(change_coordinates==T){
  #
  #
  test <- r_sst_m
  extent_raster <- c(-180,180,-90,90)
  r <- setExtent(test, extent_raster, keepres=TRUE)
  
}
test2 <- shift(test,x=-180)
#shift to -180
r_test <- raster()
## can use deseaon ater for anomalies 
names(r_sst_m)<- r_sst_name
#test<- subset(r_sst_m,1:12)
rf <- writeRaster(r_sst_m, filename="sst.rst", format="IDRISI", overwrite=TRUE,bylayer=T,suffix=r_sst_name)
#rf <- writeRaster(test, filename="sst.rst", format="IDRISI", overwrite=TRUE,bylayer=T,suffix=r_sst_name[1:12])


########################## END OF SCRIPT #######################

### This looks like the bug was fixed on 07/08/2015
#data(vdendool)
## claculate 2 leading modes
#nh_modes <- eot(x = vdendool, y = NULL, n = 2,
#                standardised = FALSE,
#                verbose = TRUE)
#plot(nh_modes, y = 1, show.bp = TRUE)
#plot(nh_modes, y = 2, show.bp = TRUE)
