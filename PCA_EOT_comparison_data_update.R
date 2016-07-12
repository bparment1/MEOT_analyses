#################################    EOT AND PCA ROTATION  #######################################
########################### SPACE-TIME VARIABILITY  #############################################
#This script examines the 1982-2016 dataset downloaded from NOAA:
#http://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.html
#
#The goal is to run EOT, MEOT and PCA on the updated SST dataset using the IDRISI and the R package "remote".
#
#AUTHOR: Benoit Parmentier                                                                       #
#DATE CREATED:07/11/2016 
#DATE MODIFIED: 07/11/2016
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
library(XML)                               # HTML funcitons

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

convert_to_numeric <-function(x){
  if(class(x)=="character"){
    x<-as.numeric(x)
  }
  ##
  if(class(x)=="factor"){
    x<-as.numeric(as.character(x))
  }
  return(x)
}

infile1_function <- file.path("/home/bparmentier/Google Drive/Papers_writing_MEOT/R_scripts/",
                             "EOT_PCA_rotation_functions_01092016.R")
source(infile1_function)

#############################################
######## Parameters and arguments  ########

CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84 # CONST 2
proj_str<- CRS_WGS84 #param 2
CRS_reg <- CRS_WGS84 # PARAM 3

file_format <- ".rst" #PARAM 4
NA_value <- -9999 #PARAM5
NA_value_SST <- 32767
NA_flag_val <- NA_value #PARAM6
out_suffix <-"NEST_prism_07112016" #output suffix for the files and ouptu folder #PARAM 7
create_out_dir_param=TRUE #PARAM8
num_cores <- 4 #PARAM 9

#station_data_fname <- file.path("/home/bparmentier/Google Drive/NEST/", "MHB_data_2006-2015.csv") #PARAM 11

#years_to_process <- 2003:2016
years_to_process <- 1982:2016
#start_date <- "2012-01-01" #PARAM 12
#end_date <- "2012-12-31" #PARAM 13 #should process by year!!!
var_name <- "sst" #PARAM 14, Name of variable of interest: bacteria measurement (DMR data)
scaling <- 1/0.0099999998

r_mask_filename <- "/home/bparmentier/Google Drive/Papers_writing_MEOT/MEOT_paper/SST_data_update_1982_2015/lsmask.nc"

out_suffix <- "sst_old_data_pca_07072016"
inDir <- "/home/bparmentier/Google Drive/Papers_writing_MEOT/000_EOT/EOT_MEOT/data/Old_data"
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


#scrping tables from html pages
#http://stackoverflow.com/questions/1395528/scraping-html-tables-into-r-data-frames-using-the-xml-package

### OLD DATA: 1982-2007
#make this a function
theurl <- "file:///home/bparmentier/Google%20Drive/Papers_writing_MEOT/000_EOT/EOT_MEOT/data/Old_data/PCA_old_S-MODE_stn_cn_20comp.html"
#l_tables <- readHTMLTable(theurl)
l_tables <- readHTMLTable(theurl,header=T) #list of tables extracted from the html documents
n.rows <- unlist(lapply(l_tables, function(t) dim(t)[1]))
variance_df1 <- l_tables[[1]]
components_df1 <- l_tables[[2]]
rownames(components_df1) <- as.character(components_df1[,1])
components_df1 <- components_df1[,-1]

#reformat
test<-lapply(components_df1,FUN=convert_to_numeric)
test2<- do.call(cbind,test)
test2 <- as.data.frame(test2)
names(test2)
names(test2)<- paste0("cp_old_",1:20)
components_df1 <- test2
  
### OLD DATA: 1982-2015
## table
theurl <- "file:///home/bparmentier/Google%20Drive/Papers_writing_MEOT/000_EOT/EOT_MEOT/data/New_data/S-MODE1_newdata_1982_2007.html"
l_tables <- readHTMLTable(theurl,header=T) #list of tables extracted from the html documents
n.rows <- unlist(lapply(l_tables, function(t) dim(t)[1]))
variance_df2 <- l_tables[[1]]
components_df2 <- l_tables[[2]]
rownames(components_df2) <- as.character(components_df2[,1])
components_df2 <- components_df2[,-1]

#reformat
test<-lapply(components_df2,FUN=convert_to_numeric)
test2<- do.call(cbind,test)
test2 <- as.data.frame(test2)
names(test2)
names(test2)<- paste0("cp_new_",1:20)
components_df2 <- test2

###

#cor_series_fun(ts1,ts2,fig=F,out_suffix)
#debug(cor_series_fun)
test <- cor_series_fun(ts1=components_df1,ts2=components_df2,fig=F,out_suffix)

barplot(as.numeric(diag(as.matrix(test))),
        ylim=c(-1,1),
        names.arg=1:20) #Okay very similar results
   
#############



########################## END OF SCRIPT

#data(vdendool) #data of 36 cols and 14 rows! very small

## claculate 2 leading modes
#nh_modes <- eot(x = vdendool, y = NULL, n = 2,
#                standardised = FALSE,
#                verbose = TRUE)
#plot(nh_modes, y = 1, show.bp = TRUE)
#plot(nh_modes, y = 2, show.bp = TRUE)
