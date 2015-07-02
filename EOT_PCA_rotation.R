#################################    EOT AND PCA ROTATION  #######################################
########################### SPACE-TIME VARIABILITY  #############################################
#This script analyzes  SST data.
#The goal is to compare EOT and PCA with rotation.
#
#AUTHOR: Benoit Parmentier                                                                       #
#DATE CREATED: 07/02/2015 
#DATE MODIFIED: 07/02/2015
#
#PROJECT: Land transitions from Remote Sensing
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
library(colorRamps)                         # contains matlab.like palette
library(lsr)

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

cv_test_fun <- function(x,y) {
  #modified from:
  #http://www.r-bloggers.com/example-8-39-calculating-cramers-v/
  chisq_val <- chisq.test(x, y, correct=FALSE)
  CV = sqrt(chisq_val$statistic /
              (length(x) * (min(length(unique(x)),length(unique(y))) - 1)))
  print.noquote("Cramér V / Phi:")
  CV <- as.numeric(CV)
  cv_obj<- list(CV=CV,chisq_val=chisq_val)
  return(cv_obj)
}

#############################################
######## Parameters and arguments  ########

inDir <- "J:/Benoit/Data/MEOT_analyses_02202015"
SST_dir <- "SST_1982_2007"
mask_fname <- "mask_rgf_1_1.tif"
out_suffix <-"_eot_pca_07022015"

lf_sst <- list.files(path=file.path(in_dir,SST_dir),pattern=".rst$",full.names=T)


#proj_str<-"+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

infile1<-"SAODI-01-1854_06-2011_test.asc"             #GHCN shapefile containing variables for modeling 2010                 
infile2<-"SAODI-01-1854_06-2011.csv"                     #List of 10 dates for the regression
infile3<-"MEOT_MSSA_Telcon_indices_12112012.xlsx"                        #LST dates name
#infile4<-"mask_rgf_1_1.rst"
infile_pca <-"pca_IDRISI_03302013_T-Mode_DIF.xlsx"
infile_s_mode_pca  <- "pca_IDRISI_04152013_S-Mode_DIF.xlsx"
npc<-20 # number of pca to produce...put this at the beginning of the script...
lag_window<-13
telind<-c("PNA","NAO","TNA","TSA","SAOD","MEI","PDO","AO","AAO","AMM","AMOsm","QBO")
mode_list_MEOT<-c("MEOT1","MEOT3", "MEOT4","MEOT7","MEOT10","MEOT15","MEOT16")
mode_list_PCA<-c("MSSA1","MSSA2","MSSA3","MSSA4","MSSA5","MSSA6")    
#mode_list_PCA<-paste("MSSA",1:15,sep="")
out_prefix<-"pca_test_04152013"

setwd(inDir)

#outDir <- "/Users/benoitparmentier/Dropbox/Data/Dissertation_paper2_04142012"
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

### STEP 0: READ IN DATASETS RELATED TO TELECONNECTION AND PREVIOUS LOADINGS

SAODI<-read.table(infile2,sep=",", header=TRUE)
#Prepare data to write out in a textfile
s_SAODI2 <- subset(SAODI, year>1981 & year<2008) #Subset rows that correspond to the conditions
s_SAODI2<- subset(s_SAODI2, select=-year) #Remove the column year
in_SAODI2<-as.vector(t(s_SAODI2))   #This transform the data frame into a one colum
write.table(in_SAODI2,file=paste("SAOD_index_1981_2007",out_prefix,".txt",sep=""),sep=",")

r_sst <- stack(lf_sst)
r_mask <- raster(file.path(in_dir,mask_fname))
r_mask_NA <- r_mask
NAvalue(r_mask_NA) <- 0
#NAvalues()
levelplot(r_sst,layer=1)
levelplot(r_mask,layer=1)

r_sst_m <- mask(r_sst,r_mask_NA,filename=paste("sst_",out_suffix,".tif",sep=""),overwrite=T)
levelplot(r_sst_m,layer=1)

#Multiply layer by weight

lat_coord<-coordinates(r_mask_NA)[,2]
w_rast<- r_mask_NA
values(w_rast)<-cos(lat_coord*pi/180) #use area in raster package for weight of a cell??
w_rast_m<-mask(w_rast,r_mask_NA)
levelplot(w_rast_m)

r_sst_m <- r_sst_m*w_rast #Maybe use overlay to avoid putting in memory?
writeRaster(r_sst_m,filename="ANOM_SST_1982_2007_weighted.tif",overwrite=TRUE)
rm(r_sst)
r_sst <-brick("ANOM_SST_1982_2007_weighted.tif")

##############################################################
#### STEP 2: get data ready for PCA by lagging data frame...

SST_sgdf<-as(r_sst,"SpatialGridDataFrame")
#SST_df<-as.data.frame(SST_rast) #drop the coordinates x and y
SST_df<-as.data.frame(SST_sgdf)
write.table(SST_df,"sst_f_weighted.txt",sep=",")
rm(SST_sgdf)

#nv<-312
SST_xy<-SST_df[,c("s1","s2")]


### DO PCA WITH SST_df: T-mode cross product from a standardized dataset...(i.e. correlation matrix)
