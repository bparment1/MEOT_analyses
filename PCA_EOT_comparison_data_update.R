#################################    EOT AND PCA ROTATION  #######################################
########################### SPACE-TIME VARIABILITY  #############################################
#This script examines the 1982-2016 dataset downloaded from NOAA:
#http://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.html
#
#The goal is to run EOT, MEOT and PCA on the updated SST dataset using the IDRISI and the R package "remote".
#
#AUTHOR: Benoit Parmentier                                                                       #
#DATE CREATED:07/11/2016 
#DATE MODIFIED: 07/16/2016
#
#PROJECT: MEOT/EOT climate variability extraction
#
##################################################################################################
#
###Loading r library and packages

library(raster)                            # loading the raster package
library(gtools)                            # loading ...
library(sp)                                # spatial objects in R
library(gplots)                            # 
library(rgdal)                             # gdal driver for R
library(RColorBrewer)                      # color scheme, palettes used for plotting
library(gdata)                             # read different format (including .xlsx)
library(plotrix)                           # plot options and functions including plotCI
library(rasterVis)                         # raster visualization
library(gridExtra)                         # graphic package
library(latticeExtra)                      # graphic package
library(colorRamps)                        # contains matlab.like palette
library(lsr)                               #
library(psych)                             # PCA
library(GPArotation)                       # PCA rotation
library(zoo)                               # Time series object and functions
library(xts)                               # Time series object and functions
library(remote)                            # EOT implementation in R/cpp
library(XML)                               # HTML funcitons
library(readODS)                           # read open data spreadsheet format
library(plyr)                              # contains "rename","revalue" and other useful functions

#################################################
###### Functions  used in the script  ##########

infile1_function <- file.path("/home/bparmentier/Google Drive/Papers_writing_MEOT/R_scripts/",
                             "PCA_EOT_comparison_data_update_function_07152016.R")
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

out_suffix <- "eot_pca_1982_2015_anom_07152016"
create_out_dir_param=TRUE #PARAM8
num_cores <- 4 #PARAM 9

years_to_process <- 1982:2015
#start_date <- "2012-01-01" #PARAM 12
#end_date <- "2012-12-31" #PARAM 13 #should process by year!!!
var_name <- "sst" #PARAM 14, Name of variable of interest: bacteria measurement (DMR data)
scaling <- 1/0.0099999998

r_mask_filename <- "/home/bparmentier/Google Drive/Papers_writing_MEOT/EOT_paper/data/lsmasked_0_180.rst"

in_dir <- "/home/bparmentier/Google Drive/Papers_writing_MEOT/EOT_paper/EOT82_15_July12run/sst_msk_0_180_1982_2015_anom/components"
setwd(in_dir)

#out_dir <- "/Users/benoitparmentier/Dropbox/Data/Dissertation_paper2_04142012"
#out_dir <- "/home/bparmentier/Google Drive/Papers_writing_MEOT/EOT_paper"
out_dir <- in_dir

create_out_dir_param = TRUE

#Create output directory

if(create_out_dir_param==TRUE){  
  out_dir <- create_dir_fun(out_dir,out_suffix)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

SST_dir <- "SST_1982_2015"
eot_dir <- "/home/bparmentier/Google Drive/Papers_writing_MEOT/workdir_terrset_08282015/anom_sst_1982_2007/components"
#mask_fname <- "mask_rgf_1_1.tif"
#eot_fname1 <- "eot_std_s7_test__EOT_Center_Std.avl"
#out_suffix <-"_eot_pca_12272015"
C#RS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84 # CONST 2

lf_sst <- list.files(path=file.path(in_dir,SST_dir),pattern=".rst$",full.names=T)

pca_fname1 <- file.path(in_dir,"sst_msk_0_180_1982_2015_anom_PCA_Center_Std_S-Mode_DIF.ods") #ods file with loadings and variance
eot_fname1 <- file.path(in_dir,"eot_sst_msk_0_180_1982_2015_anom_EOT_Center_Std.ods")

indices_fname <- "/home/bparmentier/Google Drive/Papers_writing_MEOT/EOT_paper/data/indices_new.xls"

telind<-c("PNA","NAO","TNA","TSA","SAOD","MEI","PDO","AO","AAO","AMM","AMOsm","QBO")
#mode_list_MEOT<-c("MEOT1","MEOT3", "MEOT4","MEOT7","MEOT10","MEOT15","MEOT16")
#mode_list_MEOT<-paste("MEOT",1:25,sep="")
#c("MEOT1","MEOT3", "MEOT4","MEOT7","MEOT10","MEOT15","MEOT16")

#mode_list_PCA<-c("MSSA1","MSSA2","MSSA3","MSSA4","MSSA5","MSSA6")    
#mode_list_PCA<-paste("MSSA",1:15,sep="")

lag_window<-13

########################################################
##############  Start of th script  ##############

### PART 0: READ IN DATASETS RELATED TO TELECONNECTION AND PREVIOUS LOADINGS
#http://geog.uoregon.edu/bartlein/courses/geog607/Rmd/netCDF_01.htm


#scrping tables from html pages
#http://stackoverflow.com/questions/1395528/scraping-html-tables-into-r-data-frames-using-the-xml-package

### GET DATA: 1982-2015 EOT and PCA

pca_dat <- read_ods(pca_fname1,sheet="pca_loadings")
#test <- read_ods(eot_fname1,sheet="pca_loadings")
names(pca_dat)[1] <- "fnames"
names(pca_dat)[2:ncol(pca_dat)] <- paste0("pc_",1:20)
test <- subset(pca_dat,select=paste0("pc_",1:20))
test<-lapply(test,FUN=convert_to_numeric)
test<- do.call(cbind,test)
pca_dat <- as.data.frame(test)

eot_dat <- read_ods(eot_fname1,sheet="eot_loadings")
#test <- read_ods(eot_fname1,sheet="pca_loadings")
test <- subset(eot_dat,select=paste0("eot",1:20))
test<-lapply(test,FUN=convert_to_numeric)
test<- do.call(cbind,test)
eot_dat <- as.data.frame(test)

#Import results teleconneciton indices processed by Elsa 

dat<-read.xls(indices_fname, sheet="1982-2015")
names(dat)[1:2] <- c("year","month")
dat <- rename(dat, c("QBO_30_original"="QBO")) #from plyr package

dat_indices <- subset(dat,select=telind)

##Not elegant but works for now...
test<-lapply(dat_indices ,FUN=convert_to_numeric)
test<- do.call(cbind,test)
dat_indices  <- as.data.frame(test)

#Creating time series objects
d_ts<-ts(data=dat,start=c(1982,1), end=c(2015,12), frequency=12, names=names(dat))
colnames(d_ts)
d_z<-as.zoo(d_ts)  #### THIS IS THE TIME SERIES OBJECT USED LATER ON
d_xts<-as.xts(d_ts)

## find all 7th of the month between two dates, the last being a 7th.
st <- as.Date("1982-1-15")
en <- as.Date("2015-12-15")
dseq <- seq(st, en, by="month") #Creating monthly date sequence to create a time series from a data frame
d_z2<-zoo(dat,dseq)
time(d_z)  #no time stamp??
time(d_z2) 

### NOW CARRY OUT CORRELATION BETWEEN INDICES AND COMPONENTS

#cor_series_fun(ts1,ts2,fig=F,out_suffix)
#debug(cor_series_fun)
cor_pca_df <- cor_series_fun(ts1=pca_dat,ts2=dat_indices,fig=F,out_suffix)
write.table(cor_pca_df ,file=file.path(out_dir,paste0("cor_pca_df","_",out_suffix,".txt")),sep=",")

barplot(as.numeric(diag(as.matrix(cor_pca_df))),
        ylim=c(-1,1),
        names.arg=1:20) #Okay very similar results

cor_eot_df <- cor_series_fun(ts1=eot_dat,ts2=dat_indices,fig=F,out_suffix)
write.table(cor_eot_df ,file=file.path(out_dir,paste0("cor_eot_df","_",out_suffix,".txt")),sep=",")

barplot(as.numeric(diag(as.matrix(cor_eot_df))),
        ylim=c(-1,1),
        names.arg=1:20) #Okay very similar results

### generate barplot for each components + indices?
### Generate maps+ temporal loadings figures for each (pca and eot!!)
### Generate maps+ barplots?

########################## END OF SCRIPT ###############################################

