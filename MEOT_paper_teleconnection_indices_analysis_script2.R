####################################    MEOT-MSSA PAPER   #######################################
############################      SCRIPT 2- MEOT          #######################################
#This script carries out a MSSA with varimax rotation for a timeseries of 312 images.            #
#Note that spatial patterns from MEOT and MSSA components are not analyzed in this script       #                 
#AUTHOR: Benoit Parmentier                                                                      #
#DATE: 03/09/2013            
#Version: 4
#PROJECT: Clark Labs Climate predction- MEOT/MSSA paper                                         #
#################################################################################################

###Loading R library and packages                                                      
library(gtools)                              # loading some useful tools 
library(mgcv)                                # GAM package by Simon Wood
library(sp)                                  # Spatial pacakge with class definition by Bivand et al.
library(spdep)                               # Spatial pacakge with methods and spatial stat. by Bivand et al.
library(rgdal)                               # GDAL wrapper for R, spatial utilities
library(gstat)                               # Kriging and co-kriging by Pebesma et al.
library(fields)                              # NCAR Spatial Interpolation methods such as kriging, splines
library(raster)                              # Hijmans et al. package for raster processing
library(foreign)                             # Library for format exchange (e.g. dbf,spss,sas etc.)
library(gdata)                               # various tools with xls reading
library(xts)                                 # basic package for time series analysis
library(zoo)                                 # basic package for time series analysis
#library(forecast)                            # package containing ARIMA procedures
library(rasterVis)
library(psych)
### Parameters and argument

infile1<-"SAODI-01-1854_06-2011_test.asc"             #GHCN shapefile containing variables for modeling 2010                 
infile2<-"SAODI-01-1854_06-2011.csv"                     #List of 10 dates for the regression
#infile2<-"list_365_dates_04212012.txt"
infile3<-"MEOT_MSSA_Telcon_indices_08062012.xlsx"                        #LST dates name
infile3<-"MEOT_MSSA_Telcon_indices_12112012.xlsx"                        #LST dates name
infile4<-"mask_rgf_1_1.rst"

path<-"/Users/benoitparmentier/Dropbox/Data/MEOT_paper/MEOT12272012/MEOT_working_dir_10232012/MEOT_analysis_R_10102012"
in_path<-"/Users/benoitparmentier/Dropbox/Data/MEOT_paper/MEOT12272012/MEOT_working_dir_10232012/MEOT10232011/anom_sst_1982_2007/components"
in_path<-"/Users/benoitparmentier/Dropbox/Data/MEOT_paper/MEOT12272012/MEOT_working_dir_10232012/SST_1982_2007"
setwd(in_path)

#path<-"/Users/benoitparmentier/Documents/DATA/Benoit/Clark_University/Paper_writings/MSSA_BNP/"
#on MAC:
#path<-"/Users/benoitparmentier/Dropbox/Data/MEOT_paper/MEOT12272012/MEOT_working_dir_10232012/MEOT_analysis_R_10102012"
#path<-"/Users/benoitparmentier/Documents/DATA/Benoit/Clark_University/Paper_writings/MSSA_BNP/MEOT_analysis_R_10102012"
# on Atlas:
#path<-"/home/parmentier/Data/MEOT12272012/MEOT_working_dir_10232012/MEOT_analysis_R_10102012"
#path_raster<-"/Users/benoitparmentier/Dropbox/Data/MEOT_paper/MEOT12272012/MEOT_working_dir_10232012/MSSA_EEOT_04_29_09"
telind<-c("PNA","NAO","TNA","TSA","SAOD","MEI","PDO","AO","AAO","AMM","AMOsm","QBO")
mode_list_MEOT<-c("MEOT1","MEOT3", "MEOT4","MEOT7","MEOT10","MEOT15","MEOT16")
#mode_list_MEOT<-paste("MEOT",1:25,sep="")
#c("MEOT1","MEOT3", "MEOT4","MEOT7","MEOT10","MEOT15","MEOT16")

out_prefix<-"MEOT_paper_02102013"

lf<-mixedsort(list.files(pattern="ANOM_SST.*.rst"))
SST_s<-stack(lf)
SST1<-raster(lf[1])
lag_window<-13

mask_land<-raster(infile4)
mask_land_NA<-mask_land
mask_land_NA[mask_land_NA==0]<-NA
SST_rast<-mask(SST_s,mask_land_NA,filename="ANOM_SST_1982_2007.tif",overwrite=TRUE)
class(SST_rast)

#Need to add weighting scheme base on lat long...
#SST_sgdf<-as(SST_rast,"SpatialGridDataFrame")
SST_df<-as.data.frame(SST_sgdf) #drop the coordinates x and y
rm(SST_sgdf)
cor.matrix.1<-cor(SST_df,use="complete")  
#cor.matrix.1<-layerStats(SST_rast,"pearson")

pca.2 <- principal(r = cor.matrix.1, nfactors = 20, residuals = FALSE, rotate = "varimax")
pca_varimax<-principal(r=SST_df, nfactors = 20, residuals = FALSE, rotate = "varimax",scores=TRUE)
#components are stored in the pca_varimax$scores
plot(pca.2)

### TO DO ... format the data for MSSA create a function for the window of any size ...use lag???
#unction ot lage files
lag_nb<-13
#lf_l1<-lag(lf,k=lag_nb)

function(lag_window,list_lf){
  list_lf<-vector("list",lag_window)
  for (j in 1:lag_window){
    index1<-j
    index2<-length(lf)-(lag_window-j)
    list_lf[[j]]<- lf[index1:index2]
  }
  n_tl<-length(lf)-lag_window+1
  list_t<-vector("list",n_tl)
  list_tmp<-vector("list",lag_window)
  for (i in 1:n_tl){
    for (j in 1:lag_window){
      list_tmp[[j]]<-(list_lf[[j]][i])
    }
    list_t[[i]]<-list_tmp
  }
  list_obj<-c(list_lf,list_t)
  names(list_obj)<-c("list_lf","list_t")
  return(list_obj)
}

t1<-stack(list_t[[1]])
levelplot(t1)
df1<-as.data.frame(t1)
df1_v<-as.vector(t(df1))   #This transform the data frame into a one colum

t2<-stack(list_t[[2]])
levelplot(t2)
df2<-as.data.frame(t2)
df2_v<-as.vector(t(df2))   #This transform the data frame into a one colum
head(df2)
test<-rbind(df1[,1],df1[,2])
#MSSA1_l<-lag(d_z$MSSA1,k=lag_nb)

SST1
ncol_concat<-ncol(SST1)*lag_window
test<-raster(nrow=nrow(SST1),ncol=ncol_concat)
#MSSA3_l<-lag(d_z$MSSA3,k=lag_nb)
projection(test)<-NA
extent(test)<-c(1,ncol_concat+1,1,181)

j=lag_number
projection(t1_1)<-NA
j=1

t1_1<-subset(t1,1)
t1_2<-subset(t1,2)
for (j in 1:lag_window){
  extent(t1_1)<-c(1+(360*j),361+(360*j),1,181)
}

mosaic together
#change extent for each of layer

