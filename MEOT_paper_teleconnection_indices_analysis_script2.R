####################################    MEOT-MSSA PAPER   #######################################
############################      SCRIPT 2- MEOT          #######################################
#This script carries out a MSSA with varimax rotation for a timeseries of 312 SST images.       #
#Note that spatial patterns from MEOT and MSSA components are not analyzed in this script       #                 
#AUTHOR: Benoit Parmentier                                                                      #
#DATE: 03/13/2013            
#Version: 2
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
library(GPArotation)
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

out_prefix<-"MEOT_paper_03142013"

lf<-mixedsort(list.files(pattern="ANOM_SST.*.rst"))
SST_s<-stack(lf)
SST1<-raster(lf[1])
lag_window<-13

mask_land<-raster(infile4)
mask_land_NA<-mask_land
mask_land_NA[mask_land_NA==0]<-NA
SST_rast<-mask(SST_s,mask_land_NA,filename="ANOM_SST_1982_2007.tif",overwrite=TRUE)
class(SST_rast)

#Multiply layer by weight

lat_coord<-coordinates(SST_rast)[,2]
w_rast<-mask_land
values(w_rast)<-cos(lat_coord) #use area in raster package for weight of a cell??
w_rast<-mask(w_rast,mask_land_NA)

SST_rast<-SST_rast*w_rast #Maybe use overlay to avoid putting in memory?
writeRaster(SST_rast,filename="ANOM_SST_1982_2007_weighted.tif",overwrite=TRUE)
rm(SST_rast)
SST_rast<-brick("ANOM_SST_1982_2007_weighted.tif")

SST1_m<-subset(SST_rast,1)
plot(stack(SST1,SST1_m))
rm(SST1,SST1_m)
#Need to add weighting scheme base on lat long...
SST_sgdf<-as(SST_rast,"SpatialGridDataFrame")
#SST_df<-as.data.frame(SST_rast) #drop the coordinates x and y
SST_df<-as.data.frame(SST_sgdf)
rm(SST_sgdf)

### TO DO ... format the data for MSSA create a function for the window of any size ...use lag???
#unction ot lage files

lag_grouping<-function(lag_window,list_lf){
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
  list_obj<-list(list_lf,list_t)
  names(list_obj)<-c("list_lf","list_t")
  return(list_obj)
}

nt<- length(lf)-lag_window+1 #time dimension
nl<- lag_window #lag dimenstion 
ns<- ncell(SST1)
dimension<-list(nt,nl,ns)
ns_NA<-nrow(SST_df)
lag_list_obj<-lag_grouping(lag_window,lf) 
#list_lf contains gorupings of files for every lag (13 groups if lag_window is 13)
#list_t contains gorupings of files for every lag (13 groups if lag_window is 13)
list_lf<-lag_list_obj$list_lf
list_lag_df<-vector("list",length(lag_window))
names(SST_df)<-c(lf,"s1","s2")
new_names<-paste("t",1:nt,sep="_")
for (j in 1:lag_window){
  list_lag<-list_lf[[j]]
  lag_names<-c(list_lag,"s1","s2")
  list_lag_df[[j]]<-SST_df[lag_names]
  names(list_lag_df[[j]])<-c(new_names,"s1","s2")
  list_lag_df[[j]]$lag<-rep(j,ns_NA)
  #standardize in time? the rows...
}

SST_lag_data_df<-do.call(rbind,list_lag_df)
rm(list_lag_df)
save(SST_lag_data_df,file= paste("SST_lag_data_df_",out_prefix,".RData",sep=""))
#scale on the on the row to standardize pixel value in time?

#load("~/Dropbox/Data/MEOT_paper/MEOT12272012/MEOT_working_dir_10232012/SST_1982_2007/SST_lag_data_df_MEOT_paper_03122013.RData")

X<-as.matrix(SST_lag_data_df[,1:nt])
Xt<-t(X)
Xt<-scale(Xt)
X<-t(Xt)
cX<-Xt%*%X

pca_varimax<-principal(r=cX, nfactors = 20, residuals = FALSE, rotate = "varimax")
Xpc <-predict(pca_varimax,X)

#pca.2 <- principal(r = cor.matrix.1, nfactors = 20, residuals = FALSE, rotate = "varimax")
#pca_varimax<-principal(r=SST_df, nfactors = 20, residuals = FALSE, rotate = "varimax",scores=TRUE)
#pca_varimax<-principal(r=SST_lag_data_df[,1:nt], nfactors = 20, residuals = FALSE, rotate = "varimax",scores=TRUE)

save(pca_varimax,file= paste("pca_varimax_SST_lag_data_df_",out_prefix,".RData",sep=""))
save(X,file= paste("SST_lag_data_matrix_nt_",out_prefix,".RData",sep=""))
save(Xpc,file= paste("SST_lag_data_matrix_components_",out_prefix,".RData",sep=""))
save(X_pc_data,file= paste("SST_lag_data_components_",out_prefix,".RData",sep=""))

plot(pca_varimax)

### STEP 4: get lag images for MSSA/PCA

dim(Xpc)
SST_xy_lag<-SST_lag_data_df[,nt:ncol(SST_lag_data_df)]
X_pc_data<-as.data.frame(Xpc)
X_pc_data<-cbind(X_pc_data,SST_xy_lag)

coordinates(X_pc_data)<-X_pc_data[,c("s1","s2")]
#loop through k and j with j being lag index
tmp<-subset(X_pc_data,X_pc_data$lag==1 & X_pc_data[,k)
pc1_l1<-rasterize(tmp,SST1_m)??


### END OF SCRIPT...
#pc1_l2...
#################################
# nt<- length(lf)-lag_window+1 #time dimension
# nl<- lag_window #lag dimenstion 
# ns<- ncell(SST1)
# dimension<-list(nt,ntl,ns)
# list_lf<-lag_list_obj$list_lf
# list_lag_df<-vector("list",length(lag_window))
# 
# array_preparation_fun<-function(dimension,SST_rast){
# 
#   #fill in
#   nt<- length(lf)-lag_window+1 #time dimension
#   nl<- lag_window #lag dimenstion 
#   ns<- ncell(SST1)
#   
#   SST_array<-array(dim=c(ns,nt,nl))
#   
#   for (j in 1:lag_window){
#     list_lag<-list_lf[[j]]
#     #x<-SST_sgdf[list_lag]
#     SST_array[,,j]<-as.matrix(SST_df[list_lag])
#     #x<-as.matrix(x)
#     #SST_array[,,j]<-SST_sgdf[list_lag]
#   }
# }
# SST_df<-as.data.frame(SST_sgdf)
# 
# nt<- length(lf)-lag_window+1 #time dimension
# nl<- lag_window #lag dimenstion 
# ns<- ncell(SST1)
# 
# 
# dim(SST_array[,,1]) # slice of 64800 by 300 or ns*nt for lag 1
# layerNames(SST_rast)<-lf
# SST_sgdf<-as(SST_rast,"SpatialGridDataFrame")
# SST_df<-as.data.frame(SST_rast) #drop the coordinates x and y
# 
# #312*64800/(300*64800*13)=0.08 0r about 12.5 more pixels...
# #64800*13
# #need to add x,y,w column where w=weight which is a function of lat...
# 
# #SST_array[,,1]<- as.data.frame(stack(lag_list_obj$list_lf[[1]]))
# subset(SST_sgdf,SST_sgdf)
# 
# 
# t1<-stack(list_t[[1]])
# levelplot(t1)
# df1<-as.data.frame(t1)
# df1_v<-as.vector(t(df1))   #This transform the data frame into a one colum
# 
# t2<-stack(list_t[[2]])
# levelplot(t2)
# df2<-as.data.frame(t2)
# df2_v<-as.vector(t(df2))   #This transform the data frame into a one colum
# head(df2)
# test<-rbind(df1[,1],df1[,2])
# #MSSA1_l<-lag(d_z$MSSA1,k=lag_nb)
# 
# SST1
# ncol_concat<-ncol(SST1)*lag_window
# test<-raster(nrow=nrow(SST1),ncol=ncol_concat)
# #MSSA3_l<-lag(d_z$MSSA3,k=lag_nb)
# projection(test)<-NA
# extent(test)<-c(1,ncol_concat+1,1,181)
# 
# j=lag_number
# projection(t1_1)<-NA
# j=1
# 
# t1_1<-subset(t1,1)
# t1_2<-subset(t1,2)
# #Use the shitf function from the raster package???
# for (j in 1:lag_window){
#   extent(t1_1)<-c(1+(360*j),361+(360*j),1,181)
# }
# 
# #mosaic together
# #change extent for each of layer

#https://stat.ethz.ch/pipermail/r-help/2004-February/046067.html
#loadings(fa)
#new <- varimax(fa$loadings, normalize = FALSE)
#class(new$loadings) <- "loadings"
#loadings(new)

# ML with varimax rotation:
#fact(x=life,method="norm",rotation="varimax",maxfactors=3)
