####################################    MEOT-MSSA PAPER   #######################################
############################      SCRIPT 2- MEOT          #######################################
#This script carries out a MSSA with varimax rotation for a timeseries of 312 SST images.       #
#Note that spatial patterns from MEOT and MSSA components are not analyzed in this script       #                 
#AUTHOR: Benoit Parmentier                                                                      #
#DATE: 03/22/2013            
#Version: 3
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

## Functions used in the script

### TO DO ... format the data for MSSA create a function for the window of any size ...use lag???
#unction ot lage files

lag_grouping<-function(lag_window,list_lf){
  #This function creates groups of variables corresponding to lag based on a list of files.
  #It requires 2 inputs:                                           
  # 1) lag_window: sed
  # 2) list_lf: list of names of variables from the a time series
  #The output is a list of four shapefile names produced by the function:
  # 1) list_obj with two lists:
  #  -list_lf contains gorupings of files for every lag (13 groups if lag_window is 13)
  #  -list_t contains gorupings of files for every lag (13 groups if lag_window is 13)
  
  
  #AUTHOR: Benoit Parmentier                                                                       
  #DATE: 03/22/2013                                                                                 
  #PROJECT: Clark Labs Climate predction- MEOT/MSSA paper     
  #Comments and TODO
  #
  ##################################################################################################
  
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

### Parameters and argument

infile1<-"SAODI-01-1854_06-2011_test.asc"             #GHCN shapefile containing variables for modeling 2010                 
infile2<-"SAODI-01-1854_06-2011.csv"                     #List of 10 dates for the regression
#infile2<-"list_365_dates_04212012.txt"
infile3<-"MEOT_MSSA_Telcon_indices_08062012.xlsx"                        #LST dates name
infile3<-"MEOT_MSSA_Telcon_indices_12112012.xlsx"                        #LST dates name
infile4<-"mask_rgf_1_1.rst"

# on Benoit Mac
#in_path<-"/Users/benoitparmentier/Dropbox/Data/MEOT_paper/MEOT12272012/MEOT_working_dir_03102013"

# on Atlas:
in_path<-"/home/parmentier/Data/MEOT12272012/MEOT_working_dir_03102013/"
npc<-20 # number of pca to produce...put this at the beginning of the script...

setwd(in_path)

#on MAC:
#path<-"/Users/benoitparmentier/Documents/DATA/Benoit/Clark_University/Paper_writings/MSSA_BNP/"
# on Atlas:
#path<-"/home/parmentier/Data/MEOT12272012/MEOT_working_dir_03102013/"
#path_raster<-"/Users/benoitparmentier/Dropbox/Data/MEOT_paper/MEOT12272012/MEOT_working_dir_10232012/MSSA_EEOT_04_29_09"
telind<-c("PNA","NAO","TNA","TSA","SAOD","MEI","PDO","AO","AAO","AMM","AMOsm","QBO")
mode_list_MEOT<-c("MEOT1","MEOT3", "MEOT4","MEOT7","MEOT10","MEOT15","MEOT16")
#mode_list_MEOT<-paste("MEOT",1:25,sep="")
#c("MEOT1","MEOT3", "MEOT4","MEOT7","MEOT10","MEOT15","MEOT16")

out_prefix<-"MEOT_paper_03222013"

########################################

#dates <-readLines(paste(path,"/",infile2, sep=""))
SAODI<-read.table(infile2,sep=",", header=TRUE)

#Prepare data to write out in a textfile
s_SAODI2 <- subset(SAODI, year>1981 & year<2008) #Subset rows that correspond to the conditions
s_SAODI2<- subset(s_SAODI2, select=-year) #Remove the column year
in_SAODI2<-as.vector(t(s_SAODI2))   #This transform the data frame into a one colum
write.table(in_SAODI2,file=paste("SAOD_index_1981_2007",out_prefix,".txt",sep=""),sep=",")

#Prepare data for cross lag correlation analysis and write out results
s_SAODI <- subset(SAODI, year>1981 & year<2007) #Subset rows that correspond to the conditions
s_SAODI<- subset(s_SAODI, select=-year) #Remove the column year
in_SAODI<-as.vector(t(s_SAODI))   #This transform the data frame into a one colum
write.table(in_SAODI,file=paste("SAOD_index_1981_2006",out_prefix,".txt",sep=""),sep=",")

#Import results from MEOT and MSSA analyses with teleconneciton indices
dat<-read.xls(infile3, sheet=1)
tail(dat$Date_label)
dat$SAOD<-in_SAODI  #Adding the SAOD index to all the MEOT/MSSA results in the data frame

cor(dat$SAOD,dat$TSA)  #Correlation between indices...
cor(dat$SAOD,dat$TNA)
cor(dat$SAOD,dat$AMM)

#Creating time series objects
d_ts<-ts(data=dat,start=c(1982,1), end=c(2006,12), frequency=12, names=names(dat))
colnames(d_ts)
d_z<-as.zoo(d_ts)  #### THIS IS THE TIME SERIES OBJECT USED LATER ON
d_xts<-as.xts(d_ts)

## find all 7th of the month between two dates, the last being a 7th.
st <- as.Date("1982-1-15")
en <- as.Date("2006-12-15")
dseq <- seq(st, en, by="month") #Creating monthly date sequence to create a time series from a data frame
d_z2<-zoo(dat,dseq)
time(d_z)  #no time stamp??
time(d_z2)

##############################################################
#### STEP 1: load data, apply latitude weight

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
values(w_rast)<-cos(lat_coord*pi/180) #use area in raster package for weight of a cell??
w_rast_m<-mask(w_rast,mask_land_NA)

SST_rast<-SST_rast*w_rast #Maybe use overlay to avoid putting in memory?
writeRaster(SST_rast,filename="ANOM_SST_1982_2007_weighted.tif",overwrite=TRUE)
rm(SST_rast)
SST_rast<-brick("ANOM_SST_1982_2007_weighted.tif")
#Do all in one step using overlay
#SST_rast<-overlay(x=SST_rast,y=w_rast,filename="ANOM_SST_1982_2007_weighted.tif",overwrite=TRUE,
#                   fun=function(x,y) {(return(x*y)})

SST1_m<-subset(SST_rast,1)
plot(stack(SST1,SST1_m))
#rm(SST1,SST1_m)
#Need to add weighting scheme base on lat long...

##############################################################
#### STEP 2: get data ready for PCA by lagging data frame...

SST_sgdf<-as(SST_rast,"SpatialGridDataFrame")
#SST_df<-as.data.frame(SST_rast) #drop the coordinates x and y
SST_df<-as.data.frame(SST_sgdf)
rm(SST_sgdf)

nt<- length(lf)-lag_window+1 #time dimension
nl<- lag_window #lag dimenstion 
ns<- ncell(SST1)
dimension<-list(nt,nl,ns)
ns_NA<-nrow(SST_df) #42098 if removing land areas from image
lag_list_obj<-lag_grouping(lag_window,lf) 
#list_lf contains gorupings of files for every lag (13 groups if lag_window is 13)
#list_t contains gorupings of files for every lag (13 groups if lag_window is 13)
list_lf<-lag_list_obj$list_lf
list_lag_df<-vector("list",length(lag_window))
names(SST_df)<-c(lf,"s1","s2")
new_names<-paste("t",1:nt,sep="_")

#Assign new names before doing rbind to avoid problems
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

##############################################################
#### STEP 3: prepare covariance matrix for PCA 

#load("~/Dropbox/Data/MEOT_paper/MEOT12272012/MEOT_working_dir_10232012/SST_1982_2007/SST_lag_data_df_MEOT_paper_03122013.RData")

X<-as.matrix(SST_lag_data_df[,1:nt]) #drop x, y and lag column
Xt<-t(X)
Xt<-scale(Xt)
X<-t(Xt)
cX<-Xt%*%X
#cX<-cX/ncol(X)

pca<-principal(r=cX, nfactors = npc, residuals = FALSE, covar=TRUE,rotate = "none")
#pca_test<-principal(r=cX, nfactors = npc, residuals = FALSE, covar=TRUE,scores=TURE,rotate = "none")

pca_varimax<-principal(r=cX, nfactors = npc, residuals = FALSE, covar=TRUE,rotate = "varimax")

Xpcv <-predict(pca_varimax,X) #with varimax
Xpc <-predict(pca,X) #without varimax

#Note that covar=TRUE ensures that the square matrix r is not standardized!!

save(pca,file= paste("pca_SST_lag_data_df_",out_prefix,".RData",sep=""))
save(pca_varimax,file= paste("pca_varimax_SST_lag_data_df_",out_prefix,".RData",sep=""))
save(X,file= paste("SST_lag_data_matrix_nt_",out_prefix,".RData",sep=""))
save(Xpc,file= paste("SST_lag_data_matrix_components_",out_prefix,".RData",sep=""))
save(Xpcv,file= paste("SST_lag_data_matrix_components_varimax_",out_prefix,".RData",sep=""))

plot(pca_varimax)

##############################################################
#### STEP 4: get lag images for MSSA/PCA

dim(Xpc) #547274*300??
SST_xy_lag<-SST_lag_data_df[,(nt+1):ncol(SST_lag_data_df)]
X_pc_data<-as.data.frame(Xpc) #Add coordinates and lag to pc scores
X_pc_data<-Xpc[,1:npc] #Add coordinates and lag to pc scores
X_pc_data<-cbind(X_pc_data,SST_xy_lag) #
coordinates(X_pc_data)<-X_pc_data[,c("s1","s2")] #promote to spdf
save(X_pc_data,file= paste("SST_lag_data_components_",out_prefix,".RData",sep=""))

tmp_names<-paste("PC",1:npc,sep="")
names(X_pc_data)<-tmp_names
#coordinates(X_pc_data)<-X_pc_data[,c("s1","s2")] #promote to spdf
#X_pc_data<-X_pc_data[,c("s1","s2")] #promote to spdf

#pc1_l1<-rasterize(tmp,SST1_m)??
#loop through k and j with j being lag index and k pc index

#for (j in 1:lag_window){
#  tmp<-subset(X_pc_data,X_pc_data$lag==j)
#  for (k in 1:npc){
#    pc_name<-names(tmp[k])
#    pc_lag<-rasterize(tmp,SST1_m,pc_name,fun=min,overwrite=TRUE,
#                      filename=paste("pc_component_",k,"_",j,"_",out_prefix,".rst",sep=""))
#  }
#}

####
pca_to_raster_fun<-function(pc_spdf,ref_raster=SSTm1,lag_window,out_prefix){
  #Input arguments:
  #pc_spdf: must include x,y and lag in the last 3 columns!!!
  npc<-ncol(pc_spdf)-3
  pc_scores_lf<-vector("list",npc)
  list_lag_pc<-vector("list",lag_window)
  for (k in 1:npc){
    tmp<-pc_spdf[,k]
    pc_name<-names(pc_spdf)[k]
    for (j in 1:lag_window){
      pc_scores<-subset(tmp,tmp$lag==j)
      raster_name<-paste("pc_component_",k,"_",j,"_",out_prefix,".rst",sep="")
      pc_lag<-rasterize(pc_scores,ref_raster,pc_name,fun=min,overwrite=TRUE,
                        filename=raster_name)
      list_lag_pc[[j]]<-raster_name
    }
    pc_scores_lf[[k]]<-list_lag_pc
  } 
  return(pc_scores_lf)
}
X_test<-X_pc_data[,c(1,2,21,22,23)]
list_pc<-pca_to_raster_fun(X_test,ref_raster=SST1_m,lag_window,out_prefix)

pc_scores_lf<-pca_to_raster_fun(X_pc_data,ref_raster=SST1_m,lag_window,out_prefix)
  
#k=1
#pc1_lf<-mixedsort(list.files(pattern=paste("pc_component_",k,"_",".*.","_",out_prefix,".rst$",sep="")))
#pc1_l1<-rasterize(tmp,SST1_m,"PC1",fun=min)
#mssa1<-stack(pc_scores_lf)
mssa1<-stack(pc_scores_lf[[1]])
plot(mssa1)

##############################################################
#### STEP 5: quick analysis of results

pca$loadings[,1] #extract first component...

dat<-read.xls(infile3, sheet=1)
tail(dat$Date_label)
dat$SAOD<-in_SAODI  #Adding the SAOD index to all the MEOT/MSSA results in the data frame

cor(dat$MSSA1,pca$loadings[,1])  #Correlation between indices...
pca_loadings_mat<-pca$loadings[,]
#mssa_loadings_mat<-as.matrix(dat[,macth("MSSA*",names(dat))])

cor(dat$SAOD,dat$TNA)
cor(dat$SAOD,dat$AMM)

###############   END OF SCRIPT   ###################
#####################################################

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
