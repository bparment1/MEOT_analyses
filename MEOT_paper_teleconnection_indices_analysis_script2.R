####################################    MEOT-MSSA PAPER   #######################################
############################      SCRIPT 2- MEOT          #######################################
#This script carries out a MSSA with varimax rotation for a timeseries of 312 SST images.       #
#Note that spatial patterns from MEOT and MSSA components are not analyzed in this script       #                 
#AUTHOR: Benoit Parmentier                                                                      #
#DATE: 03/25/2013            
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

load_obj <- function(f) 
{
  env <- new.env()
  nm <- load(f, env)[1]  
  env[[nm]]
}

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
####
pca_to_raster_fun<-function(pc_spdf,ref_raster=SSTm1,lag_window,out_prefix){
  #Input arguments:
  #pc_spdf: must include x,y and lag in the last 3 columns!!!
  npc<-ncol(pc_spdf)-3
  pc_scores_lf<-vector("list",npc)
  list_lag_pc<-vector("list",lag_window)
  for (k in 1:npc){
    tmp<-pc_spdf[,k]
    tmp$lag<-pc_spdf$lag
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

### Parameters and argument

infile1<-"SAODI-01-1854_06-2011_test.asc"             #GHCN shapefile containing variables for modeling 2010                 
infile2<-"SAODI-01-1854_06-2011.csv"                     #List of 10 dates for the regression
infile3<-"MEOT_MSSA_Telcon_indices_12112012.xlsx"                        #LST dates name
infile4<-"mask_rgf_1_1.rst"

npc<-20 # number of pca to produce...put this at the beginning of the script...
lag_window<-13
telind<-c("PNA","NAO","TNA","TSA","SAOD","MEI","PDO","AO","AAO","AMM","AMOsm","QBO")
mode_list_MEOT<-c("MEOT1","MEOT3", "MEOT4","MEOT7","MEOT10","MEOT15","MEOT16")
mode_list_PCA<-c("MSSA1","MSSA2","MSSA3","MSSA4","MSSA5","MSSA6")    
#mode_list_PCA<-paste("MSSA",1:15,sep="")
out_prefix<-"MEOT_paper_03222013"

# on Benoit Mac
in_path<-"/Users/benoitparmentier/Dropbox/Data/MEOT_paper/MEOT12272012/MEOT_working_dir_03102013"
# on Atlas:
#in_path<-"/home/parmentier/Data/MEOT12272012/MEOT_working_dir_03102013/"

setwd(in_path)

########################################

### STEP 0: READ IN DATASETS RELATED TO TELECONNECTION AND PREVIOUS LOADINGS

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

#Creating time series objects
d_ts<-ts(data=dat,start=c(1982,1), end=c(2006,12), frequency=12, names=names(dat))
colnames(d_ts)
d_z<-as.zoo(d_ts)  #### THIS IS THE TIME SERIES OBJECT USED LATER ON

##############################################################
#### STEP 1: load image time series data, apply latitude weight for area

lf<-mixedsort(list.files(pattern="ANOM_SST.*.rst"))
SST_s<-stack(lf)
SST1<-raster(lf[1])

mask_land<-raster(infile4)
mask_land_NA<-mask_land
mask_land_NA[mask_land_NA==0]<-NA
SST_rast<-mask(SST_s,mask_land_NA,filename="ANOM_SST_1982_2007.tif",overwrite=TRUE)
class(SST_rast)

#Multiply layer by weight

lat_coord<-coordinates(mask_land)[,2]
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
#### STEP 4: get lag images for MSSA/PCA from predited object

dim(Xpc) #547274*300??
SST_xy_lag<-SST_lag_data_df[,(nt+1):ncol(SST_lag_data_df)]
X_pc_data<-as.data.frame(Xpc) #Add coordinates and lag to pc scores
X_pc_data<-Xpc[,1:npc] #Add coordinates and lag to pc scores
X_pc_data<-cbind(X_pc_data,SST_xy_lag) #
coordinates(X_pc_data)<-X_pc_data[,c("s1","s2")] #promote to spdf
save(X_pc_data,file= paste("SST_lag_data_components_",out_prefix,".RData",sep=""))

tmp_names<-c(paste("PC",1:npc,sep=""),"s1","s2","lag")
names(X_pc_data)<-tmp_names

dim(Xpcv) #547274*300??
SST_xy_lag<-SST_lag_data_df[,(nt+1):ncol(SST_lag_data_df)]
X_pcv_data<-as.data.frame(Xpcv) #Add coordinates and lag to pc scores
X_pcv_data<-Xpcv[,1:npc] #Add coordinates and lag to pc scores
X_pcv_data<-cbind(X_pcv_data,SST_xy_lag) #
coordinates(X_pcv_data)<-X_pcv_data[,c("s1","s2")] #promote to spdf
out_prefix_pcv<-paste("varimax_",out_prefix,sep="")

save(X_pcv_data,file= paste("SST_lag_data_components_",out_prefix_pcv,".RData",sep=""))

tmp_names<-c(paste("PC",1:npc,sep=""),"s1","s2","lag")
names(X_pcv_data)<-tmp_names
out_prefix_pcv<-paste("varimax_",out_prefix,sep="")

X_test<-X_pcv_data[,c(1,2,21,22,23)]
#X_test<-X_pc_data[,c(1,2,21,22,23)]
#pc_scores_lf<-pca_to_raster_fun(X_test,ref_raster=SST1_m,lag_window,out_prefix)
pc_scores_lf<-pca_to_raster_fun(X_test,ref_raster=SST1_m,lag_window,out_prefix_pcv)
mssa1<-stack(pc_scores_lf[[1]])
plot(mssa1)

pc_scores_lf<-pca_to_raster_fun(X_pc_data,ref_raster=SST1_m,lag_window,out_prefix)
pcv_scores_lf<-pca_to_raster_fun(X_pcv_data,ref_raster=SST1_m,lag_window,out_prefix_pcv)

##############################################################
#### STEP 5: quick analysis of results

pca<-load_obj(paste("pca_SST_lag_data_df_",out_prefix,".RData",sep=""))
pca_varimax<-load_obj(paste("pca_varimax_SST_lag_data_df_",out_prefix,".RData",sep=""))

pca_loadings<-pca$loadings[,1:npc] #extract first component...
pca_v_loadings<-pca_varimax$loadings[,1:npc] #extract first component...
tmp_names<-c(paste("PC",1:npc,sep=""))
colnames(pca_loadings)<-tmp_names
colnames(pca_v_loadings)<-tmp_names

dat<-read.xls(infile3, sheet=1)
tail(dat$Date_label)
dat$SAOD<-in_SAODI  #Adding the SAOD index to all the MEOT/MSSA results in the data frame

cor(dat$MSSA1,pca$loadings[,1])  #Correlation between indices...
mssa_names<-paste("MSSA",1:15,sep="")
mssa_mat<-as.matrix(dat[,mssa_names])
telcon<-c("PNA","NAO","TNA","TSA","SAOD","MEI","PDO","AO","AAO","AMM","QBO")

tel_mat<-as.matrix(dat[,telcon])
#mssa_loadings_mat<-as.matrix(dat[,macth("MSSA*",names(dat))])
cor_mat<-cor(mssa_mat,pca_loadings)
diag(cor_mat)
image(cor_mat)
cor_mat<-cor(tel_mat,pca_loadings)
#cor(mssa_mat$MSSA1,pca_loadings$PC1)
j<-3
cor(mssa_mat[,j],pca_loadings[,j])
diag(cor_mat)

##With varimax components...
cor_mat<-cor(tel_mat,pca_v_loadings)
diag(cor_mat)

j<-1

plot(mssa_mat[,j],type="l")
par(new=TRUE)
#lines(pca_loadings[,j],type="l",col="red")
#plot(pca_loadings[,j],type="l",col="red")
plot(pca_v_loadings[,j],type="l",col="red")
cor(pca_v)
#Check spatial patterns
k=3
j=".?." #only one letter
raster_name<-paste("pc_component_",k,"_",j,"_",out_prefix,".rst$",sep="")
pc_scores_lf_tmp<-mixedsort(list.files(pattern=raster_name))
mssa<-stack(pc_scores_lf_tmp)
#quartz()
plot(mssa)
meot_rast<-mssa
temp.colors <- colorRampPalette(c('blue', 'lightgoldenrodyellow', 'red'))
layerNames(meot_rast)<-paste("Lag", 0:12,sep=" ")
title_plot<-paste("MSSA", "spatial sequence",sep=" ")
meot_rast_m<-mask(meot_rast,mask_land)
#levelplot(meot_rast_m,col.regions=temp.colors,par.settings = list(axis.text = list(font = 2, cex = 1)))
levelplot(meot_rast_m,main=title_plot, ylab=NULL,xlab=NULL,,par.settings = list(axis.text = list(font = 2, cex = 1.5),
          par.main.text=list(font=2,cex=2.2),strip.background=list(col="white")),par.strip.text=list(font=2,cex=1.5),
          col.regions=temp.colors,at=seq(-6,6,by=0.02))

s_range<-c(minValue(meot_rast),maxValue(meot_rast)) #stack min and max

####Check spatial patterns for MSSA varimax
k=3
j=".?." #only one letter
raster_name<-paste("pc_component_",k,"_",j,"_",out_prefix_pcv,".rst$",sep="")
pc_scores_lf_tmp<-mixedsort(list.files(pattern=raster_name))
mssa<-stack(pc_scores_lf_tmp)
#quartz()
plot(mssa)
meot_rast<-mssa
temp.colors <- colorRampPalette(c('blue', 'lightgoldenrodyellow', 'red'))
layerNames(meot_rast)<-paste("Lag", 0:12,sep=" ")
title_plot<-paste("MSSA", "spatial sequence",sep=" ")
meot_rast_m<-mask(meot_rast,mask_land)
#levelplot(meot_rast_m,col.regions=temp.colors,par.settings = list(axis.text = list(font = 2, cex = 1)))
levelplot(meot_rast_m,main=title_plot, ylab=NULL,xlab=NULL,,par.settings = list(axis.text = list(font = 2, cex = 1.5),
                                                                                par.main.text=list(font=2,cex=2.2),strip.background=list(col="white")),par.strip.text=list(font=2,cex=1.5),
          col.regions=temp.colors,at=seq(-6,6,by=0.02))

s_range<-c(minValue(meot_rast),maxValue(meot_rast)) #stack min and max

##Checking eigenvalues:
sum(pca$values)/(300-1)
sum(values(w_rast_m),na.rm=TRUE)*13
sum(values(mask_land),na.rm=TRUE)*13

#####FUNCTION
# NOW RUN ANALYSIS FOR PCA AND CREATE FIGURES..
#PCA tables and figs: now create a table with maximum correlation and corresponding lag.
#Table3. Maximum absolute value for lag cross correlations between climate indices and MSSA modes (lag in parenthesis).

#Creating time series objects

pca_loadings<-pca$loadings[,1:npc] #extract first component...
pca_v_loadings<-pca_varimax$loadings[,1:npc] #extract first component...
tmp_names<-c(paste("PCn",1:npc,sep="")) #PCA without rotation
colnames(pca_loadings)<-tmp_names
tmp_names<-c(paste("PCv",1:npc,sep="")) #PCA with varimax rotation
colnames(pca_v_loadings)<-tmp_names
dat2<-cbind(dat,as.data.frame(pca_loadings))
dat2<-cbind(dat2,as.data.frame(pca_v_loadings))
d_ts<-ts(data=dat2,start=c(1982,1), end=c(2006,12), frequency=12, names=names(dat2))
colnames(d_ts)
d_z2<-as.zoo(d_ts)  #### THIS IS THE TIME SERIES OBJECT USED LATER ON

telind
mode_list<-tmp_names<-c(paste("PCn",1:npc,sep="")) #PCA with varimax rotation
out_prefix_n<-paste("PCn_",out_prefix,sep="")
#telind<-mode_list_PCA : if cross cor desired
pcn_obj<-crosscor_lag_analysis_fun(telind,mode_list,d_z2,lag_window,fig=TRUE,out_prefix_n)
mssa_obj<-crosscor_lag_analysis_fun(telind,mode_list_PCA,d_z2,lag_window,fig=FALSE,out_prefix_n)

#debug(crosscor_lag_analysis_fun)

crosscor_lag_analysis_fun<-function(telind,mode_list,d_z,lag_window,fig,out_prefix){
  #This function crosss correlates between two sets of time series given some lag window.
  #Arguments:
  #1)telind: time series 1 as character vector
  #2)modelist: time series 2 as character vector
  #3)d_z: zoo object 
  #4)lag_window:
  #5)fig:
  
  lag_table_ext<-matrix(data=NA,nrow=length(telind),ncol=length(mode_list))
  lag_table_lag<-matrix(data=NA,nrow=length(telind),ncol=length(mode_list))
  lag_table_text<-matrix(data=NA,nrow=length(telind),ncol=length(mode_list)) #Formatted table used in the paper
  #lag_cross_cor_PCA<-vector("list",length(mode_list))
  lag_m<-seq(-1*lag_window,lag_window,1)
  #retain ccf!!!
  #list_ccf_lag_table
    
  #lag_cross_cor_PCA_m<-array(data=NA,nrow=length(lag_m),ncol=length(mode_list))
  for (i in 1:length(telind)){
    telindex<-telind[i]
    pos1<-match(telindex,names(d_z))
    #retain ccf!!!
    for (j in 1:length(mode_list)){
      mode_n<-mode_list[j]
      pos2<-match(mode_n,names(d_z))
      ccf_obj<-ccf(d_z[,pos1],d_z[,pos2], lag=lag_window)  #Note that ccf does not take
      
      lag_m<-seq(-1*lag_window,lag_window,1)
      ccf_obj$lag[,1,1]<-lag_m  #replacing lag values because continuous
      
      if (fig=="TRUE"){
        plot_name<-paste(telindex, "and", mode_n,"lag analysis",sep="_")
        png(paste(plot_name,"_",out_prefix,".png", sep=""))
        plot(ccf_obj, main= paste(telindex, "and", mode_n,"lag analysis",sep=" "), ylab="Cross-correlation",
             xlab="Lag (month)", ylim=c(-1,1))
        dev.off()
      }
      
      ######### NOW FIND THE m
      absext <-max(abs(ccf_obj$acf)) # maximum of the extremum
      pos<-match(absext,ccf_obj$acf) #find the position and lag, if NA it means it was negative
      if (is.na(pos)) {
        pos<-match(absext*-1,ccf_obj$acf)
        absext<-absext*-1   #recover the sign
      } 
      absext_lag<-ccf_obj$lag[pos,1,1] #This is the lag corresponding to the maximum absolute value
      
      lag_table_ext[i,j]<-absext
      lag_table_lag[i,j]<-absext_lag
      #number<-format(absext,digits=3)
      ext<-round(absext,digits=3)
      element<-paste(ext," (",absext_lag,")",sep="")
      lag_table_text[i,j]<-element
      
      ##Keep ccf lag somewhere
    }
  }
  
  lag_table_ext<-as.data.frame(lag_table_ext)
  names(lag_table_ext)<-mode_list
  rownames(lag_table_ext)<-telind
  file_name<-paste("lag_table_extremum_window_", lag_window,"_",out_prefix,".txt",sep="")
  write.table(lag_table_ext,file=file_name,sep=",")
  
  lag_table_lag<-as.data.frame(lag_table_lag)
  names(lag_table_lag)<-mode_list
  rownames(lag_table_lag)<-telind
  file_name<-paste("lag_table_lag_extremum_window_", lag_window,"_",out_prefix,".txt",sep="")
  write.table(lag_table_lag,file=file_name,sep=",")
  
  lag_table_text<-as.data.frame(lag_table_text)
  names(lag_table_text)<-mode_list
  rownames(lag_table_text)<-telind
  file_name<-paste("lag_table_lag_ext_text", lag_window,"_",out_prefix,".txt",sep="")
  write.table(lag_table_text,file=file_name,sep=",")
  
  #create return object
  
  crosscor_obj<-list(lag_table_ext,lag_table_lag,lag_table_text)
  names(crosscor_obj)<-c("extremum","lag_ext","text")
  file_name<-paste("crosscor_obj_lag_analysis_", lag_window,"_",out_prefix,".RData",sep="")
  save(crosscor_obj,file=file_name)
  
  return(crosscor_obj)
}

##############
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
