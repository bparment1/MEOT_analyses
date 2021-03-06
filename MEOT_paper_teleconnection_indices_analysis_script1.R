####################################    MEOT-MSSA PAPER   #######################################
############################      SCRIPT 1- MEOT          #######################################
#This script reads SAODI index and other telectonnection indices...                             #
#Cross-corrrelation and autocorrelation functions are generated to examine results of MEOT      #
#and MSSA analyses generated at Clark Labs.                                                     #
#Note that spatial patterns from MEOT and MSSA components are not analyzed in this script       #                 #
#AUTHOR: Benoit Parmentier                                                                      #
#DATE: 12/21/2012            
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
### Parameters and argument

infile1<-"SAODI-01-1854_06-2011_test.asc"             #GHCN shapefile containing variables for modeling 2010                 
infile2<-"SAODI-01-1854_06-2011.csv"                     #List of 10 dates for the regression
#infile2<-"list_365_dates_04212012.txt"
infile3<-"MEOT_MSSA_Telcon_indices_08062012.xlsx"                        #LST dates name
infile3<-"MEOT_MSSA_Telcon_indices_12112012.xlsx"                        #LST dates name

#path<-"/Users/benoitparmentier/Documents/DATA/Benoit/Clark_University/Paper_writings/MSSA_BNP/"
#on MAC:
#path<-"/Users/benoitparmentier/Dropbox/Data/MEOT_paper/MEOT12272012/MEOT_working_dir_10232012/MEOT_analysis_R_10102012"
#path<-"/Users/benoitparmentier/Documents/DATA/Benoit/Clark_University/Paper_writings/MSSA_BNP/MEOT_analysis_R_10102012"
# on Atlas:
path<-"/home/parmentier/Data/MEOT12272012/MEOT_working_dir_10232012/MEOT_analysis_R_10102012"
#path_raster<-"/Users/benoitparmentier/Dropbox/Data/MEOT_paper/MEOT12272012/MEOT_working_dir_10232012/MSSA_EEOT_04_29_09"


setwd(path)

out_prefix<-"MEOT_paper_01152013b_"
telind<-c("PNA","NAO","TNA","TSA","SAOD","MEI","PDO","AO","AAO","AMM","AMOsm","QBO")
mode_list_MEOT<-c("MEOT1","MEOT3", "MEOT4","MEOT7","MEOT10","MEOT15","MEOT16")
#mode_list_MEOT<-paste("MEOT",1:25,sep="")
#c("MEOT1","MEOT3", "MEOT4","MEOT7","MEOT10","MEOT15","MEOT16")

mode_list_PCA<-c("MSSA1","MSSA2","MSSA3","MSSA4","MSSA5","MSSA6")    
#mode_list_PCA<-paste("MSSA",1:15,sep="")

lag_window<-13

############### START OF SCRIPT  #####################

#Importing SAOD index...Note that it was reformatted from the first file in infile1

#dates <-readLines(paste(path,"/",infile2, sep=""))
SAODI<-read.table(paste(path,"/",infile2,sep=""),sep=",", header=TRUE)

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

######### START THE ANALYSIS

## PART I: EXPLORATORY ANALYSES WITH COMPARISON; comparison to other indices...

# AMM AND SAOD INDEX

tmp<-ccf(d_z$SAOD,d_z$AMM, lag=13)  #Note that ccf does not take
lag_m<- -13:13
tmp$lag[,1,1]<-lag_m  #replacign lag values because continuous
plot(tmp, main="Lag cross-correlation between SAOD and AMM",
     ylab="Cross-correlation",
     xlab="Lag (month)", ylim=c(-1,1),xlim=c(-13,13))
#place the tick tat the right place...                  
absext <-max(abs(tmp$acf)) # maximum of the extremum
pos<-match(absext,tmp$acf) #find the position and lag
absext_lag<-tmp$lag[pos,1,1] #This is the lag corresponding to the maximum absolute value
         
# TSA AND SAOD INDICES
         
tmp<-ccf(d_z$SAOD,d_z$TSA, lag=13)  #Note that ccf does not take
lag_m<- -13:13
tmp$lag[,1,1]<-lag_m  #replacign lag values because continuous
plot(tmp, main="Lag cross-correlation between SAOD and TSA",
     ylab="Cross-correlation",
     xlab="Lag (month)", ylim=c(-1,1))
         
absext <-max(abs(tmp$acf)) # maximum of the extremum
pos<-match(absext,tmp$acf) #find the position and lag,... NOT WORKING BECAUSE IT IS NEGATIVE
absext_lag<-tmp$lag[pos,1,1] #This is the lag corresponding to the maximum absolute value

# TNA AND SAOD INDICES
         
tmp<-ccf(d_z$SAOD,d_z$TNA, lag=13)  #Note that ccf does not take
lag_m<- -13:13
tmp$lag[,1,1]<-lag_m  #replacign lag values because continuous
plot(tmp, main="Lag cross-correlation between SAOD and TNA",
     ylab="Cross-correlation",
     xlab="Lag (month)", ylim=c(-1,1))
         
absext <-max(abs(tmp$acf)) # maximum of the extremum: you could use range!!! and match the position afterwords..
pos<-match(absext,tmp$acf) #find the position and lag, if NA it means it was negative
if (is.na(pos)) {
  pos<-match(absext*-1,tmp$acf)
} 
absext_lag<-tmp$lag[pos,1,1] #This is the lag corresponding to the maximum absolute value

# MEOT10/MEOT15 QUADRATURE...Figure??
         
tmp<-ccf(d_z$MEOT10,d_z$MEOT15, lag=13)  #Note that ccf does not take
lag_m<- -13:13
tmp$lag[,1,1]<-lag_m  #replacign lag values because continuous
plot(tmp, main="SAOD and MEOT10 lag analysis", ylab="Cross-correlation",
     xlab="Lag (month)", ylim=c(-1,1))

absext <-max(abs(tmp$acf)) # maximum of the extremum
pos<-match(absext,tmp$acf) #find the position and lag, if NA it means it was negative
if (is.na(pos)) {
  pos<-match(absext*-1,tmp$acf)
} 
absext_lag<-tmp$lag[pos,1,1] #This is the lag corresponding to the maximum absolute value


tmp<-ccf(d_z$SAOD,d_z$MEOT15, lag=13)  #Note that ccf does not take
lag_m<- -13:13
tmp$lag[,1,1]<-lag_m  #replacign lag values because continuous
plot(tmp, main="SAOD and MEOT15 lag analysis",ylab="Cross-correlation",
     xlab="Lag (month)", ylim=c(-1,1))
         
absext <-max(abs(tmp$acf)) # maximum of the extremum
pos<-match(absext,tmp$acf) #find the position and lag, if NA it means it was negative
if (is.na(pos)) {
pos<-match(absext*-1,tmp$acf)
} 
absext_lag<-tmp$lag[pos,1,1] #This is the lag corresponding to the maximum absolute value

### MEOT7/MEOT16 AND SAOD INDICES
tmp<-ccf(d_z$MEOT7,d_z$MEOT21, lag=13)  #Note that ccf does not take
tmp<-ccf(d_z$AMM,d_z$MEOT7, lag=13)  #Note that ccf does not take

tmp<-ccf(d_z$MEOT7,d_z$MEOT16, lag=13)  #Note that ccf does not take
lag_m<- -13:13
lag_m<-seq(-1*lag_window,lag_window,1)

tmp$lag[,1,1]<-lag_m  #replacign lag values because continuous
plot(tmp, main="SAOD and MEOT7 lag analysis", ylab="Cross-correlation",
     xlab="Lag (month)", ylim=c(-1,1))
labels<-c(-13,-11,-9,-7,-5,-3,0,3,5,7,9,11,13)

#using option h crates a thin stick like plot...

absext <-max(abs(tmp$acf)) # maximum of the extremum
pos<-match(absext,tmp$acf) #find the position and lag, if NA it means it was negative
if (is.na(pos)) {
  pos<-match(absext*-1,tmp$acf)
} 
absext_lag<-tmp$lag[pos,1,1] #This is the lag corresponding to the maximum absolute value
X11(6,10)
plot(d_z$MEOT7,col="blue",ylim=c(-1,1))
par(new=T)
plot(d_z$MEOT16,col="green", ylim=c(-1,1),axes=F)

tmp<-ccf(d_z$SAOD,d_z$MEOT16, lag=13)  #Note that ccf does not take
lag_m<- -13:13
tmp$lag[,1,1]<-lag_m  #replacign lag values because continuous
plot(tmp, main="SAOD and MEOT16 lag analysis",ylab="Cross-correlation",
     xlab="Lag (month)", ylim=c(-1,1))

absext <-max(abs(tmp$acf)) # maximum of the extremum
pos<-match(absext,tmp$acf) #find the position and lag, if NA it means it was negative
if (is.na(pos)) {
  pos<-match(absext*-1,tmp$acf)
} 
absext_lag<-tmp$lag[pos,1,1] #This is the lag corresponding to the maximum absolute value

## MEOT1/MEOT4 Lag function

tmp<-ccf(d_z$MEOT1,d_z$MEOT4, lag=13)  #Note that ccf does not take
lag_m<- -13:13
tmp$lag[,1,1]<-lag_m  #replacign lag values because continuous
plot(tmp, main="SAOD and MEOT10 lag analysis", ylab="Cross-correlation",
     xlab="Lag (month)", ylim=c(-1,1))

absext <-max(abs(tmp$acf)) # maximum of the extremum
pos<-match(absext,tmp$acf) #find the position and lag, if NA it means it was negative
if (is.na(pos)) {
  pos<-match(absext*-1,tmp$acf)
} 
absext_lag<-tmp$lag[pos,1,1] #This is the lag corresponding to the maximum absolute value

plot(d_z$MEOT1,col="blue",ylim=c(-1,1))
lines(d_z$MEOT4,col="green",ylim=c(-1,1))

tmp<-ccf(d_z$MEOT10,d_z$MEOT11, lag=13)  #Note that ccf does not take
tmp<-ccf(d_z$SAOD,d_z$MEOT7, lag=13)  #Note that ccf does not take

###############MET10/MEOT15
tmp<-ccf(d_z$MEOT10,d_z$MEOT15, lag=13)  #Note that ccf does not take
lag_m<- -13:13
tmp$lag[,1,1]<-lag_m  #replacign lag values because continuous
plot(tmp, main="SAOD and MEOT10 lag analysis", ylab="Cross-correlation",
     xlab="Lag (month)", ylim=c(-1,1))

absext <-max(abs(tmp$acf)) # maximum of the extremum
pos<-match(absext,tmp$acf) #find the position and lag, if NA it means it was negative
if (is.na(pos)) {
  pos<-match(absext*-1,tmp$acf)
} 
absext_lag<-tmp$lag[pos,1,1] #This is the lag corresponding to the maximum absolute value

## MEOT1/MEOT3 Lag function

tmp<-ccf(d_z$MEOT1,(d_z$MEOT3), lag=13)  #Note that ccf does not take
lag_m<- -13:13
lag_m<-seq(-1*lag_window,lag_window,1)

tmp$lag[,1,1]<-lag_m  #replacign lag values because continuous
plot(tmp, main="MEOT1 and MEOT3 lag analysis", ylab="Cross-correlation",
     xlab="Lag (month)", ylim=c(-1,1))
labels<-c(-13,-11,-9,-7,-5,-3,0,3,5,7,9,11,13)

absext <-max(abs(tmp$acf)) # maximum of the extremum
pos<-match(absext,tmp$acf) #find the position and lag, if NA it means it was negative
if (is.na(pos)) {
  pos<-match(absext*-1,tmp$acf)
} 
absext_lag<-tmp$lag[pos,1,1] #This is the lag corresponding to the maximum absolute value

lag_nb<--11
MSSA1_l<-lag(d_z$MSSA1,k=lag_nb)
MSSA3_l<-lag(d_z$MSSA3,k=lag_nb)

MSSA3<-d_z$MSSA3
MSSA1<-d_z$MSSA1

cor(MEOT1_l6,MEOT3)
plot(MSSA3,col="green")
lines(MSSA1_l,col="blue")  #The figure shows the series in sync
lines(MSSA1, col="red")
### note that positive lag means 6 months ahead!!!

## verify results: lag MEOT1 by 6 i.e. start half a year earlier...
lag_nb<-6
MEOT1_l6<-lag(d_z$MEOT1,k=lag_nb)
MEOT3_l6<-lag(d_z$MEOT1,k=lag_nb)

MEOT3_t6<-d_z$MEOT3[1:(length(d_z$MEOT3)-lag_nb)]
MEOT3<-d_z$MEOT3
MEOT1<-d_z$MEOT1

cor(MEOT1_l6,MEOT3)
plot(MEOT3,col="green")
lines(MEOT1_l6,col="blue")  #The figure shows the series in sync
lines(MEOT1, col="red")
### note that positive lag means 6 months ahead!!!

##########################################################################
### PART II : ANALYSIS FOR PAPER -PREPARING TABLEs: table 2 and table 3
#MEOT tables and figs
#Table2. Maximum absolute value for lag cross correlations between climate indices and MEOTs (lag in parenthesis).

mode_list<-mode_list_MEOT    
#telind<-mode_list_MEOT
lag_table_ext<-matrix(data=NA,nrow=length(telind),ncol=length(mode_list))
lag_table_lag<-matrix(data=NA,nrow=length(telind),ncol=length(mode_list))
lag_table_text<-matrix(data=NA,nrow=length(telind),ncol=length(mode_list))

#Make this loop a function!!! and call it on a list with mode_list and teleconnection indices...
for (i in 1:length(telind)){
  telindex<-telind[i]
  pos1<-match(telindex,names(d_z))
  for (j in 1:length(mode_list)){
    mode_n<-mode_list[j]
    pos2<-match(mode_n,names(d_z))
    ccf_obj<-ccf(d_z[,pos1],d_z[,pos2], lag=13)  #Note that ccf does not take
    lag_m<-seq(-1*lag_window,lag_window,1)
    ccf_obj$lag[,1,1]<-lag_m  #replacign lag values because continuous
    #X11(type="cairo") #Cairo because it is macos...?
    plot_name<-paste(telindex, "and", mode_n,"lag analysis",sep="_")
    png(paste(plot_name,"_",out_prefix,".png", sep=""))
    #plot(ccf_obj, main= paste(telindex, "and", mode_n,"lag analysis",sep=" "), ylab="Cross-correlation",
    #     xlab="Lag (month)", ylim=c(-1,1))
    plot(ccf_obj, main= paste(telindex, "and", mode_n,"lag analysis",sep=" "), ylab="Cross-correlation",
         xlab="Lag (month)", ylim=c(-1,1))
    #plot_name<-paste(telindex, "and", mode_n,"lag analysis",sep="_")
    #savePlot(paste(plot_name,"_",out_prefix,".png", sep=""), type="png")
    dev.off()
    absext <-max(abs(ccf_obj$acf)) # maximum of the extremum
    pos<-match(absext,ccf_obj$acf) #find the position and lag, if NA it means it was negative
    if (is.na(pos)) {
      pos<-match(absext*-1,ccf_obj$acf)
      absext<-absext*-1
      
    } 
    absext_lag<-ccf_obj$lag[pos,1,1] #This is the lag corresponding to the maximum absolute value
    
    lag_table_ext[i,j]<-absext
    lag_table_lag[i,j]<-absext_lag
    #number<-format(absext,digits=3)
    ext<-round(absext,digits=3)
    ext<-format(ext,digits=3)
    element<-paste(ext," (",absext_lag,")",sep="")
    lag_table_text[i,j]<-element
  }
}

lag_table_ext<-as.data.frame(lag_table_ext)
names(lag_table_ext)<-mode_list
rownames(lag_table_ext)<-telind
file_name<-paste("MEOT_lag_table_extremum_window_", lag_window,"_",out_prefix,".txt",sep="")
write.table(lag_table_ext,file=file_name,sep=",")

lag_table_lag<-as.data.frame(lag_table_lag)
names(lag_table_lag)<-mode_list
rownames(lag_table_lag)<-telind
file_name<-paste("MEOT_lag_table_lag_extremum_window_", lag_window,"_",out_prefix,".txt",sep="")
write.table(lag_table_lag,file=file_name,sep=",")

lag_table_text<-as.data.frame(lag_table_text)
names(lag_table_text)<-mode_list
rownames(lag_table_text)<-telind
file_name<-paste("MEOT_lag_table_lag_ext_text", lag_window,"_",out_prefix,".txt",sep="")
write.table(lag_table_text,file=file_name,sep=",")

meot_obj<-list(lag_table_ext,lag_table_lag,lag_table_text)
names(meot_obj)<-c("extremum","lag_ext","text")
file_name<-paste("meot_obj_lag_analysis_", lag_window,"_",out_prefix,".RData",sep="")
#file_name<-paste("meot_obj_lag_cross_correlation_MEOT", lag_window,"_",out_prefix,".RData",sep="")
save(meot_obj,file=file_name)

# NOW RUN ANALYSIS FOR PCA AND CREATE FIGURES..
#PCA tables and figs: now create a table with maximum correlation and corresponding lag.
#Table3. Maximum absolute value for lag cross correlations between climate indices and MSSA modes (lag in parenthesis).


mode_list<-mode_list_PCA
#telind<-mode_list_PCA : if cross cor desired

lag_table_ext<-matrix(data=NA,nrow=length(telind),ncol=length(mode_list))
lag_table_lag<-matrix(data=NA,nrow=length(telind),ncol=length(mode_list))
lag_table_text<-matrix(data=NA,nrow=length(telind),ncol=length(mode_list)) #Formatted table used in the paper
lag_cross_cor_PCA<-vector("list",length(mode_list))
lag_m<-seq(-1*lag_window,lag_window,1)

#lag_cross_cor_PCA_m<-array(data=NA,nrow=length(lag_m),ncol=length(mode_list))
for (i in 1:length(telind)){
  telindex<-telind[i]
  pos1<-match(telindex,names(d_z))
  for (j in 1:length(mode_list)){
    mode_n<-mode_list[j]
    pos2<-match(mode_n,names(d_z))
    #if(mode_n=="MSSA3"){
    #  ccf_obj<-ccf(d_z[,pos1],(d_z[,pos2]*-1), lag=13)  #Addded on 11/06
    #} 
    #if(mode_n!="MSSA3"){
    #  ccf_obj<-ccf(d_z[,pos1],d_z[,pos2], lag=13)  #Addded on 11/06
    #}
    ccf_obj<-ccf(d_z[,pos1],d_z[,pos2], lag=13)  #Note that ccf does not take

    lag_m<-seq(-1*lag_window,lag_window,1)
    ccf_obj$lag[,1,1]<-lag_m  #replacign lag values because continuous
    plot_name<-paste(telindex, "and", mode_n,"lag analysis",sep="_")
    png(paste(plot_name,"_",out_prefix,".png", sep=""))
    #plot(ccf_obj, main= paste(telindex, "and", mode_n,"lag analysis",sep=" "), ylab="Cross-correlation",
    #     xlab="Lag (month)", ylim=c(-1,1))
    plot(ccf_obj, main= paste(telindex, "and", mode_n,"lag analysis",sep=" "), ylab="Cross-correlation",
         xlab="Lag (month)", ylim=c(-1,1))
    #plot_name<-paste(telindex, "and", mode_n,"lag analysis",sep="_")
    #savePlot(paste(plot_name,"_",out_prefix,".png", sep=""), type="png")
    dev.off()
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
  }
}

lag_table_ext<-as.data.frame(lag_table_ext)
names(lag_table_ext)<-mode_list
rownames(lag_table_ext)<-telind
file_name<-paste("MSSA_lag_table_extremum_window_", lag_window,"_",out_prefix,".txt",sep="")
write.table(lag_table_ext,file=file_name,sep=",")

lag_table_lag<-as.data.frame(lag_table_lag)
names(lag_table_lag)<-mode_list
rownames(lag_table_lag)<-telind
file_name<-paste("MSSA_lag_table_lag_extremum_window_", lag_window,"_",out_prefix,".txt",sep="")
write.table(lag_table_lag,file=file_name,sep=",")

lag_table_text<-as.data.frame(lag_table_text)
names(lag_table_text)<-mode_list
rownames(lag_table_text)<-telind
file_name<-paste("MSSA_lag_table_lag_ext_text", lag_window,"_",out_prefix,".txt",sep="")
write.table(lag_table_text,file=file_name,sep=",")

#lag_cross_cor_PCA_m<-as.data.frame(lag_table_text)
#names(lag_cross_cor_PCA_m)<-mode_list
#rownames(lag_cross_cor_PCA_m)<-
#file_name<-paste("MSSA_lag_table_lag_ext_text", lag_window,"_",out_prefix,".txt",sep="")
#write.table(lag_table_text,file=file_name,sep=",")

mssa_obj<-list(lag_table_ext,lag_table_lag,lag_table_text)
names(mssa_obj)<-c("extremum","lag_ext","text")
file_name<-paste("mssa_obj_lag_analysis_", lag_window,"_",out_prefix,".RData",sep="")
#file_name<-paste("mssa_obj_lag_cross_cor_", lag_window,"_",out_prefix,".RData",sep="")
save(mssa_obj,file=file_name)

# formated figure for lag correlatioon analysis (bar plot)

### PART III: formated figure for lag correlation analysis (bar plot)

#MEOT 1 and MEOT3

list_fig_MEOT<-vector("list",3)
list_fig_MEOT[[1]]<- c("MEOT1","MEOT3","Figure_5_paper_MEOT1_MEOT3_barplot_lag_crosscorrelation")
list_fig_MEOT[[2]]<- c("MEOT7","MEOT16","Figure_8_paper_MEOT7_MEOT16_barplot_lag_crosscorrelation")
list_fig_MEOT[[3]]<- c("MEOT10","MEOT15","Figure_11_paper_MEOT10_MEOT15_barplot_lag_crosscorrelation")

#start the loop Make it a function...
for (i in 1:length(list_fig_MEOT)){
  pos1<-match(list_fig_MEOT[[i]][[1]],names(meot_obj$extremum))
  pos2<-match(list_fig_MEOT[[i]][[2]],names(meot_obj$extremum))
  names_ind<-c("PNA","NAO","TNA","TSA","SAOD","MEI","PDO","AO","AAO","AMM","AMO","QBO")
  data_plot<-cbind(as.vector(meot_obj$extremum[pos1]),as.vector(meot_obj$extremum[pos2]))
  
  #First plot without shading...
  plot_name<- paste(list_fig_MEOT[[i]][[3]],"_",out_prefix,".png",sep="")
  png(plot_name)
  heights<-as.matrix(t(data_plot))
  barplot(heights,     #data to plot
          #main=paste(names(data_plot)[pos1]," and ",names(data_plot)[pos2],sep=""),
          names.arg=names_ind,cex.names=0.8,   #names of the teleconnections indices and size of fonts of axis labes
          beside=TRUE,                         # see two barplots for comparisons...
          xlab="Teleconnection Indices",       # font.lab is 2 to make the font bold
          ylab="Lag crosscorrelation",font.lab=2,
          col=c("blue","red"),
          #density=c(0,10),
          #add=TRUE,
          ylim=c(-1,1))
  grid(nx=12,ny=10)
  legend("topright",legend=c(list_fig_MEOT[[i]][[1]],list_fig_MEOT[[i]][[2]]),pt.cex=0.9,fill=c("blue","red"),bty="n")
  #no box around legend
  box()
  dev.off()
  
  #Second plot with shading...
  #Now add shading to plot
  plot_name<- paste(list_fig_MEOT[[i]][[3]],"_",out_prefix,"shading.png",sep="")
  png(plot_name)
  barplot(heights,     #data to plot
          #main=paste(names(data_plot)[pos1]," and ",names(data_plot)[pos2],sep=""),
          names.arg=names_ind,cex.names=0.8,   #names of the teleconnections indices and size of fonts of axis labes
          beside=TRUE,                         # see two barplots for comparisons...
          xlab="Teleconnection Indices",       # font.lab is 2 to make the font bold
          ylab="Lag crosscorrelation",font.lab=2,
          col=c("blue","red"),
          #density=c(0,10),
          #add=TRUE,
          ylim=c(-1,1))
  grid(nx=12,ny=10)
  #legend("topright",legend=c(list_fig_MEOT[[i]][[1]],list_fig_MEOT[[i]][[2]]), cex=0.8, col=c("blue","red"),
   #      pch=15)
  barplot(heights,     #data to plot
          #main=paste(names(data_plot)[pos1]," and ",names(data_plot)[pos2],sep=""),
          names.arg=names_ind,cex.names=0.8,   #names of the teleconnections indices and size of fonts of axis labes
          beside=TRUE,                         # see two barplots for comparisons...
          xlab="Teleconnection Indices",       # font.lab is 2 to make the font bold
          ylab="Lag crosscorrelation",font.lab=2,
          col=c("blue","black"),
          density=c(0,10),
          add=TRUE,
          ylim=c(-1,1))
  #grid(nx=12,ny=10)
  legend("topright",legend=c(list_fig_MEOT[[i]][[1]],list_fig_MEOT[[i]][[2]]), pt.cex=1,fill=c("blue","red"),bty="n")
         #density=c(0,10))
  #legend("topright",legend=c(list_fig_MEOT[[i]][[1]],list_fig_MEOT[[i]][[2]]), pt.cex=1,fill=c("blue","red"),
  #density=c(0,100),add=TRUE)
  box()
  dev.off()
  
}

#######
######## BARPLOT FOR MSSA...#############

list_fig_MSSA<-vector("list",3)
list_fig_MSSA[[1]]<- c("MSSA1","MSSA3","Figure_18_paper_MSSA1_MSSA3_barplot_lag_crosscorrelation")
list_fig_MSSA[[2]]<- c("MSSA2","MSSA4","Figure_19_paper_MSSA2_MSSA4_barplot_lag_crosscorrelation")
list_fig_MSSA[[3]]<- c("MSSA5","MSSA6","Figure_19_paper_MSSA5_MSSA6_barplot_lag_crosscorrelation")

for (i in 1:length(list_fig_MSSA)){
  pos1<-match(list_fig_MSSA[[i]][[1]],names(mssa_obj$extremum))
  pos2<-match(list_fig_MSSA[[i]][[2]],names(mssa_obj$extremum))
  plot_name<- paste(list_fig_MSSA[[i]][[3]],"_",out_prefix,".png",sep="")
  names_ind<-c("PNA","NAO","TNA","TSA","SAOD","MEI","PDO","AO","AAO","AMM","AMO","QBO")
  data_plot<-cbind(as.vector(mssa_obj$extremum[pos1]),as.vector(mssa_obj$extremum[pos2]))
  
  png(plot_name)
  heights<-as.matrix(t(data_plot))
  barplot(heights,     #data to plot
          #main=paste(names(data_plot)[pos1]," and ",names(data_plot)[pos2],sep=""),
          names.arg=names_ind,cex.names=0.8,   #names of the teleconnections indices and size of fonts of axis labes
          beside=TRUE,                         # see two barplots for comparisons...
          xlab="Teleconnection Indices",       # font.lab is 2 to make the font bold
          ylab="Lag crosscorrelation",font.lab=2,
          col=c("blue","red"), ylim=c(-1,1))
  grid(nx=12,ny=10)      
  legend("topright",legend=c(list_fig_MSSA[[i]][[1]],list_fig_MSSA[[i]][[2]]), cex=0.9, fill=c("blue","red"),bty="n")
  box()
  dev.off()
  
  plot_name<- paste(list_fig_MSSA[[i]][[3]],"_",out_prefix,"shading.png",sep="")
  png(plot_name)
  barplot(heights,     #data to plot
          #main=paste(names(data_plot)[pos1]," and ",names(data_plot)[pos2],sep=""),
          names.arg=names_ind,cex.names=0.8,   #names of the teleconnections indices and size of fonts of axis labes
          beside=TRUE,                         # see two barplots for comparisons...
          xlab="Teleconnection Indices",       # font.lab is 2 to make the font bold
          ylab="Lag crosscorrelation",font.lab=2,
          col=c("blue","red"),
          #density=c(0,10),
          #add=TRUE,
          ylim=c(-1,1))
  grid(nx=12,ny=10)
  legend("topright",legend=c(list_fig_MSSA[[i]][[1]],list_fig_MSSA[[i]][[2]]), cex=0.9, fill=c("blue","red"),bty="n")
  #legend("topright",legend=c(list_fig_MEOT[[i]][[1]],list_fig_MEOT[[i]][[2]]), cex=0.8, col=c("blue","red"),
  #      pch=15)
  barplot(heights,     #data to plot
          #main=paste(names(data_plot)[pos1]," and ",names(data_plot)[pos2],sep=""),
          names.arg=names_ind,cex.names=0.8,   #names of the teleconnections indices and size of fonts of axis labes
          beside=TRUE,                         # see two barplots for comparisons...
          xlab="Teleconnection Indices",       # font.lab is 2 to make the font bold
          ylab="Lag crosscorrelation",font.lab=2,
          col=c("blue","black"),
          density=c(0,10),
          add=TRUE,
          ylim=c(-1,1))
  #grid(nx=12,ny=10)
  legend("topright",legend=c(list_fig_MSSA[[i]][[1]],list_fig_MSSA[[i]][[2]]), cex=0.9, fill=c("blue","red"),bty="n")
  
}

##################### CREATING MEOT SPATIAL SEQUENCE-PLOT ON 12/09/2012 ON ATLAS: Revised fig for paper...  #################
################# PLOT EEOT1
# On Mac Benoit
path_data<-"/Users/benoitparmentier/Dropbox/Data/MEOT_paper/MSSA_paper/Data_paper/MEOT_working_dir_10232012/MSSA_EEOT_04_29_09"
# On Atlas:
path_data<-"/home/parmentier/Data/MEOT12272012/MEOT_working_dir_10232012/MSSA_EEOT_04_29_09"
setwd(path_data)
lf_list<-vector("list",6)
lf_list[[1]]<-list.files(pattern="lag.*_sst_anom_LM_Partial_R_1.rst")
lf_list[[2]]<-list.files(pattern="lag.*_sst_anom_LM_Partial_R_3.rst")
lf_list[[3]]<-list.files(pattern="lag.*_sst_anom_LM_Partial_R_7.rst")
lf_list[[4]]<-list.files(pattern="lag.*_sst_anom_LM_Partial_R_16.rst")
lf_list[[5]]<-list.files(pattern="lag.*_sst_anom_LM_Partial_R_10.rst")
lf_list[[6]]<-list.files(pattern="lag.*_sst_anom_LM_Partial_R_15.rst")
#lf_list[[6]]<-list.files(pattern="lag.*_sst_anom_LM_Partial_R_4.rst")
#list_meot[[]]
#col.breaks <- pretty(s.range, n=100)
#lab.breaks <- pretty(s.range, n=5)

temp.colors <- colorRampPalette(c('blue', 'lightgoldenrodyellow', 'red'))

X11(width=24,height=12)   
#X11(width=24,height=24)
meot_names<-c("MEOT1","MEOT3","MEOT7","MEOT16","MEOT10","MEOT15")
mask_land<-raster("mask_rgf_1_1.rst")
mask_land[mask_land==0]<-NA
#NOT WORKING IN LOOP BECAUSE IT TAKES TOO MUCH TIME...
#X11(width=30,height=9)    
#MEOT 1 and MEOT3

list_fig_MEOT<-vector("list",3)
list_fig_MEOT[[1]]<- c("MEOT1","MEOT3","Figure_3_paper_MEOT1_MEOT3_sequence_spatial_pattern")
list_fig_MEOT[[2]]<- c("MEOT7","MEOT16","Figure_6_paper_MEOT7_MEOT16_sequence_spatial_pattern")
list_fig_MEOT[[3]]<- c("MEOT10","MEOT15","Figure_9_paper_MEOT10_MEOT15_sequence_spatial_pattern")
setwd(path_data)
out_prefix<-"MEOT_paper_01152012b"
for (j in 1:length(lf_list)){
  j=j+1
  #Sys.sleep(.3) #added to control plot time
  plot_name<-meot_names[j]
  #pdf(width=24,height=12,onefile=F,file=paste(plot_name,"test_",out_prefix,".pdf", sep=""))
  lf<-mixedsort(lf_list[[j]])  #Use mixedsort instead of "sort" to take into account the 
  meot_rast<-stack(lf)
  layerNames(meot_rast)<-paste("Lag", 0:12,sep=" ")
  title_plot<-paste(meot_names[j], "spatial sequence",sep=" ")
  meot_rast_m<-mask(meot_rast,mask_land)
  #levelplot(meot_rast_m,col.regions=temp.colors,par.settings = list(axis.text = list(font = 2, cex = 1)))
  levelplot(meot_rast_m,main=title_plot, ylab=NULL,xlab=NULL,,par.settings = list(axis.text = list(font = 2, cex = 1.5),
                      par.main.text=list(font=2,cex=2.2),strip.background=list(col="white")),par.strip.text=list(font=2,cex=1.5),
                      col.regions=temp.colors,at=seq(-1,1,by=0.02))
  #dev.off()
  #levelplot(meot_rast_m,col.regions=temp.colors,scales=list(x=list(rot=100, cex=2))) #control the axis not legend
  savePlot(filename=paste(plot_name,"test_",out_prefix,".tiff", sep=""), type="tiff")
  #Sys.sleep(.0) #Once plot drawned and saved, it can be deactivated 
}
dev.off()

#### DO THIS FOR MSSA
#path_data<-"/Users/benoitparmentier/Dropbox/Data/MEOT_paper/MSSA_paper/Data_paper/MEOT_working_dir_10232012/MSSA_EEOT_04_29_09"
# On Atlas:
path_data<-"/home/parmentier/Data/MEOT12272012/MEOT_working_dir_10232012/MSSA_EEOT_04_29_09"
#path_data<-path_raster
setwd(path_data)
lf_list<-vector("list",2)

#lf_list[[1]]<-list.files(pattern="lag.*_sst_anom_LM_Partial_R_1.rst")
lf_list[[1]]<-list.files(pattern="anom_sst_1982_2007_MSSA_Center_Std_mssa.*_S-Mode_CompLoadingImg_1.rst")
lf_list[[2]]<-list.files(pattern="anom_sst_1982_2007_MSSA_Center_Std_mssa.*_S-Mode_CompLoadingImg_3.rst")
#lf_list[[1]]<-list.files(pattern="anom_sst_1982_2007_MSSA_Center_Std_mssa.*_S-Mode_CompLoadingImg_2.rst")
#lf_list[[1]]<-list.files(pattern="anom_sst_1982_2007_MSSA_Center_Std_mssa.*_S-Mode_CompLoadingImg_7.rst")
#lf_list[[2]]<-list.files(pattern="anom_sst_1982_2007_MSSA_Center_Std_mssa.*_S-Mode_CompLoadingImg_8.rst")

temp.colors <- colorRampPalette(c('blue', 'lightgoldenrodyellow', 'red'))

#"anom_sst_1982_2007_MSSA_Center_Std_mssa3_S-Mode_CompLoadingImg_1.rst"
X11(width=24,height=12)    
#quartz(width=24,height=12)    

mask_land<-raster("mask_rgf_1_1.rst")
mask_land[mask_land==0]<-NA

mssa_names<-c("MSSA1","MSSA3")
#out_prefix<-"MEOT_paper_12212012b_"
for (j in 1:length(lf_list)){
  j=j+1
  #Sys.sleep(.3) #added to control plot time
  plot_name<-mssa_names[j]
  #pdf(width=24,height=12,onefile=F,file=paste(plot_name,"test_",out_prefix,".pdf", sep=""))

  lf<-mixedsort(lf_list[[j]])  #Use mixedsort instead of "sort" to take into account the 
  meot_rast<-stack(lf)
  layerNames(meot_rast)<-paste("Lag", 0:12,sep=" ")
  if (mssa_names[j]=="MSSA5"){
    meot_rast<-meot_rast*-1
  }
  meot_rast_m<-mask(meot_rast,mask_land)
  title_plot<-paste(mssa_names[j], "spatial sequence",sep=" ")
  #levelplot(meot_rast_m,col.regions=temp.colors)
  #levelplot(meot_rast_m,col.regions=temp.colors,cex.labels=4)
  levelplot(meot_rast_m,main=title_plot, ylab=NULL,xlab=NULL,,par.settings = list(axis.text = list(font = 2, cex = 1.5),
            par.main.text=list(font=2,cex=2.2),strip.background=list(col="white")),par.strip.text=list(font=2,cex=1.5),
            col.regions=temp.colors,at=seq(-1,1,by=0.02))
  plot_name<-mssa_names[j]
  #dev.off()
  savePlot(paste(plot_name,"test_",out_prefix,".tiff", sep=""), type="tiff")
  #Sys.sleep(.0) #Once plot drawned and saved, it can be deactivated
}
dev.off()

################################ COMBINED FIG FOR TEMPORAL PROFILES FOR LOADINGS ############################
#idx <- seq(as.Date('1982-01-15'), as.Date('2006-12-15'), by='month')
list_fig_MEOT<-vector("list",4)
list_fig_MEOT[[1]]<- c("MEOT1","MEOT3","Figure_4_paper_MEOT1_MEOT3_temporal_sequence_pattern")
list_fig_MEOT[[2]]<- c("MEOT7","MEOT16","Figure_7_paper_MEOT7_MEOT16_temporal_sequence_pattern")
list_fig_MEOT[[3]]<- c("MEOT10","MEOT15","Figure_10_paper_MEOT10_MEOT15_temporal_sequence_pattern")
list_fig_MEOT[[4]]<- c("MSSA1","MSSA3","Figure_13_paper_MSSA1_MSSA3_temporal_sequence_pattern")

MEOT_quadratures<-c("MEOT1","MEOT3","MEOT7","MEOT16","MEOT10","MEOT15","MSSA1","MSSA3")
#MSSA_quadratures<-c("MSSA1","MSSA3")

dat_subset<-subset(dat,select=MEOT_quadratures)

idx <- seq(as.Date('1982-01-15'), as.Date('2007-12-15'), by='12 month')  #Create a date object for labels...
datelabels<-as.character(1:length(idx))
for (i in 1:length(idx)){
  date_proc<-idx[i]
  month<-strftime(date_proc, "%b")          # extract current month of the date being processed
  day<-strftime(date_proc, "%d")
  year<-strftime(date_proc, "%y")  #Use y instead of Y for 2 digits
  datelabels[i]<-paste(year,sep="")
}

#MSSA_quadratures<-c("MSSA1","MSSA3")
#dat_subset<-subset(dat,select=MSSA_quadratures)

X11(width=10,height=14)
par(mfrow=c(2,1))

for (k in 1:4){
  MEOTa<- list_fig_MEOT[[k]][[1]]
  MEOTb<- list_fig_MEOT[[k]][[2]]
  #if (MEOTb=="MSSA3"){
  #  dat_subset[[MEOTb]]<-dat_subset[[MEOTb]]*-1
  #}
  y_range<-range(dat_subset[[MEOTa]],dat_subset[[MEOTb]])
  plot(dat_subset[[MEOTa]],type="l",col="blue",ylim=y_range,axes=FALSE,ylab="MEOT mode",xlab="Time (month)")
  lines(dat_subset[[MEOTb]],tybe="b",lty="dashed",lwd=1.2,col="darkgreen",axes=FALSE)
  breaks_lab<-seq(1,312,by=12)
  axis(side=2)
  #axis(1,at=breaks_lab, labels=datelabels) #reduce number of labels to Jan and June
  axis(side=1,las=1,
       at=breaks_lab,labels=datelabels, cex=1.5) #reduce number of labels to Jan and June
  box()
  if (MEOTb!="MSSA3"){
    legend("topleft",legend=c(MEOTa,MEOTb), cex=0.8, col=c("blue","darkgreen"),
           lty=c(1,2),lwd=2)  #lwd=line width
  }
  if (MEOTb=="MSSA3"){
    legend("topright",legend=c(MEOTa,MEOTb), cex=0.8, col=c("blue","darkgreen"),
           lty=c(1,2),lwd=2)  #lwd=line width
  }

  
  title(paste("Temporal profiles for", MEOTa, "and", MEOTb,sep=" "))
  
  #add second plot
  pos1<-match(MEOTa,names(d_z))
  pos2<-match(MEOTb,names(d_z))
  #if (MEOTb=="MSSA3"){
  #  ccf_obj<-ccf(d_z[,pos1],d_z[,pos2]*-1, lag=13,plot=FALSE)  #Note that MSSA3 needs to be reversed...11/06
  #}
  ccf_obj<-ccf(d_z[,pos1],d_z[,pos2], lag=13,plot=FALSE)  #Note that ccf does not take
  lag_m<-seq(-1*lag_window,lag_window,1)
  ccf_obj$lag[,1,1]<-lag_m  #replacign lag values because continuous
  
  #X11(type="cairo") #Cairo because it is macos...?
  #plot_name<-paste("crosscorrelation", MEOTa, "and", MEOTb,"lag_analysis",sep="_")#replace by list fig naem
  #png(paste(plot_name,"_",out_prefix,".png", sep=""))
  #plot(ccf_obj, main= paste(telindex, "and", mode_n,"lag analysis",sep=" "), ylab="Cross-correlation",
  #     xlab="Lag (month)", ylim=c(-1,1))
  plot(ccf_obj, main= paste(MEOTa, "and", MEOTb,"lag analysis",sep=" "), ylab="Cross-correlation",
       xlab="Lag (month)", ylim=c(-1,1),
       xaxt="n",lwd="2") #xaxt="n" do not display x axis while yaxt="n" means do not display y axis
  label_ccf<-seq(-10,10,by=1)
  label_ccf<-c(-13,label_ccf,13)
  axis(1,at=label_ccf,label=label_ccf,cex.axis=1)
  
  #plot(ccf_obj,type="b",axes=false)
  #axis(1,at=lag_m,label=lag_m)
  
  savePlot(paste(list_fig_MEOT[[k]][[3]],"_",out_prefix,".tiff", sep=""), type="tiff")  
  
}


r1 <- raster("MEOT1test_MEOT_paper_12092012b.tiff")
r2 <- raster("MEOT3test_MEOT_paper_12092012b.tiff")
ymin(r2)<-901
ymax(r2)<-180
r3<-merge(r1,r2)
writeRaster(r3,"test.tiff")

dev.off()
#X11(type="cairo
### Add this code...
lag_cor<-ccf_obj$acf[,,1]  #This is to access cross-cor
plot(lag_cor*-1)


#### PLOTTING EXPLAINED VARIANCE ###
setwd(path)
dat_variance<-read.xls(infile3, sheet=2)
mssa_var<-dat_variance[,2]
meot_var<-dat_variance[,3]
x_mode<-dat_variance[,1]

plot_name<-paste("mssa_meot_explained_variance",sep="_")
png(paste(plot_name,"_",out_prefix,".png", sep=""))
#X11()
y_range<-range(c(mssa_var,meot_var))
plot(mssa_var~x_mode,ylim=y_range,type="n",ylab="Variance explained in percentage",
     xlab="Mode of variability",
     data=dat_variance) #Prepare for plotting
lines(mssa_var~x_mode,col="red",lty=2,type="b",data=dat_variance)  #lty=2 is "dashed"
points(mssa_var~x_mode,col="red",pch=22,type="b",data=dat_variance)
lines(meot_var~x_mode,col="blue",lty=1,type="b",data=dat_variance)
points(meot_var~x_mode,col="blue",pch=1,type="b",data=dat_variance)
legend("topright",legend=c("MSSA variance", "MEOT variance"), 
       cex=0.9, col=c("red","blue"),
       lty=c(2,1),pch=c(22,1))
#title(paste("Prediction tmax difference and land cover ",sep=""))
dev.off()
#################################

# barplot(heights, names.arg=names_ind, axes=FALSE, axisnames=FALSE,
#         beside=TRUE,
#         col=c("blue","red"), ylim=c(-1,1))
# axis(1,at=1:length(names_ind),labels=names_ind,cex.axis=1)
# axis(2)
# mp<-barplot(heights,axes= FALSE,axisnames= FALSE,
#         beside=TRUE,
#         col=c("blue","red"), ylim=c(-1,1))


################# END OF SCRIPT #######################

# 
# 
# for (j in 1:length(lf_list)){
#   j=j+1
#   lf<-mixedsort(lf_list[[j]])  #Use mixedsort instead of "sort" to take into account the 
#   meot_rast<-stack(lf)
#   layerNames(meot_rast)<-paste("Lag", 0:12,sep=" ")
#   meot_rast_m<-mask(meot_rast,mask_land)
#   #levelplot(meot_rast_m,col.regions=temp.colors)
#   levelplot(meot_rast_m,col.regions=temp.colors,cex.labels=2)
#   
#   plot_name<-meot_names[j]
#   savePlot(paste(plot_name,"test_",out_prefix,".tiff", sep=""), type="tiff")
#   
# }
# dev.off()
#dev.copy(tiff,paste(plot_name,"_",out_prefix,".tiff", sep=""), res=75,antialias="none")

}
#dev.off()
