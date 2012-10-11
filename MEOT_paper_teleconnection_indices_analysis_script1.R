####################################    MEOT-MSSA PAPER   #######################################
############################      SCRIPT 1- MEOT          #######################################
#This script reads SAODI index and other telectonnection indices...                             #
#Cross-corrrelation and autocorrelation functions are generated to examine results of MEOT      #
#and MSSA analyses generated at Clark Labs.                                                     #
#Note that spatial patterns from MEOT and MSSA components are not analyzed in this script       #                 #
#AUTHOR: Benoit Parmentier                                                                      #
#DATE: 06/19/2012                                                                               #
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

### Parameters and argument

infile1<-"SAODI-01-1854_06-2011_test.asc"             #GHCN shapefile containing variables for modeling 2010                 
infile2<-"SAODI-01-1854_06-2011.csv"                     #List of 10 dates for the regression
#infile2<-"list_365_dates_04212012.txt"
infile3<-"MEOT_MSSA_Telcon_indices_08062012.xlsx"                        #LST dates name

#path<-"/Users/benoitparmentier/Documents/DATA/Benoit/Clark_University/Paper_writings/MSSA_BNP/"
path<-"/Users/benoitparmentier/Documents/DATA/Benoit/Clark_University/Paper_writings/MSSA_BNP/MEOT_analysis_R_10102012"

setwd(path)

out_prefix<-"MEOT_paper_10102012_"
telind<-c("PNA","NAO","TNA","TSA","SAOD","MEI","PDO","AO","AAO","AMM","AMOsm","QBO")
mode_list_MEOT<-c("MEOT1","MEOT3","MEOT7","MEOT10","MEOT15","MEOT16")
mode_list_PCA<-c("MSSA1","MSSA2","MSSA3","MSSA4","MSSA5","MSSA6")    

lag_window<-13

## START OF SCRIPT

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

### START THE ANALYSIS
plot(d_z$MEOT1)
acf(d_z$MEOT1) #THis is working, d_z is a zoo object
pacf(d_z$MEOT1)
ccf(d_z$MEOT1,d_z$MEI)
tmp<-ccf(d_z$MEOT1,d_z$MEI)
lag_m<- -21:21    #Creting a sequence of -21 ro 21
tmp$lag[,1,1]<-lag_m
plot(tmp)

# MEOT1 and MEI

tmp<-ccf(d_z$MEOT1,d_z$MEI, lag=13 ) #Note that ccf does not take
lag_m<- -13:13
tmp$lag[,1,1]<-lag_m  #replacign lag values because continuous
plot(tmp, main="Lag cross-correlation between MEOT1 and MEI")

absext <-max(abs(tmp$acf)) # maximum of the extremum
pos<-match(absext,tmp$acf) #find the position and lag
absext_lag<-tmp$lag[pos,1,1] #This is the lag corresponding to the maximum absolute value

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
         
absext <-max(abs(tmp$acf)) # maximum of the extremum
pos<-match(absext,tmp$acf) #find the position and lag, if NA it means it was negative
if (is.na(pos)) {
  pos<-match(absext*-1,tmp$acf)
} 
absext_lag<-tmp$lag[pos,1,1] #This is the lag corresponding to the maximum absolute value

# MEOT10/MEOT15 AND SAOD INDICES
         
tmp<-ccf(d_z$SAOD,d_z$MEOT10, lag=13)  #Note that ccf does not take
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

tmp<-ccf(d_z$SAOD,d_z$MEOT7, lag=13)  #Note that ccf does not take
lag_m<- -13:13
lag_m<-seq(-1*lag_window,lag_window,1)

tmp$lag[,1,1]<-lag_m  #replacign lag values because continuous
plot(tmp, main="SAOD and MEOT7 lag analysis", ylab="Cross-correlation",
     xlab="Lag (month)", ylim=c(-1,1))

absext <-max(abs(tmp$acf)) # maximum of the extremum
pos<-match(absext,tmp$acf) #find the position and lag, if NA it means it was negative
if (is.na(pos)) {
  pos<-match(absext*-1,tmp$acf)
} 
absext_lag<-tmp$lag[pos,1,1] #This is the lag corresponding to the maximum absolute value

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

### PREPARE TABLE FOR THE PAPER...

#telind<-c("PNA","NAO","TNA","TSA","SAOD","MEI","PDO","AO","AAO","AMM","AMOunsm","QBO")
mode_list<-mode_list_MEOT    

lag_table_ext<-matrix(data=NA,nrow=length(telind),ncol=length(mode_list))
lag_table_lag<-matrix(data=NA,nrow=length(telind),ncol=length(mode_list))
lag_table_text<-matrix(data=NA,nrow=length(telind),ncol=length(mode_list))

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
write.table(lag_table_ext,file="lag_table_extremum_window_13.txt",sep=",")

lag_table_lag<-as.data.frame(lag_table_lag)
names(lag_table_lag)<-mode_list
rownames(lag_table_lag)<-telind
write.table(lag_table_lag,file="lag_table_extremum_window_13.txt",sep=",")

lag_table_text<-as.data.frame(lag_table_text)
names(lag_table_text)<-mode_list
rownames(lag_table_text)<-telind
write.table(lag_table_text,file="test.txt",sep=",")

# NOW RUN ANALYSIS FOR PCA
mode_list<-mode_list_PCA
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
write.table(lag_table_ext,file="lag_table_extremum_window_13.txt",sep=",")

lag_table_lag<-as.data.frame(lag_table_lag)
names(lag_table_lag)<-mode_list
rownames(lag_table_lag)<-telind
write.table(lag_table_lag,file="lag_table_extremum_window_13.txt",sep=",")

lag_table_text<-as.data.frame(lag_table_text)
names(lag_table_text)<-mode_list
rownames(lag_table_text)<-telind
write.table(lag_table_text,file="test.txt",sep=",")


################# END OF SCRIPT #######################