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

path<-"/Users/benoitparmentier/Documents/DATA/Benoit/Clark_University/Paper_writings/MSSA_BNP"
setwd(path)

## START OF SCRIPT

#Importing SAOD index...

#dates <-readLines(paste(path,"/",infile2, sep=""))
SAODI<-read.table(paste(path,"/",infile1,sep=""),sep=" ", header=FALSE, skip=1) #This does not work
SAODI<-read.table(paste(path,"/",infile1,sep=""),sep=" ", header=FALSE)
SAODI<-read.table(paste(path,"/",infile2,sep=""),sep=",", header=TRUE)

s_SAODI <- subset(SAODI, year>1981 & year<2007)
s_SAODI<- subset(s_SAODI, select=-year)
in_SAODI<-as.vector(t(s_SAODI))   #This transform the data frame into a one colum

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

d_z<-as.zoo(d_ts)
d_xts<-as.xts(d_ts)

## find all 7th of the month between two dates, the last being a 7th.
st <- as.Date("1982-1-15")
en <- as.Date("2006-12-15")
dseq <- seq(st, en, by="month") #Creating monthly date sequence to create a time series from a data frame
d_z2<-zoo(dat,dseq)
time(d_z)  #no time stamp??
time(d_z2) 

plot(d_z2$MEOT1)
plot(d_z$MEOT1)

acf(d_z2$MEOT1)  #not working??
pacf(d_z2$MEOT1)

acf(d_z$MEOT1)
pacf(d_z$MEOT1)
ccf(d_z$MEOT1,d_z$MEI)
tmp<-ccf(d_z$MEOT1,d_z$MEI)
lag_m<- -21:21
tmp$lag[,1,1]<-lag_m
plot(tmp)

# MEOT1 and MEI

tmp<-ccf(d_z$MEOT1,d_z$MEI, lag=13  #Note that ccf does not take
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
plot(tmp, main="Lag cross-correlation between SAOD and AMM")
                  
absext <-max(abs(tmp$acf)) # maximum of the extremum
pos<-match(absext,tmp$acf) #find the position and lag
absext_lag<-tmp$lag[pos,1,1] #This is the lag corresponding to the maximum absolute value
         
# TSA AND SAOD INDICES
         
tmp<-ccf(d_z$SAOD,d_z$TSA, lag=13)  #Note that ccf does not take
lag_m<- -13:13
tmp$lag[,1,1]<-lag_m  #replacign lag values because continuous
plot(tmp, main="Lag cross-correlation between SAOD and AMM")
         
absext <-max(abs(tmp$acf)) # maximum of the extremum
pos<-match(absext,tmp$acf) #find the position and lag
absext_lag<-tmp$lag[pos,1,1] #This is the lag corresponding to the maximum absolute value

# TNA AND SAOD INDICES
         
tmp<-ccf(d_z$SAOD,d_z$TNA, lag=13)  #Note that ccf does not take
lag_m<- -13:13
tmp$lag[,1,1]<-lag_m  #replacign lag values because continuous
plot(tmp, main="Lag cross-correlation between SAOD and TNA")
         
absext <-max(abs(tmp$acf)) # maximum of the extremum
pos<-match(absext,tmp$acf) #find the position and lag, if NA it means it was negative
if (is.na(pos)) {
  pos<-match(absext*-1,tmp$acf)
} 
absext_lag<-tmp$lag[pos,1,1] #This is the lag corresponding to the maximum absolute value
         