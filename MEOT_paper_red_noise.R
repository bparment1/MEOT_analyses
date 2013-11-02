####################################    MEOT-MSSA PAPER   #######################################
############################      SCRIPT 3- MEOT          #######################################
#This script reads SAODI index and other telectonnection indices...                             #
#It examines indices against background red noise.                                                   #
#Note that spatial patterns from MEOT and MSSA components are not analyzed in this script       #                 #
#AUTHOR: Benoit Parmentier                                                                      #
#DATE: 11/02/2013           
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
library(forecast)                            # package containing ARIMA procedures
library(rasterVis)
library(plotrix)
library(matrixStats) 
### Parameters and argument

infile1<-"SAODI-01-1854_06-2011_test.asc"             #GHCN shapefile containing variables for modeling 2010                 
infile2<-"SAODI-01-1854_06-2011.csv"                     #List of 10 dates for the regression
#infile2<-"list_365_dates_04212012.txt"
infile3<-"MEOT_MSSA_Telcon_indices_08062012.xlsx"                        #LST dates name
infile3<-"MEOT_MSSA_Telcon_indices_12112012.xlsx"                        #LST dates name

#path<-"/Users/benoitparmentier/Documents/DATA/Benoit/Clark_University/Paper_writings/MSSA_BNP/"
#on MAC:
path<-"/Users/benoitparmentier/Dropbox/Data/MEOT_paper/MEOT12272012/MEOT_working_dir_10232012/MEOT_analysis_R_10102012"
#path<-"/Users/benoitparmentier/Documents/DATA/Benoit/Clark_University/Paper_writings/MSSA_BNP/MEOT_analysis_R_10102012"
# on Atlas:
#path<-"/home/parmentier/Data/MEOT12272012/MEOT_working_dir_10232012/MEOT_analysis_R_10102012"
#path_raster<-"/Users/benoitparmentier/Dropbox/Data/MEOT_paper/MEOT12272012/MEOT_working_dir_10232012/MSSA_EEOT_04_29_09"


setwd(path)

out_prefix<-"MEOT_paper_11012013b_"
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

# correlation amongh all telectonnections...

########## CREATING RED  NOISE ###########
no_sim<- 100
no_t<-nrow(dat)
ar1_coeff<-0.5

test_AR1 <- arima.sim(n=no_t, list(ar=c(ar1_coeff)))
stest<-spectrum(test_AR1)

ar1_sim_fun<-function(i,no_t,ar1_coeff,ref_ts){
  
  #test_AR1 <- arima.sim(n=no_t, list(ar=c(ar1_coef)))
  #d_ts<-ts(data=test_AR1,start=c(1982,1), end=c(2006,12), frequency=12, names="AR1_test")
  #colnames(d_ts)
  #d_AR1<-as.zoo(d_ts)  #### THIS IS THE TIME SERIES OBJECT USED LATER ON
  
  #d_AR1_spec<-spectrum(d_AR1,plot=FALSE)
  
  ts_AR1 <- arima.sim(n=no_t, list(ar=c(ar1_coeff)))
  s1<-spectrum(ts_AR1,plot=FALSE)
  
  power_spec <- s1$spec
  names(power_spec)<-paste("s_",i,sep="")
  return(as.data.frame(power_spec))
}

l_ar1<-lapply(1:no_sim,FUN=ar1_sim_fun,no_t=no_t,ar1_coeff=ar1_coeff)

lf<-do.call(cbind,l_ar1)
lf$freq <-stest$freq
mean_power<-rowMeans(lf[,1:100])
sd_power<-rowSds(lf[,1:100])

plot(lf$freq,mean_power,type="b")

y <- mean_power
x<- lf$freq
no <- no_sim
y_sd <- sd_power
ciw   <- qt(0.975, no,) * y_sd / sqrt(no)
plotCI(y=y, x=x, uiw=ciw, col="red", main=paste(" Spectrum for REd noise"), barcol="blue", lwd=1,
       ylab="mean power (A^2)", xlab="frequency")
lines(y~x)
#ciw   <- qt(0.975, no) * y_sd / sqrt(no)

#if(i==1){
#  plotCI(y=y, x=x, uiw=ciw, col=col_t[i], main=paste(" Comparison of ",metric_name," in ",mod_name,sep=""), barcol="blue", lwd=1,
#         ylab="RMSE (C)", xlab=xlab_text)
#  lines(y~x, col=col_t[i])

#}else{
#  lines(y~x, col=col_t[i])
#}

hist(as.numeric(lf[1,])) #note that the spectrum is chi-sqqured because it is a variance!!!!!
df<-no_t-1
dchisq(as.numeric(lf[1,1]), df, ncp=0, log = FALSE)


coverage<-0.95
tail <- (1 - coverage)
#df <- spec.obj$df
df<-no_t-1
upper.quantile <- 1 - tail * pchisq(df, df, lower.tail = FALSE)
lower.quantile <- tail * pchisq(df, df)
1/(qchisq(c(upper.quantile, lower.quantile), df)/df)

#Get the CI for each spectrum??


#calculate the diff between ref time series in spectrum and the actual time
#series...

#install.packages("portes")
#BoxPierce
#test power spectrum  in R
arima_mod <- arima(d_z$MEOT1, order = c(1,0,0))
arima_mod <- arima(d_z$MEI, order = c(1,0,0))
arima_mod <- acf(d_z$MEI)
arima_mod <- acf(d_z$MEOT1)
arima_mod <- pacf(d_z$MEOT1)
arima_mod <- pacf(d_z$MEI)

arima_mod <- arima(d_z$MEOT1, order = c(1,0,0))
no_sim<- 100
no_t<-nrow(dat)
ar1_coef<-arima_mod$coef[1]
frequency(d_z$MEOT1) #12

test_AR1 <- arima.sim(n=no_t, list(ar=c(ar1_coef)))
stest<-spectrum(test_AR1)
plot(stest$spec~(1/stest$freq),type="b")

l_ar1<-lapply(1:no_sim,FUN=ar1_sim_fun,no_t=no_t,ar1_coeff=ar1_coef)
lf<-do.call(cbind,l_ar1)
lf$freq <-stest$freq
mean_power<-rowMeans(lf[,1:100])
sd_power<-rowSds(lf[,1:100])

plot(1/lf$freq,mean_power,type="b")
plot(lf$freq,mean_power,type="b")

y <- mean_power
x<- lf$freq
no <- no_sim
y_sd <- sd_power
ciw   <- qt(0.975, no,) * y_sd / sqrt(no)
plotCI(y=y, x=x, uiw=ciw, col="red", main=paste(" Spectrum for REd noise"), barcol="blue", lwd=1,
       ylab="mean power (A^2)", xlab="frequency")
lines(y~x)
#ciw   <- qt(0.975, no) * y_sd / sqrt(no)

#if(i==1){
#  plotCI(y=y, x=x, uiw=ciw, col=col_t[i], main=paste(" Comparison of ",metric_name," in ",mod_name,sep=""), barcol="blue", lwd=1,
#         ylab="RMSE (C)", xlab=xlab_text)
#  lines(y~x, col=col_t[i])

#}else{
#  lines(y~x, col=col_t[i])
#}

hist(as.numeric(lf[1,])) #note that the spectrum is chi-sqqured because it is a variance!!!!!
df<-no_t-1
dchisq(as.numeric(lf[1,1]), df, ncp=0, log = FALSE)

test_AR1 <- arima.sim(n=no_t, list(ar=c(ar1_coef)))
d_ts<-ts(data=test_AR1,start=c(1982,1), end=c(2006,12), frequency=12, names="AR1_test")
colnames(d_ts)
d_AR1<-as.zoo(d_ts)  #### THIS IS THE TIME SERIES OBJECT USED LATER ON

d_AR1_spec<-spectrum(d_AR1,plot=TRUE)
stest<-spectrum(d_z$MEOT1,plot=FALSE)
stest<-spectrum(d_z$MEOT1,plot=TRUE)
stest$freq
lines(stest$spec~stest$freq)
plot(stest$spec~stest$freq,type="h")

plot(d_AR1_spec$spec ~ d_AR1_spec$freq,type="b",col="red")
lines(stest$spec~stest$freq,type="b")

plot(d_AR1_spec$spec ~ 1/d_AR1_spec$freq,type="b",col="red")
lines(stest$spec~1/stest$freq,type="b")

######## PLOTTING ON THE SAME SCALE TEH SPECTRUM

#Fit a model for MEOT1 index
arima_mod <- arima(d_z$MEOT1, order = c(1,0,0))
arima_mod <- arima(as.vector(d_z$MEOT1), order = c(1,0,0))
no_sim<- 10000
no_t<-nrow(dat)
ar1_coef<-arima_mod$coef[1]
frequency(d_z$MEOT1) #12

test_AR1 <- arima.sim(n=no_t, list(ar=c(ar1_coef)))
test_AR1 <- arima.sim(list(ar=c(ar1_coef)),n=no_t,sd=sd(d_z$MEOT1)) #innov=rnorm(n=250, mean=0, sd=0.1)

d_AR1_spec<-spectrum(test_AR1,plot=TRUE)
stelind<-spectrum(as.vector(d_z$MEOT1),plot=FALSE)
#stest<-spectrum(d_z$MEOT1,plot=TRUE)
plot(d_AR1_spec$spec ~ d_AR1_spec$freq,type="b",col="red")
lines(stelind$spec~stelind$freq,type="b")

l_ar1<-lapply(1:no_sim,FUN=ar1_sim_fun,no_t=no_t,ar1_coeff=ar1_coef)
lf<-do.call(cbind,l_ar1)
lf$freq <-stest$freq
mean_power<-rowMeans(lf[,1:no_sim])
sd_power<-rowSds(lf[,1:no_sim])

plot(mean_power ~ lf$freq,type="b",col="red")
lines(stelind$spec~stelind$freq,type="b")

#Now add CI to spectrum
y <- mean_power
x<- lf$freq
no <- no_sim
y_sd <- sd_power

m_dat<- as.data.frame(t(as.matrix(lf[,1:no_sim])))

#var.interval <- function(data, conf.level = 0.95)
l_CI<- lapply(m_dat,FUN=sd.interval,conf.level = 0.95)
ci_bands<-do.call(rbind,l_CI)

#ciw   <- qt(0.975, no,) * y_sd / sqrt(no)
plotCI(y=y, x=x, uiw=ci_bands[,1],liw=ci_bands[,2], col="red", main=paste(" Spectrum for REd noise"), barcol="blue", lwd=1,
       ylab="mean power (A^2)", xlab="frequency")
lines(y~x)
lines(stelind$spec~stelind$freq,type="b")
stest<-spectrum(as.vector(d_z$MEI))
telind<-c("PNA","NAO","TNA","TSA","SAOD","MEI","PDO","AO","AAO","AMM","AMOsm","QBO")
mode_list_MEOT<-c("MEOT1","MEOT3", "MEOT4","MEOT7","MEOT10","MEOT15","MEOT16")

coef_list<-vector("list",length(telind))
for (i in 1:length(telind)){
  telindex<-telind[i]
  pos1<-match(telindex,names(d_z))
  arima_mod <- arima(d_z[,pos1], order = c(1,0,0))
  coef_list[i] <- arima_mod$coef[1]
  names(coef_list)[i] <- telind[i]  
}

unlist(coef_list)

## DO CORRELATION MATRIX FOR ALL INDICES AS WELL AS SPECTRUM TO FIND DOMINANT FREQ


tmp<-ccf(d_z$SAOD,d_z$AMM, lag=13)  #Note that ccf does not take
lag_m<- -13:13
tmp$lag[,1,1]<-lag_m  #replacign lag values because continuous

ys1<-sine_structure_fun(1:300,12,0,2,0)
ys2<-sine_structure_fun(1:300,6,0,2,0)
ys3<-sine_structure_fun(1:300,150,0,2,0)

y_test <- ys1 + ys2 + ys3
cpgram(y_test)
plot(y_test,type="b")
y_test_spec <-spectrum(y_test)
y_test_spec <-spectrum(y_test,log=c("no")) #no log 

periodgram(y_test)
plot(y_test_spec$spec~1/y_test_spec$freq,type="h")
f1<- ((y_test_spec$freq*150*4)/300)
plot(y_test_spec$spec~ f1,type="h")
sine_structure_fun <-function(x,T,phase,a,b){
  #Create sine for a one dimensional series
  #Note that sine function uses radian unit.
  #a=amplitude
  #b=mean or amplitude 0 of the series
  #T= stands for period definition
  #phase=phase angle (in radian!!)
  
  y <- a*sin((x*pi/T)+ phase) + b
}

sd.interval <- function(data, conf.level = 0.95) {
       df = length(data) - 1
       chilower = qchisq((1 - conf.level)/2, df)
       chiupper = qchisq((1 - conf.level)/2, df, lower.tail = FALSE)
       v = sd(data)
       ci_range <- as.data.frame(cbind(df * v/chiupper, df * v/chilower)) 
}

#http://svn.r-project.org/R/trunk/src/library/stats/R/spectrum.R
spec.ci <- function (spec.obj, coverage = 0.95)
{
  ## A utility function for plot.spec which calculates the confidence
  ## interval (centred around zero). We use a conditional argument to
  ## ensure that the ci always contains zero.
  
  if (coverage < 0 || coverage >= 1)
    stop("coverage probability out of range [0,1)")
  tail <- (1 - coverage)
  df <- spec.obj$df
  upper.quantile <- 1 - tail * pchisq(df, df, lower.tail = FALSE)
  lower.quantile <- tail * pchisq(df, df)
  1/(qchisq(c(upper.quantile, lower.quantile), df)/df)
}

