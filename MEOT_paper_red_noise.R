####################################    MEOT-MSSA PAPER   #######################################
############################      SCRIPT 3- MEOT          #######################################
#This script reads SAODI index and other telectonnection indices...                            
#It examines indices against background red noise.                                                  
#Note that spatial patterns from MEOT and MSSA components are not analyzed in this script                      
#AUTHOR: Benoit Parmentier                                                                      
#DATE: 11/07/2013           
#Version: 4
#PROJECT: Clark Labs Climate predction- MEOT/MSSA paper                                         
#
#################################################################################################

###Loading R library and packages                                                      
#
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
#
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

lag_window<- 13

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

ar1_sim_fun<-function(i,no_t,ar1_coeff,ref_ts,sd_val=NULL){
  #no_t: length of time series
  #ar1_coeff: coef fitted from AR model
  
  if(is.null(sd_val)){
    ts_AR1 <- arima.sim(n=no_t, list(ar=c(ar1_coeff)))
  }else{
    ts_AR1 <- arima.sim(n=no_t, list(ar=c(ar1_coeff)),sd=sd_val)
  }
  
  s1<-spectrum(ts_AR1,plot=FALSE)
  
  power_spec <- s1$spec
  names(power_spec)<-paste("s_",i,sep="")
  return(as.data.frame(power_spec))
}

sine_structure_fun <-function(x,T,phase,a,b){
  #Create sine for a one dimensional series
  #Note that sine function uses radian unit.
  #a=amplitude
  #b=mean or amplitude 0 of the series
  #T= stands for period definition
  #phase=phase angle (in radian!!)
#  y <- a*sin((x*pi/T)+ phase) + b
  
  y <- a*sin((x*2*pi/T)+ phase) + b
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

######## PLOTTING ON THE SAME SCALE TEH SPECTRUM

#Fit a model for MEOT1 index
arima_mod <- arima(d_z$MEOT1, order = c(1,0,0))
arima_mod <- arima(as.numeric(d_z$MEOT1), order = c(1,0,0))
arima_mod_AMM<- arima(as.vector(d_z$AMM), order = c(1,0,0))

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
set.seed(100)
l_ar1<-lapply(1:no_sim,FUN=ar1_sim_fun,no_t=no_t,ar1_coeff=ar1_coef)
lf<-do.call(cbind,l_ar1)
lf$freq <-d_AR1_spec$freq
mean_power<-rowMeans(lf[,1:no_sim])
sd_power<-rowSds(lf[,1:no_sim])
median_power <- median(m_dat[,1])
hist(log(lf[,1:1]))
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

tt<-t.test(m_dat[,1]-y[1])

##test based on slope??
plot(1:150,as.vector(log(m_dat[1,])))
points(1:150,as.vector(log(m_dat[2,])),col="red")
##test based on slope??
plot(log(1:150),as.vector(log(m_dat[1,])))
points(log(1:150),as.vector(log(m_dat[2,])),col="red")
##test based on slope??
plot(log(1:150,base=10),as.vector(log(m_dat[1,])))
points(log(1:150,base=10),as.vector(log(m_dat[2,])),col="red")

y_f<- as.numeric(log(m_dat[1,]))
x_f<-log(1:150,base=10)
x_f<-lf$freq
d_f <- data.frame(x_f,y_f)
lm_m1<-lm(y_f~x_f,data=d_f)
summary(lm_m1)

get_coef<-function(y_f,x_f){
  d_f <- data.frame(as.numeric(log(x_f)),log(as.numeric(y_f)))
  lm_m1<-lm(y_f~x_f,data=d_f)
  #summary(lm_m1)
  coef(lm_m1)[2]
}

m_dat2 <- as.data.frame(t(as.matrix(m_dat)))
tt<-as.numeric(lapply(m_dat2,FUN=get_coef,x_f=1:150))

hist(as.numeric(tt))
val_telind<-get_coef(as.numeric(stelind$spec),stelind$freq)

m <- mean(tt)
s <- sd(tt)
n <- length(tt)
error <- qt(0.975,df=n-1)*s/sqrt(n)
left <- m-error
right <- m+error
ci<-cbind(left,right)
ci
val_telind-m/s

#pnorm(val_telind-m/s)

### Ok CI should be based on the log transform!!!

m_datlog<-log(m_dat) #log transform all columns...
mean_power_log<-colMeans(m_datlog)
sd_power_log <-colSds(m_datlog)
hist(m_datlog[,1])
plot(mean_power_log ~ loglf$freq),type="b",col="red")

stelind<-spectrum(as.vector(d_z$MEOT1),plot=FALSE)
lines(log(stelind$spec)~(stelind$freq),type="b")
stelind2<-spectrum(as.numeric(d_z$MEOT1)-mean(as.numeric(d_z$MEOT1)),plot=FALSE)
lines(log(stelind2$spec)~log(stelind2$freq),type="b",col="green")

y <- mean_power_log
x<- log(lf$freq)
x<- (lf$freq)
no <- no_sim

ciw   <- qt(0.975, no,) * y_sd / sqrt(no)
plotCI(y=y, x=x, uiw=ciw,liw=ciw, col="red", main=paste(" Spectrum for REd noise"), barcol="blue", lwd=1,
       ylab="log mean power (A^2)", xlab="log frequency",ylim=c(-8,8))
lines(y~x)
lines(log(stelind$spec)~log(stelind$freq),type="b",col="black")
lines(log(stelind$spec)~(stelind$freq),type="b")
lines(log(stelind2$spec)~(stelind2$freq),type="b",col="blue")
#lines(log(stelind$spec)~log(stelind$freq),type="b",col="black")
#plot(log(stelind2$spec)~(stelind2$freq),type="b")
lines(log(stelind2$spec)~(stelind2$freq),type="b")

m_datlog<-log(m_dat)
l_logCI<- lapply(m_datlog,FUN=sd.interval,conf.level = 0.95)

ci_log_bands<-do.call(rbind,l_logCI)


#ciw   <- qt(0.975, no,) * y_sd / sqrt(no)
plotCI(y=y, x=x, uiw=ci_log_bands[,1],liw=ci_log_bands[,2], col="red", main=paste(" Spectrum for REd noise"), barcol="blue", lwd=1,
       ylab="mean power (A^2)", xlab="frequency")
lines(y~x)
lines(log(stelind$spec)~stelind$freq,type="b")

## DO CORRELATION MATRIX FOR ALL INDICES AS WELL AS SPECTRUM TO FIND DOMINANT FREQ

stest<-spectrum(as.vector(d_z$MEI))
telind<-c("PNA","NAO","TNA","TSA","SAOD","MEI","PDO","AO","AAO","AMM","AMOsm","QBO")
mode_list_MEOT<-c("MEOT1","MEOT3", "MEOT4","MEOT7","MEOT10","MEOT15","MEOT16")

mat_telind <- as.matrix(subset(d_z,select=telind))
cor_x<-round(cor(mat_telind),digit=2)

coef_list<-vector("list",length(telind))
for (i in 1:length(telind)){
  telindex<-telind[i]
  pos1<-match(telindex,names(d_z))
  arima_mod <- arima(d_z[,pos1], order = c(1,0,0))
  coef_list[i] <- arima_mod$coef[1]
  names(coef_list)[i] <- telind[i]  
}

unlist(coef_list)
print(cbind(telind,round(as.numeric(coef_list),digit=3)))

d_z_telind <- subset(d_z,stelind) #matix of of correlation between indices

tmp<-ccf(d_z$SAOD,d_z$AMM, lag=13)  #Note that ccf does not take
lag_m<- -13:13
tmp$lag[,1,1]<-lag_m  #replacign lag values because continuous

x<- 1:300
ys1<-sine_structure_fun(x,12,0,2,0) 
ys2<-sine_structure_fun(x,24,0,2,0)
ys3<-sine_structure_fun(x,300,0,2,0) #fundamental harmonic
ys4<- 0.04*x 
ys5<-rnorm(300,sd=0.2)
y_test <- ys1 + ys2 + ys3 + ys4 + ys5

#cpgram(y_test)
plot(y_test,type="b")
y_test_spec <-spectrum(y_test)
y_test_spec <-spectrum(y_test,log=c("no")) #no log 
plot(ts(cbind(ys1,ys2,ys3,ys4,y_test)))
#periodgram(y_test)
plot(y_test_spec$spec~1/y_test_spec$freq,type="h")
f1<- ((y_test_spec$freq*150*4)/300)
plot(y_test_spec$spec~ f1,type="h")

## Simulate red noise:


