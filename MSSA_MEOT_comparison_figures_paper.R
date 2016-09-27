#################################    MEOT AND MSSA Paper  #######################################
########################### SPACE-TIME VARIABILITY  #############################################
#This script examines generates figures to compare MEOT and MSSA results for the paper.
#
#AUTHOR: Benoit Parmentier                                                                       #
#DATE CREATED:09/25/2016 
#DATE MODIFIED: 09/27/2016
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
library(plyr)  
library(lubridate)

#################################################
###### Functions  used in the script  ##########

infile1_function <- file.path("/home/bparmentier/Google Drive/Papers_writing_MEOT/R_scripts/",
                             "PCA_EOT_comparison_data_update_function_09252016.R")
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

out_suffix <- "_figures_mssa_meot_comp_anom_09252016"
create_out_dir_param=TRUE #PARAM8
num_cores <- 4 #PARAM 9

lag_window <- 13
years_to_process <- 1982:2007
var_name <- "sst" #PARAM 14, Name of variable of interest: bacteria measurement (DMR data)
#scaling <- 1/0.0099999998

r_mask_filename <- "/home/bparmentier/Google Drive/Papers_writing_MEOT/EOT_paper/data/lsmasked_0_180.rst"

#in_dir <- "/home/bparmentier/Google Drive/Papers_writing_MEOT/EOT_paper/EOT82_15_July12run/sst_msk_0_180_1982_2015_anom/components"
#in_dir <- "/home/bparmentier/Google Drive/Papers_writing_MEOT/MEOT_paper/SST_data_update_1982_2015/MSSA"
in_dir <- "/home/bparmentier/Google Drive/Papers_writing_MEOT/MEOT_paper/SST_data_update_1982_2015/"

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

SST_dir <- "/home/bparmentier/Google Drive/Papers_writing_MEOT/MEOT_paper/data/SST_anom_1982_2007"

mssa_dir1 <- "/home/bparmentier/Google Drive/Papers_writing_MEOT/MEOT_paper/SST_data_update_1982_2015/MSSA/anom_sst_1982_2007/components"

#original data:
#~/Dropbox/Data/MEOT_paper/MSSA_paper/Data_paper/MEOT_working_dir_10232012/MSSA_EEOT_04_29_09

meot_dir1 <- "/home/bparmentier/Google Drive/Papers_writing_MEOT/MEOT_paper/MSSA_EEOT_04_29_09"

no_comp <- 20 #length(comp_str) #number of components in the analysis
list_no_comp <- 1:no_comp

##############
#### Get MEOT and MSSA loadings

mssa_fname1 <- file.path(mssa_dir1,"anom_sst_1982_2007_MSSA_Center_Std_S-Mode_DIF.ods") #old data 82-2007

mssa1_lf <- lapply(1:no_comp,FUN=function(i){pattern_str=paste0(".*.CompLoadingImg_",i,".rst$")
                       ;mixedsort(list.files(path=mssa_dir1,pattern=pattern_str,full.names=T))})

meot1_lf <- lapply(1:no_comp,FUN=function(i){pattern_str=paste0("lag.*_sst_anom_LM_Partial_R_",i,".rst$")
                       ;mixedsort(list.files(path=meot_dir1,pattern=pattern_str,full.names=T))})

lf_sst <- mixedsort(list.files(SST_dir,pattern=".rst$",full.names=T))

##### read in other information

## Contains teleconnection indices assembled
#indices_fname <- "/home/bparmentier/Google Drive/Papers_writing_MEOT/EOT_paper/data/indices_new.xls"
indices_fname <- "/home/bparmentier/Google Drive/Papers_writing_MEOT/MEOT_paper/data/climate_indices_updated_09232016.ods"

infile1<-"SAODI-01-1854_06-2011_test.asc"             #GHCN shapefile containing variables for modeling 2010                 
infile2<-"SAODI-01-1854_06-2011.csv"                     #List of 10 dates for the regression
#infile2<-"list_365_dates_04212012.txt"
#infile3<-"MEOT_MSSA_Telcon_indices_08062012.xlsx"                        #LST dates name
infile3<-"MEOT_MSSA_Telcon_indices_12112012.xlsx"                        #LST dates name

# on bpy50:
path<-"/home/bparmentier/Google Drive/Papers_writing_MEOT/MEOT_analyses_02202015"
#path_raster<-"/Users/benoitparmentier/Dropbox/Data/MEOT_paper/MEOT12272012/MEOT_working_dir_10232012/MSSA_EEOT_04_29_09"

#telind <- c("PNA","NAO","TNA","TSA","SAOD","MEI","PDO","AO","AAO","AMM","AMOsm","QBO")
telind <- c("PNA","NAO","TNA","TSA","SAOD","MEI","PDO","AO","AAO","AMM","AMOunsm","QBO","Nino3.4","Modoki","Glob_LOT")

########################################################
##############  Start of th script  ##############

#########################
#### PART 1: READ IN DATASETS RELATED TO TELECONNECTION AND PREVIOUS LOADINGS


## Use old indices
### Indices used in the previous study

#Importing SAOD index...Note that it was reformatted from the first file in infile1

#dates <-readLines(paste(path,"/",infile2, sep=""))
SAODI<-read.table(paste(path,"/",infile2,sep=""),sep=",", header=TRUE)

#Prepare data to write out in a textfile
s_SAODI2 <- subset(SAODI, year>1981 & year<2008) #Subset rows that correspond to the conditions
s_SAODI2<- subset(s_SAODI2, select=-year) #Remove the column year
in_SAODI2<-as.vector(t(s_SAODI2))   #This transform the data frame into a one colum
write.table(in_SAODI2,file=paste("SAOD_index_1981_2007",out_suffix,".txt",sep=""),sep=",")

#Prepare data for cross lag correlation analysis and write out results
s_SAODI <- subset(SAODI, year>1981 & year<2007) #Subset rows that correspond to the conditions
s_SAODI<- subset(s_SAODI, select=-year) #Remove the column year
in_SAODI<-as.vector(t(s_SAODI))   #This transform the data frame into a one colum
write.table(in_SAODI,file=paste("SAOD_index_1981_2006",out_suffix,".txt",sep=""),sep=",")

#Import results from MEOT and MSSA analyses with teleconnection indices
dat<-read.xls(file.path(path,infile3), sheet=1)
tail(dat$Date_label)
dat$SAOD<-in_SAODI  #Adding the SAOD index to all the MEOT/MSSA results in the data frame
write.table(dat,"Telcon_indices_with_SAOD_12112012.txt",sep=",")

#write.table(dat,file.path(dirname(indices_fname),"Telcon_indices_with_SAOD_12112012.txt"),sep=",")
   
## find all 7th of the month between two dates, the last being a 7th.
st <- as.Date("1982-1-15")
en <- as.Date("2015-12-15")
dseq <- seq(st, en, by="month") #Creating monthly date sequence to create a time series from a data frame

dat_old_indices <- dat
dat_old_indices$AMOsm <- as.numeric(as.character(dat_old_indices$AMOsm))

dat_old_indices <- rename(dat_old_indices, c("ElNinoModoki"="Modoki")) #from plyr package

#Import results teleconneciton indices processed by Elsa 

#dat<-read.xls(indices_fname, sheet="1982-2015")
dat<-read_ods(indices_fname, sheet="1982-2015") #read in updated indices
names(dat)[1:2] <- c("year","month") #missing names for first two columns
dat <- rename(dat, c("QBO_30_original"="QBO")) #from plyr package

dat$AMOunsm <- c(dat_old_indices$AMOunsm,rep(NA,length=length(dat$AMOsm)-length(dat_old_indices$AMOsm)))

dat_indices <- subset(dat,select=telind)
#View(dat_indices)

##Not elegant but works for now...
test<-lapply(dat_indices ,FUN=convert_to_numeric)
test<- do.call(cbind,test)
dat_indices  <- as.data.frame(test)
dat_indices[dat_indices==-999] <- NA

#Creating time series objects
d_ts<-ts(data=dat,start=c(1982,1), end=c(2015,12), frequency=12, names=names(dat))
colnames(d_ts)
d_z<-as.zoo(d_ts)  #### THIS IS THE TIME SERIES OBJECT USED LATER ON
d_xts<-as.xts(d_ts)

## find all 7th of the month between two dates, the last being a 7th.
st <- as.Date("1982-1-15")
en <- as.Date("2007-12-15")
dseq1 <- seq(st, en, by="month") #Creating monthly date sequence to create a time series from a data frame

dat1 <- dat_indices[1:300,]

d_z1 <- zoo(dat1,dseq1)

indices_dat_dz <- zoo(dat_indices,dseq) #create zoo object from data.frame and date sequence object
indices_dat_dz1 <- zoo(dat1,dseq1)

### add back from old data

dat_old_indices$Glob_LOT <- dat_indices$Glob_LOT[1:300]
dat_old_indices$AMOsm <- as.numeric(dat_old_indices$AMOsm)

old_indices_dat_dz <- zoo(dat_old_indices,as.Date(dseq[1:300])) #create zoo object from data.frame and date sequence object
old_dat_dz<- subset(old_indices_dat_dz,select=telind)

###########################
### GET DATA: 1982-2015 EOT and PCA S mode

mssa1_obj_dat <- read_comp_results(data_filename=mssa_fname1,sheet_name="scores",dseq[1:300]) #1982-2007: old data
plot(mssa1_obj_dat$dat_dz,main="MSSA1 components for old data in the 1982-2007 time period",xlab="Year (Monthly time steps)")

#### Plot all the teleconnection indices and components loadings from EOT and PCA S mode
plot(indices_dat_dz,main="12 teleconnection indices considered")
plot(old_dat_dz)
plot(indices_dat_dz1)

plot(dat_old_indices$AMOsm ~ dseq[1:300],type="l",ylim=c(-1.5,0.7),main="AMO indices comparison")
#lines(dat1$AMOsm ~ dseq[1:300],type="l",col="red")
#nes(liindices_dat_dz1
lines(dat_old_indices$AMOunsm ~ dseq[1:300],type="l",col="green")
lines(dat1$Glob_LOT ~ dseq[1:300],type="l",col="green")

legend("topleft",c("AMOsm_old","AMOsm_new","AMOunsm_old"),col=c("black","red","green"),bty="n")

###########################
#### PART 1: NOW CARRY OUT CORRELATION BETWEEN COMPONENTS AND INDICES

### Generate (tables) Figures 3 through Figure 6

#Figure 3: cross corr MEOT to find quadratures
#Figure 4: cross corr MSSA to find quadratures
#Figure 5: MEOT-indices to identify/interpret modes 
#Figure 12: MSSA-indices to identify/interpret modes 

#################
### Compare MSSA1 (with old data 1982-2007) with teleconenction
#using old indices to check the results

mode_list <- paste0("comp_",1:no_comp)

comp_dat_dz <- mssa1_obj_dat$dat_dz
old_dat_dz<- subset(old_indices_dat_dz,select=telind)

dat_dz <- merge(old_dat_dz, comp_dat_dz, all = FALSE)

#old_dat_dz <- indices_dat_dcross_lag_telind_mssa1_obj$textz[13:408,]
mode_list <- paste0("comp_",1:20) #selected MSSA
out_suffix_str <- paste0("mssa_old_1982_2007_",out_suffix)

cross_lag_telind_mssa1_obj <- crosscor_lag_analysis_fun(telind,
                                   mode_list,
                                   d_z=dat_dz,
                                   lag_window=lag_window,
                                   fig=F,
                                   out_suffix=out_suffix_str)
cross_lag_telind_mssa1_obj$extremum #use this output to examine results?
cross_lag_telind_mssa1_obj$text

### Compare MSSA1 with MSSA1 and in cross-correlation to find quadratures
#using old indices to check the results

mode_list <- paste0("comp_",1:no_comp)

comp_dat_dz <- mssa1_obj_dat$dat_dz
#old_dat_dz <- subset(old_indices_dat_dz,select=telind)

#dat_dz <- merge(comp_dat_dz, comp_dat_dz, all = FALSE)
#old_dat_dz <- indices_dat_dz[13:408,]
mode_list <- paste0("comp_",1:20) #selected MSSA
out_suffix_str <- paste0("cross_cor_mssa1_mssa1_old_1982_2007_",out_suffix)

#debug(crosscor_lag_analysis_fun)
cross_lag_corr_mssa1_mssa1_obj <- crosscor_lag_analysis_fun(mode_list,
                                   mode_list,
                                   d_z=dat_dz,
                                   lag_window=lag_window,
                                   fig=F,
                                   out_suffix=out_suffix_str)

cross_lag_corr_mssa1_mssa1_obj$extremum #use this output to examine results?
cross_lag_telind_mssa1_obj$text

test_extremum <- as.data.frame((cross_lag_corr_mssa1_mssa1_obj$extremum))
test_extremum2  <- ifelse(test_extremum >= 0.5, 1, 0) #find quadrature
colSums((test_extremum2))

#write.table(cor_eot_pca_df ,file=file.path(out_dir,paste0("correlation_between_eot_and_pca_cor_eot_df","_",out_suffix,".txt")),sep=",")

#debug(generate_barplot_fun)
##Comparing cross-correlation between MSSA1 and telind as well as MEOT and telind

################
#### MEOT OLD DATA

### Correlation between MSSA1 (1982-2007 old data) and MEOT (1982-2015 new data)

mode_list <- paste0("MEOT",1:20)
out_suffix_str <- paste0("meot_old_1982_2007_",out_suffix)
cross_lag_telind_meot_obj <- crosscor_lag_analysis_fun(telind,
                                  mode_list,
                                  d_z=old_indices_dat_dz,
                                  lag_window=lag_window,
                                  fig=F,
                                  out_suffix=out_suffix_str)

cross_lag_telind_meot_obj$extremum

### Compare MEOT1 with MEOT1 and in cross-correlation to find quadratures
#using old indices to check the results

mode_list <- paste0("MEOT",1:20)
out_suffix_str <- paste0("cross_cor_meot1_meot1_old_1982_2007_",out_suffix)

#debug(crosscor_lag_analysis_fun)
cross_lag_corr_meot1_meot1_obj <- crosscor_lag_analysis_fun(mode_list,
                                   mode_list,
                                   d_z=old_indices_dat_dz,
                                   lag_window=lag_window,
                                   fig=F,
                                   out_suffix=out_suffix_str)

cross_lag_corr_meot1_meot1_obj$extremum #use this output to examine results?
cross_lag_corr_meot1_meot1_obj$text

test_extremum <- as.data.frame((cross_lag_corr_meot1_meot1_obj$extremum ))
test_extremum2  <- ifelse(test_extremum >= 0.5, 1, 0) #find quadrature

#################### PART 3: Generate spatial patterns figures ##############


#Figure 6. Temporal map sequences for MEOT1 (top) and MEOT3 (bottom) illustrate 
#the unfolding of a typical ENSO event in series of a 13 successive phases. 
#MEOT3 shows the buildup to El Nino while MEOT1 shows strengthening of El Nino conditions.

###############
#### Plot spatial pattern using score maps

#### MEOT
r_mask <- raster(r_mask_filename)
plot(r_mask)

#r_SST <- stack(lf_sst)
#r_SST_m  <- mask(r_SST,mask=r_mask)
#plot(r_SST_m,y=1:12)
#levelplot(r_SST_m)

infile1_function <- file.path("/home/bparmentier/Google Drive/Papers_writing_MEOT/R_scripts/",
                              "PCA_EOT_comparison_data_update_function_09252016.R")
source(infile1_function)

#### MEOT 1: 1982-2007
#temp.colors <- colorRampPalette(c('blue', 'lightgoldenrodyellow', 'red'))
names_stack_lf <- paste0("MEOT",1:no_comp)
i<-1
#lf_eot <- mssa3_lf[[i]]
#undebug(plot_lag_components)
#out_suffix_str <- 
z_range <- c(-1,1)

temp.colors <- colorRampPalette(c('blue', 'lightgoldenrodyellow', 'red'))
title_plot_<-paste(name_stack , "spatial sequence",sep=" ")

plot_lag_components(1,lf=meot1_lf,
                    lag_window = lag_window, 
                    r_mask = r_mask,
                    z_range=z_range,
                    out_suffix=out_suffix,
                    out_dir=out_dir,
                    name_lf=names_stack_lf)
#,
#                    title_plot_str=title_plot)

test_lf <- lapply(1:20,FUN=plot_lag_components,lf=meot1_lf,lag_window = lag_window, r_mask = r_mask,z_range=z_range,out_suffix=out_suffix,out_dir=out_dir,name_lf=names_stack_lf)
#/home/bparmentier/Dropbox/Data/MEOT_paper/MSSA_paper/Data_paper/MEOT_working_dir_10232012/MEOT10232011/anom_sst_1982_2007/components

#### MSSA 1: 1982-2007

names_stack_lf <- paste0("MSSA",1:no_comp)
i<-1
#lf_eot <- mssa3_lf[[i]]
#undebug(plot_lag_components)
#out_suffix_str <- 
z_range <- c(-.3,0.3)
plot_lag_components(1,lf=mssa1_lf,lag_window = lag_window, r_mask = r_mask,z_range=z_range,out_suffix,out_dir,name_lf=names_stack_lf)

test_lf <- lapply(1:20,FUN=plot_lag_components,lf=mssa1_lf,lag_window = lag_window, r_mask = r_mask,z_range=z_range,out_suffix=out_suffix,out_dir=out_dir,name_lf=names_stack_lf)


#################### PART 4: Generate temporal profiles patterns figures ##############

#idx <- seq(as.Date('1982-01-15'), as.Date('2006-12-15'), by='month')
list_fig_MEOT<-vector("list",4)
list_fig_MEOT[[1]]<- c("MEOT1","MEOT3","Figure_4_paper_MEOT1_MEOT3_temporal_sequence_pattern")
list_fig_MEOT[[2]]<- c("MEOT7","MEOT16","Figure_7_paper_MEOT7_MEOT16_temporal_sequence_pattern")
list_fig_MEOT[[3]]<- c("MEOT10","MEOT15","Figure_10_paper_MEOT10_MEOT15_temporal_sequence_pattern")
list_fig_MEOT[[4]]<- c("MSSA1","MSSA3","Figure_13_paper_MSSA1_MSSA3_temporal_sequence_pattern")

##Define the pairs
#MEOT_quadratures<-c("MEOT1","MEOT3","MEOT7","MEOT16","MEOT10","MEOT15","MSSA1","MSSA3")
MEOT_quadratures <- c("MEOT1","MEOT3","MEOT7","MEOT16","MEOT10","MEOT15","MSSA1","MSSA3","MSSA7-MSSA8","MSSA13-MSSA14","MSSA16-MSSA17")

## can also have options for three indices together
MEOT_quadratures <- c("MEOT1,MEOT3",
                      "MEOT7,MEOT16",
                      "MEOT10,MEOT15",
                      "MSSA1,MSSA3",
                      "MSSA7,MSSA8",
                      "MSSA1,MSSA14",
                      "MSSA16,MSSA17")


dat <- old_indices_dat_dz
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

list_quadratures <- strsplit(MEOT_quadratures,",")
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

#################### PART 5: Generate barplots of cross-correlation figures ##############

### Old data: meot and mssa comparison of cross corr with teleconnection indices

### Select quadrature and generate the cross-corr

#Do this for MEOT1-MEOT3,MEOT7-MEOT16,MEOT10-MEOT15

#subset()
out_suffix_str <- paste0("meot_mssa_old",out_suffix)
lf_barplot_comparison <- generate_barplot_comparison_fun(
  df1=cross_lag_telind_mssa1_obj$extremum,
  df2=cross_lag_telind_meot_obj$extremum,out_suffix=out_suffix_str,
  col_palette=NULL,out_dir=NULL)

#Do this for MSSA1-MSSA3,MSSA7-MSSA8,MSSA13-MSSA14,MSSA16-MSSA17

out_suffix_str <- paste0("meot_mssa_old",out_suffix)
lf_barplot_comparison <- generate_barplot_comparison_fun(
  df1=cross_lag_telind_mssa1_obj$extremum,
  df2=cross_lag_telind_meot_obj$extremum,out_suffix=out_suffix_str,
  col_palette=NULL,out_dir=NULL)

########################## END OF SCRIPT #########################

