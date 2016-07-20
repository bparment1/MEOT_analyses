#################################    EOT AND PCA ROTATION  #######################################
########################### SPACE-TIME VARIABILITY  #############################################
#This script examines the 1982-2016 dataset downloaded from NOAA:
#http://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.html
#
#The goal is to run EOT, MEOT and PCA on the updated SST dataset using the IDRISI and the R package "remote".
#
#AUTHOR: Benoit Parmentier                                                                       #
#DATE CREATED:07/11/2016 
#DATE MODIFIED: 07/20/2016
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
library(plyr)                              # contains "rename","revalue" and other useful functions

#################################################
###### Functions  used in the script  ##########

infile1_function <- file.path("/home/bparmentier/Google Drive/Papers_writing_MEOT/R_scripts/",
                             "PCA_EOT_comparison_data_update_function_07202016.R")
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

out_suffix <- "eot_pca_1982_2015_anom_07152016"
create_out_dir_param=TRUE #PARAM8
num_cores <- 4 #PARAM 9

years_to_process <- 1982:2015
#start_date <- "2012-01-01" #PARAM 12
#end_date <- "2012-12-31" #PARAM 13 #should process by year!!!
var_name <- "sst" #PARAM 14, Name of variable of interest: bacteria measurement (DMR data)
#scaling <- 1/0.0099999998

r_mask_filename <- "/home/bparmentier/Google Drive/Papers_writing_MEOT/EOT_paper/data/lsmasked_0_180.rst"

in_dir <- "/home/bparmentier/Google Drive/Papers_writing_MEOT/EOT_paper/EOT82_15_July12run/sst_msk_0_180_1982_2015_anom/components"
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

SST_dir <- "SST_1982_2015"
eot_dir <- "/home/bparmentier/Google Drive/Papers_writing_MEOT/EOT_paper/EOT82_15_July12run/sst_msk_0_180_1982_2015_anom/components"
pca_dir <- "/home/bparmentier/Google Drive/Papers_writing_MEOT/EOT_paper/EOT82_15_July12run/sst_msk_0_180_1982_2015_anom/components"

lf_sst <- list.files(path=file.path(in_dir,SST_dir),pattern=".rst$",full.names=T)
### Get EOT and PCA raster images here!!
lf_eot <- mixedsort(list.files(path=eot_dir,pattern="eot_sst_msk_0_180_1982_2015_anom_EOT_Center_Std_.*._R.rst$",full.names=T))
lf_pca <- mixedsort(list.files(path=pca_dir,pattern="sst_msk_0_180_1982_2015_anom_PCA_Center_Std_S-Mode_CompLoadingImg_.*.rst$",full.names=T))

pca_fname1 <- file.path(in_dir,"sst_msk_0_180_1982_2015_anom_PCA_Center_Std_S-Mode_DIF.ods") #ods file with loadings and variance
eot_fname1 <- file.path(in_dir,"eot_sst_msk_0_180_1982_2015_anom_EOT_Center_Std.ods")

## Contains teleconnection indices assembled
indices_fname <- "/home/bparmentier/Google Drive/Papers_writing_MEOT/EOT_paper/data/indices_new.xls"

telind<-c("PNA","NAO","TNA","TSA","SAOD","MEI","PDO","AO","AAO","AMM","AMOsm","QBO")

#lag_window<-13 #no lag window here

########################################################
##############  Start of th script  ##############

#########################
#### PART 0: READ IN DATASETS RELATED TO TELECONNECTION AND PREVIOUS LOADINGS

#Import results teleconneciton indices processed by Elsa 

dat<-read.xls(indices_fname, sheet="1982-2015")
names(dat)[1:2] <- c("year","month")
dat <- rename(dat, c("QBO_30_original"="QBO")) #from plyr package

dat_indices <- subset(dat,select=telind)
View(dat_indices)
     
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
en <- as.Date("2015-12-15")
dseq <- seq(st, en, by="month") #Creating monthly date sequence to create a time series from a data frame
d_z2<-zoo(dat,dseq)
time(d_z)  #no time stamp??
time(d_z2) 

indices_dat_dz <- zoo(dat_indices,dseq) #create zoo object from data.frame and date sequence object

### GET DATA: 1982-2015 EOT and PCA S mode

pca_dat <- read_ods(pca_fname1,sheet="pca_loadings")
names(pca_dat)[1] <- "fnames"
names(pca_dat)[2:ncol(pca_dat)] <- paste0("pc_",1:20)
test <- subset(pca_dat,select=paste0("pc_",1:20))
test<-lapply(test,FUN=convert_to_numeric)
test<- do.call(cbind,test)
pca_dat <- as.data.frame(test)
pca_dat_dz <- zoo(pca_dat,dseq) #create zoo object from data.frame and date sequence object
time(pca_dat_dz)  #no time stamp??

eot_dat <- read_ods(eot_fname1,sheet="eot_loadings")
#test <- read_ods(eot_fname1,sheet="pca_loadings")
test <- subset(eot_dat,select=paste0("eot",1:20))
test<-lapply(test,FUN=convert_to_numeric)
test<- do.call(cbind,test)
eot_dat <- as.data.frame(test)
eot_dat_dz <- zoo(eot_dat,dseq) #create zoo object from data.frame and date sequence object

#### Plot all the teleconnection indices and components loadings from EOT and PCA S mode
plot(indices_dat_dz)
plot(eot_dat_dz)
plot(pca_dat_dz)

###########################
#### PART 1: NOW CARRY OUT CORRELATION BETWEEN INDICES AND COMPONENTS

#cor_series_fun(ts1,ts2,fig=F,out_suffix)
#debug(cor_series_fun)
#infile1_function <- file.path("/home/bparmentier/Google Drive/Papers_writing_MEOT/R_scripts/",
#                              "PCA_EOT_comparison_data_update_function_07192016.R")
#source(infile1_function)

cor_pca_df <- cor_series_fun(ts1=dat_indices,ts2=pca_dat,fig=F,out_suffix)

write.table(cor_pca_df ,file=file.path(out_dir,paste0("cor_pca_df","_",out_suffix,".txt")),sep=",")

cor_eot_df <- cor_series_fun(ts1=dat_indices,ts2=eot_dat,fig=F,out_suffix)
#test <-lapply(cor_eot_df,FUN=convert_to_numeric)
#test<- do.call(cbind,test)
#cor_eot_df <- as.data.frame(test)

write.table(cor_eot_df ,file=file.path(out_dir,paste0("cor_eot_df","_",out_suffix,".txt")),sep=",")

### Correlation between EOT and PCA

cor_eot_pca_df <- cor_series_fun(ts1=pca_dat,ts2=eot_dat,fig=F,out_suffix)

write.table(cor_eot_pca_df ,file=file.path(out_dir,paste0("correlation_between_eot_and_pca_cor_eot_df","_",out_suffix,".txt")),sep=",")

############################
#### PART 2: Analyses and plots of results 

### generate barplot for each components + indices?
### Generate maps+ temporal loadings figures for each (pca and eot!!)
### Generate maps+ barplots?
## Get scree plots from reading variance and compare results!!!
##

###############
#### Plot spatial pattern using score maps

#### EOT
r_mask <- raster(r_mask_filename)
plot(r_mask)

eot_s <- stack(lf_eot)
eot_s <- mask(eot_s,r_mask)
names(eot_s) <- paste0("eot_",1:20)
plot(eot_s,col=matlab.like(256))


res_pix <- 600 #it will be time 2
res_pix <- 500
col_mfrow <- 5
row_mfrow <- 4

png_file_name <- file.path(out_dir,paste0("figure_spatial_pattern_levelplot_eot_stack_","comp_1_20_",out_suffix,".png"))
png(png_file_name,width=col_mfrow*res_pix,height=row_mfrow*res_pix)

#par(mfrow=c(row_mfrow,col_mfrow))
#layout=c(3, 2)
layout_m <- c(col_mfrow,row_mfrow)

temp.colors <- colorRampPalette(c('blue', 'lightgoldenrodyellow', 'red'))
title_plot<-paste("EOT", "spatial sequence",sep=" ")
#levelplot(meot_rast_m,col.regions=temp.colors)
#levelplot(meot_rast_m,col.regions=temp.colors,cex.labels=4)
levelplot(eot_s,main=title_plot, layout=layout_m,
          ylab=NULL,xlab=NULL,
          par.settings = list(axis.text = list(font = 2, cex = 1.5),
                              par.main.text=list(font=2,cex=2.2),strip.background=list(col="white")),par.strip.text=list(font=2,cex=1.5),
          col.regions=temp.colors,at=seq(-1,1,by=0.02))
dev.off()

#### PCA

pca_s <- stack(lf_pca)
pca_s <- mask(pca_s,r_mask)
names(pca_s) <- paste0("pca_",1:20)
plot(pca_s,col=matlab.like(256))
levelplot(pca_s,re=matlab.like(256))

res_pix <- 600 #it will be time 2
res_pix <- 500
col_mfrow <- 5
row_mfrow <- 4

png_file_name <- file.path(out_dir,paste0("figure_spatial_pattern_levelplot_pca_stack_","comp_1_20_",out_suffix,".png"))
png(png_file_name,width=col_mfrow*res_pix,height=row_mfrow*res_pix)

#par(mfrow=c(row_mfrow,col_mfrow))
#layout=c(3, 2)
layout_m <- c(col_mfrow,row_mfrow)

temp.colors <- colorRampPalette(c('blue', 'lightgoldenrodyellow', 'red'))
title_plot<-paste("PCA", "spatial sequence",sep=" ")
#levelplot(meot_rast_m,col.regions=temp.colors)
#levelplot(meot_rast_m,col.regions=temp.colors,cex.labels=4)
levelplot(pca_s,main=title_plot, layout=layout_m,
          ylab=NULL,xlab=NULL,
          par.settings = list(axis.text = list(font = 2, cex = 1.5),
                              par.main.text=list(font=2,cex=2.2),strip.background=list(col="white")),par.strip.text=list(font=2,cex=1.5),
          col.regions=temp.colors,at=seq(-1,1,by=0.02))
dev.off()

## can show there is more autocorrelation in EOTs?
## can show there is less variance in EOTs?

###############
#### Plot temporal profiles patterns using loadings 

### Fix correlation function

#plot_name<- paste(list_fig_MSSA[[i]][[3]],"_",out_prefix,".png",sep="")
#names_ind <- c("PNA","NAO","TNA","TSA","SAOD","MEI","PDO","AO","AAO","AMM","AMO","QBO")
#names_ind <- telind

#data_plot<-cbind(as.vector(mssa_obj$extremum[pos1]),as.vector(mssa_obj$extremum[pos2]))

### Generate for comparison between EOT and PCA??

#debug(generate_barplot_fun)
lf_barplot_comparison <- generate_barplot_comparison_fun(df1=cor_eot_df,df2=cor_pca_df,out_suffix=out_suffix ,col_palette=NULL,out_dir=NULL)

###############
#### Generate maps+ temporal loadings figures for each (pca and eot!!)

#Does it have more regional coherence compared to PCA hence it may work more for predicting specific
#variables locally?

test < - cor(as.matrix(pca_dat),as.matrix(eot_dat),use="na.or.complete")

##############
#### SCREE PLOT AND VARIANCE

#Get scree plots from reading variance and compare results!!!


########################## END OF SCRIPT ###############################################

### Make this a function now:

# i<-1 #right now works on columns only
# names_ref <- names(cor_eot_df)[i]
# names_ind <- rownames(cor_eot_df)
# 
# data_plot <- as.numeric((cor_eot_df[,i])) # make sure it is numeric
# png(plot_name)
# heights<-as.matrix(t(data_plot))
# #barplot(heights)
# barplot(heights,     #data to plot
#         #main=paste(names(data_plot)[pos1]," and ",names(data_plot)[pos2],sep=""),
#         main= names_ref,
#         names.arg=names_ind,cex.names=0.8,   #names of the teleconnections indices and size of fonts of axis labes
#         #beside=TRUE,                         # see two barplots for comparisons...
#         xlab="Teleconnection Indices",       # font.lab is 2 to make the font bold
#         ylab="Correlation",font.lab=2,
#         #col=c("blue","red"), ylim=c(-1,1))
#         col=c("blue"), ylim=c(-1,1))
# grid(nx=12,ny=10)      
# 
# #legend("topright",legend=c(list_fig_MSSA[[i]][[1]],list_fig_MSSA[[i]][[2]]), cex=0.9, fill=c("blue","red"),bty="n")
# legend("topright",legend=c(names_ref), cex=0.9, fill=c("blue"),bty="n")
# 
# box()
# dev.off()
