########################################  Generation of space time dataset #######################################
########################################### Read QC flags in R for MODIS #####################################
#This script generates synthetic space time data set in R. 
#AUTHORS: Neeti Neeti and Benoit Parmentier                                                                       
#DATE: 08/15/2013                                                                                
#PROJECT: Clark Labs, MEOT, time series                                  
###################################################################################################

###Loading R library and packages                                                      
#library(gtools)                                        # loading some useful tools 

library(sp)
library(raster)
library(rgdal)
library(vegan) 

### Functions

##Function to generate spatial image

generate_raster_region <-function(j,list_param){
  #This....
  #Arguments: no of image, no of rows and columns, projection, range of value to be used, distribution, output prefix, and output directory
  #Output: raster image in IDRISI format
  
  ## Parse input arguments
  
  j <- list_param$j
  col <- list_param$col
  row <- list_param$row
  proj_str <- list_param$proj_str
  range_val <- list_param$range_val
  distribution_val_i <- list_param$distribution_val_i
  out_suffix <- list_param$out_suffix
  out_dir <-list_param$out_dir
  
  ## Start #
  
  if(distribution_val_i=="none"){
    pix_values <- rep(1,row*col)
  }
  if(distribution_val_i=="unif"){
    pix_values <- runif(row*col)
  }
  
  r<-raster(nrows=row, ncols=col,crs=proj_str,xmn=1, xmx = col, ymn=1, ymx=row)
  r <- setValues(r,pix_values)
  
  #Write out the raster defining the region
  raster_name<-paste("reg_raster","_",row, "r_", col,"c_",out_suffix, sep="")
  writeRaster(r, filename=raster_name,NAflag=-999,bylayer=FALSE,bandorder="BSQ",overwrite=TRUE)  #Writing the data in a raster file format...  
  #return raster defining the extent and basis layer!!!
  return(r)
}

#### Parameters and arguments

#in_dir<- "/Users/benoitparmentier/Dropbox/MEOT"
in_dir <- "/Users/benoitparmentier/Dropbox/Data/MEOT_paper/synthetic_dataset_MEOT_revisions"
out_dir<- "/Users/benoitparmentier/Dropbox/Data/MEOT_paper/synthetic_dataset_MEOT_revisions"
setwd(out_dir)

col <- 12 
row <- 12
proj_str <- NA
range_val <- c(0,1)
distribution_val_i <-"none" #c("unif","norm","none") 
#distribution_val_i <-"unif" #c("unif","norm","none") 
out_suffix <- "gen_08112013.rst"
j <- 1 #number of images

############# BEGIN SCRIPT #############

#step 1: generate raster image/region

list_param_gen_raster<-list(j,col,row,proj_str,range_val,distribution_val_i,out_suffix, out_dir)
names(list_param_gen_raster)<-c("j","col","row","proj_str","range_val","distribution_val_i","out_suffix","out_dir")

r <- generate_raster_region(1,list_param_gen_raster)
plot(r)

#step 2: generate spatial structure

#1. constant structure
#2. spatial trend
#3. unipolar spatial mode
#4. bipolar spatial mode
#5. moving spatial mode (4 to 6 states)

range_val <- c(1,2)
type_spatialstructure <-c("constant","trend","unipolar", "bipolar", "moving")
out_suffix <- "spstr_08112013.rst"
out_dir <- "/Users/benoitparmentier/Dropbox/MEOT"
j <- 1 #number of images

list_param_gen_raster<-list(j,col,row,proj_str,range_val,type_spatialstructure,out_suffix, out_dir)
names(list_param_gen_raster)<-c("j","col","row","proj_str","range_val","type_spatialstructure","out_suffix","out_dir")

#Function to generate spatial structure
adding_sptatial_structure  <-function(j,list_param){
  #This....
  #Arguments: no of image, no of rows and columns, projection, range of value to be used, type of spatial pattern, output prefix, and output directory
  #Output: raster image in IDRISI format  
  ## Parse input arguments
  
  j <- list_param$j
  col <- list_param$col
  row <- list_param$row
  proj_str <- list_param$proj_str
  range_val <- list_param$range_val
  type_spatialstructure <- list_param$type_spatialstructure
  out_suffix <- list_param$out_suffix
  out_dir <-list_param$out_dir
    
  ### Start #####
  
  spstr<-raster(nrows=row, ncols=col,crs=proj_str,xmn=1, xmx = col, ymn=1, ymx=row)
 
 if(type_spatialstructure=="constant"){
    pix_values <- min(range_val)
    spstr <- setValues(spstr,pix_values)
    spstr[1:(row/2), 1:(col/2)] <- max(range_val)
    spstr[(row/2)+1:row, (col/2)+1:col] <- max(range_val)
 } 
 if(type_spatialstructure=="trend"){
    pix_values <- min(range_val)
    spstr <- setValues(spstr,pix_values)
    x.coord <- rep(1:col, each=col)
    y.coord <- rep(1:row, times=row)
    xy <- data.frame(x.coord, y.coord)
    xy.dist <- dist(xy)
    pcnm.axes <- pcnm(xy.dist)$vectors
    z.value <- pcnm.axes[,1]*100 + rnorm(row*col, min(range_val), max(range_val))
    spstr[] <- z.value
 }
  
 if(type_spatialstructure=="trend"){
    pix_values <- min(range_val)
    spstr <- setValues(spstr,pix_values)
    x.coord <- rep(1:col, each=col)
    y.coord <- rep(1:row, times=row)
    xy <- data.frame(x.coord, y.coord)
    xy.dist <- dist(xy)
    pcnm.axes <- pcnm(xy.dist)$vectors
    z.value <- pcnm.axes[,1]*100 + rnorm(row*col, min(range_val), max(range_val))
    spstr[] <- z.value
  } 
  
 if(type_Spatialstructure%in%c("trend_h","trend_v"){
   u <-xFromCol(r,colnr=1:col)
   a<-2
   b<- -1
   ux <-  a*u + b  
   ux <-rep(ux,time=row)  
   r1<-setValues(r,ux)  
   if(type_Spatialstructure==("trend_v"){
     r1 <-t(r1)
   }
 }
    
 if(type_Spatialstructure%in%c("periodic_h","periodic_v"){
   u <-xFromCol(r,colnr=1:col)
   a<-2
   b<- 0
   f<-0.125
   ux <-  a*sin(u*pi/(f*col)) + b
   ux <-rep(ux,time=row)  
   r1<-setValues(r,ux)  
   if(type_Spatialstructure==("trend_v"){
      r1 <-t(r1)
   }
 }
       
 raster_name<-paste(raster,"-",row, "r_", col,"c_",out_suffix, sep="")
 writeRaster(r, filename=raster_name,NAflag=-999,bylayer=FALSE,bandorder="BSQ",overwrite=TRUE)
 return()
}

#Generate moving structure?


#Step 3: generate temporal structure

#Function ....
adding_temporal_structure <-function(...)
{...}

nt <- 312 #month units conventially
#temp_index <-""...
temp_period <- c(12) #

#1. constant structure 
#2. temporal trend
#3. period
#4. temp_index : can be ENSO, random time series, red noise etc...
#5. moving in time: two sine functions in quadratures

#Function to combine time and space structure...
combine_space_time_structure <-function(...)
{...}
#multiply at every pixel by temporal structure...
#if more than 1 raster image, muliply in sequence matching time step.

#1.  
#2. temporal trend
#3. period
#4. temp_index : can be ENSO, random time series, red noise etc...
#5. moving in time: two sine functions in quadratures

########## END OF SCRIPT ##############