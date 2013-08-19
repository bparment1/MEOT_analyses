########################################       Generation of space-time dataset       #######################################
########################################### For Testing MSSA-MEOT and S-T mode PCA #####################################
#This script generates synthetic space time data set in R.
#Output format can be chosen between IDRISI raster and tif.
#...
#AUTHORS:Benoit Parmentier and Neeti Neeti                                                 
#DATE: 08/19/2013                                                                                
#PROJECT: Clark Labs, MEOT, time series                                  
###################################################################################################

###Loading R library and packages                                                      
#library(gtools)                                        # loading some useful tools 

library(sp)
library(raster)
library(rasterVis)
library(rgdal)
library(vegan) 

##### Functions used in this script #####

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

#Function to generate spatial structure
adding_sptatial_structure  <-function(j,list_param){
  
  #This....
  #Arguments: no of image, no of rows and columns, projection, range of value to be used, type of spatial pattern, output prefix, and output directory
  #Output: raster image in IDRISI format  
  #TO DO: modify for memory efficiency, write out individual raster images if necessary. 
  #Functions used in this code:
  
  sine_structure_fun <-function(x,T,phase,a,b){
    #Create sine for a one dimensional series
    #Note that sine function uses radian unit.
    #a=amplitude
    #b=mean or amplitude 0 of the series
    #T= stands for period definition
    #phase=phase angle (in radian!!)
    
    y <- a*sin((x*pi/T)+ phase) + b
  }
  
  ## Parse input arguments
  
  j <- list_param$j
  r <- list_param$r
  proj_str <- list_param$proj_str
  range_val <- list_param$range_val
  NA_flag_val <- list_param$NA_flag_val
  file_format <- list_param$file_format #".tif" or ".rst" for the time being
  #type_spatialstructure <- list_param$type_spatialstructure
  out_suffix <- list_param$out_suffix
  out_dir <-list_param$out_dir
  r <- list_param$r
  row <- nrow(r)  
  col<- ncol(r)
  
  ### Start #####
  
  type_spatialstructure <- character (length=15) #empty char vector to store names of spatial structures/raster layers 
  #Generate first spatial pattern: constant square   
  spstr<- r  
  type_spatialstructure[1] <- "constant_sqr"
  pix_values <- min(range_val)
  spstr <- setValues(spstr,pix_values)
  spstr[1:(row/2), 1:(col/2)] <- max(range_val)
  spstr[(row/2)+1:row, (col/2)+1:col] <- max(range_val)
  r1 <- spstr 
  
  #Generate spatial pattern 2: sine diagonal   
  type_spatialstructure[2] <- "sine_diag1" 
  pix_values <- min(range_val)
  spstr <- setValues(spstr,pix_values)
  x.coord <- rep(1:col, each=col)
  y.coord <- rep(1:row, times=row)
  xy <- data.frame(x.coord, y.coord)
  xy.dist <- dist(xy)
  pcnm.axes <- pcnm(xy.dist)$vectors
  z.value <- pcnm.axes[,1]*100 + rnorm(row*col, min(range_val), max(range_val))
  spstr[] <- z.value
  r2 <- spstr
  
  #Generate spatial pattern 3: trend in x    
  type_spatialstructure[3] <-"trend_x"
  u <-xFromCol(r,colnr=1:col)
  a<- 2 #slope
  b<- -1 #intercept
  ux <-  a*u + b  
  ux <-rep(ux,time=row)  
  r3<-setValues(r,ux)  
  
  #Generate spatial pattern 4: trend in y    
  type_spatialstructure[4] <-"trend_y"
  r4 <-t(r3)
  
  #Generate spatial pattern 5:     
  type_spatialstructure[5] <- "periodic_x1"
  u <- xFromCol(r,colnr=1:col)
  a<- 2 #amplitude in this case
  b<- 0
  T<- col
  phase <- 0
  ux <- sine_structure_fun(u,T,phase,a,b)
  ux <-rep(ux,time=row)  
  r5<-setValues(r,ux)  
  plot(r5)
  
  #Generate spatial pattern6:  
  type_spatialstructure[6] <- "periodic_y1"
  r6 <-t(r5)
  
  #Generate spatial pattern7:  
  type_spatialstructure[7] <- "periodic_xy1"
  r7 <- r5 + r6
  
  #Generate spatial pattern 8:     
  type_spatialstructure[8] <- "periodic_x2"
  u <- xFromCol(r,colnr=1:col)
  a<- 2 #amplitude in this case
  b<- 0
  T<- round(col/2)
  phase <- 0  
  ux <- sine_structure_fun(u,T,phase,a,b)
  ux <-rep(ux,time=row)  
  r8 <-setValues(r,ux)  

  #Generate spatial pattern 9:  
  type_spatialstructure[9] <- "periodic_y2"
  r9 <- t(r8)
  
  #Generate spatial pattern 10:  
  type_spatialstructure[10] <- "periodic_xy2"
  r10 <- r8 + r9
    
  #Generate spatial pattern 11:     
  type_spatialstructure[11] <- "periodic_x3"
  u <- xFromCol(r,colnr=1:col)
  a<- 2 #amplitude in this case, should be an argument of the function
  b<- 0
  T<- round(col/3)
  phase <- 0  
  ux <- sine_structure_fun(u,T,phase,a,b)
  ux <-rep(ux,time=row)  
  r11 <-setValues(r,ux)  
  
  #Generate spatial pattern 12:  
  type_spatialstructure[12] <- "periodic_y3"
  r12 <- t(r11)
  
  #Generate spatial pattern 13:  
  type_spatialstructure[13] <- "periodic_xy3"
  r13 <- r11 + r12
  
  #Generate spatial pattern 14:  
  type_spatialstructure[14] <- "constant" #This may contain randomness!!
  r14 <-r 
  
  #Generate spatial pattern 15: trend in xy, diagonal   
  type_spatialstructure[15] <-"trend_xy"
  r15 <- r3 + r4
  
  #TO DO later - Generate user defined patter: using a,b,phase and other inputs
  
  #Prepare stack to return object
  r_spat <-stack(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15)
  #names(r_spat) <- type_Spatialstructure
  layerNames(r_spat) <- type_spatialstructure
  
  #Write out spatial patterns
  if(file_format==".rst"){
    for (i in 1:nlayers(r_spat)){
      rast <-subset(r_spat,i)
      raster_name<-paste("raster","_",row, "r_", col,"c_","r_spat_s_",i,"_",
                         type_spatialstructure[i],"_",out_suffix,file_format, sep="")
      writeRaster(rast, NAflag=NA_flag_val,filename=raster_name,bylayer=FALSE,bandorder="BSQ",overwrite=TRUE)
    }
  }
  #now write out stack...
  if(file_format==".tif"){
  raster_name<-paste("raster","_",row, "r_", col,"c_","r_spat_s","_",out_suffix,file_format, sep="")
  writeRaster(r_spat, NAflag=NA_flag_val,filename=raster_name,bylayer=TRUE,bandorder="BSQ",overwrite=TRUE)   
  }

  return(r_spat)
}

#Function ....
adding_temporal_structure <-function(list_param){
  #This....
  #Arguments: no of image, no of rows and columns, projection, range of value to be used, type of spatial pattern, output prefix, and output directory
  #Output: raster image in IDRISI format  
  
  ##Functions used:
  
  temp_structure_fun <-function(x,T,phase){
    y <-sin((x*pi/T)+ phase)
  }
  
  ## Parse input arguments
  #j <- list_param$j
  nt <- list_param$nt
  phase <- list_param$phase
  temp_periods <- list_param$temp_periods
  temp_period_quadrature <- list_param$temp_period_quadrature
  
  ### Start #####
  
  #type_temporalstructure <- character (length=6)
  x<- as.numeric(1:nt)
  
  # Generate temporal signal with periodic signal
  
  y1_list<-vector("list",length=length(temp_periods))
  for (i in 1:length(temp_periods)){
    T<-temp_periods[i]
    y1_list[[i]] <- temp_structure_fun(x,T,phase)
  }
  names(y1_list) <- paste("t_period",temp_periods,sep="_")
  # Generate periodic Quadrature:
  
  T <- temp_period_quadrature
  y2_list<-vector("list",length=length(temp_periods))
  y2_list[[1]] <- temp_structure_fun(x,T,pi/2)
  y2_list[[2]] <- temp_structure_fun(x,T,0)
  
  y2<-y2_list[[1]] + y2_list[[2]]
  
  # Generate temporal trends
  a<-2
  b<- -1
  y3 <-  a*x + b  
  
  # Generate temporal temporal randomness
  y4 <- runif(nt) 
  
  y5<- rnorm(nt)
  
  #Prepare return object
  
  temp_patterns_ts_obj <- list(y1_list,y2,y3,y4,y5)
  names(temp_patterns_ts_obj) <- c("periodic","quadrature","trend","unif","norm") 
  return(temp_patterns_ts_obj)
}

#generate stack given a raster image, this replicates the original image...
generate_stack_raster<-function(r,n_stack){
  pix_vals <-getValues(r)
  list_r<-vector("list",length=n_stack)
  for (i in 1:n_stack){
    list_r[[i]]<-r
  }
  r_stack<-stack(list_r)
}

#Function to combine time and space structure, given a stack of value and temporal index
combine_space_time_structure <-function(temp_pattern,r_stack,out_suffix,NA_flag_val,file_format){
  #modify later for efficiency if needed, this is a small image!!!
  rast_ref <- raster(r_stack,1)
  r_list <-vector("list",length=nlayers(r_stack))
  names(r_list) <- paste("ts",1:nlayers(r_stack),sep="_") #use names in new raster package
  for (k in 1:nlayers(r_stack)){
    pix_vals <-getValues(subset(r_stack,k))
    r_comb <- setValues(rast_ref,pix_vals*temp_pattern[k]) 
    r_list[[k]] <- r_comb
    if(file_format==".rst"){
      raster_name<-paste("raster","_",row, "r_", col,"c_","ts_",k,"_",out_suffix,file_format, sep="")
      writeRaster(r_comb, NAflag=NA_flag_val,filename=raster_name,bylayer=FALSE,bandorder="BSQ",overwrite=TRUE)
    }
  }
  r_temp <- stack(r_list)
  #now write out stack...
  if(file_format==".tif"){
    raster_name<-paste("raster","_",row, "r_", col,"c_","ts_",k,"_",out_suffix,file_format, sep="")
    #writeRaster(r_comb, NAflag=NA_flag_val,filename=raster_name,bylayer=FALSE,bandorder="BSQ",overwrite=TRUE)
    writeRaster(r_temp, NAflag=NA_flag_val,filename=raster_name,bylayer=TRUE,bandorder="BSQ",overwrite=TRUE)   
  }
  return(r_temp)
}

#Add moving spatial pattern function here...
#...
#combine sine in x and y in small window
#Create sequence with position on a circle...(if 12, divide circle in 12 location)
#combine sequence with r stack with specific time event...

##### Parameters and arguments #####

#in_dir<- "/Users/benoitparmentier/Dropbox/MEOT"
in_dir <- "/Users/benoitparmentier/Dropbox/Data/MEOT_paper/synthetic_dataset_MEOT_revisions"
out_dir<- "/Users/benoitparmentier/Dropbox/Data/MEOT_paper/synthetic_dataset_MEOT_revisions"
setwd(out_dir)
range_val <- c(1,2)
#type_spatialstructure <-c("constant","trend","unipolar", "bipolar", "moving")
#out_suffix <- "spstr_08112013"
#out_dir <- "/Users/benoitparmentier/Dropbox/MEOT"
j <- 1 #number of images

col <- 12 
row <- 12
proj_str <- NA
range_val <- c(0,1)
distribution_val_i <-"none" #c("unif","norm","none") 
#distribution_val_i <-"unif" #c("unif","norm","none") 
file_format <-".rst"
NA_flag_val <- -9999
out_suffix <- "gen_08192013"
j <- 1 #number of images

############# BEGIN SCRIPT #############

#step 1: generate raster image/region

list_param_gen_raster<-list(j,col,row,proj_str,range_val,distribution_val_i,out_suffix, out_dir)
names(list_param_gen_raster)<-c("j","col","row","proj_str","range_val","distribution_val_i","out_suffix","out_dir")

r <- generate_raster_region(1,list_param_gen_raster)
plot(r)

#step 2: generate spatial structure

#note that r is the raste image!!
list_param_adding_spatial_structure <-list(j,r,proj_str,range_val,NA_flag_val,file_format,out_suffix, out_dir)
names(list_param_adding_spatial_structure)<-c("j","r","proj_str","range_val","NA_flag_val","file_format","out_suffix","out_dir")

r_spat_s <-adding_sptatial_structure(1,list_param_adding_spatial_structure)

#Generate moving structure?
#Add separate function for this...? yes...coming soon

#Step 3: generate temporal structure

#Also add option for temp_index

nt <- 312 #month units conventially
#temp_index <-""...
temp_periods <- c(12) #
temp_periods <- c(156,78,39,12)
temp_period_quadrature <- c(12)
phase<-pi/2

list_param_adding_temp_str <- list(nt,temp_periods,temp_period_quadrature,phase)
names(list_param_adding_temp_str) <- list("nt","temp_periods","temp_period_quadrature","phase")

#debug(adding_temporal_structure)
temp_patterns_ts_obj <- adding_temporal_structure(list_param_adding_temp_str)
#to-do later-transform object into data.frame for ease of manipulation an visualization!!!
y<-temp_patterns_ts_obj$periodic$t_period_156
plot(y,type="b")
abline(h=0, v=0, col = "gray60")

temp_pattern <- as.numeric(2*y)

#Replicate same image n times into a stack

r_stack <- generate_stack_raster(r,nt)

###############################
### NOW create a few datasets for the MEOT and MSSA/PCA diagnostics...
###############################

###############################
#Create first test dataset:
ts1<-temp_patterns_ts_obj$periodic$t_period_156 #1 full cyle
ts2<-temp_patterns_ts_obj$periodic$t_period_12 #12 months periods...25 cycles
ts3<-temp_patterns_ts_obj$trend #linear  trend  of slope 1

temp_pattern <- 2*ts2 
plot(temp_pattern,type="b")

#Prepare spatial pattern, choose from r_spat_s
r<-1*subset(r_spat_s,"constant")
r_stack<-generate_stack_raster(r,nt)

#Prepare to generate images combining time and space
out_suffix_s <- paste("r_ts1_",out_suffix,sep="")
r_ts1 <- combine_space_time_structure(temp_pattern,r_stack,out_suffix_s,NA_flag_val,file_format)

levelplot((subset(r_ts1,c(1,78,156,234,312))))
plot((subset(r_ts1,c(1,78,156,234,312))))

#########################################
##### Create second test dataset:

#Prepare temporal pattern
ts1<-temp_patterns_ts_obj$periodic$t_period_156 #1 full cyle
ts2<-temp_patterns_ts_obj$periodic$t_period_12 #12 months periods...25 cycles
ts3<-temp_patterns_ts_obj$trend #linear  trend  of slope 1
temp_pattern <- 2*ts1 + 2*ts2 
plot(temp_pattern,type="b")

#Prepare spatial pattern
r<-1*subset(r_spat_s,"constant_sqr")
r_stack<-generate_stack_raster(r,nt)

#Prepare to generate images combining time and space
out_suffix_s <- paste("r_ts2_",out_suffix,sep="")
r_ts2 <- combine_space_time_structure(temp_pattern,r_stack,out_suffix_s,NA_flag_val,file_format)

#Checking the raster time series generated:
levelplot((subset(r_ts2,c(1,78,156,234,312))))
levelplot(subset(r_ts2,c(1,178)))

########################################
##### Create third test dataset:

#Prepare temporal pattern
ts1<-temp_patterns_ts_obj$periodic$t_period_156 #1 full cyle
ts2<-temp_patterns_ts_obj$periodic$t_period_12 #12 months periods...25 cycles
temp_pattern <- 2*ts1 + 2*ts2
plot(temp_pattern,type="b")

#Prepare spatial pattern
r<-2*subset(r_spat_s,"sine_diag1")
r_stack<-generate_stack_raster(r,nt) #Can be modified for efficiency, save to ".tif"

#Prepare to generate images combining time and space
out_suffix_s <- paste("r_ts3_",out_suffix,sep="")
r_ts3 <- combine_space_time_structure(temp_pattern,r_stack,out_suffix_s,NA_flag_val,file_format)
  
levelplot((subset(r_ts3,c(1,78,156,234,312))))
levelplot(subset(r_ts3,c(1,178)))

#######################################
##### Create fourth test dataset:

#Use periodic spatial pattern
#Prepare temporal pattern
ts1<-temp_patterns_ts_obj$periodic$t_period_156 #1 full cyle
ts2<-temp_patterns_ts_obj$periodic$t_period_12 #12 months periods...25 cycles
ts4<- temp_patterns_ts_obj$norm #random normal pattern
temp_pattern <- 2*ts1 + 2*ts2 + 0.25*ts4
plot(temp_pattern,type="b")

#Prepare spatial pattern
r<-2*subset(r_spat_s,"periodic_xy3") # can combine basic spatial structures!!
r_stack<-generate_stack_raster(r,nt)

#Prepare to generate images combining time and space
out_suffix_s <- paste("r_ts4_",out_suffix,sep="")
r_ts4 <- combine_space_time_structure(temp_pattern,r_stack,out_suffix_s,NA_flag_val,file_format)

levelplot((subset(r_ts4,c(1,78,156,234,312))))
levelplot(subset(r_ts4,c(1,178)))

######################################
##### Create fifth test dataset:

#Use moving spatial pattern

### Create rgf for IDRISI to speed up the building of time series...

#out_suffix_s <- paste("r_ts4_",out_suffix,sep="")
#list.files(pattern=paste(".*.",out_suffix,sep=""),full.names=TRUE)

########## END OF SCRIPT ##############