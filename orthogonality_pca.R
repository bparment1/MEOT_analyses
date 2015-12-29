#################################    EOT AND PCA ROTATION  #######################################
########################### SPACE-TIME VARIABILITY  #############################################
#This script generates orthogonal datasets in time.
#The goal is to compare EOT and PCA with rotation and without using the synthetic dataset.
#This is part of the Time Series analysis project started at IDRISI.
#
#AUTHOR: Benoit Parmentier                                                                       #
#DATE CREATED: 12/28/2015 
#DATE MODIFIED: 12/28/2015
#
#PROJECT: MEOT/EOT climate variability extraction
#
##################################################################################################
#
###Loading r library and packages

library(raster)                             # loading the raster package
library(gtools)                             # loading ...
library(sp)                                 # spatial objects in R
library(gplots)                             # 
library(rgdal)                              # gdal driver for R
library(RColorBrewer)                       # color scheme, palettes used for plotting
library(gdata)                              # read different format (including .xlsx)
library(plotrix)                            # plot options and functions including plotCI
library(rasterVis)                          # raster visualization
library(colorRamps)                         # contains matlab.like palette
library(lsr)
library(psych)
library(GPArotation)
library(zoo)
library(xts)

#################################################
###### Functions  used in the script  ##########

create_dir_fun <- function(outDir,out_suffix){
  #if out_suffix is not null then append out_suffix string
  if(!is.null(out_suffix)){
    out_name <- paste("output_",out_suffix,sep="")
    outDir <- file.path(outDir,out_name)
  }
  #create if does not exists
  if(!file.exists(outDir)){
    dir.create(outDir)
  }
  return(outDir)
}

load_obj <- function(f){
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

generate_orthogonal_and_correlation_vectors_fun <- function(seed_number=100,n=100,rho=0.5){
  #
  #
  #
  
  #seed_number <- seed_number
  #n     <- n                   # length of vector
  #rho   <- rho                   # desired correlation = cos(angle)
  
  #n     <- 20                    # length of vector
  #rho   <- 0.6                   # desired correlation = cos(angle)
  theta <- acos(rho)             # corresponding angle
  set.seed(0)
  
  x1    <- rnorm(n, 1, 1)        # fixed given data
  set.seed(seed_number)
  
  #x2    <- rnorm(n, 2, 0.5)      # new random data
  x2    <- rnorm(n, 1, 0.5)      # new random data
  X     <- cbind(x1, x2)         # matrix
  Xctr  <- scale(X, center=TRUE, scale=FALSE)   # centered columns (mean 0)
  
  Id   <- diag(n)                               # identity matrix
  Q    <- qr.Q(qr(Xctr[ , 1, drop=FALSE]))      # QR-decomposition, just matrix Q
  P    <- tcrossprod(Q)          # = Q Q'       # projection onto space defined by x1
  x2o  <- (Id-P) %*% Xctr[ , 2]                 # x2ctr made orthogonal to x1ctr
  Xc2  <- cbind(Xctr[ , 1], x2o)                # bind to matrix
  Y    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  # scale columns to length 1
  
  x <- Y[ , 2] + (1 / tan(theta)) * Y[ , 1]     # final new vector
  cor(x1, x)                                    # check correlation = rho
  #df_data <- as.data.frame(cbind(x1,x))
  df_data <- as.data.frame(cbind(Xctr[,1],x*10))
  return(df_data)
}
# check correlation = rho

infile1_function <- file.path("/home/bparmentier/Google Drive/Papers_writing_MEOT/R_scripts/",
                              "EOT_PCA_rotation_functions_12272015.R")
functions_generation_space_time_datasets <- "generation_of_space_time_dataset_functions_12282015.R"
#script_path<-"/Users/benoitparmentier/Dropbox/Data/MEOT_paper/R_scripts_12312013" #path to script
script_path <- "/home/bparmentier/Google Drive/Papers_writing_MEOT/R_scripts/"
source(file.path(script_path,functions_generation_space_time_datasets)) #source all functions used in this script.
source(infile1_function)

############# Parameters and arguments ####################

inDir <- "/home/bparmentier/Google Drive/Papers_writing_MEOT/MEOT_analyses_02202015"

j <- 100 #number of images
col <- 10 # number of columns in the images 
row <- 10# number of rows in the images

proj_str <- NA #projection string (use PROJ4 convention or "NA")
range_val <- c(0,1)
distribution_val_i <-"none" #c("unif","norm","none") 
#distribution_val_i <-"unif" #c("unif","norm","none")
sd_val<-1
mean_val<-0
seed_nb<-100
file_format <-".rst"
NA_flag_val <- -9999
out_suffix <- "gen_12282015"

setwd(inDir)

#outDir <- "/Users/benoitparmentier/Dropbox/Data/Dissertation_paper2_04142012"
outDir <- inDir

create_outDir_param = TRUE

#Create output directory

if(create_outDir_param==TRUE){  
  outDir <- create_dir_fun(outDir,out_suffix)
  setwd(outDir)
}else{
  setwd(outDir) #use previoulsy defined directory
}
out_dir <- outDir
#########################################################
#################### BEGIN SCRIPT ########################

##### PART I : Prepare input temporal and spatial patterns #####

#step 1: generate reference raster image/region

list_param_gen_raster<-list(j,col,row,proj_str,range_val,NA_flag_val,distribution_val_i,out_suffix, out_dir,seed_nb)
names(list_param_gen_raster)<-c("j","col","row","proj_str","range_val","NA_flag_val","distribution_val_i","out_suffix","out_dir","seed_nb")

#debug(generate_raster_region)
r <- generate_raster_region(1,list_param_gen_raster)
plot(r)

#step 2: generate spatial structure
#step 2: generate reference spatial patterns

#j <- 100 #number of images
#col <- 10 # number of columns in the images 
#row <- 10 # number of rows in the images

proj_str <- NA #projection string (use PROJ4 convention or "NA")
range_val <- c(0,1)
distribution_val_i <-"none" #c("unif","norm","none") 
#distribution_val_i <-"unif" #c("unif","norm","none") 
#file_format <-".rst"
#NA_flag_val <- -9999
#out_suffix <- "gen_11252013"

#note that r is the raste image!!
list_param_adding_spatial_structure <-list(j,r,proj_str,range_val,NA_flag_val,file_format,out_suffix, out_dir)
names(list_param_adding_spatial_structure)<-c("j","r","proj_str","range_val","NA_flag_val","file_format","out_suffix","out_dir")

#debug(adding_sptatial_structure)
r_spat_s <-adding_sptatial_structure(1,list_param_adding_spatial_structure)
plot(r_spat_s) #examine available pattern
r_trend_x <- (subset(r_spat_s,"trend_x")) #examine available pattern

#Generate moving structure? Not in this case...
#no moving structure
#Step 3: generate temporal structure

#there are currently 5 datasets generated...

list_temp_pattern_produced <- vector("list",length=5)

### Generate orhogonal solution...


list_seed_nb <- c(100,101)
nt<- 100
rho_val <-0
#debug()
test2<- lapply(list_seed_nb,FUN=generate_orthogonal_and_correlation_vectors_fun,n=nt,rho=rho_val)

test<- lapply(1:length(list_seed_nb),FUN=generate_orthogonal_and_correlation_vectors_fun,n=nt,rho=rho_val)

plot(test2[[1]][,1],type="l")
lines(test2[[1]][,2],type="l",col="red")
plot(test2[[1]][,2],type="l",col="red")
cor(test2[[1]][,2],test2[[1]][,1])
cor(test2[[1]][,2],test2[[2]][,2])
plot(test2[[1]][,2],test2[[1]][,1])


list_seed_nb <- 100:109
nt<- 100
rho_val <-0

test<- lapply(list_seed_nb,FUN=generate_orthogonal_and_correlation_vectors_fun,n=nt,rho=rho_val)
list_df<- (lapply(1:length(test),FUN=function(i,x){as.data.frame(x[[i]][,2])},x=test))
df <- (do.call(cbind,list_df)) 
names(df) <- paste("v",1:10,sep="")
df$ref <- test[[1]][,1]
cor(df$v1,df$v2)
cor(df$v1,df$v3)
cor(df$v1,df$v4)
cor(df$v2,df$v4)

cor(df)
### Experimental setup: generate a set of 100 orthogonal vectors with and without random noise and run EOT and PCA for them.

#this represents hundred time steps, make it 10x10 pixels, so you need 100 time series!!
#Generate the maps then do PCA, PCA varimax and EOT and check how well it correlates back with the original variables

#Now that we have 10 variables orthogonal in time...produce map using the function provided!!
#Nothing generated in this particular dataset (see other versions)

##For each time series generate a different dataset with merging with the spatial structure
##When the 10 raster/image time series have been changed then combine all of them together by addition and then multiplication?
## Try EOT and PCA as a third step

#######################################################
######## PART II : Generate space-time datasets using input temporal and spatial patterns =
#NOW create a few datasets for the MEOT and MSSA/PCA diagnostics...


## Combine temporal and spatial structure here...


################# END OF SCRIPT ############################

#Orthogonal Vectors
#Two vectors u and v whose dot product is u·v=0 (i.e., the vectors are perpendicular) are said to be orthogonal. 
#In three-space, three vectors can be mutually perpendicular.

# k <- 5
# tstMat <- array(runif(k), dim=c(k,k))
# tstOrth <- qr.Q(qr(tstMat))
# t(tstOrth)%*%tstOrth
# 
# #http://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variable
# n     <- 20                    # length of vector
# rho   <- 0.6                   # desired correlation = cos(angle)
# theta <- acos(rho)             # corresponding angle
# x1    <- rnorm(n, 1, 1)        # fixed given data
# x2    <- rnorm(n, 2, 0.5)      # new random data
# X     <- cbind(x1, x2)         # matrix
# Xctr  <- scale(X, center=TRUE, scale=FALSE)   # centered columns (mean 0)
# 
# Id   <- diag(n)                               # identity matrix
# Q    <- qr.Q(qr(Xctr[ , 1, drop=FALSE]))      # QR-decomposition, just matrix Q
# P    <- tcrossprod(Q)          # = Q Q'       # projection onto space defined by x1
# x2o  <- (Id-P) %*% Xctr[ , 2]                 # x2ctr made orthogonal to x1ctr
# Xc2  <- cbind(Xctr[ , 1], x2o)                # bind to matrix
# Y    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  # scale columns to length 1
# 
# x <- Y[ , 2] + (1 / tan(theta)) * Y[ , 1]     # final new vector
# cor(x1, x)                                    # check correlation = rho
# 
# 
# n     <- 20                    # length of vector
# rho   <- 0.6                   # desired correlation = cos(angle)
# theta <- acos(rho)             # corresponding angle
# x1    <- rnorm(n, 1, 1)        # fixed given data
# x2    <- rnorm(n, 2, 0.5)      # new random data
# X     <- cbind(x1, x2)         # matrix
# Xctr  <- scale(X, center=TRUE, scale=FALSE)   # centered columns (mean 0)
# 
# Id   <- diag(n)                               # identity matrix
# Q    <- qr.Q(qr(Xctr[ , 1, drop=FALSE]))      # QR-decomposition, just matrix Q
# P    <- tcrossprod(Q)          # = Q Q'       # projection onto space defined by x1
# x2o  <- (Id-P) %*% Xctr[ , 2]                 # x2ctr made orthogonal to x1ctr
# Xc2  <- cbind(Xctr[ , 1], x2o)                # bind to matrix
# Y    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  # scale columns to length 1

