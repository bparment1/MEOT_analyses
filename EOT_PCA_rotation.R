#################################    EOT AND PCA ROTATION  #######################################
########################### SPACE-TIME VARIABILITY  #############################################
#This script analyzes  SST data.
#The goal is to compare EOT and PCA with rotation and without
#This is part of the Time Series analysis project started at IDRISI.
#
#AUTHOR: Benoit Parmentier                                                                       #
#DATE CREATED: 07/02/2015 
#DATE MODIFIED: 09/05/2015
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

cv_test_fun <- function(x,y) {
  #modified from:
  #http://www.r-bloggers.com/example-8-39-calculating-cramers-v/
  chisq_val <- chisq.test(x, y, correct=FALSE)
  CV = sqrt(chisq_val$statistic /
              (length(x) * (min(length(unique(x)),length(unique(y))) - 1)))
  print.noquote("Cram?r V / Phi:")
  CV <- as.numeric(CV)
  cv_obj<- list(CV=CV,chisq_val=chisq_val)
  return(cv_obj)
}


pca_to_raster_fun<-function(pc_spdf,ref_raster,NA_val,out_prefix){
  #Input arguments:
  #pc_spdf: spdf with scores components and
  #         must include x,y in the last 2 columns
  #ref_raster: reference raster giving the output extent and NA if present
  #NA_val<- values assigned to NA pixels, if "NA", use default
  
  npc<-ncol(pc_spdf)-2 #number of components 
  pc_scores_lf<-vector("list",npc) 
  for (k in 1:npc){
    pc_scores<-pc_spdf[,k] #vector containing scores for components k
    pc_name<-names(pc_spdf)[k] #name of the current component
    raster_name<-paste("pc_component_",k,"_",out_prefix,".rst",sep="")
    
    if (NA_val=="NA"){
      pc_lag<-rasterize(pc_scores,ref_raster,pc_name,fun=min,overwrite=TRUE,
                        filename=raster_name)
    }
    #Use provided NA value for the background
    if (NA_val!="NA"){
      pc_lag<-rasterize(pc_scores,ref_raster,pc_name,background=NA_val,fun=min,overwrite=TRUE,
                        filename=raster_name)
    }
    pc_scores_lf[[k]]<-raster_name
  } 
  return(pc_scores_lf)
}
#############################################
######## Parameters and arguments  ########

#inDir <- "J:/Benoit/Data/MEOT_analyses_02202015" #IDRISI 690-2 computer
inDir <- "/home/parmentier/Data/MEOT12272012/Papers_writing_MEOT/MEOT_analyses_02202015" #Atlas
inDir <- "/home/bparmentier/Google Drive/Papers_writing_MEOT/MEOT_analyses_02202015" #bpy50
SST_dir <- "SST_1982_2007"
eot_dir <- "/home/bparmentier/Google Drive/Papers_writing_MEOT/workdir_terrset_08282015/anom_sst_1982_2007/components"
mask_fname <- "mask_rgf_1_1.tif"
eot_fname1 <- "eot_std_s7_test__EOT_Center_Std.avl"
out_suffix <-"_eot_pca_07022015"
CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84 # CONST 2

lf_sst <- list.files(path=file.path(inDir,SST_dir),pattern=".rst$",full.names=T)


#proj_str<-"+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

infile1<-"SAODI-01-1854_06-2011_test.asc"             #GHCN shapefile containing variables for modeling 2010                 
infile2<-"SAODI-01-1854_06-2011.csv"                     #List of 10 dates for the regression
infile3<-"MEOT_MSSA_Telcon_indices_12112012.xlsx"                        #LST dates name
#infile4<-"mask_rgf_1_1.rst"
infile_pca_IDRISI <-"pca_IDRISI_03302013_T-Mode_DIF.xlsx"

infile_pca_spss <- "spss_pca_output_08082015.xlsx"

infile_s_mode_pca  <- "pca_IDRISI_04152013_S-Mode_DIF.xlsx"
npc <- 10 # number of pca components to produce...put this at the beginning of the script...
lag_window <-13 #not in use here...
telind<-c("PNA","NAO","TNA","TSA","SAOD","MEI","PDO","AO","AAO","AMM","AMOsm","QBO")
mode_list_MEOT<-c("MEOT1","MEOT3", "MEOT4","MEOT7","MEOT10","MEOT15","MEOT16")
mode_list_PCA<-c("MSSA1","MSSA2","MSSA3","MSSA4","MSSA5","MSSA6")    
#mode_list_PCA<-paste("MSSA",1:15,sep="")
#out_prefix<-"pca_test_04152013"

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

########################################################
##############  Start of th script  ##############

### PART 0: READ IN DATASETS RELATED TO TELECONNECTION AND PREVIOUS LOADINGS

SAODI<-read.table(file.path(inDir,infile2),sep=",", header=TRUE)
#Prepare data to write out in a textfile
s_SAODI2 <- subset(SAODI, year>1981 & year<2008) #Subset rows that correspond to the conditions
s_SAODI2<- subset(s_SAODI2, select=-year) #Remove the column year
in_SAODI2<-as.vector(t(s_SAODI2))   #This transform the data frame into a one colum
write.table(in_SAODI2,file=paste("SAOD_index_1981_2007",out_suffix,".txt",sep=""),sep=",")
 
## read in eot info:

#eot_avl_files <- mixedsort(list.files(eot_dir,"*.avl",full.names=TRUE))
#eot_avl <- read.table(eot_avl_files[11])
#df_eot_f <- lapply(eot_avl_files,FUN=function(x){read.table(x)})
eot_s7_df <- read.table(file.path(eot_dir, eot_fname1))
n_eot <- (ncol(eot_s7_df)) -1
names(eot_s7_df) <- c("time",paste("eot_",1:n_eot,sep=""))

#later on read eot_s1_df: eot using agg 1 ie no aggregation for easier comparison to EOF
##
r_sst <- stack(lf_sst)
projection(r_sst) <- CRS_WGS84 #assign projection
r_mask <- raster(file.path(inDir,mask_fname))
projection(r_mask) <- CRS_WGS84
r_mask_NA <- r_mask
NAvalue(r_mask_NA) <- 0
#NAvalues()
levelplot(r_sst,layer=1)
levelplot(r_mask,layer=1)
plot(r_sst,y=1,col=matlab.like(255))

r_sst_m <- mask(r_sst,r_mask_NA,filename=paste("sst_",out_suffix,".tif",sep=""),overwrite=T)
levelplot(r_sst_m,layer=1)

#Multiply layer by weight

lat_coord<-coordinates(r_mask_NA)[,2]
w_rast<- r_mask_NA
values(w_rast)<-cos(lat_coord*pi/180) #use area in raster package for weight of a cell??
w_rast_m<-mask(w_rast,r_mask_NA)
levelplot(w_rast_m)

r_sst_m <- r_sst_m*w_rast #Maybe use overlay to avoid putting in memory?
writeRaster(r_sst_m,filename="ANOM_SST_1982_2007_weighted.tif",overwrite=TRUE)
rm(r_sst)
r_sst <-brick("ANOM_SST_1982_2007_weighted.tif")

SST1_m <- subset(r_sst,1)
plot(SST1_m,col=matlab.like(255))  

##############################################################
#### STEP 2: get data ready for PCA by lagging data frame...

SST_sgdf<-as(r_sst,"SpatialGridDataFrame")
#SST_df<-as.data.frame(SST_rast) #drop the coordinates x and y
SST_df<-as.data.frame(SST_sgdf)
write.table(SST_df,file.path(outDir,"sst_f_weighted.txt"),sep=",",row.names=F)
rm(SST_sgdf)

#nv<-312
SST_xy<-SST_df[,c("s1","s2")]


### DO PCA WITH SST_df: T-mode cross product from a standardized dataset...(i.e. correlation matrix)

dim(SST_df)

############# PART 1: CARRY PCA IN T AND S mode with/without varimax rotation ###################

### DO PCA WITH SST_df: T-mode cross product from a standardized dataset...(i.e. correlation matrix)

A <- SST_df[1:(ncol(SST_df)-2)] #drop s1 and s2 which contain coordinates
dim(A) #n*c or 42098*312
At <- scale(A) #center and reduce using sd and mean, T mode standardized
sd(At[,2]) #this must be equal to 1 since At is standardized by columns
mean(At[,2]) #this must be equal to 0 since At is standardized by columns
names_var<-paste("t",1:ncol(At),sep="_") #Rename columns 
colnames(At)<-names_var
#perform PCA with matrix operation...
cAt<-t(At)%*%At #T mode!! --> (c*n)by(n*c) should get c*c matrix or 312*312 
diag(cAt) #equal 1 since it is a correlation, only if divided by n-1!!!
cAt<-cAt/(nrow(At)-1)
diag(cAt)
#cAt<-cAt/nrow(At) #reduce by number of row...!!! THIS IS NOT NECESSARY SINCE At has been already based on standardized??
#note that cAt is standardized and has been divided by n, so it is equivalent to a correlation matrix
Et<-eigen(cAt) #this creates a list of two objects
Et$vectors #312*312, eigenvectors...
Et$values #1*312, eigenvalues...

PC_scores <- At%*%Et$vectors # to obtain the scores, multiply the eigenvector matrix by the data matrix
PC_scores <- as.matrix(A)%*%Et$vectors # to

###COMPARING SCORES FROM FROM "eigen" and principal function
## Note that the scores are equal only if the square root of eigenvalues has been used to reduce (normalize) the Eigen Vectors!!!


##Cross product PCA T mode
#pca_SST<-principal(r=cAt, nfactors = npc, residuals = FALSE, covar=TRUE,rotate = "none",scores=T)
pca_SST_t_mode_unrotated <-principal(r=cAt, nfactors = npc, residuals = FALSE, 
                                     covar=TRUE,rotate = "none",
                                     scores=TRUE,
                                     missing=FALSE,
                                     oblique.scores=TRUE,
                                     method="regression") #psych package
class(unclass(pca_SST_t_mode_unrotated$loadings)) # extract the matrix of ??
head(unclass(pca_SST_t_mode_unrotated$loadings)) # extract the matrix of ??

pca_SST_t_mode_varimax <-principal(r=cAt, nfactors = npc, residuals = FALSE,rotate="varimax",
                                   n.obs=NA, covar=FALSE, scores=FALSE,missing=FALSE,
                                   impute="median",oblique.scores=TRUE,method="regression")
principal_t_mode_df_unrotated <- as.data.frame(unclass(pca_SST_t_mode_unrotated$loadings)) # extract the matrix of ??
principal_t_mode_df_varimax <- as.data.frame(unclass(pca_SST_t_mode_varimax$loadings)) # extract the matrix of ??

(nrow(At)-1)*ncol(At)  #Sum of the diagonal is 42097*312
sum(Et$value)   #equal to 312
sum(diag(cAt))  #
sum(pca_SST_t_mode$value)
mean(Et$vectors)  #This is equal to zero since there was centering!!
1/var(Et$vectors)[1]
var(Et$vectors)[1]
var(Et$vectors)[2]
var(Et$vectors)[3]

pc_scores <- At%*%Et$vectors # to
pc_scores <- as.matrix(A)%*%Et$vectors # to

#To obtain de scores from the principal object, use predict with the input matrix?
class(pca_SST_t_mode_unrotated)
pc_t_principal_scores1 <- predict(pca_SST_t_mode_unrotated,A)
pc_t_principal_scores2 <- predict(pca_SST_t_mode_unrotated,At)
pc_t_principal_scores1[,1]-pc_t_principal_scores2[,1]

cor(pc_t_principal_scores[,1],pc_scores[,1])
plot(pc_t_principal_scores[,1],pc_scores[,1])

plot(principal_t_mode_df_unrotated[,1],Et$vectors[,1])
cor(principal_t_mode_df_unrotated[,1],Et$vectors[,1])# ok correlation of 1!!

diff <- principal_t_mode_df_unrotated[,1] - Et$vectors[,1]
sum(diff)
(principal_t_mode_df_unrotated[,1])
mean(Et$vectors[,1])


cor(PC_scores[,1],pca_scores[,1]) 
plot(PC_scores[,1],pca_scores[,1]) 
PC_scores[1:10,1]/pca_scores[1:10,1]
sqrt(Et$values[1])
(PC_scores[1:10,1]/sqrt(Et$values[1]))/pca_scores[1:10,1]
plot(PC_scores[,1]/sqrt(Et$values[1]),pca_scores[,1]) 

mean(PC_scores[,1])
sd(PC_scores[,1])
#The regression weights are found from the inverse of the correlation matrix times the component loadings. 
#This has the result that the component scores are standard scores (mean=0, sd = 1) of the standardized 
#input.
mean(pc_t_principal_scores[,1]) # mean of 0
sd(pc_t_principal_scores[,1])  # sd of 1
#This has been standardize!!! Principal standardizes the values!!

PC_scores <- At%*%Et$vectors # to
PC_scores <- as.matrix(A)%*%Et$vectors # to
mean(PC_scores[,1])
sd(PC_scores[,1])


cor(PC_scores[,1],pca_scores[,1]) 
plot(PC_scores[,1],pca_scores[,1]) 
PC_scores[1:10,1]/pca_scores[1:10,1]
sqrt(Et$values[1])
(PC_scores[1:10,1]/sqrt(Et$values[1]))/pca_scores[1:10,1]
plot(PC_scores[,1]/sqrt(Et$values[1]),pca_scores[,1]) 

mean(pca_scores[,1]) #This has been standardize!!! Principal standardizes the values!!
sd(pca_scores[,1])

####################################################################
### NOW DO PCA WITH S-mode

## DO PCA WITH SST_df: S-mode cross product but standardized pixel-wise
#PC_scores[1:10,1]/pca_scores[1:10,1]
A <- SST_df[1:(ncol(SST_df)-2)] #remove x,y fields to obtain a data.frame with 42098*312 dimension
Atr<-t(A)
dim(Atr)# 312 *42098 rows*columns
Ats<-scale(Atr) #center time series...(observations/time profiles or pixels in images)
As<-t(Ats) #transpose back
cAs<-Ats%*%As #get cross product matrix wihtout centering or scaling in the var dimensions
cAs<-cAs/nrow(cAs)
Es<-eigen(cAs) #cross product S-mode
diag(cAs)

##Cross product PCA S mode
#pca_SST<-principal(r=cAs, nfactors = npc, residuals = FALSE, covar=TRUE,rotate = "none")
pca_SST_s_mode_unrotated <-principal(r=cAs, nfactors = npc, residuals = FALSE, covar=TRUE,rotate = "none")
unclass(pca_SST_s_mode_unrotated$loadings) # extract the matrix of ??

pca_SST_s_mode_varimax <-principal(r=cAs, nfactors = npc, residuals = FALSE,rotate="varimax",
                                   n.obs=NA, covar=FALSE, scores=FALSE,missing=FALSE,
                                   impute="median",oblique.scores=TRUE,method="regression")

principal_s_mode_df_unrotated<- as.data.frame(unclass(pca_SST_s_mode_unrotated$loadings)) # extract the matrix of ??
principal_s_mode_df_varimax<- as.data.frame(unclass(pca_SST_s_mode_varimax$loadings)) # extract the matrix of ??

sum(Es$value)
sum(diag(cAs)) #sum of eigenvalue is equal to the number of variable...
sum(pca_SST_s_mode_unrotated$value)

##Cross product PCA T mode
#pca_SST<-principal(r=cAs, nfactors = npc, residuals = FALSE, covar=TRUE,rotate = "none")
#pca_SST_s_mode<-principal(r=cAs, nfactors = npc, residuals = FALSE, covar=TRUE,rotate = "none")
#unclass(pca_SST_s_mode$loadings) # extract the matrix of ??

sum(Es$value)
sum(diag(cAs)) #sum of eigenvalue is equal to the number of variable...
#sum(pca_SST_s_mode$value)

#### PART 2: Assessing and comparing all PCA methods

###### CREATE SCORES IMAGES USING PREDICTED SCORES FROM PRINCIPAL
pca_scores <- predict(pca_SST_t_mode,A)  #generate scores from original matrix and object
pca_scores_spdf<-cbind(as.data.frame(pca_scores),SST_xy) #add coordinates
coordinates(pca_scores_spdf)<-pca_scores_spdf[,c("s1","s2")] #promote to spdf

#Now generate raster images...
pca_scores_lf<-pca_to_raster_fun(pca_scores_spdf,SST1_m,-9999,out_suffix)

pca_SST_s<-stack(pca_scores_lf)
NAvalue(pca_SST_s)<- -9999
#plot(subset(pca_IDRISI_s,1))


### Comparison Loadings and scores between R, SPSS and IDRISI

## read in IDRISI and SPSS results
pca_IDRISI_covar<-read.xls(file.path(inDir,infile_pca_IDRISI), sheet=1)
pca_IDRISI_cor<-read.xls(file.path(inDir,infile_pca_IDRISI), sheet=2)
pca_IDRISI_eigenvalues<-read.xls(file.path(inDir,infile_pca_IDRISI), sheet=3)
pca_IDRISI_eigenvectors<-read.xls(file.path(inDir,infile_pca_IDRISI), sheet=4)
pca_IDRISI_loadings<-read.xls(file.path(inDir,infile_pca_IDRISI), sheet=5)

spss_data_loadings_pca_unrotated <-read.xls(file.path(inDir,infile_pca_spss), sheet="loadings_pca_unrotated") #this is T-mode using cor matrix
spss_data_loadings_pca_varimax <-read.xls(file.path(inDir,infile_pca_spss), sheet="loadings_pca_rotated") #this is T-mode using cor matrix, varimax rotation

plot(spss_data_loadings_pca_varimax[["pc10"]],type="l")
lines(spss_data_loadings_pca_unrotated[["pc10"]],col="red")
cor(spss_data_loadings_pca_varimax[["pc10"]],spss_data_loadings_pca_unrotated[["pc10"]])
cor(spss_data_loadings_pca_varimax[["pc5"]],spss_data_loadings_pca_unrotated[["pc5"]])

plot(spss_data_loadings_pca_unrotated$pc1,unclass(pca_SST_t_mode$loadings)[,1])
cor(spss_data_loadings_pca_unrotated$pc1,unclass(pca_SST_t_mode$loadings)[,1])

plot(pca_IDRISI_loadings$CMP.1,unclass(pca_SST_t_mode$loadings)[,1])
cor(pca_IDRISI_loadings$CMP.1,unclass(pca_SST_t_mode$loadings)[,1])

pca_SST_t_mode_varimax 

head(principal_df_varimax)
plot(principal_df_varimax[["PC1"]],type="l")
lines(principal_df_unrotated[["PC1"]],col="red")
lines(spss_data_loadings_pca_varimax[["pc1"]],col="blue")

plot(principal_s_mode_df_unrotated[["PC1"]],col="red",type="l")
lines(principal_s_mode_df_varimax[["PC1"]],type="l",col="darkgreen")

### Comparison between SPSS and R, T mode unrotated!!
plot(principal_t_mode_df_unrotated[["PC1"]],col="red",type="l")
lines(principal_t_mode_df_varimax[["PC1"]],type="l",col="darkgreen")
lines(spss_data_loadings_pca_unrotated[["pc1"]],col="blue")

plot(principal_t_mode_df_unrotated[["PC1"]],spss_data_loadings_pca_unrotated[["pc1"]])
cor(principal_t_mode_df_unrotated[["PC1"]],spss_data_loadings_pca_unrotated[["pc1"]])
cor(principal_t_mode_df_unrotated[["PC2"]],spss_data_loadings_pca_unrotated[["pc2"]])
cor(principal_t_mode_df_unrotated[["PC3"]],spss_data_loadings_pca_unrotated[["pc3"]])

plot(principal_t_mode_df_varimax[["PC1"]],spss_data_loadings_pca_varimax[["pc1"]])
cor(principal_t_mode_df_varimax[["PC1"]],spss_data_loadings_pca_varimax[["pc1"]])
cor(principal_t_mode_df_varimax[["PC2"]],spss_data_loadings_pca_varimax[["pc2"]])
cor(principal_t_mode_df_varimax[["PC3"]],spss_data_loadings_pca_varimax[["pc3"]])
cor(principal_t_mode_df_varimax[["PC4"]],spss_data_loadings_pca_varimax[["pc4"]])

####### Comparison between EOT and rotation...
#ok now we should do this for different truncation...if one can show that the corr
#increases when the truncation decreases (more pc) then we will be ok...

run_pca_fun(mode="T",rotation_opt="none",df_val,matrix_val=NULL,npc=1,loadings=TRUE,scores_opt=TRUE){
  #Arguments
  #mode: T, s mode or other mode , note that this option is not in use at this moment, use matrix_val
  #rotation_opt: which option to use to rotate pc
  #matrix_val
  #if(!is.null(matrix_val)){
  #  as.matrix(df_val)
  #}

  ##Cross product PCA S mode
  #pca_SST<-principal(r=cAs, nfactors = npc, residuals = FALSE, covar=TRUE,rotate = "none")
  pca_principal_obj <-principal(r=matrix_val, nfactors = npc, residuals = FALSE, 
                                       covar=TRUE, #use covar option...
                                       rotate = rotation_opt,
                                       scores=scores_opt,
                                       oblique.scores=T,
                                       method="regression")
  
  unclass(pca_SST_s_mode_unrotated$loadings) # extract the matrix of ??
  
  pca_SST_s_mode_varimax <-principal(r=cAs, nfactors = npc, residuals = FALSE,rotate="varimax",
                                     n.obs=NA, covar=FALSE, scores=FALSE,missing=FALSE,
                                     impute="median",oblique.scores=TRUE,method="regression")
  
  principal_loadings_df <- as.data.frame(unclass(pca_principal_obj$loadings)) # extract the matrix of ??
  #Add scores...
  #
}

#diff could also be due to spatial average of 7 in EOT
#make a function for comparison...

#Also compare to PCA...

cor(eot_s7_df$eot_1,principal_s_mode_df_unrotated[["PC1"]])
cor(eot_s7_df$eot_1,principal_s_mode_df_varimax[["PC1"]])

cor(eot_s7_df$eot_2,principal_s_mode_df_unrotated[["PC2"]])
cor(eot_s7_df$eot_2,principal_s_mode_df_varimax[["PC2"]])

cor(eot_s7_df$eot_3,principal_s_mode_df_unrotated[["PC3"]])
cor(eot_s7_df$eot_3,principal_s_mode_df_varimax[["PC3"]])

cor(eot_s7_df$eot_4,principal_s_mode_df_unrotated[["PC4"]])
cor(eot_s7_df$eot_4,principal_s_mode_df_varimax[["PC4"]])

cor(eot_s7_df$eot_5,principal_s_mode_df_varimax[["PC5"]])
cor(eot_s7_df$eot_5,principal_s_mode_df_varimax[["PC5"]])

#function perform correlation as a function

lines(spss_data_loadings_pca_varimax[["pc1"]],col="blue")
plot(principal_s_mode_df_unrotated[["PC1"]],col="red",type="l")

plot(spss_data_loadings_pca_varimax[["pc10"]],type="l")
lines(spss_data_loadings_pca_unrotated[["pc10"]],col="red")

#################################################################
########################## END OF SCRIPT #########################
#################################################################

#cAt<-A%*%At
#Ae<-eigen(cAt)
#principal(r, nfactors = 1, residuals = FALSE,rotate="varimax",n.obs=NA, covar=FALSE, scores=FALSE,missing=FALSE,
#          impute="median",oblique.scores=TRUE,method="regression",...)



##Correlation PCA T mode
#pca_SST<-principal(r=cAt, nfactors = npc, residuals = FALSE, covar=TRUE,rotate = "none")

#pca_SST_cor_t_mode<-principal(r=corAt, nfactors = npc, residuals = FALSE, covar=TRUE,rotate = "none")
#unclass(pca_SST_t_mode$loadings) # extract the matrix of ??
# 
# sum(Et$value)
# sum(diag(cAt))
# sum(pca_SST_t_mode$value)
# 
# head(unclass(pca_SST_t_mode$loadings)[,1:npc])
# head(Et$vectors[,1:npc]) #not equal...because not standardized...
# Estd <-  Et$vectors%*% diag(sqrt(Et$values))
# head(Estd)
# head(unclass(pca_SST$loadings)[,1:npc])
# 
# #cAt<-A%*%At
# #Ae<-eigen(cAt)
# #principal(r, nfactors = 1, residuals = FALSE,rotate="varimax",n.obs=NA, covar=FALSE, scores=FALSE,missing=FALSE,
# #          impute="median",oblique.scores=TRUE,method="regression",...)
# 
# 
# ##############
# ### 
# require(psych) 
# data(mtcars) 
# rawd <- mtcars[,c("am","carb","cyl","disp","drat","gear","hp","mpg")] #subset the original data.frame
# 
# ## compute acp 
# .PC <- princomp(~am+carb+cyl+disp+drat+gear+hp+mpg, cor=TRUE, data=mtcars) 
# pca <- principal(rawd, nfactors = 8, residuals = T, rotate="none", scores=T) 
# 
# ## eigenvectors of the correlation matrix of raw data 
# eigens <- eigen(cor(rawd)) #eigen value decomposition
# eigens$vectors #matrix of the eigenvectors, eigens is a list
# unclass(loadings(.PC))  # component 'loadings' in princomp parlance are the eigenvectors!!!
# 
# ## correlation matrix between raw data and unrotated component scores 
# ## 'loadings' in ?principal parlance and 'component matrix' in SPSS --> loadings in princomp are the eigenvectors 
# # components scores can be obtained by multiplying the original data by the eigenvectors...
# eigens$vectors %*% diag(sqrt(eigens$values)) ##standardizing the eigenvectors...(it is already centered)
# cor(cbind(rawd, .PC$scores))  #loadings are the correlation between the PC scores and the origina variables
# unclass(pca$loadings) # extract the matrix of ??
# 
# ## extract un-rotated scores of principal components 
# head(scale(rawd) %*% eigens$vectors) # app, but very similar results 
# head(.PC$scores) 
# head(pca$scores) # scale'd, and obtained via regression on scale'd data 
# head(factor.scores(scale(rawd), 
#                    unclass(pca$loadings))) # same scores as from ?principal 
# #for differeneces between ?princomp and ?principal scores 
# #see last paragraph of Details in ?principal 
# rm(PC_std_scores)

