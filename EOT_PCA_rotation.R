#################################    EOT AND PCA ROTATION  #######################################
########################### SPACE-TIME VARIABILITY  #############################################
#This script analyzes  SST data.
#The goal is to compare EOT and PCA with rotation.
#
#AUTHOR: Benoit Parmentier                                                                       #
#DATE CREATED: 07/02/2015 
#DATE MODIFIED: 08/15/2015
#
#PROJECT: Land transitions from Remote Sensing
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

#inDir <- "J:/Benoit/Data/MEOT_analyses_02202015"
inDir <- "/home/parmentier/Data/MEOT12272012/Papers_writing_MEOT/MEOT_analyses_02202015"

SST_dir <- "SST_1982_2007"
mask_fname <- "mask_rgf_1_1.tif"
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
npc <- 10 # number of pca to produce...put this at the beginning of the script...
lag_window <-13
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

### STEP 0: READ IN DATASETS RELATED TO TELECONNECTION AND PREVIOUS LOADINGS

#SAODI<-read.table(infile2,sep=",", header=TRUE)
#Prepare data to write out in a textfile
#s_SAODI2 <- subset(SAODI, year>1981 & year<2008) #Subset rows that correspond to the conditions
#s_SAODI2<- subset(s_SAODI2, select=-year) #Remove the column year
#in_SAODI2<-as.vector(t(s_SAODI2))   #This transform the data frame into a one colum
#write.table(in_SAODI2,file=paste("SAOD_index_1981_2007",out_prefix,".txt",sep=""),sep=",")
 
r_sst <- stack(lf_sst)
projection(r_sst) <- CRS_WGS84 #assign projection
r_mask <- raster(file.path(inDir,mask_fname))
projection(r_mask) <- CRS_WGS84
r_mask_NA <- r_mask
NAvalue(r_mask_NA) <- 0
#NAvalues()
levelplot(r_sst,layer=1)
levelplot(r_mask,layer=1)

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

### DO PCA WITH SST_df: T-mode cross product from a standardized dataset...(i.e. correlation matrix)

A<-SST_df[1:(ncol(SST_df)-2)] #drop s1 and s2 which contain coordinates
dim(A) #n*c or 42098*312
At<-scale(A) #center and reduce using sd and mean
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
Et<-eigen(cAt)

##Cross product PCA T mode
#pca_SST<-principal(r=cAt, nfactors = npc, residuals = FALSE, covar=TRUE,rotate = "none")
pca_SST_t_mode<-principal(r=cAt, nfactors = npc, residuals = FALSE, covar=TRUE,rotate = "none") #psych package
unclass(pca_SST_t_mode$loadings) # extract the matrix of ??
head(unclass(pca_SST_t_mode$loadings)) # extract the matrix of ??
#<- as.data.frame(unclass(pca_SST_t_mode$loadings))) # extract the matrix of ??

pca_SST_t_mode_varimax <-principal(r, nfactors = npc, residuals = FALSE,rotate="varimax",n.obs=NA, covar=FALSE, scores=FALSE,missing=FALSE,
          impute="median",oblique.scores=TRUE,method="regression",...)

(nrow(At)-1)*ncol(At)  #Sum of the diagonal is 42097*312
sum(Et$value)   #equal to 312
sum(diag(cAt))  #
sum(pca_SST_t_mode$value)
mean(Et$vectors)  #This is equal to zero since there was centering!!
1/var(Et$vectors)[1]
var(Et$vectors)[1]
var(Et$vectors)[2]
var(Et$vectors)[3]

## read in IDRISI and SPSS results
pca_IDRISI_covar<-read.xls(file.path(inDir,infile_pca_IDRISI), sheet=1)
pca_IDRISI_cor<-read.xls(file.path(inDir,infile_pca_IDRISI), sheet=2)
pca_IDRISI_eigenvalues<-read.xls(file.path(inDir,infile_pca_IDRISI), sheet=3)
pca_IDRISI_eigenvectors<-read.xls(file.path(inDir,infile_pca_IDRISI), sheet=4)
pca_IDRISI_loadings<-read.xls(file.path(inDir,infile_pca_IDRISI), sheet=5)

spss_data_loadings_pca_unrotated <-read.xls(file.path(inDir,infile_pca_spss), sheet="loadings_pca_unrotated") #this is T-mode using cor matrix
spss_data_loadings_pca_varimax <-read.xls(file.path(inDir,infile_pca_spss), sheet="loadings_pca_rotated") #this is T-mode using cor matrix, varimax rotation

 
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


plot(spss_data_loadings_pca_unrotated$pc1,unclass(pca_SST_t_mode$loadings)[,1])
cor(spss_data_loadings_pca_unrotated$pc1,unclass(pca_SST_t_mode$loadings)[,1])

plot(pca_IDRISI_loadings$CMP.1,unclass(pca_SST_t_mode$loadings)[,1])
cor(pca_IDRISI_loadings$CMP.1,unclass(pca_SST_t_mode$loadings)[,1])


