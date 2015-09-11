#################################    EOT AND PCA ROTATION  #######################################
########################### SPACE-TIME VARIABILITY  #############################################
#This script analyzes  SST data.
#The goal is to compare EOT and PCA with rotation and without
#This is part of the Time Series analysis project started at IDRISI.
#
#AUTHOR: Benoit Parmentier                                                                       #
#DATE CREATED: 07/02/2015 
#DATE MODIFIED: 09/15/2015
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

infile1_function <- file.path("/home/bparmentier/Google Drive/Papers_writing_MEOT","EOT_PCA_rotation_functions_09102015.R")
source(infile1_function)

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
idrisi_dir <- "~/Google Drive/Papers_writing_MEOT/MEOT_analyses_02202015/workdir_terrset_08282015/anom_sst_1982_2007/components"
infile_idrisi_pca_s_mode <- "idrisi_pca_s_mode_PCA_Center_Std_S-Mode_DIF.xls"
infile_idrisi_pca_t_mode <- "idrisi_pca_t_mode_PCA_Center_Std_T-Mode_DIF.xls"

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

SAODI <- read.table(file.path(inDir,infile2),sep=",", header=TRUE)
#Prepare data to write out in a textfile
s_SAODI2 <- subset(SAODI, year>1981 & year<2008) #Subset rows that correspond to the conditions
s_SAODI2 <- subset(s_SAODI2, select=-year) #Remove the column year
in_SAODI2 <- as.vector(t(s_SAODI2))   #This transform the data frame into a one colum
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
At <- scale(A) #center and reduce using sd and mean, T mode standardized,also makes it a matrix?
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

pc_scores <- At%*%Et$vectors # to obtain the scores, multiply the eigenvector matrix by the data matrix
PC_scores <- as.matrix(A)%*%Et$vectors # to

sd(PC_scores[,1])


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
sum(pca_SST_t_mode_unrotated$values)
mean(Et$vectors)  #This is equal to zero since there was centering!!
1/var(Et$vectors)[1]
var(Et$vectors)[1]
var(Et$vectors)[2]
var(Et$vectors)[3]

mean(PC_scores[,1])
sd(PC_scores[,1])
#The regression weights are found from the inverse of the correlation matrix times the component loadings. 
#This has the result that the component scores are standard scores (mean=0, sd = 1) of the standardized 
#input.
mean(pc_t_principal_scores[,1]) # mean of 0
sd(pc_t_principal_scores[,1])  # sd of 1
#This has been standardize!!! Principal standardizes the values!!

pc_scores1 <- At%*%Et$vectors # to
pc_scores2 <- as.matrix(A)%*%Et$vectors # to

loadings_Et1 <- unlist(lapply(1:ncol(pc_scores1),FUN=function(i){cor(pc_scores1[,1],At[,i])}))
loadings_Et2 <- unlist(lapply(1:ncol(pc_scores1),FUN=function(i){cor(pc_scores1[,2],At[,i])}))


plot(pc_scores1-pc_scores2)
sd(pc_scores1[,1])
sd(pc_scores1[,2])
sd(A[,1])
sd(pc_scores2[,1]) #this is not standardized!!
sd(pc_scores2[,2])


#To obtain de scores from the principal object, use predict with the input matrix?
class(pca_SST_t_mode_unrotated)
pc_t_principal_scores1 <- predict(pca_SST_t_mode_unrotated,A)
pc_t_principal_scores2 <- predict(pca_SST_t_mode_unrotated,At)
plot(pc_t_principal_scores1[,1]-pc_t_principal_scores2[,1])

cor(pc_t_principal_scores1[,1],pc_scores[,1])
plot(pc_t_principal_scores1[,1],pc_scores[,1])

plot(principal_t_mode_df_unrotated[,1],Et$vectors[,1]) #ok eigenvectors are the same but not the scores..
cor(principal_t_mode_df_unrotated[,1],Et$vectors[,1])# ok correlation of 1!!

diff <- principal_t_mode_df_unrotated[,1] - Et$vectors[,1]
sum(diff)
(principal_t_mode_df_unrotated[,1])
mean(Et$vectors[,1])

cor(pc_t_principal_scores1[,1],pc_scores[,1])
plot(pc_t_principal_scores1[,1],pc_scores[,1])

pc_scores[1:10,1]/pc_t_principal_scores1[1:10,1]
sqrt(Et$values[1])
sd(pc_scores[,1]/sqrt(Et$values[1]))
pc_t_principal_scores1[1:10,1]
plot(pc_scores[,1]/sqrt(Et$values[1]),pc_t_principal_scores1[,1]) 


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

Es_std <- Es$vectors %*% diag(sqrt(Es$values)) ##standardizing the eigenvectors...(it is already centered)

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

plot(Es$vectors[,1],principal_s_mode_df_unrotated[,1])
plot(Es$vectors[,2],principal_s_mode_df_unrotated[,2])
cor(Es$vectors[,1],principal_s_mode_df_unrotated[,1])
cor(Es$vectors[,2],principal_s_mode_df_unrotated[,2])

plot(Es$vectors[,1],type="l")
lines(Es$vectors[,2],type="l",col="blue")
cor(Es$vectors[,1],Es$vectors[,2])

#### PART 2: Assessing and comparing all PCA methods

###### CREATE SCORES IMAGES USING PREDICTED SCORES FROM PRINCIPAL
pca_s_mode_scores <- predict(pca_SST_s_mode_unrotated,At)  #generate scores from original matrix and object
pca_s_mode_scores_spdf<-cbind(as.data.frame(pca_s_mode_scores),SST_xy) #add coordinates
coordinates(pca_s_mode_scores_spdf) <- pca_s_mode_scores_spdf[,c("s1","s2")] #promote to spdf

#Now generate raster images...
pca_s_mode_scores_spdf_lf <- pca_to_raster_fun(pca_s_mode_scores_spdf,SST1_m,-9999, out_suffix) #add a mask image??
r_s_mode_pca <- stack(pca_s_mode_scores_spdf_lf)

pca_t_mode_scores <- predict(pca_SST_t_mode_unrotated,A)  #generate scores from original matrix and object
pca_t_mode_scores_spdf<-cbind(as.data.frame(pca_t_mode_scores),SST_xy) #add coordinates
coordinates(pca_t_mode_scores_spdf) <- pca_t_mode_scores_spdf[,c("s1","s2")] #promote to spdf

#Now generate raster images...
out_suffix_str <- paste("_t_mode_",out_suffix,sep="")
pca_t_mode_scores_scores_lf <- pca_to_raster_fun(pca_t_mode_scores_spdf,SST1_m,NA_val=-9999,out_suffix_str) #add a mask image??
r_t_mode_pca  <- stack(pca_t_mode_scores_scores_lf)

pca_SST_s <- stack(r_s_mode_pca)

NAvalue(r_s_mode_pca )<- -9999
NAvalue(r_t_mode_pca )<- -9999

plot(r_s_mode_pca,1:4)
plot(r_t_mode_pca,1:4)

pca_scores[1,]

lapply(cor(pca_scores[,1],A[,1])

#plot(subset(pca_IDRISI_s,1))

idrisi_dir <- "~/Google Drive/Papers_writing_MEOT/MEOT_analyses_02202015/workdir_terrset_08282015/anom_sst_1982_2007/components"
infile_idrisi_pca_s_mode <- "idrisi_pca_s_mode_PCA_Center_Std_S-Mode_DIF.xls"
infile_idrisi_pca_t_mode <- "idrisi_pca_t_mode_PCA_Center_Std_T-Mode_DIF.xls"
#infile_idrisi1_pca_t_mode <- file.path("/home/bparmentier/Google Drive/Papers_writing_MEOT/MEOT_analyses_02202015","pca_IDRISI_03302013_T-Mode_DIF.xlsx")
#infile_idrisi1_pca_s_mode <- file.path("/home/bparmentier/Google Drive/Papers_writing_MEOT/MEOT_analyses_02202015","pca_IDRISI_03302013_S-Mode_DIF.xlsx")

#idrisi_s_mode_loadings1 <- read.xls(infile_idrisi1_pca_s_mode, sheet = "loadings", header = TRUE)
#idrisi_s_mode_loadings1 <- read.xls(infile_idrisi1_pca_s_mode, sheet = "loadings", header = TRUE)

lf_s_mode_idrisi <- mixedsort(list.files(path=idrisi_dir,
                                         pattern="idrisi_pca_s_mode_PCA_Center_Std_S-Mode_EigenVec_.*.rst$",
                                         full.names=T))
r_s_eig <- stack(lf_s_mode_idrisi)
lf_s_mode_idrisi <- mixedsort(list.files(path=idrisi_dir,
                                         pattern="idrisi_pca_s_mode_PCA_Center_Std_S-Mode_CompLoadingImg_.*.rst$",
                                         full.names=T))
r_s_loadings <- stack(lf_s_mode_idrisi)

plot(r_s_mode_pca,1:4)
plot(r_t_mode_pca,1:4)
plot(r_s_eig,1:4)
plot(r_s_loadings,1:4)

idrisi_s_mode_loadings <- read.xls (file.path(idrisi_dir,infile_idrisi_pca_s_mode), sheet = "loadings", header = TRUE)
idrisi_t_mode_loadings <- read.xls (file.path(idrisi_dir,infile_idrisi_pca_t_mode), sheet = "loadings", header = TRUE)

plot(idrisi_s_mode_loadings[,2],Es$vectors[,1])
plot(idrisi_s_mode_loadings[,3],Es$vectors[,2])
plot(idrisi_s_mode_loadings[,4],Es$vectors[,3])

plot(idrisi_t_mode_loadings[,2],Et$vectors[,1])
plot(idrisi_t_mode_loadings[,3],Et$vectors[,2])
plot(idrisi_t_mode_loadings[,4],Et$vectors[,3])

plot(idrisi_t_mode_loadings[,2],loadings_Et1)
plot(idrisi_t_mode_loadings[,3],loadings_Et2)
#plot(idrisi_s_mode_loadings[,4],Es$vectors[,3])
cor(idrisi_t_mode_loadings[,2],loadings_Et1)
cor(idrisi_t_mode_loadings[,3],loadings_Et2)

cor(idrisi_s_mode_loadings[,2],Es$vectors[,1])
cor(idrisi_s_mode_loadings[,3],Es$vectors[,2])
cor(idrisi_s_mode_loadings[,4],Es$vectors[,3])

cor(idrisi_t_mode_loadings[,2],Et$vectors[,1])
cor(idrisi_t_mode_loadings[,3],Et$vectors[,2])
cor(idrisi_t_mode_loadings[,4],Et$vectors[,3])

cor(idrisi_t_mode_loadings[,2],principal_s_mode_df_unrotated[,1])
cor(idrisi_t_mode_loadings[,3],principal_s_mode_df_unrotated[,2])
cor(idrisi_t_mode_loadings[,4],principal_s_mode_df_unrotated[,3])

####### Comparison between EOT and rotation...
#ok now we should do this for different truncation...if one can show that the corr
#increases when the truncation decreases (more pc) then we will be ok...

#Make this a function...

#no_component <- 10
data_matrix <- A # data matrix before cross product
cross_matrix <- cAs #matrix used in the PCA
rotation_opt <- "none" #rotation used in PCA
fig_opt <- FALSE #gengerate figures?
out_suffix
#debug(run_pca_fun)
ref_data <- eot_s7_df[,-1]
#test <- run_pca_fun(A,mode="T",rotation_opt=rotation_opt,matrix_val=cAt,npc=no_component,loadings=TRUE,scores_opt=TRUE) 
list_no_component <- c(10,20,30,40,50,60,70,80,90,100)
n_cores <- 4

#l_loadings <- lapply(1:length(list_pca),function(i){list_pca[[i]]$loadings})

#debug(comparison_pca_eot_fun)
comp_pca_eot_s_unrotated_obj <- comparison_pca_eot_fun(list_no_component,data_matrix,rotation_opt,
                                           fig_opt,n_cores,out_suffix)
  

#diff could also be due to spatial average of 7 in EOT
#make a function for comparison...


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

