############################      SCRIPT 2 MEOT-PCA          #######################################
#This script carries out a PCA with varimax rotation for a timeseries of 312 SST images.       #
#Note that spatial patterns from MEOT and MSSA components are not analyzed in this script       #                 
#AUTHOR: Benoit Parmentier                                                                      #
#DATE: 04/16/2013            
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
#library(forecast)                            # package containing ARIMA procedures
library(rasterVis)
library(psych)
library(GPArotation)

## Functions used in the script

### TO DO ... format the data for MSSA create a function for the window of any size ...use lag???
#unction ot lage files

load_obj <- function(f) 
{
  env <- new.env()
  nm <- load(f, env)[1]  
  env[[nm]]
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

### Parameters and argument

infile1<-"SAODI-01-1854_06-2011_test.asc"             #GHCN shapefile containing variables for modeling 2010                 
infile2<-"SAODI-01-1854_06-2011.csv"                     #List of 10 dates for the regression
infile3<-"MEOT_MSSA_Telcon_indices_12112012.xlsx"                        #LST dates name
infile4<-"mask_rgf_1_1.rst"
infile_pca <-"pca_IDRISI_03302013_T-Mode_DIF.xlsx"
infile_s_mode_pca  <- "pca_IDRISI_04152013_S-Mode_DIF.xlsx"
npc<-20 # number of pca to produce...put this at the beginning of the script...
lag_window<-13
telind<-c("PNA","NAO","TNA","TSA","SAOD","MEI","PDO","AO","AAO","AMM","AMOsm","QBO")
mode_list_MEOT<-c("MEOT1","MEOT3", "MEOT4","MEOT7","MEOT10","MEOT15","MEOT16")
mode_list_PCA<-c("MSSA1","MSSA2","MSSA3","MSSA4","MSSA5","MSSA6")    
#mode_list_PCA<-paste("MSSA",1:15,sep="")
out_prefix<-"pca_test_04152013"
#out_prefix<-"pca_test_03302013"

# on Benoit Mac
#in_path<-"/Users/benoitparmentier/Dropbox/Data/MEOT_paper/MEOT12272012/MEOT_working_dir_03102013"
# on Atlas:
in_path<-"/home/parmentier/Data/MEOT12272012/MEOT_working_dir_03102013/"

setwd(in_path)

########################################

### STEP 0: READ IN DATASETS RELATED TO TELECONNECTION AND PREVIOUS LOADINGS

SAODI<-read.table(infile2,sep=",", header=TRUE)
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

#Creating time series objects
d_ts<-ts(data=dat,start=c(1982,1), end=c(2006,12), frequency=12, names=names(dat))
colnames(d_ts)
d_z<-as.zoo(d_ts)  #### THIS IS THE TIME SERIES OBJECT USED LATER ON

##############################################################
#### STEP 1: load image time series data, apply latitude weight for area

write.table(in_SAODI,file=paste("SAOD_index_1981_2006",out_prefix,".txt",sep=""),sep=",")

#Import results from MEOT and MSSA analyses with teleconneciton indices
dat<-read.xls(infile3, sheet=1)
tail(dat$Date_label)
dat$SAOD<-in_SAODI  #Adding the SAOD index to all the MEOT/MSSA results in the data frame

#Creating time series objects
d_ts<-ts(data=dat,start=c(1982,1), end=c(2006,12), frequency=12, names=names(dat))
colnames(d_ts)
d_z<-as.zoo(d_ts)  #### THIS IS THE TIME SERIES OBJECT USED LATER ON

##############################################################
#### STEP 1: load image time series data, apply latitude weight for area

lf<-mixedsort(list.files(pattern="ANOM_SST.*.rst"))
SST_s<-stack(lf)
SST1<-raster(lf[1])

mask_land<-raster(infile4)
mask_land_NA<-mask_land
mask_land_NA[mask_land_NA==0]<-NA
SST_rast<-mask(SST_s,mask_land_NA,filename="ANOM_SST_1982_2007.tif",overwrite=TRUE)
class(SST_rast)

#Multiply layer by weight

lat_coord<-coordinates(mask_land)[,2]
w_rast<-mask_land
values(w_rast)<-cos(lat_coord*pi/180) #use area in raster package for weight of a cell??
w_rast_m<-mask(w_rast,mask_land_NA)

SST_rast<-SST_rast*w_rast #Maybe use overlay to avoid putting in memory?
writeRaster(SST_rast,filename="ANOM_SST_1982_2007_weighted.tif",overwrite=TRUE)
rm(SST_rast)
SST_rast<-brick("ANOM_SST_1982_2007_weighted.tif")
#Do all in one step using overlay
#SST_rast<-overlay(x=SST_rast,y=w_rast,filename="ANOM_SST_1982_2007_weighted.tif",overwrite=TRUE,
#                   fun=function(x,y) {(return(x*y)})

SST1_m<-subset(SST_rast,1)
plot(stack(SST1,SST1_m))

##############################################################
#### STEP 2: get data ready for PCA by lagging data frame...

SST_sgdf<-as(SST_rast,"SpatialGridDataFrame")
#SST_df<-as.data.frame(SST_rast) #drop the coordinates x and y
SST_df<-as.data.frame(SST_sgdf)
rm(SST_sgdf)

#nv<-312
SST_xy<-SST_df[,c("s1","s2")]

### DO PCA WITH SST_df: T-mode cross product from a standardized dataset...(i.e. correlation matrix)

A<-SST_df[1:(ncol(SST_df)-2)] #n*c or 42098*312
At<-scale(A) #center and reduce using sd and mean
sd(At[,2]) #this must be equal to 1 since At is standardized by columns
mean(At[,2]) #this must be equal to 0 since At is standardized by columns
names_var<-paste("t",1:ncol(At),sep="_")
colnames(At)<-names_var
cAt<-t(At)%*%At #T mode!! --> (c*n)by(n*c) should get c*c matrix or 312*312 
diag(cAt) #equal 1 since it is a correlation, only if divided by n-1!!!
cAt<-cAt/(nrow(At)-1)
diag(cAt)
#cAt<-cAt/nrow(At) #reduce by number of row...!!! THIS IS NOT NECESSARY SINCE At has been already based on standardized??
#note that cAt is standardized and has been divided by n, so it is equivalent to a correlation matrix
Et<-eigen(cAt)

##Cross product PCA T mode
#pca_SST<-principal(r=cAt, nfactors = npc, residuals = FALSE, covar=TRUE,rotate = "none")
pca_SST_t_mode<-principal(r=cAt, nfactors = npc, residuals = FALSE, covar=TRUE,rotate = "none")
unclass(pca_SST_t_mode$loadings) # extract the matrix of ??

(nrow(At)-1)*ncol(At)  #Sum of the diagonal is 42097*312
sum(Et$value)   
sum(diag(cAt))
sum(pca_SST_t_mode$value)
mean(Et$vectors)  #This is equal to zero since there was centering!!
1/var(Et$vectors)[1]
var(Et$vectors)[1]
var(Et$vectors)[2]
var(Et$vectors)[3]

#
pca_IDRISI_covar<-read.xls(infile_pca, sheet=1)
pca_IDRISI_cor<-read.xls(infile_pca, sheet=2)
pca_IDRISI_eigenvalues<-read.xls(infile_pca, sheet=3)
pca_IDRISI_eigenvectors<-read.xls(infile_pca, sheet=4)
pca_IDRISI_loadings<-read.xls(infile_pca, sheet=5)

###### CREATE SCORES IMAGES USING PREDICTED SCORES FROM PRINCIPAL
pca_scores<-predict(pca_SST_t_mode,A)  #generate scores from original matrix and object
pca_scores_spdf<-cbind(as.data.frame(pca_scores),SST_xy) #add coordinates
coordinates(pca_scores_spdf)<-pca_scores_spdf[,c("s1","s2")] #promote to spdf

#Now generate raster images...
pca_scores_lf<-pca_to_raster_fun(pca_scores_spdf,SST1_m,-9999,out_prefix)

pca_SST_s<-stack(pca_scores_lf)
NAvalue(pca_SST_s)<- -9999
#plot(subset(pca_IDRISI_s,1))

pca_IDRISI_s<-stack(mixedsort(list.files(pattern="pca_IDRISI_03302013_T-Mode_Cmp.*.RST$")))
pca_IDRISI_s<-mask(pca_IDRISI_s,mask_land_NA)
k=3
plot(subset(pca_IDRISI_s,k))
plot(subset(pca_SST_s,k))
plot(stack(subset(pca_SST_s,k),subset(pca_IDRISI_s,k)))
plot(subset(pca_SST_s,k),subset(pca_IDRISI_s,k))

## EXAMINE CORRELATION MAP SCORES BETWEEN PCA FROM R AND FROM IDRISI...

mat1_df<-as.data.frame(as(pca_SST_s,"SpatialGridDataFrame"))
mat2_df<-as.data.frame(as(pca_IDRISI_s,"SpatialGridDataFrame"))
cormat<-cor(mat1_df,mat2_df)
diag(cormat) #This is the diagonal of correlation matrix for IDRISI scores and R scores

# EXAMINE LOADINGS???

pca_IDRISI_loadings$CMP.1
#
#mssa_names <-paste("MSSA",1:15,sep="")
#mssa_mat <- d_z[,mssa_names]
pc_mat <-unclass(pca_SST_t_mode$loadings)[,1:15]
cor(pc_mat,pca_IDRISI_loadings[,1:15])
plot(diag(cor(pc_mat,pca_IDRISI_loadings[,1:15])),type="h")
### COMPARE MSSA IDRISI AND R components...
for (j in 1:15){
  #png
  j<-j+1
  plot(pc_mat[,j],type="l")
  par(new=TRUE)
  plot(pca_IDRISI_loadings[,j],type="l",col="red",axes=FALSE)
  axis(4) 
  legend("topright",col=c("black","red"),
         legend=c(colnames(pc_mat)[j],colnames(pca_IDRISI_loadings)[j]),
         lwd=1.5,lty=c(1,1))
  title(paste(colnames(pc_mat)[j]," and ",colnames(pca_IDRISI_loadings)[j]))
  #dev.off()
}

####################################################################
### NOW DO PCA WITH S-mode

## DO PCA WITH SST_df: S-mode cross product
#PC_scores[1:10,1]/pca_scores[1:10,1]
A<-SST_df[1:(ncol(SST_df)-2)] #remove x,y fields to obtain a data.frame with 42098*312 dimension
Atr<-t(A)
Ats<-scale(Atr) #center time series...
As<-t(Ats)
cAs<-Ats%*%As
cAs<-cAs/nrow(cAs)
Es<-eigen(cAs) #cross product S-mode
diag(cAs)

##Cross product PCA T mode
#pca_SST<-principal(r=cAs, nfactors = npc, residuals = FALSE, covar=TRUE,rotate = "none")
pca_SST_s_mode<-principal(r=cAs, nfactors = npc, residuals = FALSE, covar=TRUE,rotate = "none")
unclass(pca_SST_s_mode$loadings) # extract the matrix of ??

sum(Es$value)
sum(diag(cAs)) #sum of eigenvalue is equal to the number of variable...
sum(pca_SST_s_mode$value)

## NOW predict raster images...using eigen method

PC_scores_e_s_mode<-As%*%Es$vectors # this generates 312 component scores....

PC_scores_e_s_mode_scores_spdf<-cbind(as.data.frame(PC_scores_e_s_mode),SST_xy) #add coordinates
coordinates(PC_scores_e_s_mode_scores_spdf)<-PC_scores_e_s_mode_scores_spdf[,c("s1","s2")] #promote to spdf
#Now generate raster images...
out_prefix_e_s_mode<- paste("e_s_mode_",out_prefix,sep="")
PC_scores_e_s_mode_scores_lf<-pca_to_raster_fun(PC_scores_e_s_mode_scores_spdf,SST1_m,-9999,out_prefix_e_s_mode)

Es_prin<-princomp(cor = FALSE,scores=FALSE,covmat=cAs) #Ok this gives the same as Es$vectors
Es_prin_v <-varimax(Es_prin$loadings[,1:20],normalize= FALSE,eps= 1e-3)


## NOW predict raster images...using principal method

pca_s_mode_scores<-predict(pca_SST_s_mode,A)  #generate scores from original matrix and object
pca_s_mode_scores_spdf<-cbind(as.data.frame(pca_s_mode_scores),SST_xy) #add coordinates
coordinates(pca_s_mode_scores_spdf)<-pca_s_mode_scores_spdf[,c("s1","s2")] #promote to spdf

#Now generate raster images...
out_prefix_s_mode<- paste("s_mode_",out_prefix,sep="")
pca_s_mode_scores_lf<-pca_to_raster_fun(pca_s_mode_scores_spdf,SST1_m,-9999,out_prefix_s_mode)

pca_SST_s_mode_rast <-stack(pca_s_mode_scores_lf) #raster stack of images...
NAvalue(pca_SST_s_mode_rast)<- -9999
j<-1
plot(raster(pca_SST_s_mode_rast,j))
## NOW read loadings and eigenvalues from PCA in IDRISI to compare 

pca_s_mode_IDRISI_loadings<-read.xls(infile_s_mode_pca, sheet=1) 
#note scores are loadings in S-mode
#pca_s_mode_IDRISI_eigenvalues<-read.xls(infile_s_mode_pca, sheet=3)
pca_s_mode_IDRISI_loadings
pc_mat <-unclass(pca_SST_s_mode$loadings)[,1:15]
cor(pc_mat,pca_IDRISI_loadings[,1:15])
plot(diag(cor(pc_mat,pca_s_mode_IDRISI_loadings[,1:15])),type="h")
### COMPARE MSSA IDRISI AND R components...
for (j in 1:15){
  #png
  j<-j+1
  plot(pc_mat[,j],type="l")
  par(new=TRUE)
  plot(pca_IDRISI_loadings[,j],type="l",col="red",axes=FALSE)
  axis(4) 
  legend("topright",col=c("black","red"),
         legend=c(colnames(pc_mat)[j],colnames(pca_IDRISI_loadings)[j]),
         lwd=1.5,lty=c(1,1))
  title(paste(colnames(pc_mat)[j]," and ",colnames(pca_IDRISI_loadings)[j]))
  #dev.off()
}


## NOW load raster scores from PCA in IDRISI
pca_s_mode_IDRISI_s<-stack(mixedsort(list.files(pattern="pca_IDRISI_04152013_S-Mode_Cmp.*.RST$")))
pca_IDRISI_04152013_S-Mode_CompLoadingImg_.*.rst
pca_s_mode_IDRISI_s<-mask(pca_s_mode_IDRISI_s,mask_land_NA)
k=3
plot(subset(pca_IDRISI_s,k))
plot(subset(pca_SST_s,k))
plot(stack(subset(pca_SST_s,k),subset(pca_IDRISI_s,k)))
plot(subset(pca_SST_s,k),subset(pca_IDRISI_s,k))

####### PART 3 COMPARING METHODS TO PRODUCE COMPONENTS USING PRINCIPAL AND LINEAR ALGEBRA

###COMPARING EIGENVECTORS FROM "eigen" and principal function; Tmode

head(unclass(pca_SST_t_mode$loadings)[,1:npc])
head(Et$vectors[,1:npc]) #not equal...because not standardized...
cor(unclass(pca_SST_t_mode$loadings)[,1:npc],Et$vectors[,1:npc])
diag(cor(unclass(pca_SST_t_mode$loadings)[,1:npc],Et$vectors[,1:npc]))
plot(unclass(pca_SST_t_mode$loadings)[,1],Et$vectors[,1])
unclass(pca_SST_t_mode$loadings)[,1]/Et$vectors[,1]
##NOTE THAT THE RATIO IS 6.989 which is the square root of eigenvalue 1!!!
(unclass(pca_SST_t_mode$loadings)[,1]/Et$vectors[,1])/sqrt(Et$value[1])

Estd <-  Et$vectors%*% diag(sqrt(Et$values))  # Assign length to each eigenvectors??
head(Estd[,1])
head(unclass(pca_SST_t_mode$loadings)[,1])

###COMPARIING SCORES FROM FROM "eigen" and principal function
## Note that the scores are equal only if the square root of eigenvalues has been used to reduce (normalize) the Eigen Vectors!!!

PC_scores<-At%*%Et$vectors # to
PC_scores<-as.matrix(A)%*%Et$vectors # to

pca_scores<-predict(pca_SST_t_mode,A)

cor(PC_scores[,1],pca_scores[,1]) 
plot(PC_scores[,1],pca_scores[,1]) 
PC_scores[1:10,1]/pca_scores[1:10,1]
sqrt(Et$values[1])
(PC_scores[1:10,1]/sqrt(Et$values[1]))/pca_scores[1:10,1]
plot(PC_scores[,1]/sqrt(Et$values[1]),pca_scores[,1]) 

mean(PC_scores[,1])
sd(PC_scores[,1])
mean(pca_scores[,1]) #This has been standardize!!! Principal standardizes the values!!
sd(pca_scores[,1])

##SO PCA scores from principal are standardized...!!! i.e. they have been divided by the sqrt of the eigenvalues???

#Now reduced version
PC_std_scores<-At%*%Estd
plot(PC_std_scores[,1],pca_scores[,1])
PC_std_scores[1:10,1]/pca_scores[1:10,1]

(PC_std_scores[1:10,1]/(Et$values[1]))/pca_scores[1:10,1]


#varimax(x, normalize = TRUE, eps = 1e-5)

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