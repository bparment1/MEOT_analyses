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



### DO PCA WITH SST_df: T-mode cross product from a standardized dataset...(i.e. correlation matrix)

A<-SST_df[1:(ncol(SST_df)-2)] #n*c or 42098*312
At<-scale(A) #center and reduce using sd and mean
sd(At[,2]) #this must be equal to 1 since At is standardized by columns
mean(At[,2]) #this must be equal to 0 since At is standardized by columns

cAt<-t(At)%*%At #T mode!! --> (c*n)by(n*c) should get c*c matrix or 312*312 
diag(cAt) #equal 1 since it is a correlation
#cAt<-cAt/nrow(At) #reduce by number of row...!!! THIS IS NOT NECESSARY SINCE At has been already based on standardized??
#note that cAt is standardized and has been divided by n, so it is equivalent to a correlation matrix
Et<-eigen(cAt)

##Cross product PCA T mode
#pca_SST<-principal(r=cAt, nfactors = npc, residuals = FALSE, covar=TRUE,rotate = "none")
pca_SST_t_mode<-principal(r=cAt, nfactors = npc, residuals = FALSE, covar=TRUE,rotate = "none")
unclass(pca_SST_t_mode$loadings) # extract the matrix of ??

sum(Et$value)
sum(diag(cAt))
sum(pca_SST_t_mode$value)

head(unclass(pca_SST_t_mode$loadings)[,1:npc])
head(Et$vectors[,1:npc]) #not equal...because not standardized...
Estd <-  Et$vectors%*% diag(sqrt(Et$values))
head(Estd)
head(unclass(pca_SST$loadings)[,1:npc])

##### Note that the scores are equal not equal but have correlation 1 and the ratio is the square root of eigenvalues...
PC_scores<-At%*%Et$vectors
pca_scores<-predict(pca_SST_t_mode,A)

cor(PC_scores[,1],pca_scores[,1]) 
plot(PC_scores[,1],pca_scores[,1]) 
PC_scores[1:10,1]/pca_scores[1:10,1]
sqrt(Et$values[1])

## DO PCA WITH SST_df: S-mode cross product

A<-SST_df[1:(ncol(SST_df)-2)] #remove x,y fields to obtain a data.frame with 42098*312 dimension
At<-t(A)
At<-scale(At)
A<-t(At)
cAs<-At%*%A
cAs<-cAs/nrow(cAs)
Es<-eigen(cAs) #cross product S-mode
diag(cAs)
##Cross product PCA T mode
pca_SST<-principal(r=cAs, nfactors = npc, residuals = FALSE, covar=TRUE,rotate = "none")

sum(Es$value)
sum(diag(cAs))
sum(pca_SST$value)



#cAt<-A%*%At
#Ae<-eigen(cAt)
#principal(r, nfactors = 1, residuals = FALSE,rotate="varimax",n.obs=NA, covar=FALSE, scores=FALSE,missing=FALSE,
#          impute="median",oblique.scores=TRUE,method="regression",...)



##Correlation PCA T mode
#pca_SST<-principal(r=cAt, nfactors = npc, residuals = FALSE, covar=TRUE,rotate = "none")

pca_SST_cor_t_mode<-principal(r=corAt, nfactors = npc, residuals = FALSE, covar=TRUE,rotate = "none")
unclass(pca_SST_t_mode$loadings) # extract the matrix of ??

sum(Et$value)
sum(diag(cAt))
sum(pca_SST_t_mode$value)

head(unclass(pca_SST_t_mode$loadings)[,1:npc])
head(Et$vectors[,1:npc]) #not equal...because not standardized...
Estd <-  Et$vectors%*% diag(sqrt(Et$values))
head(Estd)
head(unclass(pca_SST$loadings)[,1:npc])

#cAt<-A%*%At
#Ae<-eigen(cAt)
#principal(r, nfactors = 1, residuals = FALSE,rotate="varimax",n.obs=NA, covar=FALSE, scores=FALSE,missing=FALSE,
#          impute="median",oblique.scores=TRUE,method="regression",...)


##############
### 
require(psych) 
data(mtcars) 
rawd <- mtcars[,c("am","carb","cyl","disp","drat","gear","hp","mpg")] #subset the original data.frame

## compute acp 
.PC <- princomp(~am+carb+cyl+disp+drat+gear+hp+mpg, cor=TRUE, data=mtcars) 
pca <- principal(rawd, nfactors = 8, residuals = T, rotate="none", scores=T) 

## eigenvectors of the correlation matrix of raw data 
eigens <- eigen(cor(rawd)) #eigen value decomposition
eigens$vectors #matrix of the eigenvectors, eigens is a list
unclass(loadings(.PC))  # component 'loadings' in princomp parlance are the eigenvectors!!!

## correlation matrix between raw data and unrotated component scores 
## 'loadings' in ?principal parlance and 'component matrix' in SPSS --> loadings in princomp are the eigenvectors 
# components scores can be obtained by multiplying the original data by the eigenvectors...
eigens$vectors %*% diag(sqrt(eigens$values)) ##standardizing the eigenvectors...(it is already centered)
cor(cbind(rawd, .PC$scores))  #loadings are the correlation between the PC scores and the origina variables
unclass(pca$loadings) # extract the matrix of ??

## extract un-rotated scores of principal components 
head(scale(rawd) %*% eigens$vectors) # app, but very similar results 
head(.PC$scores) 
head(pca$scores) # scale'd, and obtained via regression on scale'd data 
head(factor.scores(scale(rawd), 
                   unclass(pca$loadings))) # same scores as from ?principal 
#for differeneces between ?princomp and ?principal scores 
#see last paragraph of Details in ?principal 
