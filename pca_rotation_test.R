### DO PCA WITH SST_df: S-mode cross product

A<-SST_df[1:(ncol(SST_df)-2)]
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

### DO PCA WITH SST_df: T-mode cross product

A<-SST_df[1:(ncol(SST_df)-2)]
At<-scale(A)
cAt<-t(At)%*%At #T mode!!
cAt<-cAt/nrow(At)
Et<-eigen(cAt)

##Cross product PCA T mode
pca_SST<-principal(r=cAt, nfactors = npc, residuals = FALSE, covar=TRUE,rotate = "none")

sum(Et$value)
sum(diag(cAt))
sum(pca_SST$value)

#cAt<-A%*%At
#Ae<-eigen(cAt)
pca_SST_t_mode<-principal(r=cAt, nfactors = npc, residuals = FALSE, covar=TRUE,rotate = "none")
principal(r, nfactors = 1, residuals = FALSE,rotate="varimax",n.obs=NA, covar=FALSE, scores=FALSE,missing=FALSE,
          impute="median",oblique.scores=TRUE,method="regression",...)

##############
### 
require(psych) 
data(mtcars) 
rawd <- mtcars[,c("am","carb","cyl","disp","drat","gear","hp","mpg")] 

## compute acp 
.PC <- princomp(~am+carb+cyl+disp+drat+gear+hp+mpg, cor=TRUE, data=mtcars) 
pca <- principal(rawd, nfactors = 8, residuals = T, rotate="none", scores=T) 

## eigenvectors of the correlation matrix of raw data 
eigens <- eigen(cor(rawd)) 
eigens$vectors 
unclass(loadings(.PC))  # component 'loadings' in ?princomp parlance 
#not sure if ?principal outputs this 

## correlation matrix between raw data and unrotated component scores 
## 'loadings' in ?principal parlance and 'component matrix' in SPSS 
eigens$vectors %*% diag(sqrt(eigens$values)) 
cor(cbind(rawd, .PC$scores)) 
unclass(pca$loadings) 

## extract un-rotated scores of principal components 
head(scale(rawd) %*% eigens$vectors) # app, but very similar results 
head(.PC$scores) 
head(pca$scores) # scale'd, and obtained via regression on scale'd data 
head(factor.scores(scale(rawd), 
                   unclass(pca$loadings))) # same scores as from ?principal 
#for differeneces between ?princomp and ?principal scores 
#see last paragraph of Details in ?principal 
