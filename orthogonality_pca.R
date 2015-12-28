
#Orthogonal Vectors
#Two vectors u and v whose dot product is u·v=0 (i.e., the vectors are perpendicular) are said to be orthogonal. 
#In three-space, three vectors can be mutually perpendicular.

k <- 5
tstMat <- array(runif(k), dim=c(k,k))
tstOrth <- qr.Q(qr(tstMat))
t(tstOrth)%*%tstOrth

#http://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variable
n     <- 20                    # length of vector
rho   <- 0.6                   # desired correlation = cos(angle)
theta <- acos(rho)             # corresponding angle
x1    <- rnorm(n, 1, 1)        # fixed given data
x2    <- rnorm(n, 2, 0.5)      # new random data
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


n     <- 20                    # length of vector
rho   <- 0.6                   # desired correlation = cos(angle)
theta <- acos(rho)             # corresponding angle
x1    <- rnorm(n, 1, 1)        # fixed given data
x2    <- rnorm(n, 2, 0.5)      # new random data
X     <- cbind(x1, x2)         # matrix
Xctr  <- scale(X, center=TRUE, scale=FALSE)   # centered columns (mean 0)

Id   <- diag(n)                               # identity matrix
Q    <- qr.Q(qr(Xctr[ , 1, drop=FALSE]))      # QR-decomposition, just matrix Q
P    <- tcrossprod(Q)          # = Q Q'       # projection onto space defined by x1
x2o  <- (Id-P) %*% Xctr[ , 2]                 # x2ctr made orthogonal to x1ctr
Xc2  <- cbind(Xctr[ , 1], x2o)                # bind to matrix
Y    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  # scale columns to length 1

### Generate orhogonal solution...

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

list_seed_nb <- c(100,101)
nt<- 100
rho_val <-0
debug()
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

