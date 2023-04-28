#Bootstrap implementation used in this case is based on Residual Resampling
# Steps involved
# (1) estimation of the initial estimation of y.diff
# (2) use the residuals from (1) to calculate the restricted residuals based on S.tilde
# (3) draw a random sample from the residuals with replacement with the lenght of the original residuals - sample() available in R
# (4) test regression with the pseudo sample y*, which was created with the resampled residuals
# (5) repeat (3) and (4) n.times and save each resulting t-stat from the repetitions
# 
# Bootstrap procedure needs to be done for values of y*, but also the minimum LM stat 
# is calculated for each series, so the test has to be applied to all replicated y* time series.
# Bootstrap the distribution of the t-statistic based on the determined break dates

#-----------------------------------------------------------------------------
# Set the number of replications
#-----------------------------------------------------------------------------

n.sim <- 1000
set.seed(123)

#-----------------------------------------------------------------------------
# Initialize the parallelization framework
#-----------------------------------------------------------------------------
# 
library(foreach)
library(doSNOW)
library(parallel)
source("LeeStrazicichUnitRootTestParallelization.R")
load("Methodology data2.RData")


#######################################
# IQRs -----------------------------------
######################################

cl <- makeCluster(max(1, detectCores() - 2))

registerDoSNOW(cl)
#Define variable y, which is the variable to be analysed

myVariable <- as.numeric(IQR_FSZ)

lagmatrix <- function(x,max.lag){
  embed(c(rep(NA,max.lag),x),max.lag+1)
}
#Add diffmatrix function
diffmatrix <- function(x,max.diff = 1,max.lag = 1){
        #Add if condition to make it possible to differentiate between matrix and vector                  
        if(is.vector(x) == TRUE ){
                myx <- embed(c(rep(NA,max.lag),diff(x,max.lag,max.diff)),max.diff)
                colnames(myx) <- paste("v1.d",max.diff, sep=".")
                return(myx)
        }
        
        else if(is.matrix(x) == TRUE){
                myx <- rbind(matrix(rep(NA,max.lag), ncol = ncol(x)), matrix(diff(x,max.lag,max.diff), ncol = ncol(x)))
                mycolnames <- colnames(x)
                colnames(myx) <- paste(mycolnames,"d",max.diff, sep=".")
                return(myx)
        }
        #if matrix the result should be 0, if only a vector it should be 1
        else if(as.integer(is.null(ncol(x))) == 0 ){
                myx <- rbind(matrix(rep(NA,max.lag), ncol = ncol(x)), matrix(diff(x,max.lag,max.diff), ncol = ncol(x)))
                mycolnames <- colnames(x)
                colnames(myx) <- paste(mycolnames,"d",max.diff, sep=".")
                return(myx)
                
        }
}

n <- length(myVariable)
trend <- 1:n
breaks <- 2
model <- "crash"
myVariable.diff <- diffmatrix(myVariable, max.diff = 1, max.lag = 1)
# Run the test once to determine the values for the break dates
# Definition of all the options needed

test.results.y <-  ur.ls.bootstrap(y=myVariable , model = model, breaks = breaks , method = "GTOS",pn = 0.1, critval = "bootstrap", print.results = "print")


# Define variable Z, which is either Dt or DTt, depending on the model used
# Create Dt and DTt, depending on the breaks, found in the initial application of the LS test

#Define break dates
myBreak1 <- test.results.y[[1]]$`First break`
myBreak2 <- test.results.y[[1]]$`Second break`
  

Dt1 <-  as.matrix(cbind(trend, trend >= (myBreak1 + 1)))

#       Dummy with break in intercept and in trend

DTt1 <- as.matrix(cbind(Dt1, c(rep(0, myBreak1), 1:(n - myBreak1))))
colnames(Dt1) <- c("Trend","D")
colnames(DTt1) <- c("Trend","D","DTt")
Dt <- cbind(Dt1)
DTt <- cbind(DTt1)


Dt2 <-  as.matrix(trend >= (myBreak2 + 1))
DTt2 <- as.matrix(cbind(Dt2, c(rep(0, myBreak2), 1:(n - myBreak2))))
colnames(Dt2) <- c("D2")
colnames(DTt2) <- c("D2","DTt2")
#print(paste("Break2: ",myBreak2, sep = ""))

#Combine all Dummies into one big matrix to make it easier to include in the regressions
if(breaks == 1)
{
Dt <- cbind(Dt1)
DTt <- cbind(DTt1)
} else if(breaks == 2){
  
  Dt <- cbind(Dt1, Dt2)
  DTt <- cbind(DTt1, DTt2)
}

if(model == "crash"){
  Z <- Dt
} else if(model == "break"){
    
  Z <- DTt
}

Z.diff <- diffmatrix(Z, max.diff = 1, max.lag = 1)
#Matrix of coefficients, which omits the NA values
myZcoef <- na.omit(coef(lm(myVariable.diff ~ Z.diff)))

# delta_1 from equation (7) of Chou 2007, which is needed for the manual calculation of S.tilde according to Chou 2007
# equivalent results to the calculation with the cumulated sum
delta_1 <- as.numeric(myVariable[1] - (Z[1,] %*% myZcoef))

#Result is equivalent to the original S.tilde
myS.tilde <- myVariable - delta_1 - Z %*% myZcoef
if(model == "crash"){
y.star <- delta_1 + as.vector((Dt %*% myZcoef)) + sample(myS.tilde, size = length(myS.tilde), replace = TRUE)
} else if(model == "break"){
  y.star <- delta_1 + as.vector((DTt %*% myZcoef)) + sample(myS.tilde, size = length(myS.tilde), replace = TRUE)
  
}


#Loop for the actual boot strapping procedure

boot.y <- foreach(i=1:n.sim, .combine = rbind, .packages = 'foreach') %dopar%{
  y.star <- delta_1 + as.vector((Z %*% myZcoef)) + sample(myS.tilde, size = length(myS.tilde), replace = TRUE)
  #y.star.matrix[,i] <- y.star
  #plotl(y.star)
  test.stat <- ur.ls.bootstrap(y=y.star , model = model, breaks = breaks , lags = 10, method = "Fixed",pn = 0.1, critval = "bootstrap", print.results = "silent")
  if(breaks == 2){
  return(list(c(unlist(test.stat[[1]]$`t-stat`)), c(unlist(test.stat[[1]]$`First break`)) ,c(unlist(test.stat[[1]]$`Second break`)),  c(as.numeric(test.stat[[1]]$Runtime, units = "mins"))))
  } else if (breaks == 1){
    return(c(c(unlist(test.stat[[1]]$`t-stat`)), c(unlist(test.stat[[1]]$`First break`)) , c(as.numeric(test.stat[[1]]$Runtime, units = "mins"))))
    }
}


if(breaks == 2){
      colnames(boot.y) <- c("t-stat","First break", "Second break", "Runtime")
} else if (breaks == 1){
      colnames(boot.y) <- c("t-stat","First break",  "Runtime")
}


#Calculate the critical values, based on the distribution of the test statistic
hist(unlist(boot.y[,"t-stat"]), freq = FALSE)
quantile(unlist(boot.y[,"t-stat"]), c(.01, .025 ,.05, .1, .25, .5, .75, .95, .975,.99))
quantile(unlist(boot.y[,"First break"]), c(.01, .025 ,.05, .1, .25, .5, .75, .95, .975,.99))
quantile(unlist(boot.y[,"Second break"]), c(.01, .025 ,.05, .1, .25, .5, .75, .95, .975,.99))

hist(unlist(boot.y[,"Second break"]), freq = F)
hist(unlist(boot.y[,"First break"]))

IQR_BOOTS_RESULTS_crash <- test.results.y
IQR_BOOTS_INTERVALS_crash <- boot.y




#######################################
# MPDs -----------------------------------
######################################

cl <- makeCluster(max(1, detectCores() - 2))

registerDoSNOW(cl)
#Define variable y, which is the variable to be analysed

myVariable <- as.vector(as.numeric(na.omit(MPD_MAFSZ)))

lagmatrix <- function(x,max.lag){
  embed(c(rep(NA,max.lag),x),max.lag+1)
}
#Add diffmatrix function
diffmatrix <- function(x,max.diff = 1,max.lag = 1){
  #Add if condition to make it possible to differentiate between matrix and vector                  
  if(is.vector(x) == TRUE ){
    myx <- embed(c(rep(NA,max.lag),diff(x,max.lag,max.diff)),max.diff)
    colnames(myx) <- paste("v1.d",max.diff, sep=".")
    return(myx)
  }
  
  else if(is.matrix(x) == TRUE){
    myx <- rbind(matrix(rep(NA,max.lag), ncol = ncol(x)), matrix(diff(x,max.lag,max.diff), ncol = ncol(x)))
    mycolnames <- colnames(x)
    colnames(myx) <- paste(mycolnames,"d",max.diff, sep=".")
    return(myx)
  }
  #if matrix the result should be 0, if only a vector it should be 1
  else if(as.integer(is.null(ncol(x))) == 0 ){
    myx <- rbind(matrix(rep(NA,max.lag), ncol = ncol(x)), matrix(diff(x,max.lag,max.diff), ncol = ncol(x)))
    mycolnames <- colnames(x)
    colnames(myx) <- paste(mycolnames,"d",max.diff, sep=".")
    return(myx)
    
  }
}

n <- length(myVariable)
trend <- 1:n
breaks <- 2
model <- "break"
myVariable.diff <- diffmatrix(myVariable, max.diff = 1, max.lag = 1)
# Run the test once to determine the values for the break dates
# Definition of all the options needed

test.results.y <-  ur.ls.bootstrap(y=myVariable , model = model, breaks = breaks, lags = 11, method = "Fixed",pn = 0.1, critval = "bootstrap", print.results = "print")


# Define variable Z, which is either Dt or DTt, depending on the model used
# Create Dt and DTt, depending on the breaks, found in the initial application of the LS test


#Define break dates
myBreak1 <- test.results.y[[1]]$`First break`
myBreak2 <- test.results.y[[1]]$`Second break`


Dt1 <-  as.matrix(cbind(trend, trend >= (myBreak1 + 1)))

#       Dummy with break in intercept and in trend

DTt1 <- as.matrix(cbind(Dt1, c(rep(0, myBreak1), 1:(n - myBreak1))))
colnames(Dt1) <- c("Trend","D")
colnames(DTt1) <- c("Trend","D","DTt")
Dt <- cbind(Dt1)
DTt <- cbind(DTt1)


Dt2 <-  as.matrix(trend >= (myBreak2 + 1))
DTt2 <- as.matrix(cbind(Dt2, c(rep(0, myBreak2), 1:(n - myBreak2))))
colnames(Dt2) <- c("D2")
colnames(DTt2) <- c("D2","DTt2")
#print(paste("Break2: ",myBreak2, sep = ""))

#Combine all Dummies into one big matrix to make it easier to include in the regressions
if(breaks == 1)
{
  Dt <- cbind(Dt1)
  DTt <- cbind(DTt1)
} else if(breaks == 2){
  
  Dt <- cbind(Dt1, Dt2)
  DTt <- cbind(DTt1, DTt2)
}

if(model == "crash"){
  Z <- Dt
} else if(model == "break"){
  
  Z <- DTt
}

Z.diff <- diffmatrix(Z, max.diff = 1, max.lag = 1)
#Matrix of coefficients, which omits the NA values
myZcoef <- na.omit(coef(lm(myVariable.diff ~ Z.diff)))

# delta_1 from equation (7) of Chou 2007, which is needed for the manual calculation of S.tilde according to Chou 2007
# equivalent results to the calculation with the cumulated sum
delta_1 <- as.numeric(myVariable[1] - (Z[1,] %*% myZcoef))

#Result is equivalent to the original S.tilde
myS.tilde <- myVariable - delta_1 - Z %*% myZcoef
if(model == "crash"){
  y.star <- delta_1 + as.vector((Dt %*% myZcoef)) + sample(myS.tilde, size = length(myS.tilde), replace = TRUE)
} else if(model == "break"){
  y.star <- delta_1 + as.vector((DTt %*% myZcoef)) + sample(myS.tilde, size = length(myS.tilde), replace = TRUE)
  
}


#Loop for the actual boot strapping procedure

boot.y <- foreach(i=1:n.sim, .combine = rbind, .packages = 'foreach') %dopar%{
  y.star <- delta_1 + as.vector((Z %*% myZcoef)) + sample(myS.tilde, size = length(myS.tilde), replace = TRUE)
  #y.star.matrix[,i] <- y.star
  #plotl(y.star)
  test.stat <- ur.ls.bootstrap(y=y.star , model = model, breaks = breaks, lag = 11, method = "Fixed",pn = 0.1, critval = "bootstrap", print.results = "silent")
  if(breaks == 2){
    return(list(c(unlist(test.stat[[1]]$`t-stat`)), c(unlist(test.stat[[1]]$`First break`)) ,c(unlist(test.stat[[1]]$`Second break`)),  c(as.numeric(test.stat[[1]]$Runtime, units = "mins"))))
  } else if (breaks == 1){
    return(c(c(unlist(test.stat[[1]]$`t-stat`)), c(unlist(test.stat[[1]]$`First break`)) , c(as.numeric(test.stat[[1]]$Runtime, units = "mins"))))
  }
}


if(breaks == 2){
  colnames(boot.y) <- c("t-stat","First break", "Second break", "Runtime")
} else if (breaks == 1){
  colnames(boot.y) <- c("t-stat","First break",  "Runtime")
}



MPD_BOOTS_RESULTS <- test.results.y
MPD_BOOTS_INTERVALS <- boot.y
#Calculate the critical values, based on the distribution of the test statistic
hist(unlist(MPD_BOOTS_INTERVALS[,"t-stat"]), freq = FALSE)
quantile(unlist(MPD_BOOTS_INTERVALS[,"t-stat"]), c(.01, .025 ,.05, .1, .25, .5, .75, .95, .975,.99))
quantile(unlist(MPD_BOOTS_INTERVALS[,"First break"]), c(.01, .025 ,.05, .1, .25, .5, .75, .95, .975,.99))
quantile(unlist(MPD_BOOTS_INTERVALS[,"Second break"]), c(.01, .025 ,.05, .1, .25, .5, .75, .95, .975,.99))

hist(unlist(MPD_BOOTS_INTERVALS[,"Second break"]), freq = F)
hist(unlist(MPD_BOOTS_INTERVALS[,"First break"]), freq = F)







#######################################
# MPMs -----------------------------------
######################################

cl <- makeCluster(max(1, detectCores() - 2))

registerDoSNOW(cl)
#Define variable y, which is the variable to be analysed

myVariable <- as.vector(as.numeric(na.omit(MPM_MAFSZ)))

lagmatrix <- function(x,max.lag){
  embed(c(rep(NA,max.lag),x),max.lag+1)
}
#Add diffmatrix function
diffmatrix <- function(x,max.diff = 1,max.lag = 1){
  #Add if condition to make it possible to differentiate between matrix and vector                  
  if(is.vector(x) == TRUE ){
    myx <- embed(c(rep(NA,max.lag),diff(x,max.lag,max.diff)),max.diff)
    colnames(myx) <- paste("v1.d",max.diff, sep=".")
    return(myx)
  }
  
  else if(is.matrix(x) == TRUE){
    myx <- rbind(matrix(rep(NA,max.lag), ncol = ncol(x)), matrix(diff(x,max.lag,max.diff), ncol = ncol(x)))
    mycolnames <- colnames(x)
    colnames(myx) <- paste(mycolnames,"d",max.diff, sep=".")
    return(myx)
  }
  #if matrix the result should be 0, if only a vector it should be 1
  else if(as.integer(is.null(ncol(x))) == 0 ){
    myx <- rbind(matrix(rep(NA,max.lag), ncol = ncol(x)), matrix(diff(x,max.lag,max.diff), ncol = ncol(x)))
    mycolnames <- colnames(x)
    colnames(myx) <- paste(mycolnames,"d",max.diff, sep=".")
    return(myx)
    
  }
}

n <- length(myVariable)
trend <- 1:n
breaks <- 2
model <- "crash"
myVariable.diff <- diffmatrix(myVariable, max.diff = 1, max.lag = 1)
# Run the test once to determine the values for the break dates
# Definition of all the options needed

test.results.y <-  ur.ls.bootstrap(y=myVariable , model = model, breaks = breaks, method = "GTOS",pn = 0.1, critval = "bootstrap", print.results = "print")


# Define variable Z, which is either Dt or DTt, depending on the model used
# Create Dt and DTt, depending on the breaks, found in the initial application of the LS test

#Define break dates
myBreak1 <- test.results.y[[1]]$`First break`
myBreak2 <- test.results.y[[1]]$`Second break`


Dt1 <-  as.matrix(cbind(trend, trend >= (myBreak1 + 1)))

#       Dummy with break in intercept and in trend

DTt1 <- as.matrix(cbind(Dt1, c(rep(0, myBreak1), 1:(n - myBreak1))))
colnames(Dt1) <- c("Trend","D")
colnames(DTt1) <- c("Trend","D","DTt")
Dt <- cbind(Dt1)
DTt <- cbind(DTt1)


Dt2 <-  as.matrix(trend >= (myBreak2 + 1))
DTt2 <- as.matrix(cbind(Dt2, c(rep(0, myBreak2), 1:(n - myBreak2))))
colnames(Dt2) <- c("D2")
colnames(DTt2) <- c("D2","DTt2")
#print(paste("Break2: ",myBreak2, sep = ""))

#Combine all Dummies into one big matrix to make it easier to include in the regressions
if(breaks == 1)
{
  Dt <- cbind(Dt1)
  DTt <- cbind(DTt1)
} else if(breaks == 2){
  
  Dt <- cbind(Dt1, Dt2)
  DTt <- cbind(DTt1, DTt2)
}

if(model == "crash"){
  Z <- Dt
} else if(model == "break"){
  
  Z <- DTt
}

Z.diff <- diffmatrix(Z, max.diff = 1, max.lag = 1)
#Matrix of coefficients, which omits the NA values
myZcoef <- na.omit(coef(lm(myVariable.diff ~ Z.diff)))

# delta_1 from equation (7) of Chou 2007, which is needed for the manual calculation of S.tilde according to Chou 2007
# equivalent results to the calculation with the cumulated sum
delta_1 <- as.numeric(myVariable[1] - (Z[1,] %*% myZcoef))

#Result is equivalent to the original S.tilde
myS.tilde <- myVariable - delta_1 - Z %*% myZcoef
if(model == "crash"){
  y.star <- delta_1 + as.vector((Dt %*% myZcoef)) + sample(myS.tilde, size = length(myS.tilde), replace = TRUE)
} else if(model == "break"){
  y.star <- delta_1 + as.vector((DTt %*% myZcoef)) + sample(myS.tilde, size = length(myS.tilde), replace = TRUE)
  
}


#Loop for the actual boot strapping procedure

boot.y <- foreach(i=1:n.sim, .combine = rbind, .packages = 'foreach') %dopar%{
  y.star <- delta_1 + as.vector((Z %*% myZcoef)) + sample(myS.tilde, size = length(myS.tilde), replace = TRUE)
  #y.star.matrix[,i] <- y.star
  #plotl(y.star)
  test.stat <- ur.ls.bootstrap(y=y.star , model = model, breaks = breaks, lags = 13, method = "Fixed",pn = 0.1, critval = "bootstrap", print.results = "silent")
  if(breaks == 2){
    return(list(c(unlist(test.stat[[1]]$`t-stat`)), c(unlist(test.stat[[1]]$`First break`)) ,c(unlist(test.stat[[1]]$`Second break`)),  c(as.numeric(test.stat[[1]]$Runtime, units = "mins"))))
  } else if (breaks == 1){
    return(c(c(unlist(test.stat[[1]]$`t-stat`)), c(unlist(test.stat[[1]]$`First break`)) , c(as.numeric(test.stat[[1]]$Runtime, units = "mins"))))
  }
}


if(breaks == 2){
  colnames(boot.y) <- c("t-stat","First break", "Second break", "Runtime")
} else if (breaks == 1){
  colnames(boot.y) <- c("t-stat","First break",  "Runtime")
}

MPM_BOOTS_RESULTS_crash <- test.results.y
MPM_BOOTS_INTERVALS_crash <- boot.y


#Calculate the critical values, based on the distribution of the test statistic
hist(unlist(MPM_BOOTS_INTERVALS[,"t-stat"]), freq = FALSE)
quantile(unlist(MPM_BOOTS_INTERVALS[,"t-stat"]), c(.01, 0.025, .05, .1, .25, .5, .75, .95, .975, .99))
quantile(unlist(MPM_BOOTS_INTERVALS[,"First break"]), c(.01, 0.25, .05, .1, .25, .5, .75, .95, .975, .99))
quantile(unlist(MPM_BOOTS_INTERVALS[,"Second break"]), c(.01, 0.25, .05, .1, .25, .5, .75, .95, .975, .99))

hist(unlist(MPM_BOOTS_INTERVALS[,"Second break"]), freq = F)
hist(unlist(MPM_BOOTS_INTERVALS[,"First break"]), freq = F)
