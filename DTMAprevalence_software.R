##  Copyright 2023 Annamaria Guolo (University of Padova) 
##  Permission to use, copy, modify and distribute this software and
##  its documentation, for any purpose and without fee, is hereby granted,
##  provided that:
##  1) this copyright notice appears in all copies and
##  2) the source is acknowledged through a citation to the paper
##     Guolo A. (2022). Approximate likelihood and pseudo-likelihood inference in meta-analysis of diagnostic accuracy studies accounting for disease prevalence and study's design. Submitted.
##  The Authors make no representation about the suitability of this software
##  for any purpose.  It is provided "as is", without express or implied warranty

## Parameter vector
## theta= c(1 = mu.eta, 2 = mu.xi, 3 = mu.gamma,
##          4 = var.eta, 5 = var.xi, 6 = var.gamma,
##          7 = rho.etaxi, 8 = rho.etagamma, 9 = rho.xigamma)

library(mvtnorm)
library(numDeriv)
library(far)
library(statmod)
library(nlme)
## nodes and weights for 21 points Gauss-Hermite quadrature
objGH <- gauss.quad(21, 'hermite') 

## this.data: dataset with information about
## eta.obs, xi.obs, gamma.obs, var.eta, var.xi, var.gamma, design (cc: case-control, cohort)

## data.bin: dataset with information about
## TP, P, FP, FN, TN, N, design (cc: case-control, cohort)
## if you want, order the data so that cc and cohort studies are consecutive
## sequence <- order(this.data$design, decreasing=TRUE)
## this.data <- this.data[sequence,]

## ========================================
## APPROXIMATE MODEL: LIKELIHOOD FUNCTION
## maximumm likelihood estimation
## ========================================
lik.fn <- function(theta, this.data){
    n <- NROW(this.data)
    data.cc <- this.data[this.data$design=='cc',]
    data.cohort <- this.data[this.data$design=='cohort',]
    n.cohort <- nrow(data.cohort)
    n.cc <- n - n.cohort
    mu.eta <- theta[1]
    mu.xi <- theta[2]
    mu.gamma <- theta[3]
    cond <- (theta[7] > 1 | theta[7] < -1 | theta[8] > 1 | theta[8] < -1 | theta[9] > 1 | theta[9] < -1 )
    if(theta[4] <= 0 | theta[5] <= 0 | theta[6] <= 0 | cond==TRUE)
        return(NA)
    sigma2.eta <- theta[4]
    sigma2.xi <- theta[5]
    sigma2.gamma <- theta[6]
    rho.etaxi <- theta[7]
    rho.etagamma <- theta[8]
    rho.xigamma <- theta[9]
    f <- 0
    if(n.cohort >= 1){
        for(i in 1:n.cohort){
            eta.i <- data.cohort[i,'eta.obs']
            xi.i <- data.cohort[i,'xi.obs']
            gamma.i <- data.cohort[i,'gamma.obs']
            var.eta <- sigma2.eta + data.cohort[i,'var.eta'] 
            var.xi <- sigma2.xi + data.cohort[i,'var.xi']    
            var.gamma <- sigma2.gamma + data.cohort[i,'var.gamma']    
            cov.etaxi <- rho.etaxi*sqrt(var.eta*var.xi)
            cov.etagamma <- rho.etagamma*sqrt(var.eta*var.gamma)
            cov.xigamma <- rho.xigamma*sqrt(var.xi*var.gamma)
            var.matrix <- matrix(c(var.eta, cov.etaxi, cov.etagamma,
                                   cov.etaxi, var.xi, cov.xigamma,
                                   cov.etagamma, cov.xigamma, var.gamma),
                                 ncol=3, nrow=3)
            f <- f + dmvnorm(c(data.cohort[i,1], data.cohort[i,2], data.cohort[i,3]),
                             mean=c(mu.eta, mu.xi, mu.gamma),
                             sigma=var.matrix,
                             log=TRUE)
        }
    }
    if(n.cc >= 1){
        for(i in 1:n.cc){
            eta.i <- data.cc[i,'eta.obs']
            xi.i <- data.cc[i,'xi.obs']
            var.eta <- sigma2.eta + data.cc[i,'var.eta'] 
            var.xi <- sigma2.xi + data.cc[i,'var.xi']    
            cov.etaxi <- rho.etaxi*sqrt(var.eta*var.xi)
            var.matrix <- matrix(c(var.eta, cov.etaxi, 
                                   cov.etaxi, var.xi),
                                 ncol=2, nrow=2)
            f <- f + dmvnorm(c(data.cc[i,1], data.cc[i,2]),
                             mean=c(mu.eta, mu.xi),
                             sigma=var.matrix,
                             log=TRUE)
        }    
    }
        return(f)
}

## maximize the likelihood function
lik <- function(this.data, theta.start=theta.start, method='Nelder-Mead', gr=NULL){
    ris <- optim(theta.start, lik.fn, this.data=this.data, control=list(maxit=5000, fnscale=-1), hessian=TRUE, method=method)
    ## MLE
    theta.hat <- ris$par
    theta.var <- solve(-ris$hessian)
    ## compute the sandwich s.e.
    theta.var.sandwich <- sand.approximate(theta=theta.hat, this.data=this.data)
    return(list(theta.hat=theta.hat,
                theta.se.hessian=sqrt(diag((theta.var))),
                convergence=ris$convergence,
                var.matrix.hessian=theta.var,
                theta.se.sandwich=sqrt(diag((theta.var.sandwich))),
                var.matrix.sandwich=theta.var.sandwich))
}

## ==========================================
## APPROXIMATE MODEL: LIKELIHOOD FUNCTION
## restricted maximumm likelihood estimation
## ==========================================
## maximize w.r.t. to variance/covariance components
reml.var.to.optim <- function(var.par, this.data, fixed.par){
    ## fixed-effects parameters are fixed at fixed.par
    n <- NROW(this.data)
    data.cc <- this.data[this.data$design=='cc',]
    data.cohort <- this.data[this.data$design=='cohort',]
    n.cohort <- nrow(data.cohort)
    n.cc <- n - n.cohort  
    mu.eta <- fixed.par[1]
    mu.xi <- fixed.par[2]
    mu.gamma <- fixed.par[3]
    sigma2.eta <- var.par[1]
    sigma2.xi <- var.par[2]
    sigma2.gamma <- var.par[3]
    rho.etaxi <- var.par[4]
    rho.etagamma <- var.par[5]
    rho.xigamma <- var.par[6]
    cond <- (rho.etaxi > 1 | rho.etaxi < -1 | rho.etagamma > 1 | rho.etagamma < -1 | rho.xigamma > 1 | rho.xigamma < -1)
    if( sigma2.eta <= 0 | sigma2.xi <= 0 | sigma2.gamma <= 0 | cond==TRUE)
        return(NA)
    f <- 0
    for(i in 1:n.cc){
        eta.i <- data.cc[i,1]
        xi.i <- data.cc[i,2]
        r.i <- c(eta.i, xi.i) - c(mu.eta, mu.xi)
        var.eta <- sigma2.eta + data.cc[i,4] 
        var.xi <- sigma2.xi + data.cc[i,5]
        cov.etaxi <- rho.etaxi*sqrt(var.eta*var.xi)
        var.matrix <- matrix(c(var.eta, cov.etaxi,
                               cov.etaxi, var.xi), ncol=2, nrow=2)
        f <- f + dmvnorm(r.i,
                         mean=c(0, 0),
                         sigma=var.matrix,
                         log=TRUE)
    }
    for(i in 1:n.cohort){
        eta.i <- data.cohort[i,1]
        xi.i <- data.cohort[i,2]
        gamma.i <- data.cohort[i,3]
        r.i <- c(eta.i, xi.i, gamma.i) - c(mu.eta, mu.xi, mu.gamma)
        var.eta <- sigma2.eta + data.cohort[i,4] 
        var.xi <- sigma2.xi + data.cohort[i,5]
        var.gamma <- sigma2.gamma + data.cohort[i,6]
        cov.etaxi <- rho.etaxi*sqrt(var.eta*var.xi)
        cov.etagamma <- rho.etagamma*sqrt(var.eta*var.gamma)
        cov.xigamma <- rho.xigamma*sqrt(var.xi*var.gamma)
        var.matrix <- matrix(c(var.eta, cov.etaxi, cov.etagamma,
                               cov.etaxi, var.xi, cov.xigamma,
                               cov.etagamma, cov.xigamma, var.gamma), ncol=3, nrow=3)
        f <- f + dmvnorm(r.i,
                         mean=c(0, 0, 0),
                         sigma=var.matrix,
                         log=TRUE)
    }
    return(f)
}

## maximize w.r.t. to variance/covariance components
reml.fixed.to.optim <- function(fixed.par, this.data, var.par){
    ## variance/covariance parameters are fixed at var.par
    n <- NROW(this.data)
    data.cc <- this.data[data$design=='cc',]
    data.cohort <- this.data[data$design=='cohort',]
    n.cohort <- nrow(data.cohort)
    n.cc <- n - n.cohort  
    mu.eta <- fixed.par[1]
    mu.xi <- fixed.par[2]
    mu.gamma <- fixed.par[3]
    sigma2.eta <- var.par[1]
    sigma2.xi <- var.par[2]
    sigma2.gamma <- var.par[3]
    rho.etaxi <- var.par[4]
    rho.etagamma <- var.par[5]
    rho.xigamma <- var.par[6]   
    cond <- (rho.etaxi > 1 | rho.etaxi < -1 | rho.etagamma > 1 | rho.etagamma < -1 | rho.xigamma > 1 | rho.xigamma < -1)
    if( sigma2.eta <= 0 | sigma2.xi <= 0 | sigma2.gamma <= 0 | cond==TRUE)
        return(NA)
    f <- 0
    for(i in 1:n.cohort){
        eta.i <- data.cohort[i,1]
        xi.i <- data.cohort[i,2]
        gamma.i <- data.cohort[i,3]
        r.i <- c(eta.i, xi.i, gamma.i) - c(mu.eta, mu.xi, mu.gamma)
        var.eta <- sigma2.eta + data.cohort[i,4] 
        var.xi <- sigma2.xi + data.cohort[i,5]
        var.gamma <- sigma2.gamma + data.cohort[i,6]
        cov.etaxi <- rho.etaxi*sqrt(var.eta*var.xi)
        cov.etagamma <- rho.etagamma*sqrt(var.eta*var.gamma)
        cov.xigamma <- rho.xigamma*sqrt(var.xi*var.gamma)
        var.matrix <- matrix(c(var.eta, cov.etaxi, cov.etagamma,
                               cov.etaxi, var.xi, cov.xigamma,
                               cov.etagamma, cov.xigamma, var.gamma), ncol=3, nrow=3)
        f <- f + dmvnorm(r.i,
                         mean=c(0, 0, 0),
                         sigma=var.matrix,
                         log=TRUE)
    }
    for(i in 1:n.cc){
        eta.i <- data.cc[i,1]
        xi.i <- data.cc[i,2]
        r.i <- c(eta.i, xi.i) - c(mu.eta, mu.xi)
        var.eta <- sigma2.eta + data.cc[i,4] 
        var.xi <- sigma2.xi + data.cc[i,5]
        cov.etaxi <- rho.etaxi*sqrt(var.eta*var.xi)
        var.matrix <- matrix(c(var.eta, cov.etaxi,  
                               cov.etaxi, var.xi), ncol=2, nrow=2)
        f <- f + dmvnorm(r.i,
                         mean=c(0, 0),
                         sigma=var.matrix,
                         log=TRUE)
    }
  return(f)
}

## restricted maximum likelihood
lik.reml <- function(this.data, eps=0.0001, method='Nelder-Mead'){
  fixed.start <- c(mean(this.data[,1]), mean(this.data[,2]), mean(this.data[,3], na.rm=TRUE))
  var.start <- c(var(this.data[,1]), var(this.data[,2]), var(this.data[,3], na.rm=TRUE), cor(this.data[,1], this.data[,2]), cor(this.data[this.data$design=='cohort',1], this.data[this.data$design=='cohort',3]), cor(this.data[this.data$design=='cohort',2], this.data[this.data$design=='cohort',3]))
  var.new <- optim(var.start, reml.var.to.optim, this.data=this.data, fixed.par=fixed.start, control=list(fnscale=-1), method=method)$par
  ## iterate until convergence or until the number of iterations is less than 100
  conv <- step <- 0
  while(conv < 1 & step < 100){
    var.start <- var.new
    fixed.new <- optim(fixed.start, reml.fixed.to.optim, this.data=this.data, var.par=var.start, control=list(fnscale=-1), method=method)$par
    var.new <- optim(var.start, reml.var.to.optim, this.data=this.data, fixed.par=fixed.new, control=list(fnscale=-1), method=method)$par
    difference <- abs(var.new-var.start)
    if(any(difference[1]<=eps))
        conv <- 1
    else
        step <- step+1
  }
  theta.hat <- c(fixed.new, var.new)
  theta.var <- try(solve(-fdHess(theta.hat, fn.lik, this.data=this.data)$Hessian), silent=TRUE)
  if(any(is.na(theta.var)))
      theta.var <- matrix(NA, ncol=length(theta.hat), nrow=length(theta.hat))
  ## compute the sandwich s.e.
  theta.var.sandwich <- sand.approximate(theta=theta.hat, this.data=this.data)
  return(list(theta.hat=theta.hat,
              theta.se.hessian=sqrt(diag(theta.var)),
              var.matrix.hessian = theta.var,
              theta.se.sandwich = sqrt(diag(theta.var.sandwich)),
              var.matrix.sandwich = theta.var.sandwich))
}


## ==========================================
## EXACT TRIVARIATE MODEL: LIKELIHOOD FUNCTION
## maximum likelihood estimation
## ==========================================
#### implemented with some .C code
## data.bin includes TP, TN, FP, FN, P, N
## objGH includes nodes and weights
## for numerical integration
## Fisher's transformation of the correlations
lik.exact.fn <- function(theta.fisher, data.bin, objGH){
    n <- NROW(data.bin)
    theta <- theta.fisher
    theta[7:9] <- tanh(theta.fisher[7:9])    
    data.cc <- data.bin[data.bin$design=='cc',]
    data.cohort <- data.bin[data.bin$design=='cohort',]
    n.cohort <- nrow(data.cohort)
    n.cc <- n - n.cohort
    cond <- (theta[7] > 1 | theta[7] < -1 | theta[8] > 1 | theta[8] < -1 | theta[9] > 1 | theta[9] < -1)
    if(theta[4] <= 0 | theta[5] <= 0 | theta[6] <= 0 | cond==TRUE)
        return(NA)
    nodes1 <- nodes2 <- nodes3 <- objGH$nodes
    weights1 <- weights2 <- weights3 <- objGH$weights
    num <- length(nodes1)
    ## call to .C code
    f <- .C("int_normal",
            as.double(theta),
            as.integer(n.cohort),
            as.integer(n.cc),            
            as.double(data.cohort[,'P']),
            as.double(data.cohort[,'TP']),
            as.double(data.cohort[,'N']),
            as.double(data.cohort[,'TN']),
            as.double(data.cohort[,'P'] + data.cohort[,'N']),            
            as.double(data.cc[,'P']),
            as.double(data.cc[,'TP']),
            as.double(data.cc[,'N']),
            as.double(data.cc[,'TN']),            
            as.double(nodes1),
            as.double(nodes2),
            as.double(nodes3),
            as.double(weights1),
            as.double(weights2),
            as.double(weights3),
            as.integer(num),
            lik=as.double(0.0))$lik
    return(f)
}

## maximize the likelihood function
lik.exact <- function(theta.fisher, data.bin, objGH, method='Nelder-Mead', hessian=TRUE, maxit=maxit){
    ris <- optim(theta.fisher, lik.exact.fn, data.bin=data.bin, objGH=objGH, control=list(maxit=maxit, fnscale=-1), method=method, hessian=hessian)
    theta.hat <- ris$par
    if(hessian)
        theta.var.hessian <- solve(-ris$hessian)
    else
        theta.var.hessian <- solve(-fdHess(theta.hat, lik.exact.fn, data.bin=data.bin, objGH=objGH)$Hessian)
    ## compute the sandwich s.e.
    theta.var.sandwich <- sand.exact(theta=theta.hat, data.bin=data.bin, objGH=objGH)
    return(list(theta.hat=theta.hat,
                theta.se.hessian=sqrt(diag(theta.var.hessian)),
                theta.se.sandwich=sqrt(diag(theta.var.sandwich)),
                convergence=ris$convergence, ## check the convergence
                var.matrix.hessian=theta.var.hessian,
                var.matrix.sandwich=theta.var.sandwich))
}

## ==================================
## PSEUDO-LIKELIHOOD APPROACH
## all the correlations equal to zero
## ==================================
#### implemented with some .C code
## data.bin includes TP, TN, FP, FN, P, N
## objGH includes nodes and weights
## for numerical integration
plik.fn <- function(theta, data.bin, objGH){
    n <- NROW(data.bin)
    data.cc <- data.bin[data.bin$design=='cc',]
    data.cohort <- data.bin[data.bin$design=='cohort',]
    n.cohort <- nrow(data.cohort)
    n.cc <- n - n.cohort
    if( theta[4] <= 0 | theta[5] <= 0 | theta[6] <= 0)
        return(NA)
    nodes1 <- nodes2 <- nodes3 <- objGH$nodes
    weights1 <- weights2 <- weights3 <- objGH$weights
    num <- length(nodes1)
    f <- .C("int_pliknormal",
            as.double(theta),
            as.integer(n.cohort),
            as.integer(n.cc),           
            as.double(data.cohort[,'P']),
            as.double(data.cohort[,'TP']),
            as.double(data.cohort[,'N']),
            as.double(data.cohort[,'TN']), 
            as.double(data.cohort[,'P'] + data.cohort[,'N']),            
            as.double(data.cc[,'P']),
            as.double(data.cc[,'TP']),
            as.double(data.cc[,'N']),
            as.double(data.cc[,'TN']),             
            as.double(nodes1),
            as.double(nodes2),
            as.double(nodes3),
            as.double(weights1),
            as.double(weights2),
            as.double(weights3),
            as.integer(num),
            clik=as.double(0.0))$clik
    return(f)
}
## maximize the pseudo-likelihood function
plik <- function(theta, data.bin, objGH, method='Nelder-Mead', hessian=TRUE, maxit=maxit, two.corr=FALSE){
    theta.start.cl <- theta[1:6]
    ris <- optim(theta.start.cl, plik.fn, data.bin=data.bin, objGH=objGH, control=list(maxit=maxit, fnscale=-1), method=method, hessian=hessian)
    theta.hat <- ris$par
    theta.var.hessian <- solve(-ris$hessian)
    theta.var.sandwich <- sand.pseudo(theta=theta.hat, data.bin=data.bin, objGH=objGH, two.corr=two.corr)
    return(list(theta.hat=theta.hat,
                theta.se.hessian=sqrt(diag(theta.var.hessian)),
                theta.se.sandwich=sqrt(diag(theta.var.sandwich)),
                convergence=ris$convergence, ## check the convergence
                var.matrix.hessian=theta.var.hessian,
                var.matrix.sandwich=theta.var.sandwich))
}

## ==============================================
## PSEUDO-LIKELIHOOD APPROACH
## correlations with the prevalence equal to zero
## ==============================================
#### implemented with some .C code
## data.bin includes TP, TN, FP, FN, P, N
## objGH includes nodes and weights
## for numerical integration
plik.fn2 <- function(theta, data.bin, objGH){
    n <- NROW(data.bin)
    data.cc <- data.bin[data.bin$design=='cc',]
    data.cohort <- data.bin[data.bin$design=='cohort',]
    n.cohort <- nrow(data.cohort)
    n.cc <- n - n.cohort
    if( theta[4] <= 0 | theta[5] <= 0 | theta[6] <= 0 | theta[7] > 1 | theta[7] < -1)
        return(NA)
    nodes1 <- nodes2 <- nodes3 <- objGH$nodes
    weights1 <- weights2 <- weights3 <- objGH$weights
    num <- length(nodes1)
    f <- .C("int_pliknormal2",
            as.double(theta),
            as.integer(n.cohort),
            as.integer(n.cc),           
            as.double(data.cohort[,'P']),
            as.double(data.cohort[,'TP']),
            as.double(data.cohort[,'N']),
            as.double(data.cohort[,'TN']), 
            as.double(data.cohort[,'P'] + data.cohort[,'N']),            
            as.double(data.cc[,'P']),
            as.double(data.cc[,'TP']),
            as.double(data.cc[,'N']),
            as.double(data.cc[,'TN']),          
            as.double(nodes1),
            as.double(nodes2),
            as.double(nodes3),
            as.double(weights1),
            as.double(weights2),
            as.double(weights3),
            as.integer(num),
            clikdue=as.double(0.0))$clikdue
    return(f)
}

## maximize the likelihood function
plik2 <- function(theta, data.bin, objGH, method='Nelder-Mead', hessian=TRUE, maxit=maxit, two.corr=TRUE){
    theta.start.cl <- theta[1:7]
    ris <- optim(theta.start.cl, plik.fn2, data.bin=data.bin, objGH=objGH, control=list(maxit=maxit, fnscale=-1), method=method, hessian=hessian)
    theta.hat <- ris$par
    theta.var.hessian <- solve(-ris$hessian)
    theta.var.sandwich <- sand.pseudo(theta=theta.hat, data.bin=data.bin, objGH=objGH, two.corr=two.corr)
    return(list(theta.hat=theta.hat,
                theta.se.hessian=sqrt(diag(theta.var.hessian)),
                theta.se.sandwich=sqrt(diag(theta.var.sandwich)),
                convergence=ris$convergence, ## check the convergence
                var.matrix.hessian=theta.var.hessian,
                var.matrix.sandwich=theta.var.sandwich))
}

######################################################
## standard error using sandwich for exact likelihood
######################################################
sand.exact <- function(theta, data.bin, objGH){
    a <- fdHess(theta, lik.exact.fn, data.bin=data.bin[1,], objGH=objGH)
    values.gradient <- a$gradient   
    G <- values.gradient%*%t(values.gradient)
    H <- a$Hessian          
    for(i in 2:nrow(data.bin)){
        a <- fdHess(theta, lik.exact.fn, data.bin=data.bin[i,], objGH=objGH)
        values.gradient <- a$gradient  
        G <- G + values.gradient%*%t(values.gradient)
        H <- H + a$Hessian       
    }
    return(solve(H)%*%G%*%solve(H))
}

######################################################
## standard error using sandwich for pseudo-likelihood
######################################################
sand.pseudo <- function(theta, data.bin, objGH, two.corr=FALSE){
    if(!two.corr){
        a <- fdHess(theta, plik.fn, data.bin=data.bin[1,], objGH=objGH)
        values.gradient <- a$gradient 
        G <- values.gradient%*%t(values.gradient)
        H <- a$Hessian  
        for(i in 2:nrow(data.bin)){
            a <- fdHess(theta, plik.fn, data.bin=data.bin[i,], objGH=objGH)
            values.gradient <- a$gradient  
            G <- G + values.gradient%*%t(values.gradient)
            H <- H + a$Hessian       
        }
    }    
    if(two.corr){
        a <- fdHess(theta, plik.fn2, data.bin=data.bin[1,], objGH=objGH)
        values.gradient <- a$gradient 
        G <- values.gradient%*%t(values.gradient)
        H <- a$Hessian         
        for(i in 2:nrow(data.bin)){
            a <- fdHess(theta, plik.fn2, data.bin=data.bin[i,], objGH=objGH)
            values.gradient <- a$gradient 
            G <- G + values.gradient%*%t(values.gradient)
            H <- H + a$Hessian     
        }
    }
    return(solve(H)%*%G%*%solve(H))
}

###########################################################
## standard error using sandwich for approximate likelihood
###########################################################
sand.approximate <- function(theta, this.data){
    a <- fdHess(theta, lik.fn, this.data=this.data[1,])
    values.gradient <- a$gradient
    G <- values.gradient%*%t(values.gradient)
    H <- a$Hessian      
    for(i in 2:nrow(this.data)){
        a <- fdHess(theta, lik.fn, this.data=this.data[i,])
        values.gradient <- a$gradient
        G <- G + values.gradient%*%t(values.gradient)
        H <- H + a$Hessian   
    }
    return(solve(H)%*%G%*%solve(H))
}

