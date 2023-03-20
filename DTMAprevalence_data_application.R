##  Copyright 2023 Annamaria Guolo (University of Padova) 
##  Permission to use, copy, modify and distribute this software and
##  its documentation, for any purpose and without fee, is hereby granted,
##  provided that:
##  1) this copyright notice appears in all copies and
##  2) the source is acknowledged through a citation to the paper
##     Guolo A. (2023). Hierarchical multinomial processing tree models for meta-analysis of diagnostic accuracy studies. arXiv.
##  The Authors make no representation about the suitability of this software
##  for any purpose.  It is provided "as is", without express or implied warranty
###########################################

## APPLICATION

###======================================
### EXAMPLE: beta-D-glucane data
## Karageorgopoulos et al. (CID, 2011)
## data in Hoeyer & Kuss (SIM, 2015)
###======================================

rm(list=ls())
source('DTMAprevalence_software.R')
dyn.load('DTMAprevalence_auxiliary_functions.so')
p <- 9 ## number of paramaters in theta

TP <- c(18,37,59,12,71,15,72,7,14,6,18,12,15,37)
FP <- c(19,2,124,14,36,12,28,3,13,3,15,44,1,26)
FN <- c(4,21,39,11,26,15,4,1,2,5,2,6,4,4)
TN <- c(89,18,635,44,86,101,147,26,11,122,248,139,25,135)
design <- c('cohort','cc','cohort','cohort','cc','cohort','cc','cc','cc', 'cohort','cohort','cohort','cc','cohort')
n <- length(TP)
N <- TN + FP
P <- TP + FN
data.bin <- data.frame(TP, P, FP, N, TN, design)
data.bin <- data.bin[order(data.bin$design, decreasing=TRUE),]

colnames(data.bin) <- c('TP', 'P', 'FP', 'N', 'TN', 'design')  
sp.obs <- se.obs <- prev.obs <- eta.obs <- xi.obs <- gamma.obs <- var.eta <- var.xi <- var.gamma <- rep(NA, n)
for(i in 1:n){
    sp.obs[i] <- 1-data.bin$FP[i]/data.bin$N[i] ## estimated specificity= TN/N
    se.obs[i] <- data.bin$TP[i]/data.bin$P[i]  ## estimated sensitivity= TP/P
    var.eta[i] <- 1/(data.bin$TP[i]) + 1/(data.bin$P[i]-data.bin$TP[i])
    var.xi[i] <- 1/(data.bin$TN[i]) + 1/(data.bin$N[i]-data.bin$TN[i])
    ## apply continuity correction if and where necessary
    if(sp.obs[i]==1 | sp.obs[i]==0 | is.nan(sp.obs[i])){
        sp.obs[i] <- 1 - (data.bin$FP[i]+0.5)/(data.bin$N[i]+1)
        var.xi[i] <- 1/(data.bin$TN[i]+0.5) + 1/(data.bin$N[i]-data.bin$TN[i]+0.5)
    }
    if(se.obs[i]==1 | se.obs[i]==0 | is.nan(se.obs[i])){
        se.obs[i] <- (data.bin$TP[i]+0.5)/(data.bin$P[i]+1)
        var.eta[i] <- 1/(data.bin$TP[i]+0.5) + 1/(data.bin$P[i]-data.bin$TP[i]+0.5)
    }
    eta.obs[i] <- log(se.obs[i] / (1-se.obs[i]))
    xi.obs[i] <- log(sp.obs[i] / (1-sp.obs[i] ))   
}

data.bin.cohort <- data.bin[data.bin$design=='cohort',]
n.cohort <- nrow(data.bin.cohort)

for(i in 1:n.cohort){
    prev.obs[i] <- data.bin.cohort$P[i]/(data.bin.cohort$P[i] + data.bin.cohort$N[i])
    var.gamma[i] <- 1/(data.bin.cohort$P[i]) + 1/(data.bin.cohort$N[i])
    if(prev.obs[i] == 0){
        prev.obs[i] <- (dati.tmpbin.cohort$P[i]+0.5)/(dati.tmpbin.cohort$P[i]+dati.bin.cohort$N[i]+1)
        var.gamma[i] <- 1/(dati.tmpbin.cohort$P[i]+0.5) + 1/(dati.bin.cohort$N[i]+0.5)
    }
    gamma.obs[i] <- log(prev.obs[i] / (1-prev.obs[i]))
}

this.data <- data.frame(eta.obs, xi.obs, gamma.obs, var.eta, var.xi, var.gamma, data.bin$design)
colnames(this.data) <- c('eta.obs', 'xi.obs', 'gamma.obs', 'var.eta', 'var.xi', 'var.gamma', 'design')

theta.start <- c(mean(this.data[,1]), mean(this.data[,2]), mean(this.data[,3], na.rm=TRUE),
                 var(this.data[,1]), var(this.data[,2]), var(this.data[,3], na.rm=TRUE),
                 cor(this.data[,1], this.data[,2]), cor(this.data[1:n.cohort,1], this.data[1:n.cohort,3]),
                 cor(this.data[1:n.cohort,2], this.data[1:n.cohort,3]))

## approximate likelihood
ris.lik.approximate <- lik(this.data, theta.start=theta.start)


##pseudo-likelihood 1
ris.plik <- plik(theta.start, data.bin=data.bin, objGH=objGH, maxit=1000)

##pseudo-likelihood 2
ris.plik2 <- plik2(theta.start, data.bin=data.bin, objGH=objGH, maxit=1000)

## exact likelihood using Fisher's transformation
## starting point: maximum approximate likelihood estimate
theta.start.fisher <- ris.lik.approximate$theta.hat 
theta.start.fisher[7:9] <- atanh(theta.start[7:9])
ris.lik.exact.fisher <- lik.exact(theta.start.fisher, data.bin=data.bin, objGH=objGH, hessian=FALSE, maxit=5000)


##======================================
##=========end of the example===========
##======================================
