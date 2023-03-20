/*  Copyright 2022 Annamaria Guolo (University of Padova) 
*  Permission to use, copy, modify and distribute this software and
*   its documentation, for any purpose and without fee, is hereby granted,
*   provided that:
##  1) this copyright notice appears in all copies and
*   2) the source is acknowledged through a citation to the paper
*      Guolo A. (2022). Approximate likelihood and pseudo-likelihood inference in meta-analysis of diagnostic accuracy studies accounting for disease prevalence and study's design. Submitted.
*   The Authors make no representation about the suitability of this software
*   for any purpose.  It is provided "as is", without express or implied warranty

*  Parameter vector
*  theta= c(1 = mu.eta, 2 = mu.xi, 3 = mu.gamma,
*           4 = var.eta, 5 = var.xi, 6 = var.gamma,
*           7 = rho.etaxi, 8 = rho.etagamma, 9 = rho.xigamma)
*/
#include <R.h>
#include <Rmath.h>
#include <R_ext/RS.h>
#include <R_ext/Applic.h>


void int_normal(const double *theta,
		const int *ncohort,
		const int *ncc,
		const double *n1cohort, 
		const double *x1cohort,
		const double *n0cohort,
		const double *x0cohort,		
		const double *ntotcohort,
		const double *n1cc, 
		const double *x1cc,
		const double *n0cc,
		const double *x0cc,		
		const double *nodes1,
		const double *nodes2,
		const double *nodes3,
		const double *weights1,
		const double *weights2,
		const double *weights3,
		const int *num,
		double *lik)
{
  int i,l,j,k;
  double ftmp, px0, px1, pp, px0cc, px1cc, covetaxi, covetapi, covxipi, mubar1, mubar2, var1, var2, sigma12, rho, zeta;
  
  covetaxi = theta[6]*sqrt(theta[3]*theta[4]);
  covetapi = theta[7]*sqrt(theta[3]*theta[5]);
  covxipi = theta[8]*sqrt(theta[4]*theta[5]);
  
  for(i=0; i<(*ncohort); i++){
    ftmp = 0.0;
    for(l=0; l<(*num); l++){
      for(j=0; j<(*num); j++){
	      for(k=0; k<(*num); k++){
		
		px0 = exp(nodes2[l])/(1+exp(nodes2[l]));
		px1 = exp(nodes1[j])/(1+exp(nodes1[j]));
		pp = exp(nodes3[k])/(1+exp(nodes3[k]));
	       
		mubar1 = theta[0] + covetapi*(nodes3[k]-theta[2])/theta[5];
		mubar2 = theta[1] + covxipi*(nodes3[k]-theta[2])/theta[5];
		
		var1 = theta[3] - covetapi*covetapi/theta[5];
		var2 = theta[4] - covxipi*covxipi/theta[5];
		sigma12 = covetaxi - covetapi*covxipi/theta[5];
		
		rho = sigma12/sqrt(var1*var2);
		zeta = (nodes1[j]-mubar1)*(nodes1[j]-mubar1)/var1 +
		  (nodes2[l]-mubar2)*(nodes2[l]-mubar2)/var2 -
		  2*rho*(nodes1[j]-mubar1)*(nodes2[l]-mubar2)/sqrt(var1*var2);
		
		ftmp += dbinom(x0cohort[i], n0cohort[i], px0, 0) *
		  dbinom(x1cohort[i], n1cohort[i], px1, 0) *
		  dbinom(x1cohort[i], ntotcohort[i], pp, 0) *
		  exp(-zeta/(2*(1-rho*rho)))/(2*M_PI*sqrt(var1*var2*(1-rho*rho))) *
		  dnorm(nodes3[k], theta[2], sqrt(theta[5]), 0) *
		  exp(nodes1[j]*nodes1[j]) * exp(nodes2[l]*nodes2[l]) * exp(nodes3[k]*nodes3[k])*
		  weights1[j] * weights2[l] * weights3[k];
	      }
      }
    }
    (*lik) += log(ftmp);
  }
  
  for(i=0; i<(*ncc); i++){
    ftmp = 0.0;
    for(l=0; l<(*num); l++){
      for(j=0; j<(*num); j++){
	
	px0cc = exp(nodes2[l])/(1+exp(nodes2[l]));
	px1cc = exp(nodes1[j])/(1+exp(nodes1[j]));
	ftmp += dbinom(x0cc[i], n0cc[i], px0cc, 0) * dbinom(x1cc[i], n1cc[i], px1cc, 0) *
	  dnorm(nodes2[l], theta[1], sqrt(theta[4]), 0) *
	  dnorm(nodes1[j], theta[0]+theta[6]*sqrt(theta[3]/theta[4])*(nodes2[l]-theta[1]), sqrt(theta[3]*(1-theta[6]*theta[6])), 0) *
	  exp(nodes1[j]*nodes1[j]) * exp(nodes2[l]*nodes2[l]) * weights1[j] * weights2[l];
      }
    }
    (*lik) += log(ftmp);
  }
}




void int_pliknormal(const double *theta,
		const int *ncohort,
		const int *ncc,
		const double *n1cohort, 
		const double *x1cohort,
		const double *n0cohort,
		const double *x0cohort,		
		const double *ntotcohort,
		const double *n1cc, 
		const double *x1cc,
		const double *n0cc,
		const double *x0cc,		
		const double *nodes1,
		const double *nodes2,
		const double *nodes3,
		const double *weights1,
		const double *weights2,
		const double *weights3,
		const int *num,
		double *clik)
{
  int i,l,j,k;
  double ftmp, px0, px1, pp, px0cc, px1cc, covetaxi, covetapi, covxipi, mubar1, mubar2, var1, var2, sigma12, rho, zeta;
  
  covetaxi = 0;
  covetapi = 0;
  covxipi = 0;
  
  for(i=0; i<(*ncohort); i++){
    ftmp = 0.0;
    for(l=0; l<(*num); l++){
      for(j=0; j<(*num); j++){
	      for(k=0; k<(*num); k++){
		
		px0 = exp(nodes2[l])/(1+exp(nodes2[l]));
		px1 = exp(nodes1[j])/(1+exp(nodes1[j]));
		pp = exp(nodes3[k])/(1+exp(nodes3[k]));
	       
		mubar1 = theta[0] + covetapi*(nodes3[k]-theta[2])/theta[5];
		mubar2 = theta[1] + covxipi*(nodes3[k]-theta[2])/theta[5];
		
		var1 = theta[3] - covetapi*covetapi/theta[5];
		var2 = theta[4] - covxipi*covxipi/theta[5];
		sigma12 = covetaxi - covetapi*covxipi/theta[5];
		
		rho = sigma12/sqrt(var1*var2);
		zeta = (nodes1[j]-mubar1)*(nodes1[j]-mubar1)/var1 +
		  (nodes2[l]-mubar2)*(nodes2[l]-mubar2)/var2 -
		  2*rho*(nodes1[j]-mubar1)*(nodes2[l]-mubar2)/sqrt(var1*var2);
		
		ftmp += dbinom(x0cohort[i], n0cohort[i], px0, 0) * dbinom(x1cohort[i], n1cohort[i], px1, 0) * dbinom(x1cohort[i], ntotcohort[i], pp, 0) *
		  exp(-zeta/(2*(1-rho*rho)))/(2*M_PI*sqrt(var1*var2*(1-rho*rho))) *
		  dnorm(nodes3[k], theta[2], sqrt(theta[5]), 0) *
		  exp(nodes1[j]*nodes1[j]) * exp(nodes2[l]*nodes2[l]) * exp(nodes3[k]*nodes3[k])*
		  weights1[j] * weights2[l] * weights3[k];
	      }
      }
    }
    (*clik) += log(ftmp);
  }
  
  for(i=0; i<(*ncc); i++){
    ftmp = 0.0;
    for(l=0; l<(*num); l++){
      for(j=0; j<(*num); j++){
	
	px0cc = exp(nodes2[l])/(1+exp(nodes2[l]));
	px1cc = exp(nodes1[j])/(1+exp(nodes1[j]));
	ftmp += dbinom(x0cc[i], n0cc[i], px0cc, 0) * dbinom(x1cc[i], n1cc[i], px1cc, 0) *
	  dnorm(nodes2[l], theta[1], sqrt(theta[4]), 0) *
	  dnorm(nodes1[j], theta[0], sqrt(theta[3]), 0) *
	  exp(nodes1[j]*nodes1[j]) * exp(nodes2[l]*nodes2[l]) * weights1[j] * weights2[l];
      }
    }
    (*clik) += log(ftmp);
  }
}



void int_pliknormal2(const double *theta,
		const int *ncohort,
		const int *ncc,
		const double *n1cohort, 
		const double *x1cohort,
		const double *n0cohort,
		const double *x0cohort,		
		const double *ntotcohort,
		const double *n1cc, 
		const double *x1cc,
		const double *n0cc,
		const double *x0cc,		
		const double *nodes1,
		const double *nodes2,
		const double *nodes3,
		const double *weights1,
		const double *weights2,
		const double *weights3,
		const int *num,
		double *clikdue)
{
  int i,l,j,k;
  double ftmp, px0, px1, pp, px0cc, px1cc, covetaxi, covetapi, covxipi, mubar1, mubar2, var1, var2, sigma12, rho, zeta;
  
  covetaxi = theta[6]*sqrt(theta[3]*theta[4]);
  covetapi = 0;
  covxipi = 0;
  
  for(i=0; i<(*ncohort); i++){
    ftmp = 0.0;
    for(l=0; l<(*num); l++){
      for(j=0; j<(*num); j++){
	      for(k=0; k<(*num); k++){
		
		px0 = exp(nodes2[l])/(1+exp(nodes2[l]));
		px1 = exp(nodes1[j])/(1+exp(nodes1[j]));
		pp = exp(nodes3[k])/(1+exp(nodes3[k]));
	       
		mubar1 = theta[0] + covetapi*(nodes3[k]-theta[2])/theta[5];
		mubar2 = theta[1] + covxipi*(nodes3[k]-theta[2])/theta[5];
		
		var1 = theta[3] - covetapi*covetapi/theta[5];
		var2 = theta[4] - covxipi*covxipi/theta[5];
		sigma12 = covetaxi - covetapi*covxipi/theta[5];
		
		rho = sigma12/sqrt(var1*var2);
		zeta = (nodes1[j]-mubar1)*(nodes1[j]-mubar1)/var1 +
		  (nodes2[l]-mubar2)*(nodes2[l]-mubar2)/var2 -
		  2*rho*(nodes1[j]-mubar1)*(nodes2[l]-mubar2)/sqrt(var1*var2);
		
		ftmp += dbinom(x0cohort[i], n0cohort[i], px0, 0) * dbinom(x1cohort[i], n1cohort[i], px1, 0) * dbinom(x1cohort[i], ntotcohort[i], pp, 0) *
		  exp(-zeta/(2*(1-rho*rho)))/(2*M_PI*sqrt(var1*var2*(1-rho*rho))) *
		  dnorm(nodes3[k], theta[2], sqrt(theta[5]), 0) *
		  exp(nodes1[j]*nodes1[j]) * exp(nodes2[l]*nodes2[l]) * exp(nodes3[k]*nodes3[k])*
		  weights1[j] * weights2[l] * weights3[k];
	      }
      }
    }
    (*clikdue) += log(ftmp);
  }
  
  for(i=0; i<(*ncc); i++){
    ftmp = 0.0;
    for(l=0; l<(*num); l++){
      for(j=0; j<(*num); j++){
	
	px0cc = exp(nodes2[l])/(1+exp(nodes2[l]));
	px1cc = exp(nodes1[j])/(1+exp(nodes1[j]));
	ftmp += dbinom(x0cc[i], n0cc[i], px0cc, 0) * dbinom(x1cc[i], n1cc[i], px1cc, 0) *
	  dnorm(nodes2[l], theta[1], sqrt(theta[4]), 0) *
	  dnorm(nodes1[j], theta[0], sqrt(theta[3]*(1-theta[6]*theta[6])), 0) *
	  exp(nodes1[j]*nodes1[j]) * exp(nodes2[l]*nodes2[l]) * weights1[j] * weights2[l];
      }
    }
    (*clikdue) += log(ftmp);
  }
}

