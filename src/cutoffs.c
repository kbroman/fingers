/**********************************************************************
 * 
 * cutoffs.c
 *
 * copyright (c) 2001, Karl W Broman
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions to calculate cutoffs for the hierarchical clustering
 * in the Fingers program.
 *
 * Contains: cutoff, cutoff_s, cutoff_lr, 
 *           cutoff_llr, cutoff_llr_s, cutoff_llr_meansib, cutoff_llr_lr
 *
 * ARGUMENTS: nmar     = number of markers
 *            freq     = vector of band allele frequencies
 *            quantile = quantile to use as cutoff
 *            cutoff   = (the return value)
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include "llrdist.h"
#include "cutoffs.h"

/* cutoff: (prop'n mismatches)
   use a quantile from the distr'n of dist between unrelateds */
void cutoff(int *nmar, double *freq, double *quantile, 
	    double *cutoff)
{
  int i;
  double p, q, mean, sd;

  mean=sd=0.0;
  for(i=0; i<*nmar; i++) {
    q = 1.0 - freq[i];
    p = 2.0*q*q*(1.0-q*q); /* prob of mismatch */
    mean += p;
    sd += p*(1.0-p);
  }
  mean /= (double)(*nmar);
  sd = sqrt(sd)/(double)(*nmar);
  /* qnorm comes from R : last 2 args are "lower tail" and "give log" */
  *cutoff = qnorm(*quantile, mean, sd, 1, 0);
}


/* cutoff_s: (prop'n mismatches)
   use a quantile from the distr'n of dist between sibings */
void cutoff_s(int *nmar, double *freq, double *quantile, 
	      double *cutoff)
{
  int i;
  double p, q, mean, sd;

  mean=sd=0.0;
  for(i=0; i<*nmar; i++) {
    p = freq[i]; q = 1.0 - p;
    p = (1.0 - p*p - 2.0*p*p*p*q - 3.5*p*p*q*q - 2.0*p*q*q*q -
	 q*q*q*q);
    mean += p;
    sd += p*(1.0-p);
  }
  mean /= (double)(*nmar);
  sd = sqrt(sd)/(double)(*nmar);
  /* qnorm comes from R : last 2 args are "lower tail" and "give log" */
  *cutoff = qnorm(*quantile, mean, sd, 1, 0);
}


/* cutoff_lr: (prop'n mismatches)
   use value of LR comparing unrelateds to siblings */
void cutoff_lr(int *nmar, double *freq, double *value, 
	       double *cutoff)
{
  int i;
  double p, q, ps, pu, meanu, sdu, means, sds;

  meanu=sdu=means=sds=0.0;
  for(i=0; i<*nmar; i++) {
    p=freq[i]; q = 1.0 - p;
    pu = 2.0*q*q*(1.0-q*q); /* prob of mismatch */
    ps = (1.0 - p*p - 2.0*p*p*p*q - 3.5*p*p*q*q - 2.0*p*q*q*q -
	 q*q*q*q);
    meanu += pu; means += ps;
    sdu += pu*(1.0-pu); sds += ps*(1.0-ps);
  }
  meanu /= (double)(*nmar);
  sdu = sqrt(sdu)/(double)(*nmar);
  means /= (double)(*nmar);
  sds = sqrt(sds)/(double)(*nmar);

  /* dnorm comes from R : last 2 args are "lower tail" and "give log" */
  *cutoff = find_lr(meanu,sdu,means,sds,*value);
}


/* cutoff_llr: (Karl's LLR distance)
   uses a quantile from the distr'n of dist between unrelateds */
void cutoff_llr(int *nmar, double *freq, double *quantile, 
		double *cutoff)
{
  int i, j;
  double mean, sd, meantemp, vartemp, d[3], pr[3], p, q;

  mean=sd=0.0;
  for(i=0; i<*nmar; i++) {
    p = freq[i]; q = 1.0 - p;
    calc_llr_dist(d, p);
    pr[0] = q*q*q*q;  /* prob neither has band */
    pr[1] = 2.0*q*q*(1.0-q*q); /* prob of mismatch */
    pr[2] = 1.0-pr[0]-pr[1];
    meantemp = vartemp = 0.0;
    for(j=0; j<3; j++) {
      meantemp += pr[j]*d[j];
      vartemp += pr[j]*d[j]*d[j];
    }
    mean += meantemp;
    sd += (vartemp - meantemp*meantemp);
  }
  mean /= (double)(*nmar);
  sd = sqrt(sd)/(double)(*nmar);
  /* qnorm comes from R : last 2 args are "lower tail" and "give log" */
  *cutoff = qnorm(*quantile, mean, sd, 1, 0);
}


/* cutoff_llr_s: (Karl's LLR distance)
   uses a quantile from the distr'n of dist between siblings */
void cutoff_llr_s(int *nmar, double *freq, double *quantile, 
		  double *cutoff)
{
  int i, j;
  double mean, sd, meantemp, vartemp, d[3], pr[3], p, q;

  mean=sd=0.0;
  for(i=0; i<*nmar; i++) {
    p = freq[i]; q = 1.0 - p;
    calc_llr_dist(d, p);
    pr[0] = 0.25*p*p*q*q+p*q*q*q+q*q*q*q;  /* prob neither has band */
    pr[2] = p*p+2.0*p*p*p*q+3.25*p*p*q*q+p*q*q*q; /* prob both have band */
    pr[1] = 1.0-pr[0]-pr[2]; /* prob of mismatch */
    meantemp = vartemp = 0.0;
    for(j=0; j<3; j++) {
      meantemp += pr[j]*d[j];
      vartemp += pr[j]*d[j]*d[j];
    }
    mean += meantemp;
    sd += (vartemp - meantemp*meantemp);
  }
  mean /= (double)(*nmar);
  sd = sqrt(sd)/(double)(*nmar);
  /* qnorm comes from R : last 2 args are "lower tail" and "give log" */
  *cutoff = qnorm(*quantile, mean, sd, 1, 0);
}


/* cutoff_llr_meansib: (Karl's LLR distance)
   uses the mean distance between siblings */
void cutoff_llr_meansib(int *nmar, double *freq, double *cutoff)
{
  int i, j;
  double mean, meantemp, d[3], pr[3], p, q;

  mean=0.0;
  for(i=0; i<*nmar; i++) {
    p = freq[i]; q = 1.0 - p;
    calc_llr_dist(d, p);
    pr[0] = 0.25*p*p*q*q+p*q*q*q+q*q*q*q;  /* prob neither has band */
    pr[2] = p*p+2.0*p*p*p*q+3.25*p*p*q*q+p*q*q*q; /* prob both have band */
    pr[1] = 1.0-pr[0]-pr[2]; /* prob of mismatch */
    meantemp = 0.0;
    for(j=0; j<3; j++) 
      meantemp += pr[j]*d[j];
    mean += meantemp;
  }
  mean /= (double)(*nmar);
  /* qnorm comes from R : last 2 args are "lower tail" and "give log" */
  *cutoff = mean;
}


/* cutoff_llr_lr: (Karl's LLR distance)
   use value of LR comparing unrelateds to siblings */
void cutoff_llr_lr(int *nmar, double *freq, double *value,
		   double *cutoff)
{
  int i, j;
  double meanu, sdu, meanutemp, varutemp;
  double means, sds, meanstemp, varstemp;
  double d[3], pru[3], prs[3], p, q;

  meanu=sdu=means=sds=0.0;
  for(i=0; i<*nmar; i++) {
    p = freq[i]; q = 1.0 - p;
    calc_llr_dist(d, p);
    pru[0] = q*q*q*q;  /* prob neither has band */
    pru[1] = 2.0*q*q*(1.0-q*q); /* prob of mismatch */
    pru[2] = 1.0-pru[0]-pru[1];

    prs[0] = 0.25*p*p*q*q+p*q*q*q+q*q*q*q;  /* prob neither has band */
    prs[2] = p*p+2.0*p*p*p*q+3.25*p*p*q*q+p*q*q*q; /* prob both have band */
    prs[1] = 1.0-prs[0]-prs[2]; /* prob of mismatch */

    meanutemp = varutemp = meanstemp = varstemp = 0.0;
    for(j=0; j<3; j++) {
      meanutemp += pru[j]*d[j];
      varutemp += pru[j]*d[j]*d[j];
      meanstemp += prs[j]*d[j];
      varstemp += prs[j]*d[j]*d[j];
    }
    meanu += meanutemp;
    sdu += (varutemp - meanutemp*meanutemp);
    means += meanstemp;
    sds += (varstemp - meanstemp*meanstemp);
  }
  means /= (double)(*nmar);
  meanu /= (double)(*nmar);
  sds = sqrt(sds)/(double)(*nmar);
  sdu = sqrt(sdu)/(double)(*nmar);

  *cutoff = find_lr(meanu,sdu,means,sds,*value);
}


#define TOL    1e-6
#define MAXIT  1000

/* find point x at which dnorm(x,m1,s1)/dnorm(x,m2,s2) = value */
/* using a simple bisection method */
double find_lr(double m1, double s1, double m2, double s2, double value)
{
  int i;
  double lo=m2, hi=m1, mid=(lo+hi)/2.0;
  double flo, fhi, fmid;

  value = log(value);

  flo = dnorm(lo,m1,s1,1)-dnorm(lo,m2,s2,1);
  fhi = dnorm(hi,m1,s1,1)-dnorm(hi,m2,s2,1);
  
  for(i=0; i<MAXIT; i++) {
    fmid = dnorm(mid,m1,s1,1)-dnorm(mid,m2,s2,1);

    if(fabs(fmid-value) < TOL) return(mid);

    if(fmid > value) { 
      hi=mid; fhi=fmid;
    }
    else {
      lo=mid; flo=fmid;
    }

    mid=(lo+hi)/2.0;
  }

  warning("Didn't converge\n");
  return(mid);
}
      


/* end of cutoffs.c */
