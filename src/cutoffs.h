/**********************************************************************
 * 
 * cutoffs.h
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

/* cutoff: (prop'n mismatches)
   use a quantile from the distr'n of dist between unrelateds */
void cutoff(int *nmar, double *freq, double *quantile, 
	    double *cutoff);

/* cutoff_s: (prop'n mismatches)
   use a quantile from the distr'n of dist between sibings */
void cutoff_s(int *nmar, double *freq, double *quantile, 
	      double *cutoff);

/* cutoff_lr: (prop'n mismatches)
   use value of LR comparing unrelateds to siblings */
void cutoff_lr(int *nmar, double *freq, double *value, 
	       double *cutoff);

/* cutoff_llr: (Karl's LLR distance)
   uses a quantile from the distr'n of dist between unrelateds */
void cutoff_llr(int *nmar, double *freq, double *quantile, 
		double *cutoff);

/* cutoff_llr_s: (Karl's LLR distance)
   uses a quantile from the distr'n of dist between siblings */
void cutoff_llr_s(int *nmar, double *freq, double *quantile, 
		  double *cutoff);

/* cutoff_llr_meansib: (Karl's LLR distance)
   uses the mean distance between siblings */
void cutoff_llr_meansib(int *nmar, double *freq, double *cutoff);

/* cutoff_llr_lr: (Karl's LLR distance)
   use value of LR comparing unrelateds to siblings */
void cutoff_llr_lr(int *nmar, double *freq, double *value,
		   double *cutoff);

/* find point x at which dnorm(x,m1,s1)/dnorm(x,m2,s2) = value */
/* using a simple bisection method */
double find_lr(double m1, double s1, double m2, double s2, double value);

/* end of cutoffs.h */
