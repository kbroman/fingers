/**********************************************************************
 * 
 * llrdist.c
 *
 * copyright (c) 2001, Karl W Broman
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C function to calculate the distance matrix using the log likelihood
 * ratio in place of simply proportion mismatches for the Fingers 
 * program.
 *
 * Contains: llrdist_wrap, llrdist, calc_llr_dist
 *           reorg_imatrix, reorg_dmatrix
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include "llrdist.h"

#define TOL 1e-12

void llrdist_wrap(int *n_worms, int *n_mar, int *data, 
		  double *freq, double *dist)
{
  llrdist(*n_worms, *n_mar, data, freq, dist);
}



void llrdist(int n_worms, int n_mar, int *data, double *freq, double *dist)  
{
  int i, j, k, **Data;
  double **Dist, d[3];

  reorg_imatrix(n_worms, n_mar, data, &Data);
  reorg_dmatrix(n_worms, n_worms, dist, &Dist);

  for(k=0; k<n_mar; k++) {
    calc_llr_dist(d, freq[k]);

    for(i=0; i < n_worms; i++) {
      for(j=i; j < n_worms; j++) {
	if(Data[k][i]==Data[k][j]) {
	  if(Data[k][i]==0) {
	    Dist[i][j] += d[0]; /* neither has band */
	  }
	  else {
	    Dist[i][j] += d[2]; /* both have band */
	  }
	}
	else Dist[i][j] += d[1]; /* mismatch */
      }
    }
  }

  /* fill out matrix; divide by number of markers */
  for(i=0; i<n_worms; i++) 
    for(j=i; j<n_worms; j++) 
      Dist[j][i] = Dist[i][j] = Dist[i][j]/(double)n_mar;

}


void calc_llr_dist(double *d, double p)
{
  double q, a, b, c, z;
  if(p > 1.0 - TOL) p = 1.0 - TOL;
  if(p < TOL) p = TOL;

  q = 1.0 - p;

  /* neither have band */
  a = 4.0*log(q) - log(0.25*p*p*q*q+p*q*q*q+q*q*q*q);

  /* mismatch */
  b = log(2.0)+2.0*log(q)+log(1-q*q) - 
    log(1.0 - p*p - 2.0*p*p*p*q - 3.5*p*p*q*q - 2.0*p*q*q*q - q*q*q*q);

  /* both have band */
  c = 2.0*log(1.0-q*q) - log(p*p+2.0*p*p*p*q+13.0/4.0*p*p*q*q+p*q*q*q);

  /* find minimum */
  z = a;
  if(b < z) z = b;
  if(c < z) z = c;

  /* distances, with min > 0 */
  d[0] = a-z;
  d[1] = b-z;
  d[2] = c-z;

}




/**********************************************************************
 * 
 * reorg_imatrix
 *
 * Reorganize a matrix that is stored as a vector (by columns) as 
 * a double indexed matrix
 *
 * Afterwards, matrix indexed like Matrix[col][row]
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/

void reorg_imatrix(int n_row, int n_col, int *matrix, int ***Matrix)
{
  int i;

  *Matrix = (int **)R_alloc(n_col, sizeof(int *));

  (*Matrix)[0] = matrix;
  for(i=1; i< n_col; i++) 
    (*Matrix)[i] = (*Matrix)[i-1] + n_row;
}

/**********************************************************************
 * 
 * reorg_dmatrix
 *
 * Reorganize a matrix that is stored as a vector (by columns) as 
 * a double indexed matrix
 *
 * Afterwards, matrix indexed like Matrix[col][row]
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/

void reorg_dmatrix(int n_row, int n_col, double *matrix, double ***Matrix)
{
  int i;

  *Matrix = (double **)R_alloc(n_col, sizeof(double *));

  (*Matrix)[0] = matrix;
  for(i=1; i< n_col; i++) 
    (*Matrix)[i] = (*Matrix)[i-1] + n_row;
}

/* end of llrdist.c */
