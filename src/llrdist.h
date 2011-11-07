/**********************************************************************
 * 
 * llrdist.h
 *
 * copyright (c) 2001, Karl W Broman
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C function to calculate the distance matrix using the log likelihood
 * ratio in place of simply proportion mismatches for the Fingers 
 * program.
 *
 * Contains: llrdist_wrap, llrdist, llr
 *           reorg_imatrix, reorg_dmatrix
 *
 **********************************************************************/

void llrdist_wrap(int *n_worms, int *n_mar, int *data, 
		  double *freq, double *dist);

void llrdist(int n_worms, int n_mar, int *data, double *freq, double *dist);

void calc_llr_dist(double *dist, double p);

void reorg_imatrix(int n_row, int n_col, int *matrix, int ***Matrix);

void reorg_dmatrix(int n_row, int n_col, double *matrix, double ***Matrix);

/* end of llrdist.h */
