/*
 * Christian Gaser
 * $Id: niibrainmask.h 167 2011-12-27 09:28:29Z gaser $ 
 *
 */

#include <float.h>

#include "vollib.h"
#include "niilib.h"

#define NOPVE 0
#define KMEANS 1

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

#ifndef MIN
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#endif

extern double Kmeans(double *src, unsigned char *label, unsigned char *mask, int NI, int n_clusters, double *separations, int *dims, int thresh_mask, int thresh_kmeans, int iters_nu, int pve, double bias_fwhm);
extern void get_largest_cluster(unsigned char *bw, int dim[3]);

//int iters_nu = 40;
//int write_nu = 0;
int write_masked = 1;
int write_brainmask = 0;
int strip_param[4] = {2,5,1,2}; 
