/*
 * Christian Gaser
 * $Id$ 
 *
 */

#include "nifti/nifti1_io.h"
#include "nifti/nifti1_local.h"

#include <float.h>

#define NOPVE 0
#define KMEANS 1

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

#ifndef MIN
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#endif

extern double Kmeans(double *src, unsigned char *label, unsigned char *mask, int NI, int n_clusters, double *separations, int *dims, int thresh_mask, int thresh_kmeans, int iters_nu, int pve, double bias_fwhm);
extern nifti_image *read_nifti_double( const char *input_filename, double *image[]);
extern int splineSmooth( double *src, double lambda, double distance, int subsample, double *separations, int *dims);
extern void morph_dilate_uint8(unsigned char *vol, int dims[3], int niter, int th);
extern void morph_close_uint8(unsigned char *vol, int dims[3], int niter, int th);
extern void morph_open_uint8(unsigned char *vol, int dims[3], int niter, int th);
extern void get_largest_cluster(unsigned char *bw, int dim[3]);
extern write_nifti( const char *output_filename, double image[], int data_type, double slope, int dim[], double vox[], nifti_image *in_ptr);

//int iters_nu = 40;
//int write_nu = 0;
int write_masked = 1;
int write_brainmask = 0;
int strip_param[4] = {2,5,1,2}; 
