/*
 * Christian Gaser
 * $Id: vollib.h 207 2012-09-27 22:06:31Z gaser $ 
 *
 */

#ifndef _VOLLIB_H_
#define _VOLLIB_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <float.h>

#define SQRT2PI 2.506628

#define RINT(A) floor((A)+0.5)
#ifndef isfinite
#define isfinite(x) ((x) * (x) >= 0.) /* check for NaNs */
#endif

#ifndef TINY
#define TINY 1e-15 
#endif

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

#ifndef MIN
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#endif

#ifndef ROUND
#define ROUND( x ) ((long) ((x) + ( ((x) >= 0) ? 0.5 : (-0.5) ) ))
#endif

#ifndef isnan
#define isnan(a) ((a)!=(a)) 
#endif

#define index(A,B,C,DIM) ((C)*DIM[0]*DIM[1] + (B)*DIM[0] + (A))

#define GM     0
#define WM     1
#define CSF    2
#define SKULL1 3
#define SKULL2 4
#define SKULL3 5

void get_prctile(float *src, int *dims, double threshold[2], double prctile[2], int exclude_zeros);
void watershed3d(unsigned char *img, unsigned char *marker, int flag_dilate, int *dims);
void morph_erode_uint8(unsigned char *vol, int dims[3], int niter, double th);
void morph_dilate_uint8(unsigned char *vol, int dims[3], int niter, double th);
void morph_close_uint8(unsigned char *vol, int dims[3], int niter, double th);
void morph_open_uint8(unsigned char *vol, int dims[3], int niter, double th);
void morph_close_float(float *vol, int dims[3], int niter, double th);
void morph_open_float(float *vol, int dims[3], int niter, double th);
void subsample_double(double *in, double *out, int dim_in[3], int dim_out[3], int offset_in, int offset_out);
void subsample_float(float *in, float *out, int dim_in[3], int dim_out[3], int offset_in, int offset_out);
void subsample_uint8(unsigned char *in, float *out, int dim_in[3], int dim_out[3], int offset_in, int offset_out);
void vol_approx(float *vol, int dims[3], double voxelsize[3], int samp);
void smooth_double(double *vol, int dims[3], double voxelsize[3], double s[3], int use_mask);
void smooth_float(float *vol, int dims[3], double voxelsize[3], double s[3], int use_mask);
void smooth_subsample_double(double *vol, int dims[3], double voxelsize[3], double s[3], int use_mask, int samp);
void smooth_subsample_float(float *vol, int dims[3], double voxelsize[3], double s[3], int use_mask, int samp);
void initial_cleanup(unsigned char *probs, unsigned char *mask, int dims[3], double voxelsize[3], int strength, int remove_sinus);
void cleanup(unsigned char *probs, unsigned char *mask, int dims[3], double voxelsize[3], int strength, int gmwm_only);
void median3_uint8(unsigned char *D, int *dims);
void median3_float(float *D, int *dims);
void get_largest_component(unsigned char *bw, int dim[3]);
void laplace3(float *SEG, int dims[3], int maxiter);

#endif