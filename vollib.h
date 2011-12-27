/*
 * Christian Gaser
 * $Id: vollib.h 166 2011-12-27 08:43:49Z gaser $ 
 *
 */

#ifndef _VOLLIB_H_
#define _VOLLIB_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <float.h>

#define RINT(A) floor((A)+0.5)
#ifndef isfinite
#define isfinite(x) ((x) * (x) >= 0.) /* check for NaNs */
#endif

extern void
morph_erode_uint8(unsigned char *vol, int dims[3], int niter, int th);

extern void
morph_dilate_uint8(unsigned char *vol, int dims[3], int niter, int th);

extern void
morph_close_uint8(unsigned char *vol, int dims[3], int niter, int th);

extern void
morph_open_uint8(unsigned char *vol, int dims[3], int niter, int th);

extern void
morph_close_double(double *vol, int dims[3], int niter, double th);

extern void
morph_open_double(double *vol, int dims[3], int niter, double th);

extern void 
subsample_double(double *in, double *out, int dim_in[3], int dim_out[3]);

extern void 
subsample_float(float *in, float *out, int dim_in[3], int dim_out[3]);

extern void
smooth_double(double *vol, int dims[3], double separations[3], double s[3], int use_mask);

extern void
smooth_float(float *vol, int dims[3], double separations[3], double s[3], int use_mask);

extern void
smooth_subsample_double(double *vol, int dims[3], double separations[3], double s[3], int use_mask, int samp);

extern void
smooth_subsample_float(float *vol, int dims[3], double separations[3], double s[3], int use_mask, int samp);

#endif