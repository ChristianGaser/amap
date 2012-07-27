/*
 * Christian Gaser
 * $Id$ 
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

void
morph_erode_uint8(unsigned char *vol, int dims[3], int niter, int th);

void
morph_dilate_uint8(unsigned char *vol, int dims[3], int niter, int th);

void
morph_close_uint8(unsigned char *vol, int dims[3], int niter, int th);

void
morph_open_uint8(unsigned char *vol, int dims[3], int niter, int th);

void
morph_close_double(double *vol, int dims[3], int niter, double th);

void
morph_open_double(double *vol, int dims[3], int niter, double th);

void 
subsample_double(double *in, double *out, int dim_in[3], int dim_out[3], int offset_in, int offset_out);

void 
subsample_float(float *in, float *out, int dim_in[3], int dim_out[3], int offset_in, int offset_out);

void 
subsample_uint8(unsigned char *in, float *out, int dim_in[3], int dim_out[3], int offset_in, int offset_out);

void
smooth_double(double *vol, int dims[3], double separations[3], double s[3], int use_mask);

void
smooth_float(float *vol, int dims[3], double separations[3], double s[3], int use_mask);

void
smooth_subsample_double(double *vol, int dims[3], double separations[3], double s[3], int use_mask, int samp);

void
smooth_subsample_float(float *vol, int dims[3], double separations[3], double s[3], int use_mask, int samp);

#endif