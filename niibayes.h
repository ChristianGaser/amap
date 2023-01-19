/*
 * Christian Gaser
 * $Id: niibayes.h 209 2013-04-12 08:52:39Z gaser $ 
 *
 */

#ifndef _NIISKULLSTRIP_H_
#define _NIISKULLSTRIP_H_

#define RIGID 0
#define AFFINE 1

#include "_reg_resampling.h"
#include "_reg_affineTransformation.h"
#include "_reg_blockMatching.h"
#include "_reg_tools.h"

#ifdef _W32
    #include <time.h>
#endif

#define PrecisionTYPE float

typedef struct{
  char *targetImageName;
  char *sourceImageName;
  char *outputResultName;
  char *outputMaskName;
  char *targetMaskName;
  char **tpmImageName;

  int maxIteration;

  int backgroundIndex[3];
  PrecisionTYPE sourceBGValue;

  float targetSigmaValue;
  float sourceSigmaValue;
  int levelNumber;
  int level2Perform;

  int block_percent_to_use;
  int inlier_lts;
}PARAM;

typedef struct{
  bool targetImageFlag;
  bool sourceImageFlag;
  bool tpmImageFlag;
  bool levelNumberFlag;
  bool level2PerformFlag;
  bool outputResultFlag;
  bool outputMaskFlag;
  bool targetMaskFlag;

  bool maxIterationFlag;

  bool backgroundIndexFlag;
  
  bool alignCenterFlag;

  bool rigidFlag;
  bool affineFlag;

  bool targetSigmaFlag;
  bool sourceSigmaFlag;
  bool pyramidFlag;
  bool useGPUFlag;
  bool twoDimRegistration;
  bool verbose;
}FLAG;

mat44 *affineRegistration(PARAM *param, FLAG *flag);
template <class DTYPE> 
void float_to_ptr_4D2(nifti_image *image, float *vol4d, int index_nt);
void float_to_ptr_4D(nifti_image *image, float *vol4d, int index_nt);
template <class DTYPE> 
void double_to_ptr_4D2(nifti_image *image, double *vol4d, int index_nt);
void double_to_ptr_4D(nifti_image *image, double *vol4d, int index_nt);
template <class DTYPE> 
void float_to_ptr_4D2(nifti_image *image, float *vol4d, int index_nt);
void float_to_ptr_4D(nifti_image *image, float *vol4d, int index_nt);
template <class DTYPE> 
void ptr_to_uint8_4D2(nifti_image *image, unsigned char *vol4d, int index_nt);
void ptr_to_uint8_4D(nifti_image *image, unsigned char *vol4d, int index_nt);
template <class DTYPE> 
void ptr_to_double_4D2(nifti_image *image, double *vol4d, int index_nt);
void ptr_to_double_4D(nifti_image *image, double *vol4d, int index_nt);
template <class DTYPE> 
void ptr_to_float_4D2(nifti_image *image, float *vol4d, int index_nt);
void ptr_to_float_4D(nifti_image *image, float *vol4d, int index_nt);


extern "C" {
#include "niilib.h"
#include "vollib.h"
#include "Amap.h"
void Bayes(float *src, unsigned char *label, unsigned char *priors, unsigned char *prob, double *separations, int *dims, double bias_fwhm, double bias_rate, int do_warp);
}

#endif