/*
 * Christian Gaser
 * $Id: niiamap.h 213 2014-10-27 10:24:50Z gaser $ 
 *
 */

#ifndef _NIIAMAP_H_
#define _NIIAMAP_H_

char *mask_filename = NULL;
int n_pure_classes = 3;
int iters_amap = 200;
int subsample = 16;
int iters_nu = 20;
int iters_ICM = 50;
int pve = 5;
int correct_nu = 1;
int write_seg[3] = {0, 1, 0};
int write_nu = 0;
int write_label = 1;
int debug = 0;
#ifdef SPLINESMOOTH
  double bias_fwhm = 500.0;
#else
  double bias_fwhm = 60.0;
#endif
double thresh_brainmask = 0.05;
double thresh_kmeans = 0.5;
double weight_MRF = 0.15;

#endif