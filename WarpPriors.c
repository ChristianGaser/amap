/*
 * Christian Gaser
 * $Id: WarpPriors.c 217 2020-09-25 08:49:06Z gaser $ 
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "optimizer3d.h"
#include "diffeo3d.h"
#include "Amap.h"
#include "vollib.h"

struct dartel_prm {
  int rform;         /* regularization form: 0 - linear elastic energy; 1 - membrane energy; 2 - bending energy */
  double rparam[6];  /* regularization parameters */
  double lmreg;      /* LM regularization */
  int cycles;        /* number of cycles for full multi grid (FMG) */
  int its;           /* Number of relaxation iterations in each multigrid cycle */
  int k;             /* time steps for solving the PDE */
  int code;          /* objective function: 0 - sum of squares; 1 - symmetric sum of squares; 2 - multinomial */
};


void WarpPriors(unsigned char *prob, unsigned char *priors, float *flow, int *dims, int n_loops, int n_classes, int samp)
{
  int vol_samp, vol, i, j;
  int size_samp[4], size[4], area;
  double buf[6], ll[3]; 
  float *f, *g, *v, *flow1, *flow2, *scratch, max, *priors_float;
  int it, it0, it1, it_scratch, ndims4, dims_samp[3], area_samp;
   
  int code = 2;    /* multinomial */
  int rform = 0;   /* linear energy */
  double lmreg = 0.01;
  static double param[3] = {1.0, 1.0, 1.0};
    
  struct dartel_prm* prm = (struct dartel_prm*)malloc(sizeof(struct dartel_prm)*10);

  /* only use gm/wm */
  ndims4 = 2;

  /* define grid dimensions */
  for(j=0; j<3; j++) dims_samp[j] = (int) ceil((dims[j]-1)/((double) samp))+1;

  area = dims[0]*dims[1];
  vol  = dims[0]*dims[1]*dims[2];
  
  area_samp = dims_samp[0]*dims_samp[1];
  vol_samp  = dims_samp[0]*dims_samp[1]*dims_samp[2];
  
  f      = (float *)malloc(sizeof(float)*n_classes*vol_samp);
  g      = (float *)malloc(sizeof(float)*n_classes*vol_samp);
  v      = (float *)malloc(sizeof(float)*3*vol_samp);
  
  /* initialize size of subsampled data and add 4th dimension */
  for (i=0; i < 3; i++) size_samp[i] = dims_samp[i];    
  size_samp[3] = ndims4;
  for (i=0; i < 3; i++) size[i] = dims[i];    
  size[3] = ndims4;
  
  /* some entries are equal */
  for (j = 0; j < n_loops; j++) {
    for (i = 0; i < 3; i++)  prm[j].rparam[i] = param[i];
    prm[j].rform = rform;
    prm[j].cycles = 3;
    prm[j].its = 3;
    prm[j].code = code;
    prm[j].lmreg = lmreg;
  }

  prm[0].rparam[3] = 4.0;   prm[0].rparam[4] = 2.0;    prm[0].rparam[5] = 1e-6; prm[0].k = 0; 
  prm[1].rparam[3] = 2.0;   prm[1].rparam[4] = 1.0;    prm[1].rparam[5] = 1e-6; prm[1].k = 0; 
  prm[2].rparam[3] = 1.0;   prm[2].rparam[4] = 0.5;    prm[2].rparam[5] = 1e-6; prm[2].k = 1; 
  prm[3].rparam[3] = 0.5;   prm[3].rparam[4] = 0.25;   prm[3].rparam[5] = 1e-6; prm[3].k = 2; 
  prm[4].rparam[3] = 0.25;  prm[4].rparam[4] = 0.125;  prm[4].rparam[5] = 1e-6; prm[4].k = 4; 
  prm[5].rparam[3] = 0.125; prm[5].rparam[4] = 0.0625; prm[5].rparam[5] = 1e-6; prm[5].k = 6;

  /* use different parameters for bending energy */
  if(rform == 2) {
    prm[0].rparam[3] = 8.0;   prm[0].rparam[4] = 1e-4*prm[0].rparam[3];    prm[0].rparam[5] = 1e-4; prm[0].k = 0; 
    prm[1].rparam[3] = 4.0;   prm[1].rparam[4] = 1e-4*prm[1].rparam[3];    prm[1].rparam[5] = 1e-4; prm[1].k = 0; 
    prm[2].rparam[3] = 2.0;   prm[2].rparam[4] = 1e-4*prm[2].rparam[3];    prm[2].rparam[5] = 1e-5; prm[2].k = 1; 
    prm[3].rparam[3] = 1.0;   prm[3].rparam[4] = 1e-4*prm[3].rparam[3];    prm[3].rparam[5] = 1e-5; prm[3].k = 2; 
    prm[4].rparam[3] = 0.5;   prm[4].rparam[4] = 1e-4*prm[4].rparam[3];    prm[4].rparam[5] = 1e-6; prm[4].k = 4; 
    prm[5].rparam[3] = 0.25;  prm[5].rparam[4] = 1e-4*prm[5].rparam[3];    prm[5].rparam[5] = 1e-6; prm[5].k = 6;
  }
  
  /* subsample priors and probabilities to lower resolution */
  if (samp==1) {
      for (i = 0; i < vol*n_classes; i++) {
        f[i] = priors[i];
        g[i] = prob[i];
      }
  } else {
    for (i = 0; i < n_classes; i++) {
      subsample_uint8(priors, f, dims, dims_samp, i*vol, i*vol_samp);    
      subsample_uint8(prob  , g, dims, dims_samp, i*vol, i*vol_samp);   
    } 
  }

  /* subsample initial flow field to lower resolution */
  if (samp==1) {
      for (i = 0; i < vol*3; i++)
        v[i] = flow[i];
  } else {
    for (i = 0; i < 3; i++) {
      subsample_float(flow, v, dims, dims_samp, i*vol, i*vol_samp);    
    }
  }

  /* scale subsampled probabilities to a maximum of 0.5 */
  max = -HUGE;
  for (i=0; i < n_classes*vol_samp; i++) max = MAX(g[i], max);
  for (i=0; i < n_classes*vol_samp; i++) g[i] /= max*2.0;

  /* scale subsampled priors to a maximum of 0.5 */
  max = -HUGE;
  for (i=0; i < n_classes*vol_samp; i++) max = MAX(f[i], max);
  for (i=0; i < n_classes*vol_samp; i++) f[i] /= max*2.0;

  /* iterative warping using dartel approach */
  it = 0;
  for (it0 = 0; it0 < n_loops; it0++) {
    it_scratch = iteration_scratchsize((int *)size_samp, prm[it0].code, prm[it0].k);
    scratch = (float *)malloc(sizeof(float)*it_scratch);

    for (it1 = 0; it1 < prm[it0].its; it1++) {
      it++;
      iteration(size_samp, prm[it0].k, v, f, g, (float *)0, prm[it0].rform, prm[it0].rparam, prm[it0].lmreg, 
        prm[it0].cycles, prm[it0].its, prm[it0].code, flow, ll, scratch);              
      printf("%02d:\t%.2f\n", it, ll[0]);
      fflush(stdout);
      for (i = 0; i < 3*vol_samp; i++) v[i] = flow[i];
    }
    free(scratch);
  }
  free(f);
  free(g);
  
  /* upsample flow field */
  flow2  = (float *)malloc(sizeof(float)*3*vol);
  if (samp==1) {
      for (i = 0; i < vol*3; i++) 
        flow2[i] = v[i];
  } else {
    for (i = 0; i < 3; i++) {
      subsample_float(v, flow2, dims_samp, dims, i*vol_samp, i*vol);    
    }
  }

  free(v);
  
  /* rescale flow field */
  for (i = 0; i < vol; i++) {
    flow2[i] /= (double)dims_samp[0]/(double)dims[0]; 
    flow2[i + vol] /= (double)dims_samp[1]/(double)dims[1]; 
    flow2[i + 2*vol] /= (double)dims_samp[2]/(double)dims[2]; 
  }
  
  /* use exponentional flow */
  flow1  = (float *)malloc(sizeof(float)*3*vol);
  expdef(size, 6, -1, flow2, flow, flow1, (float *)0, (float *)0); 
  free(flow1);
  
  /* copy floating priors for sampn */
  priors_float = (float *)malloc(sizeof(float)*n_classes*vol);
  for (i = 0; i < n_classes*vol; i++) priors_float[i] = (float)priors[i];

  /* apply deformation field to priors */
  for (i = 0; i < vol; i++) {
    sampn(dims, priors_float, n_classes, vol, (double)flow[i]-1.0, (double)flow[i+vol]-1.0, (double)flow[i+(2*vol)]-1.0, buf);
    for (j = 0; j < n_classes; j++) priors[i + (j*vol)] = (unsigned char)MIN(255,ROUND(buf[j]));
  }

  /* rescue flow field */
  for (i = 0; i < 3*vol; i++) flow[i] = flow2[i];

  free(prm);
  free(flow2);
  free(priors_float);

}
