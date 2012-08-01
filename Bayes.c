/*
 * Christian Gaser
 * $Id$ 
 *
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "Amap.h"
#include "vollib.h"

/* use always 6 classes */
#define Kb 6
#define MAXK 30

void WarpPriors(unsigned char *prob, unsigned char *priors, float *flow, int *dims, int n_loops, int n_classes, int samp);

void Bayes(float *src, unsigned char *label, unsigned char *priors, unsigned char *prob, double *separations, int *dims, int correct_nu)
{
  int i, j, k, l, k1, subit, subsample_warp, kmax, k2;
  int count, subsample, masked_smoothing;
  double ll = -HUGE, llr=0.0, fwhm[3];
  float *nu, *meaninvcov, *resmean;
  double mn[MAXK], vr[MAXK], mg[MAXK], mn2[MAXK], vr2[MAXK], mg2[MAXK];
  double mnt[MAXK], vrt[MAXK];
  double mom0[MAXK], mom1[MAXK], mom2[MAXK], mgm[MAXK];
  double q[MAXK], bt[MAXK], b[MAXK], qt[MAXK];
  double tol1 = 1e-3, bias_fwhm, sum, sq, qmax, s, psum, p1, oll;
  float *flow;
  
  int area = dims[0]*dims[1];
  int vol = area*dims[2];
  int K = 0;
  int do_warp = 0;
  
  /* multiple gaussians are not yet working */
  int ngauss[6]   = {2,2,2,3,3,3};
  int iters_EM[5] = {5, 10, 10, 10, 10};

  int histo[65536], lkp[MAXK];
  double mn_thresh, mx_thresh;
  double min_src = HUGE, max_src = -HUGE;
  int cumsum[65536], order_priors[6] = {2,0,1,3,4,5};
  int n_loops = 4;

  bias_fwhm = 30.0;

  subsample_warp = (int)round(9.0/(separations[0]+separations[1]+separations[2]));
  subsample_warp = 1;
  fprintf(stderr,"subsample by factor %d\n",subsample_warp);

  if (do_warp) {
    flow  = (float *)malloc(sizeof(float)*vol*3);
    /* initialize flow field with zeros */
    for (i = 0; i < (vol*3); i++) flow[i] = 0.0;
  }
  
  if (correct_nu) {
    nu         = (float *)malloc(sizeof(float)*vol);
    meaninvcov = (float *)malloc(sizeof(float)*vol);
    resmean    = (float *)malloc(sizeof(float)*vol);
  }

  for (i=0; i<vol; i++) {
    min_src = MIN((double)src[i], min_src);
    max_src = MAX((double)src[i], max_src);
  }

  /* fill areas where no prior is defined with background */
  for (i=0; i<vol; i++) {
    sum = 0.0;
    for (j=0; j<6; j++)
      sum += (double)priors[i+j*vol];
    if (sum == 0.0) priors[i+5*vol] = 255;
  }

  /* build histogram */
  for (i = 0; i < 65536; i++) histo[i] = 0;
  for (i=0; i<vol; i++) {
    if (src[i] == 0) continue;
    histo[(int)ROUND(65535.0*((double)src[i]-min_src)/(max_src-min_src))]++;
  }

  /* find values between 0.1% and 99.9% quartile */
  cumsum[0] = histo[0];
  for (i = 1; i < 65536; i++) cumsum[i] = cumsum[i-1] + histo[i];
  for (i = 0; i < 65536; i++) cumsum[i] = (int) ROUND(1000.0*(double)cumsum[i]/(double)cumsum[65535]);
  for (i = 0; i < 65536; i++) if (cumsum[i] >= 1) break;
  mn_thresh = (double)i/65535.0*(max_src-min_src);
  for (i = 65535; i > 0; i--) if (cumsum[i] <= 999) break;
  mx_thresh = (double)i/65535.0*(max_src-min_src);

  /* K = sum(ngauss) */
  for (k1=0; k1<Kb; k1++) K += ngauss[k1]; 

  /* build lkp */
  l = 0;
  for (k1=0; k1<Kb; k1++) {
    for (j=0; j<ngauss[k1]; j++) {
      lkp[l] = k1; 
      l++;
    }
  }

  /* initialize mean, var, mg */
  for (k=0; k<K; k++) {
    mn[k] = mx_thresh * 1.0/(double)(k+1);
    mg[k] = 1.0/(double)K;
    vr[k] = mx_thresh*mx_thresh + TINY;
  }
          
  count = 0;
  /* start with a few EM iterations and after nu-corection use more iterations */
  for (j=0; j < 3; j++) {
    for (subit=0; subit<iters_EM[j]; subit++) {
      oll = ll;
      ll = llr;
      for (k=0; k<K; k++) mom0[k] = mom1[k] = mom2[k] = 0.0;
      for (k1=0; k1<Kb; k1++) mgm[k1] = 0.0;
    
      for (k=0; k<K; k++) {
        mg2[k] = mg[k];
        mn2[k] = mn[k];
        vr2[k] = vr[k];
      }
      
      for (i=0; i<vol; i++) {
      
        if (correct_nu) {
          meaninvcov[i] = resmean[i] = nu[i] = 0.0;
        }

        if(((double)src[i]>mn_thresh) && ((double)src[i]<mx_thresh)) {
          s = TINY;
          k2 = 0;
          for (k1=0; k1<Kb; k1++) {
            bt[k1] = (double) priors[i+(vol*k1[order_priors])]; 
            for (l=0; l<ngauss[k1]; l++) {
              s += bt[k1]*mg[k2];
              k2++;
            }
          }

          for (k1=0; k1<Kb; k1++) {
            bt[k1] /= s;
            mgm[k1] += bt[k1];
          }

          for (k=0; k<K; k++)
            q[k] = mg[k]*bt[lkp[k]]*exp(SQR((double)src[i]-mn[k])/(-2*vr[k]))/(SQRT2PI*sqrt(vr[k]));
          
          /* prepare prob for warping */
          sq = TINY;

          for (k1=0; k1<Kb; k1++)
            qt[k1] = mnt[k1] = vrt[k1] = 0.0;
          
          /* prepare prob, mean and var for 6 classes */
          for (k=0; k<K; k++) { 
            qt[lkp[k]] += q[k];
            mnt[lkp[k]] += mn[k];
            vrt[lkp[k]] += SQR(vr[k])/ngauss[lkp[k]];
            sq += qt[lkp[k]];
          }

          /* mean of variance estimated by euclidian distance */
          for (k1=0; k1<Kb; k1++) vrt[k1] = sqrt(vrt[k1]);

          if (correct_nu) {
          /* don't use all classes */
            for (k1=1; k1<Kb-1; k1++) {
              float tempf   = qt[k1]/sq/vrt[k1];
              resmean[i]    += tempf*(src[i] - mnt[k1]/ngauss[k1]);
              meaninvcov[i] += tempf;
            }
          }

          psum = 0.0;
          qmax = -HUGE;
          for (k1=0; k1<Kb; k1++) {
            p1 = qt[k1]/sq;
            /* save probabilities */
            prob[i+(vol*k1[order_priors])] = (unsigned char)ROUND(255.0*p1);
            psum += p1;
            if (p1 > qmax) { 
              qmax = p1;     
              kmax = k1 + 1;
            }   
          }

          if (psum>0) label[i] = kmax;
          else        label[i] = 0;

          ll += log10(sq);
          for (k=0; k<K; k++) {
            p1 = q[k]/sq;           mom0[k] += p1;
            p1 *= (double)src[i];   mom1[k] += p1;
            p1 *= (double)src[i];   mom2[k] += p1;
          }      
        }
      } 
      
      if (do_warp && (subit==0) && (j==1)) WarpPriors(prob, priors, flow, dims, n_loops, Kb, subsample_warp);
          
      for (k=0; k<K; k++) {
        mg[k] = (mom0[k]+TINY)/(mgm[lkp[k]]+TINY);
        mn[k] = mom1[k]/(mom0[k]+TINY);
        vr[k] = (mom2[k]-SQR(mom1[k])/mom0[k]+1e6*TINY)/(mom0[k]+TINY);
        vr[k] = MAX(vr[k],TINY);
        
        /* rescue previous values in case of nan */
#if defined(_WIN32)
        if ((_isnan(mg[k])) || (_isnan(mn[k])) || (_isnan(vr[k]))) {
#else
        if ((isnan(mg[k])) || (isnan(mn[k])) || (isnan(vr[k]))) {
#endif
          mg[k] = mg2[k];
          mn[k] = mn2[k];
          vr[k] = vr2[k];
          break;
        }
      }
      
      printf("%7.4f\b\b\b\b\b\b\b",ll/vol);    
      printf("\n");
      fflush(stdout);
      if(fabs(ll-oll)<tol1*vol) count++;
      if (count > 2) break;
    }
    	  
    if (correct_nu) {
      /* smoothing of residuals */
      for(i=0; i<3; i++) {
        if (j==0) fwhm[i] = 1.5*bias_fwhm;
        else fwhm[i] = bias_fwhm;
      }
        
      /* use subsampling for faster processing */
      subsample = subsample_warp;
      masked_smoothing = 0;

      smooth_subsample_float(meaninvcov, dims, separations, fwhm, masked_smoothing, subsample);
      smooth_subsample_float(resmean, dims, separations, fwhm, masked_smoothing, subsample);
      
      float numean = 0.0;
      for (i = 0; i < vol; i++) {
        if (meaninvcov[i] != 0) 
          nu[i] = (resmean[i]/meaninvcov[i]) - 1.0;
        else nu[i] = 0.0;
        numean += nu[i];
      }

      numean /= vol;
      for (i = 0; i < vol; i++) {
        nu[i] -= numean;
        src[i] -= nu[i];
      }
      
    }

  }

  for (i=0; i<vol; i++) {
    if (src[i] != 0) {
      for (k1=0; k1<Kb; k1++) {
        bt[k1] = (double) priors[i+(vol*k1[order_priors])]; 
        qt[k1] = 0.0;
      }
      
      for (k=0; k<K; k++) 
        b[k] = bt[lkp[k]]*mg[k];

      for (k=0; k<K; k++) { 
        p1 = exp(SQR((double)src[i]-mn[k])/(-2*vr[k]))/(SQRT2PI*sqrt(vr[k]));
        qt[lkp[k]] += p1*b[k];
      }
      
      /* ensure that sum of prob is 1 */
      s = TINY;
      for (k1=0; k1<Kb; k1++)
        s += qt[k1];

      qmax = -HUGE;
      psum = 0.0;
      for (k1=0; k1<Kb; k1++) {
        p1 = qt[k1]/s;
        /* save probabilities */
        prob[i+(vol*k1[order_priors])] = (unsigned char)ROUND(255.0*p1);
        psum += p1;
        if (p1 > qmax) { 
          qmax = p1;     
          kmax = k1 + 1;
        }   
      }
      
      /* label only if sum of all probabilities is > 0 */
      if (psum > 0)
        label[i] = kmax;
      else
        label[i] = 0;
    } else label[i] = 0;
  }

  for (k=0; k<Kb; k++) 
    printf("%g %g\n",mnt[k],sqrt(vrt[k]));
    

  if (do_warp) free(flow);
  if (correct_nu) {
    free(nu);
    free(meaninvcov);
    free(resmean);
  }
//  for (i=0; i<vol; i++) src[i] = nu[i];
  
    
  return;
}  
