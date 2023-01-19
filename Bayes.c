/*
 * Christian Gaser
 * $Id: Bayes.c 209 2013-04-12 08:52:39Z gaser $ 
 *
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "vollib.h"

/* use always 6 classes */
#define Kb 6
#define MAXK 30

void WarpPriors(unsigned char *prob, unsigned char *priors, float *flow, int *dims, int n_loops, int n_classes, int samp);

int latent();

void Bayes(float *src, unsigned char *label, unsigned char *priors, unsigned char *prob, double *separations, int *dims, double bias_fwhm, double bias_rate, int do_warp)
{
  int i, j, k, l, m, k1, iter, subit, subsample_warp, kmax, k2;
  int count, correct_nu;
  double ll = -HUGE, llr=0.0, fwhm[3], wp_reg, summom0, sumwp;
  float *nu;
  double mn[MAXK], vr[MAXK], mg[MAXK], mn2[MAXK], vr2[MAXK], mg2[MAXK], wp[Kb];
  double mnt[MAXK], vrt[MAXK];
  double mom0[MAXK], mom1[MAXK], mom2[MAXK], mgm[Kb], mgm2[Kb];
  double mm0[Kb], mm1[Kb], mm2[Kb];
  double mni[Kb], vri[Kb], mgi[Kb];
  double q[MAXK], bt[Kb], b[MAXK], qt[MAXK];
  double tol1 = 1e-4, sum, sq, qmax, s, psum, p1, oll;
  double vr0, m0, m1, m2;
  double threshold[2], prctile[2] = {0.1,99.9};
  float *flow;
  
  int area = dims[0]*dims[1];
  int vol = area*dims[2];
  int K = 0;
  
  int iters_EM[5] = {5,5,5,10,10};
  int ngauss0[Kb] = {1,1,1,1,1,1};
  int ngauss[Kb]  = {1,1,1,3,4,2};
  int lkp0[Kb]    = {1,2,3,4,5,6};
  int lkp[MAXK];
  int order_priors[Kb] = {2,0,1,3,4,5};
  int n_loops = 5;
  int max_iters = 1;
  
  wp_reg = 100; /* Bias wp towards 1/Kb */
  
  for (k1=0; k1<Kb; k1++) wp[k1] = 1.0/(double)Kb; 

  correct_nu = (bias_fwhm > 0) ? 1 : 0;
  
  subsample_warp = (int)round(9.0/(separations[0]+separations[1]+separations[2]));
  fprintf(stderr,"subsample by factor %d\n",subsample_warp);

  if (do_warp) {
    flow  = (float *)malloc(sizeof(float)*vol*3);
    /* initialize flow field with zeros */
    for (i = 0; i < (vol*3); i++) flow[i] = 0.0;
  }
  
  if (correct_nu) {
    nu = (float *)malloc(sizeof(float)*vol);
    fprintf(stderr,"Correct bias\n");
  }

  /* fill areas for last class where no prior is defined with background */
  for (i=0; i<vol; i++) {
    sum = 0.0;
    for (j=0; j<Kb; j++)
      sum += (double)priors[i+j*vol];
    if (sum == 0.0) priors[i+(Kb-1)*vol] = 255;
  }

  /* find values between 0.1% and 99.9% percentile */
  get_prctile(src, dims, threshold, prctile, 0);

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

  for (k1=0; k1<Kb; k1++) 
    mni[k1] = vri[k1] = mgi[k1] = mgm[k1] = mgm2[k1] = mm0[k] = mm1[k] = mm2[k] = 0.0;

  m0 = m1 = m2 = 0.0;
  for (i=0; i<vol; i++) {
    if(((double)src[i]>threshold[0]) && ((double)src[i]<threshold[1])) {
      m0++;
      m1 += (double)src[i];
      m2 += (double)SQR(src[i]);
      for (k1=0; k1<Kb; k1++) {
        float b = (double)priors[i+(vol*k1[order_priors])];
        mm0[k1] += b/255.0;
        mm1[k1] += b/255.0*src[i];
        mm2[k1] += SQR(src[i]);
      }      
    }
  }
  
  /* Construct a "Wishart-style prior" (vr0) */
  vr0 = (m2/m0 - SQR(m1/m0))/SQR(Kb);
  fprintf(stderr,"%g\n",vr0);
  
  double vr1 = 0;
  for (k1=0; k1<Kb; k1++) {
    mni[k1] = mm1[k1]/(mm0[k1]+TINY);
    vr1 += (mm2[k1] - mm1[k1]*mm1[k1]/mm0[k1]);
  }
  double sum_mm0 = 0.0;
  for (k1=0; k1<Kb; k1++)
    sum_mm0 += mm0[k1];
  
  /* initialize mean, var, mg */
  vr1 = (vr1+vr0)/(sum_mm0+1.0);
  for (k1=0; k1<Kb; k1++) {
    mgi[k1] = 1.0;
    vri[k1] = vr1;
  }

  /* extend the values to classes with multiple entries */
  double prev_mn = 0.0;
  for (k=0; k<K; k++) {
    mn[k] = mni[lkp[k]];
    mg[k] = mgi[lkp[k]];
    vr[k] = vri[lkp[k]];
    
    /* add some small value to prevent having the same values */
    mn[k] += 0.1*k*mn[k];
    vr[k] += 0.1*k*vr[k];
  }
  
  count = 0;
  /* start with a few EM iterations and after nu-corection use more iterations */
  for (iter=0; iter < max_iters; iter++) {
        
    for (subit=0; subit<iters_EM[iter]; subit++) {
      oll = ll;
      ll  = llr;
      for (k=0;  k<K;   k++)  mom0[k] = mom1[k] = mom2[k] = 0.0;
      for (k1=0; k1<Kb; k1++) mgm[k1] = mgm2[k1] = 0.0;
    
      for (k=0; k<K; k++) {
        mg2[k] = mg[k];
        mn2[k] = mn[k];
        vr2[k] = vr[k];
      }
      
      if (correct_nu)
        for (i=0; i<vol; i++) nu[i] = 0.0;

      for (i=0; i<vol; i++) {
      
        if(((double)src[i]>threshold[0]) && ((double)src[i]<threshold[1])) {
          s = TINY;
          k2 = 0;
          for (k1=0; k1<Kb; k1++) {
            bt[k1] = wp[k1]*(double) priors[i+(vol*k1[order_priors])]; 
            s += bt[k1];
/*            for (l=0; l<ngauss[k1]; l++) {
              s += bt[k1]*mg[k2];
              k2++;
            }
*/
          }

          for (k1=0; k1<Kb; k1++) {
            mgm[k1] += bt[k1]/s;
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

          if (correct_nu) {
          /* don't use all classes */
            for (k1=2; k1<3; k1++) {
              if ((label[i]-1) == k1)
                nu[i]    = src[i] - (mnt[k1]/ngauss[k1]);
            }
          }

          ll += log10(sq);
          for (k=0; k<K; k++) {
            p1 = q[k]/sq;           mom0[k] += p1;
            p1 *= (double)src[i];   mom1[k] += p1;
            p1 *= (double)src[i];   mom2[k] += p1;
          }      
        }
      } 
      
      if (do_warp && (subit==0) && (iter==1)) fprintf(stderr,"n_loops: %d\n",n_loops);

      if ((subit==0) && (iter>0)) {
        initial_cleanup(prob, label, dims, separations, 1, 1);
        cleanup(prob, label, dims, separations, 2, 1);
        if (do_warp)
          WarpPriors(prob, priors, flow, dims, n_loops, Kb, subsample_warp);
      }

      for (k=0; k<K; k++) {
        mgm2[lkp[k]] += mom0[k];
      }

      for (k=0; k<K; k++) {
        mg[k] = (mom0[k]+TINY)/(mgm2[lkp[k]]+TINY);
        mn[k] = mom1[k]/(mom0[k]+TINY);
        vr[k] = (mom2[k] - mom1[k]*mom1[k]/mom0[k] + vr0)/(mom0[k]+1);
        vr[k] = MAX(vr[k],TINY);
        fprintf(stderr,"%g\t%g\t%g\t%g\t%g\t%g\n",mn[k],mg[k],vr[k],mom0[k],mom1[k1],mom2[k]);
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
      fprintf(stderr,"\n");
/* weighting is not yet working */
if (1) {
      l = 0;
      sumwp = 0.0;
      for (k1=0; k1<Kb; k1++) {
        summom0 = 0.0;
        for (m=0; m<ngauss[k1]; m++) {
          summom0 += mom0[l];
          l++;
        }
        wp[k1] = (summom0 + wp_reg)/(mgm[k1] + wp_reg*Kb);
        sumwp += wp[k1];
      }
      for (k1=0; k1<Kb; k1++) {
        wp[k1] /= sumwp;
        fprintf(stderr,"%g ",wp[k1]);
      }
      fprintf(stderr,"\n");
}
        
      printf("%7.4f\b\b\b\b\b\b\b",ll/vol);    
      printf("\n");
      fflush(stdout);
      if(fabs(ll-oll)<tol1*vol) count++;
      if (count > 0) break;
    
        
//      if ((correct_nu) && (iter > 0)) {
      if ((correct_nu) && (iter > 0)) {
        /* smoothing of residuals */
        bias_fwhm /= bias_rate;
        for(i=0; i<3; i++)
          fwhm[i] = bias_fwhm;
            
        vol_approx(nu, dims, separations, 4);

        /* mean correct nu */
        float numean = 0.0;
        for (i = 0; i < vol; i++)
          numean += nu[i];

        numean /= vol;
        for (i = 0; i < vol; i++)
          src[i] -= (nu[i] - numean);
      }
    }
  }

  for (i=0; i<vol; i++) {
    if (src[i] != 0) {
      for (k1=0; k1<Kb; k1++) {
        bt[k1] = wp[k1]*(double) priors[i+(vol*k1[order_priors])]; 
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
    printf("%g %g\n",mnt[k]/ngauss[k],sqrt(vrt[k]));
    
  if (do_warp) free(flow);
  if (correct_nu) free(nu);
    
  return;
}  
