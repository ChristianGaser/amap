#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

#ifndef MIN
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#endif

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

#ifndef SQRT2PI
#define SQRT2PI 2.506628
#endif

#ifndef EPSILON
#define EPSILON 1.0e-9
#endif

void Bayes(double *src, unsigned char *label, unsigned char *priors, int niters, double *separations, int *dims, int iters_nu)
{
  int i, j, k, k1, subit;
  int z_area, y_dims;
  unsigned char *msk;
  double ll = -FLT_MAX, llr=0.0, *nu;
  
  int area = dims[0]*dims[1];
  int vol = area*dims[2];

  int Kb = 4;
  int ngauss[4] = {1,1,1,4};
  int K = 0, K2 = 0;
  
  if (iters_nu > 0) 
    nu = (double *)malloc(sizeof(double)*vol);

  // K = sum(ngauss)
  for (k1=0; k1<Kb; k1++)   K  += ngauss[k1]; 
  for (k1=0; k1<Kb-1; k1++) K2 += ngauss[k1]; 

  // lkp = [0 1 2 3 3 3 3]
  int lkp[K];
  int l = 0;
  for (k1=0; k1<Kb; k1++) {
    for (j=0; j<ngauss[k1]; j++) {
      lkp[l] = k1; 
      l++;
    }
  }

  // image maximum  
  double mx = -FLT_MAX;
  for (i=0; i<vol; i++) mx = MAX(src[i], mx);

  // initialize mean, var, mg
  double mn[K], vr[K], mg[K], mn2[K], vr2[K], mg2[K];
  for (k=0; k<K; k++) {
    mn[k] = mx * drand48();
    mg[k] = 1.0/(double)K;
    vr[k] = mx*mx + EPSILON;
  }
  
  msk = (unsigned char *)malloc(sizeof(unsigned char)*vol);
  
  int nm = 0;
  for (i=0; i<vol; i++) {
    msk[i] = ((priors[i + (3*vol)] < 254) && (src[i] != 0)) ? 1 : 0;
    nm += msk[i];
  }
  
  double mom0[K], mom1[K], mom2[K], mgm[Kb];
  double q[K], bt[Kb], b[K], qt[Kb];
  double tol1 = 1e-4;
  
  for (j=0; j<3; j++) {
    for (subit=0; subit<5; subit++) {
      double oll = ll;
      ll = llr;
      for (k=0; k<K; k++) {
        mom0[k] = 0.0;
        mom1[k] = 0.0;
        mom2[k] = 0.0;
      }
      for (k1=0; k1<Kb; k1++) mgm[k1] = 0.0;
    
      for (k=0; k<K; k++) {
        mg2[k] = mg[k];
        mn2[k] = mn[k];
        vr2[k] = vr[k];
      }
      for (i=0; i<vol; i++) {
        if (msk[i] > 0) {
          double s = EPSILON;
          for (k1=0; k1<Kb; k1++) {
            bt[k1] = (double) priors[i+(vol*k1)]; 
            for (l=0; l<ngauss[k1]; l++)      
              s += bt[k1]*mg[ngauss[k1]-1+l];
          }
          for (k1=0; k1<Kb; k1++) {
            bt[k1] /= s;
            mgm[k1] += bt[k1];
          }
          double sq = EPSILON;
          double qmax = -FLT_MAX;
          int kmax;
          for (k=0; k<K; k++) {
            q[k] = mg[k]*bt[lkp[k]]*exp(SQR(src[i]-mn[k])/(-2*vr[k]))/(SQRT2PI*sqrt(vr[k]));            
            sq += q[k];
            if (q[k] > qmax) {
              qmax = q[k];
              if (k>K2-1) kmax = 0; else kmax = k + 1;
            } 
          }
          label[i] = kmax;
          ll += log10(sq);
          for (k=0; k<K; k++) {
            double p1 = q[k]/sq; mom0[k] += p1;
            p1 *= src[i];        mom1[k] += p1;
            p1 *= src[i];        mom2[k] += p1;
          }      
        } else label[i] = 0;
      } 
      
      for (k=0; k<K; k++) {
        mg[k] = (mom0[k]+EPSILON)/(mgm[lkp[k]]+EPSILON);
        mn[k] = mom1[k]/(mom0[k]+EPSILON);
        vr[k] = (mom2[k]-SQR(mom1[k])/mom0[k]+1e6*EPSILON)/(mom0[k]+EPSILON);
        vr[k] = MAX(vr[k],EPSILON);
        // rescue previous values in case of nan
        if ((isnan(mg[k])) || (isnan(mn[k])) || (isnan(vr[k]))) {
          mg[k] = mg2[k];
          mn[k] = mn2[k];
          vr[k] = vr2[k];
          break;
        }
      }
      
      if((ll-oll)<tol1*nm) break;
      fprintf(stderr,"%g\n",ll);    
    }
    
    // only use values above the middle of the upper two cluster for nu-estimate
	double th_src = (double)((mn[Kb-2]+mn[Kb-3])/1.9);

    for (subit=0; subit<(iters_nu>0); subit++) {
      for (i = 0; i < vol; i++) {
        nu[i] = 0.0;
        // only use values above threshold where mask is defined for nu-estimate
        if ((src[i] > th_src) && (msk[i] > 0)) {
          double val_nu = src[i]/mn[label[i]-1];
          if ((isfinite(val_nu) && (val_nu < 10))) {
            nu[i] = val_nu;
          }
        }
      }
      fprintf(stderr,"Fitting splines ...\n");
      // spline estimate
      splineSmooth(nu, 0.01, 500.0, 4, separations, dims);
      
      // apply nu correction to source image
      for (i=0; i<vol; i++) {
        if (nu[i] != 0.0)
          src[i] /= nu[i];
      }
    }    
  }

  for (i=0; i<vol; i++) {
    if (src[i] != 0) {
      for (k1=0; k1<Kb; k1++) {
        bt[k1] = (double) priors[i+(vol*k1)]; 
        qt[k1] = 0.0;
      }
      double s = EPSILON;
      for (k=0; k<K; k++) { 
        b[k] = bt[lkp[k]]*mg[k];
        s += b[k];
      }
      double sq = EPSILON;
      for (k=0; k<K; k++) { 
        double p1 = exp(SQR(src[i]-mn[k])/(-2*vr[k]))/(SQRT2PI*sqrt(vr[k]));
        qt[lkp[k]] += p1*b[k]/s;
        sq += qt[lkp[k]];
      }
      double qmax = -FLT_MAX;
      int kmax;
      double psum = 0.0;
      for (k1=0; k1<3; k1++) {
        double p1 = qt[k1]/sq;
        psum += p1;
        if (p1 > qmax) { 
          qmax = p1;     
          kmax = k1 + 1;
        }   
      }
      // label only if sum of all probabilities is > 0.1
      if (psum > 0.1)
        label[i] = kmax;
      else
        label[i] = 0;
    } else label[i] = 0;
  }

  for (k=0; k<K2; k++) 
    fprintf(stderr,"%g %g\n",mn[k],sqrt(vr[k]));
    
  free(msk);    
  if (iters_nu > 0) 
    free(nu);
    
  return;
}  
