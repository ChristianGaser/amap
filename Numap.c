/*
 * Bmap.C
 *
 * Biased Maximum Apriori Segmentation
 *
 * Jagath Rajapakse and Frithjof Kruggel
 * Segmentation of MR Images with Intensity Inhomogeneities
 * Image & Vision Computing, (in press)
 *
 * comments to raja@cns.mpg.de 
 *
 */

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
#define SQRT2PI 2.506628

void MrfPrior(unsigned char *label, int nc, double *alpha, double *beta, int BG, int init, int *dims);


void Numap(double *src, unsigned char *label, unsigned char *prob, double *mean, int nc, int BG, int niters, int nflips, double *separations, int *dims)
{
  int i,j,index,x,y,z,iters;
  int z_area, y_dims;
  long histo[65536];
  double d, val;
  double beta[1],alpha[nc],s[nc],ss[nc];
  double var[nc],lvar[nc],p[nc],value;
  double mn_thresh, mx_thresh;
  
  int area = dims[0]*dims[1];
  int vol = area*dims[2];

  // initialize prior parameters
  MrfPrior(label, nc, alpha, beta, BG, 0, dims);

  double min_src = FLT_MAX, max_src = -FLT_MAX;
  for (i=0; i<vol; i++) {
    min_src = MIN(src[i], min_src);
    max_src = MAX(src[i], max_src);
  }

  // build histogram
  for (i = 0; i < 65536; i++) histo[i] = 0;
  for (i=0; i<vol; i++) {
    if ((int)label[i] < BG) continue;
    histo[(int)round(65535.0*(src[i]-min_src)/(max_src-min_src))]++;
  }

  // find values <1% and > 99% quartile
  long cumsum[65536];
  cumsum[0] = histo[0];
  for (i = 1; i < 65536; i++) cumsum[i] = cumsum[i-1] + histo[i];
  for (i = 0; i < 65536; i++) cumsum[i] = (long) round(1000.0*(double)cumsum[i]/(double)cumsum[65535]);
  for (i = 0; i < 65536; i++) if (cumsum[i] >= 10) break;
  mn_thresh = (double)i/65535.0*(max_src-min_src);
  for (i = 65535; i > 0; i--) if (cumsum[i] <= 990) break;
  mx_thresh = (double)i/65535.0*(max_src-min_src);

  // initialize means
  for (j=0; j<nc; j++) s[j]=ss[j]=0.0;
  for (i=0; i<vol; i++) {
    int lab = (int) label[i];
    if (lab < BG) continue;
    s[lab-BG] += 1.0; ss[lab-BG] += src[i];
  }
  for (j=0; j<nc; j++)
    mean[j] = s[j] > 0.0 ? ss[j]/s[j]: 0.0;

  // intitialize standard deviations
  for (j=0; j<nc; j++) ss[j] = 0.0;
  for (i=0; i<vol; i++) {
    int lab = (int) label[i];
    if (lab < BG) continue;
    ss[lab-BG] += SQR(src[i]-mean[lab-BG]);
  }
  for (j=0; j<nc; j++) 
    var[j] = s[j]>1.0 ? ss[j]/(s[j]-1.0): 1.0;
  
  // set new variables to speed up
  for (j=0; j<nc; j++) {
    lvar[j] = var[j] > 0.0 ? 0.5*log(var[j]) - log(alpha[j]): - log(alpha[j]);
    var[j] = 0.5/var[j];
  }

  // iterative condition modes
  for (iters=0; iters<=niters; iters++)  {
    int flips = 0;

    // loop over image voxels
    for (z=1; z<dims[2]-1; z++) 
      for (y=1; y<dims[1]-1; y++)
        for (x=1; x<dims[0]-1; x++)  {
	  
          index = x + y*dims[0] + z*area;

          int lab = (int) label[index];
          if (lab < BG) continue;
	  
          // loop over all classes
          int xi=BG; double dmin = FLT_MAX;
          double psum = 0.0;
          
          for (j=0; j<nc; j++) {

            // find the number of first order neighbors in the class
            int first=0;
            if (label[index-1] == j+BG) first++;
            if (label[index+1] == j+BG) first++;
            if (label[index-dims[0]] == j+BG) first++;
            if (label[index+dims[0]] == j+BG) first++;
            if (label[index-area] == j+BG) first++;
            if (label[index+area] == j+BG) first++;
	    
            d = SQR(src[index]-mean[j])*var[j]+lvar[j]-beta[0]*first;
	    
            p[j] = exp(-d)/(SQRT2PI*sqrt(var[j]));
            psum += p[j];
	    
            if (d < dmin) {xi = j; dmin = d;}
          }
          for (j=0; j<nc; j++)
            prob[(vol*j) + index] = (unsigned char)round(255*p[j]/psum);

          // if the class has changed increment flips and change the label
          if (xi+BG != lab) {flips++;  label[index] = (unsigned char) (xi+BG); }
	  
        }
        
    fprintf(stderr,"iters:%2d flips:%6d\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",iters, flips);

    if (flips <= nflips) break;


    // find tissue means
    for (j=0; j<nc; j++) {s[j]=ss[j]=0.0;}
      for (i=0; i<vol; i++) {
        int lab = (int) label[i];
        if (lab < BG) continue;
        val = src[i];
        // exclude values out of quartile 1/99%
        if ((val<mn_thresh) || (val>mx_thresh)) continue;
        s[lab-BG]  += 1;
        ss[lab-BG] += val;
      }
    for (j=0; j<nc; j++) {mean[j] = s[j]>0.0 ? ss[j]/s[j]: 0.0;}

    // find tissue variances
    for (j=0; j<nc; j++) {s[j]=ss[j]=0.0;}
    for (i=0; i<vol; i++) {
      int lab = (int) label[i];
      if (lab < BG) continue;
      val = src[i];
      // exclude values out of quartile 1/99%
      if ((val<mn_thresh) || (val>mx_thresh)) continue;
      s[lab-BG] += 1.0;
      ss[lab-BG] += SQR(val-mean[lab-BG]);
    }
    for (j=0; j<nc; j++) {
      var[j] = s[j]>0.0 ? ss[j]/s[j]: 1.0;
      lvar[j] = var[j] > 0.0 ? 0.5*log(var[j]) - log(alpha[j]): -log(alpha[j]);
      var[j] = 0.5/var[j];
    }        
  }
      
  fprintf(stderr,"\nFinal means*vars: "); 
  for (i=0; i<nc; i++) fprintf(stderr,"%3.2f*%3.2f  ",mean[i],sqrt(0.5/var[i])); 
  fprintf(stderr,"\n"); 

  return;
}  
