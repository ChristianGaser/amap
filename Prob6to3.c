
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

unsigned char * Prob6to3(double *src, unsigned char *prob, unsigned char *label0, double *mean, int BG, int *dims)
{
  int x,y,z,i,z_area,y_dims,ind,mxi;
  double w, mx;
  unsigned char *label, new_val[3];
  
  int area = dims[0]*dims[1];
  int vol = area*dims[2];
  
  label = (unsigned char*)malloc(sizeof(unsigned char)*vol);
  
  for (z = 1; z < dims[2]-1; z++) {
    z_area = z*area;
    for (y = 1; y < dims[1]-1; y++) {
      y_dims = y*dims[0];
      for (x = 1; x < dims[0]-1; x++) {
        ind = z_area + y_dims + x;

          switch(label0[ind]) {
            case 0: /* BG */
              new_val[0] = 0;
              new_val[1] = 0;
              new_val[2] = 0;
              break;
            case 1: /* CSFBG */
              w = (src[ind])/(mean[1]);
              if(w > 1.0) w = 1.0; if(w < 0.0) w = 0.0;
              new_val[0] = (unsigned char) round(255.0*w);
              new_val[1] = 0;
              new_val[2] = 0;
              break;
            case 2: /* CSF */
              new_val[0] = 255;
              new_val[1] = 0;
              new_val[2] = 0;
              break;
            case 3: /* GMCSF */
              w = (src[ind] - mean[1])/(mean[3]-mean[1]);
              if(w > 1.0) w = 1.0; if(w < 0.0) w = 0.0;
              new_val[0] = (unsigned char) round(255.0*(1-w));
              new_val[1] = (unsigned char) round(255.0*w);
              new_val[2] = 0;
              break;
            case 4: /* GM */
              new_val[0] = 0;
              new_val[1] = 255;
              new_val[2] = 0;
              break;
            case 5: /*WMGM */
              w = (src[ind] - mean[3])/(mean[5]-mean[3]);
              if(w > 1.0) w = 1.0; if(w < 0.0) w = 0.0;
              new_val[0] = 0;
              new_val[1] = (unsigned char) round(255.0*(1-w));
              new_val[2] = (unsigned char) round(255.0*w);
              break;
            case 6: /* WM */
              new_val[0] = 0;
              new_val[1] = 0;
              new_val[2] = 255;
              break;
          }

          prob[ind] = new_val[0];
          prob[vol + ind] = new_val[1];
          prob[(2*vol) + ind] = new_val[2];
        
          // get new label
          mx = -1e15;
          for (i = 0; i < 3; i++)
            if (new_val[i] > mx) {
              mx = new_val[i];
              mxi = i;
            }
          label[ind] = mxi+BG;
      }
    }
  }
  
  return(label);
}