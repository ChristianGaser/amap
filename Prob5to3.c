
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

unsigned char * Prob5to3(unsigned char *prob, unsigned char *label0, double *mean, int BG, int *dims)
{
  int x,y,z,i,z_area,y_dims,ind,mxi;
  double val[5],w[4],*out,psum,mx, tmp_val[3];
  unsigned char *label, new_val[3];
  
  int area = dims[0]*dims[1];
  int vol = area*dims[2];
  
  // calculate weights according to differences between means
  w[0] = (mean[2]-mean[1])/(mean[2]-mean[0]);
  w[1] = 1 - w[0];
  w[2] = (mean[4]-mean[3])/(mean[4]-mean[2]);
  w[3] = 1 - w[2];
  
  label = (unsigned char*)malloc(sizeof(unsigned char)*vol);
  
  for (z = 1; z < dims[2]-1; z++) {
    z_area = z*area;
    for (y = 1; y < dims[1]-1; y++) {
      y_dims = y*dims[0];
      for (x = 1; x < dims[0]-1; x++) {
        ind = z_area + y_dims + x;
        // get original 5 probabilities
        for (i = 0; i < 5; i++) val[i] = (double)prob[(i*vol) + ind];
        if (label0[ind] >= BG) {
          // calculate 3 new probabilities using weight
          tmp_val[0] = val[0] + (w[0]*val[1]);
          tmp_val[1] = val[2] + (w[2]*val[3]) + (w[1]*val[1]);
          tmp_val[2] = val[4] + (w[3]*val[3]);

          // sum of all probabilities
          psum = tmp_val[0] + tmp_val[1] + tmp_val[2];

          new_val[0] = (unsigned char) round(255.0*tmp_val[0]/psum);
          new_val[1] = (unsigned char) round(255.0*tmp_val[1]/psum);
          new_val[2] = (unsigned char) round(255.0*tmp_val[2]/psum);

          prob[ind] = new_val[0];
          prob[vol + ind] = new_val[1];
          prob[(2*vol) + ind] = new_val[2];
        
          // get new label
          mx = -1e15;
          for (i = 0; i <= 2; i++)
            if (new_val[i] > mx) {
              mx = new_val[i];
              mxi = i;
            }
          label[ind] = mxi+BG;
        }
      }
    }
  }
  
  return(label);
}