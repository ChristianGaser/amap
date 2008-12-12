#include <stdio.h>
#include <math.h>

#define MAX_NC 5
#define TH_COLOR 2

void MrfPrior(unsigned char *label, int nc, double *alpha, double *beta, int BG, int init, int *dims)
{
  int i, j, k, x, y, z;
  int fi, fj;
  long color[MAX_NC][7][7][7][7];
  long area;

  int f[MAX_NC-1], plab, zero, iBG;
  int n, z_area, y_dims;
  double XX, YY, L;
  
  area = dims[0]*dims[1];

  /* initialize configuration counts */
  for (i=0; i<nc; i++)
    for (f[0]=0; f[0]<7; f[0]++)
      for (f[1]=0; f[1]<7; f[1]++)
        for (f[2]=0; f[2]<7; f[2]++)
          for (f[3]=0; f[3]<7; f[3]++)
            color[i][f[0]][f[1]][f[2]][f[3]]=0;

  /* calculate configuration counts */
  n = 0;
  for (i=0; i<nc; i++) alpha[i] = 0.0;

  for (z=1; z<dims[2]-1; z++) {
    z_area=z*area; 
    for (y=1; y<dims[1]-1; y++) {
      y_dims = y*dims[0];
      for (x=1; x<dims[0]-1; x++) {
      
        plab = (int)label[z_area + y_dims + x];
        
        zero = plab;
        if (zero < BG) continue;
        n++;
        alpha[zero - BG] += 1.0;
        
        for (i=1; i<nc; i++) {
          f[i-1] = 0;	  
          iBG = i+BG;
          if ((int)label[z_area + y_dims + x-1] == iBG)        f[i-1]++;
          if ((int)label[z_area + y_dims + x+1] == iBG)        f[i-1]++;
          if ((int)label[z_area + ((y-1)*dims[0]) + x] == iBG) f[i-1]++;
          if ((int)label[z_area + ((y+1)*dims[0]) + x] == iBG) f[i-1]++;
          if ((int)label[((z-1)*area) + y_dims + x] == iBG)    f[i-1]++;
          if ((int)label[((z+1)*area) + y_dims + x] == iBG)    f[i-1]++;
        }
        for (i=nc; i<nc+1; i++) f[i-1] = 0;
        color[zero-BG][f[0]][f[1]][f[2]][f[3]]++;
      }
    }
  }

  /* evaluate alphas, either proportional to counts or made them equal */
  fprintf(stderr,"MRF priors: alpha ");
  for (i=0; i<nc; i++) {
    if (init == 0) alpha[i] /= n; else alpha[i] = 1.0;
    fprintf(stderr,"%3.3f ", alpha[i]);
  }

  /* compute beta */
  XX=0.0, YY=0.0;
  for (f[0]=0; f[0]<7; f[0]++)
    for (f[1]=0; f[1]<7; f[1]++)
      for (f[2]=0; f[2]<7; f[2]++)
        for (f[3]=0; f[3]<7; f[3]++)
          for (i=0; i<nc; i++)
            for (j=0; j<i; j++) {

              if (color[i][f[0]][f[1]][f[2]][f[3]] < TH_COLOR ||
                  color[j][f[0]][f[1]][f[2]][f[3]] < TH_COLOR) continue;
	      
              L = log(((double) color[i][f[0]][f[1]][f[2]][f[3]])/
                       (double) color[j][f[0]][f[1]][f[2]][f[3]]);
	      
              if (i == 0) {
                fi = 6;
                for (k=0; k<nc+1; k++) fi -= f[k];
              } else fi = f[i-1];
	      
              if (j == 0) { 
                fj = 6;
                for (k=0; k<nc+1; k++) fj -= f[k];
              } else fj = f[j-1];

              XX += (fi-fj)*(fi-fj);
              YY += L*(fi-fj);
  }
  beta[0] = XX/YY;
  fprintf(stderr,"\t beta %3.3f\n", beta[0]);
}

void MrfPrior5(unsigned char *label, int nc, double *alpha, double *beta, int BG, int init, int *dims)
{
  int i, j, k, x, y, z;
  int fi, fj;
  long color[MAX_NC][7][7][7][7][7][7];
  long area;

  int f[MAX_NC-1], plab, zero, iBG;
  int n, z_area, y_dims;
  double XX, YY, L;
  
  area = dims[0]*dims[1];

  /* initialize configuration counts */
  for (i=0; i<nc; i++)
    for (f[0]=0; f[0]<7; f[0]++)
      for (f[1]=0; f[1]<7; f[1]++)
        for (f[2]=0; f[2]<7; f[2]++)
          for (f[3]=0; f[3]<7; f[3]++)
            for (f[4]=0; f[4]<7; f[4]++)
              for (f[5]=0; f[5]<7; f[5]++)
                color[i][f[0]][f[1]][f[2]][f[3]][f[4]][f[5]]=0;

  /* calculate configuration counts */
  n = 0;
  for (i=0; i<nc; i++) alpha[i] = 0.0;

  for (z=1; z<dims[2]-1; z++) {
    z_area=z*area; 
    for (y=1; y<dims[1]-1; y++) {
      y_dims = y*dims[0];
      for (x=1; x<dims[0]-1; x++) {
      
        plab = (int)label[z_area + y_dims + x];
        
        zero = plab;
        if (zero < BG) continue;
        n++;
        alpha[zero - BG] += 1.0;
        
        for (i=1; i<nc; i++) {
          f[i-1] = 0;	  
          iBG = i+BG;
          if ((int)label[z_area + y_dims + x-1] == iBG)        f[i-1]++;
          if ((int)label[z_area + y_dims + x+1] == iBG)        f[i-1]++;
          if ((int)label[z_area + ((y-1)*dims[0]) + x] == iBG) f[i-1]++;
          if ((int)label[z_area + ((y+1)*dims[0]) + x] == iBG) f[i-1]++;
          if ((int)label[((z-1)*area) + y_dims + x] == iBG)    f[i-1]++;
          if ((int)label[((z+1)*area) + y_dims + x] == iBG)    f[i-1]++;
        }
        for (i=nc; i<nc+1; i++) f[i-1] = 0;
        color[zero-BG][f[0]][f[1]][f[2]][f[3]][f[4]][f[5]]++;
      }
    }
  }

  /* evaluate alphas, either proportional to counts or made them equal */
  fprintf(stderr,"MRF priors: alpha ");
  for (i=0; i<nc; i++) {
    if (init == 0) alpha[i] /= n; else alpha[i] = 1.0;
    fprintf(stderr,"%3.3f ", alpha[i]);
  }

  /* compute beta */
  XX=0.0, YY=0.0;
  for (f[0]=0; f[0]<7; f[0]++)
    for (f[1]=0; f[1]<7; f[1]++)
      for (f[2]=0; f[2]<7; f[2]++)
        for (f[3]=0; f[3]<7; f[3]++)
          for (f[4]=0; f[4]<7; f[4]++)
            for (f[5]=0; f[5]<7; f[5]++)
              for (i=0; i<nc; i++)
                for (j=0; j<i; j++) {

                  if (color[i][f[0]][f[1]][f[2]][f[3]][f[4]][f[5]] < TH_COLOR ||
                      color[j][f[0]][f[1]][f[2]][f[3]][f[4]][f[5]] < TH_COLOR) continue;
	      
                  L = log(((double) color[i][f[0]][f[1]][f[2]][f[3]][f[4]][f[5]])/
                           (double) color[j][f[0]][f[1]][f[2]][f[3]][f[4]][f[5]]);
	      
                  if (i == 0) {
                    fi = 6;
                    for (k=0; k<nc+1; k++) fi -= f[k];
                  } else fi = f[i-1];
	      
                  if (j == 0) { 
                    fj = 6;
                    for (k=0; k<nc+1; k++) fj -= f[k];
                  } else fj = f[j-1];

                  XX += (fi-fj)*(fi-fj);
                  YY += L*(fi-fj);
  }
  beta[0] = XX/YY;
  fprintf(stderr,"\t beta %3.3f\n", beta[0]);
}








