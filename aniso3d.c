#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>

#define MAX(A,B) (((A)>(B)) ? (A) : (B) )
#define MIN(A,B) (((A)<(B)) ? (A) : (B) )

/* Psi functions */
double tukeyPsi(double x, double sigma)
{
double x2;
if (x<(-sigma))
	return 0.0;

if (x>sigma)
	return 0.0;

x2 = (x/sigma);
x2 *= x2;
x2 = 1-x2;
x2 *= x2;
return (x*x2);
}

/* calculate standard deviation */
double stddev(
   double *vol, int *dims)
{
   int		x, y, z, iter;
   int		nvoxels;
   double	value, sum = 0.0, sum2 = 0.0;
   
   nvoxels = dims[0]*dims[1]*dims[2];
   
   for (z=0;z<dims[2];z++) {
     for (y=0;y<dims[1];y++) {
       for (x=0;x<dims[0];x++) {
          value = vol[z*dims[0]*dims[1] + y*dims[0] + x];
          sum  += value;
          sum2 += value*value;
       }
     }
   }
   return(sqrt((sum2-(sum*sum)/nvoxels)/(nvoxels-1)));
}


void aniso3d(double *vol, int *dims, double sigma, int Niter, double lambda)
{
	int x, y, z, iter, index; /* indexes and sizes */
	int Zind, ZindU, ZindD, Xind, XindL, XindR;
	double psi, std;
	double *prevvol;
	int Nbytes;
	int Ymax, Xmax, Zmax;
	int indexYS, indexYN, indexXL, indexXR, indexZU, indexZD;

	Nbytes = dims[0]*dims[1]*dims[2]*sizeof(double);
	Ymax = dims[0] - 1;
	Xmax = dims[1] - 1;
	Zmax = dims[2] - 1;
	prevvol = (double *)malloc(Nbytes);

  lambda = lambda/6.0; /* to do exactly the same as the m-file */
  std = stddev(vol, dims);
  sigma *= std;

	for (iter=0; iter<Niter; iter++) {
		/* copying the current vol to the previous vol */
		memcpy(prevvol, vol, Nbytes);
		for (z=0;z<dims[2];z++) {
			Zind = z*dims[0]*dims[1];
			ZindU = MAX(0, z-1)*dims[0]*dims[1];
			ZindD = MIN(Zmax, z+1)*dims[0]*dims[1];
			for (x=0;x<dims[1];x++) {
				Xind = x*dims[0];
				XindL = MAX(0, x-1)*dims[0];
				XindR = MIN(Xmax, x+1)*dims[0];
				for (y=0;y<dims[0];y++) {
					index = Zind + Xind + y;
					indexYS = Zind + Xind + MIN(Ymax, y+1);
					indexYN = Zind + Xind + MAX(0, y-1);
					indexXL = Zind + XindL + y;
					indexXR = Zind + XindR + y;
					indexZU = ZindU + Xind + y;
					indexZD = ZindD + Xind + y;
					psi = tukeyPsi(prevvol[indexXL]-prevvol[index], sigma) +
						tukeyPsi(prevvol[indexXR]-prevvol[index], sigma) +
						tukeyPsi(prevvol[indexZU]-prevvol[index], sigma) +
						tukeyPsi(prevvol[indexZD]-prevvol[index], sigma) +
						tukeyPsi(prevvol[indexYN]-prevvol[index], sigma) +
						tukeyPsi(prevvol[indexYS]-prevvol[index], sigma);
					vol[index] = prevvol[index] + lambda * psi;
				} /* end y loop */
			}  /* end x loop */
		}   /* end z loop */
	}    /* end main iterations loop */
	free(prevvol);
}
