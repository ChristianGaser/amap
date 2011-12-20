/*--------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1996, Alex P. Zijdenbos, 
              McConnell Brain Imaging Centre,
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- 
$RCSfile: MatrixSupport.cc,v $
$Revision: 1.4 $
$Author: bert $
$Date: 2004-01-29 19:42:19 $
$State: Exp $
--------------------------------------------------------------------------*/
#include "config.h"

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif 

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>		// (bert) changed from iostream.h
using namespace std;		// (bert) added
#include "dcomplex.h"
#include "MatrixSupport.h"


/*
** FFT and FHT routines
**  Copyright 1988, 1993; Ron Mayer
**  
**  fht(fz,n);
**      Does a hartley transform of "n" points in the array "fz".
**  fft(n,real,imag)
**      Does a fourier transform of "n" points of the "real" and
**      "imag" arrays.
**  ifft(n,real,imag)
**      Does an inverse fourier transform of "n" points of the "real"
**      and "imag" arrays.
**  realfft(n,real)
**      Does a real-valued fourier transform of "n" points of the
**      "real" and "imag" arrays.  The real part of the transform ends
**      up in the first half of the array and the imaginary part of the
**      transform ends up in the second half of the array.
**  realifft(n,real)
**      The inverse of the realfft() routine above.
**      
**      
** NOTE: This routine uses at least 2 patented algorithms, and may be
**       under the restrictions of a bunch of different organizations.
**       Although I wrote it completely myself; it is kind of a derivative
**       of a routine I once authored and released under the GPL, so it
**       may fall under the free software foundation's restrictions;
**       it was worked on as a Stanford Univ project, so they claim
**       some rights to it; it was further optimized at work here, so
**       I think this company claims parts of it.  The patents are
**       held by R. Bracewell (the FHT algorithm) and O. Buneman (the
**       trig generator), both at Stanford Univ.
**       If it were up to me, I'd say go do whatever you want with it;
**       but it would be polite to give credit to the following people
**       if you use this anywhere:
**           Euler     - probable inventor of the fourier transform.
**           Gauss     - probable inventor of the FFT.
**           Hartley   - probable inventor of the hartley transform.
**           Buneman   - for a really cool trig generator
**           Mayer(me) - for authoring this particular version and
**                       including all the optimizations in one package.
**       Thanks,
**       Ron Mayer; mayer@acuson.com
**
*/

#define GOOD_TRIG
char fht_version[] = "Brcwl-Hrtly-Ron-dbld";

#define SQRT2_2   0.70710678118654752440084436210484
#define SQRT2   2*0.70710678118654752440084436210484
void
fht(double *fz, int n)
{
  // double c1,s1,s2,c2,s3,c3,s4,c4;
  // double f0,g0,f1,g1,f2,g2,f3,g3;
 int i,k,k1,k2,k3,k4,kx;
 double *fi,*fn,*gi;
 TRIG_VARS;

 for (k1=1,k2=0;k1<n;k1++)
    {
     double a;
     for (k=n>>1; (!((k2^=k)&k)); k>>=1);
     if (k1>k2)
	{
	     a=fz[k1];fz[k1]=fz[k2];fz[k2]=a;
	}
    }
 for ( k=0 ; (1<<k)<n ; k++ );
 k  &= 1;
 if (k==0)
    {
	 for (fi=fz,fn=fz+n;fi<fn;fi+=4)
	    {
	     double f0,f1,f2,f3;
	     f1     = fi[0 ]-fi[1 ];
	     f0     = fi[0 ]+fi[1 ];
	     f3     = fi[2 ]-fi[3 ];
	     f2     = fi[2 ]+fi[3 ];
	     fi[2 ] = (f0-f2);	
	     fi[0 ] = (f0+f2);
	     fi[3 ] = (f1-f3);	
	     fi[1 ] = (f1+f3);
	    }
    }
 else
    {
	 for (fi=fz,fn=fz+n,gi=fi+1;fi<fn;fi+=8,gi+=8)
	    {
	     double s1,c1,s2,c2,s3,c3,s4,c4,g0,f0,f1,g1,f2,g2,f3,g3;
	     c1     = fi[0 ] - gi[0 ];
	     s1     = fi[0 ] + gi[0 ];
	     c2     = fi[2 ] - gi[2 ];
	     s2     = fi[2 ] + gi[2 ];
	     c3     = fi[4 ] - gi[4 ];
	     s3     = fi[4 ] + gi[4 ];
	     c4     = fi[6 ] - gi[6 ];
	     s4     = fi[6 ] + gi[6 ];
	     f1     = (s1 - s2);	
	     f0     = (s1 + s2);
	     g1     = (c1 - c2);	
	     g0     = (c1 + c2);
	     f3     = (s3 - s4);	
	     f2     = (s3 + s4);
	     g3     = SQRT2*c4;		
	     g2     = SQRT2*c3;
	     fi[4 ] = f0 - f2;
	     fi[0 ] = f0 + f2;
	     fi[6 ] = f1 - f3;
	     fi[2 ] = f1 + f3;
	     gi[4 ] = g0 - g2;
	     gi[0 ] = g0 + g2;
	     gi[6 ] = g1 - g3;
	     gi[2 ] = g1 + g3;
	    }
    }
 if (n<16) return;

 do
    {
     double s1,c1;
     k  += 2;
     k1  = 1  << k;
     k2  = k1 << 1;
     k4  = k2 << 1;
     k3  = k2 + k1;
     kx  = k1 >> 1;
	 fi  = fz;
	 gi  = fi + kx;
	 fn  = fz + n;
	 do
	    {
	     double g0,f0,f1,g1,f2,g2,f3,g3;
	     f1      = fi[0 ] - fi[k1];
	     f0      = fi[0 ] + fi[k1];
	     f3      = fi[k2] - fi[k3];
	     f2      = fi[k2] + fi[k3];
	     fi[k2]  = f0	  - f2;
	     fi[0 ]  = f0	  + f2;
	     fi[k3]  = f1	  - f3;
	     fi[k1]  = f1	  + f3;
	     g1      = gi[0 ] - gi[k1];
	     g0      = gi[0 ] + gi[k1];
	     g3      = SQRT2  * gi[k3];
	     g2      = SQRT2  * gi[k2];
	     gi[k2]  = g0	  - g2;
	     gi[0 ]  = g0	  + g2;
	     gi[k3]  = g1	  - g3;
	     gi[k1]  = g1	  + g3;
	     gi     += k4;
	     fi     += k4;
	    } while (fi<fn);
     TRIG_INIT(k,c1,s1);
     for (i=1;i<kx;i++)
        {
	 double c2,s2;
         TRIG_NEXT(k,c1,s1);
         c2 = c1*c1 - s1*s1;
         s2 = 2*(c1*s1);
	     fn = fz + n;
	     fi = fz +i;
	     gi = fz +k1-i;
	     do
		{
		 double a,b,g0,f0,f1,g1,f2,g2,f3,g3;
		 b       = s2*fi[k1] - c2*gi[k1];
		 a       = c2*fi[k1] + s2*gi[k1];
		 f1      = fi[0 ]    - a;
		 f0      = fi[0 ]    + a;
		 g1      = gi[0 ]    - b;
		 g0      = gi[0 ]    + b;
		 b       = s2*fi[k3] - c2*gi[k3];
		 a       = c2*fi[k3] + s2*gi[k3];
		 f3      = fi[k2]    - a;
		 f2      = fi[k2]    + a;
		 g3      = gi[k2]    - b;
		 g2      = gi[k2]    + b;
		 b       = s1*f2     - c1*g3;
		 a       = c1*f2     + s1*g3;
		 fi[k2]  = f0        - a;
		 fi[0 ]  = f0        + a;
		 gi[k3]  = g1        - b;
		 gi[k1]  = g1        + b;
		 b       = c1*g2     - s1*f3;
		 a       = s1*g2     + c1*f3;
		 gi[k2]  = g0        - a;
		 gi[0 ]  = g0        + a;
		 fi[k3]  = f1        - b;
		 fi[k1]  = f1        + b;
		 gi     += k4;
		 fi     += k4;
		} while (fi<fn);
        }
     TRIG_RESET(k,c1,s1);
    } while (k4<n);
}

void
ifft(int n, double *real, double *imag)
{
 double a,b,c,d;
 double q,r,s,t;
 int i,j,k;
 fht(real,n);
 fht(imag,n);
 for (i=1,j=n-1,k=n/2;i<k;i++,j--) {
  a = real[i]; b = real[j];  q=a+b; r=a-b;
  c = imag[i]; d = imag[j];  s=c+d; t=c-d;
  imag[i] = (s+r)*0.5;  imag[j] = (s-r)*0.5;
  real[i] = (q-t)*0.5;  real[j] = (q+t)*0.5;
 }
}

void
realfft(int n, double *real)
{
 double a,b;
 int i,j,k;
 fht(real,n);
 for (i=1,j=n-1,k=n/2;i<k;i++,j--) {
  a = real[i];
  b = real[j];
  real[j] = (a-b)*0.5;
  real[i] = (a+b)*0.5;
 }
}

void
fft(int n, double *real, double *imag)
{
 double a,b,c,d;
 double q,r,s,t;
 int i,j,k;
 for (i=1,j=n-1,k=n/2;i<k;i++,j--) {
  a = real[i]; b = real[j];  q=a+b; r=a-b;
  c = imag[i]; d = imag[j];  s=c+d; t=c-d;
  real[i] = (q+t)*.5; real[j] = (q-t)*.5;
  imag[i] = (s-r)*.5; imag[j] = (s+r)*.5;
 }
 fht(real,n);
 fht(imag,n);
}

void
realifft(int n, double *real)
{
 double a,b;
 int i,j,k;
 for (i=1,j=n-1,k=n/2;i<k;i++,j--) {
  a = real[i];
  b = real[j];
  real[j] = (a-b);
  real[i] = (a+b);
 }
 fht(real,n);
}

/*
** End of Ron Mayer's code
 */

/* This function is removed because it exists in trivial.h called from FileIO.h
double 
gauss(double mean, double std_dev)
{
   double v1 = 2.0*drand48() - 1.0;
   double v2 = 2.0*drand48() - 1.0;
   double s = v1*v1 + v2*v2;
   
   while (s >= 1.0) {
      v1 = 2.0*drand48() - 1.0;
      v2 = 2.0*drand48() - 1.0;
      s  = v1*v1 + v2*v2;
   }
   
   return mean + std_dev*v1*(::sqrt(-2.0*log(s)/s));
}
*/

void 
jacobi(double **a, unsigned n, double *d, double **v)
{
  unsigned nrot;
   int    j,iq,ip,i;
   double tresh,theta,tau,t,sm,s,h,g,c;
   double *b, *z;
   
   b = (double *) malloc((n+1)*sizeof(double));
   z = (double *) malloc((n+1)*sizeof(double));
   
   for( ip = 1 ; ip <= n ; ip++)  {
      for( iq = 1 ; iq <= n ; iq++)  {
	 v[ip][iq] = 0.0;
      }
      v[ip][ip] = 1.0;
   }
   for( ip = 1 ; ip <= n ; ip++)  {
      b[ip] = a[ip][ip];
      d[ip] = b[ip];
      z[ip] = 0.0;
   }
   nrot = 0;
   for( i = 1 ; i <= 50 ; i++)  {
      sm = 0.0;
      for( ip = 1 ; ip <= n-1 ; ip++ )  {
	 for( iq = ip+1 ; iq <= n ; iq++ )  {
	    sm = sm+fabs(a[ip][iq]);
	 }
      }
      if (sm == 0.0)
	 goto position99;
      if (i < 4)
      { tresh = 0.2*sm/(n*n);}
      else
      { tresh = 0.0;}
      for( ip = 1 ; ip <= n-1 ; ip++)  {
	 for( iq = ip+1 ; iq <= n ; iq++)  {
	    g = 100.0*fabs(a[ip][iq]);
	    if ( (i > 4) &&
		 ( (fabs(d[ip])+g) == fabs(d[ip]) ) &&
		 ( (fabs(d[iq])+g) == fabs(d[iq]) ) )
	       a[ip][iq] = 0.0;
	    else
	       if (fabs(a[ip][iq]) > tresh)  {
		  h = d[iq]-d[ip];
		  if ((fabs(h)+g) == fabs(h))
		     t = a[ip][iq]/h;
		  else {
		     theta = 0.5*h/a[ip][iq];
		     t = 1.0/(fabs(theta)+::sqrt(1.0+theta*theta));
		     if (theta < 0.0)
			t = -t;
		  }
		  c = 1.0/(::sqrt(1+ t*t ));
		  s = t*c;
		  tau = s/(1.0+c);
		  h = t*a[ip][iq];
		  z[ip] = z[ip]-h;
		  z[iq] = z[iq]+h;
		  d[ip] = d[ip]-h;
		  d[iq] = d[iq]+h;
		  a[ip][iq] = 0.0;
		  for( j = 1 ; j <= ip-1 ; j++) {
		     g = a[j][ip];
		     h = a[j][iq];
		     a[j][ip] = g-s*(h+g*tau);
		     a[j][iq] = h+s*(g-h*tau);
		  }
		  for( j = ip+1 ; j <= iq-1 ; j++)  {
		     g = a[ip][j];
		     h = a[j][iq];
		     a[ip][j] = g-s*(h+g*tau);
		     a[j][iq] = h+s*(g-h*tau);
		  }
		  for( j = iq+1 ; j <= n ; j++)  {
		     g = a[ip][j];
		     h = a[iq][j];
		     a[ip][j] = g-s*(h+g*tau);
		     a[iq][j] = h+s*(g-h*tau);
		  }
		  for( j = 1 ; j <= n ; j++)  {
		     g = v[j][ip];
		     h = v[j][iq];
		     v[j][ip] = g-s*(h+g*tau);
		     v[j][iq] = h+s*(g-h*tau);
		  }
		  nrot = nrot+1;
	       }
	 }
      }
   }
   for( ip = 1 ; ip <= n ; ip++)  {
      b[ip] = b[ip]+z[ip];
      d[ip] = b[ip];
      z[ip] = 0.0;
   }
   printf("pause in routine JACOBI\n");
   printf("50 iterations should not happen\n");
   
   getchar();
   
 position99: 
   free((char *)b);
   free((char *)z);
   return; 
}


//fft functions:
//--------------
/*****************************************************************************
*
*  void fft_basic(cbuffer,sintab,buflog,direct)
     *          This function should not be called directly.  It is called by
     *          the fft(cbuffer,len) or the ifft(cbuffer,len).
     *
     *  Input:  Takes input directly from fft(cbuffer,len) or ifft(cbuffer,len).
     *  Output: No output.
     *
     p******************************************************************************/
/* Do not touch this!!! */
#ifdef USE_COMPMAT
void fft_basic(dcomplex *cbuffer,double *sintab, int buflog, int direct)
{
   int     bufflen,halflen,quadlen;
   register int blknum,blkctr,offset,sinarg;
   double   wsin,wcos,remult,immult;
   register dcomplex *bufpt0,*bufpt1;
   
   bufflen = 1 << buflog;
   halflen = bufflen >> 1;
   quadlen = halflen >> 1;
   for(blknum = bufflen; (blkctr = --blknum) >= 0; ) {
      for(offset = 0,sinarg = buflog; --sinarg >= 0; ) {
	 offset <<= 1;
	 offset += blkctr & 1;
	 blkctr >>= 1;
      }
      if(blknum < offset) {
	 bufpt0 = cbuffer + blknum;
	 bufpt1 = cbuffer + offset;
	 remult = real(*bufpt0);
	 immult = imag(*bufpt0);
	 *bufpt0 = dcomplex(real(*bufpt1),imag(*bufpt1));
	 *bufpt1 = dcomplex(remult,immult);
      }
   }
   offset = 1;
   blknum = halflen;
   while(--buflog >= 0) {
      sinarg = halflen;
      bufpt0 = cbuffer + (bufflen + offset);
      bufpt1 = bufpt0 + offset;
      offset <<= 1;
      while((sinarg -= blknum) >= 0) {
	 if(sinarg <= quadlen) {
	    wsin = sintab[sinarg];
	    wcos = sintab[quadlen - sinarg];
	 }
	 else {
	    wsin = sintab[halflen - sinarg];
	    wcos = -(sintab[sinarg - quadlen]);
	 }
	 if(direct == 0)
	    wsin = (-wsin);
	 blkctr = blknum;
	 bufpt0 -= (bufflen + 1);
	 bufpt1 -= (bufflen + 1);
	 while(--blkctr >= 0) {
	    remult = wcos * real(*bufpt1) - wsin * imag(*bufpt1);
	    immult = wcos * imag(*bufpt1) + wsin * real(*bufpt1);
	    *bufpt1 = *bufpt0 - dcomplex(remult,immult);
	    *bufpt0 = *bufpt0 + dcomplex(remult,immult);
	    bufpt0 += offset; 
	    bufpt1 += offset;
	 }
      }
      blknum >>= 1;
   }
}
#endif

///
#ifdef USE_FCOMPMAT
void fft_basic_float(fcomplex *cbuffer,float *sintab, int buflog, int direct)
{
   int     bufflen,halflen,quadlen;
   register int blknum,blkctr,offset,sinarg;
   float   wsin,wcos,remult,immult;
   register fcomplex *bufpt0,*bufpt1;
   
   bufflen = 1 << buflog;
   halflen = bufflen >> 1;
   quadlen = halflen >> 1;
   for(blknum = bufflen; (blkctr = --blknum) >= 0; ) {
      for(offset = 0,sinarg = buflog; --sinarg >= 0; ) {
	 offset <<= 1;
	 offset += blkctr & 1;
	 blkctr >>= 1;
      }
      if(blknum < offset) {
	 bufpt0 = cbuffer + blknum;
	 bufpt1 = cbuffer + offset;
	 remult = real(*bufpt0);
	 immult = imag(*bufpt0);
	 *bufpt0 = fcomplex(real(*bufpt1),imag(*bufpt1));
	 *bufpt1 = fcomplex(remult,immult);
      }
   }
   offset = 1;
   blknum = halflen;
   while(--buflog >= 0) {
      sinarg = halflen;
      bufpt0 = cbuffer + (bufflen + offset);
      bufpt1 = bufpt0 + offset;
      offset <<= 1;
      while((sinarg -= blknum) >= 0) {
	 if(sinarg <= quadlen) {
	    wsin = sintab[sinarg];
	    wcos = sintab[quadlen - sinarg];
	 }
	 else {
	    wsin = sintab[halflen - sinarg];
	    wcos = -(sintab[sinarg - quadlen]);
	 }
	 if(direct == 0)
	    wsin = (-wsin);
	 blkctr = blknum;
	 bufpt0 -= (bufflen + 1);
	 bufpt1 -= (bufflen + 1);
	 while(--blkctr >= 0) {
	    remult = wcos * real(*bufpt1) - wsin * imag(*bufpt1);
	    immult = wcos * imag(*bufpt1) + wsin * real(*bufpt1);
	    *bufpt1 = *bufpt0 - fcomplex(remult,immult);
	    *bufpt0 = *bufpt0 + fcomplex(remult,immult);
	    bufpt0 += offset; 
	    bufpt1 += offset;
	 }
      }
      blknum >>= 1;
   }
}
#endif
// End of C code


void 
inferDimensions(unsigned long nElements, unsigned& nrows, unsigned& ncols)
{
  // Both nrows and ncols or only ncols supplied; (re)calculate nrows
  if (ncols) {
    nrows = nElements/ncols;
    return;
  }

  // Only nrows supplied; calculate ncols
  if (nrows) {
    ncols = nElements/nrows;
    return;
  }

  // Neither nrows nor ncols supplied; start with a square matrix, and 
  // reduce nrows until a match is found.
  nrows = ncols = unsigned(::sqrt(double(nElements)));

  while (nrows * ncols != nElements) {
    nrows--;
    ncols = nElements/nrows;
  }
}

void 
inferDimensions(unsigned long nElements, unsigned& nslis, unsigned& nrows, 
		unsigned& ncols)
{
  if (nslis) {
    inferDimensions(nElements/nslis, nrows, ncols);
    return;
  }

  if (nrows) {
    inferDimensions(nElements/nrows, nslis, ncols);
    return;
  }

  if (ncols) {
    inferDimensions(nElements/ncols, nslis, nrows);
    return;
  }

  // No dimension supplied; start with a cubical matrix, and 
  // reduce slis/nrows until a match is found.
  nslis = nrows = ncols = unsigned(pow(double(nElements), 1.0/3));

  while (nslis * nrows * ncols != nElements) {
    nslis--;
    inferDimensions(nElements/nslis, nrows, ncols);
  }
}

