/*
 * Christian Gaser
 * $Id: vollib.c 216 2020-09-18 07:46:33Z gaser $ 
 *
 * This code is a substantially modified version of spm_conv_vol.c 
 * from J. Ashburner
 */

#include "vollib.h"

static void 
convxy(double out[], int xdim, int ydim, double filtx[], double filty[], int fxdim, int fydim, int xoff, int yoff, double buff[])
{
  int x,y,k;
  for(y=0; y<ydim; y++)
  {
    for(x=0; x<xdim; x++)
    {
      buff[x] = out[x+y*xdim];
      if (!isfinite(buff[x]))
        buff[x] = 0.0;
    }
    for(x=0; x<xdim; x++)
    {
      double sum1 = 0.0;
      int fstart, fend;
      fstart = ((x-xoff >= xdim) ? x-xdim-xoff+1 : 0);
      fend = ((x-(xoff+fxdim) < 0) ? x-xoff+1 : fxdim);

      for(k=fstart; k<fend; k++)
        sum1 += buff[x-xoff-k]*filtx[k];
      out[x+y*xdim] = sum1;
    }
  }
  for(x=0; x<xdim; x++)
  {
    for(y=0; y<ydim; y++)
      buff[y] = out[x+y*xdim];

    for(y=0; y<ydim; y++)
    {
      double sum1 = 0.0;
      int fstart, fend;
      fstart = ((y-yoff >= ydim) ? y-ydim-yoff+1 : 0);
      fend = ((y-(yoff+fydim) < 0) ? y-yoff+1 : fydim);

      for(k=fstart; k<fend; k++)
        sum1 += buff[y-yoff-k]*filty[k];
      out[y*xdim+x] = sum1;
    }
  }
}

static void 
convxy_float(float out[], int xdim, int ydim, double filtx[], double filty[], int fxdim, int fydim, int xoff, int yoff, float buff[])
{
  int x,y,k;
  for(y=0; y<ydim; y++)
  {
    for(x=0; x<xdim; x++)
    {
      buff[x] = out[x+y*xdim];
      if (!isfinite(buff[x]))
        buff[x] = 0.0;
    }
    for(x=0; x<xdim; x++)
    {
      double sum1 = 0.0;
      int fstart, fend;
      fstart = ((x-xoff >= xdim) ? x-xdim-xoff+1 : 0);
      fend = ((x-(xoff+fxdim) < 0) ? x-xoff+1 : fxdim);

      for(k=fstart; k<fend; k++)
        sum1 += (double)buff[x-xoff-k]*filtx[k];
      out[x+y*xdim] = (float)sum1;
    }
  }
  for(x=0; x<xdim; x++)
  {
    for(y=0; y<ydim; y++)
      buff[y] = out[x+y*xdim];

    for(y=0; y<ydim; y++)
    {
      double sum1 = 0.0;
      int fstart, fend;
      fstart = ((y-yoff >= ydim) ? y-ydim-yoff+1 : 0);
      fend = ((y-(yoff+fydim) < 0) ? y-yoff+1 : fydim);

      for(k=fstart; k<fend; k++)
        sum1 += (double)buff[y-yoff-k]*filty[k];
      out[y*xdim+x] = (float)sum1;
    }
  }
}

int 
convxyz_double(double *iVol, double filtx[], double filty[], double filtz[],
  int fxdim, int fydim, int fzdim, int xoff, int yoff, int zoff,
  double *oVol, int dims[3])
{
  double *tmp, *buff, **sortedv, *obuf;
  int xy, z, y, x, k, fstart, fend, startz, endz;
  int xdim, ydim, zdim;

  xdim = dims[0];
  ydim = dims[1];
  zdim = dims[2];

  tmp = (double *)malloc(sizeof(double)*xdim*ydim*fzdim);
  buff = (double *)malloc(sizeof(double)*((ydim>xdim) ? ydim : xdim));
  sortedv = (double **)malloc(sizeof(double *)*fzdim);

  if((tmp == NULL) || (buff == NULL) || (sortedv == NULL)) {
    fprintf(stderr,"Memory allocation error\n");
    exit(EXIT_FAILURE);
  }

  startz = ((fzdim+zoff-1<0) ? fzdim+zoff-1 : 0);
  endz   = zdim+fzdim+zoff-1;

  for (z=startz; z<endz; z++)
  {
    double sum2 = 0.0;

    if (z >= 0 && z<zdim)
    {
      for (y=0;y<ydim;y++) for (x=0;x<xdim;x++)
        tmp[((z%fzdim)*xdim*ydim)+(y*xdim)+x] = iVol[(z*xdim*ydim)+(y*xdim)+x];   
      convxy(tmp+((z%fzdim)*xdim*ydim),xdim, ydim,
        filtx, filty, fxdim, fydim, xoff, yoff, buff);
    }
    if (z-fzdim-zoff+1>=0 && z-fzdim-zoff+1<zdim)
    {
      fstart = ((z >= zdim) ? z-zdim+1 : 0);
      fend = ((z-fzdim < 0) ? z+1 : fzdim);

      for(k=0; k<fzdim; k++)
      {
        int z1 = (((z-k)%fzdim)+fzdim)%fzdim;
        sortedv[k] = &(tmp[z1*xdim*ydim]);
      }

      for(k=fstart, sum2=0.0; k<fend; k++)
        sum2 += filtz[k];

      obuf = &oVol[(z-fzdim-zoff+1)*ydim*xdim];
      if (sum2)
      {
        for(xy=0; xy<xdim*ydim; xy++)
        {
          double sum1=0.0;
          for(k=fstart; k<fend; k++)
            sum1 += filtz[k]*sortedv[k][xy];

          obuf[xy] = sum1/sum2;
        }
      }
      else
        for(xy=0; xy<xdim*ydim; xy++)
          obuf[xy] = 0.0;
    }
  }
  free(tmp);
  free(buff);
  free(sortedv);
  return(0);
}

int 
convxyz_float(float *iVol, double filtx[], double filty[], double filtz[],
  int fxdim, int fydim, int fzdim, int xoff, int yoff, int zoff,
  float *oVol, int dims[3])
{
  float *tmp, *buff, **sortedv, *obuf;
  int xy, z, y, x, k, fstart, fend, startz, endz;
  int xdim, ydim, zdim;

  xdim = dims[0];
  ydim = dims[1];
  zdim = dims[2];

  tmp = (float *)malloc(sizeof(float)*xdim*ydim*fzdim);
  buff = (float *)malloc(sizeof(float)*((ydim>xdim) ? ydim : xdim));
  sortedv = (float **)malloc(sizeof(float *)*fzdim);

  if((tmp == NULL) || (buff == NULL) || (sortedv == NULL)) {
    fprintf(stderr,"Memory allocation error\n");
    exit(EXIT_FAILURE);
  }

  startz = ((fzdim+zoff-1<0) ? fzdim+zoff-1 : 0);
  endz   = zdim+fzdim+zoff-1;

  for (z=startz; z<endz; z++)
  {
    double sum2 = 0.0;

    if (z >= 0 && z<zdim)
    {
      for (y=0;y<ydim;y++) for (x=0;x<xdim;x++)
        tmp[((z%fzdim)*xdim*ydim)+(y*xdim)+x] = iVol[(z*xdim*ydim)+(y*xdim)+x];   
        convxy_float(tmp+((z%fzdim)*xdim*ydim),xdim, ydim,
            filtx, filty, fxdim, fydim, xoff, yoff, buff);
    }
    if (z-fzdim-zoff+1>=0 && z-fzdim-zoff+1<zdim)
    {
      fstart = ((z >= zdim) ? z-zdim+1 : 0);
      fend = ((z-fzdim < 0) ? z+1 : fzdim);

      for(k=0; k<fzdim; k++)
      {
        int z1 = (((z-k)%fzdim)+fzdim)%fzdim;
        sortedv[k] = &(tmp[z1*xdim*ydim]);
      }

      for(k=fstart, sum2=0.0; k<fend; k++)
        sum2 += filtz[k];

      obuf = &oVol[(z-fzdim-zoff+1)*ydim*xdim];
      if (sum2)
      {
        for(xy=0; xy<xdim*ydim; xy++)
        {
          double sum1=0.0;
          for(k=fstart; k<fend; k++)
            sum1 += filtz[k]*(double)sortedv[k][xy];

          obuf[xy] = (float)sum1/sum2;
        }
      }
      else
        for(xy=0; xy<xdim*ydim; xy++)
          obuf[xy] = 0.0;
    }
  }
  free(tmp);
  free(buff);
  free(sortedv);
  return(0);
}

int
convxyz_uint8(unsigned char *iVol, double filtx[], double filty[], double filtz[],
  int fxdim, int fydim, int fzdim, int xoff, int yoff, int zoff,
  unsigned char *oVol, int dims[3])
{
  double *tmp, *buff, **sortedv;
  int xy, z, y, x, k, fstart, fend, startz, endz;
  int xdim, ydim, zdim;
  double tmp2;
  unsigned char *obuf;

  xdim = dims[0];
  ydim = dims[1];
  zdim = dims[2];

  tmp = (double *)malloc(sizeof(double)*xdim*ydim*fzdim);
  buff = (double *)malloc(sizeof(double)*((ydim>xdim) ? ydim : xdim));
  sortedv = (double **)malloc(sizeof(double *)*fzdim);

  if((tmp == NULL) || (buff == NULL) || (sortedv == NULL)) {
    fprintf(stderr,"Memory allocation error\n");
    exit(EXIT_FAILURE);
  }

  startz = ((fzdim+zoff-1<0) ? fzdim+zoff-1 : 0);
  endz   = zdim+fzdim+zoff-1;

  for (z=startz; z<endz; z++)
  {
    double sum2 = 0.0;

    if (z >= 0 && z<zdim)
    {
      for (y=0;y<ydim;y++) for (x=0;x<xdim;x++)
        tmp[((z%fzdim)*xdim*ydim)+(y*xdim)+x] = (double)iVol[(z*xdim*ydim)+(y*xdim)+x];   
      convxy(tmp+((z%fzdim)*xdim*ydim),xdim, ydim,
        filtx, filty, fxdim, fydim, xoff, yoff, buff);
    }
    if (z-fzdim-zoff+1>=0 && z-fzdim-zoff+1<zdim)
    {
      fstart = ((z >= zdim) ? z-zdim+1 : 0);
      fend = ((z-fzdim < 0) ? z+1 : fzdim);

      for(k=0; k<fzdim; k++)
      {
        int z1 = (((z-k)%fzdim)+fzdim)%fzdim;
        sortedv[k] = &(tmp[z1*xdim*ydim]);
      }

      for(k=fstart, sum2=0.0; k<fend; k++)
        sum2 += filtz[k];

      obuf = oVol;
      obuf = &obuf[(z-fzdim-zoff+1)*ydim*xdim];
      if (sum2)
      {
        for(xy=0; xy<xdim*ydim; xy++)
        {
          double sum1=0.0;
          for(k=fstart; k<fend; k++)
            sum1 += filtz[k]*sortedv[k][xy];
          tmp2 = sum1/sum2;
          if (tmp2<0.0) tmp2 = 0.0;
          else if (tmp2>255.0) tmp2 = 255.0;
          obuf[xy] = (unsigned char)RINT(tmp2);
        }
      }
      else
        for(xy=0; xy<xdim*ydim; xy++)
          obuf[xy] = 0;
    }
  }
  free(tmp);
  free(buff);
  free(sortedv);
  return(0);
}

int 
convxyz_int16(signed short *iVol, double filtx[], double filty[], double filtz[],
  int fxdim, int fydim, int fzdim, int xoff, int yoff, int zoff,
  signed short *oVol, int dims[3])
{
  double *tmp, *buff, **sortedv;
  int xy, z, y, x, k, fstart, fend, startz, endz;
  int xdim, ydim, zdim;
  double tmp2;
  signed short *obuf;

  xdim = dims[0];
  ydim = dims[1];
  zdim = dims[2];

  tmp = (double *)malloc(sizeof(double)*xdim*ydim*fzdim);
  buff = (double *)malloc(sizeof(double)*((ydim>xdim) ? ydim : xdim));
  sortedv = (double **)malloc(sizeof(double *)*fzdim);

  if((tmp == NULL) || (buff == NULL) || (sortedv == NULL)) {
    fprintf(stderr,"Memory allocation error\n");
    exit(EXIT_FAILURE);
  }

  startz = ((fzdim+zoff-1<0) ? fzdim+zoff-1 : 0);
  endz   = zdim+fzdim+zoff-1;

  for (z=startz; z<endz; z++)
  {
    double sum2 = 0.0;

    if (z >= 0 && z<zdim)
    {
      for (y=0;y<ydim;y++) for (x=0;x<xdim;x++)
        tmp[((z%fzdim)*xdim*ydim)+(y*xdim)+x] = (double)iVol[(z*xdim*ydim)+(y*xdim)+x];   
      convxy(tmp+((z%fzdim)*xdim*ydim),xdim, ydim,
        filtx, filty, fxdim, fydim, xoff, yoff, buff);
    }
    if (z-fzdim-zoff+1>=0 && z-fzdim-zoff+1<zdim)
    {
      fstart = ((z >= zdim) ? z-zdim+1 : 0);
      fend = ((z-fzdim < 0) ? z+1 : fzdim);

      for(k=0; k<fzdim; k++)
      {
        int z1 = (((z-k)%fzdim)+fzdim)%fzdim;
        sortedv[k] = &(tmp[z1*xdim*ydim]);
      }

      for(k=fstart, sum2=0.0; k<fend; k++)
        sum2 += filtz[k];

      obuf = oVol;
      obuf = &obuf[(z-fzdim-zoff+1)*ydim*xdim];
      if (sum2)
      {
        for(xy=0; xy<xdim*ydim; xy++)
        {
          double sum1=0.0;
          for(k=fstart; k<fend; k++)
            sum1 += filtz[k]*sortedv[k][xy];
          tmp2 = sum1/sum2;
          if (tmp2<-32768.0) tmp2 = -32768.0;
          else if (tmp2>32767.0) tmp2 = 32767.0;
          obuf[xy] = (signed short)RINT(tmp2);
        }
      }
      else
        for(xy=0; xy<xdim*ydim; xy++)
          obuf[xy] = 0;
    }
  }
  free(tmp);
  free(buff);
  free(sortedv);
  return(0);
}

int 
convxyz_int32(signed int *iVol, double filtx[], double filty[], double filtz[],
  int fxdim, int fydim, int fzdim, int xoff, int yoff, int zoff,
  signed int *oVol, int dims[3])
{
  double *tmp, *buff, **sortedv;
  int xy, z, y, x, k, fstart, fend, startz, endz;
  int xdim, ydim, zdim;
  double tmp2;
  signed int *obuf;

  xdim = dims[0];
  ydim = dims[1];
  zdim = dims[2];

  tmp = (double *)malloc(sizeof(double)*xdim*ydim*fzdim);
  buff = (double *)malloc(sizeof(double)*((ydim>xdim) ? ydim : xdim));
  sortedv = (double **)malloc(sizeof(double *)*fzdim);

  if((tmp == NULL) || (buff == NULL) || (sortedv == NULL)) {
    fprintf(stderr,"Memory allocation error\n");
    exit(EXIT_FAILURE);
  }

  startz = ((fzdim+zoff-1<0) ? fzdim+zoff-1 : 0);
  endz   = zdim+fzdim+zoff-1;

  for (z=startz; z<endz; z++)
  {
    double sum2 = 0.0;

    if (z >= 0 && z<zdim)
    {
      for (y=0;y<ydim;y++) for (x=0;x<xdim;x++)
        tmp[((z%fzdim)*xdim*ydim)+(y*xdim)+x] = (double)iVol[(z*xdim*ydim)+(y*xdim)+x];   
      convxy(tmp+((z%fzdim)*xdim*ydim),xdim, ydim,
        filtx, filty, fxdim, fydim, xoff, yoff, buff);
    }
    if (z-fzdim-zoff+1>=0 && z-fzdim-zoff+1<zdim)
    {
      fstart = ((z >= zdim) ? z-zdim+1 : 0);
      fend = ((z-fzdim < 0) ? z+1 : fzdim);

      for(k=0; k<fzdim; k++)
      {
        int z1 = (((z-k)%fzdim)+fzdim)%fzdim;
        sortedv[k] = &(tmp[z1*xdim*ydim]);
      }

      for(k=fstart, sum2=0.0; k<fend; k++)
        sum2 += filtz[k];

      obuf = oVol;
      obuf = &obuf[(z-fzdim-zoff+1)*ydim*xdim];
      if (sum2)
      {
        for(xy=0; xy<xdim*ydim; xy++)
        {
          double sum1=0.0;
          for(k=fstart; k<fend; k++)
            sum1 += filtz[k]*sortedv[k][xy];
          tmp2 = sum1/sum2;
          if (tmp2<-2147483648.0) tmp2 = -2147483648.0;
          else if (tmp2>2147483647.0) tmp2 = 2147483647.0;
          obuf[xy] = (signed int)RINT(tmp2);
        }
      }
      else
        for(xy=0; xy<xdim*ydim; xy++)
          obuf[xy] = 0;
    }
  }
  free(tmp);
  free(buff);
  free(sortedv);
  return(0);
}

/* estimate minimum of A and its index in A */
void
pmin(float A[], int sA, float *minimum, int *index)
{
  int i; 
  *minimum=FLT_MAX; *index=0; /* printf("%d ",sizeof(A)/8); */
  for(i=0;i<sA;i++) {
    if ((A[i]>0) && (*minimum>A[i]))
    { 
      *minimum = A[i]; 
      *index   = i;
    }
  }
}

/* estimate x,y,z position of index i in an array size sx,sxy=sx*sy... */
void
ind2sub(int i,int *x,int *y, int *z, int sxy, int sy) {
  *z = (int)floor( (double)i / (double)sxy ) +1; 
   i = i % (sxy);
  *y = (int)floor( (double)i / (double)sy ) +1;        
  *x = i % sy + 1;
}

void
vbdist(float *V, int *dims, double *voxelsize) 
{
  
  /* main information about input data (size, dimensions, ...) */
  const int     nL = dims[0]*dims[1]*dims[2];
  const int     x  = dims[0];
  const int     y  = dims[1];
  const int     xy = x*y;

  float s1 = (float)voxelsize[0],s2 = (float)voxelsize[1],s3 = (float)voxelsize[2];
  const float   s12  = sqrt( s1*s1  + s2*s2); /* xy - voxel size */
  const float   s13  = sqrt( s1*s1  + s3*s3); /* xz - voxel size */
  const float   s23  = sqrt( s2*s2  + s3*s3); /* yz - voxel size */
  const float   s123 = sqrt(s12*s12 + s3*s3); /* xyz - voxel size */
  
  /* indices of the neighbor Ni (index distance) and euclidean distance NW */
  const int   NI[] = {  0, -1,-x+1, -x,-x-1,  -xy+1,-xy,-xy-1,  -xy+x+1,-xy+x,-xy+x-1,  -xy-x+1,-xy-x,-xy-x-1};  
  const float ND[] = {0.0, s1, s12, s2, s12,    s13, s3,  s13,     s123,  s23,   s123,     s123,  s23,   s123};
  const int   sN = sizeof(NI)/4;    
  float       DN[sN],DI[sN];
  float       DNm = FLT_MAX;
  double      mod;
  int i, n, ni, DNi = 0;

  
  /* data */
  float         *D;
  unsigned int  *I;
  
  D = (float *)malloc(sizeof(float)*nL);
  I = (unsigned int *)malloc(sizeof(unsigned int)*nL);
  
  
  /* intitialisiation */
  for (i=0;i<nL;i++) 
  {
    if (V[i]>0.5) D[i]=0.0; else D[i]=FLT_MAX; 
    I[i]=(unsigned int)i;
  }

  int u,v,w,nu,nv,nw; 
  for (i=0;i<nL;i++) 
  {
    if (D[i]>0)
    {
      ind2sub(i,&u,&v,&w,xy,x);
      
      /* read neighbor values */
      for (n=0;n<sN;n++)
      {
        ni = i + NI[n];
        ind2sub(ni,&nu,&nv,&nw,xy,x);
        if ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) ) ni=i;
        DN[n] = D[ni] + ND[n];
      }

      /* find minimum distance within the neighborhood */
      pmin(DN,sN,&DNm,&DNi);

      /* update values */
      if (DNi>0) {
        I[i] = (unsigned int)  I[i+NI[DNi]];
        D[i] = DNm; 
        ind2sub((int)I[i],&nu,&nv,&nw,xy,x); 
        D[i] = sqrt(pow((float)(u-nu)*s1,2) + pow((float)(v-nv)*s2,2) + pow((float)(w-nw)*s3,2));
      }
    }
  }
  
  for (i=nL-1;i>0;i--)
  {
    if (D[i]>0)
    {
      ind2sub(i,&u,&v,&w,xy,x);

      /* read neighbour values */
      for (n=0;n<sN;n++)
      {
        ni = i - NI[n];
        ind2sub(ni,&nu,&nv,&nw,xy,x);
        if ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) ) ni=i;
        DN[n] = D[ni] + ND[n];
      }

      /* find minimum distance within the neighborhood */
      pmin(DN,sN,&DNm,&DNi);

      /* update values */
      if (DNi>0) {
        I[i] = (unsigned int)  I[i-NI[DNi]];
        D[i] = DNm; 
        ind2sub((int)I[i],&nu,&nv,&nw,xy,x); 
        D[i] = sqrt(pow((float)(u-nu)*s1,2) + pow((float)(v-nv)*s2,2) + pow((float)(w-nw)*s3,2));
      }
    }
  }

  for (i=0;i<nL;i++) 
    V[i] = D[i];
    
  free(D);
  free(I);
}

void
distclose_uint8(unsigned char *vol, int dims[3], double voxelsize[3], int niter, double th)
{
  float *buffer;
  int i,x,y,z,j,band,dims2[3];
  unsigned char max_vol;
  
  if (niter < 1) return;

  for (i=0; i<dims[0]*dims[1]*dims[2]; i++) max_vol = MAX(max_vol,vol[i]);
  th *= (double)max_vol;

  /* add band with zeros to image to avoid clipping */  
  band = niter;
  for (i=0;i<3;i++) dims2[i] = dims[i] + 2*band;

  buffer = (float *)malloc(sizeof(float)*dims2[0]*dims2[1]*dims2[2]);

  if(buffer == NULL) {
    fprintf(stderr,"Memory allocation error\n");
    exit(EXIT_FAILURE);
  }
  
  memset(buffer,0,sizeof(float)*dims2[0]*dims2[1]*dims2[2]);
  
  /* threshold input */
  for (z=0;z<dims[2];z++) for (y=0;y<dims[1];y++) for (x=0;x<dims[0];x++) 
    buffer[index(x+band,y+band,z+band,dims2)] = (float)((double)vol[index(x,y,z,dims)]>th);
        
  vbdist(buffer, dims2, voxelsize);
  for (i=0;i<dims2[2]*dims2[1]*dims2[0];i++)
    buffer[i] = buffer[i] > (float)niter;

  vbdist(buffer, dims2, voxelsize);
  for (i=0;i<dims2[2]*dims2[1]*dims2[0];i++)
    buffer[i] = buffer[i] > (float)niter;

  /* return image */
  for (z=0;z<dims[2];z++) for (y=0;y<dims[1];y++) for (x=0;x<dims[0];x++) 
    vol[index(x,y,z,dims)] = (unsigned char)buffer[index(x+band,y+band,z+band,dims2)];
    
  free(buffer);
}

void
distclose_float(float *vol, int dims[3], double voxelsize[3], int niter, double th)
{
  float *buffer;
  int i,x,y,z,j,band,dims2[3];
  float max_vol;
  
  if (niter < 1) return;

  for (i=0; i<dims[0]*dims[1]*dims[2]; i++) max_vol = MAX(max_vol,vol[i]);
  th *= (double)max_vol;

  /* add band with zeros to image to avoid clipping */  
  band = niter;
  for (i=0;i<3;i++) dims2[i] = dims[i] + 2*band;

  buffer = (float *)malloc(sizeof(float)*dims2[0]*dims2[1]*dims2[2]);

  if(buffer == NULL) {
    fprintf(stderr,"Memory allocation error\n");
    exit(EXIT_FAILURE);
  }
  
  memset(buffer,0,sizeof(float)*dims2[0]*dims2[1]*dims2[2]);
  
  /* threshold input */
  for (z=0;z<dims[2];z++) for (y=0;y<dims[1];y++) for (x=0;x<dims[0];x++) 
    buffer[index(x+band,y+band,z+band,dims2)] = (vol[index(x,y,z,dims)]>(float)th);
        
  vbdist(buffer, dims2, voxelsize);
  for (i=0;i<dims2[2]*dims2[1]*dims2[0];i++)
    buffer[i] = buffer[i] > (float)niter;

  vbdist(buffer, dims2, voxelsize);
  for (i=0;i<dims2[2]*dims2[1]*dims2[0];i++)
    buffer[i] = buffer[i] > (float)niter;

  /* return image */
  for (z=0;z<dims[2];z++) for (y=0;y<dims[1];y++) for (x=0;x<dims[0];x++) 
    vol[index(x,y,z,dims)] = buffer[index(x+band,y+band,z+band,dims2)];
    
  free(buffer);
}

void
distopen_uint8(unsigned char *vol, int dims[3], double voxelsize[3], int niter, double th)
{
  float *buffer;
  int i,j;
  unsigned char max_vol;
  
  if (niter < 1) return;

  for (i=0; i<dims[0]*dims[1]*dims[2]; i++) max_vol = MAX(max_vol,vol[i]);
  th *= (double)max_vol;
  
  buffer = (float *)malloc(sizeof(float)*dims[0]*dims[1]*dims[2]);

  if(buffer == NULL) {
    fprintf(stderr,"Memory allocation error\n");
    exit(EXIT_FAILURE);
  }
  
  /* threshold input */
  for (i=0;i<dims[2]*dims[1]*dims[0];i++)
    buffer[i] = 1.0 - (float)((double)vol[i]>th);
        
  vbdist(buffer, dims, voxelsize);
  for (i=0;i<dims[2]*dims[1]*dims[0];i++)
    buffer[i] = buffer[i] > (float)niter;

  vbdist(buffer, dims, voxelsize);
  for (i=0;i<dims[2]*dims[1]*dims[0];i++)
    buffer[i] = buffer[i] < (float)niter;

  /* return image */
  for (i=0;i<dims[2]*dims[1]*dims[0];i++)
    vol[i] = (unsigned char)buffer[i];

  free(buffer);
}

void
distopen_float(float *vol, int dims[3], double voxelsize[3], int niter, double th)
{
  float *buffer;
  int i,j;
  float max_vol;
  
  if (niter < 1) return;

  for (i=0; i<dims[0]*dims[1]*dims[2]; i++) max_vol = MAX(max_vol,vol[i]);
  th *= (double)max_vol;
  
  buffer = (float *)malloc(sizeof(float)*dims[0]*dims[1]*dims[2]);

  if(buffer == NULL) {
    fprintf(stderr,"Memory allocation error\n");
    exit(EXIT_FAILURE);
  }
  
  /* threshold input */
  for (i=0;i<dims[2]*dims[1]*dims[0];i++)
    buffer[i] = 1.0 - ((float)vol[i]>th);
        
  vbdist(buffer, dims, voxelsize);
  for (i=0;i<dims[2]*dims[1]*dims[0];i++)
    buffer[i] = buffer[i] > (float)niter;

  vbdist(buffer, dims, voxelsize);
  for (i=0;i<dims[2]*dims[1]*dims[0];i++)
    buffer[i] = buffer[i] < (float)niter;

  /* return image */
  for (i=0;i<dims[2]*dims[1]*dims[0];i++)
    vol[i] = buffer[i];

  free(buffer);
}

void
morph_erode_uint8(unsigned char *vol, int dims[3], int niter, double th)
{
  double filt[3]={1,1,1};
  int i,j;
  unsigned char max_vol;
  
  if (niter < 1) return;

  for (i=0; i<dims[0]*dims[1]*dims[2]; i++) max_vol = MAX(max_vol,vol[i]);
  th *= (double)max_vol;
  
  /* threshold input */
  for (j=0;j<dims[2]*dims[1]*dims[0];j++)
    vol[j] = (unsigned char)((double)vol[j]>th);

  for (i=0;i<niter;i++) {
    convxyz_uint8(vol,filt,filt,filt,3,3,3,-1,-1,-1,vol,dims);
    for (j=0;j<dims[2]*dims[1]*dims[0];j++)
      vol[j] = (vol[j]>=9);
  }
}

void
morph_erode_float(float *vol, int dims[3], int niter, double th)
{
  double filt[3]={1,1,1};
  int i,j;
  float max_vol;
  
  if (niter < 1) return;

  for (i=0; i<dims[0]*dims[1]*dims[2]; i++) max_vol = MAX(max_vol,vol[i]);
  th *= (double)max_vol;

  /* threshold input */
  for (j=0;j<dims[2]*dims[1]*dims[0];j++)
    vol[j] = vol[j]>(float)th;

  for (i=0;i<niter;i++) {
    convxyz_float(vol,filt,filt,filt,3,3,3,-1,-1,-1,vol,dims);
    for (j=0;j<dims[2]*dims[1]*dims[0];j++)
      vol[j] = vol[j]>=9.0;
  }
}

void
morph_dilate_uint8(unsigned char *vol, int dims[3], int niter, double th)
{
  double filt[3]={1,1,1};
  int i,x,y,z,j,band,dims2[3];
  unsigned char max_vol;
  unsigned char *buffer;
  
  if (niter < 1) return;

  for (i=0; i<dims[0]*dims[1]*dims[2]; i++) max_vol = MAX(max_vol,vol[i]);
  th *= (double)max_vol;

  /* add band with zeros to image to avoid clipping */  
  band = niter;
  for (i=0;i<3;i++) dims2[i] = dims[i] + 2*band;

  buffer = (unsigned char *)malloc(sizeof(unsigned char)*dims2[0]*dims2[1]*dims2[2]);

  if(buffer == NULL) {
    fprintf(stderr,"Memory allocation error\n");
    exit(EXIT_FAILURE);
  }
  
  memset(buffer,0,sizeof(unsigned char)*dims2[0]*dims2[1]*dims2[2]);
  
  /* threshold input */
  for (x=0;x<dims[0];x++) for (y=0;y<dims[1];y++) for (z=0;z<dims[2];z++) 
    buffer[index(x+band,y+band,z+band,dims2)] = (unsigned char)((double)vol[index(x,y,z,dims)]>th);

  for (i=0;i<niter;i++) {
    convxyz_uint8(buffer,filt,filt,filt,3,3,3,-1,-1,-1,buffer,dims2);
    for (j=0; j<dims2[0]*dims2[1]*dims2[2]; j++) 
      buffer[j] = buffer[j]>0;
  }

  /* return image */
  for (x=0;x<dims[0];x++) for (y=0;y<dims[1];y++) for (z=0;z<dims[2];z++) 
    vol[index(x,y,z,dims)] = buffer[index(x+band,y+band,z+band,dims2)];
    
  free(buffer);
  
}

void
morph_dilate_float(float *vol, int dims[3], int niter, double th)
{
  double filt[3]={1,1,1};
  int i,x,y,z,j,band,dims2[3];
  float max_vol;
  unsigned char *buffer;
  
  if (niter < 1) return;

  for (i=0; i<dims[0]*dims[1]*dims[2]; i++) max_vol = MAX(max_vol,vol[i]);
  th *= (double)max_vol;

  /* add band with zeros to image to avoid clipping */  
  band = niter;
  for (i=0;i<3;i++) dims2[i] = dims[i] + 2*band;

  buffer = (unsigned char *)malloc(sizeof(unsigned char)*dims2[0]*dims2[1]*dims2[2]);

  if(buffer == NULL) {
    fprintf(stderr,"Memory allocation error\n");
    exit(EXIT_FAILURE);
  }
  
  memset(buffer,0,sizeof(unsigned char)*dims2[0]*dims2[1]*dims2[2]);
  
  /* threshold input */
  for (x=0;x<dims[0];x++) for (y=0;y<dims[1];y++) for (z=0;z<dims[2];z++) 
    buffer[index(x+band,y+band,z+band,dims2)] = (unsigned char)((double)vol[index(x,y,z,dims)]>th);

  for (i=0;i<niter;i++) {
    convxyz_uint8(buffer,filt,filt,filt,3,3,3,-1,-1,-1,buffer,dims2);
    for (j=0; j<dims2[0]*dims2[1]*dims2[2]; j++) 
      buffer[j] = buffer[j]>0;
  }

  /* return image */
  for (x=0;x<dims[0];x++) for (y=0;y<dims[1];y++) for (z=0;z<dims[2];z++) 
    vol[index(x,y,z,dims)] = (float)buffer[index(x+band,y+band,z+band,dims2)];
    
  free(buffer);
  
}

void
morph_close_uint8(unsigned char *vol, int dims[3], int niter, double th)
{
  double filt[3]={1,1,1};
  unsigned char *buffer;
  int i,x,y,z,j,band,dims2[3];
  unsigned char max_vol;
  
  if (niter < 1) return;

  for (i=0; i<dims[0]*dims[1]*dims[2]; i++) max_vol = MAX(max_vol,vol[i]);
  th *= (double)max_vol;

  /* add band with zeros to image to avoid clipping */  
  band = niter;
  for (i=0;i<3;i++) dims2[i] = dims[i] + 2*band;

  buffer = (unsigned char *)malloc(sizeof(unsigned char)*dims2[0]*dims2[1]*dims2[2]);

  if(buffer == NULL) {
    fprintf(stderr,"Memory allocation error\n");
    exit(EXIT_FAILURE);
  }
  
  memset(buffer,0,sizeof(unsigned char)*dims2[0]*dims2[1]*dims2[2]);
  
  /* threshold input */
  for (x=0;x<dims[0];x++) for (y=0;y<dims[1];y++) for (z=0;z<dims[2];z++) 
    buffer[index(x+band,y+band,z+band,dims2)] = (unsigned char)((double)vol[index(x,y,z,dims)]>th);

  /* dilate */
  for (i=0;i<niter;i++) {
    convxyz_uint8(buffer,filt,filt,filt,3,3,3,-1,-1,-1,buffer,dims2);
    for (j=0; j<dims2[0]*dims2[1]*dims2[2]; j++) 
      buffer[j] = (buffer[j]>0);
  }

  /* erode */
  for (i=0;i<niter;i++) {
    convxyz_uint8(buffer,filt,filt,filt,3,3,3,-1,-1,-1,buffer,dims2);
    for (j=0; j<dims2[0]*dims2[1]*dims2[2]; j++) 
      buffer[j] = (buffer[j]>=9);
  }

  /* return image */
  for (x=0;x<dims[0];x++) for (y=0;y<dims[1];y++) for (z=0;z<dims[2];z++) 
    vol[index(x,y,z,dims)] = buffer[index(x+band,y+band,z+band,dims2)];
    
  free(buffer);
}

void
morph_close_float(float *vol, int dims[3], int niter, double th)
{
  double filt[3]={1,1,1};
  unsigned char *buffer;
  int i,x,y,z,j,band,dims2[3];
  float max_vol;
  
  if (niter < 1) return;

  for (i=0; i<dims[0]*dims[1]*dims[2]; i++) max_vol = MAX(max_vol,vol[i]);
  th *= (double)max_vol;

  /* add band with zeros to image to avoid clipping */  
  band = niter;
  for (i=0;i<3;i++) dims2[i] = dims[i] + 2*band;

  buffer = (unsigned char *)malloc(sizeof(unsigned char)*dims2[0]*dims2[1]*dims2[2]);

  if(buffer == NULL) {
    fprintf(stderr,"Memory allocation error\n");
    exit(EXIT_FAILURE);
  }
  
  memset(buffer,0,sizeof(unsigned char)*dims2[0]*dims2[1]*dims2[2]);
  
  /* threshold input */
  for (x=0;x<dims[0];x++) for (y=0;y<dims[1];y++) for (z=0;z<dims[2];z++) 
    buffer[index(x+band,y+band,z+band,dims2)] = (unsigned char)((double)vol[index(x,y,z,dims)]>th);
        
  /* dilate */
  for (i=0;i<niter;i++) {
    convxyz_uint8(buffer,filt,filt,filt,3,3,3,-1,-1,-1,buffer,dims2);
    for (j=0; j<dims2[0]*dims2[1]*dims2[2]; j++) 
      buffer[j] = (buffer[j]>0);
  }

  /* erode */
  for (i=0;i<niter;i++) {
    convxyz_uint8(buffer,filt,filt,filt,3,3,3,-1,-1,-1,buffer,dims2);
    for (j=0; j<dims2[0]*dims2[1]*dims2[2]; j++) 
      buffer[j] = (buffer[j]>=9);
  }

  /* return image */
  for (x=0;x<dims[0];x++) for (y=0;y<dims[1];y++) for (z=0;z<dims[2];z++) 
    vol[index(x,y,z,dims)] = (float)buffer[index(x+band,y+band,z+band,dims2)];
    
  free(buffer);
}

void
morph_open_uint8(unsigned char *vol, int dims[3], int niter, double th)
{
  double filt[3]={1,1,1};
  int i, j;
  unsigned char max_vol;
  
  if (niter < 1) return;

  for (i=0; i<dims[0]*dims[1]*dims[2]; i++) max_vol = MAX(max_vol,vol[i]);
  th *= (double)max_vol;

  /* threshold input */
  for (j=0;j<dims[2]*dims[1]*dims[0];j++)
    vol[j] = (unsigned char)((double)vol[j]>th);

  for (i=0;i<niter;i++) {
    convxyz_uint8(vol,filt,filt,filt,3,3,3,-1,-1,-1,vol,dims);
    for (j=0;j<dims[2]*dims[1]*dims[0];j++)
      vol[j] = (vol[j]>=9);
  }

  for (i=0;i<niter;i++) {
    convxyz_uint8(vol,filt,filt,filt,3,3,3,-1,-1,-1,vol,dims);
    for (j=0; j<dims[0]*dims[1]*dims[2]; j++) 
      vol[j] = (vol[j]>0);
  }

}

void
morph_open_float(float *vol, int dims[3], int niter, double th)
{
  unsigned char *buffer;
  double filt[3]={1,1,1};
  int i, j;
  float max_vol;
  
  if (niter < 1) return;

  for (i=0; i<dims[0]*dims[1]*dims[2]; i++) max_vol = MAX(max_vol,vol[i]);
  th *= (double)max_vol;

  buffer = (unsigned char *)malloc(sizeof(unsigned char)*dims[0]*dims[1]*dims[2]);

  if(buffer == NULL) {
    fprintf(stderr,"Memory allocation error\n");
    exit(EXIT_FAILURE);
  }

  /* threshold input */
  for (j=0;j<dims[2]*dims[1]*dims[0];j++)
    buffer[j] = (unsigned char)((double)vol[j]>th);

  for (i=0;i<niter;i++) {
    convxyz_uint8(buffer,filt,filt,filt,3,3,3,-1,-1,-1,buffer,dims);
    for (j=0;j<dims[2]*dims[1]*dims[0];j++)
      buffer[j] = (buffer[j]>=9);
  }

  for (i=0;i<niter;i++) {
    convxyz_uint8(buffer,filt,filt,filt,3,3,3,-1,-1,-1,buffer,dims);
    for (j=0; j<dims[0]*dims[1]*dims[2]; j++) 
      buffer[j] = (buffer[j]>0);
  }

  for (i=0;i<dims[2]*dims[1]*dims[0];i++)
    vol[i] = (float)buffer[i];
    
  free(buffer);
}

/* First order hold resampling - trilinear interpolation */
void 
subsample_double(double *in, double *out, int dim_in[3], int dim_out[3], int offset_in, int offset_out)
{
  int i, x, y, z;
  double k111,k112,k121,k122,k211,k212,k221,k222;
  double dx1, dx2, dy1, dy2, dz1, dz2, xi, yi, zi, samp[3];
  int off1, off2, xcoord, ycoord, zcoord;

  for (i=0; i<3; i++) {
    if(dim_out[i] > dim_in[i]) samp[i] = ceil((double)dim_out[i]/(double)dim_in[i]);
    else                       samp[i] = 1.0/(ceil((double)dim_in[i]/(double)dim_out[i]));
  }
  
  for (z=0; z<dim_out[2]; z++) {
    zi = 1.0+(double)z/samp[2];
    for (y=0; y<dim_out[1]; y++) {
      yi = 1.0+(double)y/samp[1];
      for (x=0; x<dim_out[0]; x++) {
        xi = 1.0+(double)x/samp[0];
        i = z*dim_out[0]*dim_out[1] + y*dim_out[0] + x + offset_out;

        if (zi>=0 && zi<dim_in[2] && yi>=0 && yi<dim_in[1] && xi>=0 && xi<dim_in[0])  {
          xcoord = (int)floor(xi); dx1=xi-(double)xcoord; dx2=1.0-dx1;
          ycoord = (int)floor(yi); dy1=yi-(double)ycoord; dy2=1.0-dy1;
          zcoord = (int)floor(zi); dz1=zi-(double)zcoord; dz2=1.0-dz1;

          off1 = xcoord-1 + dim_in[0]*(ycoord-1 + dim_in[1]*(zcoord-1)) + offset_in;
          k222 = (double)in[off1]; k122 = (double)in[off1+1]; off2 = off1+dim_in[0];
          k212 = (double)in[off2]; k112 = (double)in[off2+1]; off1+= dim_in[0]*dim_in[1];
          k221 = (double)in[off1]; k121 = (double)in[off1+1]; off2 = off1+dim_in[0];
          k211 = (double)in[off2]; k111 = (double)in[off2+1];

          out[i] = ((((k222*dx2 + k122*dx1)*dy2 + (k212*dx2 + k112*dx1)*dy1))*dz2
                         + (((k221*dx2 + k121*dx1)*dy2 + (k211*dx2 + k111*dx1)*dy1))*dz1);
                 
        } else out[i] = 0;
      }
    }
  }
}

/* First order hold resampling - trilinear interpolation */
void subsample_uint8(unsigned char *in, float *out, int dim_in[3], int dim_out[3], int offset_in, int offset_out)
{
  int i, x, y, z;
  double k111,k112,k121,k122,k211,k212,k221,k222;
  double dx1, dx2, dy1, dy2, dz1, dz2, xi, yi, zi, samp[3];
  int off1, off2, xcoord, ycoord, zcoord;

  for (i=0; i<3; i++) {
    if(dim_out[i] > dim_in[i]) samp[i] = ceil((double)dim_out[i]/(double)dim_in[i]);
    else                       samp[i] = 1.0/(ceil((double)dim_in[i]/(double)dim_out[i]));
  }
  
  for (z=0; z<dim_out[2]; z++) {
    zi = 1.0+(double)z/samp[2];
    for (y=0; y<dim_out[1]; y++) {
      yi = 1.0+(double)y/samp[1];
      for (x=0; x<dim_out[0]; x++) {
        xi = 1.0+(double)x/samp[0];
        i = z*dim_out[0]*dim_out[1] + y*dim_out[0] + x + offset_out;

        if (zi>=0 && zi<dim_in[2] && yi>=0 && yi<dim_in[1] && xi>=0 && xi<dim_in[0])  {
          xcoord = (int)floor(xi); dx1=xi-(double)xcoord; dx2=1.0-dx1;
          ycoord = (int)floor(yi); dy1=yi-(double)ycoord; dy2=1.0-dy1;
          zcoord = (int)floor(zi); dz1=zi-(double)zcoord; dz2=1.0-dz1;

          off1 = xcoord-1 + dim_in[0]*(ycoord-1 + dim_in[1]*(zcoord-1)) + offset_in;
          k222 = (double)in[off1]; k122 = (double)in[off1+1]; off2 = off1+dim_in[0];
          k212 = (double)in[off2]; k112 = (double)in[off2+1]; off1+= dim_in[0]*dim_in[1];
          k221 = (double)in[off1]; k121 = (double)in[off1+1]; off2 = off1+dim_in[0];
          k211 = (double)in[off2]; k111 = (double)in[off2+1];

          out[i] = (float)((((k222*dx2 + k122*dx1)*dy2 + (k212*dx2 + k112*dx1)*dy1))*dz2
                         + (((k221*dx2 + k121*dx1)*dy2 + (k211*dx2 + k111*dx1)*dy1))*dz1);
                 
        } else out[i] = 0;
      }
    }
  }
}

/* First order hold resampling - trilinear interpolation */
void subsample_float(float *in, float *out, int dim_in[3], int dim_out[3], int offset_in, int offset_out)
{
  int i, x, y, z;
  double k111,k112,k121,k122,k211,k212,k221,k222;
  double dx1, dx2, dy1, dy2, dz1, dz2, xi, yi, zi, samp[3];
  int off1, off2, xcoord, ycoord, zcoord;
    
  for (i=0; i<3; i++) {
    if(dim_out[i] > dim_in[i]) samp[i] = ceil((double)dim_out[i]/(double)dim_in[i]);
    else                       samp[i] = 1.0/(ceil((double)dim_in[i]/(double)dim_out[i]));
  }
  
  for (z=0; z<dim_out[2]; z++) {
    zi = 1.0+(double)z/samp[2];
    for (y=0; y<dim_out[1]; y++) {
      yi = 1.0+(double)y/samp[1];
      for (x=0; x<dim_out[0]; x++) {
        xi = 1.0+(double)x/samp[0];
        i = z*dim_out[0]*dim_out[1] + y*dim_out[0] + x + offset_out;

        if (zi>=0 && zi<dim_in[2] && yi>=0 && yi<dim_in[1] && xi>=0 && xi<dim_in[0])  {
          xcoord = (int)floor(xi); dx1=xi-(double)xcoord; dx2=1.0-dx1;
          ycoord = (int)floor(yi); dy1=yi-(double)ycoord; dy2=1.0-dy1;
          zcoord = (int)floor(zi); dz1=zi-(double)zcoord; dz2=1.0-dz1;

          off1 = xcoord-1 + dim_in[0]*(ycoord-1 + dim_in[1]*(zcoord-1)) + offset_in;
          k222 = (double)in[off1]; k122 = (double)in[off1+1]; off2 = off1+dim_in[0];
          k212 = (double)in[off2]; k112 = (double)in[off2+1]; off1+= dim_in[0]*dim_in[1];
          k221 = (double)in[off1]; k121 = (double)in[off1+1]; off2 = off1+dim_in[0];
          k211 = (double)in[off2]; k111 = (double)in[off2+1];

          out[i] = (float)((((k222*dx2 + k122*dx1)*dy2 + (k212*dx2 + k112*dx1)*dy1))*dz2
                         + (((k221*dx2 + k121*dx1)*dy2 + (k211*dx2 + k111*dx1)*dy1))*dz1);
                 
        } else out[i] = 0;
      }
    }
  }
}

void
smooth_double(double *vol, int dims[3], double separations[3], double s0[3], int use_mask)
{
  int i;
  double xsum, ysum, zsum;
  double *x, *y, *z, s[3];
  int xyz[3], nvol, sum_mask;
  double *mask;
  unsigned char *mask2;
  
  nvol = dims[0]*dims[1]*dims[2];

  for(i=0; i<3; i++) {
    s[i] = s0[i]/separations[i];
    if(s[i] < 1.0) s[i] = 1.0;
    s[i] /= sqrt(8.0*log(2.0));
    xyz[i] = (int) RINT(6.0*s[i]);
  }
  
  x = (double *) malloc(sizeof(double)*((2*xyz[0])+1));
  y = (double *) malloc(sizeof(double)*((2*xyz[1])+1));
  z = (double *) malloc(sizeof(double)*((2*xyz[2])+1));
  
  /* build mask for masked smoothing */
  if(use_mask) {
    mask  = (double *) malloc(sizeof(double)*nvol);
    mask2 = (unsigned char *) malloc(sizeof(unsigned char)*nvol);
    sum_mask = 0;
    for(i=0; i<nvol; i++) {
      if(vol[i] == 0.0) {
        mask[i]  = 0.0;
        mask2[i] = 0;
      } else {
        mask[i]  = 1.0;
        mask2[i] = 1;
        sum_mask++;
      }
    }
  }
  
  for(i=-xyz[0]; i <= xyz[0]; i++) x[i+xyz[0]] = (double)i;
  for(i=-xyz[1]; i <= xyz[1]; i++) y[i+xyz[1]] = (double)i;
  for(i=-xyz[2]; i <= xyz[2]; i++) z[i+xyz[2]] = (double)i;
  
  xsum = 0.0; ysum = 0.0; zsum = 0.0;
  for(i=0; i < ((2*xyz[0])+1); i++) {
    x[i] = exp(-pow(x[i],2) / (2.0*pow(s[0],2)));
    xsum += x[i];
  }
  for(i=0; i < ((2*xyz[1])+1); i++) {
    y[i] = exp(-pow(y[i],2) / (2.0*pow(s[1],2)));
    ysum += y[i];
  }
  for(i=0; i < ((2*xyz[2])+1); i++) {
    z[i] = exp(-pow(z[i],2) / (2.0*pow(s[2],2)));
    zsum += z[i];
  }
  
  for(i=0; i < ((2*xyz[0])+1); i++) x[i] /= xsum;
  for(i=0; i < ((2*xyz[1])+1); i++) y[i] /= ysum;
  for(i=0; i < ((2*xyz[2])+1); i++) z[i] /= zsum;
  
  convxyz_double(vol,x,y,z,((2*xyz[0])+1),((2*xyz[1])+1),((2*xyz[2])+1),-xyz[0],-xyz[1],-xyz[2],vol,dims);
  if(use_mask) {
    /* only smooth mask if mask has values > 0 */
    if(sum_mask>0) {
      convxyz_double(mask,x,y,z,((2*xyz[0])+1),((2*xyz[1])+1),((2*xyz[2])+1),-xyz[0],-xyz[1],-xyz[2],mask,dims);
      for(i=0; i<nvol; i++) {
        if(mask2[i]>0) vol[i] /= (double)mask[i];   
        else vol[i] = 0.0; 
      }
    }
    free(mask);
    free(mask2);
  }
  
  free(x);
  free(y);
  free(z);

}

void
smooth_float(float *vol, int dims[3], double separations[3], double s0[3], int use_mask)
{
  int i;
  double xsum, ysum, zsum;
  double *x, *y, *z, s[3];
  int xyz[3], nvol, sum_mask;
  float *mask;
  unsigned char *mask2;
  
  nvol = dims[0]*dims[1]*dims[2];

  for(i=0; i<3; i++) {
    s[i] = s0[i]/separations[i];
    if(s[i] < 1.0) s[i] = 1.0;
    s[i] /= sqrt(8.0*log(2.0));
    xyz[i] = (int) RINT(6.0*s[i]);
  }
  
  x = (double *) malloc(sizeof(double)*((2*xyz[0])+1));
  y = (double *) malloc(sizeof(double)*((2*xyz[1])+1));
  z = (double *) malloc(sizeof(double)*((2*xyz[2])+1));
  
  /* build mask for masked smoothing */
  if(use_mask) {
    mask  = (float *) malloc(sizeof(float)*nvol);
    mask2 = (unsigned char *) malloc(sizeof(unsigned char)*nvol);
    sum_mask = 0;
    for(i=0; i<nvol; i++) {
      if(vol[i] == 0.0) {
        mask[i]  = 0.0;
        mask2[i] = 0;
      } else {
        mask[i]  = 1.0;
        mask2[i] = 1;
        sum_mask++;
      }
    }
  }
  
  for(i=-xyz[0]; i <= xyz[0]; i++) x[i+xyz[0]] = (double)i;
  for(i=-xyz[1]; i <= xyz[1]; i++) y[i+xyz[1]] = (double)i;
  for(i=-xyz[2]; i <= xyz[2]; i++) z[i+xyz[2]] = (double)i;
  
  xsum = 0.0; ysum = 0.0; zsum = 0.0;
  for(i=0; i < ((2*xyz[0])+1); i++) {
    x[i] = exp(-pow(x[i],2) / (2.0*pow(s[0],2)));
    xsum += x[i];
  }
  for(i=0; i < ((2*xyz[1])+1); i++) {
    y[i] = exp(-pow(y[i],2) / (2.0*pow(s[1],2)));
    ysum += y[i];
  }
  for(i=0; i < ((2*xyz[2])+1); i++) {
    z[i] = exp(-pow(z[i],2) / (2.0*pow(s[2],2)));
    zsum += z[i];
  }
  
  for(i=0; i < ((2*xyz[0])+1); i++) x[i] /= xsum;
  for(i=0; i < ((2*xyz[1])+1); i++) y[i] /= ysum;
  for(i=0; i < ((2*xyz[2])+1); i++) z[i] /= zsum;
  
  convxyz_float(vol,x,y,z,((2*xyz[0])+1),((2*xyz[1])+1),((2*xyz[2])+1),-xyz[0],-xyz[1],-xyz[2],vol,dims);
  if(use_mask) {
    /* only smooth mask if mask has values > 0 */
    if(sum_mask>0) {
      convxyz_float(mask,x,y,z,((2*xyz[0])+1),((2*xyz[1])+1),((2*xyz[2])+1),-xyz[0],-xyz[1],-xyz[2],mask,dims);
      for(i=0; i<nvol; i++) {
        if(mask2[i]>0) vol[i] /= (float)mask[i];   
        else vol[i] = 0.0; 
      }
    }
    free(mask);
    free(mask2);
  }
  
  free(x);
  free(y);
  free(z);

}

void
smooth_subsample_double(double *vol, int dims[3], double separations[3], double s[3], int use_mask, int samp)
{
  int i, nvol_samp, nvol;
  int dims_samp[3];
  double *vol_samp, separations_samp[3];
  
  /* define grid dimensions */
  for(i=0; i<3; i++) dims_samp[i] = (int) ceil((dims[i]-1)/((double) samp))+1;
  for(i=0; i<3; i++) separations_samp[i] = separations[i]*((double)dims[i]/(double)dims_samp[i]);

  nvol  = dims[0]*dims[1]*dims[2];
  nvol_samp  = dims_samp[0]*dims_samp[1]*dims_samp[2];
  vol_samp  = (double *)malloc(sizeof(double)*nvol_samp);

  subsample_double(vol, vol_samp, dims, dims_samp, 0, 0);    
  smooth_double(vol_samp, dims_samp, separations_samp, s, use_mask);
  subsample_double(vol_samp, vol, dims_samp, dims, 0, 0);    

  free(vol_samp);
}

void
smooth_subsample_float(float *vol, int dims[3], double separations[3], double s[3], int use_mask, int samp)
{
  int i, nvol_samp, nvol;
  int dims_samp[3];
  float *vol_samp;
  double separations_samp[3];
  
  /* define grid dimensions */
  for(i=0; i<3; i++) dims_samp[i] = (int) ceil((dims[i]-1)/((double) samp))+1;
  for(i=0; i<3; i++) separations_samp[i] = separations[i]*((double)dims[i]/(double)dims_samp[i]);

  nvol  = dims[0]*dims[1]*dims[2];
  nvol_samp  = dims_samp[0]*dims_samp[1]*dims_samp[2];
  vol_samp  = (float *)malloc(sizeof(float)*nvol_samp);

  subsample_float(vol, vol_samp, dims, dims_samp, 0, 0);    
  smooth_float(vol_samp, dims_samp, separations_samp, s, use_mask);
  subsample_float(vol_samp, vol, dims_samp, dims, 0, 0);    

  free(vol_samp);
}

void
initial_cleanup(unsigned char *probs, unsigned char *label, int *dims, double *voxelsize, int strength, int remove_sinus)
{
  
  double scale = 3.0/(voxelsize[0] + voxelsize[1] + voxelsize[2]);
  int vol, th, i;
  int n_initial_openings = MAX(1,round(scale*strength));
  float *sum;
  double filt[3] = {0.75, 1.0, 0.75};
  
  vol = dims[0]*dims[1]*dims[2];

  sum = (float *)malloc(sizeof(float)*vol);

  /* build a first rough mask to remove noisy parts */
  for( i = 0;  i < vol;  ++i )
    sum[i] = (float)probs[i + GM*vol] + (float)probs[i + WM*vol];

  morph_open_float(sum, dims, n_initial_openings, 0.1);
  morph_dilate_float(sum, dims, round(scale*1), 0.5);
  distclose_float(sum, dims, voxelsize, round(scale*10), 0.5);

  if(remove_sinus) {
    /* remove sinus sagittalis */
    for (i = 0; i < vol; i++)
      sum[i] = sum[i] && ( label[i] < 4 );
  }

  distclose_float(sum, dims, voxelsize, round(scale*2), 0.5);

  for( i = 0;  i < vol;  ++i )
    label[i] = (unsigned char)sum[i];
  
  free(sum);
}

void
cleanup(unsigned char *probs, unsigned char *mask, int *dims, double *voxelsize, int strength)
{
  
  double scale = 3.0/(voxelsize[0] + voxelsize[1] + voxelsize[2]);
  int niter, iter, vol, th, th_erode, th_dilate, i;
  int n_initial_openings = MAX(1,round(scale*strength));
  float *sum;
  double filt[3] = {0.75, 1.0, 0.75};
  
  niter     =  45;
  th_erode  = 153;   /* initial threshold for erosion 0.6*255.0 */
  th_dilate = (5*strength + 1)*16; /* threshold for dilation */
  
  vol = dims[0]*dims[1]*dims[2];

  sum = (float *)malloc(sizeof(float)*vol);

  /* build a first rough mask to remove noisy parts */
  for( i = 0;  i < vol;  ++i )
    sum[i] = (float)probs[i + GM*vol] + (float)probs[i + WM*vol];

  morph_open_float(sum, dims, n_initial_openings, 0.25);
  
  /* init mask with WM values that are larger than GM and CSF and threshold for erosion */
  for( i = 0;  i < vol;  ++i )
    if ((probs[i + WM*vol] > probs[i + GM*vol]) && (probs[i + WM*vol] > probs[i + CSF*vol]) && (probs[i + WM*vol] > th_erode) && (sum[i] > 0))
      mask[i] = probs[i + WM*vol];
    else mask[i] = 0;

  /* use masked WM image for all subsequent operations */
  for( i = 0;  i < vol;  ++i ) probs[i + WM*vol] = mask[i];
  
  /* mask computed from gm and wm */
  /* erosions and conditional dilations */
  for (iter=0; iter < niter; iter++) {
  
    /*  start with 2 iterations of erosions*/
    if( iter < 2 ) th = th_erode;
    else           th = th_dilate;
    
    /* mask = (mask>th).*(white+gray) */
    for( i = 0;  i < vol;  ++i ) {
      if( mask[i] > th ) {
        sum[i] = (float)probs[i + GM*vol] + (float)probs[i + WM*vol];
        mask[i] = (unsigned char)MIN(sum[i], 255.0);
      } else  mask[i] = 0;		    
    }

    /* convolve mask with filter width of 3 voxel */
    convxyz_uint8(mask,filt,filt,filt,3,3,3,-1,-1,-1,mask,dims);
  }
  
  for( i = 0;  i < vol;  ++i )
    sum[i] = (float)mask[i];

  /* use copy of mask to erode and fill holes */
  morph_erode_float(sum, dims, round(scale*4), 0.5);
  morph_close_float(sum, dims, round(scale*20), 0.5);

  /* use either original mask or new eroded and filled mask */
  for( i = 0;  i < vol;  ++i )
    mask[i] = (sum[i] > 0) ||  (mask[i] > 0);

  /* fill remaining CSF spaces */
  distclose_uint8(mask, dims, voxelsize, round(scale*4), 0.5);

  free(sum);
  
}

/* qicksort */
void
swap_uint8(unsigned char *a, unsigned char *b)
{
  float t=*a; *a=*b; *b=t;
}

void
sort_uint8(unsigned char arr[], int start, int end)
{
  if (end > start + 1)
  {
    unsigned char piv = arr[start];
    int l = start + 1, r = end;
    while (l < r)
    {
      if (arr[l] <= piv) l++;
      else swap_uint8(&arr[l], &arr[--r]);
    }
    swap_uint8(&arr[--l], &arr[start]);
    sort_uint8(arr, start, l);
    sort_uint8(arr, r, end);
  }
}


/* simple median function for uint8 */
void 
median3_uint8(unsigned char *D, int *dims)
{

  /* indices of the neighbor Ni (index distance) and euclidean distance NW */
  unsigned char NV[27];
  int i,j,k,ind,ni,x,y,z,n;
  unsigned char *M;
        
  /* output */
  M = (unsigned char *)malloc(sizeof(unsigned char)*dims[0]*dims[1]*dims[2]);

  /* filter process */
  for (z=0;z<dims[2];z++) for (y=0;y<dims[1];y++) for (x=0;x<dims[0];x++) {
    ind = index(x,y,z,dims);
    n = 0;
    /* go through all elements in a 3x3x3 box */
    for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (k=-1;k<=1;k++) {
      /* check borders */ 
      if ( ((x+i)>=0) && ((x+i)<dims[0]) && ((y+j)>=0) && ((y+j)<dims[1]) && ((z+k)>=0) && ((z+k)<dims[2])) {
        ni = index(x+i,y+j,z+k,dims);
        /* check masks and NaN or Infinities */
        if (isnan(D[ni]) || D[ni]==FLT_MAX || D[ind]==-FLT_MAX ) ni = ind;
        NV[n] = D[ni];
        n++;
      }
    }
    /* get correct n */
    n--;
    /* sort and get the median by finding the element in the middle of the sorting */
    sort_uint8(NV,0,n);
    M[ind] = NV[(int)(n/2)];
  }
   
  for (i=0;i<dims[0]*dims[1]*dims[2];i++) D[i] = M[i];
  
  free(M);
  
}

