/*
 * Christian Gaser
 * $Id$ 
 *
 */

#include <stdlib.h>
#include <math.h>
#include <memory.h>
#define RINT(A) floor((A)+0.5)

static void 
convxy(out, xdim, ydim, filtx, filty, fxdim, fydim, xoff, yoff, buff)
int xdim, ydim, fxdim, fydim, xoff, yoff;
double out[], filtx[], filty[], buff[];
{
  int x,y,k;
  for(y=0; y<ydim; y++)
  {
    for(x=0; x<xdim; x++)
    {
      buff[x] = out[x+y*xdim];
      if (!finite(buff[x]))
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


int 
convxyz_double(double *iVol, double filtx[], double filty[], double filtz[],
  int fxdim, int fydim, int fzdim, int xoff, int yoff, int zoff,
  double *oVol, int dims[3])
{
  double *tmp, *buff, **sortedv;
  int xy, z, y, x, k, fstart, fend, startz, endz;
  int xdim, ydim, zdim;

  xdim = dims[0];
  ydim = dims[1];
  zdim = dims[2];

  tmp = (double *)malloc(sizeof(double)*xdim*ydim*fzdim);
  buff = (double *)malloc(sizeof(double)*((ydim>xdim) ? ydim : xdim));
  sortedv = (double **)malloc(sizeof(double *)*fzdim);


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

      double *obuf;
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



int convxyz_uint8(unsigned char *iVol, double filtx[], double filty[], double filtz[],
  int fxdim, int fydim, int fzdim, int xoff, int yoff, int zoff,
  unsigned char *oVol, int dims[3])
{
  double *tmp, *buff, **sortedv;
  int xy, z, y, x, k, fstart, fend, startz, endz;
  int xdim, ydim, zdim;

  xdim = dims[0];
  ydim = dims[1];
  zdim = dims[2];

  tmp = (double *)malloc(sizeof(double)*xdim*ydim*fzdim);
  buff = (double *)malloc(sizeof(double)*((ydim>xdim) ? ydim : xdim));
  sortedv = (double **)malloc(sizeof(double *)*fzdim);


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

      double tmp;
      unsigned char *obuf;
      obuf = oVol;
      obuf = &obuf[(z-fzdim-zoff+1)*ydim*xdim];
      if (sum2)
      {
        for(xy=0; xy<xdim*ydim; xy++)
        {
          double sum1=0.0;
          for(k=fstart; k<fend; k++)
            sum1 += filtz[k]*sortedv[k][xy];
          tmp = sum1/sum2;
          if (tmp<0) tmp = 0;
          else if (tmp>255) tmp = 255;
          obuf[xy] = RINT(tmp);
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

int convxyz_int16(signed short *iVol, double filtx[], double filty[], double filtz[],
  int fxdim, int fydim, int fzdim, int xoff, int yoff, int zoff,
  signed short *oVol, int dims[3])
{
  double *tmp, *buff, **sortedv;
  int xy, z, y, x, k, fstart, fend, startz, endz;
  int xdim, ydim, zdim;

  xdim = dims[0];
  ydim = dims[1];
  zdim = dims[2];

  tmp = (double *)malloc(sizeof(double)*xdim*ydim*fzdim);
  buff = (double *)malloc(sizeof(double)*((ydim>xdim) ? ydim : xdim));
  sortedv = (double **)malloc(sizeof(double *)*fzdim);


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

      double tmp;
      signed short *obuf;
      obuf = oVol;
      obuf = &obuf[(z-fzdim-zoff+1)*ydim*xdim];
      if (sum2)
      {
        for(xy=0; xy<xdim*ydim; xy++)
        {
          double sum1=0.0;
          for(k=fstart; k<fend; k++)
            sum1 += filtz[k]*sortedv[k][xy];
          tmp = sum1/sum2;
          if (tmp<-32768) tmp = -32768;
          else if (tmp>32767) tmp = 32767;
          obuf[xy] = RINT(tmp);
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

int convxyz_int32(signed int *iVol, double filtx[], double filty[], double filtz[],
  int fxdim, int fydim, int fzdim, int xoff, int yoff, int zoff,
  signed int *oVol, int dims[3])
{
  double *tmp, *buff, **sortedv;
  int xy, z, y, x, k, fstart, fend, startz, endz;
  int xdim, ydim, zdim;

  xdim = dims[0];
  ydim = dims[1];
  zdim = dims[2];

  tmp = (double *)malloc(sizeof(double)*xdim*ydim*fzdim);
  buff = (double *)malloc(sizeof(double)*((ydim>xdim) ? ydim : xdim));
  sortedv = (double **)malloc(sizeof(double *)*fzdim);


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

      double tmp;
      signed int *obuf;
      obuf = oVol;
      obuf = &obuf[(z-fzdim-zoff+1)*ydim*xdim];
      if (sum2)
      {
        for(xy=0; xy<xdim*ydim; xy++)
        {
          double sum1=0.0;
          for(k=fstart; k<fend; k++)
            sum1 += filtz[k]*sortedv[k][xy];
          tmp = sum1/sum2;
          if (tmp<-2147483648.0) tmp = -2147483648.0;
          else if (tmp>2147483647.0) tmp = 2147483647.0;
          obuf[xy] = RINT(tmp);
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

void
morph_erode_uint8(unsigned char *vol, int dims[3], int niter, int th)
{
  double filt[3]={1,1,1};
  int i,j;
  
  /* threshold input */
  for (j=0;j<dims[2]*dims[1]*dims[0];j++)
    vol[j] = (vol[j]>th);

  for (i=0;i<niter;i++) {
    convxyz_uint8(vol,filt,filt,filt,3,3,3,-1,-1,-1,vol,dims);
    for (j=0;j<dims[2]*dims[1]*dims[0];j++)
      vol[j] = (vol[j]>=9);
  }
}


void
morph_dilate_uint8(unsigned char *vol, int dims[3], int niter, int th)
{
  double filt[3]={1,1,1};
  int i,x,y,z,j,band,dims2[3];
  unsigned char *buffer;

  /* add band with zeros to image to avoid clipping */  
  band = niter;
//  band = 0;
  for (i=0;i<3;i++) dims2[i] = dims[i] + 2*band;
  
  buffer = (unsigned char *)malloc(sizeof(unsigned char)*dims2[0]*dims2[1]*dims2[2]);
  memset(buffer,0,sizeof(unsigned char)*dims2[0]*dims2[1]*dims2[2]);
  
  /* threshold input */
  for (z=0;z<dims[2];z++) for (y=0;y<dims[1];y++) for (x=0;x<dims[0];x++) 
    buffer[((z+band)*dims2[0]*dims2[1])+((y+band)*dims2[0])+x+band] = (vol[(z*dims[0]*dims[1])+(y*dims[0])+x]>th);

  for (i=0;i<niter;i++) {
    convxyz_uint8(buffer,filt,filt,filt,3,3,3,-1,-1,-1,buffer,dims);
    for (j=0;j<dims2[2]*dims2[1]*dims2[0];j++)
      buffer[j] = (buffer[j]>0);
  }
  
  /* return image */
  for (z=0;z<dims[2];z++) for (y=0;y<dims[1];y++) for (x=0;x<dims[0];x++) 
    vol[(z*dims[0]*dims[1])+(y*dims[0])+x] = buffer[((z+band)*dims[0]*dims[1])+((y+band)*dims[0])+x+band];
  
  free(buffer);
}

void
morph_close_uint8(unsigned char *vol, int dims[3], int niter, int th)
{
  morph_dilate_uint8(vol, dims, niter, th);
  morph_erode_uint8(vol, dims, niter, 0);  
}

void
morph_open_uint8(unsigned char *vol, int dims[3], int niter, int th)
{
  morph_erode_uint8(vol, dims, niter, th);
  morph_dilate_uint8(vol, dims, niter, 0);
}

void
morph_close_double(double *vol, int dims[3], int niter, double th)
{
  unsigned char *buffer;
  int i;

  buffer = (unsigned char *)malloc(sizeof(unsigned char)*dims[0]*dims[1]*dims[2]);

  for (i=0;i<dims[2]*dims[1]*dims[0];i++)
    buffer[i] = (unsigned char) (vol[i] > th);
        
  morph_dilate_uint8(buffer, dims, niter, 0);
  morph_erode_uint8(buffer, dims, niter, 0);

  for (i=0;i<dims[2]*dims[1]*dims[0];i++)
    vol[i] = (double)buffer[i];
}

void
morph_open_double(double *vol, int dims[3], int niter, double th)
{
  unsigned char *buffer;
  int i;

  buffer = (unsigned char *)malloc(sizeof(unsigned char)*dims[0]*dims[1]*dims[2]);

  for (i=0;i<dims[2]*dims[1]*dims[0];i++)
    buffer[i] = (unsigned char) (vol[i] > th);
        
  morph_erode_uint8(buffer, dims, niter, 0);
  morph_dilate_uint8(buffer, dims, niter, 0);

  for (i=0;i<dims[2]*dims[1]*dims[0];i++)
    vol[i] = (double)buffer[i];
}
