#include <math.h>
#include <memory.h>
#include <stdlib.h> 
#include <stdio.h>  
#include "Amap.h"

#define CONN 18; //connectivity scheme (6 or 18)


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//% author: Keith Forbes     %
//% e-mail: keith@umpire.com %
//% tel: +27 21 674 3345     %
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int inim(int x,int y,int z,const int *dims)
{
 return ((x>=0)&&(x<dims[0]))&&((y>=0)&&(y<dims[1]))&&((z>=0)&&(z<dims[2]));
}


void dilate(unsigned char *out,unsigned char *img,const int *dims)
{
  int x,y,z,i,j,k,area;
  double val;
  
  area = dims[1]*dims[0];

  for (z=1;z<dims[2]-1;z++)
    for (y=1;y<dims[1]-1;y++)
      for (x=1;x<dims[0]-1;x++)
      {
          val = 0;
          // sum of all 26 neighbours
          for (i=-1;i<=1;i++)
            for (j=-1;j<=1;j++)
              for (k=-1;k<=1;k++)
                val += (double)img[((z+i)*area)+((y+j)*dims[0])+x+k];

        out[(z*area)+(y*dims[0])+x] = (unsigned char)(val>0);
      }
}

void mag_gradient(unsigned char *mag,unsigned char *img,const int *dims)
{
  int x,y,z,loop;
  int ind, z_area, y_dims;
  double *tmp, gx, gy, gz;

  int area = dims[1]*dims[0];
  int vol = area*dims[2];
  tmp = (double *) malloc(sizeof(double)*vol);

  for (z=1;z<dims[2]-1;z++)
  {
    z_area = z*area;
    for (y=1;y<dims[1]-1;y++)
    {
      y_dims = y*dims[0];
      for (x=1;x<dims[0]-1;x++)
      {
          ind = z_area + y_dims + x;
        gx = ((double)img[ind+1] - (double)img[ind-1])/2;
        gy = ((double)img[z_area+((y+1)*dims[0])+x] - (double)img[z_area+((y-1)*dims[0])+x])/2;
        gz = ((double)img[((z+1)*area)+y_dims+x] - (double)img[((z-1)*area)+y_dims+x])/2;
        tmp[ind] = sqrt(gx*gx + gy*gy + gz*gz);
      }
    }
  }
  double gmin =  1e15;
  double gmax = -1e15;
  for (loop=0;loop<vol;loop++)
  {
    gmax = MAX(tmp[loop], gmax);  //find max of gradient
    gmin = MIN(tmp[loop], gmin);  //find min of gradient
  }
  
  // scale data to range of 0..255
  double scale = (gmax - gmin)/255;
  for (loop=0;loop<vol;loop++)
    mag[loop] = round((tmp[loop] - gmin)/scale);

  free(tmp);
}

void watershed(unsigned char *imgout,unsigned char *g,unsigned char *marker,const int *dims
) //form the watershed transform from a gradient image and a marker image
{

  long count[256];
  long startq[256];
  long endq[256];
  int * xq[256];
  int * yq[256];
  int * zq[256];
  int * temp;
  int mmax;
  int x,y,z,area;
  int z_area, y_dims;
  int currentq=0;
  int xcurr, ycurr, zcurr;
  int ixcurr, iycurr, izcurr;
  int ind, ind0;
  int qnum, conn;
  long loop;
  long numelements;
  int i;
  static int xi[] = { 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0};
  static int yi[] = { 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1};
  static int zi[] = { 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0,-1, 1,-1, 1};

  conn = CONN;
  numelements=(long)(dims[2]*dims[1]*dims[0]);

  for (i=0;i<=255;i++)
    count[i]=0;

  mmax = -1e5;
  for (loop=0;loop<numelements;loop++)
  {
    imgout[loop]=marker[loop];  //copy the marker image into the output image
    mmax = MAX(marker[loop], mmax);  //find max of marker
  }

  for (loop=0;loop<numelements;loop++)
  {
    count[(int)g[loop]]+=1;
    count[0]+=(int)(imgout[loop]>0); //count marked pixels as zeros
  }   

  for (i=0;i<=255;i++)
  {
    xq[i] = malloc(sizeof(int)*count[i]); //create a new space
    yq[i] = malloc(sizeof(int)*count[i]); //create a new space  
    zq[i] = malloc(sizeof(int)*count[i]); //create a new space
    startq[i] =  0;
    endq[i]   = -1;
  }

  area = dims[1]*dims[0];
  
  for (z=1;z<dims[2];z++)
  {
    z_area = z*area;
    for (y=1;y<dims[1];y++)
    {
      y_dims = y*dims[0];
      for (x=1;x<dims[0];x++)
      {
        if ( imgout[z_area+y_dims+x] > 0) //if marked flood sources
        {
          if (endq[0]+1==count[0]) //if out of queue space already
          {
            count[0]*= 2; //double the number of elements available
          
            temp = malloc(sizeof(int)*count[0]); //create a new space
            memcpy(temp,xq[0],(endq[0]+1)*sizeof(int)); //copy existing data across
            free(xq[0]);
            xq[0]=temp; //point to the new larger array

            temp = malloc(sizeof(int)*count[0]); //create a new space       
            memcpy(temp,yq[0],(endq[0]+1)*sizeof(int)); //copy existing data across
            free(yq[0]);
            yq[0]=temp; //point to the new larger array

            temp = malloc(sizeof(int)*count[0]); //create a new space       
            memcpy(temp,zq[0],(endq[0]+1)*sizeof(int)); //copy existing data across
            free(zq[0]);
            zq[0]=temp; //point to the new larger array
          }
        
          endq[0]++;  //marker pixels are dealt with first
          (xq[0])[endq[0]]=x;
          (yq[0])[endq[0]]=y;
          (zq[0])[endq[0]]=z; 
        }
      }
    }
  }

  //MAIN WATERSHED LOOP

  while (currentq<=255)
  {
    if (endq[currentq]>=startq[currentq])
    {
      xcurr=(xq[currentq])[startq[currentq]];
      ycurr=(yq[currentq])[startq[currentq]];
      zcurr=(zq[currentq])[startq[currentq]];
      startq[currentq]++;

      ind0 = xcurr+(ycurr*dims[0])+(zcurr*area);
          
      for (i=0;i<conn;i++)
      {
        ixcurr = xcurr+xi[i];
        iycurr = ycurr+yi[i];
        izcurr = zcurr+zi[i];

        if (inim(ixcurr,iycurr,izcurr,dims))
        {
          ind = ixcurr+(iycurr*dims[0])+(izcurr*area);
          if (imgout[ind]==0) //not yet dealt with
          {
            imgout[ind]=imgout[ind0];

            qnum=(int)g[ind];
            if (qnum<currentq)
               qnum=currentq; //%if q is closed put on current q

            if (endq[qnum]+1==count[qnum]) //if out of queue space 
            {
              count[qnum]*=2; //double the number of elements available
                
              temp = malloc(sizeof(int)*count[qnum]); //create a new space
              memcpy(temp,xq[qnum],(endq[qnum]+1)*sizeof(int)); //copy existing data across
              free(xq[qnum]);
              xq[qnum]=temp; //point to the new larger array

              temp = malloc(sizeof(int)*count[qnum]); //create a new space
              memcpy(temp,yq[qnum],(endq[qnum]+1)*sizeof(int)); //copy existing data across
              free(yq[qnum]);
              yq[qnum]=temp; //point to the new larger array

              temp = malloc(sizeof(int)*count[qnum]); //create a new space
              memcpy(temp,zq[qnum],(endq[qnum]+1)*sizeof(int)); //copy existing data across
              free(zq[qnum]);
              zq[qnum]=temp; //point to the new larger array
            }
              
            endq[qnum]++;
            (xq[qnum])[endq[qnum]]=ixcurr;
            (yq[qnum])[endq[qnum]]=iycurr;
            (zq[qnum])[endq[qnum]]=izcurr;
          }
        }
      }
    }

    if (startq[currentq]>endq[currentq])
      currentq++;

  }
  
  // set maximum value of marker to 0
  for (loop=0;loop<numelements;loop++)
    if (imgout[loop]==(unsigned char)mmax)
      imgout[loop] = 0;
  
  for (loop=0;loop<=255;loop++)
  {
    free(xq[loop]); //clean up by freeing used space otherwise memory leak is created
    free(yq[loop]);
    free(zq[loop]);
  }
  
}

void watershed3d(unsigned char *img, unsigned char *marker, int flag_dilate, int *dims)
{
  
  unsigned char *g = (unsigned char *)malloc(sizeof(unsigned char)*dims[0]*dims[1]*dims[2]);
  unsigned char *d = (unsigned char *)malloc(sizeof(unsigned char)*dims[0]*dims[1]*dims[2]);

  mag_gradient(g, img, dims);
  if (flag_dilate != 0) {
    watershed(d, g, marker, dims);
    free(g);
    dilate(img, d, dims);
  }
  else watershed(img, g, marker, dims); 

  free(d);
}