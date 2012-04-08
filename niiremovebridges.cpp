/*
 * Christian Gaser
 * $Id$ 
 *
 */

#include "CRemoveBridges.h"

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

#ifndef ROUND
#define ROUND( x ) ((long) ((x) + ( ((x) >= 0) ? 0.5 : (-0.5) ) ))
#endif

#ifndef HUGE
#define HUGE 1e15
#endif

#define UNCHANGED  225
#define REMOVED 236
#define ADDED 245
#define ALTERNATIVE 226

extern "C"
{
#include <ParseArgv.h>
#include "nifti/nifti1_io.h"
#include "nifti/nifti1_local.h"
extern nifti_image *read_nifti_double( const char *input_filename, double *image[]);
extern int equal_image_dimensions(nifti_image *nii_ptr, nifti_image *nii_ptr2);
extern int write_nifti_double( const char *output_filename, double image[], int data_type, double slope, int dim[], double vox[], nifti_image *in_ptr);
}

double labelvalue[2] = {-1.0, 127.0};
char *changed_filename = NULL;

static ArgvInfo argTable[] = {
  {"-label", ARGV_FLOAT, (char *) 2, (char *) &labelvalue,
       "range of label to binarize image."},
  {"-changed", ARGV_STRING, (char *) 1, (char *) &changed_filename, 
       "Save image with indicated changes (unchanged 255; removed 236; added 245)."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

double EstimateKmeans(float *src, int n_classes, double *mean, int ni, int *dims, double max_src)
/* perform k-means algorithm give initial mean estimates */
{
  int i, j, j0, x, y, z, v;
  int count;
  long histo[256], lut[256], cumsum[256], vol, area;
  long z_area, y_dims;
  double diff, dmin, dx, xnorm, sum, threshold;
  unsigned char *label;

  area = dims[0]*dims[1];
  vol  = area*dims[2];

  label = (unsigned char *)malloc(sizeof(unsigned char)*vol);

  /* build intensity histogram */
  for (i = 0; i < 256; i++) histo[i] = 0;
  for (i = 0; i < vol; i++) {
    v = (int)ROUND(255.0*src[i]/max_src);
    if (v < 1) continue;
    if (v < 0) v = 0;
    if (v > 255) v = 255;	
    histo[v]++;  
  }

  /* use only value in histogram where cumsum is between 1..99% */
  cumsum[0] = histo[0];
  for (i = 1; i < 256; i++) cumsum[i] = cumsum[i-1] + histo[i];
  for (i = 0; i < 256; i++) cumsum[i] = (long) ROUND(1000.0*(double)cumsum[i]/(double)cumsum[255]);
  for (i = 0; i < 256; i++) if ((cumsum[i] <= 10) || (cumsum[i] >= 990)) histo[i] = 0;

  /* loop through */
  diff = HUGE;  count = 0;
  while (diff > 1.0 && count < ni) {

    /* assign class labels */
    for (i = 0; i < 256; i++) {
      dmin = 256.0 * 256.0;
      for (j = 0; j < n_classes; j++) {
	      dx = (double) i - mean[j];
	      dx *= dx;
	      if (dx < dmin) {
	        lut[i] = j;
	        dmin = dx;
	      }
      }
    }

    /* find the new cluster centers */
    diff = 0;
    for (i = 0; i < n_classes; i++) {
      xnorm = 0.0;
      sum = 0.0;
      for (j = 0; j < 256; j++)
	    if (lut[j] == i) {
	      xnorm += histo[j];
	      sum +=  j * histo[j];
	    }
      sum = xnorm > 0 ? sum /= xnorm : 0.0;
      dx = sum - mean[i];
      mean[i] = sum;
      dx *= dx;
      diff += dx;
    }
    count++;
  }

  /* assign final labels to voxels */
  for (i = 0; i < 256; i++) {
    dmin = HUGE;
    j0 = 0;
    for (j = 0; j < n_classes; j++) {
      if (fabs((double) i - mean[j]) < dmin) {
	      dmin = fabs((double)i - mean[j]);
	      j0 = j;
      }
    }
    lut[i] = j0;
  }
  
  lut[0] = 0;

  /* adjust for the background label */
  diff = 0;
  
  for (i = 0; i < vol; i++) {
    v = (int)ROUND(255.0*src[i]/max_src);
    if (v >= 1) {
      if (v < 0) v = 0;
      if (v > 255) v = 255;
      label[i] = (unsigned char)(lut[v] + 1);	
      diff += SQR((double)v - mean[lut[v]]);
    }
    else label[i] = 0;	
  }

  free(label);
  /* return square error */
  return(diff);
}

void Kmeans(float *src, int NI, int n_clusters, int *dims, double *mu)
{
  int i, j, l, k, x, y, z;
  double e, emin, eps;
  double max_src = -HUGE;

  long n[3];
  double mean[3];
  double var[3];
  double Mu[3];

  int val, n_classes;
  long vol, area, z_area, y_dims;

  area = dims[0]*dims[1];
  vol  = area*dims[2];
  
  int n_classes_initial = n_clusters;

  /* find maximum and mean inside mask */
  for (i = 0; i < vol; i++) {
    if (src[i]>0)
      max_src = MAX(src[i], max_src);
  }

  /* go through all sizes of cluster beginning with two clusters */
  for (n_classes=2; n_classes <= n_classes_initial; n_classes++) {

    if (n_classes == 2) {
      /* initialize for the two cluster case; */
      n[0]=0; mean[0] = 0.0; var[0] = 0.0;

      for (i = 0; i < vol; i++) {
        val = 255.0*src[i]/max_src;
        if (val < 1.0/255.0) continue;
        n[0]++;
        mean[0] += val;
        var[0]  += SQR(val);
      }
      
      Mu[0] = n[0] != 0 ? mean[0]/n[0]: 0.0;
      var[0] = n[0] > 1 ? (var[0] - n[0]*Mu[0]*Mu[0])/(n[0] - 1.0) : 1.0;
      eps = 0.5*sqrt(var[0]);
    }
    else {
      /* find the deviant (epsilon) for the node being divided */
      eps = Mu[0];
      for (i = 0; i < n_classes-2; i++)
        if (Mu[i+1] - Mu[i] < eps)
          eps = Mu[i+1] - Mu[i];
      if (255 - Mu[n_classes-2] < eps)
        eps = 255 - Mu[n_classes-2];
      eps = eps*0.5;
    }

    /* go through low order clustering */
    emin = HUGE;
    for (k = 0; k < n_classes-1; k++) {
      for (i = n_classes-1; i > k+1; i--) mean[i] = Mu[i-1];
      mean[k+1] = Mu[k] + eps;  mean[k] = Mu[k] - eps;
      for (i = 1; i < k; i++) mean[i] = Mu[i];
      e = EstimateKmeans(src, n_classes, mean, NI, dims, max_src);
      if (e < emin) {
        emin = e;
        for (i = 0; i < n_classes; i++) 
          mu[i] = mean[i];
      }
    }
    for (i = 0; i < n_classes; i++) Mu[i] = mu[i];     
  }

  e = EstimateKmeans(src, n_clusters, mu, NI, dims, max_src);
  
  fprintf(stderr,"Final means: ");
  for (i=0; i<n_clusters; i++) 
    fprintf(stderr,"%3.1f ",mu[i]); 
  fprintf(stderr,"\n");    
  return;    
}

int  main(
  int   argc,
  char  *argv[] )
{
  nifti_image *wm_ptr, *t1_ptr;
  int       nii_dimids[MAX_NII_DIMS];
  int       nii_dir[MAX_NII_DIMS];
  int       nii_map[MAX_NII_DIMS];
  unsigned long nii_lens[MAX_NII_DIMS];
  int       nii_ndims, n_classes;
  int       nifti_file_type;
  int	    x, y, z, dims[3], i, z_area, y_dims;
  long	    area, vol;
  char	    *arg_string, *extension;
  unsigned char	*wm, *wm_out;
  char	    *wm_filename, *t1_filename, *wm_out_filename;
  float    *t1;
  double    *vol_tmp, mu[3], voxelsize[3];

  /* Get arguments */
  if (ParseArgv(&argc, argv, argTable, 0) || (argc < 3)) {
    (void) fprintf(stderr, 
    "\nUsage: %s wm.nii t1.nii wm_corrected.nii\n",
                     argv[0]);
    (void) fprintf(stderr, 
      "       %s -help\n\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  
  wm_filename  = argv[1];
  t1_filename  = argv[2];
  wm_out_filename  = argv[3];
    
  /* deal with extension */
  extension = nifti_find_file_extension(wm_out_filename);
  
  /* if no valid extension was found use .nii */
  if (extension == NULL) {
    fprintf(stdout,"Use .nii as extension for %s.\n",wm_out_filename);
    (void) sprintf( extension, ".nii");;
  }

  /* read data */
  wm_ptr = read_nifti_double(wm_filename, &vol_tmp);
  
  if(wm_ptr == NULL) {
    fprintf(stderr,"Error reading %s.\n", wm_filename);
    return(EXIT_FAILURE);
  }
  
  wm  = (unsigned char *)malloc(sizeof(unsigned char)*wm_ptr->nvox);
  wm_out  = (unsigned char *)malloc(sizeof(unsigned char)*wm_ptr->nvox);
  t1 = (float *)malloc(sizeof(float)*wm_ptr->nvox);

  if((wm == NULL) || (wm_out == NULL) || (t1 == NULL)) {
    fprintf(stderr,"Memory allocation error\n");
    exit(EXIT_FAILURE);
  }

  voxelsize[0] = wm_ptr->dx;
  voxelsize[1] = wm_ptr->dy;
  voxelsize[2] = wm_ptr->dz;
  dims[0] = wm_ptr->nx;
  dims[1] = wm_ptr->ny;
  dims[2] = wm_ptr->nz;

  area = dims[0]*dims[1];
  vol  = area*dims[2];

  /* load wm image as uint8 and find bounding box */
  int xmin = 255, xmax = 0, ymin = 255, ymax = 0, zmin = 255, zmax = 0; 
  for (z = 0; z < dims[2]; z++) {
    z_area = z*area;
    for (y = 0; y < dims[1]; y++) {
      y_dims = y*dims[0];
      for (x = 0; x < dims[0]; x++) {
        i = z_area + y_dims + x;
        wm[i] = 0;
        if ((vol_tmp[i] > labelvalue[0]) && (vol_tmp[i] < labelvalue[1]))
          wm[i] = 255;
        else wm[i] = 0;
        /* find bounding box of wm image */
        if (wm[i] > 0) {
          xmax = (x > xmax) ? x : xmax;
          ymax = (y > ymax) ? y : ymax;
          zmax = (z > zmax) ? z : zmax;
          xmin = (x < xmin) ? x : xmin;
          ymin = (y < ymin) ? y : ymin;
          zmin = (z < zmin) ? z : zmin;
        }
      }
    }
  }
  
  xmin++; xmax++; ymin++; ymax++; zmin++; zmax++;
  
  fprintf(stderr,"Bounding box: x: %d-%d y: %d-%d z: %d-%d\n",xmin,xmax,ymin,ymax,zmin,zmax);

  /* read t1 */
  t1_ptr = read_nifti_double(t1_filename, &vol_tmp);
  if(t1_ptr == NULL) {
    fprintf(stderr,"Error reading %s.\n", t1_filename);
    return(EXIT_FAILURE);
  }
    
  /* check size */ 
  if (!equal_image_dimensions(wm_ptr,t1_ptr)) {   
    fprintf(stderr,"WM and T1 image have different sizes\n");
    exit(EXIT_FAILURE);
  }

  /* convert t1 inside bounding box */
  for (i = 0; i < vol; i++) t1[i] = 0;
  for (z = zmin; z < zmax; z++) {
    z_area = z*area;
    for (y = ymin; y < ymax; y++) {
      y_dims = y*dims[0];
      for (x = xmin; x < xmax; x++) {
        i = z_area + y_dims + x;
        t1[i] = (float)vol_tmp[i];
      }
    }
  }
  
  /* use K-Means to estimate averages for GM/WM */
  Kmeans(t1, 25, 3, dims, mu);
  
  CRemoveBridges *bro = new CRemoveBridges();

  /* avgWhite,threshold,avgGray */
  bro->remove(wm, t1, wm_out, dims[0],dims[1],dims[2],mu[2],((mu[2]+mu[1])/2.0),mu[1]);

  /* save file with indicated changes */
  if (changed_filename != NULL) {
    for (i = 0; i < vol; i++) {
      vol_tmp[i] = (double)wm_out[i];
      if (vol_tmp[i] == UNCHANGED) vol_tmp[i] = 255;
    }
    if(!write_nifti_double( changed_filename, vol_tmp, DT_UINT8, 1, dims, 
          voxelsize, wm_ptr))
      exit(EXIT_FAILURE);
  }
  
  int n_changes = 0;
  for (i = 0; i < vol; i++) {
  	if((wm_out[i] == ADDED) || (wm_out[i] == ALTERNATIVE)) {
  		wm_out[i] = 255;
  		n_changes++;
  	}
  	if(wm_out[i] == UNCHANGED)
  		wm_out[i] = 255;
  	if(wm_out[i] == REMOVED) {
  		wm_out[i] = 0;
  		n_changes++;
  	}
  }
  fprintf(stderr,"%d voxels changed\n",n_changes);
  
  for (i = 0; i < vol; i++) vol_tmp[i] = (double)wm_out[i];
  if(!write_nifti_double( wm_out_filename, vol_tmp, DT_UINT8, 1, dims, 
          voxelsize, wm_ptr))
    exit(EXIT_FAILURE);

  free(wm);
  free(t1);
  free(wm_out);
  free(vol_tmp);
  
  return( 0 );
}
