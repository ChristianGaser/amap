/*
 * Christian Gaser
 * $Id$ 
 *
 */

#include <ParseArgv.h>
#include "niibrainmask.h"

#include "Amap.h"
#include "niiamap.h"

static ArgvInfo argTable[] = {
  {"-iters_nu", ARGV_INT, (char *) 1, (char *) &iters_nu,
       "Number of iterations for nu correction."},
  {"-strip", ARGV_INT, (char *) 4, (char *) &strip_param,
       "Number of openings1, dilations1, openings2, and dilations2."},
  {"-write_masked", ARGV_CONSTANT, (char *) 1, (char *) &write_masked,
       "Write masked nu corrected image (default)."},
  {"-write_nu", ARGV_CONSTANT, (char *) 1, (char *) &write_nu,
       "Write unmasked nu corrected image."},
  {"-write_brainmask", ARGV_CONSTANT, (char *) 1, (char *) &write_brainmask,
       "Write mask image."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};


static int usage(void)
{
    static const char msg[] = {
        "niiamap: Segmentation with adaptive MAP\n"
        "usage: niiamap [options] in.nii out.nii\n"
    };
    fprintf(stderr, "%s", msg);
    exit(EXIT_FAILURE);
}

int
main( int argc, char **argv )
{
  /* NIFTI stuff */
  nifti_image *src_ptr;
  char      *input_filename, *output_filename, *basename, *extension;
  int       i, j, dims[3], thresh, thresh_kmeans_int;
  int		    n_classes;
  char		  buffer[1024], *str_ptr;
  unsigned char *label, *mask, *prob;
  double	  *src, *filtered, slope;
  double    offset, max_vol, min_vol, voxelsize[3];

  /* Get arguments */
  if (ParseArgv(&argc, argv, argTable, 0) || (argc < 3)) {
    (void) fprintf(stderr, 
    "\nUsage: %s [options] in.nii out.nii\n",
                     argv[0]);
    (void) fprintf(stderr, 
      "       %s -help\n\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  
  input_filename  = argv[1];
  output_filename = argv[2];

  /* deal with extension */
  extension = nifti_find_file_extension(output_filename);
  
  /* if no valid extension was found use .nii */
  if (extension == NULL) {
    fprintf(stderr,"No valid extension found for output filename %s.\n",output_filename);
    exit(EXIT_FAILURE);
  }

  /* read data and scale it to 0..255 */
  strcpy(buffer, input_filename);
  str_ptr = strrchr(buffer, '.');
  src_ptr = read_nifti_double(input_filename, &src);
  if(src_ptr == NULL) {
    fprintf(stderr,"Error reading %s.\n", input_filename);
    return(EXIT_FAILURE);
  }
  
  if(src_ptr == NULL) {
    fprintf(stderr,"Error reading %s.\n",input_filename);
    exit(EXIT_FAILURE);
  }

  mask  = (unsigned char *)malloc(sizeof(unsigned char)*src_ptr->nvox);
  label = (unsigned char *)malloc(sizeof(unsigned char)*src_ptr->nvox);
  filtered = (double *)malloc(sizeof(double)*src_ptr->nvox);
    
  if((mask == NULL) || (label == NULL)) {
    fprintf(stderr,"Memory allocation error\n");
    exit(EXIT_FAILURE);
  }

  /* get min/max */
  min_vol =  FLT_MAX; max_vol = -FLT_MAX;
  for (i = 0; i < src_ptr->nvox; i++) {
    min_vol = MIN(src[i], min_vol);
    max_vol = MAX(src[i], max_vol);
  }

  /* correct images with values < 0 */
  if (min_vol < 0) {
    for (i = 0; i < src_ptr->nvox; i++)
      src[i] = src[i] - min_vol;
    min_vol = 0.0;
  }

  /* use minimum value in the image to get mask value */
  for (i = 0; i < src_ptr->nvox; i++) {
    if (src[i] == min_vol) mask[i] = 0;
    else mask[i] = 255;
  }
     
  /* add offset to ensure that CSF values are much larger than background noise */
  offset = 0.2*max_vol;  
  for (i = 0; i < src_ptr->nvox; i++)
    if (mask[i] > 0) src[i] += offset;

  voxelsize[0] = src_ptr->dx;
  voxelsize[1] = src_ptr->dy;
  voxelsize[2] = src_ptr->dz;
  dims[0] = src_ptr->nx;
  dims[1] = src_ptr->ny;
  dims[2] = src_ptr->nz;
    
  thresh = 128;
  thresh_kmeans_int = 128;
  
  /* initially use 4 classes to consider background */
  n_classes = 4;
  max_vol = Kmeans( src, label, mask, 25, n_classes, voxelsize, dims, thresh, thresh_kmeans_int, 20, NOPVE, bias_fwhm);

  double mu[2];
  int n_voxel[2];

  /* estimate mean for the first two classes */
  for (j = 0; j < 2; j++) {
    n_voxel[j] = 0;
    for (i = 0; i < src_ptr->nvox; i++) {
      if (label[i] == j+1) {
        n_voxel[j]++;
        mu[j] += src[i];
      }
    }
    mu[j] /= (double)n_voxel[j];
  }
  
  /* set mask to zero where image has values for background */
  double mn_thresh = (mu[0]+mu[1])/3.0;
  for (i = 0; i < src_ptr->nvox; i++) {
    if (src[i] < mn_thresh) {
//      mask[i] = 0;
    }
  }
  
fprintf(stderr,"%g\n",(mu[0]+mu[1])/2.0);
  ornlm(src, filtered, 3, 1, (mu[0]+mu[1])/2.0, dims);

  /* second Kmeans with 3 classes */
  n_classes = 3;
  max_vol = Kmeans( src, label, mask, 25, n_classes, voxelsize, dims, thresh, thresh_kmeans_int, 20, KMEANS, bias_fwhm);    

  /* first rough skull-stripping */  
fprintf(stderr,".");
  morph_open_uint8(label, dims, strip_param[0], 3);

fprintf(stderr,".");
  get_largest_cluster(label, dims);
    for (i = 0; i < src_ptr->nvox; i++)
      src[i] = (double)label[i];

    (void) sprintf( buffer, "%s_brainmask%s",basename,extension); 
    
    if(!write_nifti(buffer, src, DT_UINT8, slope, 
            dims, voxelsize, src_ptr))
      exit(EXIT_FAILURE);
fprintf(stderr,":");
  morph_dilate_uint8(label, dims, strip_param[1], 0);
fprintf(stderr,".");
  morph_close_uint8(label, dims, 10, 0);
fprintf(stderr,".");
  
  /* update mask */
  for (i = 0; i < src_ptr->nvox; i++)
    mask[i] = 255*(label[i] > 0);

  n_classes = 5;
  double mean[n_classes];
  prob  = (unsigned char *)malloc(sizeof(unsigned char)*src_ptr->nvox*n_classes);

  /* use Kmeans with 5 classes */
  max_vol = Kmeans( filtered, label, mask, 25, 3, voxelsize, dims, thresh, thresh_kmeans_int, 10, NOPVE, bias_fwhm);

  /* use amap approach with PVE */
fprintf(stderr,".");
int iters_icm = 50;
  Amap( filtered, label, prob, mean, 3, 10, subsample, dims, 1, weight_MRF, voxelsize, iters_icm, offset);
fprintf(stderr,".");
  Pve5(filtered, prob, label, mean, dims);

fprintf(stderr,".");

  /* final skull-stripping */
  morph_open_uint8(label, dims, strip_param[2], 160);
fprintf(stderr,".");
  get_largest_cluster(label, dims);
fprintf(stderr,".");
  morph_dilate_uint8(label, dims, strip_param[3], 0);
fprintf(stderr,".");
  morph_close_uint8(label, dims, 10, 0);
fprintf(stderr,".");

  basename = nifti_makebasename(output_filename);

  slope = 1.0;

  /* write nu corrected image */
  if (write_nu) {
     (void) sprintf( buffer, "%s_unmasked%s",basename,extension); 

    if(!write_nifti( buffer, filtered, DT_FLOAT32, slope, dims, 
            voxelsize, src_ptr))
      exit(EXIT_FAILURE);
  }
  
  /* write masked image */
  if (write_masked) {

    for (i = 0; i < src_ptr->nvox; i++)
      src[i] = (double)label[i]*filtered[i];

    if(!write_nifti( output_filename, src, DT_FLOAT32, slope, dims, 
            voxelsize, src_ptr))
      exit(EXIT_FAILURE);
  }

  /* write brainmask */
  if (write_brainmask) {

    for (i = 0; i < src_ptr->nvox; i++)
      src[i] = (double)label[i];

    (void) sprintf( buffer, "%s_brainmask%s",basename,extension); 
    
    if(!write_nifti(buffer, src, DT_UINT8, slope, 
            dims, voxelsize, src_ptr))
      exit(EXIT_FAILURE);
    
  }
    
  free(src);
  free(label);
  free(mask);
  free(prob);
  
  return(EXIT_SUCCESS);
}