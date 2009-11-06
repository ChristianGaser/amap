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
  int		n_classes;
  char		buffer[1024], *str_ptr;
  unsigned char *label, *mask, *prob;
  double	*src, slope;
  double    max_vol, min_vol, separations[3];

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
  src_ptr = read_nifti_float(input_filename, &src);
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
     
  separations[0] = src_ptr->dx;
  separations[1] = src_ptr->dy;
  separations[2] = src_ptr->dz;
  dims[0] = src_ptr->nx;
  dims[1] = src_ptr->ny;
  dims[2] = src_ptr->nz;
    
  thresh = 128;
  thresh_kmeans_int = 128;
  
  /* initially use 4 classes to consider background */
  n_classes = 4;
  max_vol = Kmeans( src, label, mask, 25, n_classes, separations, dims, thresh, thresh_kmeans_int, iters_nu, NOPVE, 50.0);

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
      mask[i] = 0;
    }
  }
  
  /* second Kmeans with 3 classes */
  n_classes = 3;
  max_vol = Kmeans( src, label, mask, 25, n_classes, separations, dims, thresh, thresh_kmeans_int, iters_nu, KMEANS, 50.0);

  /* first rough skull-stripping */  
  morph_open_uint8(label, dims, strip_param[0], 3);
  get_largest_cluster(label, dims);
  morph_dilate_uint8(label, dims, strip_param[1], 0);
  morph_close_uint8(label, dims, 10, 0);
  
  /* update mask */
  for (i = 0; i < src_ptr->nvox; i++)
    mask[i] = 255*(label[i] > 0);

  n_classes = 6;
  double mean[n_classes];
  prob  = (unsigned char *)malloc(sizeof(unsigned char)*src_ptr->nvox*n_classes);

  /* use Kmeans with 6 classes */
  max_vol = Kmeans( src, label, mask, 25, 3, separations, dims, thresh, thresh_kmeans_int, iters_nu, 1, 50.0);

  /* use amap approach with PVE */
  Amap( src, label, prob, mean, 3, 10, 16, dims, 1, 0.5);
  Pve6(src, prob, label, mean, dims, PVELABEL);

  /* final skull-stripping */
  morph_open_uint8(label, dims, strip_param[2], 160);
  get_largest_cluster(label, dims);
  morph_dilate_uint8(label, dims, strip_param[3], 0);
  morph_close_uint8(label, dims, 10, 0);

  basename = nifti_makebasename(output_filename);

  slope = 1.0;

  /* write nu corrected image */
  if (write_nu) {
     (void) sprintf( buffer, "%s_unmasked%s",basename,extension); 

    if(!write_nifti( buffer, src, DT_FLOAT32, slope, dims, 
            separations, src_ptr))
      exit(EXIT_FAILURE);
  }
  
  /* write masked image */
  if (write_masked) {

    for (i = 0; i < src_ptr->nvox; i++)
      src[i] = (double)label[i]*src[i];

    if(!write_nifti( output_filename, src, DT_FLOAT32, slope, dims, 
            separations, src_ptr))
      exit(EXIT_FAILURE);
  }

  /* write brainmask */
  if (write_brainmask) {

    for (i = 0; i < src_ptr->nvox; i++)
      src[i] = (double)label[i];

    (void) sprintf( buffer, "%s_brainmask%s",basename,extension); 
    
    if(!write_nifti(buffer, src, DT_UINT8, slope, 
            dims, separations, src_ptr))
      exit(EXIT_FAILURE);
    
  }
    
  free(src);
  free(label);
  free(mask);
  free(prob);
  
  return(EXIT_SUCCESS);
}