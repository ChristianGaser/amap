/*
 * Christian Gaser
 * $Id: niinucorrect.c 191 2012-07-31 15:41:36Z gaser $ 
 *
 */

#include <float.h>
#include "ParseArgv.h"
#include "niilib.h"
#include "Amap.h"
#include "niinucorrect.h"
	
static ArgvInfo argTable[] = {
  {"-iters_nu", ARGV_INT, (char *) 1, (char *) &iters_nu,
       "Number of iterations for nu correction."},
  {"-thresh_kmeans", ARGV_FLOAT, (char *) 1, (char *) &thresh_kmeans,
       "Threshold for Kmeans algorithm (0..1)."},
  {"-bias", ARGV_FLOAT, (char *) 1, (char *) &bias_fwhm,
       "Bias field spline smoothing (FWHM) in mm."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};


static int usage(void)
{
    static const char msg[] = {
        "niinucorrect: Nu-correction\n"
        "usage: niiamap [options] in.nii [out.nii]\n"
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
  int       i, dims[3], thresh, thresh_kmeans_int;
  int		    count_zero;
  char		  buffer[1024];
  unsigned char *label, *mask;
  float	  *src;
  double    ratio_zeros, avg, avg8, max_src, min_src, voxelsize[3];

  /* Get arguments */
  if (ParseArgv(&argc, argv, argTable, 0) || (argc < 2)) {
    (void) fprintf(stderr, 
    "\nUsage: %s [options] in.nii [out.nii]\n",
                     argv[0]);
    (void) fprintf(stderr, 
      "       %s -help\n\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  
  input_filename  = argv[1];

  /* if not defined use original name as basename for output */
  if(argc == 3)
    output_filename = argv[2];
  else
    output_filename = argv[1];
  
  /* get basename */
  basename = nifti_makebasename(output_filename);

  /* deal with extension */
  extension = nifti_find_file_extension(output_filename);
  
  /* if no valid extension was found use .nii */
  if (extension == NULL) {
    fprintf(stdout,"Use .nii as extension for %s.\n",output_filename);
    (void) sprintf( extension, ".nii");;
  }
  
  /* read data and scale it to 0..255 */
  src_ptr = read_nifti_float(input_filename, &src, 0);
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

  /* get sure that threshold for brainmask is zero if no mask is defined */
  thresh_brainmask = 0.5;

  thresh = (int)round(255*thresh_brainmask);
  thresh_kmeans_int = (int)round(255*thresh_kmeans);

  /* get min/max/avg */
  min_src =  FLT_MAX; max_src = -FLT_MAX;
  avg = 0.0;
  for (i = 0; i < src_ptr->nvox; i++) {
    min_src = MIN((double)src[i], min_src);
    max_src = MAX((double)src[i], max_src);
    avg += (float)src[i];
  }
  avg /= (double)src_ptr->nvox;
  
  avg8 = 0.0;
  for (i = 0; i < src_ptr->nvox; i++)
    if (src[i] > (float)avg/8.0) avg8 += (double)src[i];
  avg8 /= (double)src_ptr->nvox;

  /* correct images with values < 0 */
  if (min_src < 0) {
    avg8 -= min_src;
    for (i = 0; i < src_ptr->nvox; i++)
      src[i] -= (float)min_src;
  }

  /* if no mask file is given use minimum value or zeros in the image to get mask value */
  for (i = 0; i < src_ptr->nvox; i++) {
    if ((src[i] < (float)avg8) || (src[i] > 0.95*(float)max_src)) mask[i] = 0;
    else mask[i] = 255;
  }

  count_zero = 0;
  for (i = 0; i < src_ptr->nvox; i++)
    if (mask[i] == 0) count_zero++;
  
  ratio_zeros = 100.0*(double)count_zero/(double)(src_ptr->nvox);
  if(ratio_zeros < 5)
    fprintf(stdout,"Warning: Only %g%s of the voxels outside the mask have zero values. This points to an image, which is not skull-stripped.\n", ratio_zeros,"%");
     
  voxelsize[0] = src_ptr->dx;
  voxelsize[1] = src_ptr->dy;
  voxelsize[2] = src_ptr->dz;
  dims[0] = src_ptr->nx;
  dims[1] = src_ptr->ny;
  dims[2] = src_ptr->nz;

  /* initial nu-correction works best with 6 class Kmeans approach */
  max_src = Kmeans( src, label, mask, 25, n_pure_classes, voxelsize, dims, thresh, thresh_kmeans_int, iters_nu, KMEANS, bias_fwhm);
    
  /* write nu corrected image */
  (void) sprintf( buffer, "%s_nu%s",basename,extension); 

  if(!write_nifti_float( buffer, src, src_ptr->datatype, 1.0, dims, 
          voxelsize, src_ptr))
    exit(EXIT_FAILURE);
      
  free(src);
  free(label);
  free(mask);
  
  return(EXIT_SUCCESS);
}