/*
 * Christian Gaser
 * $Id$ 
 *
 */

#include <ParseArgv.h>
#include <float.h>

#include "Amap.h"
#include "niiamap.h"

#include "nifti1/nifti1_io.h"
#include "nifti1/nifti1_local.h"

extern nifti_image *read_nifti_float( const char *input_filename, double *image[]);
	
static ArgvInfo argTable[] = {
  {"-mask", ARGV_STRING, (char *) 1, (char *) &mask_filename, 
       "Prior brainmask."},
  {"-iters", ARGV_INT, (char *) 1, (char *) &iters_amap,
       "Number of iterations to end."},
  {"-sub", ARGV_INT, (char *) 1, (char *) &subsample,
       "Subsampling for Amap approach."},
  {"-iters_nu", ARGV_INT, (char *) 1, (char *) &iters_nu,
       "Number of iterations for nu correction."},
  {"-iters_icm", ARGV_INT, (char *) 1, (char *) &iters_ICM,
       "Number of iterations for Iterative Conditional Mode (ICM)."},
  {"-no_nucorrect", ARGV_CONSTANT, (char *) 0, (char *) &correct_nu,
       "Do not use nu correction."},
  {"-mrf", ARGV_FLOAT, (char *) 1, (char *) &weight_MRF,
       "Weight of MRF prior (0..1)."},
  {"-thresh", ARGV_FLOAT, (char *) 1, (char *) &thresh_brainmask,
       "Threshold for prior brainmask (0..1)."},
  {"-thresh_kmeans", ARGV_FLOAT, (char *) 1, (char *) &thresh_kmeans,
       "Threshold for Kmeans algorithm (0..1)."},
  {"-bias", ARGV_FLOAT, (char *) 1, (char *) &bias_fwhm,
       "Bias field spline smoothing (FWHM) in mm."},
  {"-pve", ARGV_INT, (char *) 1, (char *) &pve,
       "Use Partial Volume Estimation with 5 classes (5), 6 classes (6) or do not use PVE (0)."},
  {"-write_seg", ARGV_INT, (char *) 3, (char *) &write_seg,
       "Write fuzzy segmentations as separate images. Three numbers should be given, while a '1' indicates that this tissue class should be saved. Order is CSF/GM/WM."},
  {"-write_nu", ARGV_CONSTANT, (char *) 1, (char *) &write_nu,
       "Write nu corrected image."},
  {"-write_label", ARGV_CONSTANT, (char *) 1, (char *) &write_label,
       "Write label image (default)."},
  {"-nowrite_label", ARGV_CONSTANT, (char *) 0, (char *) &write_label,
       "Do not write label image."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};


static int usage(void)
{
    static const char msg[] = {
        "niiamap: Segmentation with adaptive MAP\n"
        "usage: niiamap [options] in.nii [out.nii]\n"
    };
    fprintf(stderr, "%s", msg);
    exit(EXIT_FAILURE);
}

int
main( int argc, char **argv )
{
  /* NIFTI stuff */
  nifti_image *src_ptr, *mask_ptr;
  int       nii_dimids[MAX_NII_DIMS];
  int       nii_dir[MAX_NII_DIMS];
  int       nii_map[MAX_NII_DIMS];
  unsigned long nii_lens[MAX_NII_DIMS];
  int       nii_ndims, n_classes;
  int       nifti_file_type;
  char      *input_filename, *output_filename, *basename, *extension;
  int       i, j, dims[3], thresh, thresh_kmeans_int;
  int		    x, y, z, z_area, y_dims, count_zero;
  char		  *arg_string, buffer[1024];
  unsigned char *label, *prob, *mask, *marker, *init_mask, *priors;
  double	  *src, *buffer_vol, ratio_zeros, slope;
  double    offset, val, max_vol, min_vol, voxelsize[3];

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

  if (iters_nu <= 0)
    correct_nu = 0;

  if (iters_amap == 0) {
    for (i = 0; i < 3; i++) write_seg[i] = 0;
    fprintf(stdout,"To write segmentation you need at least one iteration.\n");
  }

  /* do not write nu corrected image if correction is not selected */
  if (!correct_nu) {
    iters_nu = -1;
    write_nu = 0;
  }
  
  if((pve != 0) && (pve != 5) && (pve != 6)) {
    fprintf(stderr,"Value for pve can be either 0 (no PVE), 5 or 6 (5 or 6 classes).\n");
    exit(EXIT_FAILURE);
  }

  switch(pve) {
  case 0:
    n_classes = 3;
    break;
  case 5:
    n_classes = 5;
    break;
  case 6:
    n_classes = 6;
    break;
  }
  
  /* read data and scale it to 0..255 */
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
  prob  = (unsigned char *)malloc(sizeof(unsigned char)*src_ptr->nvox*n_classes);
  
  if((mask == NULL) || (label == NULL) || (prob == NULL)) {
    fprintf(stderr,"Memory allocation error\n");
    exit(EXIT_FAILURE);
  }

  /* read mask and check for same size */
  if (mask_filename != NULL) {
      
    /* read volume */
    mask_ptr = read_nifti_float(mask_filename, &buffer_vol);
    if(mask_ptr == NULL) {
      fprintf(stderr,"Error reading %s.\n", mask_filename);
      return(EXIT_FAILURE);
    }
    
    /* check size */ 
    if (!equal_image_dimensions(src_ptr,mask_ptr)) {   
      fprintf(stderr,"Mask and source image have different size\n");
      exit(EXIT_FAILURE);
    }
    
    /* get min/max */  
    min_vol =  FLT_MAX; max_vol = -FLT_MAX;
    for (i = 0; i < mask_ptr->nvox; i++) {
      min_vol = MIN(buffer_vol[i], min_vol);
      max_vol = MAX(buffer_vol[i], max_vol);
    }
    /* scale mask image to a range 0..255 */
    for (i = 0; i < mask_ptr->nvox; i++) 
      mask[i] = (unsigned char) round(255*(buffer_vol[i] - min_vol)/(max_vol - min_vol));
  }

  /* get sure that threshold for brainmask is zero if no mask is defined */
  if (mask_filename == NULL) thresh_brainmask = 0.0;

  double mean[n_classes], mu[n_pure_classes], var[n_pure_classes];
  for (i = 0; i < n_pure_classes; i++)
    mu[i] = 0;

  thresh = (int)round(255*thresh_brainmask);
  thresh_kmeans_int = (int)round(255*thresh_kmeans);

  /* get min/max */
  min_vol =  FLT_MAX; max_vol = -FLT_MAX;
  for (i = 0; i < src_ptr->nvox; i++) {
    min_vol = MIN(src[i], min_vol);
    max_vol = MAX(src[i], max_vol);
  }

  /* if no mask file is given use minimum value or zeros in the image to get mask value */
  if (mask_filename == NULL) {
    for (i = 0; i < src_ptr->nvox; i++) {
      if ((src[i] == min_vol) || (src[i] == 0.0)) mask[i] = 0;
      else mask[i] = 255;
    }  
  }

  /* correct images with values < 0 */
  if (min_vol < 0) {
    for (i = 0; i < src_ptr->nvox; i++)
      src[i] = src[i] - min_vol;
  }

  /* add offset to ensure that CSF values are much larger than background noise */
  offset = 0.2*max_vol;  
  for (i = 0; i < src_ptr->nvox; i++)
    if (mask[i] > 0) src[i] += offset;

  count_zero = 0;
  for (i = 0; i < src_ptr->nvox; i++)
    if ((mask[i] == 0)) count_zero++;
  
  ratio_zeros = 100.0*(double)count_zero/(double)(src_ptr->nvox);
  if(ratio_zeros < 5)
    fprintf(stdout,"Warning: Only %g%s of the voxels outside the mask have zero values. This points to an image, which is not skull-stripped.\n", ratio_zeros,"%");
     
  voxelsize[0] = src_ptr->dx;
  voxelsize[1] = src_ptr->dy;
  voxelsize[2] = src_ptr->dz;
  dims[0] = src_ptr->nx;
  dims[1] = src_ptr->ny;
  dims[2] = src_ptr->nz;

  /* initial nu-correction works best with 6 class Kmeans approach followed by a 3 class approach */
  if (correct_nu)
    max_vol = Kmeans( src, label, mask, 25, n_pure_classes, voxelsize, dims, thresh, thresh_kmeans_int, iters_nu, KMEANS, bias_fwhm);
  
  /* final Kmeans estimation */
  max_vol = Kmeans( src, label, mask, 25, n_pure_classes, voxelsize, dims, thresh, thresh_kmeans_int, iters_nu, NOPVE, bias_fwhm);

  Amap( src, label, prob, mean, n_pure_classes, iters_amap, subsample, dims, pve, weight_MRF, voxelsize, iters_ICM, offset);

  /* PVE */
  if (pve) {
    fprintf(stdout,"Calculate Partial Volume Estimate.\n");
    if(pve==6)
      Pve6(src, prob, label, mean, dims);
    else
      Pve5(src, prob, label, mean, dims);
  }
  
  /* write nu corrected image */
  if (write_nu) {
     (void) sprintf( buffer, "%s_nu%s",basename,extension); 

    if(!write_nifti( buffer, src, DT_FLOAT32, 1.0, dims, 
            voxelsize, src_ptr))
      exit(EXIT_FAILURE);
  }
  
  /* write labeled volume */
  if (write_label) {

    /* different ranges for pve */
    if (pve) slope = 3.0/255.0;
    else slope = 1.0;

    for (i = 0; i < src_ptr->nvox; i++)
      src[i] = (double)label[i];

    (void) sprintf( buffer, "%s_seg%s",basename,extension); 

    if(!write_nifti(buffer, src, DT_UINT8, slope, 
            dims, voxelsize, src_ptr))
      exit(EXIT_FAILURE);
    
  }
  
  /* write fuzzy segmentations for each class */
  if (write_seg[0] || write_seg[1] || write_seg[2]) {
    
    slope = 1.0/255.0;

    for (j = 0; j < n_pure_classes; j++) {
      if (write_seg[j]) {
        (void) sprintf( buffer, "%s_prob%d%s",basename,j,extension); 
        
        for (i = 0; i < src_ptr->nvox; i++)
          src[i] = prob[i+(j*src_ptr->nvox)];
        
        if(!write_nifti( buffer, src, DT_UINT8, slope, 
                dims, voxelsize, src_ptr))
          exit(EXIT_FAILURE);
      
      }
    }    
  }
  
  free(src);
  free(prob);
  free(label);
  free(mask);
  
  return(EXIT_SUCCESS);
}