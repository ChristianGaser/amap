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
  {"-iters", ARGV_INT, (char *) 1, (char *) &Niters,
       "Number of iterations to end."},
  {"-sub", ARGV_INT, (char *) 1, (char *) &subsample,
       "Subsampling for Amap approach."},
  {"-iters_nu", ARGV_INT, (char *) 1, (char *) &iters_nu,
       "Number of iterations for nu correction."},
  {"-no_nucorrect", ARGV_CONSTANT, (char *) 0, (char *) &correct_nu,
       "Do not use nu correction."},
  {"-thresh", ARGV_FLOAT, (char *) 1, (char *) &thresh_brainmask,
       "Threshold for prior brainmask (0..1)."},
  {"-thresh_kmeans", ARGV_FLOAT, (char *) 1, (char *) &thresh_kmeans,
       "Threshold for Kmeans algorithm (0..1)."},
  {"-pve", ARGV_INT, (char *) 1, (char *) &pve,
       "Use Partial Volume Estimation with marginalized likelihood estimation (1) or Kmeans initialization (2) or do not use PVE (0)."},
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
        "usage: niiamap [options] in.nii out.nii\n"
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
  int       nii_ndims;
  int       nifti_file_type;
  char      *input_filename, *output_filename, *basename, *extension;
  int       i, j, dims[3], thresh, thresh_kmeans_int;
  int		x, y, z, z_area, y_dims, count_zero;
  char		*axis_order[3] = { MIzspace, MIyspace, MIxspace };
  char		*arg_string, buffer[1024], *str_ptr;
  unsigned char *label, *prob, *mask, *marker, *init_mask, *priors;
  double	*src, *vol_double, *buffer_vol, ratio_zeros, slope;
  double    val, max_vol, min_vol, separations[3];

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

  if (iters_nu <= 0)
    correct_nu = 0;

  if (Niters == 0) {
    for (i = 0; i < 3; i++) write_seg[i] = 0;
  }

  if (correct_nu)
    fprintf(stdout,"Nu correction.\n");
  else {
    iters_nu = -1;
  }

  /* do not write nu corrected image if correction is not selected */
  if (!correct_nu) write_nu = 0;

  /* read data and scale it to 0..255 */
  strcpy(buffer, input_filename);
  str_ptr = strrchr(buffer, '.');
  src_ptr = read_nifti_float(input_filename, &src);
  
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

  double mean[n_classes], mu[n_pure_classes];
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

  /* correct images with values < 0 */
  if (min_vol < 0) {
    for (i = 0; i < src_ptr->nvox; i++)
      src[i] = src[i] - min_vol;
    min_vol = 0.0;
  }

  /* if no mask file is given use minimum value in the image to get mask value */
  if (mask_filename == NULL) {
    for (i = 0; i < src_ptr->nvox; i++) {
      if (src[i] == min_vol) mask[i] = 0;
      else mask[i] = 255;
    }  
  }

  count_zero = 0;
  for (i = 0; i < src_ptr->nvox; i++)
    if ((mask[i] == 0)) count_zero++;
  
  ratio_zeros = 100.0*(double)count_zero/(double)(src_ptr->nvox);
  if(ratio_zeros < 5)
    fprintf(stdout,"Warning: Only %g%s of the voxels outside the mask have zero values. This points to an image, which is not skull-stripped.\n", ratio_zeros,"%");
     
  separations[0] = src_ptr->dx;
  separations[1] = src_ptr->dy;
  separations[2] = src_ptr->dz;
  dims[0] = src_ptr->nx;
  dims[1] = src_ptr->ny;
  dims[2] = src_ptr->nz;
    
  /* initial nu-correction works best with 5 class Kmeans approach followed by a 3 class approach */
  max_vol = Kmeans( src, label, mask, 25, n_pure_classes, separations, dims, thresh, thresh_kmeans_int, iters_nu, KMEANS);
  max_vol = Kmeans( src, label, mask, 25, n_pure_classes, separations, dims, thresh, thresh_kmeans_int, iters_nu, NOPVE);

  /* final Kmeans estimation if nu-correction was selected */
  if (correct_nu)
    max_vol = Kmeans( src, label, mask, 25, n_pure_classes, separations, dims, thresh, thresh_kmeans_int, iters_nu, pve);

  if (Niters > 0)
    Amap( src, label, prob, mean, n_pure_classes, Niters, subsample, dims, pve);

  /* PVE */
  if (pve) {
    fprintf(stdout,"Calculate Partial Volume Estimate.\n");
    Pve6(src, prob, label, mean, dims, PVELABEL);
  }
  
  basename = nifti_makebasename(output_filename);

  /* write nu corrected image */
  if (write_nu) {
     (void) sprintf( buffer, "%s_nu%s",basename,extension); 

    write_nifti( buffer, src, DT_FLOAT32, 1.0, dims, separations, src_ptr);
  }
  
  /* write labeled volume */
  if (write_label) {

    /* different ranges for pve */
    if (pve) slope = 3.0/255.0;
    else slope = 1.0;

    for (i = 0; i < src_ptr->nvox; i++)
      src[i] = (double)label[i];

    write_nifti(output_filename, src, DT_UINT8, slope, dims, separations, src_ptr);
    
  }
  
  /* write fuzzy segmentations for each class */
  if (write_seg[0] || write_seg[1] || write_seg[2]) {
    
    slope = 1.0/255.0;

    for (j = 0; j<n_pure_classes; j++) {
      if (write_seg[j]) {
        (void) sprintf( buffer, "%s_seg%d%s",basename,j,extension); 
        
        for (i = 0; i < src_ptr->nvox; i++)
          src[i] = prob[i+(j*src_ptr->nvox)];
        
        write_nifti( buffer, src, DT_UINT8, slope, dims, separations, src_ptr);
      
      }
    }    
  }
  
  free(src);
  free(prob);
  free(label);
  free(mask);
}