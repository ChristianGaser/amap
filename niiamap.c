#include <ParseArgv.h>
#include <limits.h>
#include <float.h>
#include <time_stamp.h>

#include <bicpl.h>

#include "Amap.h"
#include "niiamap.h"

#include "nifti1/nifti1_io.h"
#include "nifti1/nifti1_local.h"

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
  {"-write_fuzzy", ARGV_CONSTANT, (char *) 1, (char *) &write_fuzzy,
       "Write fuzzy segmentations as separate images."},
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
        "usage: niiamap ...\n"
    };
    fprintf(stderr, "%s", msg);
    return (-1);
}


nifti_image
*read_nifti_float( const char *input_filename, double *image[])
{
  nifti_image *src_ptr;
  double tmp;
  int i;
  
  src_ptr = nifti_image_read(input_filename, 1);
  *image = (double *)malloc(sizeof(double)*src_ptr->nvox);
  
  for (i = 0; i < src_ptr->nvox; i++) {
      switch (src_ptr->datatype) {
      case DT_INT8:
        tmp = (double) ((signed char *)src_ptr->data)[i];
        break;
      case DT_UINT8:
        tmp = (double) ((unsigned char *)src_ptr->data)[i];
        break;
      case DT_INT16:
        tmp = (double) ((signed short *)src_ptr->data)[i];
        break;
      case DT_UINT16:
        tmp = (double) ((unsigned short *)src_ptr->data)[i];
        break;
      case DT_INT32:
        tmp = (double) ((signed int *)src_ptr->data)[i];
        break;
      case DT_UINT32:
        tmp = (double) ((unsigned int *)src_ptr->data)[i];
        break;
      case DT_FLOAT32:
        tmp = (double) ((float *)src_ptr->data)[i];
        break;
      case DT_FLOAT64:
        tmp = (double) ((double *)src_ptr->data)[i];
        break;
      default:
        fprintf(stderr,"Unknown datatype\n");
        return(NULL);
        break;
      }
      /* check whether scaling is needed */
      if (src_ptr->scl_slope == 0)
        (*image)[i] = tmp;
      else
        (*image)[i] = (src_ptr->scl_slope * tmp) + src_ptr->scl_inter;
    }  
    free(src_ptr->data);
    
    return(src_ptr);
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
  char      *input_filename, *output_filename, *extension;
  int       i, j, dims[3], thresh, thresh_kmeans_int;
  int		x, y, z, z_area, y_dims, count_zero;
  char		*axis_order[3] = { MIzspace, MIyspace, MIxspace };
  char		*arg_string, buffer[1024], *str_ptr;
  unsigned char *label, *prob, *mask, *marker, *init_mask, *priors;
  double	*src, *prevsrc, *tmp_vol, ratio_zeros;
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
    return(-1);
  }

  if (iters_nu <= 0)
    correct_nu = 0;

  if (Niters == 0)
    write_fuzzy = 0;

  if (correct_nu)
    fprintf(stdout,"Nu correction.\n");
  else {
    fprintf(stdout,".\n");
    iters_nu = -1;
  }

  /* do not write nu corrected image if correction is not selected */
  if (!correct_nu) write_nu = 0;

  /* read data and scale it to 0..255 */
  strcpy(buffer, input_filename);
  str_ptr = strrchr(buffer, '.');
  src_ptr = read_nifti_float(input_filename, &src);

  mask  = (unsigned char *)malloc(sizeof(unsigned char)*src_ptr->nvox);
  label = (unsigned char *)malloc(sizeof(unsigned char)*src_ptr->nvox);
  mask  = (unsigned char *)malloc(sizeof(unsigned char)*src_ptr->nvox);
  prob  = (unsigned char *)malloc(sizeof(unsigned char)*src_ptr->nvox*n_classes);

  /* read mask and check for same size */
  if (mask_filename != NULL) {
      
    /* read volume */
    mask_ptr = read_nifti_float(mask_filename, &tmp_vol);
    
    /* check size */    
    if ((mask_ptr->nx != src_ptr->nx) || (mask_ptr->ny != src_ptr->ny) || (mask_ptr->nz != src_ptr->nz) ||
        (mask_ptr->dx != src_ptr->dx) || (mask_ptr->dy != src_ptr->dy) || (mask_ptr->dz != src_ptr->dz)) {
      fprintf(stderr,"Mask and source image have different size\n");
      return(-1);
    }
    
    /* get min/max */  
    min_vol =  FLT_MAX; max_vol = -FLT_MAX;
    for (i = 0; i < mask_ptr->nvox; i++) {
      min_vol = MIN(tmp_vol[i], min_vol);
      max_vol = MAX(tmp_vol[i], max_vol);
    }
    /* scale mask image to a range 0..255 */
    for (i = 0; i < mask_ptr->nvox; i++) 
      mask[i] = (unsigned char) round(255*(tmp_vol[i] - min_vol)/(max_vol - min_vol));
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
    fprintf(stderr,"Warning: Only %g%s of the voxels outside the mask have zero values. This points to an image, which is not skull-stripped.\n", ratio_zeros,"%");
     
  separations[0] = src_ptr->dx;
  separations[1] = src_ptr->dy;
  separations[2] = src_ptr->dz;
  dims[0] = src_ptr->nx;
  dims[1] = src_ptr->ny;
  dims[2] = src_ptr->nz;
    
  /* initial nu-correction works best with 5 class Kmeans approach followed by a 3 class approach */
  max_vol = Kmeans( src, label, mask, 25, n_pure_classes, separations, dims, thresh, thresh_kmeans_int, iters_nu, KMEANS);
  max_vol = Kmeans( src, label, mask, 25, n_pure_classes, separations, dims, thresh, thresh_kmeans_int, iters_nu, NOPVE);

  /* final Kmeans estimation */
  max_vol = Kmeans( src, label, mask, 25, n_pure_classes, separations, dims, thresh, thresh_kmeans_int, iters_nu, pve);

  /* PVE */
  if (pve) {
    fprintf(stdout,"Calculate Partial Volume Estimate.\n");
    Pve6(src, prob, label, mean, dims, PVELABEL);
  }
  
  src_ptr->nifti_type = 1;
  
  if (!strcmp(extension,".img") == 1) {
    src_ptr->nifti_type = 2;
  }
  
  if (!strcmp(extension,".hdr") == 1) {
    src_ptr->nifti_type = 2;
    strcpy(extension,".img");
  }
  
  output_filename = nifti_makebasename(output_filename);

  /* write labeled volume */
  if (write_label) {
    src_ptr->datatype = DT_UINT8;
    src_ptr->nbyper = 1;
    
    /* different ranges for pve */
    if (pve) src_ptr->scl_slope = 3.0/255.0;
    else src_ptr->scl_slope = 1.0;
    
    src_ptr->scl_inter = 0;

    src_ptr->data = NULL;
    src_ptr->data = (unsigned char *)malloc(sizeof(unsigned char)*src_ptr->nvox);
    memcpy(src_ptr->data, label, sizeof(unsigned char)*src_ptr->nvox);

    (void) sprintf( buffer, "%s%s",output_filename,extension);

    src_ptr->iname = NULL;
    src_ptr->iname = malloc(strlen(buffer));
    strcpy(src_ptr->iname, buffer);

    if ((src_ptr->nifti_type == 0) || (src_ptr->nifti_type == 2)) {
      (void) sprintf( buffer, "%s%s",output_filename,".hdr");
    }
    src_ptr->fname = NULL;
    src_ptr->fname = malloc(strlen(buffer));
    strcpy(src_ptr->fname, buffer);

    nifti_image_write(src_ptr);

    free(src_ptr->data);
  }
  
  /* write fuzzy segmentations for each class */
  if (write_fuzzy) {
    src_ptr->datatype = DT_UINT8;
    src_ptr->nbyper = 1;
    src_ptr->scl_slope = 1.0/255.0;
    src_ptr->scl_inter = 0;

    src_ptr->data = NULL;
    src_ptr->data = (unsigned char *)malloc(sizeof(unsigned char)*src_ptr->nvox*n_classes);

    for (j = 0; j<n_pure_classes; j++) {
    
      (void) sprintf( buffer, "%s_seg%d%s",output_filename,j,extension);
      
      for (i = 0; i < src_ptr->nvox; i++)
        ((unsigned char *)src_ptr->data)[i] = prob[i+(j*src_ptr->nvox)];
      
      src_ptr->iname = NULL;
      src_ptr->iname = malloc(strlen(buffer));
      strcpy(src_ptr->iname, buffer);

      if ((src_ptr->nifti_type == 0) || (src_ptr->nifti_type == 2)) {
        (void) sprintf( buffer, "%s_seg%d%s",output_filename,j,".hdr");
      }
      src_ptr->fname = NULL;
      src_ptr->fname = malloc(strlen(buffer));
      strcpy(src_ptr->fname, buffer);

      nifti_image_write(src_ptr);
    }
    
    free(src_ptr->data);
  }

  /* write nu corrected image */
  if (write_nu) {
    src_ptr->datatype = DT_FLOAT32;
    src_ptr->nbyper = 4;
    src_ptr->scl_slope = 1;
    src_ptr->scl_inter = 0;

    src_ptr->data = NULL;
    src_ptr->data = (float *)malloc(sizeof(float)*src_ptr->nvox);
      for (i = 0; i < src_ptr->nvox; i++)
        ((float *)src_ptr->data)[i] = src[i];

    (void) sprintf( buffer, "%s_nu%s",output_filename,extension);

    src_ptr->iname = NULL;
    src_ptr->iname = malloc(strlen(buffer));
    strcpy(src_ptr->iname, buffer);

    if ((src_ptr->nifti_type == 0) || (src_ptr->nifti_type == 2)) {
      (void) sprintf( buffer, "%s_nu%s",output_filename,".hdr");
    }
    src_ptr->fname = NULL;
    src_ptr->fname = malloc(strlen(buffer));
    strcpy(src_ptr->fname, buffer);

    nifti_image_write(src_ptr);

    free(src_ptr->data);
  }

//  nifti_image_free(mask_ptr);
//  nifti_image_free(src_ptr);
  
  free(src);
  free(prob);
  free(label);
  free(mask);

}