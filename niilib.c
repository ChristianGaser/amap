/*
 * Christian Gaser
 * $Id$ 
 *
 */

#include <time_stamp.h>
#include <stdlib.h>

#include "nifti1/nifti1_io.h"
#include "nifti1/nifti1_local.h"

int
equal_image_dimensions(nifti_image *nii_ptr, nifti_image *nii_ptr2) {

  if((nii_ptr->nx != nii_ptr2->nx) ||
     (nii_ptr->ny != nii_ptr2->ny) ||
     (nii_ptr->nz != nii_ptr2->nz) ||
     (nii_ptr->dx != nii_ptr2->dx) ||
     (nii_ptr->dy != nii_ptr2->dy) ||
     (nii_ptr->dz != nii_ptr2->dz)) {
    fprintf(stderr,"Error: Image %s and image %s differ.\n",nii_ptr->fname, nii_ptr2->fname);
    return(0);    
  }
  return(1);
}

void
init_nifti_header(nifti_image *nii_ptr)
{
  int i, j;

  nii_ptr->ndim = 0;

  nii_ptr->nx = nii_ptr->ny = nii_ptr->nz = nii_ptr->nt = nii_ptr->nu = 
    nii_ptr->nv = nii_ptr->nw = 0;

  for (i = 0; i < MAX_NII_DIMS; i++) {
    nii_ptr->dim[i] = 1;
  }

  nii_ptr->nvox = 0;
  nii_ptr->nbyper = 0;
  nii_ptr->datatype = DT_UNKNOWN;

  nii_ptr->dx = nii_ptr->dy = nii_ptr->dz = nii_ptr->dt = nii_ptr->du = 
    nii_ptr->dv = nii_ptr->dw = 0.0;
  
  for (i = 0; i < MAX_NII_DIMS; i++) {
    nii_ptr->pixdim[i] = 0.0;
  }

  nii_ptr->num_ext = 0;
  nii_ptr->scl_slope = 0.0;
  nii_ptr->scl_inter = 0.0;
  nii_ptr->cal_min = 0.0;
  nii_ptr->cal_max = 0.0;

  nii_ptr->qform_code = NIFTI_XFORM_UNKNOWN;
  nii_ptr->sform_code = NIFTI_XFORM_UNKNOWN;

  nii_ptr->freq_dim = 0;
  nii_ptr->phase_dim = 0;
  nii_ptr->slice_dim = 0;

  nii_ptr->slice_code = 0;
  nii_ptr->slice_start = 0;
  nii_ptr->slice_end = 0;
  nii_ptr->slice_duration = 0.0;

  nii_ptr->quatern_b = 0.0;
  nii_ptr->quatern_c = 0.0;
  nii_ptr->quatern_d = 0.0;
  nii_ptr->qoffset_x = 0.0;
  nii_ptr->qoffset_y = 0.0;
  nii_ptr->qoffset_z = 0.0;
  nii_ptr->qfac = 0.0;

  nii_ptr->toffset = 0.0;

  nii_ptr->xyz_units = NIFTI_UNITS_MM; /* Default spatial units */
  nii_ptr->time_units = NIFTI_UNITS_SEC; /* Default time units */

  nii_ptr->nifti_type = FT_ANALYZE;
  nii_ptr->intent_code = 0;
  nii_ptr->intent_p1 = 0.0;
  nii_ptr->intent_p2 = 0.0;
  nii_ptr->intent_p3 = 0.0;
  memset(nii_ptr->intent_name, 0, sizeof (nii_ptr->intent_name));

  memset(nii_ptr->descrip, 0, sizeof (nii_ptr->descrip));

  memset(nii_ptr->aux_file, 0, sizeof (nii_ptr->aux_file));
  
  nii_ptr->fname = NULL;
  nii_ptr->iname = NULL;
  nii_ptr->iname_offset = 0;
  nii_ptr->swapsize = 0;
  nii_ptr->byteorder = 0;
  nii_ptr->data = NULL;

  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      nii_ptr->qto_xyz.m[i][j] = 0.0;
      nii_ptr->qto_ijk.m[i][j] = 0.0;
      nii_ptr->sto_xyz.m[i][j] = 0.0;
      nii_ptr->sto_ijk.m[i][j] = 0.0;
    }
  }
}  

int
write_nifti( const char *output_filename, double image[], int data_type, double slope, int dim[], double vox[], nifti_image *in_ptr)
{
  nifti_image *nii_ptr, nii_rec;
  char *extension, buffer[1024];
  int i;
  
  if((data_type != DT_UINT8) && 
     (data_type != DT_FLOAT32) && 
     (data_type != DT_INT8) && 
     (data_type != DT_INT16) && 
     (data_type != DT_INT32) && 
     (data_type != DT_FLOAT32)) {
    fprintf(stderr,"Datatype %d not supported to write data.\n",data_type);
    return(0);
  }
  
  if(in_ptr == NULL) {
    nii_ptr = &nii_rec;
    init_nifti_header(nii_ptr);
  } else nii_ptr = in_ptr;
  
  extension = nifti_find_file_extension(output_filename);
  
  /* if no valid extension was found use .nii */
  if (extension == NULL) {
    fprintf(stderr,"No valid extension found for output filename %s.\n",output_filename);
    return(0);
  }

  nii_ptr->nifti_type = 1;
  
  if (!strcmp(extension,".img") == 1) {
    nii_ptr->nifti_type = 2;
  }
  
  if (!strcmp(extension,".hdr") == 1) {
    nii_ptr->nifti_type = 2;
    strcpy(extension,".img");
  }

  output_filename = nifti_makebasename(output_filename);

  nii_ptr->nx = dim[0];
  nii_ptr->ny = dim[1];
  nii_ptr->nz = dim[2];
  nii_ptr->pixdim[0] = vox[0];
  nii_ptr->pixdim[1] = vox[1];
  nii_ptr->pixdim[2] = vox[2];
  nii_ptr->dx = vox[0];
  nii_ptr->dy = vox[1];
  nii_ptr->dz = vox[2];

  nii_ptr->data = NULL;
  nii_ptr->scl_slope = slope;
  nii_ptr->scl_inter = 0;
  nii_ptr->nvox = dim[0]*dim[1]*dim[2];
  nii_ptr->nifti_type = 1;
  nii_ptr->ndim = 3;

  nii_ptr->datatype = data_type;
  switch (data_type) {
  case DT_UINT8:
    nii_ptr->nbyper = 1;
    nii_ptr->data = (unsigned char *)malloc(sizeof(unsigned char)*nii_ptr->nvox);
    /* check for memory */
    if(nii_ptr->data == NULL) {
      fprintf(stderr,"Memory allocation error\n");
      return(0);
    }
    for (i = 0; i < nii_ptr->nvox; i++)
      ((unsigned char *)nii_ptr->data)[i] = image[i];
    break;
  case DT_INT8:
    nii_ptr->nbyper = 1;
    nii_ptr->data = (signed char *)malloc(sizeof(signed char)*nii_ptr->nvox);
    /* check for memory */
    if(nii_ptr->data == NULL) {
      fprintf(stderr,"Memory allocation error\n");
      return(0);
    }
    for (i = 0; i < nii_ptr->nvox; i++)
      ((signed char *)nii_ptr->data)[i] = image[i];
    break;
  case DT_INT16:
    nii_ptr->nbyper = 2;
    nii_ptr->data = (signed short *)malloc(sizeof(signed short)*nii_ptr->nvox);
    /* check for memory */
    if(nii_ptr->data == NULL) {
      fprintf(stderr,"Memory allocation error\n");
      return(0);
    }
    for (i = 0; i < nii_ptr->nvox; i++)
      ((signed short *)nii_ptr->data)[i] = image[i];
    break;
  case DT_INT32:
    nii_ptr->nbyper = 4;
    nii_ptr->data = (signed int *)malloc(sizeof(signed int)*nii_ptr->nvox);
    /* check for memory */
    if(nii_ptr->data == NULL) {
      fprintf(stderr,"Memory allocation error\n");
      return(0);
    }
    for (i = 0; i < nii_ptr->nvox; i++)
      ((signed int *)nii_ptr->data)[i] = image[i];
    break;
  case DT_INT64:
    nii_ptr->nbyper = 8;
    nii_ptr->data = (long long *)malloc(sizeof(long long)*nii_ptr->nvox);
    /* check for memory */
    if(nii_ptr->data == NULL) {
      fprintf(stderr,"Memory allocation error\n");
      return(0);
    }
    for (i = 0; i < nii_ptr->nvox; i++)
      ((long long *)nii_ptr->data)[i] = image[i];
    break;
  case DT_FLOAT32:
    nii_ptr->nbyper = 4;
    nii_ptr->data = (float *)malloc(sizeof(float)*nii_ptr->nvox);
    /* check for memory */
    if(nii_ptr->data == NULL) {
      fprintf(stderr,"Memory allocation error\n");
      return(0);
    }
    for (i = 0; i < nii_ptr->nvox; i++)
      ((float *)nii_ptr->data)[i] = image[i];
    break;
  }
    
  (void) sprintf( buffer, "%s%s",output_filename,extension);

  nii_ptr->iname = NULL;
  nii_ptr->iname = malloc(strlen(buffer));
  strcpy(nii_ptr->iname, buffer);

  if ((nii_ptr->nifti_type == 0) || (nii_ptr->nifti_type == 2)) {
    (void) sprintf( buffer, "%s%s",output_filename,".hdr");
  }
  nii_ptr->fname = NULL;
  nii_ptr->fname = malloc(strlen(buffer));
  strcpy(nii_ptr->fname, buffer);

  nifti_image_write(nii_ptr);

  free(nii_ptr->data);
  
  return(1);
}

nifti_image
*read_nifti_float( const char *input_filename, double *image[])
{
  nifti_image *nii_ptr;
  double tmp;
  int i;
  
  nii_ptr = nifti_image_read(input_filename, 1);
  if(nii_ptr == NULL) {
    fprintf(stderr,"read_nifti_float: Error reading %s.\n", input_filename);
    return(NULL);
  }
  
  /* read as double format */
  *image = (double *)malloc(sizeof(double)*nii_ptr->nvox);
  
  /* check for memory */
  if(image == NULL) {
    fprintf(stderr,"read_nifti_float: Memory allocation error\n");
    return(NULL);
  }
  
  for (i = 0; i < nii_ptr->nvox; i++) {
      switch (nii_ptr->datatype) {
      case DT_INT8:
        tmp = (double) ((signed char *)nii_ptr->data)[i];
        break;
      case DT_UINT8:
        tmp = (double) ((unsigned char *)nii_ptr->data)[i];
        break;
      case DT_INT16:
        tmp = (double) ((signed short *)nii_ptr->data)[i];
        break;
      case DT_UINT16:
        tmp = (double) ((unsigned short *)nii_ptr->data)[i];
        break;
      case DT_INT32:
        tmp = (double) ((signed int *)nii_ptr->data)[i];
        break;
      case DT_UINT32:
        tmp = (double) ((unsigned int *)nii_ptr->data)[i];
        break;
      case DT_FLOAT32:
        tmp = (double) ((float *)nii_ptr->data)[i];
        break;
      case DT_FLOAT64:
        tmp = (double) ((double *)nii_ptr->data)[i];
        break;
      default:
        fprintf(stderr,"read_nifti_float: Unknown datatype\n");
        return(NULL);
        break;
      }
      /* check whether scaling is needed */
      if (nii_ptr->scl_slope == 0)
        (*image)[i] = tmp;
      else
        (*image)[i] = (nii_ptr->scl_slope * tmp) + nii_ptr->scl_inter;
    }  
    free(nii_ptr->data);
    
    return(nii_ptr);
}
