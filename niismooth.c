/*
 * Christian Gaser
 * $Id$ 
 *
 */

#include <ParseArgv.h>
#include <float.h>
#include <stdlib.h>

#include "nifti1/nifti1_io.h"
#include "nifti1/nifti1_local.h"

extern nifti_image *read_nifti_float( const char *input_filename, double *image[]);
extern equal_image_dimensions(nifti_image *nii_ptr, nifti_image *nii_ptr2);
extern write_nifti( const char *output_filename, double image[], int data_type, double slope, int dim[], double vox[], nifti_image *in_ptr);

double fwhm = 8;

static ArgvInfo argTable[] = {
  {"-fwhm", ARGV_FLOAT, (char *) 1, (char *) &fwhm, 
       "FWHM."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

/* Main program */

int main(int argc, char *argv[])
{
  char *infile, *outfile;
  int i, j, dims[3];
  double *input, separations[3];
  nifti_image *nii_ptr;
  double filt[3]={1,1,1};

  /* Get arguments */
  if  (ParseArgv(&argc, argv, argTable, 0) ||(argc < 2)) {
   (void) fprintf(stderr, 
   "\nUsage: %s -fwhm fwhm_in_mm <in.nii> <out.nii>\n",
        argv[0]);
   (void) fprintf(stderr, 
   "    %s -help\n\n", argv[0]);
   exit(EXIT_FAILURE);
  }
  
  infile  = argv[1];
  outfile = argv[2];
  
  /* read first image to get image parameters */
  nii_ptr = read_nifti_float(infile, &input);
  if(nii_ptr == NULL) {
    fprintf(stderr,"Error reading %s.\n", infile);
    return(EXIT_FAILURE);
  }

  separations[0] = nii_ptr->dx;
  separations[1] = nii_ptr->dy;
  separations[2] = nii_ptr->dz;
  dims[0] = nii_ptr->nx;
  dims[1] = nii_ptr->ny;
  dims[2] = nii_ptr->nz;
  
  morph_open_double(input, dims, 1, 0.098436);
//  morph_dilate_double(input, dims, 2, 0.5);
  morph_close_double(input, dims, 10, 0.5);

  if (!write_nifti( outfile, input, DT_FLOAT32, 1.0, dims, separations, nii_ptr)) 
    exit(EXIT_FAILURE);

  free(input);
  
  return(EXIT_SUCCESS);

}