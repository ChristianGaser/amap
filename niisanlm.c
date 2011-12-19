/*
 * Christian Gaser
 * $Id$ 
 *
 */

#include <float.h>
#include <stdlib.h>

#include "nifti/nifti1_io.h"
#include "nifti/nifti1_local.h"

extern nifti_image *read_nifti_double( const char *input_filename, double *image[]);
extern write_nifti( const char *output_filename, double image[], int data_type, double slope, int dim[], double vox[], nifti_image *in_ptr);

/* Main program */

int main(int argc, char *argv[])
{
  char *infile, *outfile;
  int i, j, dims[3];
  double *input, separations[3];
  nifti_image *nii_ptr;
  
  if(argc < 3)
  {
    fprintf(stderr,"\n\
Usage: %s input.nii output.nii\n\n\
	Spatial adaptive non-local means denoising filter.\n\n", argv[0]);
    return( 1 );
  }
  
  infile  = argv[1];
  outfile = argv[2];
  
  /* read first image to get image parameters */
  nii_ptr = read_nifti_double(infile, &input);
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
  
  
  anlm(input, 3, 1, 1, dims);

  if (!write_nifti( outfile, input, DT_FLOAT32, 1.0, dims, separations, nii_ptr)) 
    exit(EXIT_FAILURE);

  free(input);
  
  return(EXIT_SUCCESS);

}