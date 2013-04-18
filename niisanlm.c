/*
 * Christian Gaser
 * $Id$ 
 *
 */

#include <float.h>
#include <stdlib.h>
#if !defined(_WIN32)
  #include <libgen.h>
#endif

#include "niilib.h"

/* Main program */

int main(int argc, char *argv[])
{
  char *infile, outfile[1024];
  int i, j, dims[3];
  float *input;
  double separations[3];
  nifti_image *nii_ptr;
  
  if(argc < 2)
  {
    fprintf(stderr,"\n\
Usage: %s input.nii output.nii\n\n\
	Spatial adaptive non-local means denoising filter.\n\n", argv[0]);
    return( 1 );
  }
  
  infile  = argv[1];

  /* if not defined use original name as basename for output */
  if(argc == 3)
      (void) sprintf(outfile, "%s", argv[2]); 
  else {
    #if !defined(_WIN32)
      (void) sprintf(outfile, "%s/n%s", dirname(infile), basename(infile)); 
    #else
      fprintf(stderr,"\n\
Usage: %s input.nii output.nii\n\n\
	Spatial adaptive non-local means denoising filter.\n\n", argv[0]);
      return( 1 );
    #endif
  }
  
  /* read first image to get image parameters */
  nii_ptr = read_nifti_float(infile, &input, 0);
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
  
  
  sanlm(input, 3, 1, 1, dims);

  if (!write_nifti_float( outfile, input, DT_FLOAT32, 1.0, dims, separations, nii_ptr)) 
    exit(EXIT_FAILURE);

  free(input);
  
  return(EXIT_SUCCESS);

}