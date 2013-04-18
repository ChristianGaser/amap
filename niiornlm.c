/*
 * Christian Gaser
 * $Id$ 
 *
 */

#include <float.h>
#include <stdlib.h>

#include "ParseArgv.h"
#include "niilib.h"

double h;

static ArgvInfo argTable[] = {
  {"-h", ARGV_FLOAT, (char *) 1, (char *) &h, 
       "Noise parameter h."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

/* Main program */

int main(int argc, char *argv[])
{
  char *infile, *outfile;
  int i, j, dims[3];
  double *input, *output, separations[3];
  nifti_image *nii_ptr;

  /* Get arguments */
  if  (ParseArgv(&argc, argv, argTable, 0) ||(argc < 2)) {
   (void) fprintf(stderr, 
   "\nUsage: %s -h smoothing_size in.nii out.nii\n",
        argv[0]);
   (void) fprintf(stderr, 
   "    %s -help\n\n", argv[0]);
   exit(EXIT_FAILURE);
  }
  
  infile  = argv[1];
  outfile = argv[2];
  
  /* read first image to get image parameters */
  nii_ptr = read_nifti_double(infile, &input, 0);
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
  
  output = (double *)malloc(sizeof(double)*dims[0]*dims[1]*dims[2]);
  
  ornlm(input, output, 3, 1, h, dims);

  if (!write_nifti_double( outfile, output, DT_FLOAT32, 1.0, dims, separations, nii_ptr)) 
    exit(EXIT_FAILURE);

  free(input);
  
  return(EXIT_SUCCESS);

}