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

/* Main program */

int main(int argc, char *argv[])
{
  char **infiles, *outfile;
  int i, j, nfiles, dims[3];
  double *avg, *input, separations[3];
  nifti_image *nii_ptr, *nii_ptr2;

  /* Get arguments */
  if (argc < 2) {
   (void) fprintf(stderr, 
   "\nUsage: %s [options] [<in1.nii> ...] <out.nii>\n",
        argv[0]);
   (void) fprintf(stderr, 
   "    %s -help\n\n", argv[0]);
   exit(EXIT_FAILURE);
  }
  
  nfiles = argc - 2;
  infiles = &argv[1];
  outfile = argv[argc-1];

  /* Make sure that we have something to process */
  if (nfiles == 0) {
   (void) fprintf(stderr, "Error: No input files specified\n");
   exit(EXIT_FAILURE);
  }
  
  /* read first image to get image parameters */
  nii_ptr = read_nifti_float(infiles[0], &input);
  fprintf(stdout,"%3d: %s\n",0, infiles[0]);

  separations[0] = nii_ptr->dx;
  separations[1] = nii_ptr->dy;
  separations[2] = nii_ptr->dz;
  dims[0] = nii_ptr->nx;
  dims[1] = nii_ptr->ny;
  dims[2] = nii_ptr->nz;

  /* prepare average image */
  avg  = (double *)malloc(sizeof(double)*nii_ptr->nvox);
  for (i=0; i<nii_ptr->nvox; i++) 
    avg[i] = input[i]/(double)nfiles;
  free(input);
  
  /* read remaining images and check for image parameters */
  for (i=1; i<nfiles; i++) {
    fprintf(stdout,"%3d: %s\n",i, infiles[i]);
    nii_ptr2 = read_nifti_float(infiles[i], &input);
    
    /* check for dimensions */
    if(!equal_image_dimensions(nii_ptr,nii_ptr2))
      return(0);
    
    /* calculate average */
    for (j=0; j<nii_ptr->nvox; j++) 
      avg[j] = avg[j] + (input[j]/(double)nfiles);
    free(input);
  }

  if (!write_nifti( outfile, avg, DT_FLOAT32, 1.0, dims, separations, nii_ptr)) 
    exit(EXIT_FAILURE);

  free(avg);

}