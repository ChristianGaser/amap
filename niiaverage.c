/*
 * Christian Gaser
 * $Id$ 
 *
 */

#include <ParseArgv.h>
#include <float.h>

#include "nifti1/nifti1_io.h"
#include "nifti1/nifti1_local.h"

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
   (void) fprintf(stderr, "No input files specified\n");
   exit(EXIT_FAILURE);
  }

  nii_ptr = read_nifti_float(infiles[0], &input);

  separations[0] = nii_ptr->dx;
  separations[1] = nii_ptr->dy;
  separations[2] = nii_ptr->dz;
  dims[0] = nii_ptr->nx;
  dims[1] = nii_ptr->ny;
  dims[2] = nii_ptr->nz;

  avg  = (double *)malloc(sizeof(double)*nii_ptr->nvox);
  for (i=0; i<nii_ptr->nvox; i++) 
    avg[i] = input[i]/(double)nfiles;
  
  for (i=1; i<nfiles; i++) {
fprintf(stderr,"%s\n",infiles[i]);
    nii_ptr2 = read_nifti_float(infiles[i], &input);
    /* check for dimensions */
    for (j=0; j<nii_ptr->nvox; j++) 
      avg[j] = avg[j] + (input[j]/(double)nfiles);
  }

  write_nifti( outfile, avg, DT_FLOAT32, 1.0, dims, separations, nii_ptr);

}