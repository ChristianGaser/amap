/*
 * Christian Gaser
 * $Id$ 
 *
 */

#include "ParseArgv.h"
#include "niilib.h"

#include <float.h>
#include <stdlib.h>

extern int ParseArgv(int *argcPtr, char **argv, ArgvInfo *argTable, int flags);

char *std_filename = NULL;

static ArgvInfo argTable[] = {
  {"-std", ARGV_STRING, (char *) 1, (char *) &std_filename, 
       "Write standard deviation."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

/* Main program */

int main(int argc, char *argv[])
{
  char **infiles, *outfile;
  int i, j, nfiles, dims[3];
  double *avg, *sum_squares, *input, separations[3];
  nifti_image *nii_ptr, *nii_ptr2;

  /* Get arguments */
  if  (ParseArgv(&argc, argv, argTable, 0) ||(argc < 2)) {
   (void) fprintf(stderr, 
   "\nUsage: %s [-std std-out.nii] [<in1.nii> ...] <out.nii>\n",
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
  nii_ptr = read_nifti_double(infiles[0], &input);
  if(nii_ptr == NULL) {
    fprintf(stderr,"Error reading %s.\n", infiles[0]);
    return(EXIT_FAILURE);
  }
  fprintf(stdout,"%3d: %s\n",1, infiles[0]);

  separations[0] = nii_ptr->dx;
  separations[1] = nii_ptr->dy;
  separations[2] = nii_ptr->dz;
  dims[0] = nii_ptr->nx;
  dims[1] = nii_ptr->ny;
  dims[2] = nii_ptr->nz;

  /* prepare average image */
  avg  = (double *)malloc(sizeof(double)*nii_ptr->nvox);
  for (j=0; j<nii_ptr->nvox; j++) 
    avg[j] = input[j]/(double)nfiles;
  /* prepare sum of squares image */
  if(std_filename != NULL) {
    sum_squares  = (double *)malloc(sizeof(double)*nii_ptr->nvox);
    for (j=0; j<nii_ptr->nvox; j++) 
      sum_squares[j] = input[j]*input[j];
  }    
  free(input);
  
  /* read remaining images and check for image parameters */
  for (i=1; i<nfiles; i++) {
    fprintf(stdout,"%3d: %s\n",i+1, infiles[i]);
    nii_ptr2 = read_nifti_double(infiles[i], &input);
    if(nii_ptr2 == NULL) {
      fprintf(stderr,"Error reading %s.\n", infiles[i]);
      return(EXIT_FAILURE);
    }
    
    /* check for dimensions */
    equal_image_dimensions(nii_ptr,nii_ptr2);
    
    /* calculate average */
    for (j=0; j<nii_ptr->nvox; j++) 
      avg[j] += (input[j]/(double)nfiles);

    /* calculate sum of squares */
    if(std_filename != NULL) {
      for (j=0; j<nii_ptr->nvox; j++) 
        sum_squares[j] += input[j]*input[j];
    }
    free(input);
  }

  if(std_filename != NULL) {
    for (j=0; j<nii_ptr->nvox; j++) 
      sum_squares[j] = sqrt(1.0/((double)nfiles-1.0)*(sum_squares[j] - (double)nfiles*avg[j]*avg[j]));
    if (!write_nifti_double( std_filename, sum_squares, DT_FLOAT32, 1.0, dims, separations, nii_ptr)) 
      exit(EXIT_FAILURE);
    free(sum_squares);
  }

  if (!write_nifti_double( outfile, avg, DT_FLOAT32, 1.0, dims, separations, nii_ptr)) 
    exit(EXIT_FAILURE);

  free(avg);
  
  return(EXIT_SUCCESS);

}
