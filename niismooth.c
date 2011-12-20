/*
 * Christian Gaser
 * $Id$ 
 *
 */

#include <ParseArgv.h>
#include <float.h>
#include <stdlib.h>
#if !defined(_WIN32)
  #include <libgen.h>
#endif

#include "nifti/nifti1_io.h"
#include "nifti/nifti1_local.h"

extern nifti_image *read_nifti_double( const char *input_filename, double *image[]);
extern int write_nifti( const char *output_filename, double image[], int data_type, double slope, int dim[], double vox[], nifti_image *in_ptr);
extern int smooth_double(double *vol, int dims[3], double separations[3], double s[3], int use_mask);

double fwhm = 8.0;
int use_mask = 0;
int verbose = 0;

static
ArgvInfo argTable[] = {
  {"-fwhm", ARGV_FLOAT, (char *) 1, (char *) &fwhm, 
       "FWHM in mm."},
  {"-mask", ARGV_CONSTANT, (char *) 1, (char *) &use_mask,
       "Use masked smoothing (default no masking)."},
  {"-v", ARGV_CONSTANT, (char *) 1, (char *) &verbose,
       "Be verbose."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};


/* Main program */
int
main(int argc, char *argv[])
{
  char *infile, *outfile;
  int i, dims[3];
  double *input, separations[3], s[3];
  nifti_image *nii_ptr;

  /* Get arguments */
#if !defined(_WIN32)
  if  (ParseArgv(&argc, argv, argTable, 0) ||(argc < 2)) {
#else
  if  (ParseArgv(&argc, argv, argTable, 0) ||(argc < 3)) {
#endif
   (void) fprintf(stderr, 
   "\nUsage: %s [options] in.nii [out.nii]\n",
        argv[0]);
   (void) fprintf(stderr, 
   "    %s -help\n\n", argv[0]);
   exit(EXIT_FAILURE);
  }
  
  infile = argv[1];
  
  if (verbose)
    fprintf(stdout,"Filtering %s with FWHM of %gmm.\n", infile, fwhm);
  
  /* read first image to get image parameters */
  nii_ptr = read_nifti_double(infile, &input);
  if(nii_ptr == NULL) {
    fprintf(stderr,"Error reading %s.\n", infile);
    return(EXIT_FAILURE);
  }

  /* only allow isotropic filtering */
  for(i=0; i<3; i++) s[i] = fwhm;
  
  separations[0] = nii_ptr->dx;
  separations[1] = nii_ptr->dy;
  separations[2] = nii_ptr->dz;
  dims[0] = nii_ptr->nx;
  dims[1] = nii_ptr->ny;
  dims[2] = nii_ptr->nz;
  
  smooth_double(input, dims, separations, s, use_mask);
  
  /* if not defined use original name as basename for output */
  if(argc == 3)
    outfile = argv[2];
  else {
    outfile = argv[1];
    #if !defined(_WIN32)
      (void) sprintf(outfile, "%s/s%g%s", dirname(outfile), fwhm, basename(outfile)); 
    #endif
  }

  /* write data using same data type and rescale */
  if (!write_nifti( outfile, input, nii_ptr->datatype, 0.0, dims, separations, nii_ptr)) 
    exit(EXIT_FAILURE);
  
  return(EXIT_SUCCESS);

}
