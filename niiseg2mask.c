/*
 * Christian Gaser
 * $Id$ 
 *
 */

#include "ParseArgv.h"
#include "niilib.h"
#include "vollib.h"


static int usage(void)
{
    static const char msg[] = {
        "niiseg2mask: Create brainmask from segmentations\n"
        "usage: niiseg2mask GM WM CSF out.nii\n"
    };
    fprintf(stderr, "%s", msg);
    exit(EXIT_FAILURE);
}

int
main( int argc, char **argv )
{
  /* NIFTI stuff */
  nifti_image *src_ptr;
  char      *csf_filename, *gm_filename, *wm_filename;
  char      *basename, *extension, buffer[1024], filename[1024];
  int       i, j, dims[3], vol;
  unsigned char *probs, *mask;
  double    voxelsize[3];
  float     *src, *tmp;

  /* Get arguments */
  if (argc < 4) {
    (void) fprintf(stderr, 
    "\nUsage: %s GM WM CSF out.nii\n",
                     argv[0]);
    (void) fprintf(stderr, 
      "       %s -help\n\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  
  gm_filename  = argv[1];
  wm_filename  = argv[2];
  csf_filename = argv[3];

  /* deal with extension */
  extension = nifti_find_file_extension(gm_filename);
  
  /* if no valid extension was found use .nii */
  if (extension == NULL) {
    fprintf(stderr,"No valid extension found for output filename %s.\n",gm_filename);
    exit(EXIT_FAILURE);
  }

  /* read GM */
  src_ptr = read_nifti_float(gm_filename, &src, 0);
  if(src_ptr == NULL) {
    fprintf(stderr,"Error reading %s.\n", gm_filename);
    return(EXIT_FAILURE);
  }
  
  for (int j=0; j<3; j++) {
      voxelsize[j] = src_ptr->pixdim[j+1];
      dims[j] = src_ptr->dim[j+1];
  }

  vol = src_ptr->nvox;
  probs  = (unsigned char *)malloc(sizeof(unsigned char)*vol*3);
    
  if(probs == NULL) {
    fprintf(stderr,"Memory allocation error\n");
    exit(EXIT_FAILURE);
  }

  for (i = 0; i < src_ptr->nvox; i++)
    probs[i] = (unsigned char) round(255.0*src[i]);

  /* read WM */
  src_ptr = read_nifti_float(wm_filename, &src, 0);
  if(src_ptr == NULL) {
    fprintf(stderr,"Error reading %s.\n", wm_filename);
    return(EXIT_FAILURE);
  }
  
  for (i = 0; i < src_ptr->nvox; i++)
    probs[i + vol] = (unsigned char) round(255.0*src[i]);

  /* read CSF */
  src_ptr = read_nifti_float(csf_filename, &src, 0);
  if(src_ptr == NULL) {
    fprintf(stderr,"Error reading %s.\n", csf_filename);
    return(EXIT_FAILURE);
  }
  
  for (i = 0; i < vol; i++)
    probs[i + 2*vol] = (unsigned char) round(255.0*src[i]);
        
  double slope = 1.0/255.0;
  tmp = (float *)malloc(sizeof(float)*vol);
  mask = (unsigned char *)malloc(sizeof(unsigned char)*vol);

  cleanup(probs, mask, dims, voxelsize, 2, 1);


  for (int j = 0; j < 3; j++) {
    for (int i = 0; i < vol; i++) tmp[i] = (float)probs[i + j*vol];

    strcpy( buffer, argv[j+1]);
    basename = nifti_makebasename(buffer);
    (void) sprintf( filename, "%s_cleanup%s", basename, extension); 
    if(!write_nifti_float(filename, tmp, NIFTI_TYPE_UINT8, slope, 
                dims, voxelsize, src_ptr))
          exit(EXIT_FAILURE);
  }
  
  free(src);
  free(probs);
  
  return(EXIT_SUCCESS);
}