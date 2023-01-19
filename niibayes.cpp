
#include "niibayes.h"

void Usage(char *exec)
{
  fprintf(stderr,"niibayes\n");
  fprintf(stderr,"Usage:\t%s -target <targetImageName> -source <sourceImageName> -tpm <TPM>.\n",exec);
  fprintf(stderr,"\tSee the help for more details (-h).\n");
  return;
}


int main(int argc, char **argv)
{
  
  PARAM *param = (PARAM *)calloc(1,sizeof(PARAM));
  FLAG *flag = (FLAG *)calloc(1,sizeof(FLAG));
  float *tpm;
  int remove_sinus;
  char *basename, *extension, buffer[1024], filename[1024];
  double bias_fwhm = 100;
  double bias_rate = 0.9;
  
  flag->affineFlag = 1;
  flag->rigidFlag = 1;
  param->block_percent_to_use = 50;
  param->inlier_lts = 50;
  flag->alignCenterFlag = 1;
  flag->verbose = 0;
  
  param->tpmImageName = (char **)malloc(6);

  
  /* read the input parameter */
  for(int i=1;i<argc;i++){
    if(strcmp(argv[i], "-help")==0 || strcmp(argv[i], "-Help")==0 ||
       strcmp(argv[i], "-HELP")==0 || strcmp(argv[i], "-h")==0 ||
       strcmp(argv[i], "--h")==0 || strcmp(argv[i], "--help")==0){
      Usage(argv[0]);
      return 0;
    }
    else if(strcmp(argv[i], "-target") == 0){
      param->targetImageName=argv[++i];
      flag->targetImageFlag=1;
    }
    else if(strcmp(argv[i], "-source") == 0){
      param->sourceImageName=argv[++i];
      flag->sourceImageFlag=1;
    }
    else if(strcmp(argv[i], "-tpm") == 0){
      param->tpmImageName[0]=argv[++i];
      param->tpmImageName[1]=argv[++i];
      flag->tpmImageFlag=1;
    }
    else if(strcmp(argv[i], "-bias-fwhm") == 0){
      bias_fwhm=(float)(atof(argv[++i]));
    }
    else if(strcmp(argv[i], "-bias-rate") == 0){
      bias_rate=(float)(atof(argv[++i]));
    }
    else{
      fprintf(stderr,"Err:\tParameter %s unknown.\n",argv[i]);
      Usage(argv[0]);
      return 1;
    }
  }
  
  if(!flag->targetImageFlag || !flag->sourceImageFlag || !flag->tpmImageFlag){
    fprintf(stderr,"Error:\tThe target, source, and tpm images have to be defined.\n");
    Usage(argv[0]);
    return 1;
  }
  
  /* if no output is defined use source image name */
  if(!flag->outputMaskFlag)
    strcpy( buffer, param->sourceImageName);
  else
    strcpy( buffer, param->outputMaskName);

  /* get basename */
  basename = nifti_makebasename(buffer);

  /* deal with extension */
  extension = nifti_find_file_extension(buffer);

  /* if no valid extension was found use .nii */
  if (extension == NULL) {
    fprintf(stdout,"Use .nii as extension for %s.\n",buffer);
    strcpy(extension, ".nii");
  }

  if(!flag->levelNumberFlag) param->levelNumber=3;
  
  /* Read the maximum number of iteration */
  if(!flag->maxIterationFlag) param->maxIteration=5;

  if(!flag->level2PerformFlag) param->level2Perform=param->levelNumber;

  param->level2Perform=param->level2Perform<param->levelNumber?param->level2Perform:param->levelNumber;
  
  /* Read the target and source images */
  nifti_image *sourceImage = nifti_image_read(param->sourceImageName,true);
  if(sourceImage->data == NULL){
    fprintf(stderr, "** ERROR Error when reading the source image: %s\n", param->sourceImageName);
    return 0;
  }
  reg_changeDatatype<PrecisionTYPE>(sourceImage);

  nifti_image *tpmImage = read_nifti_float(param->tpmImageName[0],&tpm,true);
  if(tpmImage->data == NULL){
    fprintf(stderr, "** ERROR Error when reading the source image: %s\n", param->tpmImageName[0]);
    return 0;
  }
  reg_changeDatatype<PrecisionTYPE>(tpmImage);

  /* check number of classes */
  if(tpmImage->nt != 6){
    fprintf(stderr, "** ERROR TPM image has %d instead of 6 classes: %s\n", tpmImage->nt, param->tpmImageName[0]);
    return 0;
  }

  nifti_image *targetImage = nifti_image_read(param->targetImageName,true);   
  if(targetImage->data == NULL){
    fprintf(stderr, "** ERROR Error when reading the target image: %s\n", param->targetImageName);
    return 0;
  }
  reg_changeDatatype<PrecisionTYPE>(targetImage);

  nifti_image *positionFieldImage = nifti_copy_nim_info(sourceImage);
  positionFieldImage->dim[0]=positionFieldImage->ndim=5;
  positionFieldImage->dim[1]=positionFieldImage->nx=sourceImage->nx;
  positionFieldImage->dim[2]=positionFieldImage->ny=sourceImage->ny;
  positionFieldImage->dim[3]=positionFieldImage->nz=sourceImage->nz;
  positionFieldImage->dim[4]=positionFieldImage->nt=1;
  positionFieldImage->dim[5]=positionFieldImage->nu=3;
  positionFieldImage->dim[6]=positionFieldImage->nv=1;
  positionFieldImage->dim[7]=positionFieldImage->nw=1;
  positionFieldImage->pixdim[4]=positionFieldImage->dt=1.0;
  positionFieldImage->pixdim[5]=positionFieldImage->du=1.0;
  positionFieldImage->pixdim[6]=positionFieldImage->dv=1.0;
  positionFieldImage->pixdim[7]=positionFieldImage->dw=1.0;
  positionFieldImage->nvox=positionFieldImage->nx*positionFieldImage->ny*positionFieldImage->nz*positionFieldImage->nt*positionFieldImage->nu;
  if(sizeof(PrecisionTYPE)==4) positionFieldImage->datatype = NIFTI_TYPE_FLOAT32;
  else positionFieldImage->datatype = NIFTI_TYPE_FLOAT64;
  positionFieldImage->nbyper = sizeof(PrecisionTYPE);
  positionFieldImage->data = (void *)calloc(positionFieldImage->nvox, positionFieldImage->nbyper);

  mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));
  
  /* get affine transformation and its inverse */
  affineTransformation = affineRegistration(param, flag);
  *affineTransformation = nifti_mat44_inverse(*affineTransformation);

  reg_affine_positionField(affineTransformation,
            sourceImage,
            positionFieldImage);

  /* allocate the result image */
  nifti_image *resultImage = nifti_copy_nim_info(sourceImage);
  resultImage->cal_min = 0;
  resultImage->cal_max = 0;
  resultImage->scl_slope = 1.0;
  resultImage->scl_inter = 0.0;
  resultImage->datatype = tpmImage->datatype;
  resultImage->nbyper = tpmImage->nbyper;
  resultImage->data = (void *)malloc(resultImage->nvox*resultImage->nbyper);

  unsigned char *priors = (unsigned char *)malloc(sizeof(unsigned char)*sourceImage->nvox*tpmImage->nt);

  for (int j=0; j<tpmImage->nt; j++) {
  
    float_to_ptr_4D(tpmImage, tpm, j);            

    reg_resampleSourceImage<PrecisionTYPE>(sourceImage,
          tpmImage,
          resultImage,
          positionFieldImage,
          NULL,
          1,
          0);
          
    ptr_to_uint8_4D(resultImage, priors, j);

  }

  free(positionFieldImage);
  free(resultImage);
  free(targetImage);

  unsigned char *label = (unsigned char *)malloc(sizeof(unsigned char)*sourceImage->nvox);
  unsigned char *probs = (unsigned char *)malloc(sizeof(unsigned char)*sourceImage->nvox*tpmImage->nt);
  float           *src = (float *)malloc(sizeof(float)*sourceImage->nvox);
  
  ptr_to_float_4D(sourceImage, src, 0);

  double voxelsize[3];
  int dims[3];
  
  for (int j=0; j<3; j++) {
    voxelsize[j] = sourceImage->pixdim[j+1];
    dims[j] = sourceImage->dim[j+1];
  }
  
  double slope = 1.0/255.0;
  int vol = dims[0]*dims[1]*dims[2];

  int do_warp = 1;
  int correct_nu = 1;

  Bayes(src, label, priors, probs, voxelsize, dims, bias_fwhm, bias_rate, do_warp);

  int cleanup_strength = 2;
  remove_sinus = 1;
//  initial_cleanup(probs, label, dims, voxelsize, cleanup_strength, remove_sinus);
//  cleanup(probs, label, dims, voxelsize, 2, 1);
    
  float *tmp = (float *)malloc(sizeof(float)*sourceImage->nvox);
  
  for (int j=0; j<3; j++) {
    (void) sprintf( filename, "%s_prob%d%s",basename, j+1,extension); 
//    (void) sprintf( filename, "%s_prob%d_%g_%g%s",basename, j+1,bias_fwhm,bias_rate,extension); 
    for (int i = 0; i < sourceImage->nvox; i++) tmp[i] = (float)probs[i + j*vol];
    if(!write_nifti_float(filename, tmp, NIFTI_TYPE_UINT8, slope, 
                dims, voxelsize, sourceImage))
      exit(EXIT_FAILURE);
  }

  (void) sprintf( filename, "%s_bias%s", basename, extension); 
  if(!write_nifti_float(filename, src, NIFTI_TYPE_FLOAT32, 1.0, 
              dims, voxelsize, sourceImage))
    exit(EXIT_FAILURE);

if (0) {
  for (int i = 0; i < sourceImage->nvox; i++) tmp[i] = (float)label[i];
  (void) sprintf( filename, "%s_label%s", basename, extension); 
  if(!write_nifti_float(filename, tmp, NIFTI_TYPE_UINT8, slope, 
              dims, voxelsize, sourceImage))
    exit(EXIT_FAILURE);
    
  for (int i = 0; i < sourceImage->nvox; i++) tmp[i] = (float)priors[i] + (float)priors[i + vol] + (float)priors[i + 2*vol];
  
  (void) sprintf( filename, "%s_mask%s",basename, extension); 
  if(!write_nifti_float(filename, tmp, NIFTI_TYPE_UINT8, slope, 
              dims, voxelsize, sourceImage))
    exit(EXIT_FAILURE);
      
}
  free(tmp);
  free(tpmImage);
  free(sourceImage);
  free(label);
  free(src);
  free(probs);
  free(priors);

  return 0;
}

