/*
 * Christian Gaser
 * $Id: affineRegistration.cpp 191 2012-07-31 15:41:36Z gaser $ 
 *
 */

#include "niibayes.h"

#define CONVERGENCE_EPS 0.00001
#define NDEBUG 1

bool reg_test_convergence(mat44 *mat)
{
  bool convergence=true;
  if((fabsf(mat->m[0][0])-1.0f)>CONVERGENCE_EPS) convergence=false;
  if((fabsf(mat->m[1][1])-1.0f)>CONVERGENCE_EPS) convergence=false;
  if((fabsf(mat->m[2][2])-1.0f)>CONVERGENCE_EPS) convergence=false;

  if((fabsf(mat->m[0][1])-0.0f)>CONVERGENCE_EPS) convergence=false;
  if((fabsf(mat->m[0][2])-0.0f)>CONVERGENCE_EPS) convergence=false;
  if((fabsf(mat->m[0][3])-0.0f)>CONVERGENCE_EPS) convergence=false;

  if((fabsf(mat->m[1][0])-0.0f)>CONVERGENCE_EPS) convergence=false;
  if((fabsf(mat->m[1][2])-0.0f)>CONVERGENCE_EPS) convergence=false;
  if((fabsf(mat->m[1][3])-0.0f)>CONVERGENCE_EPS) convergence=false;

  if((fabsf(mat->m[2][0])-0.0f)>CONVERGENCE_EPS) convergence=false;
  if((fabsf(mat->m[2][1])-0.0f)>CONVERGENCE_EPS) convergence=false;
  if((fabsf(mat->m[2][3])-0.0f)>CONVERGENCE_EPS) convergence=false;

  return convergence;
}

mat44 *affineRegistration(PARAM *param, FLAG *flag)
{

  /* Read the target and source images */
  nifti_image *targetHeader = nifti_image_read(param->targetImageName,false);
  if(targetHeader == NULL){
    fprintf(stderr,"** ERROR Error when reading the target image: %s\n",param->targetImageName);
    return NULL;
  }
  nifti_image *sourceHeader = nifti_image_read(param->sourceImageName,false);
  if(sourceHeader == NULL){
    fprintf(stderr,"** ERROR Error when reading the source image: %s\n",param->sourceImageName);
    return NULL;
  }

  /* Flag for 2D registration */
    if(sourceHeader->nz==1 || targetHeader->nz==1){
        flag->twoDimRegistration=1;
    }

  /* Check the source background index */
  if(!flag->backgroundIndexFlag) param->sourceBGValue = 0.0;
  else{
    if(param->backgroundIndex[0] < 0 || param->backgroundIndex[1] < 0 || param->backgroundIndex[2] < 0 
       || param->backgroundIndex[0] >= sourceHeader->dim[1] || param->backgroundIndex[1] >= sourceHeader->dim[2] || param->backgroundIndex[2] >= sourceHeader->dim[3]){
      fprintf(stderr,"The specified index (%i %i %i) for background does not belong to the source image (out of bondary)\n",
          param->backgroundIndex[0], param->backgroundIndex[1], param->backgroundIndex[2]);
      return NULL;
    }
  }
    
  /* Read the input affine tranformation is defined otherwise assign it to identity */
  mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));
  affineTransformation->m[0][0]=1.0;
  affineTransformation->m[1][1]=1.0;
  affineTransformation->m[2][2]=1.0;
  affineTransformation->m[3][3]=1.0;
  if(flag->alignCenterFlag){
    mat44 *sourceMatrix;
    if(sourceHeader->sform_code>0)
      sourceMatrix = &(sourceHeader->sto_xyz);
    else sourceMatrix = &(sourceHeader->qto_xyz);
    mat44 *targetMatrix;
    if(targetHeader->sform_code>0)
      targetMatrix = &(targetHeader->sto_xyz);
    else targetMatrix = &(targetHeader->qto_xyz);
    float sourceCenter[3];
    sourceCenter[0]=(float)(sourceHeader->nx)/2.0f;
    sourceCenter[1]=(float)(sourceHeader->ny)/2.0f;
    sourceCenter[2]=(float)(sourceHeader->nz)/2.0f;
    float targetCenter[3];
    targetCenter[0]=(float)(targetHeader->nx)/2.0f;
    targetCenter[1]=(float)(targetHeader->ny)/2.0f;
    targetCenter[2]=(float)(targetHeader->nz)/2.0f;
    float sourceRealPosition[3]; reg_mat44_mul(sourceMatrix, sourceCenter, sourceRealPosition);
    float targetRealPosition[3]; reg_mat44_mul(targetMatrix, targetCenter, targetRealPosition);
    affineTransformation->m[0][3]=sourceRealPosition[0]-targetRealPosition[0];
    affineTransformation->m[1][3]=sourceRealPosition[1]-targetRealPosition[1];
    affineTransformation->m[2][3]=sourceRealPosition[2]-targetRealPosition[2];
  }
  
    /* read and binarise the target mask image */
    nifti_image *targetMaskImage=NULL;
    if(flag->targetMaskFlag){
        targetMaskImage = nifti_image_read(param->targetMaskName,true);
        if(targetMaskImage == NULL){
            fprintf(stderr,"* ERROR Error when reading the target mask image: %s\n",param->targetMaskName);
            return NULL;
        }
        /* check the dimension */
        for(int i=1; i<=targetHeader->dim[0]; i++){
            if(targetHeader->dim[i]!=targetMaskImage->dim[i]){
                fprintf(stderr,"* ERROR The target image and its mask do not have the same dimension\n");
                return NULL;
            }
        }
        reg_tool_binarise_image(targetMaskImage);
    }

  /* *********************************** */
  /* DISPLAY THE REGISTRATION PARAMETERS */
  /* *********************************** */
  if (flag->verbose ){
    printf("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    printf("\n\n");
    printf("Parameters\n");
    printf("Target image name: %s\n",targetHeader->fname);
    printf("\t%ix%ix%i voxels\n",targetHeader->nx,targetHeader->ny,targetHeader->nz);
    printf("\t%gx%gx%g mm\n",targetHeader->dx,targetHeader->dy,targetHeader->dz);
    printf("Source image name: %s\n",sourceHeader->fname);
    printf("\t%ix%ix%i voxels\n",sourceHeader->nx,sourceHeader->ny,sourceHeader->nz);
    printf("\t%gx%gx%g mm\n",sourceHeader->dx,sourceHeader->dy,sourceHeader->dz);
    printf("Maximum iteration number: %i ",param->maxIteration);
    printf("(%i during the first level)\n",2*param->maxIteration);
    printf("Percentage of blocks: %i %%",param->block_percent_to_use);
    printf(" (100%% during the first level)\n");
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n\n");
  }
    
  /* ********************** */
  /* START THE REGISTRATION */
  /* ********************** */
  
  for(int level=0; level<param->level2Perform; level++){
    /* Read the target and source image */
    nifti_image *targetImage = nifti_image_read(param->targetImageName,true);
    
    if(targetImage->data == NULL){
      fprintf(stderr, "** ERROR Error when reading the target image: %s\n", param->targetImageName);
      return NULL;
    }
    reg_changeDatatype<PrecisionTYPE>(targetImage);
    nifti_image *sourceImage = nifti_image_read(param->sourceImageName,true);
    if(sourceImage->data == NULL){
      fprintf(stderr, "** ERROR Error when reading the source image: %s\n", param->sourceImageName);
      return NULL;
    }
    reg_changeDatatype<PrecisionTYPE>(sourceImage);

    // Twice more iterations are performed during the first level
    // All the blocks are used during the first level
    int maxNumberOfIterationToPerform=param->maxIteration;
    int percentageOfBlockToUse=param->block_percent_to_use;
    if(level==0){
        maxNumberOfIterationToPerform*=2;
        percentageOfBlockToUse=100;
    }

    /* declare the target mask array */
    int *targetMask;
    int activeVoxelNumber=0;

    /* downsample the input images if appropriate */
    nifti_image *tempMaskImage=NULL;
    if(flag->targetMaskFlag){
        tempMaskImage = nifti_copy_nim_info(targetMaskImage);
        tempMaskImage->data = (void *)malloc(tempMaskImage->nvox * tempMaskImage->nbyper);
        memcpy(tempMaskImage->data, targetMaskImage->data, tempMaskImage->nvox*tempMaskImage->nbyper);
    }

    for(int l=level; l<param->levelNumber-1; l++){
        int ratio = (int)powf(2.0f,param->levelNumber-param->levelNumber+l+1.0f);

        bool sourceDownsampleAxis[8]={true,true,true,true,true,true,true,true};
        if((sourceHeader->nx/ratio) < 32) sourceDownsampleAxis[1]=false;
        if((sourceHeader->ny/ratio) < 32) sourceDownsampleAxis[2]=false;
        if((sourceHeader->nz/ratio) < 32) sourceDownsampleAxis[3]=false;
        reg_downsampleImage<PrecisionTYPE>(sourceImage, 1, sourceDownsampleAxis);

        bool targetDownsampleAxis[8]={true,true,true,true,true,true,true,true};
        if((targetHeader->nx/ratio) < 32) targetDownsampleAxis[1]=false;
        if((targetHeader->ny/ratio) < 32) targetDownsampleAxis[2]=false;
        if((targetHeader->nz/ratio) < 32) targetDownsampleAxis[3]=false;
        reg_downsampleImage<PrecisionTYPE>(targetImage, 1, targetDownsampleAxis);

        if(flag->targetMaskFlag){
            reg_downsampleImage<PrecisionTYPE>(tempMaskImage, 0, targetDownsampleAxis);
        }
    }
    targetMask = (int *)malloc(targetImage->nvox*sizeof(int));
    if(flag->targetMaskFlag){
        reg_tool_binaryImage2int(tempMaskImage, targetMask, activeVoxelNumber);
        nifti_image_free(tempMaskImage);
    }
    else{
        for(unsigned int i=0; i<targetImage->nvox; i++)
            targetMask[i]=i;
        activeVoxelNumber=targetImage->nvox;
    }


    /* smooth the input image if appropriate */
    if(flag->targetSigmaFlag){
            bool smoothAxis[8]={true,true,true,true,true,true,true,true};
      reg_gaussianSmoothing<PrecisionTYPE>(targetImage, param->targetSigmaValue, smoothAxis);
        }
    if(flag->sourceSigmaFlag){
            bool smoothAxis[8]={true,true,true,true,true,true,true,true};
      reg_gaussianSmoothing<PrecisionTYPE>(sourceImage, param->sourceSigmaValue, smoothAxis);
        }
    
    /* allocate the deformation Field image */
    nifti_image *positionFieldImage = nifti_copy_nim_info(targetImage);
    positionFieldImage->dim[0]=positionFieldImage->ndim=5;
    positionFieldImage->dim[1]=positionFieldImage->nx=targetImage->nx;
    positionFieldImage->dim[2]=positionFieldImage->ny=targetImage->ny;
    positionFieldImage->dim[3]=positionFieldImage->nz=targetImage->nz;
    positionFieldImage->dim[4]=positionFieldImage->nt=1;positionFieldImage->pixdim[4]=positionFieldImage->dt=1.0;
    if(flag->twoDimRegistration) positionFieldImage->dim[5]=positionFieldImage->nu=2;
    else positionFieldImage->dim[5]=positionFieldImage->nu=3;
    positionFieldImage->pixdim[5]=positionFieldImage->du=1.0;
    positionFieldImage->dim[6]=positionFieldImage->nv=1;positionFieldImage->pixdim[6]=positionFieldImage->dv=1.0;
    positionFieldImage->dim[7]=positionFieldImage->nw=1;positionFieldImage->pixdim[7]=positionFieldImage->dw=1.0;
    positionFieldImage->nvox=positionFieldImage->nx*positionFieldImage->ny*positionFieldImage->nz*positionFieldImage->nt*positionFieldImage->nu;
    if(sizeof(PrecisionTYPE)==4) positionFieldImage->datatype = NIFTI_TYPE_FLOAT32;
    else positionFieldImage->datatype = NIFTI_TYPE_FLOAT64;
    positionFieldImage->nbyper = sizeof(PrecisionTYPE);
    positionFieldImage->data = (void *)calloc(positionFieldImage->nvox, positionFieldImage->nbyper);
    
    /* allocate the result image */
    nifti_image *resultImage = nifti_copy_nim_info(targetImage);
    resultImage->datatype = sourceImage->datatype;
    resultImage->nbyper = sourceImage->nbyper;
    resultImage->data = (void *)calloc(resultImage->nvox, resultImage->nbyper);

    /* Set the padding value */
    if(flag->backgroundIndexFlag){
      int index[3];
      index[0]=param->backgroundIndex[0];
      index[1]=param->backgroundIndex[1];
      index[2]=param->backgroundIndex[2];
      if(flag->pyramidFlag){
        for(int l=level; l<param->levelNumber-1; l++){
          index[0] /= 2;
          index[1] /= 2;
          index[2] /= 2;
        }
      }
      param->sourceBGValue = (float)(reg_tool_GetIntensityValue(sourceImage, index));
    }
    else param->sourceBGValue = 0;

    /* initialise the block matching */
    _reg_blockMatchingParam blockMatchingParams;
        initialise_block_matching_method(   targetImage,
                                            &blockMatchingParams,
                                            percentageOfBlockToUse,    // percentage of block kept
                                            param->inlier_lts,              // percentage of inlier in the optimisation process
                                            targetMask,
                                            flag->useGPUFlag);

    mat44 updateAffineMatrix;


    /* Display some parameters specific to the current level */
    if (flag->verbose ){
      printf("Current level %i / %i\n", level+1, param->levelNumber);
      printf("Target image size: \t%ix%ix%i voxels\t%gx%gx%g mm\n",
             targetImage->nx, targetImage->ny, targetImage->nz, targetImage->dx, targetImage->dy, targetImage->dz);
      printf("Source image size: \t%ix%ix%i voxels\t%gx%gx%g mm\n",
             sourceImage->nx, sourceImage->ny, sourceImage->nz, sourceImage->dx, sourceImage->dy, sourceImage->dz);
          if(flag->twoDimRegistration)
              printf("Block size = [4 4 1]\n");
          else printf("Block size = [4 4 4]\n");
      printf("Block number = [%i %i %i]\n", blockMatchingParams.blockNumber[0],
        blockMatchingParams.blockNumber[1], blockMatchingParams.blockNumber[2]);
    }
#ifndef NDEBUG
    if(targetImage->sform_code>0)
      reg_mat44_disp(&targetImage->sto_xyz, (char *)"[DEBUG] Target image matrix (sform sto_xyz)");
    else reg_mat44_disp(&targetImage->qto_xyz, (char *)"[DEBUG] Target image matrix (qform qto_xyz)");
    if(sourceImage->sform_code>0)
      reg_mat44_disp(&sourceImage->sto_xyz, (char *)"[DEBUG] Source image matrix (sform sto_xyz)");
    else reg_mat44_disp(&sourceImage->qto_xyz, (char *)"[DEBUG] Source image matrix (qform qto_xyz)");
#endif
    if (flag->verbose ){
      printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
      reg_mat44_disp(affineTransformation, (char *)"Initial affine transformation:");
    }

    /* ****************** */
    /* Rigid registration */
    /* ****************** */
    int iteration=0;

    if((flag->rigidFlag && !flag->affineFlag) || (flag->affineFlag && flag->rigidFlag && level==0)){
      while(iteration<maxNumberOfIterationToPerform){
        /* Compute the affine transformation deformation field */
          reg_affine_positionField( affineTransformation,
                                                targetImage,
                                                positionFieldImage);
          /* Resample the source image */
          reg_resampleSourceImage<PrecisionTYPE>( targetImage,
                                                            sourceImage,
                                                            resultImage,
                                                            positionFieldImage,
                                                            targetMask,
                                                            1,
                                                            param->sourceBGValue);
          /* Compute the correspondances between blocks */
                    block_matching_method<PrecisionTYPE>(   targetImage,
                                                            resultImage,
                                                            &blockMatchingParams,
                                                            targetMask);
          /* update  the affine transformation matrix */
          optimize( &blockMatchingParams,
                  &updateAffineMatrix,
                  RIGID);
        // the affine transformation is updated
        *affineTransformation = reg_mat44_mul( affineTransformation, &(updateAffineMatrix));
#ifndef NDEBUG
        printf("[DEBUG] -Rigid- iteration %i - ",iteration);
        reg_mat44_disp(&updateAffineMatrix, (char *)"[DEBUG] updateMatrix");
        reg_mat44_disp(affineTransformation, (char *)"[DEBUG] updated affine");
#endif
  
        if(reg_test_convergence(&updateAffineMatrix)) break;
        iteration++;
      }
    }

    /* ******************* */
    /* Affine registration */
    /* ******************* */
    iteration=0;
    if(flag->affineFlag){
      while(iteration<maxNumberOfIterationToPerform){
        /* Compute the affine transformation deformation field */
          reg_affine_positionField( affineTransformation,
                  targetImage,
                  positionFieldImage);
          /* Resample the source image */
          reg_resampleSourceImage<PrecisionTYPE>( targetImage,
                    sourceImage,
                    resultImage,
                    positionFieldImage,
                                        targetMask,
                    1,
                    param->sourceBGValue);
          /* Compute the correspondances between blocks */
          block_matching_method<PrecisionTYPE>( targetImage,
                    resultImage,
                    &blockMatchingParams,
                                        targetMask);
          /* update  the affine transformation matrix */
          optimize( &blockMatchingParams,
              &updateAffineMatrix,
              AFFINE);
  
        // the affine transformation is updated
        *affineTransformation = reg_mat44_mul( affineTransformation, &(updateAffineMatrix));
#ifndef NDEBUG
        printf("[DEBUG] -Affine- iteration %i - ",iteration);
        reg_mat44_disp(&updateAffineMatrix, (char *)"[DEBUG] updateMatrix");
        reg_mat44_disp(affineTransformation, (char *)"[DEBUG] updated affine");
#endif
        if(reg_test_convergence(&updateAffineMatrix)) break;
        iteration++;
      }
    }

    free(targetMask);

    if(level==(param->level2Perform-1)){
      /* ****************** */
      /* OUTPUT THE RESULTS */
      /* ****************** */

      if(param->level2Perform != param->levelNumber){
          if(positionFieldImage->data)free(positionFieldImage->data);
          positionFieldImage->dim[1]=positionFieldImage->nx=targetHeader->nx;
          positionFieldImage->dim[2]=positionFieldImage->ny=targetHeader->ny;
          positionFieldImage->dim[3]=positionFieldImage->nz=targetHeader->nz;
          positionFieldImage->dim[4]=positionFieldImage->nt=1;positionFieldImage->pixdim[4]=positionFieldImage->dt=1.0;
          if(flag->twoDimRegistration)
              positionFieldImage->dim[5]=positionFieldImage->nu=2;
          else positionFieldImage->dim[5]=positionFieldImage->nu=3;
          positionFieldImage->pixdim[5]=positionFieldImage->du=1.0;
          positionFieldImage->dim[6]=positionFieldImage->nv=1;positionFieldImage->pixdim[6]=positionFieldImage->dv=1.0;
          positionFieldImage->dim[7]=positionFieldImage->nw=1;positionFieldImage->pixdim[7]=positionFieldImage->dw=1.0;
          positionFieldImage->nvox=positionFieldImage->nx*positionFieldImage->ny*positionFieldImage->nz*positionFieldImage->nt*positionFieldImage->nu;
          positionFieldImage->data = (void *)calloc(positionFieldImage->nvox, positionFieldImage->nbyper);
      }
      
      /* The corresponding deformation field is evaluated and saved */
      reg_affine_positionField( affineTransformation,
              targetHeader,
              positionFieldImage);
      
      /* The result image is resampled using a cubic spline interpolation */
      nifti_image_free(sourceImage);
      sourceImage = nifti_image_read(param->sourceImageName,true); // reload the source image with the correct intensity values
      nifti_image_free(resultImage);
      resultImage = nifti_copy_nim_info(targetHeader);
      resultImage->cal_min=sourceImage->cal_min;
      resultImage->cal_max=sourceImage->cal_max;
      resultImage->scl_slope=sourceImage->scl_slope;
      resultImage->scl_inter=sourceImage->scl_inter;
      resultImage->datatype = sourceImage->datatype;
      resultImage->nbyper = sourceImage->nbyper;
      resultImage->data = (void *)calloc(resultImage->nvox, resultImage->nbyper);
      reg_resampleSourceImage<double>(targetHeader,
              sourceImage,
              resultImage,
              positionFieldImage,
                            NULL,
              3,
              param->sourceBGValue);
      if(flag->outputResultFlag) {
          nifti_set_filenames(resultImage, param->outputResultName, 0, 0);
          nifti_image_write(resultImage);
      }

    }
    nifti_image_free(positionFieldImage);
    nifti_image_free(resultImage);
    nifti_image_free(targetImage);
    nifti_image_free(sourceImage);
    if(flag->verbose)
      reg_mat44_disp(affineTransformation, (char *)"Final affine transformation:");
#ifndef NDEBUG
    mat33 tempMat;
    for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
        tempMat.m[i][j] = affineTransformation->m[i][j];
      }
    }
    printf("[DEBUG] Matrix determinant %g\n", nifti_mat33_determ  (tempMat));
#endif
    if(flag->verbose)
      printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
        
  }

  nifti_image_free(targetHeader);
  nifti_image_free(sourceHeader);

  free(flag);
  free(param);

  return affineTransformation;
}
