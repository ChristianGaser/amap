/*
 *  reg_resample.cpp
 *
 *
 *  Created by Marc Modat on 18/05/2009.
 *  Copyright (c) 2009, University College London. All rights reserved.
 *  Centre for Medical Image Computing (CMIC)
 *  See the LICENSE.txt file in the nifty_reg root folder
 *
 */

#ifndef _MM_RESAMPLE_CPP
#define _MM_RESAMPLE_CPP

#include "_reg_resampling.h"
#include "_reg_affineTransformation.h"
#include "_reg_bspline.h"
#include "_reg_bspline_comp.h"
#include "_reg_tools.h"

#define PrecisionTYPE float

typedef struct{
	char *targetImageName;
	char *sourceImageName;
	char *affineMatrixName;
	char *inputCPPName;
	char *inputVelocityFieldName;
	char *outputPosName;
	char *outputDispName;
	char *outputEuclDispName;
	char *outputResultName;
	char *outputBlankName;
    char *outputJacobianName;
    char *outputJacobianMatrixName;
    char *outputCPPName;
	PrecisionTYPE sourceBGValue;
}PARAM;
typedef struct{
	bool targetImageFlag;
	bool sourceImageFlag;
	bool affineMatrixFlag;
	bool affineFlirtFlag;
	bool inputCPPFlag;
    bool inputVelocityFieldFlag;
	bool outputDispFlag;
	bool outputEuclDispFlag;
	bool outputPosFlag;
	bool outputFullDefFlag;
	bool outputResultFlag;
	bool outputBlankFlag;
	bool outputJacobianFlag;
	bool outputJacobianMatrixFlag;
    bool NNInterpolationFlag;
    bool TRIInterpolationFlag;
    bool outputCPPFlag;
}FLAG;


void PetitUsage(char *exec)
{
	fprintf(stderr,"Usage:\t%s -target <targetImageName> -source <sourceImageName> [OPTIONS].\n",exec);
	fprintf(stderr,"\tSee the help for more details (-h).\n");
	return;
}
void Usage(char *exec)
{
	printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
	printf("Usage:\t%s -target <filename> -source <filename> [OPTIONS].\n",exec);
	printf("\t-target <filename>\tFilename of the target image (mandatory)\n");
	printf("\t-source <filename>\tFilename of the source image (mandatory)\n\n");

	printf("* * OPTIONS * *\n");
    printf("\t*\tOnly one of the following tranformation is taken into account\n");
    printf("\t-aff <filename>\t\tFilename which contains an affine transformation (Affine*Target=Source)\n");
    printf("\t-affFlirt <filename>\t\tFilename which contains a radiological flirt affine transformation\n");
    printf("\t-cpp <filename>\t\tFilename of control point grid image\n");
    printf("\t-vel <filename>\t\tFilename of the velocity field image\n\n");

    printf("\t*\tThere are no limit for the required output number from the following\n");
    printf("\t-result <filename> \tFilename of the resampled image [none]\n");
	printf("\t-blank <filename> \tFilename of the resampled blank grid [none]\n");
	printf("\t-jac <filename> \tFilename of the Jacobian map image [none]\n");
	printf("\t-jacM <filename> \tFilename of the Jacobian matrix image [none]\n");
	printf("\t-opf <filename>\t\tFilename of the position field image\n");
	printf("\t-odf <filename>\t\tFilename of the displacement field image\n");
	printf("\t-oed <filename>\t\tFilename of the euclidian displacement image\n\n");
    printf("\t-ocp <filename> \tFilename of the control point grid image\n");

    printf("\t*\tOthers\n");
	printf("\t-NN \t\t\tUse a Nearest Neighbor interpolation for the source resampling (cubic spline by default)\n");
	printf("\t-TRI \t\t\tUse a Trilinear interpolation for the source resampling (cubic spline by default)\n");
	printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
	return;
}

int main(int argc, char **argv)
{
	PARAM *param = (PARAM *)calloc(1,sizeof(PARAM));
	FLAG *flag = (FLAG *)calloc(1,sizeof(FLAG));
	
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
		else if(strcmp(argv[i], "-aff") == 0){
			param->affineMatrixName=argv[++i];
			flag->affineMatrixFlag=1;
		}
		else if(strcmp(argv[i], "-affFlirt") == 0){
			param->affineMatrixName=argv[++i];
			flag->affineMatrixFlag=1;
			flag->affineFlirtFlag=1;
		}
		else if(strcmp(argv[i], "-result") == 0){
			param->outputResultName=argv[++i];
			flag->outputResultFlag=1;
		}
		else if(strcmp(argv[i], "-jac") == 0){
			param->outputJacobianName=argv[++i];
			flag->outputJacobianFlag=1;
		}
		else if(strcmp(argv[i], "-vel") == 0){
			param->inputVelocityFieldName=argv[++i];
            flag->inputVelocityFieldFlag=1;
		}
		else if(strcmp(argv[i], "-jacM") == 0){
			param->outputJacobianMatrixName=argv[++i];
			flag->outputJacobianMatrixFlag=1;
		}
		else if(strcmp(argv[i], "-cpp") == 0){
			param->inputCPPName=argv[++i];
			flag->inputCPPFlag=1;
		}
		else if(strcmp(argv[i], "-odf") == 0){
			param->outputDispName=argv[++i];
			flag->outputDispFlag=1;
		}
		else if(strcmp(argv[i], "-oed") == 0){
			param->outputEuclDispName=argv[++i];
			flag->outputEuclDispFlag=1;
		}
		else if(strcmp(argv[i], "-opf") == 0){
			param->outputPosName=argv[++i];
			flag->outputPosFlag=1;
		}
        else if(strcmp(argv[i], "-ocp") == 0){
            param->outputCPPName=argv[++i];
            flag->outputCPPFlag=1;
        }
		else if(strcmp(argv[i], "-NN") == 0){
			flag->NNInterpolationFlag=1;
		}
		else if(strcmp(argv[i], "-TRI") == 0){
			flag->TRIInterpolationFlag=1;
		}
		else if(strcmp(argv[i], "-blank") == 0){
			param->outputBlankName=argv[++i];
			flag->outputBlankFlag=1;
		}
		else{
			fprintf(stderr,"Err:\tParameter %s unknown.\n",argv[i]);
			PetitUsage(argv[0]);
			return 1;
		}
	}
	
	if(!flag->targetImageFlag || !flag->sourceImageFlag){
		fprintf(stderr,"Err:\tThe target and the source image have both to be defined.\n");
		PetitUsage(argv[0]);
		return 1;
	}
	
    /* Check the number of input images */
    if( ((unsigned int)flag->affineMatrixFlag
        + (unsigned int)flag->affineFlirtFlag 
        + (unsigned int)flag->inputCPPFlag
        + (unsigned int)flag->inputVelocityFieldFlag) > 1){
        fprintf(stderr,"Err:\tOnly one input transformation has to be assigned.\n");
        PetitUsage(argv[0]);
    }

	/* Read the target image */
	nifti_image *targetImage = nifti_image_read(param->targetImageName,true);
	if(targetImage == NULL){
		fprintf(stderr,"** ERROR Error when reading the target image: %s\n",param->targetImageName);
		return 1;
	}
	
	/* Read the source image */
	nifti_image *sourceImage = nifti_image_read(param->sourceImageName,true);
	if(sourceImage == NULL){
		fprintf(stderr,"** ERROR Error when reading the source image: %s\n",param->sourceImageName);
		return 1;
	}
	
	
	/* *********************************** */
	/* DISPLAY THE REGISTRATION PARAMETERS */
	/* *********************************** */
	printf("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
	printf("Command line:\n");
	for(int i=0;i<argc;i++) printf(" %s", argv[i]);
	printf("\n\n");
	printf("Parameters\n");
	printf("Target image name: %s\n",targetImage->fname);
	printf("\t%ix%ix%i voxels\n",targetImage->nx,targetImage->ny,targetImage->nz);
	printf("\t%gx%gx%g mm\n",targetImage->dx,targetImage->dy,targetImage->dz);
	printf("Source image name: %s\n",sourceImage->fname);
	printf("\t%ix%ix%i voxels\n",sourceImage->nx,sourceImage->ny,sourceImage->nz);
	printf("\t%gx%gx%g mm\n",sourceImage->dx,sourceImage->dy,sourceImage->dz);
	printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n\n");
	
	/* ********************** */
	/*  START THE RESAMPLING  */
	/* ********************** */
	
	/* allocate the position field image if necessary */
	bool positionFieldNeeded=false;
	if(	flag->outputResultFlag ||
		flag->outputBlankFlag ||
		flag->outputDispFlag ||
		flag->outputEuclDispFlag ||
		flag->outputPosFlag)
		positionFieldNeeded=true;
	nifti_image *positionFieldImage=NULL;
	if(positionFieldNeeded==true){
		positionFieldImage = nifti_copy_nim_info(targetImage);
		positionFieldImage->cal_min=0;
		positionFieldImage->cal_max=0;
        positionFieldImage->scl_slope = 1.0f;
        positionFieldImage->scl_inter = 0.0f;
		positionFieldImage->dim[0]=positionFieldImage->ndim=5;
		positionFieldImage->dim[1]=positionFieldImage->nx=targetImage->nx;
		positionFieldImage->dim[2]=positionFieldImage->ny=targetImage->ny;
		positionFieldImage->dim[3]=positionFieldImage->nz=targetImage->nz;
		positionFieldImage->dim[4]=positionFieldImage->nt=1;positionFieldImage->pixdim[4]=positionFieldImage->dt=1.0;
		if(targetImage->nz>1)
			positionFieldImage->dim[5]=positionFieldImage->nu=3;
		else positionFieldImage->dim[5]=positionFieldImage->nu=2;
		positionFieldImage->pixdim[5]=positionFieldImage->du=1.0;
		positionFieldImage->dim[6]=positionFieldImage->nv=1;positionFieldImage->pixdim[6]=positionFieldImage->dv=1.0;
		positionFieldImage->dim[7]=positionFieldImage->nw=1;positionFieldImage->pixdim[7]=positionFieldImage->dw=1.0;
		positionFieldImage->nvox=positionFieldImage->nx*positionFieldImage->ny*positionFieldImage->nz*positionFieldImage->nt*positionFieldImage->nu;
		if(sizeof(PrecisionTYPE)==4) positionFieldImage->datatype = NIFTI_TYPE_FLOAT32;
		else positionFieldImage->datatype = NIFTI_TYPE_FLOAT64;
		positionFieldImage->nbyper = sizeof(PrecisionTYPE);
		positionFieldImage->data = (void *)calloc(positionFieldImage->nvox, positionFieldImage->nbyper);
	}

    bool controlPointPositionNeeded=false;
    nifti_image *velocityFieldImage = NULL;
    if(flag->inputVelocityFieldFlag){
        velocityFieldImage = nifti_image_read(param->inputVelocityFieldName,true);
        if(velocityFieldImage == NULL){
            fprintf(stderr,"** ERROR Error when reading the velocity field image: %s\n",param->inputVelocityFieldName);
            return 1;
        }
        if(positionFieldNeeded==true || flag->outputCPPFlag)
            controlPointPositionNeeded=true;

        if(flag->outputJacobianFlag){
            nifti_image *jacobianImage = nifti_copy_nim_info(targetImage);
            jacobianImage->scl_slope = 1.0f;
            jacobianImage->scl_inter = 0.0f;
            jacobianImage->cal_min=0;
            jacobianImage->cal_max=0;
            jacobianImage->datatype = NIFTI_TYPE_FLOAT32;
            jacobianImage->nbyper = sizeof(float);
            jacobianImage->data = (void *)calloc(jacobianImage->nvox, jacobianImage->nbyper);
            nifti_set_filenames(jacobianImage, param->outputJacobianName, 0, 0);
            reg_bspline_GetJacobianMapFromVelocityField(velocityFieldImage,
                                                        jacobianImage);
            nifti_image_write(jacobianImage);
            nifti_image_free(jacobianImage);
            printf("Jacobian map image has been saved from the velocity field: %s\n", param->outputJacobianName);
        }
    }

	if(controlPointPositionNeeded==true || flag->inputCPPFlag){
        nifti_image *controlPointImage = NULL;
        if(flag->inputCPPFlag){
		    /* Read the CPP image */
		    controlPointImage = nifti_image_read(param->inputCPPName,true);
		    if(controlPointImage == NULL){
			    fprintf(stderr,"** ERROR Error when reading the cpp image: %s\n",param->inputCPPName);
			    return 1;
		    }
        }
        else{
            /* generate the CPP image from the velocity field */
            controlPointImage = nifti_copy_nim_info(velocityFieldImage);
            controlPointImage->scl_slope = 1.0f;
            controlPointImage->scl_inter = 0.0f;
            controlPointImage->cal_min=0;
            controlPointImage->cal_max=0;
            controlPointImage->datatype = NIFTI_TYPE_FLOAT32;
            controlPointImage->nbyper = sizeof(float);
            controlPointImage->data = (void *)calloc(controlPointImage->nvox, controlPointImage->nbyper);
            reg_spline_scaling_squaring(velocityFieldImage,
                                        controlPointImage);
            if(flag->outputCPPFlag){
                nifti_set_filenames(controlPointImage, param->outputCPPName, 0, 0);
                nifti_image_write(controlPointImage);
            }
        }
		/* apply the cubic spline interpolation to generate the position field */
		if(positionFieldNeeded==true){
			reg_bspline<PrecisionTYPE>(	controlPointImage,
							            targetImage,
							            positionFieldImage,
                                        NULL,
							            0); // new df
		}
		/* Generate the jacobian map */
		if(flag->outputJacobianFlag && !flag->inputVelocityFieldFlag){

			nifti_image *jacobianImage = nifti_copy_nim_info(targetImage);
            jacobianImage->cal_min=0;
            jacobianImage->cal_max=0;
            jacobianImage->scl_slope = 1.0f;
            jacobianImage->scl_inter = 0.0f;
			jacobianImage->datatype = NIFTI_TYPE_FLOAT32;
			jacobianImage->nbyper = sizeof(float);
			jacobianImage->data = (void *)calloc(jacobianImage->nvox, jacobianImage->nbyper);
			nifti_set_filenames(jacobianImage, param->outputJacobianName, 0, 0);
			reg_bspline_GetJacobianMap(	controlPointImage,
									   jacobianImage);
			nifti_image_write(jacobianImage);
			nifti_image_free(jacobianImage);
			printf("Jacobian map image has been saved: %s\n", param->outputJacobianName);
		}
		/* Generate the jacobian matrix image */
		if(flag->outputJacobianMatrixFlag && !flag->inputVelocityFieldFlag){
			nifti_image *jacobianImage = nifti_copy_nim_info(targetImage);
			jacobianImage->cal_min=0;
			jacobianImage->cal_max=0;
            jacobianImage->scl_slope = 1.0f;
            jacobianImage->scl_inter = 0.0f;
			jacobianImage->dim[0]=jacobianImage->ndim=5;
			jacobianImage->dim[1]=jacobianImage->nx=targetImage->nx;
			jacobianImage->dim[2]=jacobianImage->ny=targetImage->ny;
			jacobianImage->dim[3]=jacobianImage->nz=targetImage->nz;
			jacobianImage->dim[4]=jacobianImage->nt=1;jacobianImage->pixdim[4]=jacobianImage->dt=1.0;
			jacobianImage->dim[5]=jacobianImage->nu=controlPointImage->nu*controlPointImage->nu;
            jacobianImage->pixdim[5]=jacobianImage->du=1.0;
			jacobianImage->dim[6]=jacobianImage->nv=1;jacobianImage->pixdim[6]=jacobianImage->dv=1.0;
			jacobianImage->dim[7]=jacobianImage->nw=1;jacobianImage->pixdim[7]=jacobianImage->dw=1.0;
			jacobianImage->nvox=jacobianImage->nx*jacobianImage->ny*jacobianImage->nz*jacobianImage->nt*jacobianImage->nu;
			jacobianImage->datatype = NIFTI_TYPE_FLOAT32;
			jacobianImage->nbyper = sizeof(float);
			jacobianImage->data = (void *)calloc(jacobianImage->nvox, jacobianImage->nbyper);
			nifti_set_filenames(jacobianImage, param->outputJacobianMatrixName, 0, 0);
	
			reg_bspline_GetJacobianMatrix(	controlPointImage,
							jacobianImage);
	
			nifti_image_write(jacobianImage);
			nifti_image_free(jacobianImage);
			printf("Jacobian matrix image has been saved: %s\n", param->outputJacobianMatrixName);
		}
		nifti_image_free(controlPointImage);
	}
	else{
		mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));
		affineTransformation->m[0][0]=1.0f;
		affineTransformation->m[1][1]=1.0f;
		affineTransformation->m[2][2]=1.0f;
		affineTransformation->m[3][3]=1.0f;
		if(flag->affineMatrixFlag){
			reg_tool_ReadAffineFile(	affineTransformation,
							targetImage,
							sourceImage,
							param->affineMatrixName,
							flag->affineFlirtFlag);
		}
		if(positionFieldNeeded==true){
			reg_affine_positionField(	affineTransformation,
							            targetImage,
							            positionFieldImage);
		}
		free(affineTransformation);
	}
    nifti_image_free(velocityFieldImage);

	/* Resample the source image */
    if(flag->outputResultFlag){
        int inter=3;
        if(flag->TRIInterpolationFlag){
            inter=1;
        }
        else if(flag->NNInterpolationFlag){
            inter=0;
        }
        nifti_image *resultImage = nifti_copy_nim_info(targetImage);
        resultImage->cal_min=sourceImage->cal_min;
        resultImage->cal_max=sourceImage->cal_max;
        resultImage->scl_slope=sourceImage->scl_slope;
        resultImage->scl_inter=sourceImage->scl_inter;
        resultImage->datatype = sourceImage->datatype;
        resultImage->nbyper = sourceImage->nbyper;
        resultImage->data = (void *)calloc(resultImage->nvox, resultImage->nbyper);
        reg_resampleSourceImage<double>(targetImage,
                                        sourceImage,
                                        resultImage,
                                        positionFieldImage,
                                        NULL,
                                        inter,
                                        param->sourceBGValue);
        nifti_set_filenames(resultImage, param->outputResultName, 0, 0);
        nifti_image_write(resultImage);
        nifti_image_free(resultImage);
        printf("Resampled image has been saved: %s\n", param->outputResultName);
    }

	/* Resample a blank grid image */
	if(flag->outputBlankFlag){
		nifti_image *gridImage = nifti_copy_nim_info(sourceImage);
		gridImage->cal_min=0;
		gridImage->cal_max=255;
		gridImage->datatype = NIFTI_TYPE_UINT8;
		gridImage->nbyper = sizeof(unsigned char);
		gridImage->data = (void *)calloc(gridImage->nvox, gridImage->nbyper);
		unsigned char *gridImageValuePtr = static_cast<unsigned char *>(gridImage->data);
		for(int z=0; z<gridImage->nz;z++){
			for(int y=0; y<gridImage->ny;y++){
				for(int x=0; x<gridImage->nx;x++){
					if(targetImage->nz>1){
						if( x/10==(float)x/10.0 || y/10==(float)y/10.0 || z/10==(float)z/10.0)
							*gridImageValuePtr = 255;
					}
					else{
						if( x/10==(float)x/10.0 || x==targetImage->nx-1 || y/10==(float)y/10.0 || y==targetImage->ny-1)
							*gridImageValuePtr = 255;
					}
					gridImageValuePtr++;
				}
			}
		}
		
		nifti_image *resultImage = nifti_copy_nim_info(targetImage);
		resultImage->cal_min=gridImage->cal_min;
		resultImage->cal_max=gridImage->cal_max;
        resultImage->scl_slope = 1.0f;
        resultImage->scl_inter = 0.0f;
		resultImage->datatype = gridImage->datatype;
		resultImage->nbyper = gridImage->nbyper;
		resultImage->data = (void *)calloc(resultImage->nvox, resultImage->nbyper);
		reg_resampleSourceImage<double>(	targetImage,
						                    gridImage,
						                    resultImage,
						                    positionFieldImage,
                                            NULL,
						                    1,
						                    0);
		nifti_image_free(gridImage);
		nifti_set_filenames(resultImage, param->outputBlankName, 0, 0);
		nifti_image_write(resultImage);
		nifti_image_free(resultImage);
		printf("Resampled grid image has been saved: %s\n", param->outputBlankName);
	}

	/* Output the position field */
	if(flag->outputPosFlag){
		nifti_set_filenames(positionFieldImage, param->outputPosName, 0, 0);
		nifti_image_write(positionFieldImage);
		printf("Position field image has been saved: %s\n", param->outputPosName);
	}

	/* Output the displacement field */
	if(flag->outputDispFlag || flag->outputEuclDispFlag ){
		nifti_image *displacementFieldImage = nifti_copy_nim_info(positionFieldImage);
        displacementFieldImage->scl_slope = 1.0f;
        displacementFieldImage->scl_inter = 0.0f;
		displacementFieldImage->data = (void *)calloc(displacementFieldImage->nvox, displacementFieldImage->nbyper);
		memcpy(displacementFieldImage->data, positionFieldImage->data, displacementFieldImage->nvox*displacementFieldImage->nbyper);
		
		nifti_image *euclidianDispImage=NULL;
		PrecisionTYPE *euclidianDispPtr=NULL;
		if(flag->outputEuclDispFlag){
			euclidianDispImage = nifti_copy_nim_info(targetImage);
			if(sizeof(PrecisionTYPE)==4) euclidianDispImage->datatype = NIFTI_TYPE_FLOAT32;
			else euclidianDispImage->datatype = NIFTI_TYPE_FLOAT64;
			euclidianDispImage->nbyper=sizeof(PrecisionTYPE);
			euclidianDispImage->scl_slope = 1.0f;
			euclidianDispImage->scl_inter = 0.0f;
			euclidianDispImage->data = (void *)calloc(euclidianDispImage->nvox, euclidianDispImage->nbyper);
			nifti_set_filenames(euclidianDispImage, param->outputEuclDispName, 0, 0);
			euclidianDispPtr=static_cast<PrecisionTYPE *>(euclidianDispImage->data);
		}
		
		mat44 *target_voxel_2_real=NULL;
		if(targetImage->sform_code)
			target_voxel_2_real=&targetImage->sto_xyz;
		else target_voxel_2_real=&targetImage->qto_xyz;
		
		if(targetImage->nz>1){
			PrecisionTYPE *fullDefPtrX=static_cast<PrecisionTYPE *>(displacementFieldImage->data);
			PrecisionTYPE *fullDefPtrY=&fullDefPtrX[targetImage->nvox];
			PrecisionTYPE *fullDefPtrZ=&fullDefPtrY[targetImage->nvox];
			PrecisionTYPE position[3];
			for(int z=0; z<displacementFieldImage->nz; z++){
				for(int y=0; y<displacementFieldImage->ny; y++){
					for(int x=0; x<displacementFieldImage->nx; x++){
						position[0]=x*target_voxel_2_real->m[0][0] + y*target_voxel_2_real->m[0][1]
							+ z*target_voxel_2_real->m[0][2] + target_voxel_2_real->m[0][3];
						position[1]=x*target_voxel_2_real->m[1][0] + y*target_voxel_2_real->m[1][1]
							+ z*target_voxel_2_real->m[1][2] + target_voxel_2_real->m[1][3];
						position[2]=x*target_voxel_2_real->m[2][0] + y*target_voxel_2_real->m[2][1]
							+ z*target_voxel_2_real->m[2][2] + target_voxel_2_real->m[2][3];
						*fullDefPtrX -= position[0];
						*fullDefPtrY -= position[1];
						*fullDefPtrZ -= position[2];
						
						if(flag->outputEuclDispFlag){
							*euclidianDispPtr++ = sqrt((*fullDefPtrX * (*fullDefPtrX))+
								(*fullDefPtrY * (*fullDefPtrY))+(*fullDefPtrZ * (*fullDefPtrZ)));
						}
						
						fullDefPtrX++;
						fullDefPtrY++;
						fullDefPtrZ++;
					}
				}
			}
		}
		else{
			PrecisionTYPE *fullDefPtrX=static_cast<PrecisionTYPE *>(displacementFieldImage->data);
			PrecisionTYPE *fullDefPtrY=&fullDefPtrX[targetImage->nvox];
			PrecisionTYPE position[3];
			for(int y=0; y<displacementFieldImage->ny; y++){
				for(int x=0; x<displacementFieldImage->nx; x++){
					position[0]=x*target_voxel_2_real->m[0][0] + y*target_voxel_2_real->m[0][1] + target_voxel_2_real->m[0][3];
					position[1]=x*target_voxel_2_real->m[1][0] + y*target_voxel_2_real->m[1][1] + target_voxel_2_real->m[1][3];
					*fullDefPtrX -= position[0];
					*fullDefPtrY -= position[1];
					
					if(flag->outputEuclDispFlag){
						*euclidianDispPtr++ = sqrt((*fullDefPtrX * (*fullDefPtrX))+
												   (*fullDefPtrY * (*fullDefPtrY)));
					}
					
					fullDefPtrX++;
					fullDefPtrY++;
				}
			}
			
		}
		if(flag->outputEuclDispFlag){
			nifti_image_write(euclidianDispImage);
			printf("Euclidian displacement image has been saved: %s\n", param->outputEuclDispName);
			nifti_image_free(euclidianDispImage);
		}
		if(flag->outputDispFlag){
			nifti_set_filenames(displacementFieldImage, param->outputDispName, 0, 0);
			nifti_image_write(displacementFieldImage);
			printf("Displacement field image has been saved: %s\n", param->outputDispName);
		}
		nifti_image_free(displacementFieldImage);
	}


	if(positionFieldNeeded==true) nifti_image_free(positionFieldImage);
	nifti_image_free(targetImage);
	nifti_image_free(sourceImage);
	
	free(flag);
	free(param);
	
	return 0;
}

#endif
