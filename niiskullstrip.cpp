
#include "niiskullstrip.h"

void PetitUsage(char *exec)
{
	fprintf(stderr,"Aladin - Seb.\n");
	fprintf(stderr,"Usage:\t%s -target <targetImageName> -source <sourceImageName> [OPTIONS].\n",exec);
	fprintf(stderr,"\tSee the help for more details (-h).\n");
	return;
}
void Usage(char *exec)
{
	printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
	printf("Block Matching algorithm for global registration.\n");
	printf("Based on Ourselin et al., \"Reconstructing a 3D structure from serial histological sections\",\n");
    printf("Image and Vision Computing, 2001\n");
	printf("This code has been written by Marc Modat (m.modat@ucl.ac.uk) and Pankaj Daga,\n");
	printf("for any comment, please contact them.\n");
	printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
	printf("Usage:\t%s -target <filename> -source <filename> [OPTIONS].\n",exec);
	printf("\t-target <filename>\tFilename of the target image (mandatory)\n");
	printf("\t-source <filename>\tFilename of the source image (mandatory)\n");
	printf("* * OPTIONS * *\n");
	printf("\t-rigOnly\t\tTo perform a rigid registration only (rigid+affine by default)\n");
	printf("\t-affDirect\t\tDirectly optimize 12 DoF affine [default is rigid initially then affine]\n");
	printf("\t-inaff <filename>\tFilename which contains an input affine transformation (Affine*Target=Source) [none]\n");
	printf("\t-affFlirt <filename>\tFilename which contains an input affine transformation from Flirt [none]\n");
    printf("\t-tmask <filename>\tFilename of a mask image in the target space\n");
	printf("\t-result <filename>\tFilename of the resampled image [outputResult.nii]\n");
	printf("\t-maxit <int>\t\tNumber of iteration per level [5]\n");
	printf("\t-smooT <float>\t\tSmooth the target image using the specified sigma (mm) [0]\n");
	printf("\t-smooS <float>\t\tSmooth the source image using the specified sigma (mm) [0]\n");
	printf("\t-ln <int>\t\tNumber of level to perform [3]\n");
	printf("\t-lp <int>\t\tOnly perform the first levels [ln]\n");
	
	printf("\t-nac\t\t\tUse the nifti header origins to initialise the translation\n");
	
	printf("\t-bgi <int> <int> <int>\tForce the background value during\n\t\t\t\tresampling to have the same value as this voxel in the source image [none]\n");

	printf("\t-%%v <int>\t\tPercentage of block to use [50]\n");
	printf("\t-%%i <int>\t\tPercentage of inlier for the LTS [50]\n");
	printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
	return;
}

int main(int argc, char **argv)
{
	
	PARAM *param = (PARAM *)calloc(1,sizeof(PARAM));
	FLAG *flag = (FLAG *)calloc(1,sizeof(FLAG));

	flag->affineFlag=1;
	flag->rigidFlag=1;
	param->block_percent_to_use=50;
	param->inlier_lts=50;
    flag->alignCenterFlag=1;

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
		    for(int j=0; j<6; j++) 
		        param->tpmImageName[j]=argv[++i];
	        flag->tpmImageFlag=1;
		}
		else if(strcmp(argv[i], "-result") == 0){
			param->outputResultName=argv[++i];
			flag->outputResultFlag=1;
		}
		else if(strcmp(argv[i], "-smooT") == 0){
			param->targetSigmaValue=(float)(atof(argv[++i]));
			flag->targetSigmaFlag=1;
		}
		else if(strcmp(argv[i], "-smooS") == 0){
			param->sourceSigmaValue=(float)(atof(argv[++i]));
			flag->sourceSigmaFlag=1;
		}
		else if(strcmp(argv[i], "-rigOnly") == 0){
			flag->affineFlag=0;
		}
		else if(strcmp(argv[i], "-affDirect") == 0){
			flag->rigidFlag=0;
		}
		else if(strcmp(argv[i], "-nac") == 0){
			flag->alignCenterFlag=0;
		}
		else if(strcmp(argv[i], "-bgi") == 0){
			param->backgroundIndex[0]=atoi(argv[++i]);
			param->backgroundIndex[1]=atoi(argv[++i]);
			param->backgroundIndex[2]=atoi(argv[++i]);
			flag->backgroundIndexFlag=1;
		}
		else if(strcmp(argv[i], "-%v") == 0){
			param->block_percent_to_use=atoi(argv[++i]);
		}
		else if(strcmp(argv[i], "-%i") == 0){
			param->inlier_lts=atoi(argv[++i]);
		}
		else{
			fprintf(stderr,"Err:\tParameter %s unknown.\n",argv[i]);
			PetitUsage(argv[0]);
			return 1;
		}
	}
	
	if(!flag->targetImageFlag || !flag->sourceImageFlag){
		fprintf(stderr,"Err:\tThe target and the source image have to be defined.\n");
		PetitUsage(argv[0]);
		return 1;
	}
	
	param->levelNumber=3;
//	param->maxIteration=5;
	param->maxIteration=1;
	param->level2Perform=param->levelNumber;
	param->level2Perform=param->level2Perform<param->levelNumber?param->level2Perform:param->levelNumber;
	
	/* Read the target and source images */
	nifti_image *targetHeader = nifti_image_read(param->targetImageName,false);
	if(targetHeader == NULL){
		fprintf(stderr,"** ERROR Error when reading the target image: %s\n",param->targetImageName);
		return NULL;
	}

    nifti_image *sourceImage = nifti_image_read(param->sourceImageName,true);
    if(sourceImage->data == NULL){
    	fprintf(stderr, "** ERROR Error when reading the source image: %s\n", param->sourceImageName);
    	return NULL;
    }

	nifti_image **tpmImage = (nifti_image **)calloc(6,sizeof(nifti_image));
    for (int j=0; j<6; j++) {
        tpmImage[j] = nifti_image_read(param->tpmImageName[j],true);
        if(tpmImage[j]->data == NULL){
    	    fprintf(stderr, "** ERROR Error when reading the source image: %s\n", param->tpmImageName[j]);
    	    return NULL;
        }
    }
    
    nifti_image *targetImage = nifti_image_read(param->targetImageName,true);		
    if(targetImage->data == NULL){
			fprintf(stderr, "** ERROR Error when reading the target image: %s\n", param->targetImageName);
			return NULL;
    }

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

    reg_affine_positionField(affineRegistration(param, flag),
							targetHeader,
							positionFieldImage);

    /* allocate the result image */
    unsigned char *priors = (unsigned char *)malloc(sizeof(unsigned char)*sourceImage->nvox*6);
//	nifti_image **resultImage = (nifti_image **)calloc(6,sizeof(nifti_image));
    for (int j=0; j<6; j++) {
    
        nifti_image *resultImage = nifti_copy_nim_info(targetImage);
        resultImage = nifti_copy_nim_info(targetHeader);
        resultImage->cal_min=tpmImage[0]->cal_min;
        resultImage->cal_max=tpmImage[0]->cal_max;
        resultImage->scl_slope=tpmImage[0]->scl_slope;
        resultImage->scl_inter=tpmImage[0]->scl_inter;
        resultImage->datatype = DT_UINT8;
        resultImage->nbyper = 1;
        resultImage->data = (void *)calloc(resultImage->nvox, resultImage->nbyper);
        reg_resampleSourceImage<double>(targetHeader,
							tpmImage[j],
							resultImage,
							positionFieldImage,
                            NULL,
							3,
							param->sourceBGValue);
							
        for(int i=0; i<sourceImage->nvox;i++) 
            priors[i+j*sourceImage->nvox] = resultImage->data[i];

    }

    unsigned char *label = (unsigned char *)malloc(sizeof(unsigned char)*sourceImage->nvox);
    double *src   = (double *)malloc(sizeof(double)*sourceImage->nvox);
    
//    Bayes( sourceImage, label, priors, 100, separations, dims, 10);

//    nifti_set_filenames(resultImage[0], "test2.nii", 0, 0);
//	nifti_image_write(resultImage[0]);

	return 0;
}

