
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
	printf("\t-affDirect\t\tDirectly optimize 12 DoF affine [default is rigid initially then affine]\n");
    printf("\t-tmask <filename>\tFilename of a mask image in the target space\n");
	printf("\t-maxit <int>\t\tNumber of iteration per level [5]\n");
	printf("\t-smooT <float>\t\tSmooth the target image using the specified sigma (mm) [0]\n");
	printf("\t-smooS <float>\t\tSmooth the source image using the specified sigma (mm) [0]\n");
	printf("\t-ln <int>\t\tNumber of level to perform [3]\n");
	printf("\t-lp <int>\t\tOnly perform the first levels [ln]\n");
	
	printf("\t-noOrigin\t\t\tDo not use the nifti header origins to initialise the translation\n");
	
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
	float *tpm;
	
	flag->affineFlag = 1;
	flag->rigidFlag = 1;
	param->block_percent_to_use = 50;
	param->inlier_lts = 50;
    flag->alignCenterFlag = 0;

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
			param->tpmImageName=argv[++i];
			flag->tpmImageFlag=1;
		}
		else if(strcmp(argv[i], "-smooT") == 0){
			param->targetSigmaValue=(float)(atof(argv[++i]));
			flag->targetSigmaFlag=1;
		}
		else if(strcmp(argv[i], "-smooS") == 0){
			param->sourceSigmaValue=(float)(atof(argv[++i]));
			flag->sourceSigmaFlag=1;
		}
		else if(strcmp(argv[i], "-affDirect") == 0){
			flag->rigidFlag=0;
		}
		else if(strcmp(argv[i], "-noOrigin") == 0){
			flag->alignCenterFlag=1;
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
	
	if(!flag->targetImageFlag || !flag->sourceImageFlag || !flag->tpmImageFlag){
		fprintf(stderr,"Err:\tThe target, source, and tpm images have to be defined.\n");
		PetitUsage(argv[0]);
		return 1;
	}
	
	param->levelNumber=3;
	param->maxIteration=5;
//	param->maxIteration=1;
	param->level2Perform=param->levelNumber;
	param->level2Perform=param->level2Perform<param->levelNumber?param->level2Perform:param->levelNumber;
	
	/* Read the target and source images */
    nifti_image *sourceImage = nifti_image_read(param->sourceImageName,true);
    if(sourceImage->data == NULL){
    	fprintf(stderr, "** ERROR Error when reading the source image: %s\n", param->sourceImageName);
    	return NULL;
    }

    nifti_image *tpmImage = read_nifti_float(param->tpmImageName,&tpm);
    if(tpmImage->data == NULL){
    	fprintf(stderr, "** ERROR Error when reading the source image: %s\n", param->tpmImageName);
    	return NULL;
    }

    /* check number of classes */
    if(tpmImage->nt != 6){
    	fprintf(stderr, "** ERROR TPM image has %d instead of 6 classes: %s\n", tpmImage->nt, param->tpmImageName);
    	return NULL;
    }

    /* ensure that nvox is that of a 3D image */
    tpmImage->nvox  /= tpmImage->nt;
    tpmImage->dim[0] = tpmImage->ndim = 4;

    nifti_image *targetImage = nifti_image_read(param->targetImageName,true);		
    if(targetImage->data == NULL){
			fprintf(stderr, "** ERROR Error when reading the target image: %s\n", param->targetImageName);
			return NULL;
    }

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
    resultImage->data = (void *)calloc(resultImage->nvox, resultImage->nbyper);

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
    double *src          = (double *)malloc(sizeof(double)*sourceImage->nvox);
    
    ptr_to_double_4D(sourceImage, src, 0);

    double voxelsize[3];
    int dims[3];
    
    for (int j=0; j<3; j++) {
        voxelsize[j] = sourceImage->pixdim[j+1];
        dims[j] = sourceImage->dim[j+1];
    }
    
    double slope = 1.0;

    Bayes(src, label, priors, probs, voxelsize, dims, 1);
    
    int cleanup_strength = 1;
    double scale = 3.0/(sourceImage->dx + sourceImage->dy + sourceImage->dz);
    
    cleanup(probs, label, dims, cleanup_strength, scale);

    for (int i = 0; i < sourceImage->nvox; i++)
      src[i] = (double)priors[i];

    if(!write_nifti_double("test0.nii", src, NIFTI_TYPE_FLOAT32, slope, 
            dims, voxelsize, sourceImage))
      exit(EXIT_FAILURE);

    free(tpmImage);
    free(sourceImage);
    free(label);
    free(src);
    free(probs);
    free(priors);

	return 0;
}

