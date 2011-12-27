
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
    nifti_image *resultImage;

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
	param->maxIteration=5;
	param->level2Perform=param->levelNumber;
	param->level2Perform=param->level2Perform<param->levelNumber?param->level2Perform:param->levelNumber;
	
	resultImage = affineRegistration(param, flag);
	
	if (resultImage == NULL) {
		fprintf(stderr,"Error\n");
		return 1;
	}
	
	nifti_set_filenames(resultImage, param->outputResultName, 0, 0);
	nifti_image_write(resultImage);

	return 0;
}

