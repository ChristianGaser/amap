#define RIGID 0
#define AFFINE 1

#include "_reg_resampling.h"
#include "_reg_affineTransformation.h"
#include "_reg_blockMatching.h"
#include "_reg_tools.h"

#ifdef _WINDOWS
    #include <time.h>
#endif

#define PrecisionTYPE float

typedef struct{
	char *targetImageName;
	char *sourceImageName;
	char *outputResultName;

	int maxIteration;

	int backgroundIndex[3];
	PrecisionTYPE sourceBGValue;

	float targetSigmaValue;
	float sourceSigmaValue;
	int levelNumber;
	int level2Perform;

	int block_percent_to_use;
	int inlier_lts;
}PARAM;

typedef struct{
	bool targetImageFlag;
	bool sourceImageFlag;
	bool levelNumberFlag;
	bool level2PerformFlag;
	bool outputResultFlag;

	bool maxIterationFlag;

	bool backgroundIndexFlag;
	
	bool alignCenterFlag;

	bool rigidFlag;
	bool affineFlag;

	bool targetSigmaFlag;
	bool sourceSigmaFlag;
	bool pyramidFlag;
	bool useGPUFlag;
    bool twoDimRegistration;
}FLAG;

extern nifti_image *affineRegistration(PARAM *param, FLAG *flag);
