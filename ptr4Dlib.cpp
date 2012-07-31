/*
 * Christian Gaser
 * $Id$ 
 *
 */

#include "niiskullstrip.h"

template <class DTYPE>
void float_to_ptr_4D2(nifti_image *image, float *vol4d, int index_nt)
{
    DTYPE *imagePtr = static_cast<DTYPE *>(image->data);
	DTYPE currentMin=0;
	DTYPE currentMax=0;
	switch(image->datatype){
		case NIFTI_TYPE_UINT8:
			currentMin=(DTYPE)std::numeric_limits<unsigned char>::max();
			currentMax=0;
			break;
		case NIFTI_TYPE_INT8:
			currentMin=(DTYPE)std::numeric_limits<char>::max();
			currentMax=(DTYPE)-std::numeric_limits<char>::max();
			break;
		case NIFTI_TYPE_UINT16:
			currentMin=(DTYPE)std::numeric_limits<unsigned short>::max();
			currentMax=0;
			break;
		case NIFTI_TYPE_INT16:
			currentMin=(DTYPE)std::numeric_limits<short>::max();
			currentMax=-(DTYPE)std::numeric_limits<short>::max();
			break;
		case NIFTI_TYPE_UINT32:
			currentMin=(DTYPE)std::numeric_limits<unsigned int>::max();
			currentMax=0;
			break;
		case NIFTI_TYPE_INT32:
			currentMin=(DTYPE)std::numeric_limits<int>::max();
			currentMax=-(DTYPE)std::numeric_limits<int>::max();
			break;
		case NIFTI_TYPE_FLOAT32:
			currentMin=(DTYPE)std::numeric_limits<float>::max();
			currentMax=-(DTYPE)std::numeric_limits<float>::max();
			break;
		case NIFTI_TYPE_FLOAT64:
			currentMin=(DTYPE)std::numeric_limits<double>::max();
			currentMax=-(DTYPE)std::numeric_limits<double>::max();
			break;
	}

    if(image->scl_slope==0) image->scl_slope=1.0f;

	imagePtr = static_cast<DTYPE *>(image->data);
	for(unsigned int index=0; index<image->nvox; index++){
		double value = ((double)vol4d[index+(int)index_nt*(int)image->nvox] - image->scl_inter) / image->scl_slope;
		
		*imagePtr++=(DTYPE)value;
	}

}

void 
float_to_ptr_4D(nifti_image *image, float *vol4d, int index_nt)
{
	switch(image->datatype){
		case NIFTI_TYPE_UINT8:
			float_to_ptr_4D2<unsigned char>(image, vol4d, index_nt);
			break;
		case NIFTI_TYPE_INT8:
			float_to_ptr_4D2<char>(image, vol4d, index_nt);
			break;
		case NIFTI_TYPE_UINT16:
			float_to_ptr_4D2<unsigned short>(image, vol4d, index_nt);
			break;
		case NIFTI_TYPE_INT16:
			float_to_ptr_4D2<short>(image, vol4d, index_nt);
			break;
		case NIFTI_TYPE_UINT32:
			float_to_ptr_4D2<unsigned int>(image, vol4d, index_nt);
			break;
		case NIFTI_TYPE_INT32:
			float_to_ptr_4D2<int>(image, vol4d, index_nt);
			break;
		case NIFTI_TYPE_FLOAT32:
			float_to_ptr_4D2<float>(image, vol4d, index_nt);
			break;
		case NIFTI_TYPE_FLOAT64:
			float_to_ptr_4D2<double>(image, vol4d, index_nt);
			break;
		default:
			printf("err\tfloat_to_ptr_4D2\tThe image data type is not supported\n");
			return;
	}
}

template <class DTYPE>
void double_to_ptr_4D2(nifti_image *image, double *vol4d, int index_nt)
{
    DTYPE *imagePtr = static_cast<DTYPE *>(image->data);
	DTYPE currentMin=0;
	DTYPE currentMax=0;
	switch(image->datatype){
		case NIFTI_TYPE_UINT8:
			currentMin=(DTYPE)std::numeric_limits<unsigned char>::max();
			currentMax=0;
			break;
		case NIFTI_TYPE_INT8:
			currentMin=(DTYPE)std::numeric_limits<char>::max();
			currentMax=(DTYPE)-std::numeric_limits<char>::max();
			break;
		case NIFTI_TYPE_UINT16:
			currentMin=(DTYPE)std::numeric_limits<unsigned short>::max();
			currentMax=0;
			break;
		case NIFTI_TYPE_INT16:
			currentMin=(DTYPE)std::numeric_limits<short>::max();
			currentMax=-(DTYPE)std::numeric_limits<short>::max();
			break;
		case NIFTI_TYPE_UINT32:
			currentMin=(DTYPE)std::numeric_limits<unsigned int>::max();
			currentMax=0;
			break;
		case NIFTI_TYPE_INT32:
			currentMin=(DTYPE)std::numeric_limits<int>::max();
			currentMax=-(DTYPE)std::numeric_limits<int>::max();
			break;
		case NIFTI_TYPE_FLOAT32:
			currentMin=(DTYPE)std::numeric_limits<double>::max();
			currentMax=-(DTYPE)std::numeric_limits<double>::max();
			break;
		case NIFTI_TYPE_FLOAT64:
			currentMin=(DTYPE)std::numeric_limits<double>::max();
			currentMax=-(DTYPE)std::numeric_limits<double>::max();
			break;
	}

    if(image->scl_slope==0) image->scl_slope=1.0f;

	imagePtr = static_cast<DTYPE *>(image->data);
	for(unsigned int index=0; index<image->nvox; index++){
		double value = ((double)vol4d[index+(int)index_nt*(int)image->nvox] - image->scl_inter) / image->scl_slope;
		
		*imagePtr++=(DTYPE)value;
	}

}

void 
double_to_ptr_4D(nifti_image *image, double *vol4d, int index_nt)
{
	switch(image->datatype){
		case NIFTI_TYPE_UINT8:
			double_to_ptr_4D2<unsigned char>(image, vol4d, index_nt);
			break;
		case NIFTI_TYPE_INT8:
			double_to_ptr_4D2<char>(image, vol4d, index_nt);
			break;
		case NIFTI_TYPE_UINT16:
			double_to_ptr_4D2<unsigned short>(image, vol4d, index_nt);
			break;
		case NIFTI_TYPE_INT16:
			double_to_ptr_4D2<short>(image, vol4d, index_nt);
			break;
		case NIFTI_TYPE_UINT32:
			double_to_ptr_4D2<unsigned int>(image, vol4d, index_nt);
			break;
		case NIFTI_TYPE_INT32:
			double_to_ptr_4D2<int>(image, vol4d, index_nt);
			break;
		case NIFTI_TYPE_FLOAT32:
			double_to_ptr_4D2<float>(image, vol4d, index_nt);
			break;
		case NIFTI_TYPE_FLOAT64:
			double_to_ptr_4D2<double>(image, vol4d, index_nt);
			break;
		default:
			printf("err\tdouble_to_ptr_4D2\tThe image data type is not supported\n");
			return;
	}
}

template <class DTYPE>
void ptr_to_uint8_4D2(nifti_image *image, unsigned char *vol4d, int index_nt)
{
    DTYPE *imagePtr = static_cast<DTYPE *>(image->data);
	DTYPE currentMin=0;
	DTYPE currentMax=0;
	switch(image->datatype){
		case NIFTI_TYPE_UINT8:
			currentMin=(DTYPE)std::numeric_limits<unsigned char>::max();
			currentMax=0;
			break;
		case NIFTI_TYPE_INT8:
			currentMin=(DTYPE)std::numeric_limits<char>::max();
			currentMax=(DTYPE)-std::numeric_limits<char>::max();
			break;
		case NIFTI_TYPE_UINT16:
			currentMin=(DTYPE)std::numeric_limits<unsigned short>::max();
			currentMax=0;
			break;
		case NIFTI_TYPE_INT16:
			currentMin=(DTYPE)std::numeric_limits<short>::max();
			currentMax=-(DTYPE)std::numeric_limits<short>::max();
			break;
		case NIFTI_TYPE_UINT32:
			currentMin=(DTYPE)std::numeric_limits<unsigned int>::max();
			currentMax=0;
			break;
		case NIFTI_TYPE_INT32:
			currentMin=(DTYPE)std::numeric_limits<int>::max();
			currentMax=-(DTYPE)std::numeric_limits<int>::max();
			break;
		case NIFTI_TYPE_FLOAT32:
			currentMin=(DTYPE)std::numeric_limits<float>::max();
			currentMax=-(DTYPE)std::numeric_limits<float>::max();
			break;
		case NIFTI_TYPE_FLOAT64:
			currentMin=(DTYPE)std::numeric_limits<double>::max();
			currentMax=-(DTYPE)std::numeric_limits<double>::max();
			break;
	}

    if(image->scl_slope==0) image->scl_slope=1.0f;

	imagePtr = static_cast<DTYPE *>(image->data);
	for(unsigned int index=0; index<image->nvox; index++){

		double value = (double)*imagePtr * image->scl_slope + image->scl_inter;
		if (value > 1.0) {
			printf("err\tValues > 1.0 are not supported\n");
		}
	    vol4d[index+(int)index_nt*(int)image->nvox] = (unsigned char)round(value*255);
		
		*imagePtr++;
	}

}

void
ptr_to_uint8_4D(nifti_image *image, unsigned char *vol4d, int index_nt)
{
	switch(image->datatype){
		case NIFTI_TYPE_UINT8:
			ptr_to_uint8_4D2<unsigned char>(image, vol4d, index_nt);
			break;
		case NIFTI_TYPE_INT8:
			ptr_to_uint8_4D2<char>(image, vol4d, index_nt);
			break;
		case NIFTI_TYPE_UINT16:
			ptr_to_uint8_4D2<unsigned short>(image, vol4d, index_nt);
			break;
		case NIFTI_TYPE_INT16:
			ptr_to_uint8_4D2<short>(image, vol4d, index_nt);
			break;
		case NIFTI_TYPE_UINT32:
			ptr_to_uint8_4D2<unsigned int>(image, vol4d, index_nt);
			break;
		case NIFTI_TYPE_INT32:
			ptr_to_uint8_4D2<int>(image, vol4d, index_nt);
			break;
		case NIFTI_TYPE_FLOAT32:
			ptr_to_uint8_4D2<float>(image, vol4d, index_nt);
			break;
		case NIFTI_TYPE_FLOAT64:
			ptr_to_uint8_4D2<double>(image, vol4d, index_nt);
			break;
		default:
			printf("err\tptr_to_uint8_4D2\tThe image data type is not supported\n");
			return;
	}
}

template <class DTYPE>
void ptr_to_double_4D2(nifti_image *image, double *vol4d, int index_nt)
{
    DTYPE *imagePtr = static_cast<DTYPE *>(image->data);
	DTYPE currentMin=0;
	DTYPE currentMax=0;
	switch(image->datatype){
		case NIFTI_TYPE_UINT8:
			currentMin=(DTYPE)std::numeric_limits<unsigned char>::max();
			currentMax=0;
			break;
		case NIFTI_TYPE_INT8:
			currentMin=(DTYPE)std::numeric_limits<char>::max();
			currentMax=(DTYPE)-std::numeric_limits<char>::max();
			break;
		case NIFTI_TYPE_UINT16:
			currentMin=(DTYPE)std::numeric_limits<unsigned short>::max();
			currentMax=0;
			break;
		case NIFTI_TYPE_INT16:
			currentMin=(DTYPE)std::numeric_limits<short>::max();
			currentMax=-(DTYPE)std::numeric_limits<short>::max();
			break;
		case NIFTI_TYPE_UINT32:
			currentMin=(DTYPE)std::numeric_limits<unsigned int>::max();
			currentMax=0;
			break;
		case NIFTI_TYPE_INT32:
			currentMin=(DTYPE)std::numeric_limits<int>::max();
			currentMax=-(DTYPE)std::numeric_limits<int>::max();
			break;
		case NIFTI_TYPE_FLOAT32:
			currentMin=(DTYPE)std::numeric_limits<float>::max();
			currentMax=-(DTYPE)std::numeric_limits<float>::max();
			break;
		case NIFTI_TYPE_FLOAT64:
			currentMin=(DTYPE)std::numeric_limits<double>::max();
			currentMax=-(DTYPE)std::numeric_limits<double>::max();
			break;
	}

    if(image->scl_slope==0) image->scl_slope=1.0f;

	imagePtr = static_cast<DTYPE *>(image->data);
	for(unsigned int index=0; index<image->nvox; index++){

		double value = (double)*imagePtr * image->scl_slope + image->scl_inter;
	    vol4d[index+(int)index_nt*(int)image->nvox] = value;
		
		*imagePtr++;
	}

}

void
ptr_to_double_4D(nifti_image *image, double *vol4d, int index_nt) 
{

	switch(image->datatype){
		case NIFTI_TYPE_UINT8:
			ptr_to_double_4D2<unsigned char>(image, vol4d, index_nt);
			break;
		case NIFTI_TYPE_INT8:
			ptr_to_double_4D2<char>(image, vol4d, index_nt);
			break;
		case NIFTI_TYPE_UINT16:
			ptr_to_double_4D2<unsigned short>(image, vol4d, index_nt);
			break;
		case NIFTI_TYPE_INT16:
			ptr_to_double_4D2<short>(image, vol4d, index_nt);
			break;
		case NIFTI_TYPE_UINT32:
			ptr_to_double_4D2<unsigned int>(image, vol4d, index_nt);
			break;
		case NIFTI_TYPE_INT32:
			ptr_to_double_4D2<int>(image, vol4d, index_nt);
			break;
		case NIFTI_TYPE_FLOAT32:
			ptr_to_double_4D2<float>(image, vol4d, index_nt);
			break;
		case NIFTI_TYPE_FLOAT64:
			ptr_to_double_4D2<double>(image, vol4d, index_nt);
			break;
		default:
			printf("err\tptr_to_double_4D2\tThe image data type is not supported\n");
			return;
	}
}