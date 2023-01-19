/*--------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1996, John G. Sled, 
              McConnell Brain Imaging Centre,
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- 
$RCSfile: splineSmooth.cc,v $
$Revision$
$Author$
$Date$
$State: Exp $
--------------------------------------------------------------------------*/
/* ----------------------------- MNI Header -----------------------------------
@NAME       : splineSmooth.c,v
@INPUT      : 
@OUTPUT     : (none)
@RETURNS    : 
@DESCRIPTION: Tool for smoothing and extrapolating data in minc volumes
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : April 21, 1996 (John G. Sled)
@MODIFIED   : Log: splineSmooth.c,v 
 * Revision 1.2  1996/04/23  13:36:58  jgsled
 * Working version.  Problems with thin plate splines have been fixed.
 *
 * Revision 1.1  1996/04/21  23:41:50  jgsled
 * Initial version of SplineSmooth tool
 * - B spline implementation appears to work
 *
@COPYRIGHT  : 1996
---------------------------------------------------------------------------- */

#include <stdio.h>
#include <iostream>   
//using namespace std;    
using std::cerr;
#include <math.h>
#include <Matrix.h>  
#include <TBSpline.h> 
#undef ROUND
#undef SIGN

//---------------------------------------------------------------------------------
//  Implementation notes
/*
  The spline basis functions are defined in a world coordinate system aligned
  with the voxel coordinate system and sharing the same origin.

 */

//---------------------------------------------------------------------------------
// Declarations
int fitSplinesToVolumeLookup(TBSplineVolume *spline, float *src, 
                              int subsample, double *separations, int *dims);
void smoothVolumeLookup(TBSplineVolume *spline, float *src, int *dims);
//int splineSmooth( float *src, double lambda, double distance, int subsample, double *separations, int *dims);

//--------------------------------------------------------------------------------
// main program
extern "C" {
int splineSmooth( float *src, double lambda, double distance, int subsample, double *separations, int *dims)    
{
  int i;
  
  // create spline basis
  Spline *theSplines;
  
  double start[3] = { 0.0, 0.0, 0.0 };
  theSplines = new TBSplineVolume(start, separations, dims,
                                      distance, lambda, 1);
  
  // do least squares fit to data 
  if(fitSplinesToVolumeLookup((TBSplineVolume *)theSplines, src,
                               subsample, separations, dims) == TRUE) {
    // write smooth function to volume
    smoothVolumeLookup((TBSplineVolume *) theSplines, src, dims);
  } else {                               
    cerr << "Spline fit failed: No fitting is used.\n";
    for (i=0; i < dims[0]*dims[1]*dims[2]; i++) src[i] = 0.0;
  }
      
  return(0);
} 
}

//-----------------------------------------------------------------------------
// Supporting functions


int
fitSplinesToVolumeLookup(TBSplineVolume *spline, float* src,
                         int subsample, double* separations, int* dims)
{
  int i,x,y,z,j,k;
  long area, vol, z_area, y_dims;
  double value;

  area = dims[0]*dims[1];
  vol  = area*dims[2];

  for(z = 0; z < dims[2]; z += subsample) {
    z_area = z*area;
    for(y = 0; y < dims[1]; y += subsample) {
      y_dims = y*dims[0];
      for(x = 0; x < dims[0]; x += subsample) {
        value = (double)src[z_area + y_dims + x];
        if(value > 0)
          spline->addDataPoint(x,y,z, value);
      }
    }
  }

  if(spline->fit() == FALSE) // fit splines to the data
    {
      cerr << "Fatal Error: Spline fit failed.\n";
      return(FALSE);
    } else return(TRUE);
}

void 
smoothVolumeLookup(TBSplineVolume *spline, float* src, int* dims)
{
  int x,y,z;
  double value;
  long area, vol, z_area, y_dims;

  area = dims[0]*dims[1];
  vol  = area*dims[2];

  for (z = 0; z < dims[2]; z++) {
    z_area = z*area;
    for (y = 0; y < dims[1]; y++) {
      y_dims = y*dims[0];
      for (x = 0; x < dims[0]; x++) {
    value = (*spline)(x,y,z); 
    src[z_area + y_dims + x] = (float)value;
    }
    }
  }
}
