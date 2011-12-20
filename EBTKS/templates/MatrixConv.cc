/*--------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1996, Alex P. Zijdenbos, 
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
$RCSfile: MatrixConv.cc,v $
$Revision: 1.4 $
$Author: bert $
$Date: 2004-12-08 17:06:56 $
$State: Exp $
--------------------------------------------------------------------------*/
#include <config.h>
#include <iostream>
using namespace std;
#include "Matrix.h"
//#include "miscTemplateFunc.h"

//
//-------------------------// 
//

#ifdef USE_COMPMAT
template <class Type> 
Mat<dcomplex> asCompMat(const Mat<Type>& Re, const Mat<Type>& Im)
{
  if ((Im.getcols() != Re.getcols()) || (Im.getrows() != Re.getrows())) {
    cerr << "asCompMat: Re and Im matrices don't have the same dimensions; using Re only"
	 << endl;
    return asCompMat(Re);
  }

  Mat<dcomplex> cast(Re.getrows(), Re.getcols());
  dcomplex     *castPtr = (dcomplex *) cast.getEl()[0];
  const Type  *rePtr   = Re.getEl()[0];
  const Type  *imPtr   = Im.getEl()[0];
  for (unsigned i=Re.nElements(); i; i--)
    *castPtr++ = dcomplex(*rePtr++, *imPtr++);
  
  return cast;
}

template <class Type> 
Mat<dcomplex> asCompMat(const Mat<Type>& A)
{
  Mat<dcomplex> cast(A.getrows(), A.getcols());
  dcomplex     *castPtr = (dcomplex *) cast.getEl()[0];
  Type         *aPtr    = (Type *) A.getEl()[0];
  for (unsigned i=A.nElements(); i; i--)
    *castPtr++  = dcomplex(*aPtr++);
  
  return cast;
}
#endif

#ifdef USE_FCOMPMAT
template <class Type> 
Mat<fcomplex> asFcompMat(const Mat<Type>& A)
{
  Mat<fcomplex> cast(A.getrows(), A.getcols());
  fcomplex     *castPtr = (fcomplex *) cast.getEl()[0];
  const Type  *aPtr    = A.getEl()[0];
  for (unsigned i=A.nElements(); i; i--)
    *castPtr++  = fcomplex(*aPtr++);
  
  return cast;
}

template <class Type> 
Mat<fcomplex> asFcompMat(const Mat<Type>& Re, const Mat<Type>& Im)
{
  if ((Im.getcols() != Re.getcols()) || (Im.getrows() != Re.getrows())) {
    cerr << "asCompMat: Re and Im matrices don't have the same dimensions; using Re only"
	 << endl;
    return asFcompMat(Re);
  }

  Mat<fcomplex> cast(Re.getrows(), Re.getcols());
  fcomplex     *castPtr = (fcomplex *) cast.getEl()[0];
  const Type  *rePtr   = Re.getEl()[0];
  const Type  *imPtr   = Im.getEl()[0];
  for (unsigned i=Re.nElements(); i; i--)
    *castPtr++ = fcomplex(*rePtr++, *imPtr++);
  
  return cast;
}
#endif

template Mat<dcomplex> asCompMat(Mat<double> const &);
template Mat<dcomplex> asCompMat(Mat<dcomplex> const &);


#if !defined(__GNUC__) && defined(__sgi)
#include "MatrixSpec.cc"
#include "SimpleArraySpec.cc"

#ifdef USE_COMPMAT
template class Mat<dcomplex>;
template Mat<dcomplex> inv(const Mat<dcomplex> &);
template Mat<dcomplex> fft(const Mat<double> &, unsigned, unsigned);
template Mat<dcomplex> fft(const Mat<dcomplex> &, unsigned, unsigned);
template Mat<dcomplex> ifft(const Mat<double> &, unsigned, unsigned);
template Mat<dcomplex> ifft(const Mat<dcomplex> &, unsigned, unsigned);
template class SimpleArray<dcomplex>;
template class IndexStruct<dcomplex>;
#endif // USE_COMPMAT
#endif // !defined(__GNUC__) && defined(__sgi)

