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
$RCSfile: Polynomial.cc,v $
$Revision: 1.3 $
$Author: stever $
$Date: 2003-11-17 04:07:52 $
$State: Exp $
--------------------------------------------------------------------------*/
#include <config.h>
#include "Polynomial.h"

using namespace std;



//
// Constructors
//

Polynomial::Polynomial(unsigned dim, const DblArray& coef)
: _expComb(dim, coef.size()),
  _coef(coef)
{
  unsigned maxOrder = unsigned(rint(pow(coef.size(), 1.0/dim)));
  if (intPower(maxOrder, dim) != coef.size()) {
    cerr << "# polynomial coefficients (" << coef.size()
	 << ") does not match polynomial dimension (" << dim << ")" << endl;
    exit(EXIT_FAILURE);
  }

  _allExpComb(dim, maxOrder - 1);

  _nDimensions = _expComb.getrows();
  _nCoef       = _expComb.getcols();
  assert(_nDimensions && _nCoef);
}

Polynomial::Polynomial(const IntMat& expComb, const DblArray& coef)
: _expComb(expComb),
  _coef(coef)
{
  _nDimensions = _expComb.getrows();
  _nCoef       = _expComb.getcols();

  if (coef.size() != _nCoef) {
    cerr << "Size of the exponent combinations array (" << _nCoef
	 << ") does not match the # coefficients (" << coef.size() << ")" << endl;
    exit(EXIT_FAILURE);
  }
  
  assert(_nDimensions && _nCoef);
}

// 1D fitted polynomials
Polynomial::Polynomial(unsigned mvo, const FloatArray& x, const DblArray& f)
: _expComb(1, mvo + 1)
{
  _nCoef = mvo + 1;

  for (unsigned i = 0; i < _nCoef; i++)
    _expComb(0, i) = i;
  
  _nDimensions = 1;

  _fit(FlMat(x, FlMat::COLUMN), f);
}

Polynomial::Polynomial(unsigned mvo, const IntArray& x, const DblArray& f)
: _expComb(1, mvo + 1)
{
  _nCoef = mvo + 1;

  for (unsigned i = 0; i < _nCoef; i++)
    _expComb(0, i) = i;
  
  _nDimensions = 1;

  _fit(FlMat(asFloatArray(x), FlMat::COLUMN), f);
}

// 2D fitted polynomials

Polynomial::Polynomial(unsigned mvo, const FloatArray& x, const FloatArray& y, 
		       const DblArray& f)
{
  _allExpComb(2, mvo);
  _pruneExpComb(mvo);
      
  _nDimensions = _expComb.getrows();
  _nCoef       = _expComb.getcols();
  assert(_nDimensions && _nCoef);

  FlMat X(x, FlMat::COLUMN);
  
  _fit(X.appendRight(FlMat(y, FlMat::COLUMN)), f);
}

Polynomial::Polynomial(unsigned mvo, const IntArray& x, const IntArray& y, 
		       const DblArray& f)
{
  _allExpComb(2, mvo);
  _pruneExpComb(mvo);
      
  _nDimensions = _expComb.getrows();
  _nCoef       = _expComb.getcols();
  assert(_nDimensions && _nCoef);

  FlMat X(asFloatArray(x), FlMat::COLUMN);

  _fit(X.appendRight(FlMat(asFloatArray(y), FlMat::COLUMN)), f);
}

// n-D fitted polynomials
Polynomial::Polynomial(unsigned mvo, const FlMat& x, const DblArray& f)
{
  _allExpComb(x.getcols(), mvo);
  _pruneExpComb(mvo);
  
  _nDimensions = _expComb.getrows();
  _nCoef       = _expComb.getcols();
  assert(_nDimensions && _nCoef);

  _fit(x, f);
}

Polynomial::Polynomial(unsigned mvo, const IntMat& x, const DblArray& f)
{
  _nDimensions = x.getcols();

  _allExpComb(_nDimensions, mvo);
  _pruneExpComb(mvo);

  _nCoef = _expComb.getcols();
  assert(_nDimensions && _nCoef);

  _fit(asFlMat(x), f);
}

#ifdef HAVE_INTERFACE
Polynomial::Polynomial(const MRegion& xy, const DblArray& f)
{
  _allExpComb(2, 1);
  _pruneExpComb(1);

  _nDimensions = _expComb.getrows();
  _nCoef       = _expComb.getcols();
  assert(_nDimensions && _nCoef);

  Array<FloatArray> X(2);
  X[0] = asFloatArray(xy.xCoord());
  X[1] = asFloatArray(xy.yCoord());
  
  _fit(X, f);
}

Polynomial::Polynomial(unsigned mvo, const MRegion& xy, const DblArray& f)
{
  _allExpComb(2, mvo);
  _pruneExpComb(mvo);
  
  _nDimensions = _expComb.getrows();
  _nCoef       = _expComb.getcols();
  assert(_nDimensions && _nCoef);

  Array<FloatArray> X(2);
  X[0] = asFloatArray(xy.xCoord());
  X[1] = asFloatArray(xy.yCoord());
  
  _fit(X, f);
}
#endif

//
// Evaluation functions
//

double
Polynomial::operator () (float x) const
{
  if (_nDimensions != 1) {
    cerr << "Polynomial::operator (): Error: cannot evaluate a " << _nDimensions 
	 << "-dimensional polynomial with 2 coordinates." << endl;
    return 0.0;
  }

  double result = 0.0;
  const double *coefPtr  = _coef.contents();
  const int    *powerPtr = _expComb[0];
  for (unsigned i = _nCoef; i; i--)
    result += *coefPtr++ * intPower(x, *powerPtr++);
  
  return result;
}

double
Polynomial::operator () (float x, float y) const
{
  if (_nDimensions != 2) {
    cerr << "Polynomial::operator (): Error: cannot evaluate a " << _nDimensions 
	 << "-dimensional polynomial with 2 coordinates." << endl;
    return 0.0;
  }

  double result = 0.0;
  const  double *coefPtr = _coef.contents();
  const  int    *xPwrPtr = _expComb[0];
  const  int    *yPwrPtr = _expComb[1];
  for (unsigned i = _nCoef; i; i--)
    result += *coefPtr++ * intPower(x, *xPwrPtr++) * intPower(y, *yPwrPtr++);

  return result;
}

double
Polynomial::operator () (float x, float y, float z) const
{
  if (_nDimensions != 3) {
    cerr << "Polynomial::operator (): Error: cannot evaluate a " << _nDimensions 
	 << "-dimensional polynomial with 3 coordinates." << endl;
    return 0.0;
  }

  double result = 0.0;
  const  double *coefPtr = _coef.contents();
  const  int    *xPwrPtr = _expComb[0];
  const  int    *yPwrPtr = _expComb[1];
  const  int    *zPwrPtr = _expComb[2];
  for (unsigned i = 0; i < _nCoef; i++)
    result += *coefPtr++ * 
      intPower(x, *xPwrPtr++) * 
      intPower(y, *yPwrPtr++) * 
      intPower(z, *zPwrPtr++);

  return result;
}

double
Polynomial::operator () (const FloatArray& X) const
{
  if (_nDimensions != X.size()) {
    cerr << "Polynomial::operator (): Error: cannot evaluate a " << _nDimensions 
	 << "-dimensional polynomial with " << X.size() << " coordinates." << endl;
    return 0.0;
  }

  double result = 0.0;
  for (unsigned i = 0; i < _nCoef; i++) {
    double temp = _coef[i];
    const float *xPtr = X.contents();
    for (unsigned coord = 0; coord < _nDimensions; coord++)
      temp *= intPower(*xPtr++, _expComb[coord][i]);
    result += temp;
  }

  return result;
}

//
// Private functions
//

void
Polynomial::_allExpComb(unsigned dimension, unsigned maxOrder)
{
// Set up <_expComb>

  unsigned nExpComb = (unsigned) intPower(maxOrder + 1, dimension);
  _expComb = IntMat(dimension, nExpComb);
  _expComb.fill(0);

// Determine all possible combinations of exponents of the variables
// (i.e., create a number series with base <maxOrder + 1>)

  for (unsigned i = 1; i < nExpComb; i++) {
    Boolean addOne = TRUE;
    for (unsigned coord = 0; coord < dimension; coord++) {
      int *coordArray = (int *) _expComb[coord];
      if (addOne) {
	if (coordArray[i-1] == maxOrder)
	  coordArray[i] = 0;
	else {
	  coordArray[i] = coordArray[i-1] + 1;
	  addOne = FALSE;
	}
      }
      else
	coordArray[i] = coordArray[i-1];
    }
  }
}

void
Polynomial::_pruneExpComb(unsigned mvo)
{
  unsigned nDimensions = _expComb.getrows();
  unsigned nExpComb    = _expComb.getcols();

  if (!nDimensions || !nExpComb)
    return;

  unsigned dim;
  unsigned col = 0;
  for (unsigned i = 0; i < nExpComb; i++) {
    unsigned expSum = 0;
    for (dim = 0; dim < nDimensions; dim++)
      expSum += _expComb(dim, i);

    if (expSum <= mvo) {
      if (col != i) 
	for (dim = 0; dim < nDimensions; dim++)
	  _expComb(dim, col) = _expComb(dim, i);
      col++;
    }
  }

  _expComb.resize(nDimensions, col);
}

void
Polynomial::_fit(const FlMat& X, const DblArray& F)
{
  if (X.getcols() != _nDimensions) {
    cerr << "Error in polynomial::_fit: _expComb and X have different sizes" 
	 << endl;
    return;
  }

  typedef Array<DblArray> ArrayOfDblArrays;

  IntArray                nPowers(_nDimensions);
  Array<ArrayOfDblArrays> powers(_nDimensions);

  // Determine the maximum exponent of each variable

  for (unsigned coord = 0; coord < _nDimensions; coord++) {
    int nPow = nPowers[coord] = 2*_expComb.row(coord).max() + 1;
    powers[coord] = ArrayOfDblArrays(nPow);

    // Raise X to all necessary powers

    for (unsigned power = 0; power < nPow; power++)
      powers[coord][power] = asDblArray(array(X.col(coord))) ^ power;
  }

  DblMat A(_nCoef, _nCoef);
  DblMat B(_nCoef, 1);

  for (unsigned row = 0; row < _nCoef; row++) {
    unsigned exponent = _expComb[0][row];
    DblArray XrowPower(powers[(unsigned int)0][exponent]);
    for (unsigned coord = 1; coord < _nDimensions; coord++) {
      exponent = _expComb[coord][row];
      XrowPower *= powers[coord][exponent];
    }

    for (unsigned col = 0; col <= row; col++) {
      DblArray Xpower(XrowPower);

      for (unsigned coord = 0; coord < _nDimensions; coord++) {
	exponent = _expComb[coord][col];
	Xpower *= powers[coord][exponent];
      }

      A(row, col) = sum(Xpower);
      if (row != col)
	A(col, row) = A(row, col);
    }

    B(row) = sum(XrowPower * F);
  }

  DblMat coefMat(inv(A)*B);

  _coef = DblArray(_nCoef);

  for (unsigned i = 0; i < _nCoef; i++)
    _coef[i] = coefMat(i);
}
