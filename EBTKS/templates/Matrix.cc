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
$RCSfile: Matrix.cc,v $
$Revision: 1.7 $
$Author: rotor $
$Date: 2010-05-18 23:01:21 $
$State: Exp $
--------------------------------------------------------------------------*/
#include <config.h>
#include "Matrix.h"
#include "FileIO.h"  //included to load compressed data
#include <sys/stat.h>
#include "unistd.h"
#ifdef HAVE_MATLAB
#include "matlabSupport.h"
#endif
#include "miscTemplateFunc.h"

using namespace std;


#ifdef USE_COMPMAT
#ifdef __GNUC__
#pragma do_not_instantiate Mat<dcomplex>::array
#pragma do_not_instantiate Mat<dcomplex>::applyIndexFunction
#endif
#endif

#ifdef USE_FCOMPMAT
#ifdef __GNUC__
#pragma do_not_instantiate Mat<fcomplex>::array
#endif
#endif


/*********************************************************
template <class Type> const Mat<Type>::MATLAB = 0;
template <class Type> const Mat<Type>::RAW = 1;
template <class Type> const Mat<Type>::ASCII  = 2;
//template <class Type>  tells belong to template class 
*********************************************************/

#ifndef __GNUC__
template <class Type> unsigned Mat<Type>::_rangeErrorCount = 25;
template <class Type> Boolean  Mat<Type>::flushToDisk = FALSE;
#endif

//
// Some mathematical operations
//

template <class T1, class T2>
Mat<T1> operator * (const Mat<T1>& A, const Mat<T2>& B)
{
  unsigned arows = A.getrows();
  unsigned acols = A.getcols();

  unsigned brows    = B.getrows();
  unsigned bcols    = B.getcols();
  unsigned bmaxcols = B.getmaxcols();

  Mat<T1> Temp(arows, bcols);

  if (acols != brows) {
    cerr << "Mat sizes incompatible for *" << endl;
    return Temp;
  }
  
  T1        *tempPtr = (T1 *) Temp.getEl()[0];
  const T1 **rowPtr  = A.getEl();
  
  for (unsigned i = arows; i; i--) {
    const T2 *argColPtr = B.getEl()[0];
    
    for (unsigned j = bcols; j; j--) {
      const T1 *ptr    = *rowPtr;
      const T2 *argPtr = argColPtr++;
      *tempPtr = (T1) 0;
      
      for (unsigned k = acols; k; k--) {
	*tempPtr += *ptr++ * (T1) *argPtr;
	argPtr += bmaxcols;
      }
      tempPtr++;
    }
    rowPtr++;
  }
  
  return Temp;
}

template <class T1, class T2>
Mat<T1>& pmultEquals(Mat<T1>& A, const Mat<T2>& B)
{
  unsigned nrows = A.getrows();
  unsigned ncols = A.getcols();

  if (!(A.isvector() && B.isvector() && (A.length() == B.length())) &&
      ((B.getrows() != nrows) || (B.getcols() != ncols))) {
    cerr << "Matrices of incompatible sizes for pmultEquals" << endl;
    return A;
  }

  T1       *aPtr = (T1 *) A.getEl()[0];
  const T2 *bPtr = B.getEl()[0];

  for (unsigned i = nrows; i; i--)
    for (unsigned j = ncols; j; j--)
      *aPtr++ *= (T1) *bPtr++;
  
  return A;
}

template <class T1, class T2>
Mat<T1>& pdivEquals(Mat<T1>& A, const Mat<T2>& B)
{
  unsigned nrows = A.getrows();
  unsigned ncols = A.getcols();

  if (!(A.isvector() && B.isvector() && (A.length() == B.length())) &&
      ((B.getrows() != nrows) || (B.getcols() != ncols))) {
    cerr << "Matrices of incompatible sizes for pdivEquals" << endl;
    return A;
  }

  T1       *aPtr = (T1 *) A.getEl()[0];
  const T2 *bPtr = B.getEl()[0];

  for (unsigned i = nrows; i; i--)
    for (unsigned j = ncols; j; j--)
      *aPtr++ /= (T1) *bPtr++;
  
  return A;
}

/*********************************
Mat class definitions and member functions
**********************************/
//
//-------------------------// 
//

template <class Type>
Mat<Type>::Mat(const Mat<Type>& A)
{
  _rows = A._rows; _cols = A._cols;
  _maxrows = A._maxrows; _maxcols = A._maxcols;
  _el = 0;
   
  _allocateEl();

  if (_maxrows && _maxcols && _el)
    memcpy(*_el, *A._el, _maxrows*_maxcols*sizeof(Type));
}

//
//-------------------------// 
//

template <class Type>
Mat<Type>::Mat(const SimpleArray<Type>& A, char dir)
{
  unsigned N = A.size();
  
  if (dir == Mat<Type>::ROW) {
    _rows = _maxrows = 1;
    _cols = _maxcols = N;
    _el = 0;
    _allocateEl();
  }
  else {
    _rows = _maxrows = N;
    _cols = _maxcols = 1;
    _el = 0;
    _allocateEl();
  }
  
  Type *elPtr = *_el;
  for (unsigned i = 0; i < N; i++)
    *elPtr++ = A[i];
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>::Mat(unsigned nrows, unsigned ncols, Type value)
{
   _rows = _maxrows = nrows;
   _cols = _maxcols = ncols;
   _el = 0;
   
   _allocateEl();
   
   if (value != 0.0)
      fill(value);
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>::Mat(unsigned nrows, unsigned ncols)
{
   _rows = _maxrows = nrows;
   _cols = _maxcols = ncols;
   _el = 0;
   
   _allocateEl();
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>::Mat(unsigned nrows, unsigned ncols, const Type *data)
{
   _rows = _maxrows = nrows;
   _cols = _maxcols = ncols;
   _el = 0;
   
   _allocateEl();
   
   size_t  nBytesInRow = ncols*sizeof(Type);
   const Type *dataPtr = data;
   for (unsigned row = 0; row < nrows; row++) {
      memcpy(_el[row], dataPtr, nBytesInRow);
      dataPtr += ncols;
   }
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>::Mat(const char *filename, int type, const char *varname)
{
  _rows = _maxrows = 0;
  _cols = _maxcols = 0;
  _el = 0;
  
  load(filename, type, varname);
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>::~Mat()
{
  clear();
}

//
//-------------------------// 
//

template <class Type>
void
Mat<Type>::clear()
{
   if (_el) {
      #ifdef DEBUG
      cout << "Freeing allocated memory at " << _el << endl;
      #endif
      if( _el[0] ) {
        delete [] _el[0];
      }
      delete [] _el; // delete row of pointers  
      _el = 0;
   }

   _rows = _maxrows = _cols = _maxcols = 0;
}

/****************************************************************************
These overloaded operators allow accessing of individual elements within the
matrix.  By convention, the first element of the matrix is labeled element 
zero.  Since the matrix is stored in row major format, The element numbers increase going to the right and at the following row. The element at the
last row and last column corresponds to the largest possible value of the 
input argument.
*****************************************************************************/


/****************************************************************************
This accessing operator allows modification of individual elements of the 
matrix. A seperate function is used for accessing elements within
complex matrices.
*****************************************************************************/
template <class Type>
Type& 
Mat<Type>::operator () (unsigned n) 
{
  if (n >= _rows*_cols) {
    if (_rangeErrorCount) {
      cerr << "Error: index " << n << " exceeds matrix dimensions. ";
      cerr << "Changed to " << _rows*_cols - 1 << endl;
      _rangeErrorCount--;
    }
    n = _rows*_cols - 1;
  }

  unsigned int rowindex;
  unsigned int colindex;
  rowindex= (n/_cols);
  colindex= (n%_cols);
  return _el[rowindex][colindex];
}

/****************************************************************************
This accessing operator prevents accidental modification by having a constant
output argument. A seperate function is used for accessing elements within
complex matrices.
*****************************************************************************/

template <class Type>
Type
Mat<Type>::operator () (unsigned n) const 
{
  if (n >= _rows*_cols) {
    if (_rangeErrorCount) {
      cerr << "Error: index " << n << " exceeds matrix dimensions. ";
      cerr << "Changed to " << _rows*_cols - 1 << endl;
      _rangeErrorCount--;
    }
    n = _rows*_cols - 1;
  }

  unsigned int rowindex;
  unsigned int colindex;
  rowindex= (n/_cols);
  colindex= (n%_cols);
  return _el[rowindex][colindex];
}

//
//-------------------------// 
//
template <class Type>
Type& 
Mat<Type>::operator () (unsigned r, unsigned c) {
  using std::min;		// (bert) force it to find the right min()

  if ((r >= _rows) || (c >= _cols)) {
    if (_rangeErrorCount) {
      cerr << "Error: indices (" << r << ", " << c << ") exceed matrix dimensions. "
	   << "Changed to (" << min(r, _rows - 1) << ", " << min(c, _cols - 1)
	   << ")" << endl;
      _rangeErrorCount--;
    }
    r = min(r, _rows - 1);
    c = min(c, _cols - 1);
  }

  return(_el[r][c]);
}

//
// 
//
template <class Type>
Type
Mat<Type>::operator () (unsigned r, unsigned c) const {
  using std::min;		// (bert) force it to find the right min().

  if ((r >= _rows) || (c >= _cols)) {
    if (_rangeErrorCount) {
      cerr << "Error: indices (" << r << ", " << c << ") exceed matrix dimensions. "
	   << "Changed to (" << min(r, _rows - 1) << ", " << min(c, _cols - 1)
	   << ")" << endl;
      _rangeErrorCount--;
    }
    r = min(r, _rows - 1);
    c = min(c, _cols - 1);
  }

  return(_el[r][c]);
}

//
//------------------// 
//
//cropt a section of the matrix
template <class Type>
Mat<Type>
Mat<Type>::operator () (unsigned r1, unsigned r2, unsigned c1, unsigned c2) const {
  if ((r1 > r2) || (c1 > c2) || (r2 >= _rows) || (c2>= _cols)){
    cerr << "Error in cropping: improper row or column sizes." << endl;
    cerr << r1 << " to " << r2 << " and" << endl;
    cerr << c1 << " to " << c2 << endl;
    exit(1);
  }
  Mat<Type> A(r2-r1+1,c2-c1+1);
  Type *Aptr = A._el[0];
  for (unsigned i=r1; i<=r2 ; i++){
    for (unsigned j=c1; j<=c2; j++)
      *Aptr++=_el[i][j];
  }
  return(A);     
}

//
//-------------------------// 
//
template <class Type>
const Type * 
Mat<Type>::operator [] (unsigned index) const {
  if (index >= _rows) {
    if (_rangeErrorCount) {
      cerr << "Error: index " << index << " exceeds matrix dimensions. ";
      cerr << "Changed to " << _rows - 1 << endl;
      _rangeErrorCount--;
    }

    index = _rows - 1;
  }

  return _el[index];
}

//
//-------------------------// 
//
template <class Type>
void
Mat<Type>::operator () (const Mat<Type>& A)
{
  if ((A._maxrows != _maxrows) || (A._maxcols != _maxcols)) {
    _maxrows = A._maxrows;
    _maxcols = A._maxcols;
    _allocateEl();
  }

  _rows = A._rows;
  _cols = A._cols;

  if (_maxrows && _maxcols && _el)
    memcpy(*_el, *A._el, _maxrows*_maxcols*sizeof(Type));
}

//
//-------------------------// 
//


/**********************************************************************
Assignment operators for two dimensional matrix objects
***********************************************************************/

template <class Type>
Mat<Type>& 
Mat<Type>::operator = (const Mat<Type>& A)
{
  if (this == &A) 
    return *this;

  if ((A._maxrows != _maxrows) || (A._maxcols != _maxcols)) {
    _maxrows = A._maxrows;
    _maxcols = A._maxcols;
    _allocateEl();
  }
  
  _rows = A._rows;
  _cols = A._cols;

  if (_maxrows && _maxcols && _el)
    memcpy(*_el,*A._el, _maxrows*_maxcols*sizeof(Type));

  return *this;
}

template <class Type>
Mat<Type>& 
Mat<Type>::absorb(Mat<Type>& A)
{
  if (this == &A)
    return *this;

  // Delete current contents;
  if (_el) {
    if( _el[0] ) {
      delete [] _el[0];  // delete data block
    }
    delete [] _el; // delete row of pointers  
  }
  
  // Copy all from A
  _maxrows = A._maxrows;
  _maxcols = A._maxcols;
  _rows    = A._rows;
  _cols    = A._cols;
  _el      = A._el;
  
  // Empty A
  A._maxrows = A._maxcols = A._rows = A._cols = 0;
  A._el = 0;

  return *this;
}

//
//-------------------------// 
//

template <class Type>
Mat<Type>&
Mat<Type>::pad(unsigned nrows, unsigned ncols, int row, int col, Type value)
{
  if ((nrows == _rows) && (ncols == _cols) && !row && !col)
    return *this;

  char tempFile[256];
  get_temp_filename(tempFile);

  // Save the current matrix to disk
  if (flushToDisk && saveRaw(tempFile)) {
    unsigned argRows = _rows;
    unsigned argCols = _cols;

    clear();

    _rows = _maxrows = nrows;
    _cols = _maxcols = ncols;

    _allocateEl();
    fill(value);

    insert(tempFile, argRows, argCols, row, col);
  }
  else {
    // Saving fails, create a new matrix
    Mat<Type> result(nrows, ncols, value);
    result.insert(*this, row, col);
    
    this->absorb(result);
  }

  unlink(tempFile);

  return *this;
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>&
Mat<Type>::resize(unsigned nrows, unsigned ncols, int row, int col)
{
  if ((nrows == _rows) && (ncols == _cols))
    return *this;

  // Origin maintained, shrinking: no reallocation necessary
  if (!row && !col && (nrows <= _maxrows) && (ncols <= _maxcols)) {
    _rows = nrows;
    _cols = ncols;

    cerr << "This type of resizing is insecure!! Should be fixed..." << endl;

    return *this;
  }
  
  // This could be done by flushing to disk in order to save memory
  Mat<Type> T(nrows, ncols);

  Type *elPtr = T._el[0];
  for (int i = 0, i0 = row; i < nrows; i++, i0++) {
    Boolean valid = (i0 >= 0) && (i0 < _rows);
    for (int j = 0, j0 = col; j < ncols; j++, j0++)
      *elPtr++ = (valid && (j0 >= 0) && (j0 < _cols)) ? _el[i0][j0] : Type(0);
  }

  absorb(T);

  return *this;
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>&
Mat<Type>::resize()
{
  _rows = _maxrows;
  _cols = _maxcols;
  
  cerr << "This type of resizing is insecure!! Should be fixed..." << endl;

  return *this;
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>
Mat<Type>::appendRight(const Mat<Type>& A) const
{
  using std::max;		// (bert) force it to find the right max().

  Mat<Type> Temp(max(_rows, A._rows), _cols + A._cols);
  Temp.insert(*this);
  Temp.insert(A, 0, _cols);

  return Temp;
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>
Mat<Type>::appendBelow(const Mat<Type>& A) const
{
  using std::max;		// (bert) force it to use the right max()

  Mat<Type> Temp(_rows + A._rows, max(_cols, A._cols));
  Temp.insert(*this);
  Temp.insert(A, _rows, 0);

  return Temp;
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>&
Mat<Type>::swapRows(unsigned r1, unsigned r2)
{
  if (r1 == r2)
    return *this;

  if ((r1 >= _rows) || (r2 >= _rows)) {
    cerr << "Error in swapRows: improper row indices " << r1 
	 << "," << r2 << " for matrix with " << _rows << " rows" << endl;

    return *this;
  }

  Type *r1ptr = _el[r1];
  Type *r2ptr = _el[r2];

  for (unsigned i = _cols; i; i--, r1ptr++, r2ptr++)
    swap(*r1ptr, *r2ptr);

  return *this;
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>&
Mat<Type>::swapCols(unsigned c1, unsigned c2)
{
  if (c1 == c2)
    return *this;

  if ((c1 >= _cols) || (c2 >= _cols)) {
    cerr << "Error in swapCols: improper column indices " << c1 
	 << "," << c2 << " for matrix with " << _cols << " cols" << endl;

    return *this;
  }

  for (unsigned r = 0; r < _rows; r++)
    swap(_el[r][c1], _el[r][c2]);
  
  return *this;
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>
Mat<Type>::residual(unsigned row, unsigned col) const
{
  if ((_rows <= 1) || (_cols <= 1) || (row >= _rows) || (col >= _cols)) {
    cerr << "Error: residual(" << row << ", " << col << ") of " 
	 << _rows << "x" << _cols << " matrix." << endl;
    return *this;
  }

  Mat<Type> result(_rows - 1, _cols - 1);

  Type **rowPtr       = _el;
  Type **resultRowPtr = result._el;

  for (unsigned i = 0; i < _rows; i++) {
    Type *elPtr = *rowPtr++;
    if (i != row) {
      Type *resultElPtr = *resultRowPtr++;
      for (unsigned j = 0; j < _cols; j++) {
	if (j != col)
	  *resultElPtr++ = *elPtr;
	elPtr++;
      }
    }
  }
   
  return result;
}

//
//-------------------------//
//
template <class Type>
Mat<Type>&
Mat<Type>::insert(const Mat<Type>& A, int row, int col)
{
  // This could be written more CPU-efficient
  const Type **ArowPtr = (const Type **) A._el;
  int          destRow = row;
  for (unsigned i = 0; i < A._rows; i++, destRow++) {
    Boolean valid = (destRow >= 0) && (destRow < _rows);
    const Type *AelPtr  = *ArowPtr++;
    int         destCol = col;
    for (unsigned j = 0; j < A._cols; j++, destCol++, AelPtr++)
      if (valid && (destCol >= 0) && (destCol < _cols))
	_el[destRow][destCol] = *AelPtr;
  }

  return *this;
}

template <class Type>
Mat<Type>&
Mat<Type>::insert(const char *path, unsigned nrows, unsigned ncols, int row, int col)
{
  InputFile argFile(path);
  if (!argFile) {
    cerr << "Couldn't open file " << path << endl;
    return *this;
  }

  _checkMatrixDimensions(path, nrows, ncols);

  // Create a row-sized buffer
  Type *buffer = 0;
  allocateArray(ncols, buffer);
  if (!buffer) {
    cerr << "Couldn't allocate buffer" << endl;
    return *this;
  }

  int destRow = row;
  for (unsigned i = 0; i < nrows; i++, destRow++) {
    if (!(argFile.stream().read((char *) buffer, ncols*sizeof(Type)))) {
      cerr << "Error while reading file " << path << endl;
      freeArray(buffer);
      return *this;
    }
    Boolean valid = (destRow >= 0) && (destRow < _rows);
    const Type *bufferPtr = buffer;
    int         destCol = col;
    for (unsigned j = 0; j < ncols; j++, destCol++, bufferPtr++) {
      if (valid && (destCol >= 0) && (destCol < _cols))
	_el[destRow][destCol] = *bufferPtr;
    }
  }

  freeArray(buffer);
  
  return *this;
}

//
//-------------------------//
//
template <class Type>
Mat<Type>&
Mat<Type>::fill(Type value)
{
  Type *elPtr = _el[0];
  for (unsigned i = _rows ; i ; i--)
    for (unsigned j = _cols ; j ; j--) 
      *elPtr++ = value ;
  
  return *this;
}

//
//-------------------------//
//
template <class Type>
Mat<Type>&
Mat<Type>::fill(Type value, unsigned r1, unsigned r2, unsigned c1, unsigned c2)
{
  if ((r2 < r1) || (r2 >= _rows) || (c2 < c1) || (c2 >= _cols)) {
    cerr << "Error in Mat::fill: invalid row or column arguments." << endl;
    cerr << r1 << " to " << r2 << " and" << endl;
    cerr << c1 << " to " << c2 << endl;
    exit(1);
  }
    
  Type **rowPtr = _el + r1;
  for (unsigned i = r2 - r1 + 1 ; i ; i--) {
    Type *elPtr = *rowPtr++ + c1;  
    for (unsigned j = c2 - c1 + 1 ; j ; j--) 
      *elPtr++ = value ;
  }
  
  return *this;
}

//
//-------------------------//
//
template <class Type>
Mat<Type>&
Mat<Type>::fillEllips(Type value, double rowDiameter, double colDiameter)
{
  double row = double(_rows - 1)/2;
  double col = double(_cols - 1)/2;

  if (rowDiameter <= 0)
    rowDiameter = _rows;
  if (colDiameter <= 0)
    colDiameter = _cols;

  double rowRadiusSqr = SQR(rowDiameter/2);
  double colRadiusSqr = SQR(colDiameter/2);

  Type *elPtr = _el[0];
  for (unsigned i = 0; i < _rows; i++) {
    double rowDist = SQR(i - row)/rowRadiusSqr;
    for (unsigned j = 0; j < _cols; j++, elPtr++) 
      if (rowDist + SQR(j - col)/colRadiusSqr <= 1)
	*elPtr = value;
  }
  
  return *this;
}

//
//-------------------------//
//
template <class Type>
Mat<Type>&
Mat<Type>::fillEllips(Type value, double row, double col, double rowDiameter, 
		      double colDiameter)
{
  if (rowDiameter <= 0)
    rowDiameter = 2*MIN(row + 0.5, _rows - row - 0.5);
  if (colDiameter <= 0)
    colDiameter = 2*MIN(col + 0.5, _cols - col - 0.5);

  double rowRadiusSqr = SQR(rowDiameter/2);
  double colRadiusSqr = SQR(colDiameter/2);

  Type *elPtr = _el[0];
  for (unsigned i = 0; i < _rows; i++) {
    double rowDist = SQR(i - row)/rowRadiusSqr;
    for (unsigned j = 0; j < _cols; j++, elPtr++) 
      if (rowDist +  SQR(j - col)/colRadiusSqr <= 1)
	*elPtr = value;
  }
  
  return *this;
}

//
//-------------------------//
//
template <class Type>
Mat<Type>&
Mat<Type>::eye()
{
  unsigned mindim = (_rows > _cols) ? _cols : _rows;
  
  Type **dataPtr = _el;
  
  for(unsigned j=0; j< _rows;j++)
    for(unsigned k=0;k<_cols;k++)
      dataPtr[j][k]=(Type) 0;
  
  for (unsigned i=0 ; i < mindim ; i++)
    dataPtr[i][i] = (Type)1.0;
  
  return *this;
}

//
//-------------------------//
//
template <class Type>
Mat<Type>&
Mat<Type>::randuniform(double min, double max)
{
  double range = max - min;

  Type *elPtr = _el[0];
  for(unsigned i = _rows; i; i--)
    for(unsigned j = _cols; j; j--)
      *elPtr++ = Type(drand48() * range + min);
  
  return *this;
}

//
//-------------------------//
//
template <class Type>
Mat<Type>&
Mat<Type>::randnormal(double mean, double std)
{
   Type *elPtr = _el[0];
   for(unsigned i = _rows; i; i--)
     for(unsigned j = _cols ; j; j--)
       *elPtr++ = Type(gauss(mean, std));
   
   return *this;
}
//
//-------------------------//
//
template <class Type>
Mat<Type>&
Mat<Type>::hamming()
{
   //note len = _rows
   unsigned len=_rows;
   double alphaIncrement = 8.0*atan(1.0)/(len - 1);
   double alpha = 0;
   
   for(unsigned i=0 ; i < len ; i++){
      _el[i][0] = (Type)(0.54 - 0.46*::cos(alpha));
      alpha += alphaIncrement;
   }
   
   return *this;
}
//
//-------------------------//
//
template <class Type>
Mat<Type>&
Mat<Type>::hanning()
{
   unsigned len=_rows;
   double alphaIncrement = 8.0*atan(1.0)/(len-1);
   double alpha = 0;
   
   for(unsigned i=0 ; i < len; i++) {
      _el[i][0] = (Type)(0.5 - 0.5*::cos(alpha));
      alpha += alphaIncrement;
   }
   
   return *this;
}
//
//-------------------------//
//
template <class Type>
Mat<Type>&
Mat<Type>::blackman()
{
   unsigned len=_rows;
   double alphaIncrement = 8.0*atan(1.0)/(len - 1);
   double alpha = 0.0;
   
   for(unsigned i=0 ; i < len ; i++) {
      _el[i][0] = (Type)(0.42 - 0.5*::cos(alpha) + 0.08*::cos((2*alpha)));
      alpha += alphaIncrement;
   }
   
   return *this;
}


//This function creates a diagnonal matrix from the existing row or column 
//vector object

template <class Type> 
Mat<Type> 
Mat<Type>::diag() const
{ 
//Ensure that calling object is a vector
  unsigned dimensions=0;

  if(this->isvector())
   { 
//Size of new square matrix will equal the size of the vector 
    dimensions=this->length();  
   }  
  else
    {
      cerr<<"Error:calling object is not a row or column vector"<<endl;
      exit(1);
    }
   
  Mat<Type> Diagonalized_Matrix(dimensions,dimensions,(Type)0);
  for(unsigned i=0 ; i < dimensions ; i++)
//The () operator for the calling object is used to access each element in order    
    Diagonalized_Matrix(i,i)=this->operator() (i);
   
  return(Diagonalized_Matrix);
}  

//This function creates a diagnonal matrix from the existing row or column 
//vector object

/****************************************************************************
This function creates a toeplitz matrix from the two input vectors

   ie.  if c=[a b c] and r=[A B C], the result would be  a B C
                                                         b a B
                                                         c b a
*****************************************************************************/
 
template<class Type>
Mat<Type>
Mat<Type>::toeplitz(const Mat<Type>& r) const
{
//Make sure that both input arguments are vectors
  if(!(r.isvector() && this->isvector() ))
   {
     cerr<<"Error:One or both of the input arguments is/are not a vector"<<endl;
     exit(1);
   }
//temporarily store all the lengths so we can safely and quickly use them in
//math
  unsigned c_length= this->length();
  unsigned r_length=r.length();
  unsigned Row_Buffer_length=(r_length - 1 + c_length);
  Mat<Type> Row_Buffer(1,Row_Buffer_length,(Type)0);

//Fill the row buffer with the two vectors concatenated
//The vector (calling object) c is concatenated backwards to the front of 
//the row buffer  

  for(unsigned i=0;i < c_length;i++)
    Row_Buffer(i)= this->operator () (c_length - 1 - i);
  for(unsigned j=1; j < r_length; j++)
    Row_Buffer(c_length - 1 + j)=r(j);
  
  Mat<Type> Toeplitz_M(1,r_length,(Type)0);
//Get first row of toeplitz matrix  
  Toeplitz_M=Row_Buffer(0,0,c_length -1,Row_Buffer_length -1);
//Construct the matrix by rolling and appending sections of Row_Buffer
  
  for(unsigned k=2;k <= c_length;k++)
    Toeplitz_M = Toeplitz_M.appendBelow(Row_Buffer(0,0,c_length-k,Row_Buffer_length-k));
  
  return(Toeplitz_M);   
} 

/*********************************************************************/
//Display functions for two dimensional matrix objects
/*********************************************************************/

/****************************************************/
//This function displays the entire matrice
//This one does not work for <type> complex elements
/***************************************************/

template <class Type>
ostream&
Mat<Type>::display(ostream& os) const
{
  for (unsigned i=0 ; i < _rows ; i++) {
    for (unsigned j=0 ; j < _cols ; j++)
      os << _el[i][j] << " ";
    os << endl;
  }
  return os;
}

/**********************************************************************
This function will display a subsection of the matrix determined by 
input arguments. The first argument is the starting row(the first row
is row 0) followed by the ending row.  The third and fourth argument
are analagous for the start and end columns. A seperate function is
used to display matrices with complex elements
***********************************************************************/

template <class Type>
ostream&
Mat<Type>::display(ostream& os, unsigned r1, unsigned r2, unsigned c1, unsigned c2) const
{
//Make sure that the stop row and stop column are greater than the 
//start row and start column

  if ((r1 > r2) || (c1 > c2)){
    cerr << "Error in display: improper row or column sizes." << endl;
    cerr << r1 << " to " << r2 << " and" << endl;
    cerr << c1 << " to " << c2 << endl;
    exit(1);
  }

//Check to see if the stop row and stop column are exclusively less than
//the actual number of rows and columns in the matrix.  This
//prevents accidental access to an undefined memory location

  if ((r2 >=_rows) || (c2 >=_cols)) {
    cerr<<"The requested _rows or columns are not defined for this ";
    cerr<<"matrix"<<endl;
    exit(1);
  } 
  
  for (unsigned i=r1 ; i <= r2 ; i++) {
    for (unsigned j=c1 ; j <= c2 ; j++)
      os << _el[i][j] << " ";
    os << endl;
  }
  return os;
}

/**********************************************************************
This function will display a subsection of the matrix determined by 
input arguments.  The first argument is the starting row(the first row
is row 0) followed by the ending row.  The third and fourth argument
are analagous for the start and end columns.  This works or only matrices
with complex elements.
***********************************************************************/

//
//------------
//
template <class Type>
int
Mat<Type>::operator != (const Mat<Type>& Arg) const
{
  if ((_rows != Arg._rows) || (_cols != Arg._cols))
    return 1;

  Type **argRowPtr = Arg._el;
  Type **rowPtr    = _el;
  for (unsigned i=_rows ; i != 0 ; i--) {
    Type *argColPtr = *argRowPtr++;
    Type *colPtr    = *rowPtr++;
    for (unsigned j=_cols ; j != 0 ; j--)
      if (*colPtr++ != *argColPtr++)
	return 1;
  }
  
  return 0;
}

template <class Type>
Mat<Type>& 
Mat<Type>::operator += (dcomplex x)
{
  Type addend = Type(real(x));

  Type *elPtr = _el[0];
  for (unsigned i=_rows ; i != 0 ; i--)
    for (unsigned j=_cols ; j != 0 ; j--, elPtr++)
      *elPtr += addend;
  
  return *this;
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>&
Mat<Type>::operator *= (dcomplex x)
{
  double scale = real(x);

  Type *elPtr = _el[0];
  for (unsigned i=_rows ; i != 0 ; i--)
    for (unsigned j=_cols ; j != 0 ; j--, elPtr++)
      *elPtr = Type(scale * *elPtr);
  
  return *this;
}

//
//-------------------------// 
//
template <class Type>
Mat<Type> 
Mat<Type>::t() const            // Transpose operation
{
   Mat<Type> Temp(_cols, _rows);
   Type **tempDataRowPtr = Temp._el;
   for (unsigned i=0 ; i < _cols ; i++) {
      Type *tempDataColPtr = *tempDataRowPtr++;
      for (unsigned j=0 ; j < _rows ; j++)
	 *tempDataColPtr++ = _el[j][i];
   }
   
   return Temp;
}

//************************
//reverse the filter or rotate by 180
//used for convolution and morphology 
//-------------------------// 
//
template <class Type>
Mat<Type> 
Mat<Type>::rotate180() const           
{
   Mat<Type> Temp(_rows,_cols);
   unsigned kr,kc;   
   for (unsigned i=0 ; i < _rows ; i++) {
      kr = _rows -i -1;
      for (unsigned j=0 ; j < _cols ; j++){
	 kc = _cols -j -1;
	 Temp(kr,kc) = _el[i][j];
      }
   }
   
   return Temp;
}

/**************************

inv - matrix inversion - Gauss-Jordan elimination with full pivoting.

The returned matrix is always a Type matrix.  Exits and prints
error message with singular matrices or bad MATRIX structures.

***************************/
template <class Type>
Mat<Type> 
Mat<Type>::inv() const
{
//Good results if used on float or double only
// check for square matrix 
  if(_rows != _cols) {
    cerr << endl << "Mat inversion ERROR: non-square, size = " 
	 << _rows << " x " << _cols << endl;
    return Mat<Type>();
  }
   
// check pointer 
  if(!_el) {
    cerr << endl << "ERROR: No matrix for matrix invert" << endl;
    return Mat<Type>();
  }
  
// size of square matrix 
  unsigned      n = _rows;
  Mat<Type>   Ai(*this);
  Type **a = Ai._el;
  
// allocate index arrays and set to zero 
  unsigned *pivot_flag = (unsigned *) calloc(n,sizeof(unsigned));
  unsigned *swap_row = (unsigned *) calloc(n,sizeof(unsigned));
  unsigned *swap_col = (unsigned *) calloc(n,sizeof(unsigned));
  
  if(!pivot_flag || !swap_row || !swap_col) {
    cerr << "Allocation error in matrix invert" << endl;
    return Mat<Type>();
  }
  
  for(unsigned i = 0 ; i < n ; i++) {   // n iterations of pivoting 
      
// find the biggest pivot element 
    unsigned irow,icol;
    unsigned row, col;
    Type big = (Type)0;
    for(row = 0 ; row < n ; row++) {
      if(!pivot_flag[row]) {       // only unused pivots 
	for(col = 0 ; col < n ; col++) {
	  if(!pivot_flag[col]) {
	    Type abs_element = (Type)fabs((double)a[row][col]);
	    if(abs_element >= big) {
	      big = abs_element;
	      irow = row;
	      icol = col;
	    }
	  }
	}
      }
    }
    pivot_flag[icol]++;    // mark this pivot as used 
    
    // swap rows to make this diagonal the biggest absolute pivot 
    Ai.swapRows(irow, icol);

    /* Old code for row swap
    if(irow != icol) {
      for(unsigned col = 0 ; col < n ; col++) {
	Type temp = a[irow][col];
	a[irow][col] = a[icol][col];
	a[icol][col] = temp;
      }
    }
    */
    
    // store what we swaped 
    swap_row[i] = irow;
    swap_col[i] = icol;
      
    // bad news if the pivot is zero 
    if(a[icol][icol] == (Type) 0) {
      cerr << "ERROR: matrix invert: SINGULAR MATRIX" << endl;
      return Mat<Type>();
    }
    
    // divide the row by the pivot 
    //note here if the opreration is integer results are bad
    Type pivot_inverse = (Type)1/a[icol][icol];
    a[icol][icol] = (Type) 1;        // pivot = 1 to avoid round off 
    for(col = 0 ; col < n ; col++)
      a[icol][col] = a[icol][col]*pivot_inverse;
    
    // fix the other _rows by subtracting 
    for(row = 0 ; row < n ; row++) {
      if(row != icol) {
	Type temp = a[row][icol];
	a[row][icol] = (Type) 0;
	for( col = 0 ; col < n ; col++)
	  a[row][col] = a[row][col]-a[icol][col]*temp;
      }
    }
  }
  
  // fix the affect of all the swaps for final answer 
  for(int swap = n-1 ; swap >= 0 ; swap--) {
    Ai.swapCols(swap_row[swap], swap_col[swap]);

    /* Old code for column swap
    if(swap_row[swap] != swap_col[swap]) {
      for(unsigned row = 0 ; row < n ; row++) {
	Type temp = a[row][swap_row[swap]];
	a[row][swap_row[swap]] = a[row][swap_col[swap]];
	a[row][swap_col[swap]] = temp;
      }
    }
    */
  }

     
  // free up all the index arrays 
  free((char *)pivot_flag);
  free((char *)swap_row);
  free((char *)swap_col);
  
  return Ai;
}

//
//-------------------------// 
//
template <class Type>
Mat<Type> 
Mat<Type>::h() const        // Hermitian Transpose operation
{
   Mat<Type> Temp(_cols, _rows);
   Type **tempDataRowPtr = Temp._el;
   for (unsigned i=0 ; i < _cols ; i++) {
      Type *tempDataColPtr = *tempDataRowPtr++;
      for (unsigned j=0 ; j < _rows ; j++)
	 *tempDataColPtr++ = _el[j][i];
   }
   return Temp;
}

//
//-------------------------// 
//
template <class Type>
Type 
Mat<Type>::min(unsigned *row, unsigned *col) const
{
   Type min = **_el;
   unsigned    rowOfMin = 0, colOfMin = 0;
   unsigned    rowindx = 0, colindx = 0;

   Type **rowPtr = _el;
   for(unsigned i = _rows ; i != 0 ; i--) {
      Type *colPtr = *rowPtr++;
      colindx = 0;
      for(unsigned j = _cols ; j != 0 ; j--){
	 if (*colPtr < min) {
	    min = *colPtr;
	    rowOfMin = rowindx;
	    colOfMin = colindx;
	 }
	 colPtr++;
         colindx++; 
      }
      rowindx++;
   }
   
   if (row)
      *row = rowOfMin;
   if (col)
      *col = colOfMin;
   
   return(min);
}

//
//-------------------------// 
//
template <class Type>
Type 
Mat<Type>::max(unsigned *row, unsigned *col) const
{
   Type max = **_el;
   unsigned    rowOfMax = 0, colOfMax = 0;
   unsigned    rowindx = 0, colindx = 0;
   
   Type **rowPtr = _el;
   for(unsigned i = _rows ; i != 0 ; i--) {
      Type *colPtr = *rowPtr++;
      colindx = 0;
      for(unsigned j = _cols ; j != 0 ; j--){
	 if (*colPtr > max) {
	    max = *colPtr;
	    rowOfMax = rowindx;
	    colOfMax = colindx;
	 }
	 colPtr++;
         colindx++;
      }
      rowindx++;
   }
   
   if (row)
      *row = rowOfMax;
   if (col)
      *col = colOfMax;
   
   return(max);
}

//
//-------------------------// 
//
template <class Type>
Type
Mat<Type>::median(Type minVal, Type maxVal) const
{
  return array(minVal, maxVal).medianVolatile();
}

//
//-------------------------// 
//
template <class Type>
dcomplex
Mat<Type>::csum() const
{
  dcomplex result = 0;
   
  Type **rowPtr = _el;
  for(unsigned i = _rows ; i != 0 ; i--) {
    Type *colPtr = *rowPtr++;
    for(unsigned j = _cols ; j != 0 ; j--)
      result += *colPtr++;
  }
  
  return result;
}

//
//-------------------------// 
//
template <class Type>
dcomplex
Mat<Type>::csum2() const
{
  dcomplex result = 0;
   
  Type **rowPtr = _el;
  for(unsigned i = _rows ; i != 0 ; i--) {
    Type *colPtr = *rowPtr++;
    for(unsigned j = _cols ; j != 0 ; j--, colPtr++)
      result += dcomplex(*colPtr * *colPtr);
  }
  
  return result;
}

//
//-------------------------// 
//
template <class Type>
dcomplex
Mat<Type>::ctrace() const
{
   unsigned mindim = (_rows >= _cols) ? _cols : _rows;
   
   Type **dataPtr = _el;
   dcomplex temp = 0;
   for (unsigned i=0 ; i < mindim ; i++)
      temp += dataPtr[i][i];
   
   return(temp);
}

//
//-------------------------// 
//
template <class Type>
dcomplex
Mat<Type>::cdet() const
{
  dcomplex determinant = 0;

  if (!_rows || (_rows != _cols)) {
    cerr << "Error: determinant of non-square or empty matrix" << endl;
    return determinant;
  }

  if (_cols > 1) {
    int factor = 1;
    Type *elPtr = *_el;
    for(unsigned col = 0; col < _cols; col++, factor *= -1)
      determinant += dcomplex(*elPtr++) * (residual(0, col).det() * factor);
  }
  else 
    determinant = **_el;

  return determinant; 
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>&
Mat<Type>::applyElementWise(double (*function)(double))
{
  Type *elPtr = _el[0];
  
  for(unsigned i = _rows; i; i--)
    for(unsigned j = _cols; j; j--, elPtr++)
      *elPtr = (Type) function((double) *elPtr);
  
  return(*this);
}

template <class Type>
Mat<Type>&
Mat<Type>::applyElementWiseC2D(double (*function)(dcomplex))
{
   Type *elPtr = _el[0];
   
   for(unsigned i = _rows; i; i--)
     for(unsigned j = _cols; j; j--, elPtr++)
       *elPtr = (Type) function((dcomplex) *elPtr);
   
   return(*this);
}

template <class Type>
Mat<Type>&
Mat<Type>::applyElementWiseC2C(dcomplex (*function)(dcomplex))
{
   Type *elPtr = _el[0];
   
   for(unsigned i = _rows; i; i--)
     for(unsigned j = _cols; j; j--, elPtr++)
       *elPtr = Type(real(function((dcomplex) *elPtr)));
   
   return(*this);
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>&
Mat<Type>::applyIndexFunction(IndexFunction F)
{
  Type *elPtr = *_el;
  for (unsigned r = 0; r < _rows; r++)
    for (unsigned c = 0; c < _rows; c++)
      *elPtr++ = Type(F(r, c));

  return(*this);
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>&
Mat<Type>::applyIndexFunction(ComplexIndexFunction F)
{
  Type *elPtr = *_el;
  for (unsigned r = 0; r < _rows; r++)
    for (unsigned c = 0; c < _rows; c++)
      *elPtr++ = Type(asDouble(F(r, c)));

  return(*this);
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>&
Mat<Type>::exp()
{
  return applyElementWise(::exp);
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>&
Mat<Type>::log()
{
  return applyElementWise(::log);
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>&
Mat<Type>::cos()
{
  return applyElementWise(::cos);
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>&
Mat<Type>::sin()
{
  Type *elPtr = _el[0];
  
  for(unsigned i = _rows; i; i--)
    for(unsigned j = _cols; j; j--, elPtr++)
      *elPtr = (Type) ::sin(double(*elPtr));
  
  return(*this);
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>&
Mat<Type>::round()
{
  return applyElementWise(::rint);
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>&
Mat<Type>::abs()
{
  return applyElementWise(::fabs);
}

//
//-------------------------// 
//

template <class Type>
Mat<Type>&
Mat<Type>::sqrt()
{
  return applyElementWise(::sqrt);
}

template <class Type>
Mat<Type>&
Mat<Type>::pow(double exponent)
{
  Type *elPtr = _el[0];
  for(unsigned i=_rows ; i != 0 ; i--)
    for(unsigned j=_cols ; j != 0 ; j--, elPtr++)
      *elPtr = (Type)(::pow((double)*elPtr, exponent));
  
  return *this;
}

template <class Type>
Mat<Type>&
Mat<Type>::conj()
{
  return *this;
}

//
//-------------------------//
//
template <class Type>
Mat<Type>
Mat<Type>::transposeXself() const
{
  Mat<Type>  result(_cols, _cols);
  Type     **resultEl = result._el;

  for (unsigned i = 0; i < _cols; i++) {
    Type *resultPtr = resultEl[i];
    
    for (unsigned j = 0; j < i; j++) {
      Type& value = *resultPtr++;
      value = Type(0);
      for (unsigned k = 0; k < _rows; k++)
	value += _el[k][i] * _el[k][j];

      resultEl[j][i] = value;
    }
    
    Type& value = *resultPtr++;
    value = Type(0);
    for (unsigned k = 0; k < _rows; k++)
      value += _el[k][i] * _el[k][i];
  }
  
  return result;
}

//
//-------------------------//
//
template <class Type>
void 
Mat<Type>::eig(Mat<Type>& D, Mat<Type>& V) const
{
   //Mat<double> Ad = (DblMat) (this);
   
   double  *d;
   unsigned np1;
   unsigned *positions, itemp;
   double  temp;
   
   // if (!Ad.getEl()) {
   if (!_el) {
      printf("eig: invalid input matrix pointer\n");
      exit(1);
   }
   
   if (_rows != _cols) {
      cerr << "eig: matrix is not square -- " << _rows << " x " << _cols << endl;
      exit(1);
   }
   
   np1 = _rows + 1;
   Mat<double> A(np1,np1);
   Mat<double> v(np1,np1);
   d = (double *) malloc((np1)*sizeof(double));
   unsigned i, j;
   for (i=1 ; i <= _rows ; i++)
      for (j=1 ; j <= _cols ; j++)
	 A(i,j) = (double)_el[i-1][j-1];
   //A(i,j) = Ad.getEl[i-1][j-1];
   
   jacobi((double **) A.getEl(), _rows, d,(double **) v.getEl());
   
// sort eigenvalues in descending order 
   
   positions = (unsigned *) malloc( (_rows + 1) * sizeof(unsigned));
   for (i=1 ; i <= _rows ; i++)
      positions[i] = i;
   for (i=1 ; i <= _rows ; i++)
      for (j=1 ; j < _rows ; j++){
	 if (d[j] < d[j+1]) {
	    temp = d[j]; d[j] = d[j+1]; d[j+1] = temp;
	    itemp = positions[j]; 
	    positions[j] = positions[j+1]; // Keep track of original positions 
	    positions[j+1] = itemp;
	 } 
      }
   
   for (i=0 ; i < _rows ; i++)
      for (j=0 ; j < _cols ; j++) {
	 V(i,j) = (Type)v(i+1,positions[j+1]);
	 D(i,j) = (Type)0.0;
      }
   for (i=0 ; i < _rows ; i++)
      D(i,i) = (Type)d[i+1];
   
   
   free((char *)d);
   free((char *)positions);
}

//
//-------------------------// 
//
template <class Type>
Mat<Type> 
Mat<Type>::house() const
{
   if (_cols != 1) {
      cerr << "Error: input to house is not a column vector." << endl;
      exit(1);
   }
   
   unsigned n = _rows;
   Type mu = (Type) 0.0;
   Type **rowPtr = _el;
   for(unsigned i=n ; i != 0 ; i--, rowPtr++) {
      mu += (**rowPtr) * (**rowPtr);
   }
   
   mu = (Type)(::sqrt((double)mu));
   Mat<Type> v(*this);
   if (mu != (Type) 0.0){
      Type xFirst = **_el;
      Type beta = (xFirst > (Type) 0.0) ? xFirst + mu : xFirst - mu;
      
      Type **rowPtr  = _el + 1;
      Type **vRowPtr = v._el + 1;
      for(unsigned i=n-1 ; i != 0 ; i--, rowPtr++, vRowPtr++){
	 **vRowPtr = **rowPtr / beta;
      }
   }
   
   **v._el = (Type) 1.0;
   
   return(v);
}

//
//-------------------------//Not sure 
//
template <class Type>
Mat<Type> 
Mat<Type>::rowhouse(const Mat<Type>& V) const
{
   Type beta;

   if (V._cols != 1) {
      cerr << "Error: input to rowhouse is not a column vector." << endl;
      exit(1);
   }
   
   if (V._rows != _rows) {
      cerr << "Error: vector input to rowhouse is wrong length." << endl;
      exit(1);
   }
   
   Type nv = (Type) 0.0;
   Type **vRowPtr = V._el;
   for (unsigned i=V._rows ; i != 0 ; i--, vRowPtr++)
      nv += (**vRowPtr) * (**vRowPtr);
   
   if (nv == (Type) 0.0){
      cerr << "Error: vector input to rowhouse is all Zeros." << endl;
      exit(1);
   }
   beta = (Type)(-2.0)/nv;
   Mat <Type> w(_cols,1);
   w =    (*this).t();
   w =    w * V;
   w =  w * beta;
   return ( *this + V * w.t() );
}


//
//-------------------------//
//
template <class Type>
Mat<Type>
Mat<Type>::convolv2d(const Mat<Type>& filter) const
{
  unsigned    row,column, i, j;
  unsigned    dead_rows, dead_cols, kernel_rows, kernel_cols, l_cols, l_rows;
  Type   *outptr, *inptr, *filterptr;

  Mat<Type> out(_rows,_cols);
  Mat<Type> filterrev = filter.rotate180();
  kernel_rows = filter.getrows();
  kernel_cols = filter.getcols();
  dead_rows = kernel_rows/2;
  dead_cols = kernel_cols/2;
  l_rows = 2*dead_rows + _rows;
  l_cols = 2*dead_cols +_cols;
  Mat<Type> large(l_rows, l_cols);
  Mat<Type> outlarge(l_rows, l_cols);

/*copy in into center of large*/
  for(i=0; i < _rows ; i++){
    for(j=0; j < _cols ; j++)
      large._el[i+dead_rows][j+dead_cols] = _el[i][j];
  }

  for(row=0; row < large._rows -kernel_rows+1; row++){
    outptr = outlarge._el[row+dead_rows];
    for (column=0; column<large._cols-kernel_cols+1; column++){
      Type sumval = 0;
      for (i=0; i< kernel_rows; i++){
	inptr = large._el[i+row];
	inptr = inptr + column;
	filterptr = (Type *) filterrev._el[i];
	for (j=0; j< kernel_cols; j++)
	  sumval += (*inptr++) * (*filterptr++);
      }

      outptr[column+dead_cols] = sumval;
    }
  }
  
  for(i=0; i < _rows ; i++){
    for(j=0; j < _cols ; j++)
      out._el[i][j] = outlarge._el[i+dead_rows][j+dead_cols];
  }

  return(out);
}

/***************/

template <class Type>
Histogram
Mat<Type>::histogram(double minin, double maxin, unsigned n) const
{
  if (maxin <= minin) {
    minin = min();
    maxin = max();
  }

  Histogram histogram(minin, maxin, n);

  const Type **rowPtr = (const Type **) _el;
  for (unsigned i = _rows; i; i--) {
    const Type *elPtr = *rowPtr++;
    for (unsigned j = _cols; j; j--)
      histogram.add(*elPtr++);
  }
  
  return histogram;
}

/*
//this histogram works in the following way
//if there is no input in the histogram ex: histogram(), then
//it takes the defaults 0,0 and does the histogram from the smallest
//value in the image to the largest value.
//if the input is histogram(0,255) it will give the histogram
//between 0 and 255.  If the input is 10,20 it will give only 10,20
//values

template <class Type>
Mat<Type>
Mat<Type>::histogram(int minin, int maxin) const
{
  if(minin == maxin){
    minin=(int)min();
    maxin=(int)max();
  }
  int range = maxin - minin + 1;
  Mat<Type> histogram1;
  histogram1 = Zeros<Type>(1, range);
  Type *histPtr = histogram1._el[0];
  Type **rowPtr = _el;
  for (int i = _rows; i != 0; i--) {
    Type *colPtr = *rowPtr++;
    for (int j = _cols; j != 0; j--, colPtr++) {
      if ((*colPtr >= minin) && (*colPtr <= maxin)) {
	int roundedEl = (int) rint((double)*colPtr);
	(*(histPtr + roundedEl))++;
      }
    }
  }
  
  return(histogram1);
}
*/

//
//-------------------------//
//
template <class Type>
Mat<Type>&
Mat<Type>::histmod(const Histogram& hist1, const Histogram& hist2)
{
  return map(equalize(hist1, hist2));
}

//
//-------------------------//
//
template <class Type>
SimpleArray<Type>
Mat<Type>::array(Type minVal, Type maxVal) const
{
  unsigned N = 0;
  Type   **rowPtr;
  Type    *elPtr;

  if (maxVal <= minVal) {
    minVal = min();
    maxVal = max();
    N = nElements();
  }
  else {
    // Scan matrix to determine # elements in the range
    rowPtr = _el;
    for (unsigned i = _rows; i; i--) {
      elPtr = *rowPtr++;
      for (unsigned j = _cols; j; j--, elPtr++)
	if ((*elPtr >= minVal) && (*elPtr <= maxVal))
	  N++;
    }
  }
  
  SimpleArray<Type> array(N);

  if (N) {
    Type *arrayPtr = array.contents();
    rowPtr = _el;
    for (unsigned i = _rows; i; i--) {
      elPtr = *rowPtr++;
      for (unsigned j = _cols; j; j--)
	if ((*elPtr >= minVal) && (*elPtr <= maxVal))
	  *arrayPtr++ = *elPtr++;
    }
  }

  return array;
}

//
//-------------------------//
//
template <class Type>
void 
Mat<Type>::qr(Mat<Type>& R, Mat<Type>& Q) const
{
   /*****************Don't need this because R and Q get reassigned later
   if ((_rows != R._rows) && (_cols != R._cols)) {
      cerr << "Error: R matrix of incompatible size for QR." << endl;
      exit(1);
   }
   if ((_rows != Q._rows) && (_rows != Q._cols)) {
      cerr << "Error: Q matrix of incompatible size for QR." << endl;
      exit(1);
   }
   ***********************************************************************/
   
   Mat<Type> Atmp(*this);
   //Mat<Type> vtmp(Zeros<Type>(_rows,1));
   //This vector vtmp should automatically be initialized to zero
   Mat<Type> vtmp(_rows,1);
   Mat<Type> QTemp(_rows,_rows);
//This below instantiation used to be at position 12345
   Mat<Type> Temp(*this);
   unsigned mindim = (_rows > _cols) ? _cols : _rows;
   
// This can be made more efficient -- Alex, 2/4/94

   unsigned i, k;
   int      j; 
   for(j=0 ; j < mindim ; j++){
      vtmp.resize(_rows-j,1);
      
      for(i=j ; i < _rows ; i++)
	 vtmp(i-j,0) = Atmp(i,j);
      
      vtmp = vtmp.house();
//Position 12345
      Temp(Temp.rowhouse(vtmp));
      
      for(i=j ; i < _rows ; i++)
	 for(k=j ; k < _cols ; k++)
	    Atmp(i,k) = Temp(i-j,k-j);
      
      if (j < (_rows-1)){
	 for(i=j+1 ; i < _rows ; i++)
	    Atmp(i,j) = vtmp(i-j,0);
      }
      
      Temp.resize(_rows-j-1,_cols-j-1);
      for(i=0 ; i < _rows-j-1 ; i++)
	 for(k=0 ; k < _cols-j-1 ; k++)
	    Temp(i,k) = Atmp(i+j+1,k+j+1);
   }
   
   R = Mat<Type>(_rows,_cols);
   for(i=0 ; i < mindim ; i++) 
      for(k=i ; k < mindim ; k++)
	 R(i,k) = Atmp(i,k);
   
   vtmp.resize();
   Q = Mat<Type>(_rows,_rows);
   Q.eye();
   for(j = mindim - 1 ; j >= 0 ; j--){
      vtmp.resize(_rows-j,1);
      QTemp.resize(_rows-j,_rows-j);
      vtmp(0,0) = (Type) 1;
      for(i=1 ; i < _rows-j ; i++)
	 vtmp(i,0) = Atmp(i+j,j);
      for(i=0 ; i < _rows-j ; i++)
	 for(k=0 ; k < _rows-j ; k++)
	    QTemp(i,k) = Q(i+j,k+j);
      QTemp = QTemp.rowhouse(vtmp);
      for( i=0 ; i < _rows-j ; i++)
	 for(k=0 ; k < _rows-j ; k++)
	    Q(i+j,k+j) = QTemp(i,k);
   }
   vtmp.resize();
   if (_rows < _cols){
      for(j=_rows ; j < _cols ; j++){
	 for(i=0 ; i < _rows ; i++){
	    vtmp(i,0) = _el[i][j];
	 }
	 vtmp = Q.t() * vtmp;
	 for(i=0 ; i < _rows ; i++){
	    R(i,j) = vtmp(i,0);
	 }
      }
   }
}

/*********************************************************
template <class Type>
void 
Mat<Type>::qrtwo(Mat<Type>& R, Mat<Type>& Q) const
{
   cerr << "Dcomplex 3 qr not implemented" << endl;
}
**********************************************************/
//
//-------------------------//
//
/*****************************************************************************
This qr function only works with the _rows greater than or equal to the number 
of columns.  It may be possible to modify this so that it works with any size
matrix.
******************************************************************************/
template <class Type>
void 
Mat<Type>::qr(Mat<Type>& R, Mat<Type>& Q, Mat<Type>& P) const
{
// This can be made more efficient -- Alex, 2/4/94 
   
   int i, j, k, r, *piv, mindim;
   unsigned rloc, cloc;
   Type tau, dtmp;
   Mat<Type> Atmp(*this);

/*********changed 6/17/94******************   
   P = Eye<Type>(_cols);
   R = Zeros<Type>(_rows,_cols);
*******************************************/
   if(_cols>_rows){
      cerr<<"_Rows must be greater than or equal to columns"<<endl;
      exit(1);
   } 
   P=Mat<Type>(_cols,_cols);
   P.eye();
   R=Mat<Type>(_rows,_cols);

   piv = (int *)malloc(_cols * sizeof(int));
   if (!piv){
      cerr << "Error forming piv integer vector in qr." << endl;
      cerr << "_rows = " << _rows << " _cols = " << _cols << endl;
      exit(1);
   }
   Mat<Type> c(_cols,1);
   Mat<Type> v(_rows,1);
   Mat<Type> tmp(_cols,1);
   Mat<Type> Temp(_rows,_cols);
   Mat<Type> QTemp(_rows,_rows);
   Mat<Type> PTemp(_cols,_cols);
   for(j=0 ; j < _cols ; j++){
      c(j) = (Type)0.0;
      piv[j] = 0;
      for(i=0 ; i < _rows ; i++)
	 c(j) = c(j) + Atmp(i,j)*Atmp(i,j);
   }
   r = -1;
   tau = c.max(&rloc, &cloc);
   k = rloc;

   
   while(tau > 0){
      r++;
      piv[r] = k;
      //Swap the rth and kth columns of A      
      for(i=0 ; i < _rows ; i++){
	 dtmp=Atmp(i,r); Atmp(i,r)=Atmp(i,k); Atmp(i,k)=dtmp;
      } 	 
      //Swap c(r) and c(k)      
      dtmp=c(r);   c(r)=c(k);     c(k)=dtmp;

        
/**********changed 6/17/94******************      
      PTemp = Eye<Type>(_cols);
********************************************/
//Obtain an identity matrix to perform swaps to calulate the Pth permutation
//matrix

      PTemp=Mat<Type> (_cols,_cols);
      PTemp.eye();

//swap row j=r and piv[j]=k of PTemp to create the jth permutation matrix
      PTemp(r,r)=(Type)0.0; 
      PTemp(k,k)=(Type)0.0; 
      PTemp(r,k)=(Type)1.0;
      PTemp(k,r)=(Type)1.0;
//Accumulate the Permutation matrix initially, P is an identity matrix
      P = P*PTemp;
//Resize Temp and v to prepare for the rth householder       
      Temp.resize(_rows-r,1);  
      v.resize(_rows-r,1);
//Temp gets elements A(r:m,r) 
      for(i=r ; i < _rows ; i++)
	 Temp(i-r,0) = Atmp(i,r);
//v gets house(A(r:m,r))
      v = Temp.house();
//Resize Temp to allow the rowhouse 
      Temp.resize(_rows-r,_cols-r);
//Temp=A(r:m,r:n)
      for(i=r ; i < _rows ; i++)
	 for(j=r ; j < _cols ; j++)
	    Temp(i-r,j-r) = Atmp(i,j);
//Temp=rowhouse(A(r:m,r:n),v(r:m))
      Temp = Temp.rowhouse(v);
//A(r:m,r:n)=Temp
      for(i=r ; i < _rows ; i++)
	 for(j=r ; j < _cols ; j++)
	    Atmp(i,j) = Temp(i-r,j-r);
//A(r+1:m,r)=v(r+1:n)
      for(i=r+1 ; i < _rows ; i++)
	 Atmp(i,r) = v(i-r);

      for(i=r+1 ; i < _cols ; i++)
	 c(i) = c(i) - Atmp(r,i)*Atmp(r,i);

      if ( r < (_cols -1))
      {
         /******************changed 6/21/94*************************
	 tmp.resize(_cols-1-r,1);
	 for(i=r+1 ; i < _cols ; i++)
	    tmp(i-r-1,0) = c(i);
	 tau = tmp.max(&rloc, &cloc);
         ***********************************************************/
         //tau gets max(c(r+1:n))
         tau = (c(r+1,_cols-1,0,0)).max(&rloc,&cloc); 
	 k = rloc + r + 1;
      }
      else
	 tau = (Type) 0;
   }
//end while
//R gets upper triangular part of A
   for(i=0 ; i < _rows ; i++)
      for(j=i ; j < _cols ; j++)
	 R(i,j) = Atmp(i,j);
   
   mindim = r+1;
   v.resize();
  
  /**************changed 6/17/94******************   
   Q = Eye<Type>(_rows);
   ***********************************************/
   Q=Mat<Type>(_rows,_rows);
   Q.eye();

   for( j = mindim - 1 ; j >= 0 ; j--){
      v.resize(_rows-j,1);
      QTemp.resize(_rows-j,_rows-j);
      v(0,0) = (Type) 1;
      for(i=1 ; i < _rows-j ; i++)
	 v(i,0) = Atmp(i+j,j);
      for(i=0 ; i < _rows-j ; i++)
	 for(k=0 ; k < _rows-j ; k++)
	    QTemp(i,k) = Q(i+j,k+j);
      QTemp = QTemp.rowhouse(v);
      for(i=0 ; i < _rows-j ; i++)
	 for(k=0 ; k < _rows-j ; k++)
	    Q(i+j,k+j) = QTemp(i,k);
   }
   v.resize();
   if (_rows < _cols){
      for(j=_rows ; j < _cols ; j++){
	 for(i=0 ; i < _rows ; i++){
	    v(i,0) = _el[i][j];
	 }
	 v = Q.t() * v;
	 for(i=0 ; i < _rows ; i++){
	    R(i,j) = v(i,0);
	 }
      }
   }
   free((char *)piv);
}

//
//-------------------------//
//
/*****************************************************************************
This svd function only works with the _rows greater than or equal to the number 
of columns.  It may be possible to modify this so that it works with any size
matrix.
******************************************************************************/

template <class Type>
void 
Mat<Type>::svd(Mat<Type>& U, Mat<Type>& S, Mat<Type>& V) const
{
// This can be made more efficient -- Alex, 2/4/94
/************************changed 6/27/94*************************************   
   if ((U._rows != _rows) || (U._cols != _rows)){
      cerr << "Error in svd: left eigenvector input matrix of wrong size." << endl;
      cerr << _rows << "x" << _rows << endl;
      exit(1);
   }
   if ((S._rows != _rows) || (S._cols != _cols)){
      cerr << "Error in svd: singular value input matrix of wrong size." << endl;
      cerr << _rows << "x" << _rows << endl;
      exit(1);
   }
   if ((V._rows != _cols) || (V._cols != _cols)){
      cerr << "Error in svd: right eigenvector input matrix of wrong size." << endl;
      cerr << _rows << "x" << _cols << endl;
      exit(1);
   }
******************************************************************************/
   
   if(_cols>_rows){
     cerr<<"_Rows must be greater than or equal to columns"<<endl;
     exit(1);
   } 
   
   U=Mat<Type> (_rows,_rows);
   S=Mat<Type> (_rows,_cols);
   V=Mat<Type> (_cols,_cols);
   
   Mat<Type> C(_cols,_cols);
   Mat<Type> D(_cols,_cols);
   Mat<Type> R(_rows,_cols);
   Mat<Type> P(_cols,_cols);
   
   C = t() * (*this);
   C.eig(D, V);
   
   ((*this)*V).qr(R,U,P);
   
   unsigned mindim = (_rows > _cols) ? _cols : _rows;
   
   S = U.t()*(*this)*V*P;
   unsigned i, j;
   for(i=0 ; i < mindim ; i++){
     Type& sii = S(i, i);
     if (double(sii) < 0.0) {
       sii = int(fabs(double(sii)));
       for(j=0 ; j < _rows ; j++)
	 U(j,i) = -U(j,i);
     }
   }
   V = V*P;
   
   for(i=0 ; i < _rows ; i++)
     for(j=0 ; j < _cols ; j++){
       if ( i != j)
	 S(i,j) = (Type)0.0;
     }
}

//
//-------------------------//
//
template < class Type>
Mat<Type>&
Mat<Type>::clip(Type minVal, Type maxVal, Type minFill, Type maxFill)
{
  Type **rowPtr = _el;
  for(unsigned i = _rows; i; i--) {
    Type *elPtr = *rowPtr++;
    for (unsigned j = _cols; j; j--, elPtr++) {
      if (*elPtr < minVal)
	*elPtr = minFill;
      if (*elPtr > maxVal)
	*elPtr = maxFill;
    }
  }
  
  return(*this);
}

//
//-------------------------//
//
template < class Type>
Mat<Type>&
Mat<Type>::map(const ValueMap& valueMap)
{
  Type **rowPtr = _el;
  for(unsigned i = _rows; i; i--) {
    Type *elPtr = *rowPtr++;
    for (unsigned j = _cols; j; j--, elPtr++)
      *elPtr = Type(valueMap(double(*elPtr)));
  }
  
  return(*this);
}

//
//-------------------------//
//
//The scale function will scale the image between 0 and 255
//unless given a range with a minimum and maximum value
//then it finds 256 levels and adds the minimum to all values
 
//note that: ((x-min)*(Max-Min)/(max-min))+Min
//          = x*((Max-Min)/(max-min)) + ((Minmax - Maxmin)/(max-min))

template < class Type>
Mat<Type>&
Mat<Type>::scale(double minout, double maxout, double minin, double maxin)
{
  if (minin >= maxin) {
    minin = asDouble(min());
    maxin = asDouble(max());
  }

  map(LinearMap(minin, maxin, minout, maxout));
  clip(Type(minout), Type(maxout), Type(minout), Type(maxout));

  return *this;
}

template <class Type>
void 
Mat<Type>::linearinterpolate(unsigned Urow, unsigned Drow, unsigned Ucol, unsigned Dcol, Mat<Type>& y) const
{

   unsigned     i, nrows, ncols, out_rows, out_cols;
   unsigned     n, m;
   
   nrows = _rows;
   ncols = _cols;
   
   out_cols = ((ncols-1)*Urow+1)/Drow;
   out_rows = ((nrows-1)*Ucol+1)/Dcol;

/**************changed 6/23/94**********************   
   remake(y, out_rows, out_cols);
   y.fill((Type)0);
   y = Zeros<Type>(out_rows, out_cols);
   
   yi.fill((Type)0);
   yi = Zeros<Type>(nrows, out_cols);
****************************************************/
   y=Mat<Type> (out_rows,out_cols);
   Mat<Type> yi(nrows,out_cols);
   
   Mat<Type> tempr(1,ncols*Urow+1);
   
   for(i=0 ; i < nrows ; i++){
      tempr.fill((Type)0);
/*************************
      tempr = Zeros<Type>(1,ncols*Urow + 1);
**************************/
      tempr(0) = _el[i][0];
      for(n=1 ; n < ncols ; n++){
	 tempr(n*Urow) = _el[i][n];
	 for(m=1 ; m < Urow ; m++){
	    tempr((n-1)*Urow + m) = 
	       tempr((n-1)*Urow) + (tempr(n*Urow)-tempr((n-1)*Urow))*m/Urow;
	 }
      }
      for(n=0 ; n < out_cols ; n++)
	 yi(i,n) = tempr(n*Drow);
   }
   
   Mat<Type> tempc(1,nrows*Ucol + 1);
   
   for(i=0 ; i < out_cols ; i++){
/***************changed 6/23/94****************   
      tempc = Zeros<Type>(1,nrows*Ucol + 1);
**********************************************/
      tempc.fill((Type)0);
      tempc(0) = yi(0,i);
      for(n=1 ; n < nrows ; n++){
	 tempc(n*Ucol) = yi(n, i);
	 for(m=1 ; m < Ucol ; m++){
	    tempc((n-1)*Ucol + m) = 
	       tempc((n-1)*Ucol) + (tempc(n*Ucol)-tempc((n-1)*Ucol))*m/Ucol;
	 }
      }
      for(n=0 ; n < out_rows ; n++)
	 y(n,i) = tempc(n*Dcol);
   }
}

//
//-------------------------//
//

template <class Type>
Mat<Type> 
Mat<Type>::filter(Mat<Type>& B,Mat<Type>& A) const
{
   unsigned     n, i, P, Q, N;
   Type  temp;
   
   if ((B._rows != 1) && (B._cols != 1)){
      cerr << "error in filter: numerator not a vector." << endl;
      cerr << B._rows << "x" << B._cols << endl;
      exit(1);
   }
   if ((A._rows != 1) && (A._cols != 1)){
      cerr << "error in filter: denominator not a vector." << endl;
      cerr << A._rows << "x" << A._cols << endl;
      exit(1);
   }
   if ((_rows != 1) && (_cols != 1)){
      cerr << "error in filter: input not a vector." << endl;
      cerr << _rows << "x" << _cols << endl;
      exit(1);
   }
   if (A(0) == 0){
      cerr << "error in filter: first element in the denominator" << endl;
      cerr << "      must be non-zero." << endl;
      exit(1);
   }
   P = A.length();
   P=P-1;   
   Q = B.length();   
   Q=Q-1;
   N = this->length();
   temp = A(0);
   if (A(0) != 1.0){
      for (i=0 ; i <= Q ; i++)
	 B(i) = B(i)/temp;
      for (i=0 ; i <= P ; i++)
	 A(i) = A(i)/temp;
   }
   
   Mat<Type> y(_rows,_cols);
   Mat<Type> pastx(1,Q+1);   
   Mat<Type> pasty(1,P+1);
   
   for(n=0 ; n < N ; n++){
      y(n) = 0;
      pastx(0) = this->operator () (n);
      for (i=1 ; i <= P ; i++)
	 y(n) = y(n) - A(i)*pasty(i);
      for (i=0 ; i <= Q ; i++)
	 y(n) = y(n) + B(i)*pastx(i);
      pasty(0) = y(n);
      for(i=P ; i > 0 ; i--)
	 pasty(i) = pasty(i-1);
      for(i=Q ; i > 0 ; i--)
	 pastx(i) = pastx(i-1);
   }
   
   return(y);
}

/***************************morphology functions*******************************/
#ifdef USE_DBLMAT
template <class Type>
Mat<Type>
Mat<Type>::erode(const Mat<double>& strel) const
{
  unsigned strelHeight = strel.getrows();
  unsigned strelWidth  = strel.getcols();

  if (((strelHeight == 1) && (strelWidth == 1)) || !strelHeight || !strelWidth)
    return Mat<Type>(*this);

  Mat<Type> padMatrix(pad(strelHeight/2, strelWidth/2));
  Mat<Type> result(_rows, _cols);

  Type   *padMatrixPtr1 = padMatrix._el[0];
  const double *strelStart = strel.getEl()[0];
  Type   *resultPtr     = result._el[0];

  unsigned padMatrixPtr2incr = padMatrix._cols - strelWidth;
  unsigned padMatrixPtr1incr = (strelWidth/2) * 2;
     
  for (unsigned y = _rows; y != 0; y--) {
    for (unsigned x = _cols; x != 0; x--) {
      Type   *padMatrixPtr2 = padMatrixPtr1;
      double  minimum = MAXDOUBLE;
      const double *strelPtr = strelStart;
      for (register unsigned j = strelHeight; j != 0; j--) {
	for (register unsigned i = strelWidth; i != 0; i--) {
	  if (*strelPtr >= 0)
	    minimum = MIN(minimum, (double) *padMatrixPtr2 - *strelPtr);
	  strelPtr++;
	  padMatrixPtr2++;
	}
	padMatrixPtr2 += padMatrixPtr2incr;
      }
      *resultPtr = Type(minimum);
      
      resultPtr++;
      padMatrixPtr1++;
    }
    padMatrixPtr1 += padMatrixPtr1incr;
  }
  
  return result;
}

template <class Type>
Mat<Type>
Mat<Type>::dilate(const Mat<double>& strel) const
{
   unsigned strelHeight = strel.getrows();
   unsigned strelWidth  = strel.getcols();

   if (((strelHeight == 1) && (strelWidth == 1)) || !strelHeight || !strelWidth)
      return Mat<Type>(*this);

   unsigned r=0;
   unsigned c=0;

   if ((strelHeight % 2) == 0) {
      strelHeight++;
      r++;
   }
   
   if ((strelWidth % 2) == 0) {
      strelWidth++;
      c++;
   }
   
   Mat<double> newStrel(strelHeight, strelWidth, -1);
   newStrel.insert(strel.rotate180(), r, c);
   
   Mat<Type> padMatrix(pad(strelHeight/2, strelWidth/2));
   Mat<Type> result(_rows, _cols);
   
   Type   *padMatrixPtr1 = padMatrix._el[0];
   const double *strelStart    = newStrel.getEl()[0];
   Type   *resultPtr     = result._el[0];
   
   unsigned padMatrixPtr2incr = padMatrix._cols - strelWidth;
   unsigned padMatrixPtr1incr = (strelWidth/2) * 2;
   
   for (unsigned y = _rows; y != 0; y--) {
      for (unsigned x = _cols; x != 0; x--) {
	 Type   *padMatrixPtr2 = padMatrixPtr1;
	 double  maximum = -MAXDOUBLE;
	 const double *strelPtr = strelStart;
	 for (register unsigned j = strelHeight; j != 0; j--) {
	    for (register unsigned i = strelWidth; i != 0; i--) {
	       if (*strelPtr >= 0)
		  maximum = MAX(maximum, (double) *padMatrixPtr2 + *strelPtr);
	       strelPtr++;
	       padMatrixPtr2++;
	    }
	    padMatrixPtr2 += padMatrixPtr2incr;
	 }
	 *resultPtr = Type(maximum);
	 
	 resultPtr++;
	 padMatrixPtr1++;
      }
      padMatrixPtr1 += padMatrixPtr1incr;
   }
   
   return result;
}

#endif // USE_DBLMAT

/******************************Old code***************************************
template <class Type>
Mat<Type> 
Mat<Type>::filter(Mat<Type>& B, Mat<Type>& A, Mat<Type>& x)
{
   unsigned     n, i, P, Q, N;
   Type  temp;
   
   if ((B._rows != 1) && (B._cols != 1)){
      cerr << "error in filter: numerator not a vector." << endl;
      cerr << B._rows << "x" << B._cols << endl;
      exit(1);
   }
   if ((A._rows != 1) && (A._cols != 1)){
      cerr << "error in filter: denominator not a vector." << endl;
      cerr << A._rows << "x" << A._cols << endl;
      exit(1);
   }
   if ((x._rows != 1) && (x._cols != 1)){
      cerr << "error in filter: input not a vector." << endl;
      cerr << x._rows << "x" << x._cols << endl;
      exit(1);
   }
   if (A(0) == 0){
      cerr << "error in filter: first element in the denominator" << endl;
      cerr << "      must be non-zero." << endl;
      exit(1);
   }
   P = length(A) - 1;   Q = length(B) - 1;   N = length(x);
   temp = A(0);
   if (A(0) != 1.0){
      for (i=0 ; i <= Q ; i++)
	 B(i) = B(i)/temp;
      for (i=0 ; i <= P ; i++)
	 A(i) = A(i)/temp;
   }
   
   Mat<Type> y(x._rows, x._cols);
   Mat<Type> pastx(1,Q+1);   
   Mat<Type> pasty(1,P+1);   
   pastx = Zeros<Type>(1,Q+1);
   pasty = Zeros<Type>(1,P+1);
   
   for(n=0 ; n < N ; n++){
      y(n) = 0;
      pastx(0) = x(n);
      for (i=1 ; i <= P ; i++)
	 y(n) = y(n) - A(i)*pasty(i);
      for (i=0 ; i <= Q ; i++)
	 y(n) = y(n) + B(i)*pastx(i);
      pasty(0) = y(n);
      for(i=P ; i > 0 ; i--)
	 pasty(i) = pasty(i-1);
      for(i=Q ; i > 0 ; i--)
	 pastx(i) = pastx(i-1);
   }
   
   return(y);
}
****************************************************************************/



//Need to be modified 
//
//-------------------------//
//

/****************************************************************
template <class Type>
void 
Mat<Type>::radialscale(Mat<Type>& scale, Mat<Type>& y) const
{  //To recover old code look in rmold.c
   unsigned    slen, out_rows, out_cols, nrows, ncols;
   unsigned    i, j;
   Type smax, smin, fPI;
   Type xrc, xcc, yrc, ycc;
   Type xang, yang, xrad, yrad, xr, xc;
   unsigned rloc,cloc;
   
   fPI = 4.0*atan(1.0);
   slen = 360;
   nrows = _rows;   ncols = _cols;
   smin = min(scale);     smax = max(scale);
   out_rows = (unsigned)(floor(smax*nrows + 0.5));
   out_cols = (unsigned)(floor(smax*ncols + 0.5));
   y=Mat<Type>(out_rows,out_cols);
   
   xrc = (Type)(nrows)/2.0;     xcc = (Type)(ncols)/2.0;
   yrc = (Type)(out_rows)/2.0;   ycc = (Type)(out_cols)/2.0;
   
   for(i=0 ; i < out_rows ; i++)
      for(j=0 ; j < out_cols ; j++){
	 yrad = (Type)(::sqrt((i-yrc)*(i-yrc) + (j-ycc)*(j-ycc)));
	 yang = (Type)floor(atan2((j-ycc) , (i-yrc))*180/fPI);
	 if (yang < (Type)0) 
	    yang = yang + (Type)360;
	 xrad = yrad/scale((int)(yang));
	 xang = yang;
	 xr   = (Type)floor(xrc + xrad*(Type)cos(xang*fPI/(Type)180) );
	 xc   = (Type)floor(xcc + xrad*(Type)sin(xang*fPI/(Type)180));
	 if ((xr<(Type)0) || (xr>=(Type)nrows) || (xc<(Type)0) || (xc>=(Type)ncols))
	    y(i,j) = (Type)0;
	 else
	    y(i,j) = _el[(unsigned)xr][(unsigned)xc];
	 
      }
}
******************************************************************************/

//
//-------------------------//
//
//
//------------------------/
//

#ifdef HAVE_MATLAB
/*****************************************************************************/
//This function goes to the named mat file and searches for the given variable 
//Next, if possible, it reinitializes the calling object to represent the
//same matrix as that defined by the Matlab variable.  The Matlab matrix
//is stored in double format.  The function Matlab2Mat will convert the 
//Matlab data to the type of its calling object (this).  This may mean a loss
//of precision.
template <class Type>
Boolean
Mat<Type>::loadMatlab(const char *filename, const char *varname)
{
  if (!filename || !varname) {
    cerr << "Mat<Type>::loadMatlab(): no filename or variable name specified" << endl;
    return FALSE;
  }

  //Check to ensure that the file can be opened successfully for reading
  MATFile *infile = matOpen((char *) filename,"r");
  if(infile==NULL) {
    cerr << "Mat<Type>::loadMatlab(): Can't open file: " << filename << endl;
    return FALSE;
  }
   
  //Assign pointer to the Matrix object allocated from the file and ensure that
  //the named variable was found
  Matrix *MatlabMatrix = matGetMatrix(infile, (char *) varname);
  if(MatlabMatrix==NULL) {
    cerr << "Mat<Type>::loadMatlab(): Could not find: " << varname << " in file: " 
	 << filename << endl;
    //close the file and quit 
    matClose(infile);
    return FALSE;
  }
  
  //Check to see if the Matlab Matrix was full or sparse, dcomplex or real
  //if (mxIsDcomplex(MatlabMatrix) == 0)
  if (mxIsComplex(MatlabMatrix) == 0)
    *this = matlab2Mat(MatlabMatrix);
  else {
    cerr<<"Mat<Type>::loadMatlab(): The variable: "<<varname<<" in the file: "<<filename;
    cerr<<" is dcomplex.  Create a dcomplex Mat object and call its";
    cerr<<" loadMatlab function."<<endl;
    matClose(infile);
    mxFreeMatrix(MatlabMatrix);
    return FALSE;
  }
  
  //We were successful
  matClose(infile);
  mxFreeMatrix(MatlabMatrix);

  return TRUE;
}

#endif // HAVE_MATLAB
//
//-------------------------//
//
template <class Type>
Boolean 
Mat<Type>::loadRaw(const char *filename, unsigned nrows, unsigned ncols)
{
  InputFile matrixFile(filename);

  if (!matrixFile) {
    cerr << "Error in loadRaw: error opening file." << endl;
    return FALSE;
  }

  _checkMatrixDimensions(filename, nrows, ncols);
  
  if ((nrows && (nrows != _rows)) || (_cols && (ncols != _cols))) {
    _rows = _maxrows = nrows;
    _cols = _maxcols = ncols;
    _allocateEl();
  }
  
  if (!matrixFile.stream().read((char *)_el + _maxrows*sizeof(Type *),
				_maxrows*_maxcols*sizeof(Type)))
    return FALSE;
  
  return TRUE;
}

//
//-------------------------//
//
template <class Type>
Boolean
Mat<Type>::loadAscii(const char *filename)
{
  InputFile infile(filename);
  if (!infile){
    cerr << "Error in loadAsccii: error opening file." << endl;
    return FALSE;
  }
  
  if (!(infile >> _rows >> _cols))
    return FALSE;
  
  _maxrows = _rows;
  _maxcols = _cols;
  
  _allocateEl();
  
  for(unsigned i=0; i < _rows; i++) {
    for(unsigned j=0; j < _cols ; j++)
      if (!(infile >> _el[i][j]))
	return FALSE;
  }
  
  return TRUE;
}

//
//-------------------------//
//
template <class Type>
Boolean
Mat<Type>::load(const char *filename, int type, const char *varname)
{
  Boolean status = FALSE;

  switch(type) {
  case RAW:
    status = loadRaw(filename);
    break;
  case ASCII:
    status = loadAscii(filename);
    break;
  case MATLAB:
#ifdef HAVE_MATLAB
    status = loadMatlab(filename, varname);
#else
    cerr << "Can't load MATLAB files. Please compile with -DHAVE_MATLAB" << endl;
#endif
    break;
  default:
    cerr << "Unrecognized type for loading" << endl;
    break;
  }

  return status;
}

//
//-------------------------//
//
template <class Type>
Boolean
Mat<Type>::save(const char *filename, int type, const char *varname) const
{
  Boolean result;		// bert - added

  switch(type) {
  case RAW:
    result = saveRaw(filename);
    break;
  case ASCII:
    result = saveAscii(filename);
    break;
  case MATLAB:
#ifdef HAVE_MATLAB
    result = saveMatlab(filename, varname);
#else
    cerr << "Can't save MATLAB files. Please compile with -DHAVE_MATLAB" << endl;
#endif
    break;
  default:
    cerr << "Unrecognized type for saving" << endl;
    result = 0;
    break;
  }
  return (result);		// bert - return a value.
}

//
//-------------------------//
//
#ifdef HAVE_MATLAB
/*****************************************************************************/
// Saves the calling object as a Matlab readable matrix to the named file with the
// given variable name. Since Matlab stores in double format, the calling object is 
// first converted to a temporary double array. The third argument lets the user
// update a file("u") or rewrite an existing file or create a new one("w")
template <class Type>
Boolean
Mat<Type>::saveMatlab(const char *filename, const char *varname, const char *option) const
{
  if (!_rows || !_cols) {
    cerr << "Attempt to save empty matrix as matlab file. Not saved" << endl;
    return FALSE;
  }

  //Create a temporary double array
  double *temp = 0;
  allocateArray(_rows*_cols, temp);
  if (!temp) {
    cerr << "Couldn't allocate temp array; matrix not saved" << endl;
    return FALSE;
  }

  // Fill temp array columnwise (this is Matlab's storage convention)
  double *tempPtr = temp;
  for (unsigned col = 0; col < _cols; col++)
    for (unsigned row = 0; row < _rows; row++)
      *tempPtr++ = _el[row][col];
  
  // Write the matrix out
  Boolean status = ::saveMatlab(filename, varname, option, _rows, _cols, temp);

  // Free temporary array
  freeArray(temp);

  return status;
}
#endif // HAVE_MATLAB

//
//-------------------------//
//
template <class Type>
Boolean
Mat<Type>::saveRaw(const char *filename) const
{
   ofstream outfile(filename);
   if (!outfile){
      cerr << "Error in saveRaw: error opening file." << endl;
      return FALSE;
   }
   
//outfile.write((unsigned char *) &nrows, sizeof(unsigned));
//outfile.write((unsigned char *) &ncols, sizeof(unsigned));
   
   outfile.write((char *)_el + _maxrows*sizeof(Type *) , _maxrows*_maxcols*sizeof(Type));
   
   outfile.close();

   return (outfile) ? TRUE : FALSE;
}

//
//-------------------------//
//
template <class Type>
Boolean
Mat<Type>::saveAscii(const char *filename) const
{
   unsigned   i, j;
   
   ofstream outfile(filename);
   if (!outfile){
      cerr << "Error in saveAscciifile: error opening file." << endl;
      return FALSE;
   }
   
   outfile << _rows << " " << _cols << endl;
   
   for(i=0; i < _rows; i++)
   {
      for(j=0; j < _cols ; j++)
	 outfile << _el[i][j] << " ";
      outfile << endl;
   }
   
   outfile.close();

   return (outfile) ? TRUE : FALSE;
}

//
//-------------------------//
//
 

//
//-------------------------//
//

//
// Private functions start
//
template <class Type>
void
Mat<Type>::_allocateEl()
{
  if (_el) {
#ifdef DEBUG
    cout << "Freeing allocated memory at " << _el << endl;
#endif
    if( _el[0] ) {
      delete [] _el[0];  // delete data block
    }
    delete [] _el; // delete row of pointers  
  }
  _el = 0;

  unsigned int nBytes = _maxrows*_maxcols*sizeof(Type);

  if (nBytes) {
    typedef Type * TypePtr;
    _el = (Type **) new TypePtr[_maxrows];
    assert(_el);
    _el[0] = (Type *) new Type[_maxrows*_maxcols];
    assert(_el[0]);
    memset(_el[0], 0, nBytes);  //set all elements to zero

#ifdef DEBUG
    cout << "Allocated " << nBytes << " bytes at " <<  _el << endl;
#endif
   
    for (unsigned i=1 ; i < _maxrows ; i++)
      _el[i] = (Type *)(&(_el[i-1][_maxcols]));
  }
}

/*************************************************************************
This function will actually copy another matrix to a subsection of the
calling matrix
Since it modifies data internally, it has been declared private
**************************************************************************/

template <class Type> 
void 
Mat<Type>::_modify_sub_section(unsigned r1,unsigned r2,unsigned c1,unsigned c2,const Mat<Type>& B)
{
  if ((r1 > r2) || (c1 > c2) || (r2 >= _rows) || (c2>= _cols)) { // (bert)
      cerr << "Error in cropting: improper row or column sizes." << endl;
      cerr << r1 << " to " << r2 << " and" << endl;
      cerr << c1 << " to " << c2 << endl;
      exit(1);
   }
//If the selected subsection and B object are not the same size, then done   
   if((B.getrows() != (r2-r1+1)) || (B.getcols() != (c2-c1+1)))
   {
     cerr<<"Error:Input Matrix and subsection selections don't argree in size";
     cerr<<endl;
   } 
   
   
   for (unsigned i=r1; i<=r2 ; i++)
     for (unsigned j=c1; j<=c2; j++)
	   _el[i][j] = B(i-r1,j-c1);
          
}

template <class Type>
void
Mat<Type>::_checkMatrixDimensions(const char *path, unsigned& nrows,
				  unsigned& ncols) const
{
  if (Path(path).hasCompressedExtension())
    return;

  struct stat buf;
  stat(path, &buf);
  inferDimensions(buf.st_size/sizeof(Type), nrows, ncols);
}

template <class Type>
Mat<Type>&
Mat<Type>::_fft(unsigned nrows, unsigned ncols, FFTFUNC fftFunc)
{
  using std::max;		// (bert) force it to use the right max()
  using std::abs;		// (bert) and force the right abs()!

  // Verify dimensions of FFT
  if ((nrows > 1) && ((nrows < _rows) || !isPowerOf2(nrows))) {
    cerr << "Warning! Mat<Type>::fft():" << endl
	 << "  Requested # _rows for FFT (" << nrows << ") invalid;" << endl;
    nrows = max(unsigned(4), nextPowerOf2(max(nrows, _rows)));
    cerr << "  increased to " << nrows << endl;
  }

  if ((ncols > 1) && ((ncols < _cols) || !isPowerOf2(ncols))) {
    cerr << "Warning! Mat<Type>::fft():" << endl
	 << "  Requested # _cols for FFT (" << ncols << ") invalid;" << endl;
    ncols = max(unsigned(4), nextPowerOf2(max(ncols, _cols)));
    cerr << "  increased to " << ncols << endl;
  }

  Boolean doX = (ncols != 1) && (_cols > 1) ? TRUE : FALSE;
  Boolean doY = (nrows != 1) && (_rows > 1) ? TRUE : FALSE;

  if (!nrows)
    nrows = max(unsigned(4), nextPowerOf2(_rows));
  else if (nrows == 1)
    nrows = _rows;

  if (!ncols)
    ncols = max(unsigned(4), nextPowerOf2(_cols));
  else if (ncols == 1)
    ncols = _cols;

  // Pad matrix to the final FFT dimensions
  pad(nrows, ncols, (nrows - _rows)/2, (ncols - _cols)/2, 0);

  // Allocate temporary arrays
  double *real = 0;
  double *imag = 0;
  unsigned maxDim = max(doX ? _cols : 1, doY ? _rows : 1);

  allocateArray(maxDim, real);
  allocateArray(maxDim, imag);
  assert(real && imag);

  unsigned row, col;

  // Take 1D FFT in X (row) direction
  if (doX)
    for(row = 0; row < _rows; row++) {
      // Fill temporary FFT array
      double  *realPtr   = real;
      double  *imagPtr   = imag;
      Type *sourcePtr = _el[row];
      for(col = _cols; col != 0; col--) {
	dcomplex value(*sourcePtr++);
	*realPtr++ = value.real();
	*imagPtr++ = value.imag();
      }

      // Calculate 1D FFT
      fftFunc(_cols, real, imag);

      // Put results back
      realPtr   = real;
      imagPtr   = imag;
      sourcePtr = _el[row];
      for(col = _cols; col != 0; col--)
	*sourcePtr++ = Type(abs(dcomplex(*realPtr++, *imagPtr++)));
    }

  // Take 1D FFT in Y (column) direction
  if (doY)
    for(col = 0; col < _cols; col++) {
      // Fill temporary FFT array
      double  *realPtr   = real;
      double  *imagPtr   = imag;
      Type    *sourcePtr = _el[0] + col;
      for(row = _rows; row != 0; row--) {
	dcomplex value(*sourcePtr);
	*realPtr++ = value.real();
	*imagPtr++ = value.imag();
	sourcePtr += _cols;
      }

      // Calculate 1D FFT
      fftFunc(_rows, real, imag);

      // Put results back
      realPtr   = real;
      imagPtr   = imag;
      sourcePtr = _el[0] + col;
      for(row = _rows; row != 0; row--) {
	*sourcePtr = Type(abs(dcomplex(*realPtr++, *imagPtr++)));
	sourcePtr += _cols;
      }
    }

  // Free temporary arrays
  freeArray(real);
  freeArray(imag);

  return *this;
}

template <class Type>
Mat<Type>&
Mat<Type>::ifft(unsigned nrows, unsigned ncols)
{
  _fft(nrows, ncols, ::ifft);
  unsigned factor = 1;
  if (nrows != 1) factor *= _rows;
  if (ncols != 1) factor *= _cols;
  
  return *this /= factor;
}

#ifdef HAVE_MATLAB
//This function converts a Matlab matrix object to a Mat object
//The first argument is a pointer to a Matlab matrix object
//This function checks the storage format of the Matlab object
//(Full or Sparse) and does the appropriate conversions to return a Mat
//object of the calling object's type.
template <class Type>
Mat<Type>
Mat<Type>::matlab2Mat(Matrix *MatlabMatrix) const
{  
   //Make sure that pointer is not NULL
   if(MatlabMatrix==NULL)
     {
       cerr<<"The input argument points to a non existent object(NULL)"<<endl;
       exit(1);
     }
   
   //See if the object is real or dcomplex
   //if(mxIsDcomplex(MatlabMatrix))
   if(mxIsComplex(MatlabMatrix))
     {
       cerr<<"The Matlab matrix object is dcomplex.  Create a dcomplex Mat ";
       cerr<<"object first then repeat the call"<<endl;
       exit(1);
     }
   
   unsigned Matlab_rows=mxGetM(MatlabMatrix);
   unsigned Matlab_cols=mxGetN(MatlabMatrix);
   Mat<Type> Temp(Matlab_rows,Matlab_cols);
   double *MatlabDATA=mxGetPr(MatlabMatrix);
   //See if the object is Full or Sparse and do corresponding action
   if(mxIsFull(MatlabMatrix))
     { //Since the Matlab storage representation of the object
       //is exactly the transpose of the Mat representation, The object
       //is initialized with its indices reversed.
       
       for(unsigned i=0;i<Matlab_cols;i++)
          for(unsigned j=0;j<Matlab_rows;j++)
            {
             Temp._el[j][i]= (Type) *MatlabDATA;
             MatlabDATA++;
            } 
       return(Temp);
        
     }
   else if(mxIsSparse(MatlabMatrix))
     {
       //This code was derived from code form the external interface guide for 
       // Matlab pg. 7 (Copyright Math Works)
       int *ir=mxGetIr(MatlabMatrix);
       int *jc=mxGetJc(MatlabMatrix);
       int temp_index;
       for(int j=0; j<Matlab_cols;j++)
        {
          temp_index=jc[j+1];
          for(int k=jc[j];k<temp_index;k++)
            Temp._el[ir[k]][j] = (Type) MatlabDATA[k];
        }
       
       return(Temp);
     }

   return(NULL);
}

template <class Type>
Matrix *
Mat<Type>::mat2Matlab(const char *name) const
{ 
  Matrix *CopyofThis;

  //Allocate the memory for the Matlab matrix and set its name
  CopyofThis=mxCreateFull(_rows,_cols,0);
  mxSetName(CopyofThis,name);

  double *MatlabData=mxGetPr(CopyofThis);
  for(unsigned i=0;i< _cols;i++)
    for(unsigned j=0;j<_rows;j++)
     {
      *MatlabData= (double) _el[j][i];
       MatlabData++;
     } 
  
  return(CopyofThis);
}

#endif // HAVE_MATLAB

//
// Private functions end
//

template <class Type>
ostream&
operator << (ostream& os, const Mat<Type>& A)
{
  return A.display(os);
}

// 
// Type conversions
//

#ifdef USE_DBLMAT
template <class Type> Mat<double> asDblMat(const Mat<Type>& A)
{
  Mat<double> cast(A.getrows(), A.getcols());
  double     *castPtr = (double *) cast.getEl()[0];
  const Type *aPtr    = A.getEl()[0];
  for (unsigned i=A.nElements(); i; i--)
    *castPtr++ = asDouble(*aPtr++);
  
  return cast;
}
#endif // USE_DBLMAT

//
//-------------------------// 
//

#ifdef USE_FLMAT
template <class Type> Mat<float> asFlMat(const Mat<Type>& A)
{
  Mat<float>  cast(A.getrows(), A.getcols());
  float      *castPtr = (float *) cast.getEl()[0];
  const Type *aPtr    = A.getEl()[0];
  for (unsigned i=A.nElements(); i; i--)
    *castPtr++ = (float) asDouble(*aPtr++);
  
  return cast;
}
#endif

//
//-------------------------// 
//

#ifdef USE_INTMAT
template <class Type> Mat<int> asIntMat(const Mat<Type>& A)
{
  Mat<int>    cast(A.getrows(), A.getcols());
  int        *castPtr = (int *) cast.getEl()[0];
  const Type *aPtr    = A.getEl()[0];
  for (unsigned i=A.nElements(); i; i--)
    *castPtr++  = (int) asDouble(*aPtr++);

  return cast;
}
#endif

//
//-------------------------// 
//

#ifdef USE_UINTMAT
template <class Type> Mat<unsigned int> asUIntMat(const Mat<Type>& A)
{
  Mat<unsigned int> cast(A.getrows(), A.getcols());
  unsigned        *castPtr = (unsigned *) cast.getEl()[0];
  const Type      *aPtr    = A.getEl()[0];
  for (unsigned i=A.nElements(); i; i--)
    *castPtr++  = (unsigned) asDouble(*aPtr++);
  
  return cast;
}
#endif

//
//-------------------------// 
//

#ifdef USE_SHMAT
template <class Type> Mat<short> asShMat(const Mat<Type>& A)
{
  Mat<short>  cast(A.getrows(), A.getcols());
  short      *castPtr = (short *) cast.getEl()[0];
  const Type *aPtr    = A.getEl()[0];
  for (unsigned i=A.nElements(); i; i--)
    *castPtr++  = (short) asDouble(*aPtr++);
  
  return cast;
}
#endif

//
//-------------------------// 
//

#ifdef USE_USHMAT
template <class Type> Mat<unsigned short> asUShMat(const Mat<Type>& A)
{
  Mat<unsigned short>  cast(A.getrows(), A.getcols());
  unsigned short      *castPtr = (unsigned short *) cast.getEl()[0];
  const Type          *aPtr    = A.getEl()[0];
  for (unsigned i=A.nElements(); i; i--)
    *castPtr++  = (unsigned short) asDouble(*aPtr++);
  
  return cast;
}
#endif

//
//-------------------------// 
//

#ifdef USE_CHRMAT
template <class Type> Mat<char> asChrMat(const Mat<Type>& A)
{
  Mat<char>   cast(A.getrows(), A.getcols());
  char       *castPtr = (char *) cast.getEl()[0];
  const Type *aPtr    = A.getEl()[0];
  for (unsigned i=A.nElements(); i; i--)
    *castPtr++  = (char) asDouble(*aPtr++);
   
  return cast;
}
#endif

//
//-------------------------// 
//

#ifdef USE_UCHRMAT
template <class Type> Mat<unsigned char> asUChrMat(const Mat<Type>& A)
{
  Mat<unsigned char> cast(A.getrows(), A.getcols());
  unsigned char     *castPtr = (unsigned char *) cast.getEl()[0];
  const Type        *aPtr    = A.getEl()[0];
  for (unsigned i=A.nElements(); i; i--)
    *castPtr++  = (unsigned char) asDouble(*aPtr++);
  
  return cast;
}
#endif

// (bert) added this declaration here to make the compiler happier.  This
// way it doesn't see the default parameters as being listed as part of
// the function definition.  I think this was causing a compilation problem
// Jason found.
//
template<class Type>
SimpleArray<Type> array(const Mat<Type>& A, Type minVal, Type maxVal)
{
    return A.array(minVal, maxVal); 
}

//////////////////////////////////////////////////////////////////////

#ifdef __GNUC__

/* For the GNU compiler, include the explicit specializations in this file,
 * so they can be seen during class instantiation.
 */
#include "MatrixSpec.cc"

#define _INSTANTIATE_MAT(Type) \
         template class Mat<Type>;                   \
         template Mat<Type>& pmultEquals(Mat<Type>&, const Mat<Type> &); \
         template Mat<Type>& pdivEquals(Mat<Type>&, const Mat<Type> &); \
         template Mat<Type> inv(const Mat<Type> &); \
         template Mat<Type> operator*<Type>(double, const Mat<Type> &); \
         template<> unsigned Mat<Type>::_rangeErrorCount = 25; \
         template<> Boolean  Mat<Type>::flushToDisk = FALSE;

_INSTANTIATE_MAT(int);
_INSTANTIATE_MAT(float);
_INSTANTIATE_MAT(double);
template class Zeros<double>;
template class Ones<double>;
template class Eye<double>;
template ostream & operator<<(ostream &, Mat<double> const &);
template SimpleArray<double> array(Mat<double> const &, double, double);
#ifdef USE_COMPMAT
_INSTANTIATE_MAT(dcomplex);
template Mat<dcomplex> fft(const Mat<double> &, unsigned, unsigned);
template Mat<dcomplex> fft(const Mat<dcomplex> &, unsigned, unsigned);
template Mat<dcomplex> ifft(const Mat<double> &, unsigned, unsigned);
template Mat<dcomplex> ifft(const Mat<dcomplex> &, unsigned, unsigned);
#endif // USE_COMPMAT
#else
#ifdef USE_COMPMAT
template Mat<dcomplex>& pmultEquals(Mat<dcomplex>&, const Mat<dcomplex> &);
template Mat<dcomplex>& pdivEquals(Mat<dcomplex>&, const Mat<dcomplex> &);
#endif // USE_COMPMAT
template Mat<double>& pmultEquals(Mat<double>&, const Mat<double> &);
template Mat<double>& pdivEquals(Mat<double>&, const Mat<double> &);
#endif // __GNUC__ not defined
