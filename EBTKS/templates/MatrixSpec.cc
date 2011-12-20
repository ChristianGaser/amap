/* MatrixSpec.cc - Explicit specializations for the "Mat" template class.
 *
 * This file is not intended to be compiled separately, but rather to
 * be included in an "appropriate" file.  What "appropriate" means
 * seems to depend upon the compiler you're using. 
 */
#include <complex.h>

using namespace std;


#ifdef USE_COMPMAT
template <>
Mat<dcomplex>& 
Mat<dcomplex>::operator += (dcomplex addend)
{
  dcomplex *elPtr = _el[0];
  for (unsigned i=_rows ; i != 0 ; i--)
    for (unsigned j=_cols ; j != 0 ; j--, elPtr++)
      *elPtr += addend;
  
  return *this;
}
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
template <>
Mat<fcomplex>& 
Mat<fcomplex>::operator += (dcomplex addend)
{
  fcomplex *elPtr = _el[0];
  for (unsigned i=_rows ; i != 0 ; i--)
    for (unsigned j=_cols ; j != 0 ; j--, elPtr++)
      *elPtr += addend;
  
  return *this;
}
#endif // USE_FCOMPMAT

#ifdef USE_COMPMAT
template <>
Mat<dcomplex>&
Mat<dcomplex>::operator *= (dcomplex scale)
{
  dcomplex *elPtr = _el[0];
  for (unsigned i=_rows ; i != 0 ; i--)
    for (unsigned j=_cols ; j != 0 ; j--, elPtr++)
      *elPtr *= scale;
  
  return *this;
}
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
template <>
Mat<fcomplex>&
Mat<fcomplex>::operator *= (dcomplex scale)
{
  fcomplex *elPtr = _el[0];
  for (unsigned i=_rows ; i != 0 ; i--)
    for (unsigned j=_cols ; j != 0 ; j--, elPtr++)
      *elPtr *= scale;
  
  return *this;
}
#endif // USE_FCOMPMAT

#ifdef USE_COMPMAT
template <>
Mat<dcomplex> 
Mat<dcomplex>::inv() const
{
  cerr << "Mat<dcomplex>::inv() called but not implemented" << endl;
  return Mat<dcomplex>(*this);
}
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
template <>
Mat<fcomplex> 
Mat<fcomplex>::inv() const
{
  cerr << "Mat<fcomplex>::inv() called but not implemented" << endl;
  return Mat<fcomplex>(*this);
}
#endif // USE_FCOMPMAT

#ifdef USE_COMPMAT
template <>
dcomplex
Mat<dcomplex>::min(unsigned *, unsigned *) const
{
  cerr << "Mat<dcomplex>::min() called but not implemented" << endl;
  return 0;
}
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
template <>
fcomplex
Mat<fcomplex>::min(unsigned *, unsigned *) const
{
  cerr << "Mat<fcomplex>::min() called but not implemented" << endl;
  return 0;
}
#endif // USE_FCOMPMAT

#ifdef USE_COMPMAT
template <>
dcomplex
Mat<dcomplex>::max(unsigned *, unsigned *) const
{
  cerr << "Mat<dcomplex>::max() called but not implemented" << endl;
  return 0;
}
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
template <>
fcomplex
Mat<fcomplex>::max(unsigned *, unsigned *) const
{
  cerr << "Mat<fcomplex>::max() called but not implemented" << endl;
  return 0;
}
#endif // USE_FCOMPMAT

#ifdef USE_COMPMAT
template <>
dcomplex
Mat<dcomplex>::median(dcomplex, dcomplex) const
{
  cerr << "Mat<dcomplex>::median() called but not implemented" << endl;
  return 0;
}
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
template <>
fcomplex
Mat<fcomplex>::median(fcomplex, fcomplex) const
{
  cerr << "Mat<fcomplex>::median() called but not implemented" << endl;
  return 0;
}
#endif // USE_FCOMPMAT

#ifdef USE_COMPMAT
template <>
Mat<dcomplex>&
Mat<dcomplex>::applyElementWise(double (*function)(double))
{
  dcomplex *elPtr = _el[0];
  
  for(unsigned i = _rows; i; i--)
    for(unsigned j = _cols; j; j--, elPtr++)
      *elPtr = (dcomplex) function(real(*elPtr));
  
  return(*this);
}
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
template <>
Mat<fcomplex>&
Mat<fcomplex>::applyElementWise(double (*function)(double))
{
  fcomplex *elPtr = _el[0];
  
  for(unsigned i = _rows; i; i--)
    for(unsigned j = _cols; j; j--, elPtr++)
      *elPtr = (fcomplex) function(real(*elPtr));
  
  return(*this);
}
#endif // USE_FCOMPMAT

#ifdef USE_COMPMAT
template <>
Mat<dcomplex>&
Mat<dcomplex>::applyElementWiseC2C(dcomplex (*function)(dcomplex))
{
  dcomplex *elPtr = _el[0];
  
  for(unsigned i = _rows; i; i--)
    for(unsigned j = _cols; j; j--, elPtr++)
      *elPtr = function(*elPtr);
  
  return(*this);
}
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
template <>
Mat<fcomplex>&
Mat<fcomplex>::applyElementWiseC2C(dcomplex (*function)(dcomplex))
{
  fcomplex *elPtr = _el[0];
  
  for(unsigned i = _rows; i; i--)
    for(unsigned j = _cols; j; j--, elPtr++)
      *elPtr = fcomplex(function(*elPtr));
  
  return(*this);
}
#endif // USE_FCOMPMAT

#ifdef USE_COMPMAT
template <>
Mat<dcomplex>&
Mat<dcomplex>::applyIndexFunction(ComplexIndexFunction F)
{
  dcomplex *elPtr = *_el;
  for (unsigned r = 0; r < _rows; r++)
    for (unsigned c = 0; c < _rows; c++)
      *elPtr++ = dcomplex(F(r, c));

  return(*this);
}
#endif // USE_COMPMAT

#ifdef USE_COMPMAT
template <>
Mat<dcomplex>&
Mat<dcomplex>::log()
{
  using std::log;
  dcomplex *elPtr = _el[0];
  
  for(unsigned i = _rows; i; i--)
    for(unsigned j = _cols; j; j--, elPtr++)
      *elPtr = log(*elPtr);
  
  return(*this);
}
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
template <>
Mat<fcomplex>&
Mat<fcomplex>::log()
{
  using std::log;
  fcomplex *elPtr = _el[0];
  
  for(unsigned i = _rows; i; i--)
    for(unsigned j = _cols; j; j--, elPtr++)
      *elPtr = log(*elPtr);
  
  return(*this);
}
#endif // USE_FCOMPMAT

#ifdef USE_COMPMAT
template <>
Mat<dcomplex>&
Mat<dcomplex>::cos()
{
  using std::cos;
  dcomplex *elPtr = _el[0];
  
  for(unsigned i = _rows; i; i--)
    for(unsigned j = _cols; j; j--, elPtr++)
      *elPtr = cos(*elPtr);
  
  return(*this);
}
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
template <>
Mat<fcomplex>&
Mat<fcomplex>::cos()
{
  using std::cos;
  fcomplex *elPtr = _el[0];
  
  for(unsigned i = _rows; i; i--)
    for(unsigned j = _cols; j; j--, elPtr++)
      *elPtr = cos(*elPtr);
  
  return(*this);
}
#endif // USE_FCOMPMAT

#ifdef USE_COMPMAT
template <>
Mat<dcomplex>&
Mat<dcomplex>::sin()
{
  using std::sin;
  dcomplex *elPtr = _el[0];
  
  for(unsigned i = _rows; i; i--)
    for(unsigned j = _cols; j; j--, elPtr++)
      *elPtr = sin(*elPtr);
  
  return(*this);
}
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
template <>
Mat<fcomplex>&
Mat<fcomplex>::sin()
{
  using std::cos;
  fcomplex *elPtr = _el[0];
  
  for(unsigned i = _rows; i; i--)
    for(unsigned j = _cols; j; j--, elPtr++)
      *elPtr = sin(*elPtr);
  
  return(*this);
}
#endif // USE_FCOMPMAT

#ifdef USE_COMPMAT
template <>
Mat<dcomplex>&
Mat<dcomplex>::round()
{
  dcomplex *elPtr = _el[0];
   
  for(unsigned i = _rows; i; i--)
    for(unsigned j = _cols; j; j--, elPtr++)
      *elPtr = dcomplex(rint(real(*elPtr)), rint(real(*elPtr)));
  
  return(*this);
}
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
template <>
Mat<fcomplex>&
Mat<fcomplex>::round()
{
  fcomplex *elPtr = _el[0];
   
  for(unsigned i = _rows; i; i--)
    for(unsigned j = _cols; j; j--, elPtr++)
      *elPtr = fcomplex(rint(real(*elPtr)), rint(real(*elPtr)));
  
  return(*this);
}
#endif // USE_FCOMPMAT

#ifdef USE_COMPMAT
template <>
Mat<dcomplex>&
Mat<dcomplex>::abs()
{
  dcomplex *elPtr = _el[0];
  
  for(unsigned i = _rows; i; i--)
    for(unsigned j = _cols; j; j--, elPtr++)
      *elPtr = dcomplex( std::abs(*elPtr) );
  
  return(*this);
}
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
template <>
Mat<fcomplex>&
Mat<fcomplex>::abs()
{
  fcomplex *elPtr = _el[0];
  
  for(unsigned i = _rows; i; i--)
    for(unsigned j = _cols; j; j--, elPtr++)
      *elPtr = (fcomplex) ::abs(*elPtr);
  
  return(*this);
}
#endif // USE_FCOMPMAT

#ifdef USE_COMPMAT
template <>
Mat<dcomplex>&
Mat<dcomplex>::pow(double)
{
  cerr << "Mat<dcomplex>::pow() called but not implemented" << endl;
  return *this;
}
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
template <>
Mat<fcomplex>&
Mat<fcomplex>::pow(double)
{
  cerr << "Mat<fcomplex>::pow() called but not implemented" << endl;
  return *this;
}
#endif // USE_FCOMPMAT

#ifdef USE_COMPMAT
template <>
Mat<dcomplex>&
Mat<dcomplex>::conj()
{
  dcomplex *elPtr = _el[0];
  
  for(unsigned i = _rows; i; i--)
    for(unsigned j = _cols; j; j--, elPtr++)
      *elPtr = std::conj(*elPtr);
  
  return(*this);
}
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
template <>
Mat<fcomplex>&
Mat<fcomplex>::conj()
{
  fcomplex *elPtr = _el[0];
  
  for(unsigned i = _rows; i; i--)
    for(unsigned j = _cols; j; j--, elPtr++)
      *elPtr = std::conj(*elPtr);
  
  return(*this);
}
#endif // USE_FCOMPMAT

#ifdef USE_COMPMAT
template <>
void 
Mat<dcomplex>::eig(Mat<dcomplex>&, Mat<dcomplex>&) const
{
  cerr << "Mat<dcomplex>::eig() called but not implemented" << endl;
}
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
template <>
void 
Mat<fcomplex>::eig(Mat<fcomplex>&, Mat<fcomplex>&) const
{
  cerr << "Mat<fcomplex>::eig() called but not implemented" << endl;
}
#endif // USE_FCOMPMAT

#ifdef USE_COMPMAT
template <>
Mat<dcomplex> 
Mat<dcomplex>::house() const
{
  cerr << "Mat<dcomplex>::house() called but not implemented" << endl;
  return Mat<dcomplex>(*this);
}
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
template <>
Mat<fcomplex> 
Mat<fcomplex>::house() const
{
  cerr << "Mat<fcomplex>::house() called but not implemented" << endl;
  return Mat<fcomplex>(*this);
}
#endif // USE_FCOMPMAT

#ifdef USE_COMPMAT
template <>
Histogram
Mat<dcomplex>::histogram(double, double, unsigned) const
{
  cerr << "Mat<dcomplex>::histogram() called but not implemented" << endl;
  return Histogram(0);
}
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
template <>
Histogram
Mat<fcomplex>::histogram(double, double, unsigned) const
{
  cerr << "Mat<fcomplex>::histogram() called but not implemented" << endl;
  return Histogram(0);
}
#endif // USE_FCOMPMAT

#ifdef USE_COMPMAT
template <>
void 
Mat<dcomplex>::qr(Mat<dcomplex>&, Mat<dcomplex>&, Mat<dcomplex>&) const
{
  cerr << "Mat<dcomplex>::qr() called but not implemented" << endl;
}
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
template <>
void 
Mat<fcomplex>::qr(Mat<fcomplex>&, Mat<fcomplex>&, Mat<fcomplex>&) const
{
  cerr << "Mat<fcomplex>::qr() called but not implemented" << endl;
}
#endif // USE_FCOMPMAT

#ifdef USE_COMPMAT
template <>
void 
Mat<dcomplex>::svd(Mat<dcomplex>&, Mat<dcomplex>&, Mat<dcomplex>&) const
{
  cerr << "Mat<dcomplex>::svd() called but not implemented" << endl;
}
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
template <>
void 
Mat<fcomplex>::svd(Mat<fcomplex>&, Mat<fcomplex>&, Mat<fcomplex>&) const
{
  cerr << "Mat<fcomplex>::svd() called but not implemented" << endl;
}
#endif // USE_FCOMPMAT

#ifdef USE_COMPMAT
template <>
Mat<dcomplex>&
Mat<dcomplex>::clip(dcomplex, dcomplex, dcomplex, dcomplex)
{
  cerr << "Mat<dcomplex>::clip called but not implemented" << endl;
  return *this;
}
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
template <>
Mat<fcomplex>&
Mat<fcomplex>::clip(fcomplex, fcomplex, fcomplex, fcomplex)
{
  cerr << "Mat<fcomplex>::clip called but not implemented" << endl;
  return *this;
}
#endif // USE_FCOMPMAT

#ifdef USE_COMPMAT
template <>
Mat<dcomplex>&
Mat<dcomplex>::map(const ValueMap&)
{
  cerr << "Mat<dcomplex>::map called but not implemented" << endl;
  return *this;
}
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
template <>
Mat<fcomplex>&
Mat<fcomplex>::map(const ValueMap&)
{
  cerr << "Mat<fcomplex>::map called but not implemented" << endl;
  return *this;
}
#endif // USE_FCOMPMAT

#ifdef USE_COMPMAT
template <>
void 
Mat<dcomplex>::linearinterpolate(unsigned Urow, unsigned Drow, unsigned Ucol, unsigned Dcol, Mat<dcomplex>& y) const
{
  cerr << "Mat<dcomplex>::linearinterpolate not implemented" << endl;
}
#endif // USE_COMPMAT

#ifdef USE_COMPMAT
template <>
Mat<dcomplex> 
Mat<dcomplex>::filter(Mat<dcomplex>& B,Mat<dcomplex>& A) const
{
  cerr << "Mat<dcomplex>::filter not implemented" << endl;
  return B;
}
#endif // USE_COMPMAT

#ifdef USE_DBLMAT
#ifdef USE_COMPMAT
template <>
Mat<dcomplex>
Mat<dcomplex>::erode(const Mat<double>&) const
{
  cerr << "Mat<dcomplex>::erode() called but not implemented" << endl;
  return Mat<dcomplex>(*this);
}
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
template <>
Mat<fcomplex>
Mat<fcomplex>::erode(const Mat<double>&) const
{
  cerr << "Mat<fcomplex>::erode() called but not implemented" << endl;
  return Mat<fcomplex>(*this);
}
#endif // USE_FCOMPMAT

#ifdef USE_COMPMAT
template <>
Mat<dcomplex>
Mat<dcomplex>::dilate(const Mat<double>&) const
{
  cerr << "Mat<dcomplex>::dilate() called but not implemented" << endl;
  return Mat<dcomplex>(*this);
}
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
template <>
Mat<fcomplex>
Mat<fcomplex>::dilate(const Mat<double>&) const
{
  cerr << "Mat<fcomplex>::dilate() called but not implemented" << endl;
  return Mat<fcomplex>(*this);
}
#endif // USE_FCOMPMAT

#endif // USE_DBLMAT

#ifdef HAVE_MATLAB
#ifdef USE_COMPMAT
Boolean
Mat<dcomplex>::loadMatlab(const char *filename, const char *varname)
{
  if (!filename || !varname) {
    cerr << "Mat<dcomplex>::loadMatlab(): no filename or variable name specified" << endl;
    return FALSE;
  }

   //Check to ensure that the file can be opened successfully for reading
   MATFile *infile = matOpen((char *) filename,"r");
   if(infile==NULL) {
     cerr << "Mat<dcomplex>::loadMatlab(): Can't open file: "<<filename<<endl;
     return FALSE;
   }
   
   //Assign pointer to the Matrix object allocated from the file and ensure that
   //the named variable was found
   Matrix *MatlabMatrix = matGetMatrix(infile, (char *) varname);
   if(MatlabMatrix==NULL) {
     cerr << "Mat<dcomplex>::loadMatlab(): Could not find: "
	  << varname << " in file: " << filename << endl;
     //close the file and quit 
     matClose(infile);
     return FALSE;
   }
   
   //Check to see if the Matlab Matrix was full or sparse, dcomplex or real
   //if (mxIsDcomplex(MatlabMatrix) == 0) {
   if (mxIsComplex(MatlabMatrix) == 0) {
     Mat<double> Temp(_rows,_cols);
     Temp.loadMatlab(filename,varname); 
   }
   else
     *this = matlab2Mat(MatlabMatrix);         
   
   //We were successful
   matClose(infile);
   mxFreeMatrix(MatlabMatrix);

   return TRUE;
}
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
Boolean
Mat<fcomplex>::loadMatlab(const char *, const char *)
{
  cerr << "Mat<fcomplex>::loadMatlab() called but not implemented" << endl;
  return(0);     
}
#endif // USE_FCOMPMAT
#endif // HAVE_MATLAB

#ifdef USE_COMPMAT
template <>
Boolean
Mat<dcomplex>::loadAscii(const char *)
{
  cerr << "Mat<dcomplex>::loadAscii() not implemented" << endl;
  return FALSE;
}
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
template <>
Boolean
Mat<fcomplex>::loadAscii(const char *)
{
  cerr << "Mat<fcomplex>::loadAscii() not implemented" << endl;
  return FALSE;
}
#endif // USE_FCOMPMAT

#ifdef HAVE_MATLAB
#ifdef USE_COMPMAT
template <>
Boolean
Mat<dcomplex>::saveMatlab(const char *filename, const char *varname, const char *option) const
{
  if (!_rows || !_cols) {
    cerr << "Attempt to save empty matrix as matlab file. Not saved" << endl;
    return FALSE;
  }

  //Create a temporary double arrays for the real and imaginary parts
  double *real = 0;
  allocateArray(_rows*_cols, real);
  if (!real) {
    cerr << "Couldn't allocate real array; matrix not saved" << endl;
    return FALSE;
  }

  double *imag = 0;
  allocateArray(_rows*_cols, imag);
  if (!imag) {
    cerr << "Couldn't allocate imag array; matrix not saved" << endl;
    freeArray(real);
    return FALSE;
  }

  // Fill real and imag arrays columnwise (this is Matlab's storage convention)
  double *realPtr = real;
  double *imagPtr = imag;
  for (unsigned col = 0; col < _cols; col++)
    for (unsigned row = 0; row < _rows; row++) {
      dcomplex value = _el[row][col];
      *realPtr++ = value.real();
      *imagPtr++ = value.imag();
    }
  
  // Write the matrix out
  Boolean status = ::saveMatlab(filename, varname, option, _rows, _cols, real, imag);

  // Free temporary arrays
  freeArray(real);
  freeArray(imag);

  return status;
}
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
template <>
Boolean
Mat<fcomplex>::saveMatlab(const char *filename, const char *varname, const char *option) const
{
  if (!_rows || !_cols) {
    cerr << "Attempt to save empty matrix as matlab file. Not saved" << endl;
    return FALSE;
  }

  //Create a temporary double arrays for the real and imaginary parts
  double *real = 0;
  allocateArray(_rows*_cols, real);
  if (!real) {
    cerr << "Couldn't allocate real array; matrix not saved" << endl;
    return FALSE;
  }

  double *imag = 0;
  allocateArray(_rows*_cols, imag);
  if (!imag) {
    cerr << "Couldn't allocate imag array; matrix not saved" << endl;
    freeArray(real);
    return FALSE;
  }

  // Fill real and imag arrays columnwise (this is Matlab's storage convention)
  double *realPtr = real;
  double *imagPtr = imag;
  for (unsigned col = 0; col < _cols; col++)
    for (unsigned row = 0; row < _rows; row++) {
      fcomplex value = _el[row][col];
      *realPtr++ = value.real();
      *imagPtr++ = value.imag();
    }
  
  // Write the matrix out
  Boolean status = ::saveMatlab(filename, varname, option, _rows, _cols, real, imag);

  // Free temporary arrays
  freeArray(real);
  freeArray(imag);

  return status;
}
#endif // USE_FCOMPMAT
#endif // HAVE_MATLAB

#ifdef USE_COMPMAT
template <>
Boolean
Mat<dcomplex>::saveAscii(const char *) const
{
  cerr << "Mat<dcomplex>::saveAscii() not implemented" << endl;
  return FALSE;
}
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
template <>
Boolean
Mat<fcomplex>::saveAscii(const char *) const
{
  cerr << "Mat<fcomplex>::saveAscii() not implemented" << endl;
  return FALSE;
}
#endif // USE_FCOMPMAT

#ifdef USE_COMPMAT
template <>
Mat<dcomplex>&
Mat<dcomplex>::_fft(unsigned nrows, unsigned ncols, FFTFUNC fftFunc)
{
  using std::max;		// (bert) force it to use the right max()

  // Verify dimensions of FFT
  if ((nrows > 1) && ((nrows < _rows) || !isPowerOf2(nrows))) {
    cerr << "Warning! Mat<dcomplex>::fft():" << endl
	 << "  Requested # _rows for FFT (" << nrows << ") invalid;" << endl;
    nrows = max(unsigned(4), nextPowerOf2(max(nrows, _rows)));
    cerr << "  increased to " << nrows << endl;
  }

  if ((ncols > 1) && ((ncols < _cols) || !isPowerOf2(ncols))) {
    cerr << "Warning! Mat<dcomplex>::fft():" << endl
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
      dcomplex *sourcePtr = _el[row];
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
	*sourcePtr++ = dcomplex(*realPtr++, *imagPtr++);
    }

  // Take 1D FFT in Y (column) direction
  if (doY)
    for(col = 0; col < _cols; col++) {
      // Fill temporary FFT array
      double  *realPtr   = real;
      double  *imagPtr   = imag;
      dcomplex *sourcePtr = _el[0] + col;
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
	*sourcePtr = dcomplex(*realPtr++, *imagPtr++);
	sourcePtr += _cols;
      }
    }

  // Free temporary arrays
  freeArray(real);
  freeArray(imag);

  return *this;
}
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
template <>
Mat<fcomplex>&
Mat<fcomplex>::_fft(unsigned nrows, unsigned ncols, FFTFUNC fftFunc)
{
  using std::max;		// (bert) force it to use the right max()

  // Verify dimensions of FFT
  if ((nrows > 1) && ((nrows < _rows) || !isPowerOf2(nrows))) {
    cerr << "Warning! Mat<fcomplex>::fft():" << endl
	 << "  Requested # _rows for FFT (" << nrows << ") invalid;" << endl;
    nrows = max(unsigned(4), nextPowerOf2(max(nrows, _rows)));
    cerr << "  increased to " << nrows << endl;
  }

  if ((ncols > 1) && ((ncols < _cols) || !isPowerOf2(ncols))) {
    cerr << "Warning! Mat<fcomplex>::fft():" << endl
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
  unsigned maxDim = ::max(doX ? _cols : 1, doY ? _rows : 1);
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
      fcomplex *sourcePtr = _el[row];
      for(col = _cols; col != 0; col--) {
	fcomplex value(*sourcePtr++);
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
	*sourcePtr++ = fcomplex(*realPtr++, *imagPtr++);
    }

  // Take 1D FFT in Y (column) direction
  if (doY)
    for(col = 0; col < _cols; col++) {
      // Fill temporary FFT array
      double  *realPtr   = real;
      double  *imagPtr   = imag;
      fcomplex *sourcePtr = _el[0] + col;
      for(row = _rows; row != 0; row--) {
	fcomplex value(*sourcePtr);
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
	*sourcePtr = fcomplex(*realPtr++, *imagPtr++);
	sourcePtr += _cols;
      }
    }

  // Free temporary arrays
  freeArray(real);
  freeArray(imag);

  return *this;
}
#endif // USE_FCOMPMAT

#ifdef HAVE_MATLAB
#ifdef USE_COMPMAT
template <>
Mat<dcomplex>
Mat<dcomplex>::matlab2Mat(Matrix *MatlabMatrix) const
{  
  //Make sure that pointer is not NULL
  if(MatlabMatrix==NULL)
     {
       cerr<<"The input argument points to a non existent object(NULL)"<<endl;
       exit(1);
     }
   
   unsigned Matlab_rows=mxGetM(MatlabMatrix);
   unsigned Matlab_cols=mxGetN(MatlabMatrix);
   Mat<dcomplex> Temp(Matlab_rows,Matlab_cols);
   double *MatlabRealDATA=mxGetPr(MatlabMatrix);
   double *MatlabImagDATA=mxGetPi(MatlabMatrix);
   //See if the object is Full or Sparse and do corresponding action
   if(mxIsFull(MatlabMatrix))
     { //Since the Matlab storage representation of the object
       //is exactly the transpose of the Mat representation, The object
       //is initialized with its indices reversed.
       
       for(unsigned i=0;i<Matlab_cols;i++)
          for(unsigned j=0;j<Matlab_rows;j++)
            {
             Temp._el[j][i]= dcomplex( *MatlabRealDATA, *MatlabImagDATA);
             MatlabRealDATA++;
             MatlabImagDATA++;
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
       for(unsigned j=0; j<Matlab_cols;j++)
        {
          temp_index=jc[j+1];
          for(int k=jc[j];k<temp_index;k++)
            Temp._el[ir[k]][j] = dcomplex(MatlabRealDATA[k],MatlabImagDATA[k]);
        }
       
       return(Temp);
     }

   return(Temp);
}
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
template <>
Mat<fcomplex>
Mat<fcomplex>::matlab2Mat(Matrix *) const
{  
  cerr << "Mat<fcomplex>::_Matlab2Mat() called but not implemented" << endl;
  return Mat<fcomplex>(*this);
}
#endif // USE_FCOMPMAT

#ifdef USE_COMPMAT
template <>
Matrix *
Mat<dcomplex>::mat2Matlab(const char *name) const
{ 
  Matrix *CopyofThis;

  //Allocate the memory for the Matlab matrix and set its name
  CopyofThis=mxCreateFull(_rows,_cols,1);
  mxSetName(CopyofThis,name);

  double *MatlabPrData=mxGetPr(CopyofThis);
  double *MatlabPiData=mxGetPi(CopyofThis);
  for(unsigned i=0;i< _cols;i++)
    for(unsigned j=0;j<_rows;j++)
     {
      *MatlabPrData= _el[j][i].real();
       MatlabPrData++;
      *MatlabPiData= _el[j][i].imag();
       MatlabPiData++; 
     } 
  
  return(CopyofThis);
} 
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
template <>
Matrix *
Mat<fcomplex>::mat2Matlab(const char *) const
{ 
  cerr << "Mat<fcomplex>::_Mat2Matlab() called but not implemented" << endl;
  return 0;
}
#endif // USE_FCOMPMAT
#endif // HAVE_MATLAB

#ifdef USE_COMPMAT
template <class Type>
Mat<dcomplex> 
fft(const Mat<Type>& A, unsigned nrows, unsigned ncols) 
{
  return asCompMat(A).fft(nrows, ncols); 
}

template <class Type>
Mat<dcomplex> 
ifft(const Mat<Type>& A, unsigned nrows, unsigned ncols) 
{
  return asCompMat(A).ifft(nrows, ncols); 
}

#ifdef USE_DBLMAT
Mat<double>
applyElementWiseC2D(const Mat<dcomplex>& A, double (*function)(const dcomplex&))
{
  unsigned nrows = A.getrows();
  unsigned ncols = A.getcols();

  Mat<double> T(nrows, ncols);
   
  double  *TelPtr = (double *) T.getEl()[0];
  dcomplex *AelPtr = (dcomplex *) A.getEl()[0];
  
  for(unsigned i = nrows; i; i--)
    for(unsigned j = ncols ; j != 0 ; j--)
      *TelPtr++ = function(*AelPtr++);

  return T;
}

Mat<double>
arg(const Mat<dcomplex>& A)
{
  unsigned nrows = A.getrows();
  unsigned ncols = A.getcols();

  Mat<double> T(nrows, ncols);
   
  double  *TelPtr = (double *) T.getEl()[0];
  dcomplex *AelPtr = (dcomplex *) A.getEl()[0];
  
  for(unsigned i = nrows; i; i--)
    for(unsigned j = ncols ; j != 0 ; j--)
      *TelPtr++ = arg(*AelPtr++);

  return T;
}

static double
real(const dcomplex& d)
{
    return d.real();
}

static double
imag(const dcomplex& d)
{
    return d.imag();
}

Mat<double>
real(const Mat<dcomplex>& A)
{
  return applyElementWiseC2D(A, &real);
}

Mat<double>
imag(const Mat<dcomplex>& A)
{
  return applyElementWiseC2D(A, &imag);
}
#endif // USE_DBLMAT
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
template <class Type>
Mat<fcomplex> 
ffft(const Mat<Type>& A, unsigned nrows, unsigned ncols) 
{
  return asFcompMat(A).fft(nrows, ncols); 
}

template <class Type>
Mat<fcomplex> 
fifft(const Mat<Type>& A, unsigned nrows, unsigned ncols) 
{
  return asFcompMat(A).ifft(nrows, ncols); 
}
#ifdef USE_FLMAT
Mat<float>
applyElementWiseC2D(const Mat<fcomplex>& A, double (*function)(const dcomplex&))
{
  unsigned nrows = A.getrows();
  unsigned ncols = A.getcols();

  Mat<float> T(nrows, ncols);
   
  float  *TelPtr = (float *) T.getEl()[0];
  fcomplex *AelPtr = (fcomplex *) A.getEl()[0];
  
  for(unsigned i = nrows; i; i--)
    for(unsigned j = ncols ; j != 0 ; j--)
      *TelPtr++ = function(*AelPtr++);

  return T;
}

Mat<float>
arg(const Mat<fcomplex>& A)
{
  unsigned nrows = A.getrows();
  unsigned ncols = A.getcols();

  Mat<float> T(nrows, ncols);
   
  float    *TelPtr = (float *) T.getEl()[0];
  fcomplex *AelPtr = (fcomplex *) A.getEl()[0];
  
  for(unsigned i = nrows; i; i--)
    for(unsigned j = ncols ; j != 0 ; j--)
      *TelPtr++ = arg(*AelPtr++);

  return T;
}

Mat<float>
real(const Mat<fcomplex>& A)
{
  return applyElementWiseC2D(A, &dcomplex::real);
}

Mat<float>
imag(const Mat<fcomplex>& A)
{
  return applyElementWiseC2D(A, &dcomplex::imag);
}
#endif // USE_FLMAT
#endif // USE_FCOMPMAT

