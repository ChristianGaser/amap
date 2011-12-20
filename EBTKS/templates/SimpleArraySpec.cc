/* From SimpleArray.cc */

#if HAVE_FINITE
#ifndef finite
extern "C" int finite(double);
#endif /* finite() not defined (as macro) */
#define FINITE(x) finite(x)
#elif HAVE_ISFINITE
#ifndef isfinite
extern "C" int isfinite(double);
#endif /* isfinite() not defined (as macro) */
#define FINITE(x) isfinite(x)
#else
#error "Neither finite() nor isfinite() is defined on your system"
#endif /* HAVE_ISFINITE */

#ifdef USE_COMPMAT
template <>
SimpleArray<dcomplex>&
SimpleArray<dcomplex>::prune()
{
  unsigned i, j;
  for (i = 0, j = 0; i < _size; i++) {
    dcomplex value = getElConst(i);
    if (FINITE(real(value)) && FINITE(imag(value))) {
      if (i != j)
	setEl(j, value);
      j++;
    }
  }

  newSize(j);

  return *this;
}
#endif /* USE_COMPMAT */

#ifdef USE_FCOMPMAT
template<>
SimpleArray<fcomplex>&
SimpleArray<fcomplex>::prune()
{
  for (unsigned i = 0, j = 0; i < _size; i++) {
    fcomplex value = getElConst(i);
    if (FINITE(real(value)) && FINITE(imag(value))) {
      if (i != j)
	setEl(j, value);
      j++;
    }
  }

  newSize(j);

  return *this;
}
#endif

#ifdef USE_COMPMAT
template <>
SimpleArray<dcomplex>
SimpleArray<dcomplex>::round(unsigned n) const 
{
  SimpleArray<dcomplex> result(_size);

  const dcomplex *sourcePtr = _contents;
  dcomplex       *destPtr   = result.contents();
  if (n) {
    unsigned factor = (unsigned) pow(10.0, double(n));
    for (unsigned i = _size; i; i--) {
      *destPtr++ = dcomplex(ROUND(factor*(real(*sourcePtr)))/factor,
			   ROUND(factor*(imag(*sourcePtr)))/factor);
      sourcePtr++;
    }
  }
  else
    for (unsigned i = _size; i; i--) {
      *destPtr++ = dcomplex(ROUND(real(*sourcePtr)),
			   ROUND(imag(*sourcePtr)));
      sourcePtr++;
    }
  
  return result;
}
#endif /* USE_COMPMAT */

#ifdef USE_FCOMPMAT
template <>
SimpleArray<fcomplex>
SimpleArray<fcomplex>::round(unsigned n) const 
{
  SimpleArray<fcomplex> result(_size);

  const fcomplex *sourcePtr = _contents;
  fcomplex       *destPtr   = result.contents();
  if (n) {
    unsigned factor = (unsigned) pow(10.0, double(n));
    for (unsigned i = _size; i; i--) {
      *destPtr++ = fcomplex(ROUND(factor*(real(*sourcePtr)))/factor,
			   ROUND(factor*(imag(*sourcePtr)))/factor);
      sourcePtr++;
    }
  }
  else
    for (unsigned i = _size; i; i--) {
      *destPtr++ = fcomplex(ROUND(real(*sourcePtr)),
			   ROUND(imag(*sourcePtr)));
      sourcePtr++;
    }
  
  return result;
}
#endif

#ifdef USE_COMPMAT
template <>
SimpleArray<dcomplex> 
operator ^ (double base, const SimpleArray<dcomplex>& array) {
  unsigned N = array.size();
  SimpleArray<dcomplex> result(N);
  const dcomplex *sourcePtr = array.contents();
  dcomplex *resultPtr = result.contents();
  for (unsigned i = N; i != 0; i--)
    *resultPtr++ = dcomplex(pow(base, *sourcePtr++));
  return result;
} 
#endif

#ifdef USE_FCOMPMAT
template <>
SimpleArray<fcomplex> 
operator ^ (double base, const SimpleArray<fcomplex>& array) {
  unsigned N = array.size();
  SimpleArray<fcomplex> result(N);
  const fcomplex *sourcePtr = array.contents();
  fcomplex *resultPtr = result.contents();
  for (unsigned i = N; i != 0; i--)
    *resultPtr++ = fcomplex(pow(base, *sourcePtr++));
  return result;
} 
#endif

// Explicit specializations of the "abs" function for unsigned types
// (Yes, they're just identity operations.)
//
template <>
SimpleArray<unsigned int> SimpleArray<unsigned int>::abs() const
{
  return *this;
}

template <>
SimpleArray<char> SimpleArray<char>::abs() const
{
  return *this;
}

template <>
SimpleArray<unsigned char> SimpleArray<unsigned char>::abs() const
{
  return *this;
}

template <>
SimpleArray<unsigned short> SimpleArray<unsigned short>::abs() const
{
  return *this;
}
