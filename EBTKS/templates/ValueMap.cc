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
$RCSfile: ValueMap.cc,v $
$Revision: 1.4 $
$Author: stever $
$Date: 2003-11-17 04:07:52 $
$State: Exp $
--------------------------------------------------------------------------*/

#include "ValueMap.h"

using namespace std;


#ifdef __GNUC__
#define _INSTANTIATE_VALUEMAP(Type) \
  template class LUT<Type>; \
  template SimpleArray<Type>& map(SimpleArray<Type>&, const ValueMap&); \
  template SimpleArray<Type>  mapConst(const SimpleArray<Type>&, const ValueMap&); \
  template SimpleArray<Type>& map(SimpleArray<Type>&, const Array<LinearMap>&); \
  template SimpleArray<Type>  mapConst(const SimpleArray<Type>&,const Array<LinearMap>&);
#endif

#if 0				// moved to header file
/*************************
 * Abstract ValueMap class
 *************************/

ValueMap& 
ValueMap::operator () (const ValueMap& map)
{
  return concat(map); 
}

ostream&
operator << (ostream& os, const ValueMap& map)
{
  return map.print(os);
}
#endif /* 0 */
/******************
 * Linear map class
 ******************/

#if 0				// moved to header file
ValueMap& 
LinearMap::inv() 
{ 
  _factor = 1.0/_factor; 
  _offset = -_offset; 
  return *this; 
}

ValueMap& 
LinearMap::concat(const ValueMap&)
{
  cerr << "LinearMap::concat() called but not implemented" << endl;
  return *this;
}

ValueMap& 
LinearMap::concat(const LinearMap& map)
{
  _offset += _factor*map._offset; 
  _factor *= map._factor;
  return *this;
}

LinearMap&
LinearMap::operator () (double factor, double offset) 
{ 
  _factor = factor; 
  _offset = offset; 
  return *this; 
}

LinearMap& 
LinearMap::operator () (double sourceMin, double sourceMax, double destMin,
			double destMax) 
{
  _factor = (destMax - destMin)/(sourceMax - sourceMin);
  _offset = destMin - _factor*sourceMin;
  return *this; 
}

ostream&
LinearMap::print(ostream& os) const 
{
  os << "(" << _factor << ", " << _offset << ")"; 
  return os;
}
#endif /* 0 */

/********************
 * Lookup table class
 ********************/

template <class Type>
LUT<Type>::LUT(unsigned length)
  : _source(length),
    _dest(length)
{
  _source.newSize(0);
  _dest.newSize(0);
}

template <class Type>
LUT<Type>::LUT(const SimpleArray<Type>& source, const SimpleArray<Type>& dest)
  : _source(source),
    _dest(dest)
{
  if (size(_source) != size(_dest)) {
    cerr << "Warning! LUT created with different array sizes. Truncating..." << endl;
    unsigned newSize = min(size(_source), size(_dest));
    _source.newSize(newSize);
    _dest.newSize(newSize);
  }

  _sort();
}

template <class Type>
LUT<Type>::LUT(const LUT<Type>& map)
  : _source(map._source),
    _dest(map._dest)
{}

template <class Type>
LUT<Type>::~LUT()
{
  _source.destroy();
  _dest.destroy();
}

template <class Type>
LUT<Type>&
LUT<Type>::add(Type source, Type dest)
{
  unsigned len = size(_source);

  if (!len) {
    _source.append(source);
    _dest.append(dest);
    return *this;
  }

  const Type *sourcePtr = _source.contents();
  unsigned i = 0;
  while ((i < len) && (*sourcePtr < source)) {
    i++;
    sourcePtr++;
  }

  if (*sourcePtr != source) {
    _source.insert(source, i);
    _dest.insert(dest, i);
  }

  return *this;
}

template <class Type>
LUT<Type>&
LUT<Type>::operator = (const LUT<Type>& map)
{
  _source = map._source;
  _dest   = map._dest;

  return *this;
}

template <class Type>
ValueMap&
LUT<Type>::inv()
{
  SimpleArray<Type> temp(_source);
  _source = _dest;
  _dest   = temp;

  _sort();

  return *this;
}

template <class Type>
ValueMap& 
LUT<Type>::concat(const ValueMap&)
{
  cerr << "LUT<Type>::concat() called but not implemented" << endl;
  return *this;
}

template <class Type>
ValueMap&
LUT<Type>::concat(const LUT<Type>& map)
{
  Type *destPtr = _dest.contents();

  for (unsigned i = size(_source); i; i--, destPtr++)
    *destPtr = Type(map(double(*destPtr)));

  return *this;
}

template <class Type>
double
LUT<Type>::operator () (double sourceValue) const
{
  unsigned length = size(_source);

  const Type *sourcePtr = _source.contents();
  for (unsigned i = 0; i < length; i++, sourcePtr++)
    if (*sourcePtr >= sourceValue)
      if (!i || (double(*sourcePtr) - sourceValue) < (sourceValue - *(sourcePtr - 1)))
	return double(_dest[i]);
      else 
	return double(_dest[i - 1]);
  
  return double(_dest[length - 1]);
}

template <class Type>
double
LUT<Type>::reverse(double destValue) const
{
  unsigned length = size(_dest);

  if (!length)
    return 0;
  if (length < 2)
    return _source[(unsigned int)0];

  unsigned index = 0;
  const Type *destPtr = _dest.contents();
  double minDiff = fabs(destValue - *destPtr++);
  for (unsigned i = 1; i < length; i++) {
    double diff = fabs(destValue - *destPtr++);
    if (diff < minDiff) {
      index = i;
      minDiff = diff;
    }
  }
  
  return double(_source[index]);
}

template <class Type>
ostream&
LUT<Type>::print(ostream& os) const
{
  unsigned length = size(_source);
  for (unsigned i = 0; i < length; i++)
    os << _source[i] << " " << _dest[i] << endl;
  return os;
}

template <class Type>
void
LUT<Type>::_sort()
{
  SimpleArray<unsigned> indexArray(_source.qsortIndexAscending());
  _source.reorder(indexArray);
  _dest.reorder(indexArray);
}

/**************************
 * Mapping of other objects
 **************************/

template <class Type> 
SimpleArray<Type>&
map(SimpleArray<Type>& array, const ValueMap& map)
{
  Type    *elPtr = array.contents();
  unsigned N     = array.size();

  for (unsigned i = N; i; i--, elPtr++)
    *elPtr = Type(map(double(*elPtr)));

  return array;
}

template <class Type>
SimpleArray<Type> 
mapConst(const SimpleArray<Type>& array, const ValueMap& valueMap) 
{
  SimpleArray<Type> A(array); return map(A, valueMap); 
}

template <class Type> 
SimpleArray<Type>&
map(SimpleArray<Type>& array, const Array<LinearMap>& maps)
{
  unsigned N = array.size();

  if (maps.size() != N) {
    cerr << "SimpleArray::map(): bad size of maps array" << endl;
    return array;
  }

  Type            *elPtr  = array.contents();
  const LinearMap *mapPtr = maps.contents();

  for (unsigned i = N; i; i--, elPtr++)
    *elPtr = Type((*mapPtr++)(double(*elPtr)));

  return array;
}

template <class Type>
SimpleArray<Type>  
mapConst(const SimpleArray<Type>& array,const Array<LinearMap>& maps)
{
  SimpleArray<Type> A(array); return mapConst(A, maps); 
}

#ifdef __GNUC__
  _INSTANTIATE_VALUEMAP(double);
#endif
