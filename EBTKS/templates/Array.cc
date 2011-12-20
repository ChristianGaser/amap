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
$RCSfile: Array.cc,v $
$Revision: 1.3 $
$Author: bert $
$Date: 2004-12-08 16:43:44 $
$State: Exp $
--------------------------------------------------------------------------*/
#include <config.h>
#include "Array.h"
#include "dcomplex.h"
#include <assert.h>
#include <iostream>		// (bert) changed from iostream.h
#include <math.h>
//#include <string>		// (bert) changed from string.h
using namespace std;		// (bert) added

/******************
 * Array base class
 ******************/

#ifndef __GNUC__
template <class Type> unsigned Array<Type>::_arrayCtr = 0;
template <class Type> Boolean  Array<Type>::_debug = FALSE;
template <class Type> unsigned Array<Type>::_rangeErrorCount = 25;
#endif

template <class Type>
unsigned
Array<Type>::debug(Boolean on) 
{
  unsigned N = _arrayCtr;

  if (_debug = on)
    _arrayCtr = 0;
    
  return N;
}

//
// constructors for Array class
//

template <class Type>
Array<Type>::Array (unsigned sz)
{
  _self = this;
  _size = _maxSize = sz;

  if (_size) {
    _contents = new Type[_size];
    assert(_contents);
  }
  else
    _contents = 0;

  if (_debug) {
    _arrayCtr++;
    cout << "C" << _arrayCtr << ":" << long(this) << ":" << _size << " " << flush;
  }
}

template <class Type>
Array<Type>::Array (const Type& value, unsigned sz)
{
  _self = this;
  _size = _maxSize = sz;

  if (_size) {
    _contents = new Type[_size];
    assert(_contents);
    clear(value);
  }
  else
    _contents = 0;

  if (_debug) {
    _arrayCtr++;
    cout << "C" << _arrayCtr << ":" << long(this) << ":" << _size << " " << flush;
  }
}

template <class Type>
Array<Type>::Array(const Type *initArray, unsigned nElements)
{
  _self = this;
  _size = _maxSize = nElements;

  if (_size) {
    _contents = new Type[_size];
    assert(_contents);
    memcpy(_contents, initArray, _size*sizeof(Type));
  }
  else
    _contents = 0;

  if (_debug) {
    _arrayCtr++;
    cout << "C" << _arrayCtr << ":" << long(this) << ":" << _size << " " << flush;
  }
}

//
// Copy constructor
//

template <class Type>
Array<Type>::Array (const Array<Type>& array)
{
  _self     = this;
  _size     = _maxSize = 0;
  _contents = 0;

  operator = (array);

  if (_debug) {
    _arrayCtr++;
    cout << "C" << _arrayCtr << ":" << long(this) << ":" << _size << " " << flush;
  }
}

//
// destructor
//

template <class Type>
Array<Type>::~Array ()
{
  if (_debug) {
    _arrayCtr--;
    cout << "D" << _arrayCtr << ":" << long(this) << ":" << _size << " " << flush;
  }

  destroy();
}

//
// Iterator functions
//
/*
template <class Type>
void
Array<Type>::resetIterator(unsigned i)
{
  if (!_size)
    return;

  assert(i < _size);
  _itIndex = i;

  return;
}
*/
template <class Type>
void
Array<Type>::resetIterator(unsigned i) const
{
  if (!_size)
    return;

  assert(i < _size);
  _self->_itIndex = i;

  return;
}

template <class Type> 
Array<Type>
Array<Type>::operator () (unsigned nElements) const
{
  if (nElements > _size) {
    if (_rangeErrorCount) {
      cerr << "Warning! Array::operator(" << nElements 
	   << ") called with on array of size " << _size << ". Value truncated!" << endl;
      _rangeErrorCount--;
    }
    nElements = _size;
  }
  
  Array<Type> subArray(nElements);

  Type *sourcePtr = _contents;
  Type *destPtr   = subArray._contents;
  for (register unsigned i = nElements; i != 0; i--)
    *destPtr++ = *sourcePtr++;

  return(subArray);
}

template <class Type> Array<Type>
Array<Type>::operator () (unsigned start, unsigned end) const
{
  unsigned n = end - start + 1;
  if (start + n > _size) {
    if (_rangeErrorCount) {
      cerr << "Warning! Array::operator(" << start << ", " << end
	   << ") called with on array of size " << _size << ". Truncated!" << endl;
      _rangeErrorCount--;
    }
    n = _size - start;
  }

  Array<Type> subArray(n);

  Type *sourcePtr = _contents + start;
  Type *destPtr   = subArray._contents;
  for (register unsigned i = n; i != 0; i--)
    *destPtr++ = *sourcePtr++;

  return(subArray);
}

template <class Type>
Boolean 
Array<Type>::operator ! () const 
{
  return Boolean(!_size); 
}

template <class Type>
Type&
Array<Type>::first() 
{ 
  return _contents[0]; 
}

template <class Type>
Type&
Array<Type>::last()
{ 
  return _contents[_size - 1]; 
}

template <class Type> 
unsigned 
Array<Type>::size() const
{ 
  return _size; 
}

template <class Type> 
const Type *
Array<Type>::contents() const 
{ 
  return (_size) ? _contents : 0; 
}

template <class Type> 
Type *
Array<Type>::contents()
{ 
  return (_size) ? _contents : 0; 
}

template <class Type> 
Array<Type>::operator const Type *() const 
{ 
  return (_size) ? _contents : 0; 
}

template <class Type> 
Array<Type>::operator Type *()
{ 
  return (_size) ? _contents : 0; 
}

template <class Type> 
Array<Type>& 
Array<Type>::operator = (const Array<Type>& array ) 
{
  if (this == &array) return *this;

  newSize(array.size());

  resetIterator();
  array.resetIterator();

  for (unsigned i = _size; i; i--)
    (*this)++ = array++;

  return *this;
}

template <class Type> Array<Type>& 
Array<Type>::absorb(Array<Type>& array) {
  if (this == &array) return *this;

  if (_contents)
    delete [] _contents;
  _size = _maxSize = array._size;
  _contents = array._contents;

  array._size = 0;
  array._contents = 0;

  return(*this);
}

template <class Type> Array<Type>& 
Array<Type>::operator () (const Type *newContents, unsigned size) 
{
  if (size > _maxSize) {
    if (_contents)
      delete [] _contents;
    _contents = new Type[_size = _maxSize = size];
    assert(_contents);
  }
  else
    _size = size;

  register const Type *sourcePtr = newContents;
  register Type *destPtr   = _contents;
  for (register unsigned i = _size; i != 0; i--)
    *destPtr++ = *sourcePtr++;

  return(*this);
}

template <class Type> 
Type *
Array<Type>::asCarray(Type *destPtr) const
{
  if (!_size)
    return 0;

  if (!destPtr)
    destPtr = new Type[_size];

  if (destPtr) {
    const Type *sourcePtr = _contents;
    for (register unsigned i = _size; i; i--)
      *destPtr++ = *sourcePtr++;
  }

  return destPtr;
}

template <class Type>
Array<Type>&  
Array<Type>::append(const Type value) 
{
  if (_maxSize<=_size)_grow(SIZE_INCREMENT); 
  _contents[_size++]=value; 
  return *this;
}

template <class Type> 
Array<Type>& 
Array<Type>::append(const Array<Type>& array)
{
  unsigned nToAdd = array._size;

  if (!nToAdd)
    return *this;

  unsigned oldSize = _size;

  newSize(_size + nToAdd);

  const Type *sourcePtr = array._contents;
  Type *destPtr         = _contents + oldSize;
  for (register unsigned i = nToAdd; i != 0; i--)
    *destPtr++ = *sourcePtr++;

  return *this;
}

template <class Type> 
Array<Type>& 
Array<Type>::insert(const Type& value, unsigned index)
{
  if (index > _size) {
    if (_rangeErrorCount) {
      cerr << "Warning! Attempt to insert element outside range of array" << endl;
      _rangeErrorCount--;
    }
    return *this;
  }
  
  if (index == _size)
    return append(value);
  
  if (_maxSize <= _size)
    _grow(SIZE_INCREMENT);

  register Type *sourcePtr = _contents + _size - 1;
  register Type *destPtr   = sourcePtr + 1;
  for (register unsigned i = _size - index; i != 0; i--)
    *destPtr-- = *sourcePtr--;
  *destPtr = value;

  _size++;

  return *this;
}

template <class Type> 
Array<Type>& 
Array<Type>::insert(const Array<Type>& array, unsigned index)
{
  if (array._size == 0)
    return *this;

  unsigned oldSize = _size;

  newSize(_size + array._size);

  unsigned i;
  Type *sourcePtr = _contents + oldSize - 1;
  Type *destPtr   = sourcePtr + array._size;
  for (i = oldSize - index; i != 0; i--)
    *destPtr-- = *sourcePtr--;

  sourcePtr = array._contents + array._size - 1;
  for (i = array._size; i != 0; i--)
    *destPtr-- = *sourcePtr--;

  return *this;
}

template <class Type> 
Array<Type>& 
Array<Type>::replace(const Array<Type>& array, unsigned index)
{
  if (array._size == 0)
    return *this;

  if (index + array._size > _size)
    newSize(index + array._size);

  register Type *sourcePtr = array._contents;
  register Type *destPtr   = _contents + index;
  for (register unsigned i = array._size; i != 0; i--)
    *destPtr++ = *sourcePtr++;

  return *this;
}

template <class Type>
Type
Array<Type>::remove(unsigned index)
{
  if (!_size) {
    if (_rangeErrorCount) {
      _rangeErrorCount--;
      cerr << "Warning! Attempt to remove element from empty array" << endl;
    }
      
    return _contents[0];
  }

  if (index >= _size)
    _rangeError(index);

  if (index == _size - 1) {
    _size--;
    return _contents[index];
  }

  Type value(_contents[index]);

  register Type *destPtr   = _contents + index;
  register Type *sourcePtr = destPtr + 1;
  for (unsigned i = _size - index - 1; i != 0; i--)
    *destPtr++ = *sourcePtr++;
  _size--;

  return value;
}

template <class Type>
Type
Array<Type>::removeLast()
{
  if (!_size) {
    if (_rangeErrorCount) {
      _rangeErrorCount--;
      cerr << "Warning! Attempt to remove element from empty array" << endl;
    }
    return _contents[0];
  }
    
  _size--;
  
  return _contents[_size];
}

template <class Type>
Array<Type>&
Array<Type>::reorder(const Array<unsigned>& indices)
{
  Array<Type>     temp(*this);
  Type           *elPtr  = _contents;
  const unsigned *indexPtr = indices.contents();

  unsigned n = min(indices.size(), _size);

  for (unsigned i = n; i; i--, elPtr++, indexPtr++)
    if (*indexPtr < _size)
      *elPtr = temp[*indexPtr];

  return *this;
}

template <class Type>
Array<Type>&
Array<Type>::shuffle()
{
  for (unsigned i = 0; i < _size; i++) {
    unsigned j = unsigned(drand48()*_size);
    if (i != j) {
      Type temp    = _contents[i];
      _contents[i] = _contents[j];
      _contents[j] = temp;
    }
  }

  return *this;
}

template <class Type>
void
Array<Type>::clear(const Type& value)
{
  resetIterator();
  for (unsigned i = _size; i; i--)
    (*this)++ = value;
}

template <class Type>
void
Array<Type>::newSize(unsigned size)
{
  if (size == _size)
    return;

  if (size <= _maxSize) {
    _size = size;
    return;
  }

  Type *newContents = new Type[size];
  assert(newContents);

  if (_size) {
    register Type *sourcePtr = _contents;
    register Type *destPtr   = newContents;
    for (register unsigned i = _size; i != 0; i--)
      *destPtr++ = *sourcePtr++;
  }

  if (_contents)
    delete [] _contents;

  _contents = newContents;
  _size = _maxSize = size;
}

template <class Type>
Array<Type>&
Array<Type>::destroy()
{
  if (_contents) {
    delete [] _contents; 
    _contents = 0;
  }

  _size = _maxSize = 0;

  return *this;
}

template <class Type>
Array<Type>&
Array<Type>::operator >> (unsigned n)
{
  if (_size) {
    n %= _size;
    Array<Type> temp(n);

    unsigned i;
    Type *sourcePtr = _contents + _size - 1;
    Type *destPtr   = temp._contents + n - 1;
    for (i = n; i; i--)
      *destPtr-- = *sourcePtr--;

    destPtr = _contents + _size - 1;
    for (i = _size - n; i; i--)
      *destPtr-- = *sourcePtr--;

    sourcePtr = temp._contents + n - 1;
    for (i = n; i; i--)
      *destPtr-- = *sourcePtr--;
  }
  return *this;
}

template <class Type>
Array<Type>&
Array<Type>::operator << (unsigned n)
{
  if (_size) {
    n %= _size;
    Array<Type> temp(n);

    unsigned i;
    Type *sourcePtr = _contents;
    Type *destPtr   = temp._contents;
    for (i = n; i; i--)
      *destPtr++ = *sourcePtr++;

    destPtr = _contents;
    for (i = _size - n; i; i--)
      *destPtr++ = *sourcePtr++;

    sourcePtr = temp._contents;
    for (i = n; i; i--)
      *destPtr++ = *sourcePtr++;
  }
  return *this;
}

template <class Type>
Array<Type>
Array<Type>::sample(unsigned maxN) const
{
  double step = double(_size - 1)/(maxN - 1);

  if (step <= 1.0)
    return Array<Type>(*this);

  Array<Type> result(maxN);
  
  Type *destPtr = result._contents;
  double offset = 0;
    
  for (unsigned i = maxN; i; i--, offset += step)
    *destPtr++ = *(_contents + unsigned(floor(offset)));

  return result;
}

template <class Type>
Array<Type>
Array<Type>::applyElementWise(Type (*function) (Type)) const
{
  Array<Type> result(_size);

  Type *sourcePtr = _contents;
  Type *destPtr   = result._contents;

  for (unsigned i = _size; i != 0; i--)
    *destPtr++ = function(*sourcePtr++);

  return result;
}

template <class Type> 
ostream&
Array<Type>::print(ostream& os) const
{
  cerr << "Array<Type>::print(): Cannot print an Array" << endl;
  return os;
}

//
// Private functions
//

template <class Type> 
void
Array<Type>::_grow(unsigned amount)
{
  unsigned size = _size;
  newSize(_maxSize + amount);
  _size = size;
}

template <class Type> 
void
Array<Type>::_rangeError(unsigned& index) const
{
  if (_rangeErrorCount) {
    _rangeErrorCount--;
    cerr << "Corrected: index " << index << " into array of size " << _size 
	 << " !" << endl;
  }

  index = size() - 1;
}

template <class Type> 
void
Array<Type>::_notImplementedError() const
{
  cerr << "Array function called but not implemented" << endl;
  assert(0);
}

template <class Type> 
unsigned 
size(const Array<Type>& array) 
{
  return array.size(); 
}

template <class Type> 
ostream& 
operator << (ostream& os, const Array<Type>& array) 
{ 
  return array.print(os); 
}

#ifdef __GNUC__
#define _INSTANTIATE_ARRAY(Type)                       \
         template class Array<Type>;                   \
         template<> unsigned Array<Type>::_arrayCtr = 0;          \
         template<> Boolean  Array<Type>::_debug = FALSE;         \
         template<> unsigned Array<Type>::_rangeErrorCount = 25;  \
         template<> unsigned size(const Array<Type> &);

_INSTANTIATE_ARRAY(char);
_INSTANTIATE_ARRAY(unsigned char);
_INSTANTIATE_ARRAY(short);
_INSTANTIATE_ARRAY(int);
_INSTANTIATE_ARRAY(unsigned int);
_INSTANTIATE_ARRAY(unsigned short);
_INSTANTIATE_ARRAY(float);
_INSTANTIATE_ARRAY(double);
_INSTANTIATE_ARRAY(dcomplex);
#include "Path.h"
_INSTANTIATE_ARRAY(Path);
#include "ValueMap.h"
_INSTANTIATE_ARRAY(LinearMap);
#include "SimpleArray.h"
_INSTANTIATE_ARRAY(SimpleArray<unsigned> );
_INSTANTIATE_ARRAY(SimpleArray<float> );
#endif
