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
$RCSfile: CachedArray.cc,v $
$Revision: 1.10 $
$Author: bert $
$Date: 2004-12-08 17:02:18 $
$State: Exp $
--------------------------------------------------------------------------*/
#include <config.h>
#include "CachedArray.h"
#include "Histogram.h"
#include "FileIO.h"
#include <assert.h>
#include <stdio.h>
#include <unistd.h>


using namespace std;


//
// CachedArray class
//

#ifndef __GNUC__
template <class Type> const unsigned CachedArray<Type>::_DEFAULT_BLOCK_SIZE = 32768;
template <class Type> const unsigned CachedArray<Type>::_DEFAULT_N_BLOCKS   = 2;
template <class Type> unsigned CachedArray<Type>::_rangeErrorCount = 25;
#endif // __GNUC__ not defined
//
// Constructors/destructor
//

template <class Type>
CachedArray<Type>::CachedArray(unsigned size, unsigned nBlocks, unsigned blockSize)
  : SimpleArray<Type>(0)
{
  _self = this;
  _head = 0;
  _blocks = 0;
  this->_size = 0;
  _blockSize = 0;
  _nBlocks = 0;
  _itBlock = 0;
  _itBlockPtr = 0;
  this->_itIndex = 0;

  _initialize(size, nBlocks, blockSize);
  _openStream();

  resetHitRate();
}

template <class Type>
CachedArray<Type>::CachedArray(const Type *init, unsigned size, unsigned nBlocks, 
			       unsigned blockSize) 
  : SimpleArray<Type>(0)
{
  _self = this;
  _head = 0;
  _blocks = 0;
  this->_size = 0;
  _blockSize = 0;
  _nBlocks = 0;
  _itBlock = 0;
  _itBlockPtr = 0;
  this->_itIndex = 0;
  
  _initialize(size, nBlocks, blockSize);
  _openStream();

  copyFromCarray(init, size);

  resetHitRate();
}

template <class Type>
CachedArray<Type>::CachedArray(const CachedArray<Type>& array)
  : SimpleArray<Type>(0)
{
  _self = this;
  _head = 0;
  _blocks = 0;
  this->_size = 0;
  _blockSize = 0;
  _nBlocks = 0;
  _itBlock = 0;
  _itBlockPtr = 0;
  this->_itIndex = 0;

  _initialize(array._size, array._nBlocks, array._blockSize);
  _openStream();

  Array<Type>::operator = (array);
}

template <class Type>
CachedArray<Type>::CachedArray(const SimpleArray<Type>& array, unsigned nBlocks, 
			       unsigned blockSize)
  : SimpleArray<Type>(0)
{
  _self = this;
  _head = 0;
  _blocks = 0;
  this->_size = 0;
  _blockSize = 0;
  _nBlocks = 0;
  _itBlock = 0;
  _itBlockPtr = 0;
  this->_itIndex = 0;

  _initialize(array.size(), nBlocks, blockSize);
  _openStream();

  Array<Type>::operator = (array);
}

template <class Type>
CachedArray<Type>::~CachedArray()
{
  _destroy();
}

//
// Binary I/O functions
//

template <class Type>
ostream&
CachedArray<Type>::saveBinary(ostream& os, unsigned n, unsigned start) const
{
  if (start >= this->_size) {
    if (this->_size)
      if (_rangeErrorCount) {
	_rangeErrorCount--;
	cerr << "CachedArray::saveBinary: start out of range" << endl; 
      }
    
    return os;
  }
  
  if (!n)
    n = this->_size - start;
  else if (start + n > this->_size) {
    n = this->_size - start;
    if (_rangeErrorCount) {
      _rangeErrorCount--;
      cerr << "CachedArray::saveBinary: n too large; truncated" << endl;
    }
  }

  unsigned startInBlock = start % _blockSize;
  for (unsigned b = start / _blockSize; b < _maxNblocks; b++) {
    const CacheBlock<Type> *block = _read(b);
    unsigned nInBlock = MIN(n, _blockSize - startInBlock);
    //cout << "_size:" << _size << " n:" << n << " b:" << b 
    //<< " startInBlock:" << startInBlock << " nInBlock:" << nInBlock << endl;
    os.write((char *) block->_contents + startInBlock, nInBlock*sizeof(Type));
    startInBlock = 0;
    n -= nInBlock;
  }

  _revertIterator();

  return os;
}

template <class Type>
istream&
CachedArray<Type>::loadBinary(istream& is, unsigned n, unsigned start)
{
  if (!n)
    n = this->_size;

  newSize(start + n);

  if (this->_size) {
    unsigned startInBlock = start % _blockSize;
    for (unsigned b = start / _blockSize; b < _maxNblocks; b++) {
      CacheBlock<Type> *block = _read(b);
      unsigned nInBlock = MIN(n, _blockSize - startInBlock);
      
      is.read((char *) block->_contents + startInBlock, nInBlock*sizeof(Type));
      block->_changed = TRUE;
      startInBlock = 0;
      n -= nInBlock;
    }
  }

  _revertIterator();

  return is;
}

//
// Iterator functions
//
/*
template <class Type>
void
CachedArray<Type>::resetIterator(unsigned i)
{
  if (!this->_size)
    return;

  assert(i < this->_size);
  _self->_itBlock    = i / _blockSize;
  _self->_itBlockPtr = _read(_itBlock)->_contents;
  _self->_itIndex    = i - _itBlock * _blockSize;
  _self->_blocks[_itBlock]->_changed = TRUE;
}
*/
template <class Type>
void
CachedArray<Type>::resetIterator(unsigned i) const
{
  if (!this->_size)
    return;

  _self->_itBlock    = i / _blockSize;
  _self->_itBlockPtr = _read(_self->_itBlock)->_contents;
  _self->_itIndex    = i - _itBlock * _blockSize;
  // This line should be removed, and only be present in the non-const version
  // iIt is added here because g++ doesn't appear to be able to resolve the
  // const and non-const versions of resetIterator,
  _self->_blocks[_itBlock]->_changed = TRUE;

  assert((i < this->_size) && 
	 (_itBlock < _maxNblocks) && 
	 _itBlockPtr && 
	 (this->_itIndex < _blockSize));

  if (this->_debug)
    cout << "CachedArray::resetIterator:" << endl
	 << "   i:" << i 
	 << " _itIndex:" << this->_itIndex
	 << " _nBlocks:" << _nBlocks
         << " _maxNblocks:" << _maxNblocks
	 << " _blockSize:" << _blockSize 
	 << " _itBlock:" << _itBlock << endl;
}

template <class Type>
Type&
CachedArray<Type>::current()
{
  if (_self->_itIndex >= _blockSize) {
    _self->_itBlockPtr = _read(++_itBlock)->_contents;
    _self->_itIndex = 0;
    _self->_blocks[_itBlock]->_changed = TRUE;
  }

  return *(_itBlockPtr + this->_itIndex);
}

template <class Type>
const Type&
CachedArray<Type>::current() const
{
  if (_self->_itIndex >= _blockSize) {
    _self->_itBlockPtr = _read(++_self->_itBlock)->_contents;
    _self->_itIndex = 0;
  }

  return *(_itBlockPtr + this->_itIndex);
}

//
// Ascending operators
//
template <class Type>
Type&
CachedArray<Type>::operator ++()
{
  if (++(_self->_itIndex) >= _blockSize) {
    _self->_itBlockPtr   = _read(++_itBlock)->_contents;
    _self->_itIndex = 0;
    _self->_blocks[_itBlock]->_changed = TRUE;
  }

  return *(_itBlockPtr + this->_itIndex);
}

template <class Type>
const Type&
CachedArray<Type>::operator ++() const
{
  if (++(_self->_itIndex) >= _blockSize) {
    _self->_itBlockPtr   = _read(++(_self->_itBlock))->_contents;
    _self->_itIndex = 0;
  }

  return *(_itBlockPtr + this->_itIndex);
}

template <class Type>
Type&
CachedArray<Type>::operator ++(int)
{
  if (this->_itIndex >= _blockSize) {
    _self->_itBlockPtr   = _read(++(_self->_itBlock))->_contents;
    _self->_itIndex = 0;
    _self->_blocks[_itBlock]->_changed = TRUE;
  }

  return *(_itBlockPtr + (_self->_itIndex)++);
}

template <class Type>
const Type&
CachedArray<Type>::operator ++(int) const
{
  if (this->_itIndex >= _blockSize) {
    _self->_itBlockPtr = _read(++(_self->_itBlock))->_contents;
    _self->_itIndex    = 0;
  }

  return *(_itBlockPtr + (_self->_itIndex)++);
}

//
// Descending iterators
//
template <class Type>
Type&
CachedArray<Type>::operator --()
{
  if (int(--(_self->_itIndex)) < 0) {
    _self->_itBlockPtr = _read(--_itBlock)->_contents;
    _self->_itIndex = _blockSize - 1;
    _self->_blocks[_itBlock]->_changed = TRUE;
  }

  return *(_itBlockPtr + this->_itIndex);
}

template <class Type>
const Type&
CachedArray<Type>::operator --() const
{
  if (int(--(_self->_itIndex)) < 0) {
    _self->_itBlockPtr = _read(--(_self->_itBlock))->_contents;
    _self->_itIndex = _blockSize - 1;
  }

  return *(_itBlockPtr + this->_itIndex);
}

template <class Type>
Type&
CachedArray<Type>::operator --(int)
{
  if (int(this->_itIndex) < 0) {
    _self->_itBlockPtr = _read(--(_self->_itBlock))->_contents;
    _self->_itIndex = _blockSize - 1;
    _self->_blocks[_itBlock]->_changed = TRUE;
  }

  return *(_itBlockPtr + (_self->_itIndex)--);
}

template <class Type>
const Type&
CachedArray<Type>::operator --(int) const
{
  if (int(this->_itIndex) < 0) {
    _self->_itBlockPtr = _read(--(_self->_itBlock))->_contents;
    _self->_itIndex = _blockSize - 1;
  }

  return *(_itBlockPtr + (_self->_itIndex)--);
}

//
// Copy operators
//

template <class Type>
CachedArray<Type>&
CachedArray<Type>::copyFromCarray(const Type *init, unsigned size)
{
  cerr << "CachedArray::copyFromCarray() not implemented" << endl;

  return *this;
}

//
// Set functions
//
/*
template <class Type>
CachedArray<Type>::operator SimpleArray<Type> () const
{
  SimpleArray<Type> array(this->_size);
  Type             *destPtr = array.contents();

  resetIterator();
  for (unsigned i = this->_size; i; i--)
    *destPtr++ = (*this)++;
  
  return array;
}
*/
template <class Type>
void
CachedArray<Type>::newSize(unsigned size)
{
  if (!size) {
    _destroy();
    return;
  }

  if (!this->_size) {
    _initialize(size, _nBlocks ? _nBlocks : _DEFAULT_N_BLOCKS, 
		_blockSize ? _blockSize : _DEFAULT_BLOCK_SIZE);
    _openStream();
    return;
  }

  _flush();

  delete [] _blocks;
  delete _head;
  this->_size = 0;

  _initialize(size, _nBlocks, _blockSize);
  
  // Create file at requested size
  _self->_s.seekg(_maxNblocks*_blockSize*sizeof(Type));
  _self->_s.put('\0');
  assert(_s.is_open());
  
  CacheBlock<Type> *block = _head;
  for (unsigned i = 0; block; i++, block = block->_next) {
    assert(block->read(_self->_s, i));
    if (this->_debug)
      cout << "<read block " << i << " at " << long(block) << ">" << flush;
  }
}

template <class Type>
double
CachedArray<Type>::hitRate() const
{
  double hits = _hits;

  CacheBlock<Type> *block = _head;
  while (block) {
    hits += block->_nRead;
    hits += block->_nWrite;
    block = block->_next;
  }

  return hits/(hits + _misses);
}

template <class Type>
void
CachedArray<Type>::qsortAscending()
{
  if (!this->_size)
    cerr << "Warning: qsort attempted on empty CachedArray" << endl;

  _qsort(0, this->_size - 1);
}

template <class Type> 
Type
CachedArray<Type>::median() const
{
  assert(this->_size);

  if (this->_size <= _blockSize) {
    if (!_blocks[0])
      _read(0);
    return (*_blocks[0])(this->_size).medianVolatile();
  }

  CachedArray<Type> array(*this);

  return array.medianVolatile();
}

template <class Type> 
Type
CachedArray<Type>::medianVolatile()
{
  assert(this->_size);

  if (this->_size <= _blockSize) {
    if (!_blocks[0])
      _read(0);
    return (*_blocks[0])(this->_size).medianVolatile();
  }
  return _histMedian();
}

//
// Arithmetic operations
//
/*
template <class Type> 
CachedArray<Type>&
CachedArray<Type>::operator += (Type value)
{
  resetIterator();
  for (unsigned i = _size; i; i--)
    (*this)++ += value;
  
  return *this;
}

template <class Type> 
CachedArray<Type>&
CachedArray<Type>::operator += (const CachedArray<Type>& array)
{
  assert(_size == array._size);

  resetIterator();
  array.resetIterator();
  for (unsigned i = _size; i; i--)
    (*this)++ += array++;

  return *this;
}

template <class Type> 
CachedArray<Type>&
CachedArray<Type>::operator -= (const CachedArray<Type>& array)
{
  assert(_size == array._size);

  resetIterator();
  array.resetIterator();
  for (unsigned i = _size; i; i--)
    (*this)++ -= array++;
  
  return *this;
}

template <class Type> 
CachedArray<Type>&
CachedArray<Type>::operator *= (double value)
{
  resetIterator();
  for (unsigned i = _size; i; i--)
    (*this)++ *= value;
  
  return *this;
}

template <class Type> 
CachedArray<Type>&
CachedArray<Type>::operator *= (const CachedArray<Type>& array)
{
  assert(_size == array._size);

  resetIterator();
  array.resetIterator();
  for (unsigned i = _size; i; i--)
    (*this)++ *= array++;

  return *this;
}

template <class Type> 
CachedArray<Type>&
CachedArray<Type>::operator /= (const CachedArray<Type>& array)
{
  assert(this->_size == array._size);

  resetIterator();
  array.resetIterator();
  for (unsigned i = this->_size; i; i--)
    (*this)++ /= array++;
  
  return *this;
}
*/

template <class Type> 
CachedArray<Type>
CachedArray<Type>::ln() const
{
  CachedArray<Type> result(this->_size);

  resetIterator();
  result.resetIterator();
  for (unsigned i = this->_size; i; i--)
    result++ = Type(::log(double((*this)++)));

  return result;
}

template <class Type> 
CachedArray<Type>
CachedArray<Type>::log() const
{
  CachedArray<Type> result(this->_size);

  resetIterator();
  result.resetIterator();
  for (unsigned i = this->_size; i; i--)
    result++ = Type(::log10(double((*this)++)));

  return result;
}

template <class Type> 
CachedArray<Type>
CachedArray<Type>::exp() const   
{ 
  return ::exp(1.0)^(*this); 
}

template <class Type> 
CachedArray<Type>
CachedArray<Type>::exp10() const 
{ 
  return 10^(*this); 
}

template <class Type>
CachedArray<Type>
CachedArray<Type>::sample(unsigned maxN) const
{
  double step = double(this->_size - 1)/(maxN - 1);

  if (step <= 1.0)
    return CachedArray<Type>(*this);

  CachedArray<Type> result(maxN);
  result.resetIterator();

  double offset = 0;
  for (unsigned i = 0; i < maxN; i++, offset += step)
    result++ = getElConst(unsigned(::floor(offset)));

  return result;
}

template <class Type>
CachedArray<Type>
CachedArray<Type>::applyElementWise(Type (*function) (Type)) const
{
  CachedArray<Type> result(this->_size);

  resetIterator();
  result.resetIterator();
  for (unsigned i = this->_size; i; i--)
    result++ = function((*this)++);

  return result;
}

template <class Type>
CachedArray<Type>
CachedArray<Type>::map(const ValueMap& map) const
{
  CachedArray<Type> result(this->_size);

  resetIterator();
  result.resetIterator();
  for (unsigned i = this->_size; i; i--)
    result++ = Type(map((*this)++));

  return result;
}

//
// Private functions
//

template <class Type>
void
CachedArray<Type>::_initialize(unsigned size, unsigned nBlocks, unsigned blockSize)
{
  if (size) {
    assert(nBlocks && blockSize);

    this->_size       = size;
    _blockSize  = blockSize;
    _nBlocks    = nBlocks;
    _maxNblocks = unsigned(::ceil(double(size)/blockSize));

    if (_maxNblocks < nBlocks)
      _nBlocks = _maxNblocks;

    typedef CacheBlock<Type> *CBptr;
    _blocks = new CBptr[_maxNblocks];
    assert(_blocks);
  
    unsigned i;
    for (i = 0; i < _maxNblocks; i++)
      _blocks[i] = 0;

    _head = _blocks[0] = new CacheBlock<Type>(0, blockSize);
    assert(_head);
    for (i = 1; i < _nBlocks; i++)
      assert(_blocks[i] = _blocks[i-1]->addBlock(i, blockSize));
  }

  if (this->_debug) {
    cout << endl << "Created blocks:" << endl;
    for (unsigned i = 0; i < _nBlocks; i++)
      cout << "  " << long(_blocks[i]) << endl;
  }
}

template <class Type>
void
CachedArray<Type>::_openStream()
{
  if (_s.is_open())
    _s.close();

  if (this->_size) {
    char path[256];
    get_temp_filename(path);

    // Have to specify 'trunc' for some versions of the libraries, 
    // otherwise the file may not be created.
    //
    _s.open(path, ios::in | ios::out | ios::trunc);

    // Unlink the file so it will be deleted automatically when
    // closed.
    //
    unlink(path);

    assert(_s.is_open());

    // Create file at requested size
    _s.seekg(_maxNblocks*_blockSize*sizeof(Type));
    _s.put('\0');
  }
}

template <class Type>
void
CachedArray<Type>::_destroy()
{
  if (this->_size) {
    delete _head;
    _head = 0;

    delete [] _blocks;
    _blocks = 0;

    _s.close();
    _blockSize = 0;
    _maxNblocks = 0;
    this->_size = 0;

    _itBlock = 0;
    _itBlockPtr = 0;
    this->_itIndex = 0;
  }

  resetHitRate();
}

template <class Type>
void
CachedArray<Type>::_flush() const
{
  CacheBlock<Type> *block = _head;
  while (block) {
    block->write(_self->_s);
    block = block->_next;
  }
}

template <class Type>
CacheBlock<Type> *
CachedArray<Type>::_read(unsigned block) const
{
  if (this->_debug)
    cout << "<request for block " << block << ">" << flush;

  if (_blocks[block])
    return _blocks[block];

  // Locate least accessed block, first by write-access, then by read-access
  CacheBlock<Type> *readBlock = _head;
  CacheBlock<Type> *curBlock  = _head->_next;

  if (this->_debug)
    cout << "(" << long(readBlock) << ",r:" << readBlock->_nRead 
	 << ",w:" << readBlock->_nWrite << ")" << flush;

  while (curBlock) {
    if (this->_debug)
      cout << "(" << long(curBlock) << ",r:" << curBlock->_nRead 
	   << ",w:" << curBlock->_nWrite << ")" << flush;

    if ((curBlock->_nWrite && (curBlock->_nWrite < readBlock->_nWrite)) ||
	(!curBlock->_nWrite && (curBlock->_nRead < readBlock->_nRead)))
      readBlock = curBlock;

    curBlock = curBlock->_next;
  }

  // Reset the read/write counters in all blocks
  curBlock = _head;
  while (curBlock) {
    _self->_hits += curBlock->_nRead + curBlock->_nWrite;
    curBlock->_nRead = curBlock->_nWrite = 0;
    curBlock = curBlock->_next;
  }

  _self->_misses++;

  _blocks[readBlock->_ID] = 0;
    
  assert(readBlock->read(_self->_s, block));

  _blocks[block] = readBlock;

  if (this->_debug)
    cout << "<read block " << block << " at " << long(readBlock) << ">" << flush;

  return readBlock;
}

//
// Median and quicksort support functions
//

template <class Type>
void
CachedArray<Type>::_qsort(int p, int r)
{
  if (p < r) {
    int q = _partition(p, r);
    _qsort(p, q);
    _qsort(q + 1, r);
  }
}

template <class Type> 
Type
CachedArray<Type>::_histMedian(unsigned nBelow, unsigned nAbove)
{
  assert(this->_size);
  if (this->_debug)
      cout << "Begin: " << nBelow << " : " << nAbove << endl;

  if (this->_size <= _blockSize) {
    unsigned fullSize = nBelow + this->_size + nAbove;
    if (fullSize % 2)
      return _randomizedSelect(0, this->_size - 1, (fullSize+1)/2 - nBelow);
    else
      return _randomizedSelect(0, this->_size - 1, fullSize/2 - nBelow);
  }

  Type floor, ceil;
  extrema(&floor, &ceil);

  if (this->_debug)
      cout << "Floor and Ceiling: " << floor << " : " << ceil << endl;

  if (floor == ceil)
    return floor;


  Histogram hist(floor, ceil, MAX(this->_size/100, 10));
  resetIterator();
  for (unsigned i = this->_size; i; i--)
    hist.add((*this)++);

  if (this->_debug) {
      //  cout << endl << "Contents: " << *this << endl;
      //  cout << "Hist: " << hist << endl;
      cout << "[" << nBelow << ", " << nAbove << "]" << endl;
  }

  unsigned bin;
  double histMedian = hist.median(&bin, nBelow, nAbove);

  if (this->_debug)
      cout << "(" << bin << " : " << hist[bin] << " : " << histMedian << ") " << flush;

  unsigned nBelow2, nAbove2;
  removeAllNotIn(Type(hist.binStart(bin)), Type(hist.binStart(bin + 1)), 
		 &nBelow2, &nAbove2);

  if (this->_debug)
      cout << "nBelow2 : nAbove2 " << nBelow2 << " : " << nAbove2 << endl;

  return _histMedian(nBelow + nBelow2, nAbove + nAbove2);
}

// Computes the i'th order statistic (the i'th smallest element) of 
// an array.  This is the routine that does the real work, and is
// analogous to quicksort.

template <class Type>
Type
CachedArray<Type>::_randomizedSelect(int p, int r, int i)
{
// If we're only looking at one element of the array, we must have
// found the i'th order statistic.

  if (p == r)
    return getElConst(p);
   
// Partition the array slice A[p..r].  This rearranges the array so that
// all elements of A[p..q] are less than all elements of A[q+1..r].
// (Nothing fancy here, this is just a slight variation on the standard
// partition done by quicksort.)

  int q = _randomizedPartition(p, r);
  int k = q - p + 1;

// Now that we have partitioned A, its i'th order statistic is either
// the i'th order statistic of the lower partition A[p..q], or the
// (i-k)'th order stat. of the upper partition A[q+1..r] -- so 
// recursively find it.

  if (i <= k)
    return _randomizedSelect(p, q, i);
  else
    return _randomizedSelect(q+1, r, i-k);
}

template <class Type>
int
CachedArray<Type>::_randomizedPartition(int p, int r)
{
// Compute a random number between p and r

  int i = (random() / (MAXINT / (r-p+1))) + p;
  // Changed because of problems locating random() on a Sun
  //int i = int(ROUND(drand48() * (r-p+1) + p));

// Swap elements p and i

  Type temp = getElConst(p);
  setEl(p, getElConst(i));
  setEl(i, temp);

  return _partition(p, r);
}

template <class Type>
int
CachedArray<Type>::_partition(int p, int r)
{
  Type x = getElConst(p);
  int i  = p-1;
  int j  = r+1;
  
  while (1) {
    do { j--; } while (getElConst(j) > x);
    do { i++; } while (getElConst(i) < x);
    
    if (i < j) {
      Type temp = getElConst(i);
      setEl(i, getElConst(j));
      setEl(j, temp);
    }
    else
      return j;
  }      
}

template <class Type>
CachedArray<Type> 
operator ^ (double base, const CachedArray<Type>& array) {
  unsigned N = array.size();

  CachedArray<Type> result(N);

  array.resetIterator();
  result.resetIterator();
  for (unsigned i = N; i; i--)
    result++ = Type(pow(base, double(array++)));

  return result;
} 

//
// CacheBlock functions
//

template <class Type>
CacheBlock<Type>::CacheBlock(unsigned ID, unsigned size)
  : SimpleArray<Type>(size)
{
  _next = 0;
  _nBytes = size*sizeof(Type);
  _changed = FALSE;
  _ID = ID;
  _nRead = _nWrite = 0;
  _self = this;
}

template <class Type>
CacheBlock<Type>::~CacheBlock()
{
  if (_next)
    delete _next;

  _changed = FALSE;
  _ID = 0;
  _nRead = _nWrite = 0;
}

template <class Type>
CacheBlock<Type> *
CacheBlock<Type>::addBlock(unsigned ID, unsigned size)
{
  CacheBlock *previous = this;
  CacheBlock *block    = _next;

  while (block) {
    previous = block;
    block = block->_next;
  }

  block = new CacheBlock(ID, size);

  if (block)
    previous->_next = block;
    
  return block;
}

template <class Type>
Boolean
CacheBlock<Type>::read(fstream& s, unsigned ID)
{
  if (_changed) {
    if (this->_debug)
      cout << "<w" << _ID << ">" << flush;
    write(s);
  }
  if (this->_debug)
    cout << "<x" << _ID << "><r" << ID << ">" << flush;

  _ID = ID;
  _nRead = _nWrite = 0;
  _changed = FALSE;

  s.seekg(_ID*_nBytes);
  s.read((char *) this->_contents, _nBytes);

  return s ? TRUE : FALSE;
}

template <class Type>
Boolean
CacheBlock<Type>::write(fstream& s) const
{
  s.seekg(_ID*_nBytes);
  s.write((char *) this->_contents, _nBytes);

  return s ? TRUE : FALSE;
}

#ifdef __GNUC__
#define _INSTANTIATE_CACHEDARRAY(Type)                           \
  template class CachedArray<Type>;                              \
  template<> const unsigned CachedArray<Type>::_DEFAULT_BLOCK_SIZE = 32768; \
  template<> const unsigned CachedArray<Type>::_DEFAULT_N_BLOCKS   = 2;     \
  template CachedArray<Type> operator ^ (double, CachedArray<Type> const &); \
  template class CacheBlock<Type>;                               \
  template<> unsigned CachedArray<Type>::_rangeErrorCount = 25;


_INSTANTIATE_CACHEDARRAY(char);
_INSTANTIATE_CACHEDARRAY(float);
#endif // __GNUC__ defined
