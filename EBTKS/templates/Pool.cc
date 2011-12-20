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
$RCSfile: Pool.cc,v $
$Revision: 1.2 $
$Author: bert $
$Date: 2003-04-16 16:57:23 $
$State: Exp $
--------------------------------------------------------------------------*/
#include <config.h>
#include "Pool.h"
#include <assert.h>
#include <iostream>		// (bert) changed from iostream.h
using namespace std;		// (bert) added

template <class Type>
Pool<Type>::Pool(unsigned nElements)
: _elementSize(sizeof(Type)),
  _expansionSize(nElements),
  _blocks(_POOL_N_BLOCKS)
{
  assert(_elementSize >= sizeof(_Link *));
  assert(nElements);
  _head = 0;
}

template <class Type>
Pool<Type>::~Pool()
{
  ocIterator blockIt(_blocks);
  Type      *block;
  while (block = (Type *) blockIt++)
    delete [] block;
}

template <class Type>
void
Pool<Type>::_grow()
{
  Type *newBlock = new Type[_expansionSize];
  assert(newBlock);

  _blocks.add(newBlock);
  
  Type *elementPtr = newBlock;
  for (unsigned i = _expansionSize - 1; i; i--, elementPtr++)
    ((_Link *) elementPtr)->next = (_Link *) (elementPtr + 1);

// This will also do the trick; the _Link struct is really not necessary
//    *(Type **) elementPtr = elementPtr + 1;

  ((_Link *) elementPtr)->next = 0;

  _head = (_Link *) newBlock;
}

#ifdef __GNUC__
#define _INSTANTIATE_POOL(Type)                       \
         template class Pool<Type>;

#include "MPoint.h"
_INSTANTIATE_POOL(MPoint);
_INSTANTIATE_POOL(MPoint3D);
_INSTANTIATE_POOL(MWorldPoint);
#endif
