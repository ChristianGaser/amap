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
$RCSfile: OrderedCltn.cc,v $
$Revision: 1.3 $
$Author: jason $
$Date: 2004-01-19 15:38:15 $
$State: Exp $
--------------------------------------------------------------------------*/
#include <iostream>		// (bert) changed from iostream.h
using namespace std;		// (bert) added

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#include "OrderedCltn.h"

/*
 * ordered collection iterator class functions
 */

void * 
ocIterator::operator ++() { // Prefix
  if (!pOC || (++uCurIndex >= pOC->uEndIndex))
    return 0;

  return pOC->Contents[uCurIndex];
}

void * 
ocIterator::operator ++(int) { // Postfix
  if (!pOC || (uCurIndex >= pOC->uEndIndex))
    return 0;

  return pOC->Contents[uCurIndex++];
}

void * 
ocIterator::operator --() { // Prefix
  if (!pOC || ((int) --uCurIndex < 0))		// typecasting because of unsigned
    return 0;

  return pOC->Contents[uCurIndex];
}
   
void * 
ocIterator::operator --(int) { // Postfix
  if (!pOC || ((int) uCurIndex < 0))		// typecasting because of unsigned
    return 0;

  return pOC->Contents[uCurIndex--];
}
   
/*
 * ordered collection class functions
 */

// display error and exit routines

void OrderedCltn::no_mem_err () const {
      cerr << "insufficient memory for collection of size " << uCapacity << "\n";
      exit (-1);
}   


void OrderedCltn::range_err (unsigned u) const {
      cerr << "Collection range error (index " << u
              << ", range [0 - " << uEndIndex-1 << "])\n";
      exit (-1);
}


void OrderedCltn::not_impl_err () const {
      cerr << "Attempted bogus operation on OrderedCltn" << endl;
}

// default constructor (NOTE: malloc is called instead of the usual new because
// realloc will be called when the array has to resized)

OrderedCltn::OrderedCltn ( unsigned uCap ) {
      if ( uCap == 0 ) 
            uCap = 1;
      uCapacity = uCap;
      uEndIndex = 0;
      Contents = (void **)malloc ( sizeof (void *) * uCap );
      if ( Contents == NULL )
           no_mem_err ();
}


OrderedCltn::OrderedCltn ( const OrderedCltn& OC ) {
      uCapacity = OC.uCapacity;
      uEndIndex = 0;
      Contents = (void **)malloc ( sizeof (void *) * uCapacity  );
      if ( Contents == NULL )
             no_mem_err ();
      ocIterator itOC (OC);
      void *pv;
      while ( pv = itOC++ )
            add (pv);
}


OrderedCltn::~OrderedCltn () {
      free ((char *)Contents);
}


// increase the capacity of the collection

void OrderedCltn::reSize (unsigned uNewCap) {
     if ( uNewCap <= uCapacity ) 
           return;
     uCapacity = uNewCap;
     Contents = (void **)realloc ( (char *)Contents, sizeof (void *) * uNewCap );
     if ( Contents == NULL )
           no_mem_err ();
}


/*
 *   find index of a specified pointer
 */

int
OrderedCltn::indexOf (const void* Ptr) const 
{
  for (register unsigned i = 0; i < uEndIndex; i++)
    if (Contents[i] == Ptr)
      return (int) i;

  return -1;
}


/*
 *   return the total number of occurrences of an object
 */

unsigned OrderedCltn::occurrencesOf (const void* Ptr) const  {
     ocIterator i_OC ( *this );
     register unsigned u = 0;
     register void* p;
     while ( p = i_OC++ )
     if ( Ptr == p ) 
             u++;
     return u;
}   
  
/*
 *   add a pointer at the specified index
 */

unsigned OrderedCltn::add ( const void* Newptr, unsigned uIndex )  {
     if ( uIndex > uEndIndex )
            range_err (uIndex);
     if ( uEndIndex == uCapacity ) 
            reSize (uCapacity + EXPANSION_INCREMENT);
     for ( register unsigned u = uEndIndex; u > uIndex; u-- )
            Contents[u] = Contents[u-1];

     Contents[uIndex] = (void*)Newptr;
     uEndIndex++;
     return (uIndex);
}

/* 
 *  add a new pointer after the specified pointer
 */

unsigned OrderedCltn::addAfter (const void* Ptr, const void *Newptr)  {
   register unsigned uIndex;
   if ((uIndex = indexOf (Ptr)) >= uEndIndex)
      return ~0;
   return add (Newptr, uIndex+1);
}

/*
 *  add an ordered collection of pointers at the beginning
 */

const OrderedCltn& 
OrderedCltn::addAllFirst (const OrderedCltn& cltn)  {
  unsigned n2add = cltn.size();
  if (n2add > 0) {
    if (uEndIndex + n2add >= uCapacity )
      reSize(uEndIndex + n2add + EXPANSION_INCREMENT);

    unsigned i;
    void **sourcePtr = Contents + uEndIndex - 1;
    void **destPtr   = Contents + uEndIndex + n2add - 1;
    for (i = uEndIndex; i; i--)
      *destPtr-- = *sourcePtr--;
    
    sourcePtr = cltn.Contents;
    destPtr   = Contents;
    for (i = n2add; i; i--)
      *destPtr++ = *sourcePtr++;

    uEndIndex += n2add;
  }

  return cltn;
}

/*
 *  add an ordered collection of pointers at the end
 */

const OrderedCltn& OrderedCltn::addAllLast (const OrderedCltn& cltn)  {
        if ( uEndIndex+cltn.size() >= uCapacity )
                reSize(uEndIndex+cltn.size()+EXPANSION_INCREMENT);
        for (register int i=0; i<cltn.size(); i++) 
	  Contents[uEndIndex++] = cltn.Contents[i];
        return cltn;
}

// add a new pointer before the specified pointer
unsigned OrderedCltn::addBefore (const void* Ptr, const void *Newptr)
{
   register unsigned uIndex;
   if ((uIndex = indexOf (Ptr)) >= uEndIndex)
      return ~0;
   return add (Newptr, uIndex);
}


// remove the contents at the specified index from the collection

void* OrderedCltn::remove (unsigned uIndex) {
   if (uIndex >= uEndIndex)
      range_err (uIndex);
   void* Delptr = Contents[uIndex];
   uEndIndex--;
   for (register unsigned u = uIndex; u < uEndIndex; u++)
      Contents[u] = Contents[u+1];
   return (Delptr);
}

void *
OrderedCltn::remove()
{
  if (uEndIndex != 0) 
    return remove(uEndIndex - 1); 
  return 0;
}

// remove the specified pointer if it exists

unsigned OrderedCltn::remove (const void* Ptr) {
   register unsigned uIndex;
   if ((uIndex = indexOf (Ptr)) >= uEndIndex)
     return 0;
   remove (uIndex);
   return 1;    
}

OrderedCltn& OrderedCltn::operator = (const OrderedCltn& OC) {   
  free ((char *)Contents);
  uCapacity = OC.uCapacity;
  uEndIndex = 0;
  Contents = (void **)malloc (sizeof (void *) * uCapacity);
  if (Contents == 0)
    no_mem_err();
  ocIterator itOC (OC);
  void *pv;
  while (pv = itOC++)
    add (pv);
  return (*this);
}


// equivalence (NOTE: all entries must be in the same sequence)

int OrderedCltn::operator ==(const OrderedCltn& OC) const {
   if (uEndIndex != OC.uEndIndex)
      return 0;
   ocIterator itThis (*this), 
                     itOC (OC);
   register void* pv;
   while (pv = itThis++)
      if (pv != itOC++) 
            return 0;
   return 1;
}


OrderedCltn& OrderedCltn::operator &=(const OrderedCltn& OC) {
   register unsigned uIndex = uEndIndex;
   while (uIndex--)
      if (!OC.indexOf ((const void*)Contents[uIndex])) 
	 remove (uIndex);  // ... remove it from this collection
   return *this;
}

