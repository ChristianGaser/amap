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
$RCSfile: Dictionary.cc,v $
$Revision: 1.2 $
$Author: bert $
$Date: 2003-04-16 15:03:14 $
$State: Exp $
--------------------------------------------------------------------------*/
#include "Dictionary.h"
#include <assert.h>
#include <iostream>		// (bert) changed from iostream.h
using namespace std;		// (bert) added

/******************
 * Dictionary class
 ******************/

//
// Constructors
//

template <class Key, class Value> 
Dictionary<Key, Value>::Dictionary(const Dictionary<Key, Value>& dict)
: _defaultKey(dict._defaultKey),
  _defaultValue(dict._defaultValue)
{
  _initialize();

  DictIterator<Key, Value> entryIt(dict);
  DictEntry<Key, Value>   *entry;
  while (entry = entryIt++)
    (*this)[entry->key] = entry->value;
}

//
// Copy operator
//

template <class Key, class Value> 
Dictionary<Key, Value>&
Dictionary<Key, Value>::operator = (const Dictionary<Key, Value>& dict)
{
  delete _head;
  _initialize();

  _defaultKey   = dict._defaultKey;
  _defaultValue = dict._defaultValue;

  DictIterator<Key, Value> entryIt(dict);
  DictEntry<Key, Value>   *entry;
  while (entry = entryIt++)
    (*this)[entry->key] = entry->value;

  return *this;
}

//
// Get functions
//

template <class Key, class Value> 
Value *
Dictionary<Key, Value>::at(const Key& requestedKey)
{
  if (!_head)
    return 0;

// Locate key  

  DictEntry<Key, Value> *dictEntry = _head;

  while(1) {
    if (requestedKey != dictEntry->key) {
      if (!(dictEntry = dictEntry->_next))
	return 0;
    }
    else
      return &dictEntry->value;
  }
}

template <class Key, class Value> 
Value&
Dictionary<Key, Value>::operator [] (const Key& requestedKey)
{
// Empty dictionary, create a new entry

  if (!_head) {
    _current = new DictEntry<Key, Value>(requestedKey, _defaultValue);
    assert(_current);
    _size++;
    _current->_next = _current->_previous = 0;
    _head = _tail = _current;
    return _current->value;
  }

// Locate key  

  DictEntry<Key, Value> *dictEntry = _current;

  if (requestedKey < dictEntry->key) // Search backward
    while(1) {
      if (requestedKey != dictEntry->key) {
	if (requestedKey > dictEntry->key) { // Insert after dictEntry
	  _current = new DictEntry<Key, Value>(requestedKey, _defaultValue);
	  assert(_current);
	  _size++;
	  _current->_next = dictEntry->_next;
	  _current->_previous = dictEntry;
	  if (dictEntry != _tail)
	    dictEntry->_next->_previous = _current;
	  else
	    _tail = _current;
	  dictEntry->_next = _current;
	  return _current->value;
	}
      
	DictEntry<Key, Value> *previous = dictEntry->_previous;
	if (!previous) {             // Insert at the beginning
  	  _current = new DictEntry<Key, Value>(requestedKey, _defaultValue);
	  assert(_current);
	  _size++;
	  _current->_next = dictEntry;
	  _current->_previous = 0;
	  dictEntry->_previous = _current;
	  _head = _current;
	  return _current->value;
	}
      
	dictEntry = previous;
      }
      else
	return dictEntry->value; // Key found; return value
    }
  else                                // Search forward
    while(1) {
      if (requestedKey != dictEntry->key) {
	if (requestedKey < dictEntry->key) { // Insert before dictEntry
	  _current = new DictEntry<Key, Value>(requestedKey, _defaultValue);
	  assert(_current);
	  _size++;
	  _current->_previous = dictEntry->_previous;
	  _current->_next  = dictEntry;
	  if (dictEntry != _head)
	    dictEntry->_previous->_next = _current;
	  else
	    _head = _current;
	  dictEntry->_previous = _current;
	  return _current->value;
	}
	
	DictEntry<Key, Value> *next = dictEntry->_next;
	if (!next) {             // Append at the end
	  _current = new DictEntry<Key, Value>(requestedKey, _defaultValue);
	  assert(_current);
	  _size++;
	  _current->_previous = dictEntry;
	  _current->_next = 0;
	  dictEntry->_next = _current;
	  _tail = _current;
	  return _current->value;
	}
	
	dictEntry = next;
      }
      else
	return dictEntry->value;  // Key found; return value
    }
}

template <class Key, class Value> 
Value
Dictionary<Key, Value>::operator [] (const Key& requestedKey) const
{
// Empty dictionary, create a new entry

  if (!_head) {
    cerr << "Warning! Requested key " << requestedKey 
	 << " not in dictionary; returning default value" << endl;
    return _defaultValue;
  }

// Locate key  

  DictEntry<Key, Value> *dictEntry = _current;

  if (requestedKey < dictEntry->key) // Search backward
    while(1) {
      if (requestedKey != dictEntry->key) {
	if (!(dictEntry = dictEntry->_previous)) {
	  cerr << "Warning! Requested key not in dictionary; returning default value"
	       << endl;
	  return _defaultValue;
	}
      }
      else
	return dictEntry->value;
    }
  else // Search forward
    while(1) {
      if (requestedKey == dictEntry->key)  // Key found; return value
	return dictEntry->value;
      
      dictEntry = dictEntry->_next;
      
      if (!dictEntry) {
	cerr << "Warning! Requested key not in dictionary; returning default value"
	     << endl;
	return _defaultValue;
      }
    }
}

template <class Key, class Value> 
const Key&
Dictionary<Key, Value>::keyOf(const Value& value)
{
  if (!_head) {
    cerr << "Warning! Requested value not in dictionary; returning default key" << endl;
    return _defaultKey;
  }

// Locate key  

  DictEntry<Key, Value> *dictEntry = _head;

  while(1) {
    if (dictEntry->value == value)
      return dictEntry->key;
    
    dictEntry = dictEntry->_next;

    if (!dictEntry) {
      cerr << "Warning! Requested value not in dictionary; returning default key"
	   << endl;
      return _defaultKey;
    }
  }
}

template <class Key, class Value> 
Boolean
Dictionary<Key, Value>::containsKey(const Key& key) const
{
  if (!_head)
    return FALSE;

  DictEntry<Key, Value> *dictEntry = _head;

  while(1) {
    if (dictEntry->key == key)
      return TRUE;
    
    dictEntry = dictEntry->_next;

    if (!dictEntry)
      return FALSE;
  }
}

template <class Key, class Value> 
Boolean
Dictionary<Key, Value>::containsValue(const Value& value) const
{
  if (!_head)
    return FALSE;

  DictEntry<Key, Value> *dictEntry = _head;

  while(1) {
    if (dictEntry->value == value)
      return TRUE;
    
    dictEntry = dictEntry->_next;

    if (!dictEntry)
      return FALSE;
  }
}

//
// Entry removing functions
//

template <class Key, class Value> 
void
Dictionary<Key, Value>::remove(const Key& requestedKey)
{
  if (!_head)
    return;

// Locate key  

  DictEntry<Key, Value> *dictEntry = _head;

  while(1) {
    if (requestedKey == dictEntry->key) { // Key found; delete it
      dictEntry->_previous->_next = dictEntry->_next;
      dictEntry->_next->_previous = dictEntry->_previous;
      delete dictEntry;
      _size--;
    }

    dictEntry = dictEntry->_next;

    if (!dictEntry)
      return;
  }
}

//
// Operators
//
/*
template <class Key, class Value>
Dictionary<Key, Value>&
Dictionary<Key, Value>::operator += (const Dictionary<Key, Value>& dict)
{
  DictIterator<Key, Value> dictIt(dict);
  DictEntry<Key, Value>   *dictEntry;
  while (dictEntry = dictIt++)
    (*this)[dictEntry->key] += dictEntry->value;

  return *this;
}

Dictionary<void *, void *>&
Dictionary<void *, void *>::operator += (const Dictionary<void *, void *>& dict)
{
  DictIterator<Key, Value> dictIt(dict);
  DictEntry<Key, Value>   *dictEntry;
  while (dictEntry = dictIt++)
    os << dictEntry->key << ": " << dictEntry->value << endl;

  return *this;
}
*/

template <class Key, class Value>
Boolean
Dictionary<Key, Value>::operator == (const Dictionary<Key, Value>& dict) const
{
  if (_size != dict._size)
    return FALSE;

  if (!_head)
    return TRUE;

  DictIterator<Key, Value> entry1it(*this);
  DictIterator<Key, Value> entry2it(dict);
  DictEntry<Key, Value>   *entry1, *entry2;
  while ((entry1 = entry1it++) && (entry2 = entry2it++) &&
	 (entry1->key == entry2->key) && (entry1->value == entry2->value));

  if (entry1 || entry2)
    return FALSE;

  return TRUE;
}

//
// Private functions
//

template <class Key, class Value> 
void
Dictionary<Key, Value>::_notImplementedError()
{
  cerr << "Not implemented Dictionary function called" << endl;
}

template <class Key, class Value>
ostream&
operator << (ostream& os, const Dictionary<Key, Value>&)
{
  cerr << "Dictionary::operator << not implemented" << endl;
  return os;
}

/***************************
 * Dictionary iterator class
 ***************************/

template <class Key, class Value>
DictEntry<Key, Value> *
DictIterator<Key, Value>::operator ++ () // Prefix increment operator
{
  if (!_dictEntry)
    return 0;

  return _dictEntry = _dictEntry->_next;
}

template <class Key, class Value>
DictEntry<Key, Value> *
DictIterator<Key, Value>::operator ++ (int)  // Postfix increment operator
{
  if (!_dictEntry)
    return 0;

  DictEntry<Key, Value> *previous = _dictEntry;
  _dictEntry = _dictEntry->_next;

  return previous;
}

template <class Key, class Value>
DictEntry<Key, Value> *
DictIterator<Key, Value>::operator -- () // Prefix decrement operator
{
  if (!_dictEntry)
    return 0;

  return _dictEntry = _dictEntry->_previous;
}

template <class Key, class Value>
DictEntry<Key, Value> *
DictIterator<Key, Value>::operator -- (int) // Postfix decrement operator
{
  if (!_dictEntry)
    return 0;

  DictEntry<Key, Value> *next = _dictEntry;
  _dictEntry = _dictEntry->_previous;

  return next;
}

// g++ wants these declarations.... aargh. Supposedly this will be fixed in v2.8
#ifdef __GNUC__
#define _INSTANTIATE_DICTIONARY(keytype, valuetype) \
         template class DictEntry<keytype, valuetype>;  \
         template class DictIterator<keytype, valuetype>;  \
         template class Dictionary<keytype, valuetype>;

#ifdef HAVE_MATLAB
#include "mat.h"
#include "MString.h"
_INSTANTIATE_DICTIONARY(matfile *, MString);
#endif
#endif

