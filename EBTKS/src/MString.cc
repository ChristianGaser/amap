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
$RCSfile: MString.cc,v $
$Revision: 1.2 $
$Author: stever $
$Date: 2003-11-17 04:07:52 $
$State: Exp $
--------------------------------------------------------------------------*/
#include <config.h>
#include "MString.h"
#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include "OrderedCltn.h"

using namespace std;


unsigned MString::_MAX_LOAD_LENGTH = 512;

//
// Constructors/destructor
//

MString::MString(int length)
: SimpleArray<char>((unsigned) length + 1)
{
  _contents[0] = '\0';
}

MString::MString(const char *string)
: SimpleArray<char>(strlen(string) + 1)
{
  strcpy(_contents, string);
}

MString::MString(const MString& mString)
: SimpleArray<char>(mString._size)
{
  strcpy(_contents, mString._contents);
}

MString::MString(char c, unsigned n)
: SimpleArray<char>(n + 1)
{
  for (unsigned i = 0; i < n; i++)
    _contents[i] = c;

  _contents[n] = '\0';
}

MString::~MString()
{}

//
// Get functions
//

Boolean
MString::isInteger(int *value) const
{
  if (_size <= 1)
    return FALSE;

  unsigned nChars = _size;
  char *charPtr = _contents;
  if (*charPtr == '-') {
    if (_size == 2) 
      return FALSE; // A dash by itself is no integer
    else  {
      charPtr++;
      nChars--;
    }
  }

  for (unsigned i = nChars - 1; i != 0; i--, charPtr++)
    if (!isdigit(*charPtr))
      return FALSE;

  if (value)
    return sscanf(_contents, "%d", value);

  return TRUE;
}

Boolean
MString::isDouble(double *value) const 
{
  if (_size <= 1)
    return FALSE;

  char *charPtr = _contents;
  for (unsigned i = _size - 1; i != 0; i--, charPtr++)
    if (!isdigit(*charPtr) && (*charPtr != '-') && (*charPtr != '.'))
      return FALSE;

  if (value)
    return sscanf(_contents, "%lf", value);

  return TRUE;
}

Boolean
MString::contains(const MString& subString) const 
{
  return (strstr(_contents, subString._contents) != 0); 
}
 
 Boolean  
MString::contains(const char *charPtr) const 
{
  return (strstr(_contents, charPtr) != 0); 
}

//
// Conversion operators
//

MString::operator unsigned() const
{
  unsigned value = 0;
  if (!sscanf(_contents, "%u", &value))
    cerr << "Warning! Couldn't convert " << *this << " to unsigned" << endl;
  return value;
}

MString::operator int() const
{
  int value = 0;
  if (!sscanf(_contents, "%d", &value))
    cerr << "Warning! Couldn't convert " << *this << " to int" << endl;
  return value;
}

MString::operator float() const
{
  float value = 0;
  if (!sscanf(_contents, "%f", &value))
    cerr << "Warning! Couldn't convert " << *this << " to float" << endl;
  return value;
}

MString::operator double() const
{
  double value = 0;
  if (!sscanf(_contents, "%lf", &value))
    cerr << "Warning! Couldn't convert " << *this << " to double" << endl;
  return value;
}

MString&
MString::pad(char c, int n)
{
  if (n > 0)
    *this += MString(c, n);
  else if (n < 0) {
    n = ::abs(n);
    *this = MString(c, n) + *this;
  }
  
  return *this;
}

MString&
MString::pad(const MString& str, int n)
{
  if (!str) {
    cerr << "MString::pad: attempt to pad with empty string" << endl;
    return *this;
  }

  if (n > 0)
    for (unsigned i = 0; i < n; i++)
      *this += str;
  else if (n < 0) {
    MString temp(str);
    n = ::abs(n);
    for (unsigned i = 1; i < n; i++)
      temp += str;

    temp += *this;

    *this = temp;
  }

  return *this;
}

Boolean 
MString::operator < (const MString& string) const
{
//  return length() < string.length(); 
  return (strcmp(_contents, string._contents) < 0);
}

Boolean 
MString::operator > (const MString& string) const
{
//  return length() > string.length(); 
  return (strcmp(_contents, string._contents) > 0);
}

//
// Operators
//

MString&
MString::operator = (const char *string)
{
  if (!string)
    return *this;

  newSize(strlen(string) + 1);
  
  strcpy(_contents, string);

  return *this;
}

MString&
MString::operator = (const MString& mString)
{
  if (this == &mString) 
    return *this;

  newSize(mString.length() + 1);

  strcpy(_contents, mString._contents);

  return *this;
}

MString&
MString::operator = (int value)
{
  char temp[200];
  sprintf(temp, "%d", value);

  return operator = (temp);
}

MString&
MString::operator = (char value)
{
  newSize(2);
  _contents[0] = value;
  _contents[1] = '\0';

  return *this;
}

MString&
MString::operator = (double value)
{
  char temp[200];
  sprintf(temp, "%.8g", value);

  return operator = (temp);
}

MString&
MString::operator += (const char *string)
{
  if (!string)
    return *this;

  newSize(length() + strlen(string) + 1);

  strcat(_contents, string);

  return *this;
}

MString&
MString::operator += (const MString& mString)
{
  if (mString.isEmpty())
    return *this;

  newSize(length() + mString.length() + 1);

  strcat(_contents, mString._contents);

  return *this;
}

MString&
MString::operator += (char extraChar)
{
  newSize(length() + 2);

  _contents[_size - 2] = extraChar;
  _contents[_size - 1] = '\0';

  return *this;
}

MString&
MString::operator += (int value)
{
  char temp[200];
  sprintf(temp, "%d", value);

  return operator += (temp);
}

MString&
MString::operator += (double value)
{
  char temp[200];
  sprintf(temp, "%.8g", value);

  return operator += (temp);
}

MString
MString::operator + (const char *string) const
{
  MString resultString(*this);
  resultString += string;
  return resultString;
}

MString
MString::operator + (const MString& mString) const
{
  MString resultString(*this);
  resultString += mString;
  return resultString;
}

MString
MString::operator () (unsigned index) const
{
  if (index >= length())
    return MString(0);

  MString newString(length() - index);
  for (register int i = index, j = 0; i < _size; i++, j++)
    newString._contents[j] = _contents[i];

  return newString;
}

MString
MString::operator () (unsigned index, unsigned nChar) const
{
  if (!nChar)
    return MString(0);

  if (index + nChar + 1 > _size)
    nChar = _size - index + 1;

  MString newString(nChar);
  for (register int i = index, j = 0; j < nChar; i++, j++)
    newString._contents[j] = _contents[i];
  newString._contents[nChar] = '\0';

  return newString;
}

//
// Special functions
//

//
// Force a string template onto the string. Both strings will be broken
// up into tokens seperated by <separator>. The template tokens will replace
// the corresponding tokens in the original string, except when the template
// token is a '*'.
//

MString&
MString::applyTemplate(const MString& templateString, const char *separator)
{
  if (this->isEmpty() || templateString.isEmpty()) {
    delete [] _contents;
    _maxSize = _size = templateString._size;
    _contents = new char[_size];
    assert(_contents);
    strcpy(_contents, templateString._contents);
    return *this;
  }

  MString      templateCopy(templateString);
  char        *tokenPtr;
  OrderedCltn  sourceTokens, templateTokens;

// Store pointers to template tokens in <templateTokens>
// (strtok replaces the separators with NULL characters in <templateCopy>)

  tokenPtr = strtok(templateCopy._contents, separator);
  while (tokenPtr) {
    templateTokens.add((void *) tokenPtr);
    tokenPtr = strtok(NULL, separator);
  }

// Same for <sourceTokens>. Don't need to copy the source string, as it will
// be destroyed anyway.

  tokenPtr = strtok(_contents, separator);
  while (tokenPtr) {
    sourceTokens.add((void *) tokenPtr);
    tokenPtr = strtok(NULL, separator);
  }
  
// Store template tokens in the target string; replace * by the coresponding
// tokens of the source string

  int nExtraSourceTokens = sourceTokens.size() - templateTokens.size();

  char      *sourceToken, *templateToken;
  MString    targetString;
  ocIterator templateTokenIt(templateTokens);
  ocIterator sourceTokenIt(sourceTokens);

  while (templateToken = (char *) templateTokenIt++) {
    if (sourceToken = (char *) sourceTokenIt++) {
      int templateTokenLength = strlen(templateToken);
      for (int j = 0; j < templateTokenLength; j++) {
	char character = templateToken[j];
	if (character == '*') {
	  targetString += sourceToken;
	  for (int k = 1; k <= nExtraSourceTokens; k++) {
	    targetString += separator;
	    targetString += (char *) sourceTokenIt++;
	  }
	}
	else
	  targetString += character;
      }
    }
    else
      targetString += templateToken;

    targetString += separator;
  }
    
    
// Store targetString as final result; cut off the last added separator.

  delete [] _contents;
  _maxSize = _size = targetString._size - 1;
  _contents = new char [_size];
  assert(_contents);
  strncpy(_contents, targetString._contents, _size - 1);
  _contents[_size - 1] = '\0';

  return *this;
}

MString
MString::chop(unsigned n)
{
  unsigned strLength = length();
  strLength -= MIN(n, strLength);
  MString tail((*this)(strLength));
  newSize(strLength + 1);
  _contents[strLength] = '\0';

  return tail;
}

ostream&
operator << (ostream& os, const MString& mString)
{
  return (os << mString._contents);
}

istream&
operator >> (istream& is, MString& mString)
{
  if (mString._maxSize < MString::_MAX_LOAD_LENGTH) {
    if (mString._maxSize)
      delete [] mString._contents;
    mString._size = mString._maxSize = MString::_MAX_LOAD_LENGTH;
    mString._contents = new char[mString._size];
    assert(mString._contents);
  }
  else
    mString._size = mString._maxSize;
  
  mString._contents[0] = '\0';

  is.width(mString._size);

  is >> mString._contents;

  mString._size = strlen(mString._contents) + 1;

  return is;
}

/***********************
 * MStringIterator class
 ***********************/

MStringIterator::MStringIterator(const MString& mString, unsigned start)
: _separator(0)
{
  _mString = &mString;
  _index   = start;
}

MStringIterator::MStringIterator(const MString& mString, const MString& separator,
				 unsigned start)
: _separator(separator)
{
  _mString = &mString;
  _index = start;
}

MStringIterator::~MStringIterator()
{
}

MString
MStringIterator::first()
{
  _index = 0;

  if (!_separator) {
    MString firstChar(1);
    firstChar += _mString->_contents[_index];
    return firstChar;
  }

  MString token(_nextToken());

  _index = 0;

  return token;
}

MString 
MStringIterator::start(unsigned index)
{
  if (_index >= _mString->_size)
    return MString(0);

  _index = index;

  if (!_separator) {
    MString firstChar(1);
    firstChar += _mString->_contents[_index];
    return firstChar;
  }

  MString token(_nextToken());

  _index = index;

  return token;
}

char
MStringIterator::character()
{
  if (_index >= _mString->_size)
    return '\0';
  
  return _mString->_contents[_index++];
}

MString
MStringIterator::operator ++ ()
{
  if (!_separator) {
    if (_index >= _mString->_size)
      return MString(0);
    else {
      MString character(1);
      character += _mString->_contents[_index++];
      return character;
    }
  }

  return _nextToken();
}

MString
MStringIterator::_nextToken()
{
  if (_index >= _mString->_size)
    return MString(0);

  char character = _mString->_contents[_index++];

  // Skip leading separators

  while ((character != '\0') && _separator.contains(character))
    character = _mString->_contents[_index++];

  if (character == '\0')
    return MString(0);

  _index--;

  MString token;
  Boolean gotToken = FALSE;

  // Extract the token

  while (!gotToken) {
    char character = _mString->_contents[_index++];
    gotToken = ((character == '\0') || _separator.contains(character));
    if (!gotToken)
      token += character;
  }

  return token;
}
