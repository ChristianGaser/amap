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
$RCSfile: dcomplex.cc,v $
$Revision: 1.3 $
$Author: bert $
$Date: 2003-04-16 18:44:21 $
$State: Exp $
--------------------------------------------------------------------------*/

#include "dcomplex.h"
#include <iostream>		// (bert) changed from iostream.h
using namespace std;		// (bert) added

// A few functions that define (bogus) math ops for complex

static int _errorCount_dcomplex = 100;

int operator < (const dcomplex&, const dcomplex&) {
  if (_errorCount_dcomplex) {
    cerr << "Comparison of dcomplex numbers undefined" << endl;
    _errorCount_dcomplex--;
  }
  return 0;
}

int operator <= (const dcomplex&, const dcomplex&) {
  if (_errorCount_dcomplex) {
    cerr << "Comparison of dcomplex numbers undefined" << endl;
    _errorCount_dcomplex--;
  }
  return 0;
}

int operator > (const dcomplex&, const dcomplex&) {
  if (_errorCount_dcomplex) {
    cerr << "Comparison of dcomplex numbers undefined" << endl;
    _errorCount_dcomplex--;
  }
  return 0;
}

int operator >= (const dcomplex&, const dcomplex&) {
  if (_errorCount_dcomplex) {
    cerr << "Comparison of dcomplex numbers undefined" << endl;
    _errorCount_dcomplex--;
  }
  return 0;
}

int operator && (const dcomplex&, const dcomplex&) {
  if (_errorCount_dcomplex) {
    cerr << "Comparison of dcomplex numbers undefined" << endl;
    _errorCount_dcomplex--;
  }
  return 0;
}

int operator || (const dcomplex&, const dcomplex&) {
  if (_errorCount_dcomplex) {
    cerr << "Comparison of dcomplex numbers undefined" << endl;
    _errorCount_dcomplex--;
  }
  return 0;
}

