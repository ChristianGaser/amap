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
$RCSfile: MPoint.cc,v $
$Revision: 1.2 $
$Author: stever $
$Date: 2003-11-17 04:07:52 $
$State: Exp $
--------------------------------------------------------------------------*/
#include <config.h>
#include "MPoint.h"
#include <assert.h>
//#include "MRegion.h"

using namespace std;


Pool<MPoint>      MPoint::_pool; //      = Pool<MPoint>();
//Pool<MPoint3D>    MPoint3D::_pool    = Pool<MPoint3D>();
Pool<MWorldPoint> MWorldPoint::_pool; // = Pool<MWorldPoint>();

MPoint&
MPoint::magnify(int mag)
{
  if (mag)
    if (mag > 0) {
      x *= mag;
      y *= mag;
    }
    else {
      x /= -mag;
      y /= -mag;
    }

  return *this;
}

istream&
operator >> (istream& is, MPoint *&point)
{
  if (!point) {
    int x, y;
    is >> x >> y;
    point = new MPoint(x, y);
    assert(point);
  }
  else
    is >> point->x >> point->y; 

  return is;
}

istream&
operator >> (istream& is, MPoint& point)
{
  is >> point.x >> point.y; 
  return is;
}

ostream&
operator << (ostream& os, const MPoint& point)
{
  os << point.x << " " << point.y; 
  return os;
}







