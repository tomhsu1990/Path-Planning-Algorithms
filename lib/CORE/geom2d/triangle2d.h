/*****************************************************************
 * File: triangle2d.h
 * Synopsis:
 *      Basic 2-dimensional triangle geometry
 * Author: Chee and Tom Hsu (2016)
 *
 *****************************************************************
 * CORE Library Version 2.1 (Feb 2016)
 *       Chee Yap <yap@cs.nyu.edu>
 *       Tom Hsu <chhsu@cs.nyu.edu>
 *
 * Copyright (c) 1995, 1996, 1998, 1999, 2000, 2001, 2016 Exact Computation Project
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $Id: triangle2d.h,v 1.1 2006/04/03 18:55:31 exact Exp $
 *****************************************************************/

#ifndef _TRIANGLE2D_H
#define _TRIANGLE2D_H

#include <CORE/geom2d/point2d.h>
#include <CORE/geom2d/segment2d.h>
#include <CORE/geom2d/circle2d.h>

class Triangle2d : public GeomObj {

private:
 
  Point2d p1;  // the 3 points defining the triangle
  Point2d p2;
  Point2d p3;

  int orient;

  Point2d cp;

public:
  Triangle2d(const Point2d& p1, const Point2d& p2, const Point2d& p3);
   //initialized to the oriented triangle through points p1, p2, p3

  virtual ~Triangle2d();

  Triangle2d& operator=(const Triangle2d& T);

  //operations

  Point2d center();
  //return the center of the circle

  double radius();
  //returns the radius.
  //precond: the orientation of the circle is not 0

  Point2d point1() const { return p1; }
  Point2d point2() const { return p2; }
  Point2d point3() const { return p3; }

  bool is_degerate() const { return orient == 0; }
  //returns true if the defining points are collinear

  bool is_trivial() const {return p1 == p2; }
  //returns true if radius is zero

  int orientation() const { return orient; }

  int side_of(const Point2d& p) const;
  // returns -1, +1 or 0 if p lies right of, left of or on the triangle
  // respectively

  bool inside(const Point2d& p);
  //returns true if p lies inside of the triangle

  bool outside(const Point2d& p);

  bool contains(const Point2d& p) const ;
  //returns true if p lies on the triangle, false otherwise

  double distance(const Point2d& p);
  //returns the distance between p and the triangle

  double distance(const Line2d& l);
  //returns the distance between l and the triangle
  //distance from center to l minus radius

  double distance(Triangle2d& D);
  //returns the distance between this triangle and triangle D

  bool operator==(const Triangle2d& D) const ;
  bool operator!=(const Triangle2d& D) const
  { return !operator==(D); }

  friend std::ostream& operator<<(std::ostream& out, Triangle2d& t);
  friend std::istream& operator>>(std::istream& in, Triangle2d t); //?? Triangle2d &

  // decide if point p is inside triangle defined by
  // a, b, c
  static bool inside(Point2d a, Point2d b, Point2d c, Point2d p);

}; // class Triangle2d


#endif
