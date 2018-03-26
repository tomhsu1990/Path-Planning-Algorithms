/*****************************************************************
 * File: line3d.cc
 * Synopsis:
 *      Basic 3-dimensional geometry
 * Author: Shubin Zhao (shubinz@cs.nyu.edu), 2001.
 *
 *****************************************************************
 * CORE Library Version 1.4 (July 2001)
 *       Chee Yap <yap@cs.nyu.edu>
 *       Chen Li <chenli@cs.nyu.edu>
 *       Zilin Du <zilin@cs.nyu.edu>
 *
 * Copyright (c) 1995, 1996, 1998, 1999, 2000, 2001 Exact Computation Project
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $Id: line3d.cpp,v 1.1 2006/04/03 18:55:31 exact Exp $
 *****************************************************************/

#include <CORE/geom3d/line3d.h>

/************************************************************
 * *   constructors
************************************************************/

Line3d::Line3d(const Point3d &p, const Vector &v) : p0(p), p1(p + v), V(v) { }

Line3d::Line3d(const Point3d &_p0, const Point3d &_p1) : p0(_p0), p1(_p1), V(_p1 - _p0) { }

Line3d::Line3d(const Line3d &l) : p0(l.p0), p1(l.p1), V(l.V) { }

Line3d::Line3d() : p0(0.0, 0.0, 0.0), p1(0.0, 0.0, 0.0), V(0) {}
 // line passes through the origin with direction 0.

/************************************************************
 * *   Member functions 
************************************************************/

double Line3d::distance(const Point3d& q) const {
  Vector Vpq = q - p0;
  Vector u = V.cross( Vpq );
  return u.norm()/V.norm();
}

Point3d Line3d::projection(const Point3d& q) const {
  Vector Vpq = q - p0;
  double lambda = dotProduct(Vpq, V)/dotProduct(V, V);
  return p0 + lambda * V;
}

bool Line3d::contains(const Point3d& p) const {
  return V.cross( p - p0 ).isZero();
}

bool Line3d::isCoincident(const Line3d& l) const {
  return contains( l.startPt() ) && contains( l.stopPt() );
}

bool Line3d::isSkew(const Line3d& l2) const {
  //assert(intersects(l2) >= 0);
  double d = dotProduct(p0 - l2.startPt(), V.cross(l2.direction()));
  return (d != 0);
}

bool Line3d::isParallel(const Line3d& l2) const {
  return V.cross(l2.direction()).isZero();
}

int Line3d::intersects(const Line3d &l) const {
  // return the dimension of intersection
  // return -1 if not intersect
  if( !coplanar(p0, p1, l.startPt(), l.stopPt()) )
    return -1;
  else if( isCoincident( l ) ) 
    return 1;
  else if( isParallel( l ) ) 
    return -1;
  else 	// intersection is point
    return 0;    
}

GeomObj* Line3d::intersection(const Line3d &l) const {
  if( intersects( l ) == -1 ) 
    return NULL;

  if( isCoincident( l ) )
    return new Line3d(*this);
    
  Vector u = l.direction().cross(V);
  Vector w = l.direction().cross( l.startPt() - p0 );
  double lambda = dotProduct(w,u) / dotProduct(u,u);
  Vector t = lambda * V;
  return new Point3d(p0 + t);
}

std::ostream &operator<<(std::ostream &o, const Line3d &l) {
  return o << "Line3d[" << l.p0 << "," << l.p1 << ";" << l.V << "]";
}

std::istream& operator>>(std::istream& in, Line3d& l) {
  Point3d pStart, pEnd;
  std::cout << "\nInput start point(-,-,-):";
  in >> pStart;
  std::cout << "Input end point(-,-,-):";
  in >> pEnd;
  l = Line3d(pStart, pEnd); 
  return in;
}
