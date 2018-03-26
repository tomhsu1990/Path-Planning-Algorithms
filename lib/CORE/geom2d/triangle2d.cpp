/*****************************************************************
 * File: triangle2d.cpp
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
 * $Id: triangle2d.cpp,v 1.1 2006/04/03 18:55:31 exact Exp $
 *****************************************************************/

#include <CORE/geom2d/triangle2d.h>

Triangle2d::Triangle2d(const Point2d& P1, const Point2d& P2, const Point2d& P3)
    : p1(P1), p2(P2), p3(P3), cp(0, 0){
  orient = orientation2d(p1, p2, p3);
  if(orient == 0){
    //exit(1);
  }
}

Triangle2d::~Triangle2d(){
}

Triangle2d& Triangle2d::operator=(const Triangle2d& T){
  p1 = T.p1;
  p2 = T.p2;
  p3 = T.p3;
  orient = T.orient;
  cp = T.cp;
  return *this;
}

Point2d Triangle2d::center(){
  cp.setX((p1.X()+p2.X()+p3.X())/3.0f);
  cp.setY((p1.Y()+p2.Y()+p3.Y())/3.0f);
  return cp;
}

bool Triangle2d::operator==(const Triangle2d& t) const {
   if (!t.contains(p1)) return false;
   if (!t.contains(p2)) return false;
   if (!t.contains(p3)) return false;
   return true;
}

double Triangle2d::distance(const Point2d& p) {
   return 0;
}

double Triangle2d::distance(const Line2d& l) {
   return 0;
}

double Triangle2d::distance(Triangle2d& D) {
   return 0;
}

int Triangle2d::side_of(const Point2d& p) const {
   return 0;
}

bool Triangle2d::outside(const Point2d& p) {
   return (orient * side_of(p)) < 0;
}

bool Triangle2d::inside(const Point2d& p) {
   return (orient * side_of(p)) > 0;
}

bool Triangle2d::contains(const Point2d& p) const {
   return side_of(p) == 0;
}

std::ostream& operator<<(std::ostream& out, Triangle2d& t){
   out << "Triangle2d[ center=" << t.center() << "]";
   return out;
}

std::istream& operator>>(std::istream& in, Triangle2d t){
  Point2d p1, p2, p3;
  std::cout << "\nInput first point(-,-):";
  in >> p1;
  std::cout << "\nInput second point(-,-):";
  in >> p2;
  std::cout << "\nInput third point(-,-):";
  in >> p3;
  t = Triangle2d(p1, p2, p3);
  return in;
}

/*
 inside decides if a point P is Inside of the triangle
 defined by A, B, C.
*/
bool Triangle2d::inside(Point2d a, Point2d b, Point2d c, Point2d p){
  float ax, ay, bx, by, cx, cy, apx, apy, bpx, bpy, cpx, cpy;
  float cCROSSap, bCROSScp, aCROSSbp;

  ax = c.X() - b.X();  ay = c.Y() - b.Y();
  bx = a.X() - c.X();  by = a.Y() - c.Y();
  cx = b.X() - a.X();  cy = b.Y() - a.Y();
  apx= p.X() - a.X();  apy= p.Y() - a.Y();
  bpx= p.X() - b.X();  bpy= p.Y() - b.Y();
  cpx= p.X() - c.X();  cpy= p.Y() - c.Y();

  aCROSSbp = ax*bpy - ay*bpx;
  cCROSSap = cx*apy - cy*apx;
  bCROSScp = bx*cpy - by*cpx;

  return ((aCROSSbp >= 0.0f) && (bCROSScp >= 0.0f) && (cCROSSap >= 0.0f));
}
