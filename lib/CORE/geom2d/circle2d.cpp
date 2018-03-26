/*****************************************************************
 * File: circle2d.cc
 * Synopsis:
 *      Basic 2-dimensional geometry
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
 * $Id: circle2d.cpp,v 1.1 2006/04/03 18:55:31 exact Exp $
 *****************************************************************/


#include <CORE/geom2d/point2d.h>
#include <CORE/geom2d/circle2d.h>

Circle2d::Circle2d(const Point2d& P1, const Point2d& P2, const Point2d& P3)   
    : p1(P1), p2(P2), p3(P3), cp(0), rp(0)
{
   orient = orientation2d(p1, p2, p3);
   if (orient == 0 ) {
         std::cout << "Circle2d::Circle2d: non-admissable triple";
         exit(1);
      }
}

Circle2d::~Circle2d() {
   if (cp) delete cp;
   if (rp) delete rp;
}

 //center a, with b0 on circle
Circle2d::Circle2d(const Point2d& a, const Point2d& b0) : p1(b0)
{  
   p2 = b0;
   p2.rotate90(a);
   p3 = p2;
   p3.rotate90(a);
   cp = new Point2d(a);
   rp = new double(a.distance(b0));
}

Circle2d::Circle2d() : p1(0,0), p2(0,0), p3(0,0) { }

Circle2d::Circle2d(const Point2d& m, double r) 
{
   Point2d a(m.X(), m.Y()+r);
   Point2d b(m.X(), m.Y()-r);
   Point2d c(m.X()+r, m.Y());
   p1 = a;
   p2 = b;
   p3 = c;
   cp = new Point2d(m);
   rp = new double(r);
}

Circle2d& Circle2d::operator=(const Circle2d& C)
{
   p1 = C.p1;
   p2 = C.p2;
   p3 = C.p3;
   orient = C.orient;
   cp = C.cp;
   rp = C.rp;
   return *this;
} 

Point2d Circle2d::center() {
   if (cp == 0) {
      if (collinear(p1, p2, p3)) {
 	 std::cout << "Circle2d::center(): points are collinear." << std::endl;
      }
      Line2d l1 = p_bisector(p1, p2);
      Line2d l2 = p_bisector(p2, p3);
      Point2d m;
      cp = (Point2d *)l1.intersection(l2);
   }

   return *cp;
}

double Circle2d::radius() {
   if (rp == 0)
      rp = new double( center().distance(p1) );
   return *rp;
}

bool Circle2d::operator==(const Circle2d& c) const {
   if (!c.contains(p1)) return false;
   if (!c.contains(p2)) return false;
   if (!c.contains(p3)) return false;
   return true;
}

double Circle2d::distance(const Point2d& p) {
   double d = p.distance(center()) - radius();
   return d;
}

double Circle2d::distance(const Line2d& l) {
   double d = l.distance(center()) - radius();
   if (d < 0) d = 0;
   return d;
}

double Circle2d::distance(Circle2d& c) {
   double d = center().distance(c.center()) - radius() - c.radius();
   if (d < 0) d = 0;
   return d;
}

int Circle2d::side_of(const Point2d& p) const {
   double ax = p1.X();
   double ay = p1.Y();
 
   double bx = p2.X() - ax;
   double by = p2.Y() - ay;
   double bw = bx*bx + by*by;

   double cx = p3.X() - ax;
   double cy = p3.Y() - ay;
   double cw = cx*cx + cy*cy; 
     
   double d1 = cy*bw - by*cw;
   double d2 = bx*cw - cx*bw;
   double d3 = by*cx - bx*cy;

   double px = p.X() - ax;
   double py = p.Y() - ay;
   double pw = px*px + py*py;

   double D = px*d1 + py*d2 + pw*d3;

   if (D != 0)
      return (D > 0) ? 1 : -1;
   else
      return 0;
} 
 
bool Circle2d::outside(const Point2d& p) {
   return (orient * side_of(p)) < 0;
}

bool Circle2d::inside(const Point2d& p) {
   return (orient * side_of(p)) > 0;
}

bool Circle2d::contains(const Point2d& p) const {
   return side_of(p) == 0;
}

std::ostream& operator<<(std::ostream& out, Circle2d& c) {
   out << "Circle2d[ center=" << c.center() << ", radius=" << c.radius() << "]";
   return out;
}

std::istream& operator>>(std::istream& in, Circle2d c) {
  Point2d p1, p2, p3;
  std::cout << "\nInput first point(-,-):";
  in >> p1;
  std::cout << "\nInput second point(-,-):";
  in >> p2;
  std::cout << "\nInput third point(-,-):";
  in >> p3;
  c = Circle2d(p1, p2, p3);
  return in;
}

   
