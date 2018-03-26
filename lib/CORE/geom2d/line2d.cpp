/*****************************************************************
 * File: line2d.cc
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
 * $Id: line2d.cpp,v 1.2 2009/02/12 03:50:53 exact Exp $
 *****************************************************************/


#include <CORE/geom2d/line2d.h>

Line2d::Line2d(const Point2d &p, const Vector &v) 
               : p0(p), p1(p + v), V(v) { assert(p != p+v); }
Line2d::Line2d(const Point2d &p, const Point2d &q) 
                : p0(p), p1(q), V(p - q) { assert(p != q); }
Line2d::Line2d(const Line2d &l) : p0(l.p0), p1(l.p1), V(l.V) 
                { assert(l.p0 != l.p1); }

Line2d::Line2d(const double& a, const double& b, const double& c){
  	// line with equation aX+bY+c=0
	if (b != 0) {
		p0 = Point2d(0.0, -c/b);
		p1 = Point2d(1.0, (-c-a)/b);
	}
	else {
		p0 = Point2d(-c/a, 0.0);
		p1 = Point2d(-c/a, 1.0);
	}
 	V = p0 - p1;
}

Line2d::Line2d() : p0(0.0, 0.0), p1(0.0, 0.0), V(0) {}
 // line passes through the origin with direction 0.
 // improper line 

 /*************************************************************
  *   distance and others
 *************************************************************/

double Line2d::distance(Point2d p) const {
  double d = (p1.X()-p0.X())*(p0.Y()-p.Y()) -
             (p1.Y()-p0.Y())*(p0.X()-p.X());
  d /= p0.distance(p1);
  return (d < 0) ? -d : d;      
}    

Point2d Line2d::projection(const Point2d& q) const {
  Vector Vpq = q - p0;
  double lambda = dotProduct(Vpq, V)/dotProduct(V, V);
  return p0 + lambda * V;
}

double Line2d::projectionLambda(const Point2d& q) const {
  Vector Vpq = q - p0;
  Vector V1 = p1-p0;
  double lambda = dotProduct(Vpq, V1)/dotProduct(V1, V1);
  return lambda;
}

int Line2d::orientation( const Point2d& p ) const {
  double d2 = (p1.Y()-p0.Y())*(p0.X()-p.X());
  double d1 = (p1.X()-p0.X())*(p0.Y()-p.Y());
   if (d1 == d2)
      return 0;
   else
      return (d1 > d2) ? +1 : -1;   /// which is which
}

double Line2d::y_abs() const {
   if(isVertical()){
      //fprintf(stderr, "Error in computing abs for a vertical line",
        //	__FILE__, __LINE__, false);
      return 0;
   } else
      return (p0.Y() - slope()*p0.X());
}

// returning the tangent of the slope angle
double Line2d::slope() const 
{
   double dx = p0.X() - p1.X();
   double dy = p0.Y() - p1.Y();
   if(dx == 0){
     //fprintf(stderr, "Error in computing slope for a vertical line",
          //     __FILE__, __LINE__, false);
     return 0;
   } else
	return dy / dx;
}

 /*************************************************************
  *   intersection
 *************************************************************/

int Line2d::intersects(const Line2d& l) const{ 
  // return the dimension of intersection
  // return -1 if not intersect
  if( isParallel( l ) ) 
    return -1;
  else if( isCoincident( l ) ) 
    return 1;
  else 	// intersection is point
    return 0;    
}

GeomObj* Line2d::intersection(const Line2d& l) const {
   if( is_trivial() || l.is_trivial() ) 
     return NULL;

  if( isCoincident( l ) )
    return new Line2d(*this);

  if( isParallel( l ) )
    return NULL;

  double u = det(l.direction(), V);
  double w = det(l.direction(), l.startPt() - p0 );
  double lambda = w/u;
  Vector t = lambda * V;
  return new Point2d(p0 + t);
}

int orientation2d( Line2d& l, Point2d& p)
{ 
   Point2d start = l.startPt();
   Point2d end = l.stopPt(); 
  double d2 = (end.Y()-start.Y())*(start.X()-p.X());
  double d1 = (end.X()-start.X())*(start.Y()-p.Y());
   if (d1 == d2)
      return 0;
   else
      return (d1 > d2) ? +1 : -1;   /// which is which
}

std::ostream & operator<<(std::ostream &out, const Line2d &l) {
  return out << "Line2d[" << l.p0 << "===" << l.p1 << ";" << l.V << "]";
}

std::istream& operator>>(std::istream& in, Line2d& l) {
   // syntax: {[} p {===} q {]}

   Point2d p, q;
   char c;

   do in.get(c); while (isspace(c));
   if (c != '[') in.putback(c);

   in >> p;

   do in.get(c); while (isspace(c));
   while (c == '=') in.get(c);
   while (isspace(c)) in.get(c);
   in.putback(c);

   in >> q;

   do in.get(c); while (c == ' ');
   if (c != ']') in.putback(c);

   l = Line2d(p, q);
   return in;
}

Line2d p_bisector(const Point2d& p, const Point2d& q) { 
   return Line2d(p, q).rotate90(midPoint(p, q));
}
