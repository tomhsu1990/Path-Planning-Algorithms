
/*****************************************************************
 * File: point3d.cc
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
 * $Id: point3d.cpp,v 1.1 2006/04/03 18:55:31 exact Exp $
 *****************************************************************/

#include <CORE/geom3d/point3d.h>
#include <CORE/geom3d/plane3d.h>
#include <CORE/geom3d/segment3d.h>

// midPt(p, q) returns (p+q)/2:
Point3d midPt3d ( Point3d& a, Point3d& b)
{
  return Point3d((a.X()+b.X())/2,(a.Y()+b.Y())/2, (a.Z()+b.Z())/2);
}

/* returns true if points a, b, c and d are coplanar, i.e.,
 * orientation(a, b, c, d) = 0, and false otherwise. 
 */
bool coplanar(const Point3d& a, const Point3d& b, 
                     const Point3d& c, const Point3d& d)
{
   return ( orientation3d(a, b, c, d) == 0 );
}

/************************************************************
 *  CONSTRUCTORS 
 ************************************************************/

Point3d::Point3d() : x(0), y(0), z(0) { }

Point3d::Point3d(double _x, double _y, double _z) : x(_x), y(_y), z(_z) { }

Point3d::Point3d(const Point3d& p) : x(p.x), y(p.y), z(p.z) { }

Point3d::Point3d(const Vector& v) :  x(v[0]), y(v[1]), z(v[2]) { 
  if( v.dimension() > 3 ) throw Vector::RangeException();
}

/************************************************************
 *  METHODS
 ************************************************************/

Point3d& Point3d::operator=(const Point3d& p)
{
   x = p.x; y = p.y; z = p.z; return *this;
}

Point3d Point3d::operator+(const Point3d& p) const {
  return Point3d(x + p.x, y + p.y, z + p.z);
}

//Point3d Point3d::operator-(const Point3d& p) const {
//  return Point3d(x - p.x, y - p.y, z - p.z);
//}

Point3d Point3d::operator*(const double& a) const {
  return Point3d(x*a, y*a, z*a);
}

Point3d Point3d::operator+(const Vector& v) const {
  return Point3d(x + v[0], y + v[1], z + v[2]);
}

Point3d Point3d::operator-(const Vector& v) const {
  return Point3d(x - v[0], y - v[1], z - v[2]);
}

Vector Point3d::operator-(const Point3d &p) const {
  return Vector(x - p.x, y - p.y, z - p.z);
}

double Point3d::distance(const Point3d& p) const
{
   double dx = x - p.x;
   double dy = y - p.y;
   double dz = z - p.z;
   return sqrt(dx*dx + dy*dy + dz*dz);
}

bool Point3d::operator==(const Point3d& p) const
{
   return (x == p.x) && (y == p.y) && (z == p.z);
}

   // compute signed volume of a tetrahedron
double signed_volume(const Point3d& a, const Point3d& b, 
                     const Point3d& c, const Point3d& d) {
  double  a11 = a.x-d.x;
  double  a12 = a.y-d.y;
  double  a13 = a.z-d.z;
  double  a21 = b.x-d.x;
  double  a22 = b.y-d.y;
  double  a23 = b.z-d.z;
  double  a31 = c.x-d.x;
  double  a32 = c.y-d.y;
  double  a33 = c.z-d.z;

  double s = a11*(a22*a33-a23*a32) - a12*(a21*a33-a23*a31) + a13*(a21*a32-a22*a31);
  return s;
}

/* orientation3d(a, b, c, d) 
 *   computes the orientation of points a, b, c, d as the sign
 *   of the determinant
 *              | ax  ay  az 1 |
 *              | bx  by  bz 1 |
 *              | cx  cy  cz 1 |
 *              | dx  dy  dz 1 |
 *   i.e., it returns +1 if d lies in the opposite side w.r.t. the 
 *   counter-clockwise side of plane formed by a, b, c
 */
int orientation3d(const Point3d& a, const Point3d& b, 
                  const Point3d& c, const Point3d& d)
{
  double s = signed_volume(a, b, c, d);
  if( s == 0 ) return 0;
  else         return (s>0)? 1 : -1;
}


/* area(a, b, c) returns 1/2 times the determinant of orientation(a,b,c)
 * above.  This is the signed area of the triangle determined by a, b, c,
 * positive if orientation(a,b,c) > 0, and negative otherwise.  */

double volume(const Point3d& a, const Point3d& b, 
                     const Point3d& c, const Point3d& d)
{
  double signed_vol = signed_volume(a, b, c, d);
  return (signed_vol > 0)? signed_vol : -signed_vol;
}

bool Point3d::insideMultiHalfspace( Plane3d* H[], int H_n ){
  // test whether the corner is inside the n half space
  for(int i=0;i<H_n;++i){
    // apply the point to the plane equation
    if(H[i]->apply(Point3d(x, y, z)) > 0){
      return false;
    }
  }
  return true;
}

double Point3d::separationL( const Segment3d& s ) const {
  Point3d p0 = s.startPt();
  Point3d p1 = s.stopPt();
  Point3d p(x,y,z);
  double s1 = dotProduct(p1-p0, p-p0);
  double s2 = dotProduct(p1-p0, p-p1);
  if( s1*s2 <= 0)
    return s.toLine().distance(p);
  else
    return (s1>0)?p1.distance(p) : p0.distance(p);
}

std::ostream &operator<<(std::ostream &o, const Point3d& p) {
  return o << "Point(" << p.x << "," << p.y << "," << p.z << ")";
}

std::istream& operator>>(std::istream& in, Point3d& p)
{
  double x, y, z;  
  char c;

  do in.get(c); while (in && isspace(c));

  if (!in) return in;

  if (c != '(') in.putback(c);
  in >> x;

  do in.get(c); while (isspace(c));
  if (c != ',') in.putback(c);
  in >> y;

  do in.get(c); while (isspace(c));
  if (c != ',') in.putback(c);
  in >> z; 

  do in.get(c); while (c == ' ');
  if (c != ')') in.putback(c);

  p = Point3d(x, y, z);
  return in;
}



