/*****************************************************************
 * File: plane3d.cc
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
 * $Id: plane3d.cpp,v 1.2 2006/04/03 19:38:35 exact Exp $
 *****************************************************************/

#include <CORE/geom3d/segment3d.h>
#include <CORE/geom3d/plane3d.h>

/************************************************************
 * *  constructors
 ************************************************************/

  //copy constructor
Plane3d::Plane3d(const Plane3d & pl) : a(pl.A()), b(pl.B()), c(pl.C()), 
  d(pl.displacement()), n( pl.normal() ){
}  

 // plane with direction v passes through point p 
Plane3d::Plane3d(const Point3d & p, const Vector &v) : 
  a(v[0]), b(v[1]), c(v[2]), d(0.0), n(v) {
  d = - dotProduct(v, p.toVector());
}

 //plane passes through points p1, p2, p3 
Plane3d::Plane3d(const Point3d &p1, const Point3d &p2, const Point3d &p3) {
  Vector V31 = p1 - p3;
  Vector V32 = p2 - p3;
  n = V31.cross( V32 ); 
  a = n[0]; b=n[1]; c=n[2];
  d = - dotProduct(n, p3.toVector());
}

 //plane passes through point p and line l
Plane3d::Plane3d(const Point3d &p, const Line3d &l) {
  Plane3d(p, l.startPt(), l.stopPt() );
}

 //plane passes through point p and segment s
Plane3d::Plane3d(const Point3d &p, const Segment3d &s) {
  Point3d p1 = s.startPt();
  Point3d p2 = s.stopPt();
  Plane3d(p, p1, p2 );
}
 
 // plane determined by vector and displacement
Plane3d::Plane3d(const Vector& v1, double d1) : 
  a(v1[0]), b(v1[1]), c(v1[2]), d(d1), n(v1) {
}

 // plane determined by equation
Plane3d::Plane3d(double a1, double b1, double c1, double d1)
    : a(a1), b(b1), c(c1), d(d1), n(a1, b1, c1) { }

/************************************************************
 * *  member functions
 ************************************************************/

double* Plane3d::coeffients() const {

 double* params = new double[4];
  params[0] = a;
  params[0] = b;
  params[0] = c;
  params[0] = d;
  
  return params;
}
  // test coincidence
bool Plane3d::isCoincident(const Plane3d& pl) const {
  Vector Vcross = n.cross( pl.normal() );
  if( ( !Vcross.isZero() ) && (d == pl.displacement()) )
    return true;
  else
    return false;
}

bool Plane3d::isParallel(const Plane3d& pl) const {
        return n.cross( pl.normal() ).isZero();
}

bool Plane3d::isParallel(const Line3d& l) const {
        return dotProduct( n, l.direction() ) == 0;
}


bool Plane3d::contains( const Point3d& p ) const {
  return apply(p) == 0 ;
}

bool Plane3d::contains( const Line3d& l ) const {
  return ( contains( l.startPt() ) && contains( l.stopPt() ) );
}

bool Plane3d::contains( const Segment3d& s ) const {
  return ( contains( s.startPt() ) && contains( s.stopPt() ) );
}

Point3d Plane3d::projection(const Point3d& p) const {
  double lambda = apply(p) / dotProduct(n, n);
  return (p - n * lambda );
}

Line3d Plane3d::projection(const Line3d& l) const {
  return Line3d( projection( l.startPt() ), projection( l.stopPt() ) );
}

Segment3d Plane3d::projection(const Segment3d& s) const {
        return Segment3d( projection( s.startPt() ), projection( s.stopPt() ) );
}

  //distance
  // d = A*px + B*px + c*pz + d/
double Plane3d::distance( const Point3d& p ) const {
  double dist = apply(p);
  return fabs(dist)/n.norm();
}
  
double Plane3d::distance( const Line3d& l ) const {
  if( intersects( l ) ) return 0;
  else  //line l must be paralell with this plane
    return distance( l.startPt() );
}

double Plane3d::distance( const Segment3d& s ) const {
  if( intersects( s ) ) return 0.0;
  else   // minimum distance must be at end points
    return std::min( distance( s.startPt() ), distance( s.stopPt() ) );
}

 // test if the line is paralell to the plane
 // return dim of intersection
int Plane3d::intersects( const Line3d& l ) const {
  if( contains( l ) ) 
    return 1;
  else 
    return (dotProduct(n, l.direction()) != 0.0)? 0 : -1;
}

 // test if two end points lie on different side of the plane
 // return dim of intersection
int Plane3d::intersects( const Segment3d& s ) const {
  if( contains(s) ) return 1;
  
  double dist1 = apply( s.startPt() );
  double dist2 = apply( s.stopPt()  );

  if( dist1 * dist2 > 0 )
    return -1;
  else 
    return 0;
}

int Plane3d::intersects( const Plane3d& pl ) const {
  if( isCoincident( pl ) ) return 2;
  
  Vector Vcross = n.cross( pl.normal() );
  if( !Vcross.isZero() )
    return 1;
  else 
      return -1;
}

  // need improvement
GeomObj* Plane3d::intersection( const Line3d& l ) const {  
  if( intersects( l ) == -1 ) 
    return NULL;

  Point3d interPt = l.startPt() + 
            (l.direction() * 
             ( - apply( l.startPt() ) / 
               (dotProduct(n, l.direction())) ) );
               
  return new Point3d(interPt);
}

  // need improvement
GeomObj* Plane3d::intersection( const Segment3d& s ) const {
  if( contains( s ) )    return new Segment3d(s);

  GeomObj* obj = intersection( s.toLine() );
  if( obj == NULL )
    return NULL;
  else {
    Point3d* p = (Point3d *)obj;
    if( s.contains( *p ) )
      return obj;
    else
      return NULL;
  }
}


// we know the direction of intersection line is cross product of two normals
 // we just need to get an intersection point
GeomObj* Plane3d::intersection( const Plane3d& pl ) const {
  if( isParallel( pl ) ){
    return NULL;
  }
  else {  // intersection must be a line
    Vector l_direction = n.cross( pl.normal() ); //direction of intersection line
    Point3d p1;  // projection of the origin on this plane
    Point3d p2;  // projection of the origin on plane pl
    if( d == 0 )                 
      p1 = ORIGIN_3D; 
    else {
      double K = - dotProduct(n, n)/d;
      p1 = Point3d(a/K, b/K, c/K);
    }
      // compute projection of origin on pl
    if( pl.displacement() == 0 ) 
      p2 = ORIGIN_3D;
    else {
      double K = - dotProduct(pl.normal(), pl.normal())/pl.displacement();
      p2 = Point3d(pl.A()/K, pl.B()/K, pl.C()/K);
    }
        
    if( p1 == p2 )
      return new Line3d(p1, l_direction);

    Line3d line_on_this( p1, projection( p2 ));
    Line3d line_on_this2( pl.projection( p1 ), p2 );

    Point3d* pInterPt  = (Point3d *)pl.intersection( line_on_this );
    Point3d* pInterPt2 = (Point3d *)intersection( line_on_this2 );

    if(pInterPt != NULL)
      return new Line3d(*pInterPt, l_direction);
    else
      return new Line3d(*pInterPt2, l_direction);
  }
}

std::ostream& operator<<(std::ostream& o, const Plane3d& pl) {
  return o << "Plane3d [normal=" << pl.n << ", displacement=" << pl.d << "]";
}
 
std::istream& operator>>(std::istream& in, Plane3d& pl) {
  Point3d p;    
  std::cout << "\nInput normal of plane(-,-,-):";
  in >> p;  // point is in the same format as vector in 3d
  
  double disp;
  std::cout << "Input displacement of plane:";
  in >> disp;
  
  pl = Plane3d( p.toVector(), disp );
  return in;
}

