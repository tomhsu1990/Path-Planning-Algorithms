
/*****************************************************************
 * File: segment2d.cc
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
 * $Id: segment2d.cpp,v 1.1 2006/04/03 18:55:31 exact Exp $
 *****************************************************************/

#include <CORE/geom2d/segment2d.h>

/************************************************************
 * *   constructors
************************************************************/
 //ray segment
Segment2d::Segment2d(const Point2d &p, const Vector &v) : p0(p), p1(p + v) { 
  directed = false;      
  open = false;          
}

Segment2d::Segment2d(const Point2d &_p0, const Point2d &_p1) : p0(_p0), p1(_p1) { 
  directed = false; 
  open = false;     
}

Segment2d::Segment2d(const Segment2d &s) : p0(s.startPt()), p1(s.stopPt()) { 
  directed = false; 
  open = false;     
}

Segment2d::Segment2d() : p0(0.0, 0.0), p1(0.0, 0.0) {
  directed = false; 
  open = false;     
}
 // line passes through the origin with direction 0.

/************************************************************
 * *   Member functions 
************************************************************/

  // test if two segments are coincident
bool Segment2d::isCoincident( const Segment2d& s) const {
  if( directed )
    return (p0==s.startPt() && p1==s.stopPt());
  else if( (p0==s.startPt() && p1==s.stopPt()) || (p1==s.startPt() && p0==s.stopPt()) )
    return true;
  else 
    return false;
}
    
 // returns the Euclidean distance between this segment and point q
 // problem exists when segment is open
double Segment2d::distance( const Point2d& p ) const{
  double s1 = dotProduct(p1-p0, p-p0);
  double s2 = dotProduct(p1-p0, p-p1);
  if( s1*s2 <= 0)
    return toLine().distance( p );
  else
    return (s1>0)? p1.distance( p ) : p0.distance( p );
}

 // returns the point on segment closest to q
 /* NOTE:
  * problem exist when segment is open!
  * need to be solved later
  */
Point2d Segment2d::nearPt( const Point2d& p ) const {
  double s1 = dotProduct(p1-p0, p-p0);
  double s2 = dotProduct(p1-p0, p-p1);
  if( s1*s2 <= 0)
    return toLine().projection( p );
  else
    return (s1>0)? p1 : p0;
}

bool Segment2d::onSegment(Point2d r)
{
    if (r.X() <= std::max(p0.X(), p1.X()) && r.X() >= std::min(p0.X(), p1.X()) &&
        r.Y() <= std::max(p0.Y(), p1.Y()) && r.Y() >= std::min(p0.Y(), p1.Y()))
        return true;
    return false;
}

bool Segment2d::contains( const Point2d& p ) const {
  if( !toLine().contains( p ) )
    return false;
   
   // must be collinear here
  double sign = dotProduct(p-p0, p-p1);
  if( sign < 0 )
    return true;
  else if( sign==0 && !open )
    return true;
  else 
    return false;
}

 //decides whether this segment intersects t 
int Segment2d::intersects( const Line2d& l) const {  
  if( l.contains(p0) && l.contains(p1) ) // this segment is on l
    return 1;

  int s1 = l.orientation( p0 );
  int s2 = l.orientation( p1 );
  
  if(s1==0 && s2==0)  //is on
    return 1;
  else if( s1*s2>0 )
    return -1;
  else
    return 0;
}
  

  //return dim of intersection
int Segment2d::intersects( const Segment2d& s) const {
  int res = intersects( s.toLine());
  if( res == -1 )
    return -1;
  else if( res == 0 )
    if( s.intersects( toLine() ) == 0 )       return  0; //double check
    else                                      return -1;
  else {	// must be collinear
    double s1 = dotProduct( s.startPt()-p0, s.startPt()-p1 );  // o<-----*-->o
    double s2 = dotProduct( s.stopPt()-p0, s.stopPt()-p1 );  
// o<---------o<--*
    double s3 = dotProduct( p0-s.startPt(), p0-s.stopPt() );   // *<-----o-->*
    double s4 = dotProduct( p1-s.startPt(), p1-s.stopPt() );   // *<---------*<--o

    if( s1<0 || s2<0 || s3<0 || s4<0 ) //  -o---*--o----*-  OR: -o---*-----*-o-
                                       //one end point is contained
      return 1;                        // then must be a segment intersection
    else if( s1==0 && s2==0 )      // -8===========8-
      return 1;                    // coincident
    else if( s1>0 && s2>0 )    
      return -1;   //  -o--------o-*----*- or  *----*--o--------o-
    else     // one end point is coincident        -o--------8-------o-
      if( open || s.isOpen() ) return -1;
      else        return  0;
  }
}

  // return intersection point if this segment and l intersect at a single point
  // the intersection point is returned 
GeomObj* Segment2d::intersection( const Line2d& l ) const {
  if( l.isCoincident( toLine() ) )
    return new Segment2d(*this);
    
  Point2d* p = (Point2d*)l.intersection( toLine() );
  if( p!=NULL && contains( *p ) )  return p;
  else              return  NULL; 
}

  // return intersection -- segment or point
GeomObj* Segment2d::intersection( const Segment2d& s ) const {
  
  if( toLine().isCoincident( s.toLine() ) ) {
     // colinear 
     // first intersect the segment s with ray starting from p0 
    Segment2d interSegment;
    double s1 = dotProduct(s.startPt() - p0, p1-p0 );
    double s2 = dotProduct(s.stopPt()  - p0, p1-p0 );
    if( s1<0 && s2<0 )
      return NULL;
    else if(s1>=0 && s2>=0)  // copy s
      interSegment = Segment2d( s );
    else if (s1>=0)
      interSegment = Segment2d(p0, s.startPt());
    else  
      interSegment = Segment2d(p0, s.stopPt());

     // then intersect the resulting segment with ray starting from p1
    double s3 = dotProduct(interSegment.startPt() - p1, p1-p0 );
    double s4 = dotProduct(interSegment.stopPt()  - p1, p1-p0 );
    if( s3>0 && s4>0 )
      return NULL;
    else if(s3<=0 && s4>0)
      interSegment = Segment2d(interSegment.startPt(), p1);
    else if(s4<=0 && s3>0)
      interSegment = Segment2d(interSegment.stopPt() , p1);
    //else -- s3<=0 && s4<=0
    //do nothing cause intersection remains unchanged

    if( interSegment.isTrivial() )  //intersection is point
      return new Point2d( interSegment.startPt() );
    else 
      return new Segment2d( interSegment );
  }
  else{  // use line-line intersection
    Point2d* p = (Point2d *)toLine().intersection( s.toLine() );
    if( p!=NULL && contains(*p) && s.contains(*p) )    // check
      return p;
    else 
      return NULL;
  }

}


std::ostream& operator<<(std::ostream& out, const Segment2d & s) {
  return out << "Segment2d[" << s.p0 << "," << s.p1 << "; directed=" << s.directed 
  		   << "; open=" << s.open << "]";
}

std::istream& operator>>(std::istream& in, Segment2d& l) {
  Point2d pStart, pEnd;
  std::cout << "\nInput segment start point(-,-):";
  in >> pStart;
  std::cout << "Input segment end point(-,-):";
  in >> pEnd;
  l = Segment2d(pStart, pEnd); 
  return in;
}

