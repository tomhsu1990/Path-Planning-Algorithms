#ifndef _POLYGON2D_H_
#define _POLYGON2D_H_

#include <CORE/geom2d/point2d.h>
#include <CORE/geom2d/line2d.h>
#include <CORE/geom2d/segment2d.h>
#include <CORE/geom2d/triangle2d.h>
#include <CORE/linearAlgebra.h>

class Polygon2d : public GeomObj {

   /* An instance C of the data type polygon
      in the plane passing through points pts.
    */

private:
 
   std::vector<Point2d> pts;  // the points defining the polygon

   Point2d cp;
   int orient; // clockwise: 1
               // counter-clockwise: -1

public:

   Polygon2d(const std::vector<Point2d>& p);
   //initialized to the oriented polygon through points p

   Polygon2d(const Point2d& p);
   //initialized to the trivial polygon with center p

   Polygon2d();
   //initialized to the trivial polygon with center (0,0)

   Polygon2d(const Polygon2d& poly);
   //copy constructor

   virtual ~Polygon2d(); 

   Polygon2d& operator=(const Polygon2d& poly);

   //operations

   Point2d center();
   //return the center of the polygon

   void addPoint(Point2d pt) { pts.push_back(pt); }
   void addPoint(double x, double y) { pts.push_back(Point2d(x, y)); }
   std::vector<Point2d> points() const { return pts; }
   Point2d pointp(unsigned p) const { return pts[p]; }

   //bool is_degerate() const { return orient == 0; }
   //returns true if the defining points are collinear

   int orientation() const { return orient; }
   void setOrientation() { orient=clockwise(); }

   int side_of(const Point2d& p) const; 
   // returns -1, +1 or 0 if p lies right of, left of or on the circle
   // respectively

   bool inside(const Point2d& p);
   //returns true if p lies inside of the polygon

   bool outside(const Point2d& p);
   
   bool contains(const Point2d& p) const ;
   //returns true if p lies on the polygon, false otherwise

   double distance(const Point2d& p);
   //returns the distance between p and the polygon

   double distance(const Line2d& l);
   //returns the distance between l and the polygon
   //distance from center to l minus radius

   double distance(Polygon2d& P);
   //returns the distance between this polygon and polygon P

   bool operator==(const Polygon2d& D) const ;
   bool operator!=(const Polygon2d& D) const 
	{ return !operator==(D); }

   friend std::ostream& operator<<(std::ostream& out, Polygon2d& poly);
   friend std::istream& operator>>(std::istream& in, Polygon2d poly); //?? Polygon2d &

   int clockwise();

   // compute area of a contour/polygon
   float area(const std::vector<Point2d> &contour);

   // triangulate a contour/polygon, places results in STL vector
   // as series of triangles.
   bool processTriangulate(const std::vector<Point2d> &contour, std::vector<Point2d> &result);

   bool snip(const std::vector<Point2d> &contour, int u, int v, int w, int n, int *V);

}; // class Polygon2d


#endif
