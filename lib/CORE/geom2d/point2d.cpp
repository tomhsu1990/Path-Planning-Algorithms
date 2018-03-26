/*****************************************************************
 * File: point2d.ccp
 * Synopsis:
 *      Basic 2-dimensional geometry
 * Author: Shubin Zhao (shubinz@cs.nyu.edu), 2001.
 * 	   Chee Yap (yap@cs.nyu.edu), 2002.
 *
 * To Do:
 *    Generation of random and special point sets
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
 * $Id: point2d.cpp,v 1.1 2006/04/03 18:55:31 exact Exp $
 *****************************************************************/

#include "CORE/CORE.h"
#include "CORE/geom2d/point2d.h"

// returns the midpoint between a and b
Point2d midPoint(const Point2d& a, const Point2d& b)
{ return Point2d( (a.X() + b.X())/2 , (a.Y() + b.Y())/2 );
}

// returns "asymmetric center" point on the line through a, b.
// default: alpha = 0.5 (the midpoint or center)
Point2d aCenter(const Point2d& a, const Point2d& b, machine_double alpha)
{ return Point2d( (1 - alpha) * a.X() + alpha * b.X() , 
		  (1 - alpha) * a.Y() + alpha * b.Y() );
}

/* orientation: computes the orientation of points a, b, and c as the sign
                of the determinant
                | ax  ay  1 |
                | bx  by  1 |
                | cx  cy  1 |
             i.e., it returns +1 if point c lies left of the directed line
                  through a and b, 0 if a,b,and c are collinear, and
                  -1 otherwise.
 */
int orientation2d(const Point2d& a, const Point2d& b, const Point2d& c)
{
   double d1 = (a.X() - b.X()) * (a.Y() - c.Y());
   double d2 = (a.Y() - b.Y()) * (a.X() - c.X());
   if (d1 == d2)
      return 0;
   else
      return (d1 > d2) ? +1 : -1;
}

bool leftTurn(const Point2d& a, const Point2d& b, const Point2d& c)
{ return (orientation2d(a,b,c)>0); }

bool rightTurn(const Point2d& a, const Point2d& b, const Point2d& c)
{ return (orientation2d(a,b,c)<0); }

double area(const Point2d& a, const Point2d& b, const Point2d& c)
/* computes the signed area of the triangle determined by a, b, c,
   positive if orientation(a,b,c) > 0, and negative otherwise.  */
{
   return ((a.X()-b.X()) * (a.Y()-c.Y()) -
           (a.Y()-b.Y()) * (a.X()-c.X()));
}

bool collinear(const Point2d& a, const Point2d& b, const Point2d& c)
	/* returns true if points a, b, and c are collinear, i.e.,
   	orientation(a, b, c) = 0, and false otherwise. */
{
   return (a.Y()-b.Y()) * (a.X()-c.X()) ==
          (a.X()-b.X()) * (a.Y()-c.Y());
}

bool between(const Point2d& a, const Point2d& b, const Point2d& c)
	/* returns true if a,b,c are collinear and b is strictly
	 * between a and c. */
{
	return ((collinear(a,b,c) && dotProduct(a-b, c-b)<0));
}

bool betweenVar(const Point2d& a, const Point2d& b, const Point2d& c)
	/* returns true if dotProduct(a-b,c-b)<0;
	 * If a,b,c are collinear, this is equivalent
	 * to b lying strictly between a and c. */
{
	return (dotProduct(a-b, c-b)<0);
}

Point2d::Point2d() : x(0.0), y(0.0) {}
Point2d::Point2d(double _x, double _y) : x(_x), y(_y) {}
Point2d::Point2d(const Point2d &p) : x(p.x), y(p.y) {}
Point2d::Point2d(Vector v) : x(v[0]), y(v[1]) { }

Point2d& Point2d::operator=(const Point2d& p)
{
   x = p.x;
   y = p.y;
   return *this;
}

double Point2d::distance(const Point2d p) const
{
   double dx = x - p.x;
   double dy = y - p.y;
   return sqrt(dx*dx + dy*dy);
}

Vector Point2d::operator-(const Point2d &p) const
{
   return Vector(x - p.x, y-p.y);
}

Point2d Point2d::operator+(const Vector &v) const
{
   return Point2d(x+v[0], y+v[1]);
}

Point2d Point2d::rotate90(const Point2d& p)
{
   double px = p.X();
   double py = p.Y();
   double dx = X() - px;
   double dy = Y() - py;
   return Point2d(px-dy, py+dx);
}

bool Point2d::operator==(const Point2d &p) const
{
   return (x == p.x) && (y == p.y);
}

std::ostream &operator<<(std::ostream &out, const Point2d p) {
  out << "Point2d(" << p.x << "," << p.y << ")";
  return out;
}

/* OLD: 
std::istream& operator>>(std::istream& in, Point2d& p)
	// Format: ( nnn , nnn )
	// where the white spaces, "(", "," and ")" are all optional.
	// but no other non-white chars may be used.  
{
  double x, y;  // CORE disallow in >> double????
  char c;
  // char buf[256];

  do in.get(c); while (in && isspace(c));

  if (!in) return in;

  if (c != '(') in.putback(c);

  in >> x;

  do in.get(c); while (isspace(c));
  if (c != ',') in.putback(c);

  in >> y;

  do in.get(c); while (c == ' ');
  if (c != ')') in.putback(c);

  p = Point2d(x, y);
  return in;
}
*/

// NEW version of operator>> uses the auxilliary functions:
// 	startNum(char),
// 	getToChar(istream, char),
// 	getToNum(istream, char).

// startNum(c): check if char c can be the start of a literal number
  // i.e., if c occurs in string '01234567890+-.'
bool startNum(char c) {
  switch(c) {
    case '0': ;  case '1': ;  case '2': ;  case '3': ;  case '4': ;  
    case '5': ;  case '6': ;  case '7': ;  case '8': ;  case '9': ; 
    case '+': ;  case '-': ; 
    case '.': return true; break;
  }
  return false;
}//startNum(c)

// getToChar(istream, mark):
// 	-- moves input cursor to the next occurrence of mark
// 	-- in any case, it will stop if it gets to the beginning
// 	   of the next number (indicated by startNum()).
bool getToChar( std::istream& in, char mark) {
  char c;
  char buf[256];
  bool got_mark = false;
  do {
	in.get(c); 
	if (!isspace(c)) {
	     if (c == mark) {
		got_mark = true; // this is what we were looking for
	     } else if (c == '#') {	// rest of line is discarded as comment
		in.getline(buf,256);
	     } else if (startNum(c)) {
		     //IOErrorFlag = 2;	// I/O Error 2:  did not find
					// mark before next number
		in.putback(c);	// we accept this even if mark not yet found
		return false;		// Usually not considered error
	     } else {
		     //IOErrorFlag = 1;	// I/O Error 1: 
		return false;		// Unexpected character
	     }
	}
  } while (!(got_mark || in.eof()));
  return (got_mark);
}//getToChar

// getToNum(istream, mark, strict):
// 	-- moves input cursor to the beginning of the next number. 
bool getToNum( std::istream& in, char mark, bool strict=false) {
	// default value for strict is false
  // NOTES: We ignore spaces, comment lines
  // and accept at most one copy of the "mark" character.
  // [It must be exactly one copy when strict=true]
  // If so, we return true.
  // Return false if reach eof with no number, or
  // saw 2 "mark" characters, or saw ANY other unexpected chars.
  char c;
  char buf[256];
  bool got_mark = false;
  bool got_num = false;

  do {
	in.get(c); 
	if (!isspace(c)) {
	     if (c == mark) {
		if (got_mark) {
			//IOErrorFlag = 2;	// I/O Error 2:
			return false; 		// Second mark found
		} else got_mark = true; // this is what we were looking for
	     } else if (startNum(c)) {	// we may get this before finding mark
		in.putback(c);	// but we accept this (not considered error)
		got_num = true;
	     } else if (c == '#') {	// rest of line is discarded as comment
		in.getline(buf,256);
	     } else {
		     //IOErrorFlag = 1;	// I/O Error 1: 
		return false;		// Unexpected character
	     }
	}
  } while (!(got_num || in.eof()));
  return (got_num && (got_mark || !strict));
}//getToNum


// Operator>>>(in, p) :  Reads an input point which has the format
// 			( nnn , nnn )
// where the white spaces, "(", "," and ")" are all optional.
std::istream& operator>>(std::istream& in, Point2d& p)
{
	// NOTE: No other non-white chars may be used.  Whenever we see #, 
	// the rest of line is omitted, and it counts as white space.
	// E.g., "[ nnn : nnn ]" is no good but "( nnn, nnn) #xxx" is OK.
  double x, y;  

  // if (!(getToNum(in, '(') || IOErrorFlag || in.eof())) 
  if (!(getToNum(in, '(') || in.eof())) 
	return in;		// No first value
  else
  	in >> x;
  // if (!(getToNum(in, ',') ||  IOErrorFlag ||in.eof())) 
  if (!(getToNum(in, ',') || in.eof())) 
	return in;		// no second value
  else
	in >> y;
  getToChar(in, ')');	// This terminates if the start of the next number
  p = Point2d(x, y);	// 	occurs before ')'.
  return in;
} // operator>>


// readPoints(iS, pA, MaxN, N) :
//    reads a sequence of points from istream iS into Point2d array pA.
int readPoints(std::istream & iS, Point2d *pA, int MaxN=1000, int N=0){
  // DEFAULT VALUES: MaxN=1000, N=0
  // The input stream constains a sequence of 2K+1 numbers of the form
  //
  //     [NN]   ( x1 ,  y1 )  ( x2 , y2 )  ... ( xK , yK )
  //
  // The i-th point is (xi, yi). 
  //    (0) NN is optional if N is given as argument (then N is set to NN)
  //    (1) Any of the "(", "," and ")" are optional
  //    (2) Newlines, extra white spaces, '#' are all ignored.
  //    (3) Also, everything after '#' is discarded.
  // If N > MaxN, nothing is read and 0 is returned.  
  // Returns the number of points actually read, i.e, min(K, N).
  
  int i;
  // IOErrorFlag = 0;		// initially no error
  if (N == 0) {			// read NN from input
     if (getToNum(iS, ' '))
  	iS >> N;
     else
  	return 0;
  }
  if (N > MaxN) return 0;
  for (i = 0; i < N; i++){
	  if (iS.eof()) break;
	  // if (!IOErrorFlag) {
	  if (true) {
	  	iS >> pA[i];
	  } else {
		i = 0;
	        break;
	  }
  }
  return i;  // number of points read
}//readPoints(...)

