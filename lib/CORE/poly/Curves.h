/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 * 	$Source: /home/exact/cvsroot/exact/corelib2/inc/CORE/poly/Curves.h,v $
 * 	$Revision: 1.21 $ $Date: 2010/11/08 14:52:09 $
 *
 * This file is part of CORE (http://cs.nyu.edu/exact/core/); you may
 * redistribute it under the terms of the Q Public License version 1.0.
 * See the file LICENSE.QPL distributed with CORE.
 *
 * Licensees holding a valid commercial license may use this file in
 * accordance with the commercial license agreement provided with the
 * software.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * File: Curves.h
 *
 * Description: 
 * 	Two templated classes are defined here:
 *		Curve and BiPoly
 *	These classes are parametrized by the number type
 *		(called NT) which represents the
 *		domain of the coefficients of the underlying
 *		polynomials.  Standard default is NT=BigInt, but
 *		we will allow NT=int, NT=BigRat, NT=BigFloat, NT=Expr.
 *	BiPoly represents the class of bivariate polynomials,
 *		i.e.,  each BiPoly object is an element of NT[X,Y].
 *		We store each BiPoly as a list of polynomials in X.
 *	Curve represents the class of plane curves whose equation
 *		is A(X,Y)=0, for some BiPoly A(X,Y).
 *	Features:
 *		--Constructor from strings such as
 *			"3 x^2 + 7 xy^2 - 4 x + 13"
 *			"3 x^2 + 7 xy^2 = 4 x - 13"  (has an equality sign)
 *			"(x -1)^5 (y+2)^2"	     (parenthesis)
 *
 *		  Mar 2009: Originally, only integer coefficients assumed.
 *		  	We now generalize this function so that
 *		        the coefficients can now accept decimal numbers, e.g.
 *				"0.3 x ( y^2 - 1.7 xy)"
 *		--toString() methods that is the inverse of the
 *		  	above string constructors	(Nov 2010)
 *			
 *		--Basic plot functions
 *
 *	To Do:
 *	  --String constructor generalizations:
 *	  	(1) allow one "=" sign (e.g., "3 x^2 = y^2 - xy")(DONE)
 *		(2) allow general parenthesis(DONE)
 *		(3) allow X and Y (DONE)
 *	  --We should be able to read/write
 *	  	curve definitions from/to files
 *	  --Plot should be more efficient (use previous roots
 *	  	to help find the next roots, there should be
 *	  	a "plot structure" that is persistent)
 *	  --Plot should refine in both x- and y-increments.
 *	  --Plot should have some option to show the
 *	  	x- and y-axes, and to label some points.
 *	  --verticalIntersect(...) should be implemented using
 *	        Polynomial<BigFloat>, not Polynomial<Expr> for efficiency
 *	  --the plot parameters (eps,xmin,xmax,ymin,ymax) must be
 *	        made part of the Curve class (static members).
 *	        Incorporate the "setParams" method into class.
 *	  --We should allow BiPoly to be representated as a poly in X, with
 *	        coeffs. being polynomials in Y.  All we need is to store
 *	        a boolean flag to check which rep is used. This allow some
 *	        operations to be more efficient.
 *
 *  Author:  Vikram Sharma and Chee Yap (April 12, 2004)
 *  	     Jihun (Dec 2006), Philip Davidson (Dec 2006)
 *  
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 ***************************************************************************/


#ifndef CORE_CURVES_H
#define CORE_CURVES_H

#include <fstream>
#include <string>
#include <list>
#include <iostream>
#include "CORE/poly/Poly.h"

// NOTE: This IntervalT.h is only used in eval1 and eval2 methods(see below)
// (Shang Wang, Aug 2011)
#include "CORE/IntervalT.h"
// pre-define
//template <class NT>
//class IntervalT;

// NOTE: at the end of this file, we include "Curves.tcc" which has
//   is based on Sturm sequences.  An alternative is to include 
//   "CurvesDesc.tcc" which is based on Descartes method.
//   The parser for bivariate polynomials is implemented in "Curves.tcc".


using namespace std;

CORE_BEGIN_NAMESPACE

// ==================================================
// Curve Class
// ==================================================

//typedef BigInt NT;
//typedef Expr   NT;
//typedef Polynomial<NT>        PolyNT;
//typedef std::vector<Expr>	VecExpr;
//typedef std::vector<BigInt>	VecBigInt;
//typedef std::vector<NT>       VecNT;
//typedef std::vector<Polynomial<NT> >	VecPoly;

template <class NT>
class Monomial{
  //Helper class to store the coefficients for given x-deg and y-deg 
  //Used by string input routines
 public:
  NT coeff;
  int xdeg;
  int ydeg;

  Monomial(){
  }

  Monomial(NT& cf, int& dx, int& dy){
    coeff = cf;
    xdeg = dx;
    ydeg = dy;
  }

  void dump(){
    std::cout << coeff << "X^" << xdeg << " Y^" << ydeg;
  }
}; // Monomial Class


//	Class of Bivariate polynomials
//	Viewed as a polynomial in Y with
//	coefficients which are polynomials in X
template <class NT>
class BiPoly{
 private:
  //The following are used in the constructor from strings.
  //For more details see the related constructor.
  void constructFromString(string& s, char myX='x', char myY='y');
  void constructX(int n, BiPoly<NT>& P);
  void constructY(int n, BiPoly<NT>& P);
  int getnumber(const char* c, int i, unsigned int len, BiPoly<NT> & P);
  bool isint(char c);
  int getint(const char* c, int i, unsigned int len, int & n);
  int matchparen(const char* cstr, int start);
  int getbasicterm(string s, BiPoly<NT> & P);
  int getterm(string s, BiPoly<NT> & P);

 public:
  int ydeg; //Y-degree of the polynomial
  std::vector<Polynomial<NT> > coeffX; //vector of (1+ydeg) polynomials in X
	// If ydeg = d, then the polynomial is F(X,Y) =
        //   (Y^d * coeffX[d]) + (Y^{d-1} * coeffX[d-1]) +...+ (coeffX[0]).

  ////////////////////////////////////////////////////////
  //Constructors
  ////////////////////////////////////////////////////////

  ///BiPoly()
  BiPoly(); // zero polynomial

  ///BiPoly(n)
  ///  Should we make this constructor explicit?
  BiPoly(int n);// creates a BiPoly with nominal y-degree equal to n.

  ///BiPoly(NT)
  BiPoly( const NT& c ); 

  ///BiPoly(vp): construct from a vector vp of Polynomial<NT>
  BiPoly(std::vector<Polynomial<NT> >& vp); 

  ///BiPoly(n, ap): construct from an array ap of n polynomials
  /// Useful alternative to BiPoly(vp) where vp is vector of Polynomial<NT>
  BiPoly(int n, Polynomial<NT> * ap ); 

  ///BiPoly(p, flag=false):
  ///	if true, it converts polynomial p(X) into P(Y)
  /// 	if false, it creates the polynomial Y-p(X)
  ///   NOTE: this constructor is made "explicit" to avoid mysterious
  ///   bugs (as it happened to us).  Perhaps should make more of these
  ///   constructors explicit...
  explicit BiPoly(Polynomial<NT> p, bool flag=false);
  
  //BiPoly(deg, d[], C[]):
  //	Takes in a list of list of coefficients.
  //	Each cofficient list represents a polynomial in X
  //
  //  deg - ydeg of the bipoly
  //  d[] - array containing the degrees of each coefficient (i.e., X poly)
  //  C[] - list of coefficients, we use array d to select the
  //      coefficients
  BiPoly(int deg, int *d, NT *C);

  //BiPoly(String s, char myX, char myY)
  //  myX and myY are names of the two variables.
  //  Default values of myX and myY are 'x' and 'y'.
  //  The string s has the form "3 x^2 + 7 xy^2 - 4 x + 13"
  //
  //  For now, we assume no parentheses, * or =.
  
  BiPoly(const string& s, char myX='x', char myY='y');
  BiPoly(const char* s, char myX='x', char myY='y');

  // copy constructor
  BiPoly(const BiPoly<NT>&);

  //Destructor
  ~BiPoly();
  //Destructor helper
  void deleteCoeffX();


  ////////////////////////////////////////////////////////
  // METHODS
  ////////////////////////////////////////////////////////
  
  /// toString(xvar,yvar):
  /// 	Outputs a string representation of this bipoly,
  /// 	in a very natural form, starting from leading coefficients
  /// 	in Y.  This operation is the "inverse" of the bipoly constructor
  /// 	from string (inverse in the semantic sense -- the same bipoly is
  ///	constructed, though the string format might have changed. 
  ///	E.g.,  If bp is a BiPolyNT, then the following test will pass:
  ///		if (bp == BiPolyNT(bp.toString())) cout << "CORRECT" << endl;
  std::string toString(char xvar = 'x', char yvar = 'y') ;
  
  // filedump (msg, ofs, com, com2)
  // 	where msg, com, com2 are strings.
  // 	msg is an message and com, com2 are the strings
  // 	preceding each output line
  // 	(e.g., msg="BiVariate Polynomial"  and com=com2="# ")
  // This is called by the other dump functions
  void dump(std::ostream & os, std::string msg = "",
      std::string com="# ", std::string com2 = "# ") const;
  // dump(ofs, msg, com) -- dump to file
  //void dump(std::ofstream & ofs, std::string msg,
  //    std::string com="# ", std::string com2="# ") const;

  // dump(msg, com) -- dump to std output
  void dump(std::string msg="", std::string com="",
      std::string com2="") const;

  /*Cannot work with these two functions right now.
    BiPoly as per now can only handle BigInt and int since
    Expr cannot be handled by Polynomial class.*/
  
  // yPolynomial(x) 
  //   returns the polynomial (in Y) when we substitute X=x
  
  /* BiPoly<NT> yPolynomial(const Expr & x) {

    VecExpr vE;

    for (int i=0; i<= ydeg; i++) {
      vE.push_back(coeffX[i].eval(x));
    }
    
    return BiPoly<NT>(vE);
  }//yPolynomial
  */
  
  Polynomial<NT> yPolynomial(const NT & x);

  // Expr version of yPoly (temporary hack)
  Polynomial<Expr> yExprPolynomial(const Expr & x);

  // BF version of yPoly (temporary hack)
  Polynomial<BigFloat> yBFPolynomial(const BigFloat & x);

  // xPolynomial(y) 
  //   returns the polynomial (in X) when we substitute Y=y
  //   
  //   N.B. May need the
  //   		Polynomial<Expr> xExprPolynomial(Expr y)
  //   version too...
  //
  Polynomial<NT> xPolynomial(const NT & y) ;
  
  // Returns the polynmoial obtained by fixing X=x
	
  Polynomial<NT> fixX( const NT &x ) ;
	
  // Returns the polynmoial obtained by fixing Y=y
	
  Polynomial<NT> fixY( const NT &y ) ;
	
  // getYdegree()
  int getYdegree() const;
  
  // getXdegree()
  int getXdegree() const;

  // getCoeff(i)
  Polynomial<NT> getCoeff(int i) const;

  // getTrueYdegree
  int getTrueYdegree() ;

//////////////////////////////////////////////////
////  For interval evaluation, we provide THREE VERSIONS of eval.
////  	We need this in the work on Modified Miranda.
////  	-- Chee and Shang Wang (Aug2011)  
////
////  (1) eval (really should be called "eval0")
////  		based on Horner's evaluation of polynomial.
////
////  (2) eval1 -- using the mean value form
////
////  		f(I,J) = f(mx,my) + fx(I,J).I' + fy(I,J).J'
////
////  		where m(I)=mx,  m(J)=my, I' = I-mx,  J' = J-my
////
////		NOTE: fx(I,J) and fy(I,J) uses eval0.  We can rewrite eval1 in terms of eval0:
////		
////		eval1(f, I, J) = eval0(f, mx, my) + eval0(fx, I, J).I' + eval0(fy, I, J).J'.
////
////  (3) eval2 -- using the second order mean value form
////
////  		f(I,J) = f(mx,my + fx(mx,my).I' + fy(mx,my).J'
////  			 + fxx(I,J).I'^2  + fyy(I,J).J'^2  + 2.fxy(I,J).I'.J'
////
////		THIS VERSION HAS QUADRATIC CONVERGENCE!
////		NOTE: fxx(I,J), fyy(I,J) and fxy(I,J) uses eval0 (Horner's rule).
////		We can rewrite eval2 in terms of eval1 and eval0:
////
////		eval2(f, I, J) = eval0(f, mx, my) + eval1(fx, I, J).I' + eval1(fy, I, J).J'.
////
////  (3) eval3 -- using the slop form (see Stahl's thesis)
////
////  		f(I,J) = f(mx,my) + (I-mx)*eval(g(x)) + (Y-my)*eval0(h(I,J)).
////        This is based upon the slope form of a function: F(y) = F(my) + (Y-my)*G(y),
////        where G(y) is the slope function. 
////		THIS VERSION ALSO HAS QUADRATIC CONVERGENCE!
////		
	

  // templated eval(also known as eval0)
  // assume NT \subseteq T 
  // T requires *, + , =, *=, +=
  // the last two operators 
  // are only for efficiency
  // sufficient for composition
  template< class T > 
  T eval( const T& x, const T& y ) const; 
  //T eval1( const T& x, const T& y ) const; 
  //T eval2( const T& x, const T& y ) const; 

  template < class T >
  IntervalT<T> eval1( const IntervalT<T> &x, const IntervalT<T> &y ) const;

  template < class T >
  IntervalT<T> eval2( const IntervalT<T> &x, const IntervalT<T> &y ) const;
  
  template < class T >
  IntervalT<T> eval3( const IntervalT<T> &x, const IntervalT<T> &y ) const;

	

  // operator version of eval
  template< class T >
    T operator() ( const T &x, const T &y ) const { return eval( x, y ); } 

  //eval(x,y)
  Expr eval(Expr x, Expr y) ;//Evaluate the polynomial at (x,y)

  /// Five special bipolys:
  ///   REMARK: note that we return references to the static variables
  ///   representing these constants!
  /// Zero() = 0
  static const BiPoly<NT> & Zero();
  /// Unity() = 1
  static const BiPoly<NT> & Unity();
  /// NegUnity() = -1
  static const BiPoly<NT> & NegUnity();
  /// IdentityX = x
  static const BiPoly<NT> & IdentityX();
  /// IdentityY = y
  static const BiPoly<NT> & IdentityY();

  ////////////////////////////////////////////////////////
  // Polynomial arithmetic (these are all self-modifying)
  ////////////////////////////////////////////////////////
  
  // Expands the nominal y-degree to n;
  //	Returns n if nominal y-degree is changed to n
  //	Else returns -2

  int expand(int n);

  // contract() gets rid of leading zero polynomials
  //	and returns the new (true) y-degree;
  //	It returns -2 if this is a no-op

  int contract();
  
  // return true if char c is instance of a int, float or rational
  bool is_float_or_rational(char c);
  
  // Self-assignment
  BiPoly<NT> & operator=( const BiPoly<NT>& P);

  // Self-addition
  BiPoly<NT> & operator+=( const BiPoly<NT>& P);
   
  // Self-subtraction
  BiPoly<NT> & operator-=( const BiPoly<NT>& P);

  // Self-multiplication
  BiPoly<NT> & operator*=( const BiPoly<NT>& P);
  

  // Multiply by a polynomial in X
  BiPoly<NT> & mulXpoly( Polynomial<NT> & p);

  //Multiply by a constant
  BiPoly<NT> & mulScalar(const NT & c);

  // mulYpower: Multiply by Y^i (COULD be a divide if i<0)
  BiPoly<NT> & mulYpower(int s);
  
  // Divide by a polynomial in X.
  // We replace the coeffX[i] by the pseudoQuotient(coeffX[i], P)
  BiPoly<NT> & divXpoly( Polynomial<NT> & p);
  
  //Using the standard definition of pseudRemainder operation.
  //	--No optimization!
  BiPoly<NT>  pseudoRemainderY (BiPoly<NT> & Q);

  //Partial Differentiation wrt Y
  BiPoly<NT> & differentiateY();

  //Partial Differentiation wrt X
  BiPoly<NT> & differentiateX();

  //Multiple Partial Differentiation wrt X and Y
  BiPoly<NT> & differentiateXY(int m, int n);//m times wrt X and n times wrt Y

  //Discriminant of P wrt Y
  Polynomial<NT> discriminantY( BiPoly<NT> P);

  ///Discriminant of P wrt X
  Polynomial<NT> discriminantX( BiPoly<NT> P);

  //Represents the bivariate polynomial in (R[X])[Y] as a member
  //of (R[Y])[X].
  //But since our polynomials in X can only have NT coeff's thus
  // to represent the above polynomial we switch X and Y once
  // the conversion has been done.
  //NOTE: This is different from replacing X by Y which was
  //      done in the case of the constructor from a polynomial in X
  //Needed to calculate resultant wrt X.
  BiPoly<NT> & convertXpoly();

  Polynomial<NT> getCoeffY( int i ) const;
  Polynomial<NT>& getCoeffY( int i );
  
  //Set Coeffecient to the polynomial passed as a parameter
  bool setCoeff(int i, Polynomial<NT> p);

  void reverse();
  Polynomial<NT> replaceYwithX();

  //Binary-power operator
  BiPoly<NT>& pow(unsigned int n);

  //Returns a Bipoly corresponding to s, which is supposed to
  //contain as place-holders the chars 'x' and 'y'.
  BiPoly<NT> getbipoly(string s);
};//BiPoly Class

  ////////////////////////////////////////////////////////
  // Four special BiPoly's
  ////////////////////////////////////////////////////////
	  // Zero() = 0
	template < class NT >
	const BiPoly<NT> & BiPoly<NT>::Zero() {
	  static BiPoly<NT> zero;
	  return zero;
	}//Zero
	
	  // Unity() = 1
	template < class NT >
	const BiPoly<NT> & BiPoly<NT>::Unity() {
	  static BiPoly<NT> unity("1");
	  return unity;
	}//Unity
	
	  // NegUnity() = -1
	template < class NT >
	const BiPoly<NT> & BiPoly<NT>::NegUnity() {
	  static BiPoly<NT> negUnity("-1");
	  return negUnity;
	}//NegUnity
	
	  // IdentityX() = x
	template < class NT >
	const BiPoly<NT> & BiPoly<NT>::IdentityX() {
	  static BiPoly<NT> identX("x");
	  return identX;
	}//IdentityX
	
	  // IdentityY() = y
	template < class NT >
	const BiPoly<NT> & BiPoly<NT>::IdentityY() {
	  static BiPoly<NT> identY("y");
	  return identY;
	}//IdentityY

  ////////////////////////////////////////////////////////
  // Helper Functions
  ////////////////////////////////////////////////////////
//Experimental version of constructor from strings containing general 
//parentheses


// zeroPinY(P)
//	checks whether a Bi-polynomial is a zero Polynomial
template <class NT>
bool zeroPinY(BiPoly<NT> & P);

// gcd(P,Q)
//   This gcd is based upon the subresultant PRS to avoid
//   exponential coeffecient growth and gcd computations, both of which 
//   are expensive since the coefficients are polynomials

template <class NT>
BiPoly<NT> gcd( BiPoly<NT>& P ,BiPoly<NT>& Q);

// resY(P,Q):
//      Resultant of Bi-Polys P and Q w.r.t. Y.
//      So the resultant is a polynomial in X
template <class NT>
Polynomial<NT>  resY( BiPoly<NT>& P ,BiPoly<NT>& Q);

// resX(P,Q):
//      Resultant of Bi-Polys P and Q w.r.t. X.
//      So the resultant is a polynomial in Y
//	We first convert P, Q to polynomials in X. Then 
// 	call resY and then turn it back into a polynomial in Y
//	QUESTION: is this last switch really necessary???
template <class NT>
Polynomial<NT>  resX( BiPoly<NT>& P ,BiPoly<NT>& Q);

//Equality operator for BiPoly
template <class NT>
bool operator==(const BiPoly<NT>& P, const BiPoly<NT>& Q);

//Inequality operator for BiPoly
template <class NT>
bool operator!=(const BiPoly<NT>& P, const BiPoly<NT>& Q);

//Addition operator for BiPoly
template <class NT>
 BiPoly<NT> operator+(const BiPoly<NT>& P, const BiPoly<NT>& Q);

//Subtraction operator for BiPoly
template <class NT>
 BiPoly<NT> operator-(const BiPoly<NT>& P, const BiPoly<NT>& Q);

//Multiplication operator for BiPoly
template <class NT>
 BiPoly<NT> operator*(const BiPoly<NT>& P, const BiPoly<NT>& Q);

 // output operator for polynomial
 // implemented on Aug, 2011
 template <class NT>
 std::ostream& operator<<( std::ostream& P, BiPoly<NT>& Q);


  ////////////////////////////////////////////////////////
  //Curve Class
  //  	extends BiPoly Class
  ////////////////////////////////////////////////////////

template < class NT >
class Curve : public BiPoly<NT> {
public:
  // Colors for plotting curves

  const static int NumColors=7;
  static double red_comp(int i){
  	static double RED_COMP[] = {0.9, 0.8, 0.7, 0.6, 0.8, 0.8, 0.7};
	return RED_COMP[i % NumColors];
  }
  static double green_comp(int i){
  	static double GREEN_COMP[] = {0.5, 0.9, 0.3, 0.9, 0.7, 0.55, 0.95};
	return GREEN_COMP[i % NumColors];
  }
  static double blue_comp(int i){
  	static double BLUE_COMP[] = {0.8, 0.3, 0.8, 0.5, 0.4, 0.85, 0.35};
	return BLUE_COMP[i % NumColors];
  }

  Curve(); // zero polynomial
  
  //Curve(vp):
  //    construct from a vector of polynomials
  Curve(std::vector<Polynomial<NT> >& vp);
  //	  : BiPoly<NT>(vp){
  //}
  
  //Curve(p):
  //	Converts a polynomial p(X) to a BiPoly in one of two ways:
  // 	    (1) if flag is false, the result is Y-p(X) 
  // 	    (2) if flag is true, the result is p(Y) 
  //    The default is (1) because we usually want to plot the
  //        graph of the polynomial p(X)
  Curve(Polynomial<NT> p, bool flag=false);
  //	  : BiPoly<NT>(p, flag){
  //}

  //Curve(deg, d[], C[]):
  //	Takes in a list of list of coefficients.
  //	Each cofficient list represents a polynomial in X
  //
  //  deg - ydeg of the bipoly
  //  d[] - array containing the degrees of each coefficient (i.e., X poly)
  //  C[] - list of coefficients, we use array d to select the
  //      coefficients
  Curve(int deg, int *d, NT *C);
  //	  : BiPoly<NT>(deg, d, C){
  //}

  Curve(const BiPoly<NT> &P);
  //	  : BiPoly<NT>(P){
  //}

  //Curve(n) -- the nominal y-degree is n
  Curve(int n);

  //Creates a curve from a string (no parentheses, no *, no =)
  Curve(const string & s, char myX='x', char myY='y');
  Curve(const char* s, char myX='x', char myY='y');

  /////////////////////////////////////////////////////////////////////////
  // verticalIntersections(x, vecI, aprec=0, range=[1,0]):
  //    The list vecI is passed an isolating intervals for y's such that (x,y)
  //    lies on the curve, and y lies in the range.
  //    If aprec is non-zero (!), the intervals have with < 2^{-aprec}.
  //    If range=[1,0], we want roots in (-infty,+infty)
  //    Return is -2 if curve equation does not depend on Y
  //    	-1 if infinitely roots at x,
  //    	0 if no roots at x
  //    	1 otherwise

  int verticalIntersections(const BigFloat & x, BFVecInterval & vI,
		    int aprec=0, BFInterval range=INVALID_BFInterval );
  
  // TO DO: 
  // 		horizontalIntersections(...)
  
  /////////////////////////////////////////////////////////////////////////
  // plot(eps, x1, y1, x2, y2)
  //
  // 	All parameters have defaults
  //
  //    Gives the points on the curve at resolution "eps".  Currently,
  //    eps is viewed as delta-x step size (but it could change).
  //    The display is done in the rectangale 
  //    defined by [(x1, y1), (x2, y2)].
  //    The output is written into a file in the format specified
  //    by our drawcurve function (see COREPATH/ext/graphics).
  //
  //    Heuristic: the open polygonal lines end when number of roots
  //    changes...
  //
  int  plot( BigFloat eps=0.1, BigFloat x1=-1.0,
	     BigFloat y1=-1.0, BigFloat x2=1.0, BigFloat y2=1.0, int fileNo=1);


// selfIntersections():
//   this should be another member function that lists
//   all the self-intersections of a curve
//
//  template <class NT>
//  void selfIntersections(BFVecInterval &vI){
//  ...
//  }

};// Curve class


  ////////////////////////////////////////////////////////
  // BiPoly Composition (Davidson)
  ////////////////////////////////////////////////////////
  
// Bivariate = Univariate( Bivariate ) 

// standard composition, for P = c_i X^i
template < class NT >
BiPoly<NT> composeNaive ( const Polynomial<NT>& P, const BiPoly<NT>& X ) { 

  BiPoly<NT> P_X;
  int i, deg_P;
  NT ci;
 
  deg_P = P.getTrueDegree();
  
  for ( i = 0 ; i <= deg_P ; i++ ) { 
    ci = P.getCoeff(i);      
    if ( ci == 0 ) continue;  // ! requires == 0 defined for NT .. ?isZero?

    BiPoly<NT> pKi = X;         // maybe move pKi outside and reuse, if assignment clears properly
    pKi.pow(i);
    P_X += pKi.mulScalar(ci);
  }
  return P_X;
}

// Bivariate = Univariate( Bivariate ) 

// standard composition, for P = c_i X^i
template < class NT >
BiPoly<NT> composeHorner ( const Polynomial<NT>& P, const BiPoly<NT>& X ) { 

  int i, deg_P;
  NT ci;
 
  deg_P = P.getTrueDegree();
  if ( deg_P < 0 ) { return BiPoly<NT>(); } //return zero polynomial

  int cd = 0;
  ci = P.getCoeff(deg_P);

  BiPoly<NT> P_X = BiPoly<NT>( 0, &cd, &ci );

  for ( i = deg_P-1 ; i >= 0 ; --i ) { 
    P_X *= X;
    ci = P.getCoeff(i);      
    P_X += BiPoly<NT>( 0, &cd, &ci );
  }

  return P_X;

}

template < class NT >
BiPoly<NT> compose ( const Polynomial<NT>& P, const BiPoly<NT>& X ) { 
  return P(X);
}


// compose sequence c_i A^i B^{deg-i}, using c_i of P 
// used in Moebius transform

template < class NT >
BiPoly<NT> composeBinary ( const Polynomial<NT>& P, const BiPoly<NT>& A, const BiPoly<NT>& B ) { 

  BiPoly<NT> P_AB;
  int i, deg_P;
  NT ci;
  
  deg_P = P.getTrueDegree();

  for ( i = 0 ; i <= deg_P ; i++ ) { 

    ci = P.getCoeff(i);      
    if ( ci == 0 ) continue;

    BiPoly<NT> pKA = A;
    BiPoly<NT> pKB = B;
    pKA.pow(i);
    pKB.pow(deg_P-i);
    pKA *= pKB;
    P_AB += pKA.mulScalar(ci);

  }

  return P_AB;

}

// This file used to include composePoly.h. However, the functions
// in that file have since been moved into Curves.h / Poly.h. All
// compose* functions were originally written by Philip Davidson (2006).


// Bivariate from Bivariate 
template< class NT >
BiPoly<NT> composeNaive ( const BiPoly<NT>& P, const BiPoly<NT>& A, const BiPoly<NT>& B ) { 

  BiPoly<NT> P_XY;

  //P is a Bivariate Polynomial
  //  cerr << "deg(y)" << P.getYdegree() << endl;

  int deg_P = P.getYdegree();

  Polynomial<NT> ci;
  
  for ( int i = 0 ; i <= deg_P ; i++ ) { 

    ci = P.getCoeffY( i );                     //x-Polynomial coefficient
    //    cerr << "i: " << i << "->" << ci << endl;
    if ( ci.getTrueDegree() == -1 ) continue;  //zero-case

    // xCi = ci(A);
    // bivariate = univariate( bivariate )
    BiPoly<NT> xCi = composeNaive( ci, A );         //as above


    // yCi = B^i
    BiPoly<NT> yCi = B;                        //init with Y-argument
    yCi.pow(i);                                //raise to power

    // yCi = ci(A) B^i
    // multiply the two BiPoly
    yCi  *= xCi; 


    //accumulate in result;
    P_XY += yCi; 

  }  

  return P_XY;

}



// Bivariate from Bivariate 
template< class NT >
BiPoly<NT> composeHorner ( const BiPoly<NT>& P, const BiPoly<NT>& A, const BiPoly<NT>& B ) { 


  //P is a Bivariate Polynomial

  int deg_P = P.getYdegree();

  if ( deg_P < 0 ) return BiPoly<NT>(); //return zero polynomial

  Polynomial<NT> ci;
  
  std::vector< Polynomial<NT> > civ;
  civ.push_back( P.getCoeffY(deg_P) );

  BiPoly<NT> P_XY = compose( civ[0], A );

  for ( int i = deg_P - 1 ; i >= 0 ; --i ) { 

    P_XY *= B;

    civ[0] = P.getCoeffY( i );                     //x-Polynomial coefficient
    if ( civ[0].getTrueDegree() == -1 ) continue;
    P_XY += compose( civ[0], A );

  }  

  return P_XY;

}


// Bivariate from Bivariate 
template< class NT >
BiPoly<NT> compose ( const BiPoly<NT>& P, const BiPoly<NT>& A, const BiPoly<NT>& B ) { 
  return P( A, B );
}


  ////////////////////////////////////////////////////////
  // Curve helper functions
  ////////////////////////////////////////////////////////


//Xintersections(C, D, vI):
//  returns the list vI of x-ccordinates of possible intersection points.
//  Assumes that C & D are quasi-monic.(or generally aligned)
template <class NT>
void  Xintersections( Curve<NT>& P ,Curve<NT>& Q, BFVecInterval &vI);

//Yintersections(C, D, vI):
//	similar to Xintersections
template <class NT>
void  Yintersections( Curve<NT>& P ,Curve<NT>& Q, BFVecInterval &vI);

// Display Intervals
void showIntervals(char* s, BFVecInterval &vI);

// Set Display Parameters
// ...

////////////////////////////////////////////////////////
// IMPLEMENTATIONS ARE FOUND IN Curves.tcc
////////////////////////////////////////////////////////

// This class is used to cover up differences between
// the Sturm and Descartes root isolation methods. Currently,
// Sturm is used as tests showed that it was faster.
//
// NOTE(narayan): This may not be the case any longer, as
// progs/curves shows that Descartes is consistently faster
// than Sturm.
//
// NOTE(narayan): The "correct" design in terms of software
// engineering is to have Sturm<NT> and Descartes<NT> inherit
// from the same interface.
template <typename NT> class RootIsolator {
private:
  //---------------------------------------------------
  //
  // CHANGE THIS typedef IF YOU WISH TO CHANGE THE ROOT
  // ISOLATION METHOD.
  //
  //---------------------------------------------------
  //typedef Descartes<NT> IsolationMethod;
  typedef Sturm<NT> IsolationMethod;
public:
  RootIsolator(const Polynomial<NT> &poly) :
    poly_(poly),
    method_(new IsolationMethod(poly)) {
  }
  void isolateRoots(BFVecInterval &vI) {
    method_->isolateRoots(vI);
  }
  void isolateRoots(BFInterval &i, BFVecInterval &vi) {
    method_->isolateRoots(i.first, i.second, vi);
  }
  void newtonRefineAllRoots(BFVecInterval &vi, int a_prec) {
    method_->newtonRefineAllRoots(vi, a_prec);
  }
  ~RootIsolator() {
    delete method_;
  }
private:
  const Polynomial<NT> &poly_;
  IsolationMethod *method_;
};

// Curves.tcc does not concern itself with the type of root
// isolation method used any longer. All calls to isolation
// methods are directed through the class RootIsolator above.
#include <CORE/poly/Curves.tcc>

CORE_END_NAMESPACE
#endif
/*************************************************************************** */
// END
/*************************************************************************** */
