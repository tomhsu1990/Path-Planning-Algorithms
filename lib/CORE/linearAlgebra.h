/******************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2002 Exact Computation Project
 * 
 * File: LinearAlgebra.h
 * Synopsis:
 *      Linear Algebra Extension of Core Library introducing
 *              class Vector
 *              class Matrix
 *
 * Written by
 *       Shubin Zhao (shubinz@cs.nyu.edu) (2001)
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $Id: linearAlgebra.h,v 1.6 2009/02/24 21:30:13 exact Exp $
 *****************************************************************/

#ifndef CORE_LINEAR_ALGEBRA_H
#define CORE_LINEAR_ALGEBRA_H

#ifndef CORE_LEVEL
#  define CORE_LEVEL 3
#endif

#include <cstdarg>
#include <CORE/CORE.h>

class Vector;
class Matrix;

////////////////////////////////////////////////////////////////////////
//  Class Vector
//     Generic vectors
//     Operations implemented:  addition, subtraction, dot product
////////////////////////////////////////////////////////////////////////

class Vector {
private:
   int     dim;
   double* _rep;
public:
   class RangeException { };
   class ArithmeticException { };

   explicit Vector(int);
   Vector();
   Vector(double, double);
   Vector(double, double, double);
   Vector(double, double, double, double);
   Vector(const Vector&);
   Vector(int, double *);
   ~Vector();

   Vector& operator=(const Vector&);

   bool operator==(const Vector&);
   bool operator!=(const Vector&);
   const Vector& operator+=(const Vector&);
   const Vector& operator-=(const Vector&);
   const Vector& operator*=(double);

   const double& operator[](int) const;
   double& operator[](int);

   double X() const { return _rep[0]; }
   double Y() const { return _rep[1]; }
   double Z() const { return _rep[2]; }
   void setX( const double a ) { _rep[0] = a; }
   void setY( const double a ) { _rep[1] = a; }
   void setZ( const double a ) { _rep[2] = a; }
   void set( const double a, const double b, const double c ){ _rep[0] = a; _rep[1] = b; _rep[2] = c; }

   double norm() const;
   void normalize() { _rep[0] = _rep[0]/norm(); _rep[1] = _rep[1]/norm(); _rep[2] = _rep[2]/norm(); };
   double maxnorm() const;
   double infnorm() const;
   double dimension() const {return dim;}
   bool isZero() const;
   Vector cross(const Vector &v) const; 
   static Vector crossProduct(int, ...);

   friend Vector operator+(const Vector&, const Vector&);
   friend Vector operator+(const Vector&);
   friend Vector operator-(const Vector&, const Vector&);
   friend Vector operator-(const Vector&);
   double operator*(const Vector& other) {
     return _rep[0] * other._rep[0] + _rep[1] * other._rep[1] + _rep[2] * other._rep[2];
   }
   friend Vector operator*(const Vector&, double);
   friend Vector operator*(double, const Vector&);
   friend Vector operator*(const Matrix&, const Vector&);
   friend Vector operator*(const Vector&, const Matrix&);
   friend double dotProduct(const Vector&, const Vector&);

   friend double det(const Vector& u, const Vector& v);  // u,v are 2d vectors
   friend double det(const double a, const double b,
                           const double c, const double d);

   friend std::istream& operator>>(std::istream&, Vector&);
   friend std::ostream& operator<<(std::ostream&, const Vector&);
}; //vector

////////////////////////////////////////////////////////////////////////
//  Class Matrix
//     Generic matrices
//     Operations implemented:  addition, subtraction, multiplication
////////////////////////////////////////////////////////////////////////

class Matrix {
private:
   int dim1, dim2;
   double* _rep;

public:
   class RangeException { };
   class ArithmeticException { };

   explicit Matrix(int);
   Matrix(int, int);
   Matrix(int, int, double *);
   Matrix(double, double,
          double, double);
   Matrix(double, double, double,
          double, double, double,
          double, double, double);
   Matrix(const Matrix&);
   ~Matrix();

   Matrix& operator=(const Matrix&);

   bool operator==(const Matrix&);
   bool operator!=(const Matrix&);

   const Matrix& operator+=(const Matrix&);
   const Matrix& operator-=(const Matrix&);
   const Matrix& operator*=(double);

   const double& operator()(int, int) const;
   double& operator()(int, int);

   // added by chen li
   //   const Vector& row(int i) const;
   //   const Vector& col(int i) const;
   Matrix matrixAlgebraRemainder(int, int) const;
   double valueAlgebraRemainder(int, int) const;
   const Matrix& transpose();

   double determinant() const;

   int dimension_1() const { return dim1; }
   int dimension_2() const { return dim2; }

   friend Matrix operator+(const Matrix&, const Matrix&);
   friend Matrix operator-(const Matrix&, const Matrix&);
   friend Matrix operator*(const Matrix&, double);
   friend Matrix operator*(double, const Matrix&);
   friend Vector operator*(const Vector&, const Matrix&);
   friend Vector operator*(const Matrix&, const Vector&);
   friend Matrix operator*(const Matrix&, const Matrix&);
   friend Matrix transpose(const Matrix&);

   friend std::istream& operator>>(std::istream&, Matrix&);
   friend std::ostream& operator<<(std::ostream&, const Matrix&);

}; //Matrix

Matrix operator+(const Matrix&, const Matrix&);
Matrix operator-(const Matrix&, const Matrix&);
Matrix operator*(const Matrix&, double);
Matrix operator*(double, const Matrix&);
Vector operator*(const Vector&, const Matrix&);
Vector operator*(const Matrix&, const Vector&);
Matrix operator*(const Matrix&, const Matrix&);
Matrix transpose(const Matrix&);

/* 
double det(const double a, const double b,
                const double c, const double d);
// 12/11/08, chee: remove "& v" to remove overloading error with the
// 	friend declaration of det() in Vector class.
double det(const Vector u, const Vector v);  // u,v are 2d vectors
*/

std::istream& operator>>(std::istream&, Matrix&);
std::ostream& operator<<(std::ostream&, const Matrix&);

#endif

