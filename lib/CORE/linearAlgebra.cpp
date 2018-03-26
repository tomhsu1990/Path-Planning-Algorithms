//
//  File: linearAlgebra.ccp
//    	  -- class Vector implementation
//        -- class Matrix implementation
//
//  REMARK: in Core2, we implemented a templated linearAlgebra class
//      in $(COREPATH)/inc/CORE/linearAlgebraT.h.
//
//  TODO:
//  	This should introduce namespace.
//
//  Linear Algebra Extension of the CORE library, ver. 1.0
//
//    Copyright (c) 1998, 1999, 2000 Exact Computation Project
//    written by Igor Pechtchanski (pechtcha@cs.nyu.edu)
//
//
//  $Id: linearAlgebra.cpp,v 1.5 2010/06/16 15:27:54 exact Exp $
//

#ifndef CORE_LEVEL
#define CORE_LEVEL 3
#endif

#include <CORE/linearAlgebra.h>

// Note: 5/23/2010: Narayan pointed out the ZeroConst in the original
// Core, being a static, causes runtime bus error, because
// there is no guarantee that statics across compilation
// units is compiled in the right order.
// We now remove that construct.
static Vector VECTOR_ZERO_2D(0.0, 0.0);
static Vector VECTOR_ZERO_3D(0.0, 0.0, 0.0);

//////////////////////////////////////////////////////////////////////
// Vector class implementation
//////////////////////////////////////////////////////////////////////

Vector::Vector(int d) : dim(d) {
   _rep = new double[dim];
   for (int i = 0; i < dim; i++)
      _rep[i] = 0;
}

Vector::Vector() : dim(-1) {
   _rep = NULL;
}

Vector::Vector(double x, double y) : dim(2) {
   _rep = new double[dim];
   _rep[0] = x;
   _rep[1] = y;
}

Vector::Vector(double x, double y, double z) : dim(3) {
   _rep = new double[dim];
   _rep[0] = x;
   _rep[1] = y;
   _rep[2] = z;
}

// four dimensional points:
Vector::Vector(double x, double y, double z, double w) : dim(4) {
   _rep = new double[dim];
   _rep[0] = x;
   _rep[1] = y;
   _rep[2] = z;
   _rep[3] = w;
}

Vector::Vector(int d, double *element) : dim(d) {
  _rep = new double[dim];
  for (int i = 0; i<dim; i++)
    _rep[i] = element[i];
}

Vector::Vector(const Vector& v) : dim(v.dim) {
   _rep = new double[dim];
   for (int i = 0; i < dim; i++)
      _rep[i] = v._rep[i];
}

Vector::~Vector() {
   delete[] _rep;
}

Vector& Vector::operator=(const Vector& v) {
   if (dim != v.dim) {
      dim = v.dim;
      delete[] _rep;
      _rep = new double[dim];
   }
   for (int i = 0; i < dim; i++)
      _rep[i] = v._rep[i];
   return *this;
}

bool Vector::operator==(const Vector& v) {
   if (dim != v.dim) return false;
   for (int i = 0; i < dim; i++)
      if (_rep[i] != v._rep[i]) return false;
   return true;
}

bool Vector::operator!=(const Vector& v) {
   return !(*this == v);
}

const Vector& Vector::operator+=(const Vector& v) {
   if (dim != v.dim) throw ArithmeticException();
   for (int i = 0; i < dim; i++)
      _rep[i] += v._rep[i];
   return *this;
}

const Vector& Vector::operator-=(const Vector& v) {
   if (dim != v.dim) throw ArithmeticException();
   for (int i = 0; i < dim; i++)
      _rep[i] -= v._rep[i];
   return *this;
}

const Vector& Vector::operator*=(double a) {
   for (int i = 0; i < dim; i++)
      _rep[i] *= a;
   return *this;
}

const double& Vector::operator[](int i) const {
   if (i < 0 || i >= dim) throw RangeException();
   return _rep[i];
}

double& Vector::operator[](int i) {
   if (i < 0 || i >= dim) throw RangeException();
   return _rep[i];
}

double Vector::norm() const {
   return sqrt(dotProduct(*this, *this));
}

double Vector::maxnorm() const {
   double n = 0;
   //   for (int i = 0; i < dim; i++)
   //      if (fabs(_rep[i]) > n) n = fabs(_rep[i]);
   return n;
}

double Vector::infnorm() const {
   double n = 0;
   for (int i = 0; i < dim; i++) {
      n += fabs(_rep[i]);	   
   }
   return n;
}

bool Vector::isZero() const {
  for (int i=0; i<dim; i++)
    if (_rep[i] != 0) return false;
  return true;
}

Vector operator+(const Vector& a, const Vector& b) {
   Vector v = a;
   v += b;
   return v;
}

Vector operator+(const Vector& a) {
   return Vector(0,0,0) + a;
}

Vector operator-(const Vector& a, const Vector& b) {
   Vector v = a;
   v -= b;
   return v;
}

Vector operator-(const Vector& a) {
   return Vector(0,0,0) - a;
}

Vector operator*(const Vector& a, double b) {
   Vector v = a;
   v *= b;
   return v;
}

Vector operator*(double a, const Vector& b) {
   Vector v = b;
   v *= a;
   return v;
}

double dotProduct(const Vector& a, const Vector& b) {
   if (a.dim != b.dim) throw Vector::ArithmeticException();
   double d = 0;
   for (int i = 0; i < a.dim; i++)
      d += a._rep[i] * b._rep[i];
   return d;
}

Vector Vector::cross(const Vector &v) const {
   if (dim != 3){
     std::cerr << "Vector::cross( Vector &v ) -- has to be in 3-dimension\n";
     return Vector(0,0,0);
   }
   return Vector(_rep[1] * v._rep[2] - _rep[2] * v._rep[1], 
                 _rep[2] * v._rep[0] - _rep[0] * v._rep[2], 
                 _rep[0] * v._rep[1] - _rep[1] * v._rep[0]);
}


Vector Vector::crossProduct(int dim, ...) {
  // Caution: conform to the calling syntax exactly! -- Chen Li
  // there should be (dim-1) Vector pointers passed in as arguments
  va_list ap;
  // Vector *v[dim-1];
  Vector **v = new Vector* [dim-1];
  int arg_no = 0;
  va_start(ap, dim);
  // read in the arguments
  while (arg_no < dim-1) {
    v[arg_no++] = va_arg(ap, Vector*);
  }
  va_end(ap);

  Matrix m(dim, dim);
  for (int k=0; k<dim; k++) 
    m(0, k) = 1.0;

  for (int i=0; i < dim-1; i++) {
    // verify the validity of inputs
    if (v[i]->dimension() != dim) 
      throw Vector::ArithmeticException();
    for (int j=0; j < dim; j++)
      m(i+1, j) = (*v[i])[j];
  }
  
  Vector result(dim);
  for (int l=0; l<dim; l++)
    result[l] = m.valueAlgebraRemainder(0, l);
  return result;
}  


std::istream& operator>>(std::istream& i, Vector& v) {
   int dim;
   i >> dim;
   Vector w(dim);
   for (int j=0; j<dim; ++j)
     i >> w[j];

   v = w;
   return i;
}

std::ostream& operator<<(std::ostream& o, const Vector& v) {
   o << "Vector(";
   for (int i = 0; i < v.dim; i++)
      o << (i>0?",":"") << v._rep[i];
   o << ")";
   return o;
}


//////////////////////////////////////////////////////////////////////
// Matrix class implementation
//////////////////////////////////////////////////////////////////////

Matrix::Matrix(int d) : dim1(d), dim2(d) {
   _rep = new double[dim1*dim2];
   for (int i = 0; i < dim1*dim2; i++)
      _rep[i] = 0;
}

Matrix::Matrix(int d1, int d2) : dim1(d1), dim2(d2) {
   _rep = new double[dim1*dim2];
   for (int i = 0; i < dim1*dim2; i++)
      _rep[i] = 0;
}

Matrix::Matrix(double m00, double m01,
               double m10, double m11) : dim1(2), dim2(2) {
   _rep = new double[dim1*dim2];
   _rep[0] = m00; _rep[1] = m01;
   _rep[2] = m10; _rep[3] = m11;
}

Matrix::Matrix(double m00, double m01, double m02,
	       double m10, double m11, double m12,
               double m20, double m21, double m22) : dim1(3), dim2(3) {
   _rep = new double[dim1*dim2];
   _rep[0] = m00; _rep[1] = m01; _rep[2] = m02;
   _rep[3] = m10; _rep[4] = m11; _rep[5] = m12;
   _rep[6] = m20; _rep[7] = m21; _rep[8] = m22;
}

Matrix::Matrix(int m, int n, double *data) : dim1(m), dim2(n) {
  _rep = new double[dim1 * dim2];
  for (int i=0; i < dim1 * dim2; i++)
    _rep[i] = data[i];
}

Matrix::Matrix(const Matrix& m) : dim1(m.dim1), dim2(m.dim2) {
   _rep = new double[dim1*dim2];
   for (int i = 0; i < dim1*dim2; i++)
      _rep[i] = m._rep[i];
}

Matrix::~Matrix() {
   delete[] _rep;
}

Matrix& Matrix::operator=(const Matrix& m) {
  if (dim1*dim2 != m.dim1*m.dim2) {
    delete[] _rep;
    _rep = new double[m.dim1*m.dim2];
  }
  dim1 = m.dim1;
  dim2 = m.dim2;
  int i = 0;
  // Chen Li, if use for statemen it gives some strange compiler error
  // so I changed it to "equivalent" while statement. ???
  //  for (i = 0; i < dim1*dim2; i++) {
  while (i < dim1 * dim2) {
    _rep[i] = m._rep[i];
    i++;
  }
  return (*this);
}

bool Matrix::operator==(const Matrix& m) {
   if (dim1 != m.dim1 || dim2 != m.dim2) return false;
   for (int i = 0; i < dim1*dim2; i++)
      if (_rep[i] != m._rep[i]) return false;
   return true;
}

bool Matrix::operator!=(const Matrix& m) {
   return !(*this == m);
}

const Matrix& Matrix::operator+=(const Matrix& m) {
   if (dim1 != m.dim1 || dim2 != m.dim2) throw ArithmeticException();
   for (int i = 0; i < dim1*dim2; i++)
      _rep[i] += m._rep[i];
   return *this;
}

const Matrix& Matrix::operator-=(const Matrix& m) {
   if (dim1 != m.dim1 || dim2 != m.dim2) throw ArithmeticException();
   for (int i = 0; i < dim1*dim2; i++)
      _rep[i] -= m._rep[i];
   return *this;
}

const Matrix& Matrix::operator*=(double d) {
   for (int i = 0; i < dim1*dim2; i++)
      _rep[i] *= d;
   return *this;
}

const double& Matrix::operator()(int i, int j) const {
// why is this check commented out?
//   if (i < 0 || i >= dim1 || j < 0 || j >= dim2) throw RangeException();
   return _rep[i*dim2 + j];
}

double& Matrix::operator()(int i, int j) {
// why is this check commented out?
//   if (i < 0 || i >= dim1 || j < 0 || j >= dim2) throw RangeException();
   return _rep[i*dim2 + j];
}

const Matrix& Matrix::transpose() {
   if (dim1 == dim2) {
      for (int i = 0; i < dim1; i++)
         for (int j = i+1; j < dim2; j++) {
            double t = _rep[i*dim2 + j];
            _rep[i*dim2 + j] = _rep[j*dim1 + i];
            _rep[j*dim1 + i] = t;
         }
   } else {
      int d = dim1;
      dim1 = dim2;
      dim2 = d;
      double *r = _rep;
      _rep = new double[dim1*dim2];
      for (int i = 0; i < dim1; i++)
         for (int j = 0; j < dim2; j++)
            _rep[i*dim2 + j] = r[j*dim1 + i];
   }
   return *this;
}

// not right!
 // BUG to be fixed!
double Matrix::determinant() const {
  if (dim1 != dim2) throw ArithmeticException();
  Matrix A = *this;
  int i, j, k;
  for (i = 0; i < dim1 - 1; i++) {
    if (A(i,i) == 0) { // pivoting
      int p = i;
      for (int q = i + 1; q < dim1; q++) {
	if (A(q, i) != 0) {
	  p = q;
	  break;
	}
      }
      if (p == i) return 0;
      //      assert(p != i);
      // swap
      double tmp;
      for (int d = i; d < dim1; d++) {
	tmp = A(i, d);
	A(i, d) = A(p, d);
	A(p, d) = tmp;
      }
    }

    for (j = i + 1; j < dim1; j++) {
      A(j, i) /= A(i, i);
      for (k = i + 1; k < dim1; k++) {
	A(j,k) -= A(i,k) * A(j,i);
      }
    }
  }
   
  double det = 1;
  for (i = 0; i < dim1; i++)
      det *= A(i,i);
  return det;
}

// friend det function (for 2x2 determinants)
double det(const double a, const double b,
                const double c, const double d) {
        return (a*d - b*c);
}

// friend det function (u, v are 2d vectors)
double det(const Vector& u, const Vector& v) {
        assert(u.dimension()==2); assert(v.dimension()==2);
        return( u[0] * v[1] - u[1] * v[0] );
}

Matrix Matrix::matrixAlgebraRemainder(int m, int n) const {
  Matrix R(dim1-1, dim2-1);
  int pos=0;

  for (int i = 0; i < dim1; i++) {
    if (i == m) continue;
    for (int j = 0; j < dim2; j++) {
      if (j == n) continue;
      R._rep[pos++] = _rep[i*dim1+j];
    }
  }
  return R;
}

double Matrix::valueAlgebraRemainder(int m, int n) const {
  if ((m+n) % 2) { // (m+n) is odd
    return -matrixAlgebraRemainder(m, n).determinant();
  } else { // (m+n) is even
    return matrixAlgebraRemainder(m, n).determinant();
  }
}

Matrix operator+(const Matrix& a, const Matrix& b) {
   Matrix m = a;
   m += b;
   return m;
}

Matrix operator-(const Matrix& a, const Matrix& b) {
   Matrix m = a;
   m -= b;
   return m;
}

Matrix operator*(const Matrix& a, double b) {
   Matrix m = a;
   m *= b;
   return m;
}

Matrix operator*(double a, const Matrix& b) {
   Matrix m = b;
   m *= a;
   return m;
}

Vector operator*(const Vector& a, const Matrix& b) {
   if (a.dim != b.dim1) throw Matrix::ArithmeticException();
   Vector v(b.dim2);
   for (int i = 0; i < b.dim2; i++) {
      double d = 0;
      for (int j = 0; j < a.dim; j++)
         d += a._rep[j]*b._rep[j*b.dim2 + i];
      v._rep[i] = d;
   }
   return v;
}

Vector operator*(const Matrix& a, const Vector& b) {
   if (b.dim != a.dim2) throw Matrix::ArithmeticException();
   Vector v(a.dim1);
   for (int i = 0; i < a.dim1; i++) {
      double d = 0;
      for (int j = 0; j < b.dim; j++)
         d += a._rep[i*a.dim2 + j]*b._rep[j];
      v._rep[i] = d;
   }
   return v;
}

Matrix operator*(const Matrix& a, const Matrix& b) {
   if (a.dim2 != b.dim1) throw Matrix::ArithmeticException();
   Matrix m(a.dim1, b.dim2);
   for (int i = 0; i < a.dim1; i++)
      for (int j = 0; j < b.dim2; j++) {
         double d = 0;
         for (int k = 0; k < a.dim2; k++)
            d += a._rep[i*a.dim2 + k]*b._rep[k*b.dim2 + j];
         m._rep[i*b.dim2 + j] = d;
      }
   return m;
}

Matrix transpose(const Matrix& a) {
   Matrix m = a;
   m.transpose();
   return m;
}


std::istream& operator>>(std::istream& i, Matrix& m) {
   int dim;
   i >> dim;
   Matrix n(dim);
   for (int j=0; j<dim; ++j)
     for (int k=0; k<dim; ++k)
       i >> n(j,k);
   m = n;
   return i;
}

std::ostream& operator<<(std::ostream& o, const Matrix& m) {
   o << "Matrix(";
   for (int i = 0; i < m.dim1; i++) {
      o << (i>0?";":"");
      for (int j = 0; j < m.dim2; j++)
         o << (j>0?",":"") << m._rep[i*m.dim2 + j];
   }
   o << ")";
   return o;
}



