/******************************************************************
 * Core Library Version 2.0, October 2007
 * Copyright (c) 1995-2002 Exact Computation Project
 * 
 * File: LinearAlgebraT.h
 * Synopsis:
 *      Templated Linear Algebra Extension of Core Library introducing
 *              class Vector
 *              class Matrix
 *
 * Written by
 *       Jihun Yu (jihun@cs.nyu.edu) (2007)
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $Id: linearAlgebraT.h,v 1.7 2010/06/08 19:50:31 exact Exp $
 *****************************************************************/

#ifndef CORE_LINEAR_ALGEBRAT_H
#define CORE_LINEAR_ALGEBRAT_H

#include <cstdarg>
#include <CORE/CORE.h>

template <class T>
class VectorT;

template <class T>
class MatrixT;

////////////////////////////////////////////////////////////////////////
//  Class Vector
//     Generic vectors
//     Operations implemented:  addition, subtraction, dot product
////////////////////////////////////////////////////////////////////////

template<class T>
class VectorT {
private:
   // NOTE : should ideally use unsigned int for such
   // values.
   int     dim;
   T* _rep;
public:
   class RangeException { };
   class ArithmeticException { };

   explicit VectorT(int);
   VectorT();
   VectorT(T, T);
   VectorT(T, T, T);
   VectorT(const VectorT<T>&);
   VectorT(int, T *);
   ~VectorT();

   const VectorT<T>& operator=(const VectorT<T>&);

   bool operator==(const VectorT<T>&);
   bool operator!=(const VectorT<T>&);
   const VectorT<T>& operator+=(const VectorT<T>&);
   const VectorT<T>& operator-=(const VectorT<T>&);
   const VectorT<T>& operator*=(T);

   const T& operator[](int) const;
   T& operator[](int);

   T norm() const;
   T maxnorm() const;
   T infnorm() const;
   int dimension() const {return dim;}
   bool isZero() const;
   VectorT<T> cross(const VectorT<T> &v) const; 
   static VectorT<T> crossProduct(int, ...);

   template <class U>
   friend VectorT<U> operator+(const VectorT<U>&, const VectorT<U>&);
   template <class U>
   friend VectorT<U> operator-(const VectorT<U>&, const VectorT<U>&);
   template <class U>
   friend VectorT<U> operator-(const VectorT<U>&);
   template <class U>
   friend VectorT<U> operator*(const VectorT<U>&, T);
   template <class U>
   friend VectorT<U> operator*(T, const VectorT<U>&);
   template <class U>
   friend VectorT<U> operator*(const MatrixT<U>&, const VectorT<U>&);
   template <class U>
   friend VectorT<U> operator*(const VectorT<U>&, const MatrixT<U>&);
   template <class U>
   friend T dotProduct(const VectorT<U>&, const VectorT<U>&);

   template<class U>
   friend std::istream& operator>>(std::istream&, VectorT<U>&);
   template<class U>
   friend std::ostream& operator<<(std::ostream&, const VectorT<U>&);
};

////////////////////////////////////////////////////////////////////////
//  Class Matrix
//     Generic matrices
//     Operations implemented:  addition, subtraction, multiplication
////////////////////////////////////////////////////////////////////////

template<class T>
class MatrixT {
private:
   int dim1, dim2;
   T* _rep;

public:
   class RangeException { };
   class ArithmeticException { };

   explicit MatrixT(int);
   MatrixT(int, int);
   MatrixT(int, int, T *);
   MatrixT(T, T,
          T, T);
   MatrixT(T, T, T,
          T, T, T,
          T, T, T);
   MatrixT(const MatrixT<T>&);
   ~MatrixT();

   MatrixT<T>& operator=(const MatrixT<T>&);

   bool operator==(const MatrixT<T>&);
   bool operator!=(const MatrixT<T>&);
   bool isZero() const;
   bool isIdentity() const;

   const MatrixT<T>& operator+=(const MatrixT<T>&);
   const MatrixT<T>& operator-=(const MatrixT<T>&);
   const MatrixT<T>& operator*=(T);

   const T& operator()(int, int) const;
   T& operator()(int, int);

   // added by chen li
   //   const VectorT& row(int i) const;
   //   const VectorT& col(int i) const;
   MatrixT<T> matrixAlgebraRemainder(int, int) const;
   T valueAlgebraRemainder(int, int) const;
   void rowExchange(int, int);
   const MatrixT<T>& transpose();

   T determinant() const;

   // Calculates the (fraction free) inverse and determinant
   // of the matrix T based on the extended bareiss method.
   // The adjoint is written to "outputMatrix" which is assumed to be
   // an NxN matrix, and the return value is the determinant.
   T bareissInverse(MatrixT<T> *outputMatrix) const;

   int dimension_1() const { return dim1; }
   int dimension_2() const { return dim2; }

   template<class U>
   friend MatrixT<U> operator+(const MatrixT<U>&, const MatrixT<U>&);
   template<class U>
   friend MatrixT<U> operator-(const MatrixT<U>&, const MatrixT<U>&);
   template<class U>
   friend MatrixT<U> operator*(const MatrixT<U>&, U);
   template<class U>
   friend MatrixT<U> operator*(U, const MatrixT<U>&);
   template<class U>
   friend VectorT<U> operator*(const VectorT<U>&, const MatrixT<U>&);
   template<class U>
   friend VectorT<U> operator*(const MatrixT<U>&, const VectorT<U>&);
   template<class U>
   friend MatrixT<U> operator*(const MatrixT<U>&, const MatrixT<U>&);
   template<class U>
   friend MatrixT<U> transpose(const MatrixT<U>&);

   template<class U>
   friend U det(const U a, const U b,
                const U c, const U d);
   template<class U>
   friend T det(const VectorT<U>&  u, const VectorT<U> & v);  // u,v are 2d vectors
   
   template<class U>
   friend std::istream& operator>>(std::istream&, MatrixT<U>&);
   template<class U>
   friend std::ostream& operator<<(std::ostream&, const MatrixT<U>&);

}; //Matrix

// Other functions available for Matrices:
   template<class U>
   MatrixT<U> identity(int n);

//////////////////////////////////////////////////////////////////////
// Vector class implementation
//////////////////////////////////////////////////////////////////////
template <class T>
VectorT<T>::VectorT(int d) : dim(d) {
   _rep = new T[dim];
   for (int i = 0; i < dim; i++)
      _rep[i] = 0;
}

template <class T>
VectorT<T>::VectorT() : dim(-1) {
   _rep = NULL;
}

template <class T>
VectorT<T>::VectorT(T x, T y) : dim(2) {
   _rep = new T[dim];
   _rep[0] = x;
   _rep[1] = y;
}

template <class T>
VectorT<T>::VectorT(T x, T y, T z) : dim(3) {
   _rep = new T[dim];
   _rep[0] = x;
   _rep[1] = y;
   _rep[2] = z;
}

template <class T>
VectorT<T>::VectorT(int d, T *element) : dim(d) {
  _rep = new T[dim];
  for (int i = 0; i<dim; i++)
    _rep[i] = element[i];
}

template <class T>
VectorT<T>::VectorT(const VectorT<T>& v) : dim(v.dim) {
   _rep = new T[dim];
   for (int i = 0; i < dim; i++)
      _rep[i] = v._rep[i];
}

template <class T>
VectorT<T>::~VectorT() {
   delete[] _rep;
}

template <class T>
const VectorT<T>& VectorT<T>::operator=(const VectorT<T>& v) {
   if (dim != v.dim) {
      dim = v.dim;
      delete[] _rep;
      _rep = new T[dim];
   }
   for (int i = 0; i < dim; i++)
      _rep[i] = v._rep[i];
   return *this;
}

template <class T>
bool VectorT<T>::operator==(const VectorT<T>& v) {
   if (dim != v.dim) return false;
   for (int i = 0; i < dim; i++)
      if (_rep[i] != v._rep[i]) return false;
   return true;
}

template <class T>
bool VectorT<T>::operator!=(const VectorT<T>& v) {
   return !(*this == v);
}

template <class T>
const VectorT<T>& VectorT<T>::operator+=(const VectorT<T>& v) {
   if (dim != v.dim) throw ArithmeticException();
   for (int i = 0; i < dim; i++)
      _rep[i] += v._rep[i];
   return *this;
}

template <class T>
const VectorT<T>& VectorT<T>::operator-=(const VectorT<T>& v) {
   if (dim != v.dim) throw ArithmeticException();
   for (int i = 0; i < dim; i++)
      _rep[i] -= v._rep[i];
   return *this;
}

template <class T>
const VectorT<T>& VectorT<T>::operator*=(T a) {
   for (int i = 0; i < dim; i++)
      _rep[i] *= a;
   return *this;
}

template <class T>
const T& VectorT<T>::operator[](int i) const {
   if (i < 0 || i >= dim) throw RangeException();
   return _rep[i];
}

template <class T>
T& VectorT<T>::operator[](int i) {
   if (i < 0 || i >= dim) throw RangeException();
   return _rep[i];
}

template <class T>
T VectorT<T>::norm() const {
   return sqrt(dotProduct(*this, *this));
}

template <class T>
T VectorT<T>::maxnorm() const {
   T n = 0;
   //   for (int i = 0; i < dim; i++)
   //      if (fabs(_rep[i]) > n) n = fabs(_rep[i]);
   return n;
}

template <class T>
T VectorT<T>::infnorm() const {
   T n = 0;
   for (int i = 0; i < dim; i++)
      n += fabs(_rep[i]);
   return n;
}

template <class T>
bool VectorT<T>::isZero() const {
  for (int i=0; i<dim; i++)
    if (_rep[i] != 0) return false;
  return true;
}

template <class T>
VectorT<T> operator+(const VectorT<T>& a, const VectorT<T>& b) {
   VectorT<T> v = a;
   v += b;
   return v;
}

template <class T>
VectorT<T> operator-(const VectorT<T>& a, const VectorT<T>& b) {
   VectorT<T> v = a;
   v -= b;
   return v;
}

template <class T>
VectorT<T> operator-(const VectorT<T>& a) {
   return VectorT<T>(0,0,0) - a;
}

template <class T>
VectorT<T> operator*(const VectorT<T>& a, T b) {
   VectorT<T> v = a;
   v *= b;
   return v;
}

template <class T>
VectorT<T> operator*(T a, const VectorT<T>& b) {
   VectorT<T> v = b;
   v *= a;
   return v;
}

template <class T>
T dotProduct(const VectorT<T>& a, const VectorT<T>& b) {
   if (a.dim != b.dim) throw VectorT<T>::ArithmeticException();
   T d = 0;
   for (int i = 0; i < a.dim; i++)
      d += a._rep[i] * b._rep[i];
   return d;
}

template <class T>
VectorT<T> VectorT<T>::cross(const VectorT<T> &v) const {
   if (dim != 3) std::cerr << "VectorT<T>::cross( Vector &v ) -- has to be in 3-dimension\n";
   return VectorT(_rep[1] * v._rep[2] - _rep[2] * v._rep[1], 
                 _rep[2] * v._rep[0] - _rep[0] * v._rep[2], 
                 _rep[0] * v._rep[1] - _rep[1] * v._rep[0]);
}


template <class T>
VectorT<T> VectorT<T>::crossProduct(int dim, ...) {
  // Caution: conform to the calling syntax exactly! -- Chen Li
  // there should be (dim-1) Vector pointers passed in as arguments
  va_list ap;
  // Vector *v[dim-1];
  VectorT<T> **v = new VectorT<T>* [dim-1];
  int arg_no = 0;
  va_start(ap, dim);
  // read in the arguments
  while (arg_no < dim-1) {
    v[arg_no++] = va_arg(ap, VectorT<T>*);
  }
  va_end(ap);

  MatrixT<T> m(dim, dim);
  for (int k=0; k<dim; k++) 
    m(0, k) = 1.0;

  for (int i=0; i < dim-1; i++) {
    // verify the validity of inputs
    if (v[i]->dimension() != dim) 
      throw VectorT<T>::ArithmeticException();
    for (int j=0; j < dim; j++)
      m(i+1, j) = (*v[i])[j];
  }
  
  VectorT<T> result(dim);
  for (int l=0; l<dim; l++)
    result[l] = m.valueAlgebraRemainder(0, l);
  return result;
}  


template <class T>
std::istream& operator>>(std::istream& i, VectorT<T>& v) {
   int dim;
   i >> dim;
   VectorT<T> w(dim);
   for (int j=0; j<dim; ++j)
     i >> w[j];

   v = w;
   return i;
}

template <class T>
std::ostream& operator<<(std::ostream& o, const VectorT<T>& v) {
   o << "Vector(";
   for (int i = 0; i < v.dim; i++)
      o << (i>0?",":"") << v._rep[i];
   o << ")";
   return o;
}

//////////////////////////////////////////////////////////////////////
// Matrix class implementation
//////////////////////////////////////////////////////////////////////

template <class T>
MatrixT<T>::MatrixT(int d) : dim1(d), dim2(d) {
   // SHOULD ASSERT (d>0) 
   _rep = new T[dim1*dim2];
   for (int i = 0; i < dim1*dim2; i++)
      _rep[i] = 0;
}

template <class T>
MatrixT<T>::MatrixT(int d1, int d2) : dim1(d1), dim2(d2) {
   _rep = new T[dim1*dim2];
   for (int i = 0; i < dim1*dim2; i++)
      _rep[i] = 0;
}

template <class T>
MatrixT<T>::MatrixT(T m00, T m01,
               T m10, T m11) : dim1(2), dim2(2) {
   _rep = new T[dim1*dim2];
   _rep[0] = m00; _rep[1] = m01;
   _rep[2] = m10; _rep[3] = m11;
}

template <class T>
MatrixT<T>::MatrixT(T m00, T m01, T m02,
	       T m10, T m11, T m12,
               T m20, T m21, T m22) : dim1(3), dim2(3) {
   _rep = new T[dim1*dim2];
   _rep[0] = m00; _rep[1] = m01; _rep[2] = m02;
   _rep[3] = m10; _rep[4] = m11; _rep[5] = m12;
   _rep[6] = m20; _rep[7] = m21; _rep[8] = m22;
}

template <class T>
MatrixT<T>::MatrixT(int m, int n, T *data) : dim1(m), dim2(n) {
  _rep = new T[dim1 * dim2];
  for (int i=0; i < dim1 * dim2; i++)
    _rep[i] = data[i];
}

template <class T>
MatrixT<T>::MatrixT(const MatrixT<T>& m) : dim1(m.dim1), dim2(m.dim2) {
   _rep = new T[dim1*dim2];
   for (int i = 0; i < dim1*dim2; i++)
      _rep[i] = m._rep[i];
}

template <class T>
MatrixT<T>::~MatrixT() {
   delete[] _rep;
}

template <class T>
MatrixT<T>& MatrixT<T>::operator=(const MatrixT<T>& m) {
  if (dim1*dim2 != m.dim1*m.dim2) {
    delete[] _rep;
    _rep = new T[m.dim1*m.dim2];
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

template <class T>
bool MatrixT<T>::operator==(const MatrixT<T>& m) {
   if (dim1 != m.dim1 || dim2 != m.dim2) return false;
   for (int i = 0; i < dim1*dim2; i++)
      if (_rep[i] != m._rep[i]) return false;
   return true;
}

template <class T>
bool MatrixT<T>::operator!=(const MatrixT<T>& m) {
   return !(*this == m);
}

// Chee: added this for testing matrix inversion (Jun 2010)
template <class T>
bool MatrixT<T>::isZero() const {
  for (int i=0; i<dim1*dim2; i++)
    if (_rep[i] != 0) return false;
  return true;
}

// Chee: added this for testing matrix inversion (Jun 2010)
template <class T>
bool MatrixT<T>::isIdentity() const {
  for (int i=0; i<dim1; i++)
    for (int j=0; j<dim2; j++)
       if (i !=j) {
	 if (_rep[i] != 0) return false;
       } else {
	 if (_rep[i] != 1) return false;
       }
  return true;
}


template <class T>
const MatrixT<T>& MatrixT<T>::operator+=(const MatrixT<T>& m) {
   if (dim1 != m.dim1 || dim2 != m.dim2) throw ArithmeticException();
   for (int i = 0; i < dim1*dim2; i++)
      _rep[i] += m._rep[i];
   return *this;
}

template <class T>
const MatrixT<T>& MatrixT<T>::operator-=(const MatrixT<T>& m) {
   if (dim1 != m.dim1 || dim2 != m.dim2) throw ArithmeticException();
   for (int i = 0; i < dim1*dim2; i++)
      _rep[i] -= m._rep[i];
   return *this;
}

template <class T>
const MatrixT<T>& MatrixT<T>::operator*=(T d) {
   for (int i = 0; i < dim1*dim2; i++)
      _rep[i] *= d;
   return *this;
}

template <class T>
const T& MatrixT<T>::operator()(int i, int j) const {
   if (i < 0 || i >= dim1 || j < 0 || j >= dim2) throw RangeException();
   return _rep[i*dim2 + j];
}

template <class T>
T& MatrixT<T>::operator()(int i, int j) {
   if (i < 0 || i >= dim1 || j < 0 || j >= dim2) throw RangeException();
   return _rep[i*dim2 + j];
}

template <class T>
const MatrixT<T>& MatrixT<T>::transpose() {
   if (dim1 == dim2) {
      for (int i = 0; i < dim1; i++)
         for (int j = i+1; j < dim2; j++) {
            T t = _rep[i*dim2 + j];
            _rep[i*dim2 + j] = _rep[j*dim1 + i];
            _rep[j*dim1 + i] = t;
         }
   } else {
      int d = dim1;
      dim1 = dim2;
      dim2 = d;
      T *r = _rep;
      _rep = new T[dim1*dim2];
      for (int i = 0; i < dim1; i++)
         for (int j = 0; j < dim2; j++)
            _rep[i*dim2 + j] = r[j*dim1 + i];
   }
   return *this;
}

// Chee: added identity matrix constructor
template <class T>
MatrixT<T> identity(int n) {
   MatrixT<T> id(n); 
   for (int i=0; i< n; i++)
     id(i,i) = 1;
   return id;
}

namespace bareiss_algorithm {

template <typename T> bool bareiss(MatrixT<T> *matrix, const unsigned int n_steps) {
  const unsigned int n = matrix->dimension_1();

  // Initialize the pivot transformation array.
  // NOTE : Is there an unwritten rule against using STL in CORE ?
  // the vector< > class is better suited, that way i dont have to remember
  // to delete this array at all exit points. (Variable length arrays are
  // also liberally used through CORE, but they are not ISO C++) .
  unsigned int *p = new unsigned int[n];
  // Initially, p[i] = i.
  for (unsigned int i = 0; i < n; ++i) {
    p[i] = i;
  }

  MatrixT<T> &m = *matrix;

  // The
  const T m_0_0 = 1;

  for (unsigned int k = 0; k < n_steps; ++k) {
    // Calculate the pivot for this iteration.
    for (unsigned int c = k; c < n_steps + 1; ++c) {
      // This is the only case the algorithm fails, when no pivot
      // could be found (all entries of the given rows are zeroes)
      // and we cannot avoid a zero in the diagonal by any column
      // stop.
      if (c == n_steps) {
        delete[] p;
        return false;
      }

      if (m(k, p[c]) != 0) {
        const unsigned int temp = p[c];
        p[c] = p[k];
        p[k] = temp;
        break;
      }
    }

    for (unsigned int i = k + 1; i < n; ++i) {
      for (unsigned int j = k + 1; j < n; ++j) {
        const T &d = k ? m(k-1, p[k-1]) : m_0_0;
        m(i, p[j]) = (m(i, p[j]) * m(k, p[k]) - m(i, p[k]) * m(k, p[j])) / d;
      }
    }
  }

  delete [] p;
  return true;
}

}

template <typename T>
T MatrixT<T>::bareissInverse(MatrixT<T> *output) const {
  // Assert the matrix being inverted is a square
  // matrix, also assert that the output is a square matrix
  // too and moreover of the same size as this.
  assert (dim1 == dim2);
  assert(output->dim1 == output->dim2);
  assert(output->dim1 == dim1);

  // To maintain consistency with the algorithm.
  const unsigned int n = dim1;

  MatrixT<T> ext(2*n);
  const MatrixT<T> &this_m = *this;
  // This specialized loop does the following.
  //
  // Copy the matrix A to the NW corner of ext.
  // Copy I into the NE corner of ext.
  // Copy -I into the SW corner of ext.
  // Copy 0 into the SE corner of ext.
  for (unsigned int i = 0; i < n; ++i) {
    for (unsigned int j = 0; j < n; ++j) {
      // Copy the original matrix to the
      ext(i, j) = this_m(i, j);
      ext(i + n, j + n) = 0;
      if (i == j) {
        ext(i + n,j) = -1;
        ext(i, j + n) = 1;
      } else {
        ext(i + n, j) = 0;
        ext(i, j + n) = 0;
      }
    }
  }

  // Run n bareiss steps.
  // NOTE : Change this to (n-1) if you wish to run
  // n-1 steps.
  if (!bareiss_algorithm::bareiss(&ext, n)) {
    return 0;
  }

  // The adjoint is D'.
  for (unsigned int i = 0; i < n; ++i) {
    for (unsigned int j = 0; j < n; ++j) {
      (*output)(i, j) = ext(i + n, j + n);
    }
  }

  // The determinant is the lower right most entry of
  // A'
  return ext(n-1, n-1);
}

template <class T>
T MatrixT<T>::determinant() const {
  if (dim1 != dim2) throw ArithmeticException();

  const unsigned int n = dim1;

  // We need to make a copy since bareiss is in place.
  MatrixT<T> A = *this;
  if (!bareiss_algorithm::bareiss(&A, n-1)) {
    return 0;
  }

  return A(n-1, n-1);
}

// friend det function (for 2x2 determinants)
template <class T>
T det(const T a, const T b,
                const T c, const T d) {
        return (a*d - b*c);
}

// friend det function (u, v are 2d vectors)
template <class T>
T det(const VectorT<T>& u, const VectorT<T> & v) {
        assert(u.dimension()==2); assert(v.dimension()==2);
        return( u[0] * v[1] - u[1] * v[0] );
}

template <class T>
MatrixT<T> MatrixT<T>::matrixAlgebraRemainder(int m, int n) const {
  MatrixT<T> R(dim1-1, dim2-1);
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

template <class T>
T MatrixT<T>::valueAlgebraRemainder(int m, int n) const {
  if ((m+n) % 2) { // (m+n) is odd
    return -matrixAlgebraRemainder(m, n).determinant();
  } else { // (m+n) is even
    return matrixAlgebraRemainder(m, n).determinant();
  }
}

template <class T>
MatrixT<T> operator+(const MatrixT<T>& a, const MatrixT<T>& b) {
   MatrixT<T> m = a;
   m += b;
   return m;
}

template <class T>
MatrixT<T> operator-(const MatrixT<T>& a, const MatrixT<T>& b) {
   MatrixT<T> m = a;
   m -= b;
   return m;
}

template <class T>
MatrixT<T> operator*(const MatrixT<T>& a, T b) {
   MatrixT<T> m = a;
   m *= b;
   return m;
}

template <class T>
MatrixT<T> operator*(T a, const MatrixT<T>& b) {
   MatrixT<T> m = b;
   m *= a;
   return m;
}

template <class T>
VectorT<T> operator*(const VectorT<T>& a, const MatrixT<T>& b) {
   if (a.dim != b.dim1) throw typename MatrixT<T>::ArithmeticException();
   VectorT<T> v(b.dim2);
   for (int i = 0; i < b.dim2; i++) {
      T d = 0;
      for (int j = 0; j < a.dim; j++)
         d += a._rep[j]*b._rep[j*b.dim2 + i];
      v._rep[i] = d;
   }
   return v;
}

template <class T>
VectorT<T> operator*(const MatrixT<T>& a, const VectorT<T>& b) {
   if (b.dim != a.dim2) throw typename MatrixT<T>::ArithmeticException();
   VectorT<T> v(a.dim1);
   for (int i = 0; i < a.dim1; i++) {
      T d = 0;
      for (int j = 0; j < b.dim; j++)
         d += a._rep[i*a.dim2 + j]*b._rep[j];
      v._rep[i] = d;
   }
   return v;
}

template <class T>
MatrixT<T> operator*(const MatrixT<T>& a, const MatrixT<T>& b) {
   if (a.dim2 != b.dim1) throw typename MatrixT<T>::ArithmeticException();
   MatrixT<T> m(a.dim1, b.dim2);
   for (int i = 0; i < a.dim1; i++)
      for (int j = 0; j < b.dim2; j++) {
         T d = 0;
         for (int k = 0; k < a.dim2; k++)
            d += a._rep[i*a.dim2 + k]*b._rep[k*b.dim2 + j];
         m._rep[i*b.dim2 + j] = d;
      }
   return m;
}

template <class T>
MatrixT<T> transpose(const MatrixT<T>& a) {
   MatrixT<T> m = a;
   m.transpose();
   return m;
}

template <class T>
void MatrixT<T>::rowExchange(int r1, int r2) {
  T tmp;

  MatrixT<T>& A = *this;

  for(int i = 0; i < dim2; i++) {
    tmp = A(r1, i);
    A(r1, i) = A(r2, i);
    A(r2, i) = tmp;
  }
}

template <class T>
std::istream& operator>>(std::istream& i, MatrixT<T>& m) {
   int dim;
   i >> dim;
   MatrixT<T> n(dim);
   for (int j=0; j<dim; ++j)
     for (int k=0; k<dim; ++k)
       i >> n(j,k);
   m = n;
   return i;
}

template <class T>
std::ostream& operator<<(std::ostream& o, const MatrixT<T>& m) {
   o << "MatrixT(";
   for (int i = 0; i < m.dim1; i++) {
      o << (i>0?";":"") << std::endl;
      for (int j = 0; j < m.dim2; j++)
         o << (j>0?",":"") << m._rep[i*m.dim2 + j];
   }
   o << ")";
   return o;
}

#endif

