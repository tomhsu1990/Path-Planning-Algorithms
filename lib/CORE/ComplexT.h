/******************************************************************
 * File: ComplexT.h
 * Synopsis:
 *      Templated Interval class that defines the standard operations
 *      on complex numbers.
 *
 *        ComplexT<NT>
 *        PolarComplexT<NT>
 *
 *      This is used in progs/ceval, implementing the Sagraloff-Yap
 *      complex root isolation algorithm based on 8-point test.
 *
 *
 * Written by
 *       Narayan Kamath (kamath.narayan@gmail.com), June 2010 (Oxford)
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *****************************************************************/

#ifndef COMPLEX_T_H_
#define COMPLEX_T_H_

#include <iostream>

#include "SupportT.h"

// General note :
// This file provides two representations of complex numbers,
// the polar and cartesian form. It also provides a convenient
// mechanism for converting between forms. Conversion between
// forms involves trignometric operations AND WILL SUFFER FROM
// PRECISION ISSUES. (except at Level 3), so use it with care.

// Forward declaration for the Polar complex number class.
template <typename NT> class PolarComplexT;

// This class represents a complex number in its cartesian form
// (its real and imaginary parts are stored separately. This
// class defines the standard complex arithmetic operations
// {+-/*} and also the arg( ) and mod( ) operations. The arg
// operation uses the principal branch [0, 2*pi) .
template <typename NT> class ComplexT {
public:
  // Construct a complex number with its real and imaginary
  // parts.
  ComplexT(const NT &re, const NT &im) :
    re_(re),
    im_(im) {
  }
  // To provide a natural conversion between a real number
  // and a complex number.
  ComplexT(const NT &re) :
    re_(re),
    im_(0) {
  }
  // Default constructor, this is the representation of zero.
  ComplexT() : re_(0), im_(0) { }
  ComplexT(const ComplexT<NT> &rhs) :
    re_(rhs.re_),
    im_(rhs.im_) {
  }
  // Convert between polar and complex form.
  ComplexT(const PolarComplexT<NT> &rhs) :
    re_(rhs.re()),
    im_(rhs.im()) {
  }
  // We dont allocate anything on the heap.
  ~ComplexT() { }

  // Accessors for real and imaginary parts.
  const NT &re() const {
    return re_;
  }
  const NT &im() const {
    return im_;
  }
  // The argument of the complex number.
  //
  // The value of arg( ) will always lie between
  // [0, 2pi)
  const NT arg() const {
    NT arg = complex_support::atanT<NT>(im_ / re_);
    if (arg < 0) {
      arg += complex_support::piT<NT>();
    }

    return arg;
  }
  // Return |z| .
  const NT mod() const {
    return sqrt(re_ * re_ + im_ * im_);
  }

  // Setting a complex number to a real number. We provide
  // an overloaded integer version of the same function as well.
  ComplexT<NT> &operator=(const NT &r) {
    re_ = r;
    im_ = 0;
    return *this;
  }
  ComplexT<NT> &operator=(const int r) {
    re_ = r;
    im_ = 0;
    return *this;
  }

  // Standard operations on complex numbers. We currently define
  // the self modifying operators inline. The four standard
  // operations on intervals are (as with real numbers) multiplication
  // addition, subraction and division. Multiplication with a scalar
  // is treated as a special case (for efficiency, since we can
  // always "upgrade" a real number to a complex number and then
  // perform the multiplication.

  // Multiplication
  // (a + ib)(c+ id) = (ac - bd) + i(ad + bc)
  ComplexT<NT>& operator*=(const ComplexT<NT> &t) {
    NT re = re_ * t.re_ - im_ * t.im_;
    im_ = re_ * t.im_ + im_ * t.re_;
    re_ = re;

    return *this;
  }
  // Addition.
  // (a + ib) + (c + id) = (a + c) + i(b + d)
  ComplexT<NT>& operator+=(const ComplexT<NT> &s) {
    re_ += s.re_;
    im_ += s.im_;

    return *this;
  }
  // Subtraction.
  // (a + ib) - (c + id) = (a - c) + i(b - d)
  ComplexT<NT>& operator-=(const ComplexT<NT> &s) {
    re_ -= s.re_;
    im_ -= s.im_;

    return *this;
  }
  // Division.
  // (a + ib) / (c + id)
  // Multiply both numerator and denominator by (c - id) to get
  //
  // (ac + bd)/(c^2 + d^2) + i (bc - ad)/(c^2 + d^2)
  ComplexT<NT>& operator/=(const ComplexT<NT> &t) {
    NT mag_sqr = t.re_ * t.re_ + t.im_ * t.im_;
    NT re = (re_*t.re_ + im_*t.im_) / mag_sqr;
    im_ = (im_*t.re_ - re_*t.im_) / mag_sqr;
    re_ = re;

    return *this;
  }

  // Tests for equality with other complex numbers, we also
  // define the notion of equality with a real number as well.
  bool operator==(const ComplexT<NT> &rhs) const {
    return (this == &rhs) || (re_ == rhs.re_ && im_ == rhs.im_);
  }
  bool operator!=(const ComplexT<NT> &rhs) const {
    return !this->operator==(rhs);
  }
  bool operator==(const NT &rhs) const {
    return re_ == rhs && im_ == 0;
  }
  bool operator!=(const NT &rhs) const {
    return !this->operator==(rhs);
  }
private:
  NT re_;
  NT im_;
};

// The polar form of a complex number.
template <typename NT> class PolarComplexT {
public:
  // Construct a complex number with its real and imaginary
  // parts.
  PolarComplexT(const NT &mod, const NT &arg) :
    mod_(mod),
    arg_(arg) {

    while (arg_ >= 2*PI) {
      arg_ -= 2*PI;
    }

    // Remove this assert, it is very expensive.
    assert (arg_ >= 0 && arg_ < 2*PI);
  }
  // To provide a natural conversion between a real number
  // and a complex number.
  PolarComplexT(const NT &re) :
    mod_(re),
    arg_(0) {
  }
  PolarComplexT() : mod_(0), arg_(0) { }
  PolarComplexT(const PolarComplexT<NT> &rhs) :
    mod_(rhs.mod()),
    arg_(rhs.arg()) {
  }
  PolarComplexT(const ComplexT<NT> &rhs) :
    mod_(rhs.mod()),
    arg_(rhs.arg()) {

  }
  // We dont allocate anything on the heap.
  ~PolarComplexT() { }

  // Accessors for real and imaginary parts.
  const NT re() const {
    return mod_ * complex_support::cosT<NT>(arg_);
  }
  const NT im() const {
    return mod_ * complex_support::sinT<NT>(arg_);
  }
  // The argument of the complex number.
  //
  // The value of arg( ) will always lie between
  // [0, 2pi)
  const NT &arg() const {
    return arg_;
  }
  const NT &mod() const {
    return mod_;
  }

  // Setting a complex number to a real number. We provide
  // an overloaded integer version of the same function as well.
  PolarComplexT<NT> &operator=(const NT &r) {
    mod_ = r;
    arg_ = 0;
    return *this;
  }
  PolarComplexT<NT> &operator=(const int r) {
    mod_ = r;
    arg_ = 0;
    return *this;
  }

  // Standard operations on complex numbers. We currently define
  // the self modifying operators inline. The four standard
  // operations on intervals are (as with real numbers) multiplication
  // addition, subraction and division. Multiplication with a scalar
  // is treated as a special case (for efficiency, since we can
  // always "upgrade" a real number to a complex number and then
  // perform the multiplication.

  // Multiplication
  PolarComplexT<NT>& operator*=(const PolarComplexT<NT> &t) {
    mod_ *= t.mod_;
    arg_ += t.arg_;

    if (arg_ >= 2*PI) {
      arg_ -= 2*PI;
    }

    assert(arg_ >= 0 && arg_ < 2*PI);
    return *this;
  }
  // Division.
  PolarComplexT<NT>& operator/=(const PolarComplexT<NT> &t) {
    mod_ /= t.mod_;
    arg_ -= t.arg_;

    if (arg_ < 0) {
      arg_ += 2*PI;
    }

    assert(arg_ >= 0 && arg_ < 2*PI);
    return *this;
  }

  // Addition. This is very inefficient, we convert the polar
  // number into the regular form, add them , and then convert
  // them back. The same comment holds for subtraction.
  PolarComplexT<NT>& operator+=(const PolarComplexT<NT> &s) {
    // Convert both operations out of polar form.
    ComplexT<NT> a = *this;
    ComplexT<NT> b= s;
    // Perform the addition.
    a += b;
    // Convert back to polar form
    (*this) = a;

    return *this;
  }
  PolarComplexT<NT>& operator-=(const PolarComplexT<NT> &s) {
    ComplexT<NT> a = *this;
    ComplexT<NT> b= s;
    a -= b;
    (*this) = a;

    return *this;
  }

  // Tests for equality with other complex numbers, we also
  // define the notion of equality with a real number as well.
  bool operator==(const PolarComplexT<NT> &rhs) const {
    return (this == &rhs) ||
        (mod_ == rhs.mod_ && (mod_ == 0 || arg_ == rhs.arg_));
  }
  bool operator!=(const PolarComplexT<NT> &rhs) const {
    return !this->operator==(rhs);
  }
  bool operator==(const NT &rhs) const {
    return rhs == mod && arg_ == 0;
  }
  bool operator!=(const NT &rhs) const {
    return !this->operator==(rhs);
  }

private:
  // This is used frequently, so makes sense to keep it static.
  //
  // NOTE: Could be made public if people need it. At Level 2 this
  // will have the precision defined by SupportT.h.
  static const NT PI;
  NT mod_;
  NT arg_;
};

template <typename NT> const NT PolarComplexT<NT>::PI = complex_support::piT<NT>();

// For now implemented in terms of self modifying operators, we
// can make these faster if required, at the cost of some code duplication.

// Operations between cartesian complex numbers.
template <typename NT>
inline ComplexT<NT> operator-(const ComplexT<NT> &lhs, const ComplexT<NT> &rhs) {
  ComplexT<NT> ret = lhs;
  ret-=rhs;
  return ret;
}
template <typename NT>
inline ComplexT<NT> operator+(const ComplexT<NT> &lhs, const ComplexT<NT> &rhs) {
  ComplexT<NT> ret = lhs;
  ret+=rhs;
  return ret;
}
template <typename NT>
inline ComplexT<NT> operator*(const ComplexT<NT> &lhs, const ComplexT<NT> &rhs) {
  ComplexT<NT> ret = lhs;
  ret*=rhs;

  return ret;
}
template <typename NT>
inline ComplexT<NT> operator/(const ComplexT<NT> &lhs, const ComplexT<NT> &rhs) {
  ComplexT<NT> ret = lhs;
  ret/=rhs;
  return ret;
}

// Operations between polar complex numbers.
template <typename NT>
inline PolarComplexT<NT> operator-(const PolarComplexT<NT> &lhs,
                                   const PolarComplexT<NT> &rhs) {
  PolarComplexT<NT> ret = lhs;
  ret-=rhs;
  return ret;
}
template <typename NT>
inline PolarComplexT<NT> operator+(const PolarComplexT<NT> &lhs,
    const PolarComplexT<NT> &rhs) {
  PolarComplexT<NT> ret = lhs;
  ret+=rhs;
  return ret;
}
template <typename NT>
inline PolarComplexT<NT> operator*(const PolarComplexT<NT> &lhs,
    const PolarComplexT<NT> &rhs) {
  PolarComplexT<NT> ret = lhs;
  ret*=rhs;
  return ret;
}
template <typename NT>
inline PolarComplexT<NT> operator/(const PolarComplexT<NT> &lhs,
    const PolarComplexT<NT> &rhs) {
  PolarComplexT<NT> ret = lhs;
  ret/=rhs;
  return ret;
}

// Some useful operations with real numbers. Note that these operators
// are potentially ambiguous with the BigFloat operators.
template <typename NT>
inline ComplexT<NT> operator*(const ComplexT<NT> &lhs, const NT &rhs) {
  ComplexT<NT> ret = lhs;
  ret*=rhs;

  return ret;
}
template <typename NT>
inline ComplexT<NT> operator/(const ComplexT<NT> &lhs, const NT &rhs) {
  ComplexT<NT> ret(lhs.re()/rhs, lhs.im()/rhs);

  return ret;
}

// ostream operator for ComplexT
template <typename NT>
inline std::ostream& operator<<(std::ostream& os, const ComplexT<NT>& x) {
  os << "[" << x.re() << " + (" << x.im() << ")i]";
  return os;
}

// ostream operator for PolarComplexT
template <typename NT>
inline std::ostream& operator<<(std::ostream& os, const PolarComplexT<NT> & x) {
  os << "[ " << x.mod() << "*e^i(" << x.arg() << ")]";
  return os;
}

#endif  // COMPLEX_T_INL_
