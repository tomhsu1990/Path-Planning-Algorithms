/****************************************************************************
 * BigRat.h -- Big Rational number class based on mpq in GMP
 *
 * Core Library Version 2.0, March 2006
 * Copyright (c) 1995-2006 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of Core Library (http://cs.nyu.edu/exact/core); you 
 * may redistribute it under the terms of the Q Public License version 1.0.
 * See the file LICENSE.QPL distributed with Core Library.
 *
 * Licensees holding a valid commercial license may use this file in
 * accordance with the commercial license agreement provided with the
 * software.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 * WWW URL: http://cs.nyu.edu/exact/core
 * Email: exact@cs.nyu.edu
 *
 * $Id: BigRat.h,v 1.13 2010/11/23 17:58:36 exact Exp $
 ***************************************************************************/
#ifndef __CORE_BIGRAT_H__
#define __CORE_BIGRAT_H__

#include <CORE/Gmpq.h>
#include <CORE/BigInt.h>

extern std::istream& extract(std::istream &i, mpq_ptr x);

CORE_BEGIN_NAMESPACE

#ifndef CORE_DISABLE_REFCOUNTING
  typedef RcGmpq BigRatBase;
#else
  typedef Gmpq BigRatBase;
#endif

/// \class BigRat BigRat.h
/// \brief BigRat is a big rational number class based on <tt>mpq</tt> in GMP
class BigRat : public BigRatBase {
  typedef BigRatBase base_cls;
public:
  /// \name constructors and destructor
  //@{
  /// default constructor
  BigRat() {}
  /// constructor for <tt>int</tt>
  BigRat(int i) : base_cls(static_cast<long>(i)) {}
  /// constructor for <tt>unsigned int</tt>
  BigRat(unsigned int i) : base_cls(static_cast<unsigned long>(i)) {}
  /// constructor for <tt>long</tt>
  BigRat(long i) : base_cls(i) {}
  /// constructor for <tt>unsigned long</tt>
  BigRat(unsigned long i) : base_cls(i) {}
  /// constructor for <tt>double</tt>
  BigRat(double i) : base_cls(i) {}
  /// constructor for <tt>char*</tt> (no implicit conversion)
  explicit BigRat(const char* s) : base_cls(s) {}
  /// constructor for <tt>std::string</tt> (no implicit conversion)
  explicit BigRat(const std::string& s) : base_cls(s.c_str()) {}
  /// constructor for <tt>BigInt</tt>
  BigRat(const BigInt& z) : base_cls(z.mp()) {}
  /// constructor for <tt>BigInt, BigInt</tt>
  BigRat(const BigInt& num, const BigInt& den) : base_cls(num.mp(), den.mp()) {}
  //@}

public:
  /// \name assignment functions
  //@{
  /// assignment function for <tt>BigRat</tt>
  void set(const BigRat& rhs) { base_cls::operator=(rhs); }
  /// assignment function for <tt>int</tt>
  void set(int i, unsigned long den = 1)
  { set(static_cast<long>(i), den); }
  /// assignment function for <tt>unsigned int</tt>
  void set(unsigned int i, unsigned long den = 1)
  { set(static_cast<unsigned long>(i), den); }
  /// assignment function for <tt>long</tt>
  void set(long i, unsigned long den = 1)
  { mpq_set_si(mp(), i, den); }
  /// assignment function for <tt>unsigned long</tt>
  void set(unsigned long i, unsigned long den = 1)
  { mpq_set_ui(mp(), i, den); }
  /// assignment function for <tt>double</tt>
  void set(double i)
  { mpq_set_d(mp(), i); }
  /// assignment function for <tt>char*</tt>
  int set(const char* str, int base = 0) 
  { return mpq_set_str(mp(), str, base); }
  /// assignment function for <tt>std::string</tt>
  int set(const std::string& str, int base = 0)
  { return mpq_set_str(mp(), str.c_str(), base); }
  /// assignment function for <tt>BigInt</tt>
  void set(const BigInt& x)
  { mpq_set_z(mp(), x.mp()); }
  /// assignment function for <tt>BigInt, BigInt</tt>
  void set(const BigInt& num, const BigInt& den) {
    mpz_set(mpq_numref(mp()), num.mp());
    mpz_set(mpq_denref(mp()), den.mp());
    mpq_canonicalize(mp());
  }
  //@}

  /// \name arithmetic functions
  //@{
  /// addition for <tt>BigRat+BigRat</tt> 
  void add(const BigRat& x, const BigRat& y)
  { mpq_add(mp(), x.mp(), y.mp()); }
  /// subtraction for <tt>BigRat-BigRat</tt>
  void sub(const BigRat& x, const BigRat& y)
  { mpq_sub(mp(), x.mp(), y.mp()); }
  /// multiplication for <tt>BigRat*BigRat</tt>
  void mul(const BigRat& x, const BigRat& y)
  { mpq_mul(mp(), x.mp(), y.mp()); }
  /// division for <tt>BigRat/BigRat</tt> 
  void div(const BigRat& x, const BigRat& y)
  { mpq_div(mp(), x.mp(), y.mp()); }
  /// negation
  void neg(const BigRat& x)
  { mpq_neg(mp(), x.mp()); }
  /// absolute value
  void abs(const BigRat& x)
  { mpq_abs(mp(), x.mp()); }
  /// inversion (i.e., reciprocal)
  void inv(const BigRat& x)
  { mpq_inv(mp(), x.mp()); }
  void inv()
  { BigRat r; r.inv(*this); *this = r; }
  //@}
  //@}
  
  /// \name shift functions
  //@{
  /// left shift
  void mul_2exp(const BigRat& x, int y)
  { mul_2exp(x, static_cast<long>(y)); }
  /// left shift
  void mul_2exp(const BigRat& x, unsigned int y)
  { mul_2exp(x, static_cast<unsigned long>(y)); }
  /// left shift
  void mul_2exp(const BigRat& x, long y)
  { if (y>=0) mpq_mul_2exp(mp(),x.mp(),y); else mpq_div_2exp(mp(),x.mp(),-y); }
  /// left shift
  void mul_2exp(const BigRat& x, unsigned long y)
  { mpq_mul_2exp(mp(), x.mp(), y); }
  /// right shift
  void div_2exp(const BigRat& x, int y)
  { div_2exp(x, static_cast<long>(y)); }
  /// right shift
  void div_2exp(const BigRat& x, unsigned int y)
  { div_2exp(x, static_cast<unsigned long>(y)); }
  /// right shift
  void div_2exp(const BigRat& x, long y)
  { if (y>=0) mpq_div_2exp(mp(),x.mp(),y); else mpq_mul_2exp(mp(),x.mp(),-y); }
  /// right shift
  void div_2exp(const BigRat& x, unsigned long y)
  { mpq_div_2exp(mp(), x.mp(), y); }
  //@}

  /// \name comparison functions
  //@{ 
  /// compare with <tt>BigRat</tt>
  int cmp(const BigRat& x) const
  { return mpq_cmp(mp(), x.mp()); }
  /// compare with <tt>int</tt>
  int cmp(int x, unsigned long den) const
  { return mpq_cmp_si(mp(), x, den); }
  /// compare with <tt>unsigned int</tt>
  int cmp(unsigned int x, unsigned long den) const
  { return mpq_cmp_ui(mp(), x, den); }
  /// compare with <tt>long</tt>
  int cmp(long x, unsigned long den) const
  { return mpq_cmp_si(mp(), x, den); }
  /// compare with <tt>unsigned long</tt>
  int cmp(unsigned long x, unsigned long den) const
  { return mpq_cmp_ui(mp(), x, den); }
  /// compare with <tt>double</tt>
  int cmp(double x) const
  { return this->cmp(BigRat(x)); }
  /// return non-zero if equals (faster than cmp)
  int equal(const BigRat& x) const
  { return mpq_equal(mp(), x.mp()); }
  //@}
  
  /// \name conversion functions
  //@{
  /// return double value
  double get_d() const
  { return mpq_get_d(mp()); }
  /// return signed long value
  long get_si() const
  { return static_cast<long>(get_d()); }
  /// return unsigned long value
  unsigned long get_ui() const
  { return static_cast<unsigned long>(get_d()); }
  std::string get_str(int base = 10) const {
    int len = mpz_sizeinbase(mpq_numref(mp()), base) 
            + mpz_sizeinbase(mpq_denref(mp()), base) + 3;
    _gmp_alloc_cstr tmp(len);
    return std::string(mpq_get_str(tmp.str, base, mp()));
  }
  /// return BigInt value
  BigInt get_z() const
  { BigInt r; mpz_div(r.mp(), mpq_numref(mp()), mpq_denref(mp())); return r; }
  //@}
  
  /// \name miscellaneous functions
  //@{
  /// swap function 
  void swap(BigRat& other)
  { mpq_swap(mp(), other.mp()); }
  /// canonicalization
  void canonicalize() 
  { mpq_canonicalize(mp()); }
  BigRat reciprocal() const { 
    BigRat r; 
    mpz_set(mpq_numref(r.mp()), mpq_denref(mp()));  
    mpz_set(mpq_denref(r.mp()), mpq_numref(mp()));  
    return r;
  }
  /// return numerator
  BigInt numerator() const
  { return BigInt(mpq_numref(mp())); }
  /// return numerator
  BigInt num() const
  { return BigInt(mpq_numref(mp())); }
  /// return denominator
  BigInt denominator() const
  { return BigInt(mpq_denref(mp())); }
  /// return denominator
  BigInt den() const
  { return BigInt(mpq_denref(mp())); }
  /// return sign
  int sgn() const
  { return mpq_sgn(mp()); }
  /// return upper bound of MSB
  long uMSB() const
  { return num().ceillg() - den().floorlg(); }
  /// return lower bound of MSB
  long lMSB() const
  { return num().floorlg() - den().ceillg(); }
  //@}

public: // C++ operators
  /// \name unary, increment, decrement operators
  //@{
  /// unary plus operator
  BigRat operator+() const
  { return BigRat(*this); }
  /// unary negation operator
  BigRat operator-() const
  { BigRat r; r.neg(*this); return r; }
  /// prefix increment operator
  BigRat& operator++()
  { add(*this, 1); return *this; }
  /// postfix increment operator
  BigRat operator++(int)
  { BigRat r(*this); ++(*this); return r; }
  /// prefix decrement operator
  BigRat& operator--()
  { sub(*this, 1); return *this; }
  /// postfix decrement operator
  BigRat operator--(int)
  { BigRat r(*this); --(*this); return r; }
  //@}

  /// \name assignment and compound assignment operators
  //@{
  /// assignment operator for <tt>BigRat</tt>
  BigRat& operator=(const BigRat& rhs)
  { base_cls::operator=(rhs); return *this; }
  /// assignment operator for <tt>int</tt>
  BigRat& operator=(int rhs)
  { set(rhs); return *this; }
  /// assignment operator for <tt>unsigned int</tt>
  BigRat& operator=(unsigned int rhs)
  { set(rhs); return *this; }
  /// assignment operator for <tt>long</tt>
  BigRat& operator=(long rhs)
  { set(rhs); return *this; }
  /// assignment operator for <tt>unsigned long</tt>
  BigRat& operator=(unsigned long rhs)
  { set(rhs); return *this; }
  /// assignment operator for <tt>double</tt>
  BigRat& operator=(double rhs)
  { set(rhs); return *this; }
  /// assignment operator for <tt>char*</tt>
  BigRat& operator=(const char* rhs)
  { set(rhs); return *this; }
  /// assignment operator for <tt>std::string</tt>
  BigRat& operator=(const std::string& rhs)
  { set(rhs); return *this; }
  /// assignment operator for <tt>BigInt</tt>
  BigRat& operator=(const BigInt& rhs)
  { set(rhs); return *this; }

  /// compound assignment operator <tt>+=</tt>
  BigRat& operator+=(const BigRat& rhs)
  { add(*this, rhs); return *this; }
  /// compound assignment operator <tt>-=</tt>
  BigRat& operator-=(const BigRat& rhs)
  { sub(*this, rhs); return *this; }
  /// compound assignment operator <tt>*=</tt>
  BigRat& operator*=(const BigRat& rhs)
  { mul(*this, rhs); return *this; }
  /// compound assignment operator <tt>/=</tt>
  BigRat& operator/=(const BigRat& rhs)
  { div(*this, rhs); return *this; }

  /// compound assignment operator <tt><<=</tt>
  BigRat& operator<<=(int i)
  { mul_2exp(*this, i); return *this; }
  /// compound assignment operator <tt><<=</tt>
  BigRat& operator<<=(unsigned int ui)
  { mul_2exp(*this, ui); return *this; }
  /// compound assignment operator <tt><<=</tt>
  BigRat& operator<<=(long l)
  { mul_2exp(*this, l); return *this; }
  /// compound assignment operator <tt><<=</tt>
  BigRat& operator<<=(unsigned long ul)
  { mul_2exp(*this, ul); return *this; }
  /// compound assignment operator <tt>>>=</tt>
  BigRat& operator>>=(int i) 
  { div_2exp(*this, i); return *this; }
  /// compound assignment operator <tt>>>=</tt>
  BigRat& operator>>=(unsigned int ui) 
  { div_2exp(*this, ui); return *this; }
  /// compound assignment operator <tt>>>=</tt>
  BigRat& operator>>=(long l) 
  { div_2exp(*this, l); return *this; }
  /// compound assignment operator <tt>>>=</tt>
  BigRat& operator>>=(unsigned long ul) 
  { div_2exp(*this, ul); return *this; }
  //@}

#ifndef CORE_DISABLE_OLDNAMES
  /// \name back-compatiable functions
  //@{
  /// Has Exact Division
  static bool hasExactDivision() { return true; }
  /// set value from <tt>const char*</tt>
  int set_str(const char* s, int base = 0) { return set(s, base); }
  /// div2
  BigRat div2() const { BigRat r; r.div_2exp(*this, 1U); return r; }
  /// intValue
  int intValue() const { return static_cast<int>(doubleValue()); }
  /// longValue
  long longValue() const { return static_cast<long>(doubleValue()); }
  /// doubleValue
  double doubleValue() const { return get_d(); }
  /// BigIntValue
  BigInt BigIntValue() const { return get_z(); }
  //@}
#endif
};

/// \addtogroup BigRatArithmeticOperators
//@{
/// BigRat + BigRat
inline BigRat operator+(const BigRat& x, const BigRat& y)
{ BigRat r; r.add(x, y); return r; }
/// BigRat + int
inline BigRat operator+(const BigRat& x, int y)
{ return x + BigRat(y); }
/// int + BigRat
inline BigRat operator+(int x, const BigRat& y)
{ return BigRat(x) + y; }
/// BigRat + unsigned int
inline BigRat operator+(const BigRat& x, unsigned int y)
{ return x + BigRat(y); }
/// unsigned int + BigRat
inline BigRat operator+(unsigned int x, const BigRat& y)
{ return BigRat(x) + y; }
/// BigRat + long
inline BigRat operator+(const BigRat& x, long y)
{ return x + BigRat(y); }
/// long + BigRat
inline BigRat operator+(long x, const BigRat& y)
{ return BigRat(x) + y; }
/// BigRat + unsigned long
inline BigRat operator+(const BigRat& x, unsigned long y)
{ return x + BigRat(y); }
/// unsigned long + BigRat
inline BigRat operator+(unsigned long x, const BigRat& y)
{ return BigRat(x) + y; }
/// BigRat + double
inline BigRat operator+(const BigRat& x, double y)
{ return x + BigRat(y); }
/// double + BigRat
inline BigRat operator+(double x, const BigRat& y)
{ return BigRat(x) + y; }
/// BigRat + BigInt
inline BigRat operator+(const BigRat& x, const BigInt& y)
{ return x + BigRat(y); }
/// BigInt + BigRat
inline BigRat operator+(const BigInt& x, const BigRat& y)
{ return BigRat(x) + y; }

/// BigRat - BigRat
inline BigRat operator-(const BigRat& x, const BigRat& y)
{ BigRat r; r.sub(x, y); return r; }
/// BigRat - int
inline BigRat operator-(const BigRat& x, int y)
{ return x - BigRat(y); }
/// int - BigRat
inline BigRat operator-(int x, const BigRat& y)
{ return BigRat(x) - y; }
/// BigRat - unsigned int
inline BigRat operator-(const BigRat& x, unsigned int y)
{ return x - BigRat(y); }
/// unsigned int - BigRat
inline BigRat operator-(unsigned int x, const BigRat& y)
{ return BigRat(x) - y; }
/// BigRat - long
inline BigRat operator-(const BigRat& x, long y)
{ return x - BigRat(y); }
/// long - BigRat
inline BigRat operator-(long x, const BigRat& y)
{ return BigRat(x) - y; }
/// BigRat - unsigned long
inline BigRat operator-(const BigRat& x, unsigned long y)
{ return x - BigRat(y); }
/// unsigned long - BigRat
inline BigRat operator-(unsigned long x, const BigRat& y)
{ return BigRat(x) - y; }
/// BigRat - double
inline BigRat operator-(const BigRat& x, double y)
{ return x - BigRat(y); }
/// double - BigRat
inline BigRat operator-(double x, const BigRat& y)
{ return BigRat(x) - y; }
/// BigRat - BigInt
inline BigRat operator-(const BigRat& x, const BigInt& y)
{ return x - BigRat(y); }
/// BigInt - BigRat
inline BigRat operator-(const BigInt& x, const BigRat& y)
{ return BigRat(x) - y; }

/// BigRat * BigRat
inline BigRat operator*(const BigRat& x, const BigRat& y)
{ BigRat r; r.mul(x, y); return r; }
/// BigRat * int
inline BigRat operator*(const BigRat& x, int y)
{ return x * BigRat(y); }
/// int * BigRat
inline BigRat operator*(int x, const BigRat& y)
{ return BigRat(x) * y; }
/// BigRat * unsigned int
inline BigRat operator*(const BigRat& x, unsigned int y)
{ return x * BigRat(y); }
/// unsigned int * BigRat
inline BigRat operator*(unsigned int x, const BigRat& y)
{ return BigRat(x) * y; }
/// BigRat * long
inline BigRat operator*(const BigRat& x, long y)
{ return x * BigRat(y); }
/// long * BigRat
inline BigRat operator*(long x, const BigRat& y)
{ return BigRat(x) * y; }
/// BigRat * unsigned long
inline BigRat operator*(const BigRat& x, unsigned long y)
{ return x * BigRat(y); }
/// unsigned long * BigRat
inline BigRat operator*(unsigned long x, const BigRat& y)
{ return BigRat(x) * y; }
/// BigRat * double
inline BigRat operator*(const BigRat& x, double y)
{ return x * BigRat(y); }
/// double * BigRat
inline BigRat operator*(double x, const BigRat& y)
{ return BigRat(x) * y; }
/// BigRat * BigInt
inline BigRat operator*(const BigRat& x, const BigInt& y)
{ return x * BigRat(y); }
/// BigInt * BigRat
inline BigRat operator*(const BigInt& x, const BigRat& y)
{ return BigRat(x) * y; }

/// BigRat / BigRat
inline BigRat operator/(const BigRat& x, const BigRat& y)
{ BigRat r; r.div(x, y); return r; }
/// BigRat / int
inline BigRat operator/(const BigRat& x, int y)
{ return x / BigRat(y); }
/// int / BigRat
inline BigRat operator/(int x, const BigRat& y)
{ return BigRat(x) / y; }
/// BigRat / unsigned int
inline BigRat operator/(const BigRat& x, unsigned int y)
{ return x / BigRat(y); }
/// unsigned int / BigRat
inline BigRat operator/(unsigned int x, const BigRat& y)
{ return BigRat(x) / y; }
/// BigRat / long
inline BigRat operator/(const BigRat& x, long y)
{ return x / BigRat(y); }
/// long / BigRat
inline BigRat operator/(long x, const BigRat& y)
{ return BigRat(x) / y; }
/// BigRat / unsigned long
inline BigRat operator/(const BigRat& x, unsigned long y)
{ return x / BigRat(y); }
/// unsigned long / BigRat
inline BigRat operator/(unsigned long x, const BigRat& y)
{ return BigRat(x) / y; }
/// BigRat / double
inline BigRat operator/(const BigRat& x, double y)
{ return x / BigRat(y); }
/// double / BigRat
inline BigRat operator/(double x, const BigRat& y)
{ return BigRat(x) / y; }
/// BigRat / BigInt
inline BigRat operator/(const BigRat& x, const BigInt& y)
{ return x / BigRat(y); }
/// BigInt / BigRat
inline BigRat operator/(const BigInt& x, const BigRat& y)
{ return BigRat(x) / y; }
//@}

/// BigRat << int
inline BigRat operator<<(const BigRat& x, int y)
{ BigRat r; r.mul_2exp(x, y); return r; }
/// BigRat << unsigned int
inline BigRat operator<<(const BigRat& x, unsigned int y)
{ BigRat r; r.mul_2exp(x, y); return r; }
/// BigRat << long
inline BigRat operator<<(const BigRat& x, long y)
{ BigRat r; r.mul_2exp(x, y); return r; }
/// BigRat << unsigned long
inline BigRat operator<<(const BigRat& x, unsigned long y)
{ BigRat r; r.mul_2exp(x, y); return r; }
/// BigRat >> int
inline BigRat operator>>(const BigRat& x, int y)
{ BigRat r; r.div_2exp(x, y); return r; }
/// BigRat >> unsigned int
inline BigRat operator>>(const BigRat& x, unsigned int y)
{ BigRat r; r.div_2exp(x, y); return r; }
/// BigRat >> long
inline BigRat operator>>(const BigRat& x, long y)
{ BigRat r; r.div_2exp(x, y); return r; }
/// BigRat >> unsigned long
inline BigRat operator>>(const BigRat& x, unsigned long y)
{ BigRat r; r.div_2exp(x, y); return r; }
//@}

/// \addtogroup BigRatComparisonOperators
//@{
/// BigRat == BigRat
inline bool operator==(const BigRat& x, const BigRat& y)
{ return x.equal(y) != 0; }
/// BigRat == int
inline bool operator==(const BigRat& x, int y)
{ return x.cmp(y, 1UL) == 0; }
/// int == BigRat
inline bool operator==(int x, const BigRat& y)
{ return y.cmp(x, 1UL) == 0; }
/// BigRat == unsigned int
inline bool operator==(const BigRat& x, unsigned int y)
{ return x.cmp(y, 1UL) == 0; }
/// unsigned int == BigRat
inline bool operator==(unsigned int x, const BigRat& y)
{ return y.cmp(x, 1UL) == 0; }
/// BigRat == long
inline bool operator==(const BigRat& x, long y)
{ return x.cmp(y, 1UL) == 0; }
/// long == BigRat
inline bool operator==(long x, const BigRat& y)
{ return y.cmp(x, 1UL) == 0; }
/// BigRat == unsigned long
inline bool operator==(const BigRat& x, unsigned long y)
{ return x.cmp(y, 1UL) == 0; }
/// unsigned long == BigRat
inline bool operator==(unsigned long x, const BigRat& y)
{ return y.cmp(x, 1UL) == 0; }
/// BigRat == double
inline bool operator==(const BigRat& x, double y)
{ return x.cmp(y) == 0; }
/// double == BigRat
inline bool operator==(double x, const BigRat& y)
{ return y.cmp(x) == 0; }
/// BigRat == BigInt
inline bool operator==(const BigRat& x, const BigInt& y)
{ return x == BigRat(y);}
/// BigInt == BigRat
inline bool operator==(const BigInt& x, const BigRat& y)
{ return BigRat(x) ==  y;}
 
/// BigRat != BigRat
inline bool operator!=(const BigRat& x, const BigRat& y)
{ return x.equal(y) == 0; }
/// BigRat != int
inline bool operator!=(const BigRat& x, int y)
{ return x.cmp(y, 1UL) != 0; }
/// int != BigRat
inline bool operator!=(int x, const BigRat& y)
{ return y.cmp(x, 1UL) != 0; }
/// BigRat != unsigned int
inline bool operator!=(const BigRat& x, unsigned int y)
{ return x.cmp(y, 1UL) != 0; }
/// unsigned int != BigRat
inline bool operator!=(unsigned int x, const BigRat& y)
{ return y.cmp(x, 1UL) != 0; }
/// BigRat != long
inline bool operator!=(const BigRat& x, long y)
{ return x.cmp(y, 1UL) != 0; }
/// long != BigRat
inline bool operator!=(long x, const BigRat& y)
{ return y.cmp(x, 1UL) != 0; }
/// BigRat != unsigned long
inline bool operator!=(const BigRat& x, unsigned long y)
{ return x.cmp(y, 1UL) != 0; }
/// unsigned long != BigRat
inline bool operator!=(unsigned long x, const BigRat& y)
{ return y.cmp(x, 1UL) != 0; }
/// BigRat != double
inline bool operator!=(const BigRat& x, double y)
{ return x.cmp(y) != 0; }
/// double != BigRat
inline bool operator!=(double x, const BigRat& y)
{ return y.cmp(x) != 0; }
/// BigRat != BigInt
inline bool operator!=(const BigRat& x, const BigInt& y)
{ return x != BigRat(y);}
/// BigInt != BigRat
inline bool operator!=(const BigInt& x, const BigRat& y)
{ return BigRat(x) !=  y;}
 
/// BigRat >= BigRat
inline bool operator>=(const BigRat& x, const BigRat& y)
{ return x.cmp(y) >= 0; }
/// BigRat >= int
inline bool operator>=(const BigRat& x, int y)
{ return x.cmp(y, 1UL) >= 0; }
/// int >= BigRat
inline bool operator>=(int x, const BigRat& y)
{ return y.cmp(x, 1UL) <= 0; }
/// BigRat >= unsigned int
inline bool operator>=(const BigRat& x, unsigned int y)
{ return x.cmp(y, 1UL) >= 0; }
/// unsigned int >= BigRat
inline bool operator>=(unsigned int x, const BigRat& y)
{ return y.cmp(x, 1UL) <= 0; }
/// BigRat >= long
inline bool operator>=(const BigRat& x, long y)
{ return x.cmp(y, 1UL) >= 0; }
/// long >= BigRat
inline bool operator>=(long x, const BigRat& y)
{ return y.cmp(x, 1UL) <= 0; }
/// BigRat >= unsigned long
inline bool operator>=(const BigRat& x, unsigned long y)
{ return x.cmp(y, 1UL) >= 0; }
/// unsigned long >= BigRat
inline bool operator>=(unsigned long x, const BigRat& y)
{ return y.cmp(x, 1UL) <= 0; }
/// BigRat >= double
inline bool operator>=(const BigRat& x, double y)
{ return x.cmp(y) >= 0; }
/// double >= BigRat
inline bool operator>=(double x, const BigRat& y)
{ return y.cmp(x) <= 0; }
/// BigRat >= BigInt
inline bool operator>=(const BigRat& x, const BigInt& y)
{ return x >= BigRat(y);}
/// BigInt >= BigRat
inline bool operator>=(const BigInt& x, const BigRat& y)
{ return BigRat(x) >=  y;}
 
/// BigRat <= BigRat
inline bool operator<=(const BigRat& x, const BigRat& y)
{ return x.cmp(y) <= 0; }
/// BigRat <= int
inline bool operator<=(const BigRat& x, int y)
{ return x.cmp(y, 1UL) <= 0; }
/// int <= BigRat
inline bool operator<=(int x, const BigRat& y)
{ return y.cmp(x, 1UL) >= 0; }
/// BigRat <= unsigned int
inline bool operator<=(const BigRat& x, unsigned int y)
{ return x.cmp(y, 1UL) <= 0; }
/// unsigned int <= BigRat
inline bool operator<=(unsigned int x, const BigRat& y)
{ return y.cmp(x, 1UL) >= 0; }
/// BigRat <= long
inline bool operator<=(const BigRat& x, long y)
{ return x.cmp(y, 1UL) <= 0; }
/// long <= BigRat
inline bool operator<=(long x, const BigRat& y)
{ return y.cmp(x, 1UL) >= 0; }
/// BigRat <= unsigned long
inline bool operator<=(const BigRat& x, unsigned long y)
{ return x.cmp(y, 1UL) <= 0; }
/// unsigned long <= BigRat
inline bool operator<=(unsigned long x, const BigRat& y)
{ return y.cmp(x, 1UL) >= 0; }
/// BigRat <= double
inline bool operator<=(const BigRat& x, double y)
{ return x.cmp(y) <= 0; }
/// double <= BigRat
inline bool operator<=(double x, const BigRat& y)
{ return y.cmp(x) >= 0; }
/// BigRat <= BigInt
inline bool operator<=(const BigRat& x, const BigInt& y)
{ return x <= BigRat(y);}
/// BigInt <= BigRat
inline bool operator<=(const BigInt& x, const BigRat& y)
{ return BigRat(x) <=  y;}
 
/// BigRat > BigRat
inline bool operator>(const BigRat& x, const BigRat& y)
{ return x.cmp(y) > 0; }
/// BigRat > int
inline bool operator>(const BigRat& x, int y)
{ return x.cmp(y, 1UL) > 0; }
/// int > BigRat
inline bool operator>(int x, const BigRat& y)
{ return y.cmp(x, 1UL) < 0; }
/// BigRat > unsigned int
inline bool operator>(const BigRat& x, unsigned int y)
{ return x.cmp(y, 1UL) > 0; }
/// unsigned int > BigRat
inline bool operator>(unsigned int x, const BigRat& y)
{ return y.cmp(x, 1UL) < 0; }
/// BigRat > long
inline bool operator>(const BigRat& x, long y)
{ return x.cmp(y, 1UL) > 0; }
/// long > BigRat
inline bool operator>(long x, const BigRat& y)
{ return y.cmp(x, 1UL) < 0; }
/// BigRat > unsigned long
inline bool operator>(const BigRat& x, unsigned long y)
{ return x.cmp(y, 1UL) > 0; }
/// unsigned long > BigRat
inline bool operator>(unsigned long x, const BigRat& y)
{ return y.cmp(x, 1UL) < 0; }
/// BigRat > double
inline bool operator>(const BigRat& x, double y)
{ return x.cmp(y) > 0; }
/// double > BigRat
inline bool operator>(double x, const BigRat& y)
{ return y.cmp(x) < 0; }
/// BigRat > BigInt
inline bool operator>(const BigRat& x, const BigInt& y)
{ return x > BigRat(y);}
/// BigInt > BigRat
inline bool operator>(const BigInt& x, const BigRat& y)
{ return BigRat(x) >  y;}
 
/// BigRat < BigRat
inline bool operator<(const BigRat& x, const BigRat& y)
{ return x.cmp(y) < 0; }
/// BigRat < int
inline bool operator<(const BigRat& x, int y)
{ return x.cmp(y, 1UL) < 0; }
/// int < BigRat
inline bool operator<(int x, const BigRat& y)
{ return y.cmp(x, 1UL) > 0; }
/// BigRat < unsigned int
inline bool operator<(const BigRat& x, unsigned int y)
{ return x.cmp(y, 1UL) < 0; }
/// unsigned int < BigRat
inline bool operator<(unsigned int x, const BigRat& y)
{ return y.cmp(x, 1UL) > 0; }
/// BigRat < long
inline bool operator<(const BigRat& x, long y)
{ return x.cmp(y, 1UL) < 0; }
/// long < BigRat
inline bool operator<(long x, const BigRat& y)
{ return y.cmp(x, 1UL) > 0; }
/// BigRat < unsigned long
inline bool operator<(const BigRat& x, unsigned long y)
{ return x.cmp(y, 1UL) < 0; }
/// unsigned long < BigRat
inline bool operator<(unsigned long x, const BigRat& y)
{ return y.cmp(x, 1UL) > 0; }
/// BigRat < double
inline bool operator<(const BigRat& x, double y)
{ return x.cmp(y) < 0; }
/// double < BigRat
inline bool operator<(double x, const BigRat& y)
{ return y.cmp(x) > 0; }
/// BigRat < BigInt
inline bool operator<(const BigRat& x, const BigInt& y)
{ return x < BigRat(y);}
/// BigInt < BigRat
inline bool operator<(const BigInt& x, const BigRat& y)
{ return BigRat(x) <  y;}
//@}
 //@}
 
/// \addtogroup BigRatIostreamOperators 
//@{
/// istream operator for <tt>BigRat</tt>
inline std::istream& operator>>(std::istream& is, BigRat& x)
{ return ::extract(is, x.mp()); }
/// ostream operator for <tt>BigRat</tt>
inline std::ostream& operator<<(std::ostream& os, const BigRat& x)
{ return os << x.get_str(); }
//@}

#ifndef CORE_DISABLE_OLDNAMES 
/// \addtogroup BigRatBackCompatiableFunctions
//@{
/// comparison
inline int cmp(const BigRat& x, const BigRat& y) { return x.cmp(y); }
/// sign
inline int sign(const BigRat& a) { return a.sgn(); }
inline int sgn(const BigRat& a) { return a.sgn(); }
/// abs
inline BigRat abs(const BigRat& a) { BigRat r; r.abs(a); return r; }
/// neg
inline BigRat neg(const BigRat& a) { BigRat r; r.neg(a); return r; }
/// div2
inline BigRat div2(const BigRat& a) { BigRat r(a); return r.div2(); }
/// numerator
inline BigInt numerator(const BigRat& a) { return a.numerator(); }
/// denominator
inline BigInt denominator(const BigRat& a) { return a.denominator(); }
/// longValue
inline long longValue(const BigRat& a) { return a.longValue(); }
/// doubleValue
inline double doubleValue(const BigRat& a) { return a.doubleValue(); }
/// return BigInt value
inline BigInt BigIntValue(const BigRat& a) { return a.BigIntValue(); }
/// exact division
inline BigRat div_exact(const BigRat& x, const BigRat& y) 
{ BigRat z; z.div(x, y); return z; }
/// gcd
/// \brief we may want replace this by the definition of gcd of a quotient 
///        field of a UFD [Yap's book, Chap.3]
inline BigRat gcd(const BigRat&, const BigRat&)
{ return BigRat(1);  }
/// is integer ?
inline bool isInteger(const BigRat& x) { return x.den() == 1; }
// Chee: 8/8/2004: need isDivisible to compile Polynomial<BigRat>
// A trivial implementation is to return true always. But this
// caused tPolyRat to fail.
// So we follow the definition of
// Expr::isDivisible(e1, e2) which checks if e1/e2 is an integer.  
/// is divisible ?
inline bool isDivisible(const BigRat& x, const BigRat& y) 
{ BigRat r; r.div(x, y); return isInteger(r); }
/// floorlg
inline long floorLg(const BigRat& a) { return a.lMSB(); }
inline long floorlg(const BigRat& a) { return a.lMSB(); }
/// ceillg
inline long ceilLg(const BigRat& a) { return a.uMSB(); }
inline long ceillg(const BigRat& a) { return a.uMSB(); }
//@}
#endif

CORE_END_NAMESPACE

#endif /*__CORE_BIGRAT_H__*/
