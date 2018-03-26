/****************************************************************************
 * BigFloat.h -- Big Floating-point number class based on mpfr in MPFR
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
 * $Id: BigFloat.h,v 1.49 2010/11/23 17:58:36 exact Exp $
 ***************************************************************************/
#ifndef __CORE_BIGFLOAT_H__
#define __CORE_BIGFLOAT_H__

#include <CORE/Mpfr.h>
#include <CORE/CoreIO.h>
#include <CORE/BigRat.h>
#include <string>
#include <cstring>
#include <iostream>
#include <cassert>
#include <cmath>
#include <stdio.h>

// cygwin defines log2 as a macro
// we need to undefine log2 to have our own log2
#undef log2

/* Known Issues:

  1. mpfr_set_str() returns 0 if the string is a valid number, otherwise 
     returns -1. So you cannot get the exactness of the results.
 */

// Invert a round mode
#define INVERT_RND(rnd) \
((rnd == GMP_RNDU) ? GMP_RNDD : ((rnd == GMP_RNDD) ? GMP_RNDU : rnd))
// Invert sign
inline void invert_sgn(mpfr_ptr x)
{ MPFR_SIGN(x) = -MPFR_SIGN(x); }

/* subtraction */
inline int mpfr_z_sub(mpfr_ptr z, mpz_srcptr x, mpfr_srcptr y, mp_rnd_t rnd)
{ int r = -mpfr_sub_z(z, y, x, INVERT_RND(rnd)); invert_sgn(z); return r; }
inline int mpfr_q_sub(mpfr_ptr z, mpq_srcptr x, mpfr_srcptr y, mp_rnd_t rnd)
{ int r = -mpfr_sub_q(z, y, x, INVERT_RND(rnd)); invert_sgn(z); return r; }

/* remove trailing zeros (by limbs) */
void mpfr_remove_trailing_zeros(mpfr_t x);
/* C++-style input of mpfr */
std::istream& operator>> (std::istream &, mpfr_ptr);
/* convert mpfr to string using specified ndigits, base, fmt, rnd, etc.*/
std::string mpfr2str(mpfr_srcptr, size_t, int, int, rnd_t, bool, bool, bool);

CORE_BEGIN_NAMESPACE

// constant of default precision for integer, IEEE single and double
const size_t INT_PREC = sizeof(int)*8;
const size_t SINGLE_PREC = 24;
const size_t DOUBLE_PREC = 53;
 
#ifndef CORE_DISABLE_REFCOUNTING
  typedef RcMpfr BigFloatBase;
#else
  typedef Mpfr BigFloatBase;
#endif

/// \class BigFloat BigFloat.h
/// \brief BigFloat is a big floating-point number class based on MPFR
class BigFloat : public BigFloatBase {
  typedef BigFloatBase base_cls;
public: // public typedefs
  typedef BigInt ZT;
  typedef BigRat QT;
public:
  /// \name constructors (auto version)
  //@{
  /// default constructor
  BigFloat() : base_cls() {}
  /// copy constructor
  BigFloat(const BigFloat& rhs) : base_cls(rhs) {}
  /// constructor for <tt>int</tt> (use INT_PREC by default)
  BigFloat(int i, rnd_t rnd = MPFR_RND)
    : base_cls(static_cast<long>(i), INT_PREC, rnd) {}
  /// constructor for <tt>unsigned int</tt> (use INT_PREC by default)
  BigFloat(unsigned int i, rnd_t rnd = MPFR_RND)
    : base_cls(static_cast<unsigned long>(i), INT_PREC, rnd) {}
  /// constructor for <tt>long</tt> (use INT_PREC by default)
  BigFloat(long i, rnd_t rnd = MPFR_RND)
    : base_cls(i, INT_PREC, rnd) {}
  /// constructor for <tt>unsigned long</tt> (use INT_PREC by default)
  BigFloat(unsigned long i, rnd_t rnd = MPFR_RND)
    : base_cls(i, INT_PREC, rnd) {}
  /// constructor for <tt>double`</tt> (use DOUBLE_PREC by default)
  BigFloat(double i, rnd_t rnd = MPFR_RND)
    : base_cls(i, DOUBLE_PREC, rnd) {}
  /// constructor for <tt>BigInt</tt> 
  BigFloat(const BigInt& x, rnd_t rnd = MPFR_RND)
    : base_cls(x.mp(), (std::max)(count_prec(x), (prec_t)32), rnd) {}
  /// constructor for <tt>BigRat</tt> (use DOUBLE_PREC by default)
  BigFloat(const BigRat& x, rnd_t rnd = MPFR_RND)
    : base_cls(x.mp(), (std::max)(count_prec(x), (prec_t)32), rnd) {}

#ifndef CORE_LEVEL_1_NO_WRAPPERS  
BigFloat(const DoubleWrapper &wrap, rnd_t rnd = MPFR_RND)
 	: base_cls(wrap.doubleValue(), DOUBLE_PREC, rnd) {
  }
#endif

  /// constructor for <tt>char*</tt> (no implicit conversion)
  explicit BigFloat(const char* s, int base = 10, rnd_t rnd = MPFR_RND)
    : base_cls(s, base, (std::max)(count_prec(s), (prec_t)53), rnd) {}
  /// constructor for <tt>char</tt> (no implicit conversion)
  explicit BigFloat(const std::string& s, int base = 10, rnd_t rnd = MPFR_RND)
    : base_cls(s.c_str(), base,
      (std::max)(count_prec(s.c_str()), (prec_t)53), rnd) {}
  //@}

  /// \name constructors (fixed version)
  //@{
  /// constructor for <tt>short</tt> with specified precision
  BigFloat(int i, prec_t prec, rnd_t rnd = MPFR_RND)
    : base_cls(static_cast<long>(i), prec, rnd) {}
  /// constructor for <tt>unsigned short</tt> with specified precision
  BigFloat(unsigned int i, prec_t prec, rnd_t rnd = MPFR_RND)
    : base_cls(static_cast<unsigned long>(i), prec, rnd) {}
  /// constructor for <tt>long</tt> with specified precision
  BigFloat(long i, prec_t prec, rnd_t rnd = MPFR_RND)
    : base_cls(i, prec, rnd) {}
  /// constructor for <tt>unsigned long</tt> with specified precision
  BigFloat(unsigned long i, prec_t prec, rnd_t rnd = MPFR_RND)
    : base_cls(i, prec, rnd) {}
  /// constructor for <tt>double</tt> with specified precision
  BigFloat(double i, prec_t prec, rnd_t rnd = MPFR_RND)
    : base_cls(i, prec, rnd) {}
  /// constructor for <tt>BigInt</tt> with specified precision
  BigFloat(const BigInt& x, prec_t prec, rnd_t rnd = MPFR_RND)
    : base_cls(x.mp(), prec, rnd) {}
  /// constructor for <tt>BigRat</tt> with specified precision
  BigFloat(const BigRat& x, prec_t prec, rnd_t rnd = MPFR_RND)
    : base_cls(x.mp(), prec, rnd) {}
  /// constructor for <tt>BigFloat</tt> with specified precision
  BigFloat(const BigFloat& x, prec_t prec, rnd_t rnd = MPFR_RND)
    : base_cls(x.mp(), prec, rnd) {}
  /// constructor for <tt>char*</tt> with specified precision
  explicit BigFloat(const char* s, int base, prec_t prec, rnd_t rnd = MPFR_RND)
    : base_cls(s, base, prec, rnd) {}
  /// constructor for <tt>std::string</tt> with specified precision
  explicit BigFloat(const std::string& s,int base,prec_t prec,rnd_t rnd=MPFR_RND)
    : base_cls(s.c_str(), base, prec, rnd) {}

  /// constructor with value \f$i*2^e\f$ for <tt>int</tt>
  BigFloat(int i, exp_t e, prec_t prec, rnd_t rnd = MPFR_RND)
    : base_cls(static_cast<long>(i), e, prec, rnd) {}
  /// constructor with value \f$i*2^e\f$ for <tt>unsigned int</tt>
  BigFloat(unsigned int i, exp_t e, prec_t prec, rnd_t rnd = MPFR_RND)
    : base_cls(static_cast<unsigned long>(i), e, prec, rnd) {}
  /// constructor with value \f$i*2^e\f$ for <tt>long</tt>
  BigFloat(long i, exp_t e, prec_t prec, rnd_t rnd = MPFR_RND)
    : base_cls(i, e, prec, rnd) {}
  /// constructor with value \f$i*2^e\f$ for <tt>unsigned long</tt>
  BigFloat(unsigned long i, exp_t e, prec_t prec, rnd_t rnd = MPFR_RND)
    : base_cls(i, e, prec, rnd) {}
  //@}

public:
  /// \name precision accessors
  //@{
  /// return current precision (bit length of mantissa)
  prec_t get_prec() const
  { return mpfr_get_prec(mp())-1; }
  /// set current precision (bit length of mantissa). Previous value lost.
  void set_prec(prec_t prec) {
    if(get_prec() == prec) return;
    prec += 1;
    prec = (std::min)((std::max)(prec, prec_t(MPFR_PREC_MIN)), prec_t(MPFR_PREC_MAX));
    mpfr_set_prec(mp(), prec);
  }
  /// round the current precision to prec, using specified rounding mode.
  int prec_round(prec_t prec, rnd_t rnd = MPFR_RND) {
    prec += 1;
    prec = (std::min)((std::max)(prec, prec_t(MPFR_PREC_MIN)), prec_t(MPFR_PREC_MAX));
    return mpfr_prec_round(mp(), prec, rnd);
  } 
  //@}

  /// \name exponent accessors
  //@{
  /// return exponent of Mpfr
  // internal represenation of Mpfr is
  // (-1)^s * m * 2^e
  // m is the mantissa with with 0.5 < m < 1, e is the sxponent
  exp_t get_exp() const
  { return sgn() ? mpfr_get_exp(mp()) : 0; }
  //@}

public:
  /// \name assignment functions (raw version)
  //@{
  /// assignment functions for <tt>BigFloat</tt>
  int r_set(const BigFloat& rhs, rnd_t rnd = MPFR_RND)
  { return mpfr_set(mp(), rhs.mp(), rnd); }
  /// assignment functions for <tt>int</tt>
  int r_set(int i, rnd_t rnd = MPFR_RND)
  { return mpfr_set_si(mp(), i, rnd); }
  /// assignment functions for <tt>unsigned int</tt>
  int r_set(unsigned int i, rnd_t rnd = MPFR_RND)
  { return mpfr_set_ui(mp(), i, rnd); }
  /// assignment functions for <tt>long</tt>
  int r_set(long i, rnd_t rnd = MPFR_RND)
  { return mpfr_set_si(mp(), i, rnd); }
  /// assignment functions for <tt>unsigned long</tt>
  int r_set(unsigned long i, rnd_t rnd = MPFR_RND)
  { return mpfr_set_ui(mp(), i, rnd); }
  /// assignment functions for <tt>double</tt>
  int r_set(double i, rnd_t rnd = MPFR_RND)
  { return mpfr_set_d(mp(), i, rnd); }
  /// assignment functions for <tt>BigInt</tt>
  int r_set(const BigInt& x, rnd_t rnd = MPFR_RND)
  { return mpfr_set_z(mp(), x.mp(), rnd); }
  /// assignment functions for <tt>BigRat</tt>
  int r_set(const BigRat& x, rnd_t rnd = MPFR_RND) // use DOUBLE_PREC
  { return mpfr_set_q(mp(), x.mp(), rnd); }
  /// assignment functions for <tt>char*</tt>
  int r_set(const char* str, int base = 10, rnd_t rnd = MPFR_RND)
  { return mpfr_set_str(mp(), str, base, rnd); }
  /// assignment functions for <tt>std::string</tt>
  int r_set(const std::string& str, int base = 10, rnd_t rnd = MPFR_RND)
  { return mpfr_set_str(mp(), str.c_str(), base, rnd); }

  /// set value to be \f$i*2^e\f$ for <tt>int</tt>
  int r_set_2exp(int i, exp_t e, rnd_t rnd = MPFR_RND)
  { return mpfr_set_si_2exp(mp(), i, e, rnd); }
  /// set value to be \f$i*2^e\f$ for <tt>unsigned int</tt>
  int r_set_2exp(unsigned int i, exp_t e, rnd_t rnd = MPFR_RND)
  { return mpfr_set_ui_2exp(mp(), i, e, rnd); }
  /// set value to be \f$i*2^e\f$ for <tt>long</tt>
  int r_set_2exp(long i, exp_t e, rnd_t rnd = MPFR_RND)
  { return mpfr_set_si_2exp(mp(), i, e, rnd); }
  /// set value to be \f$i*2^e\f$ for <tt>unsigned long</tt>
  int r_set_2exp(unsigned long i, exp_t e, rnd_t rnd = MPFR_RND)
  { return mpfr_set_ui_2exp(mp(), i, e, rnd); }
  //@}

  /// \name assignment functions (fixed version)
  //@{
  /// assignment functions for <tt>T</tt>
  template <typename T> int set(const T& x, prec_t prec, rnd_t rnd = MPFR_RND)
  { set_prec(prec); return r_set(x, rnd); }
  /// assignment functions for <tt>char*</tt>
  int set(const char* x, int base, prec_t prec, rnd_t rnd = MPFR_RND)
  { set_prec(prec); return r_set(x, base, rnd); }
  /// assignment functions for <tt>std::string</tt>
  int set(const std::string& x, int base, prec_t prec, rnd_t rnd = MPFR_RND)
  { set_prec(prec); return r_set(x.c_str(), base, rnd); }
  /// set value to be \f$i*2^e\f$ for <tt>T</tt>
  template <typename T> 
  int set_2exp(const T& x, exp_t e, prec_t prec, rnd_t rnd = MPFR_RND)
  { set_prec(prec); return r_set_2exp(x, e, rnd); }
  //@}

  /// \name assignment functions (auto version)
  //@{
  /// assignment functions for <tt>T</tt>
  template <typename T> int set(const T& x, rnd_t rnd = MPFR_RND)
  { set_prec(count_prec(x)); return r_set(x, rnd); }
  /// assignment functions for <tt>char*</tt>
  int set(const char* x, int base = 10, rnd_t rnd = MPFR_RND)
  { set_prec(count_prec(x)); return r_set(x, base, rnd); }
  /// assignment functions for <tt>std::string</tt>
  int set(const std::string& x, int base = 10, rnd_t rnd = MPFR_RND)
  { set_prec(count_prec(x)); return r_set(x, base, rnd); }
  /// set value to be \f$i*2^e\f$ for <tt>T</tt>
  template <typename T> 
  int set_2exp(const T& x, exp_t e, rnd_t rnd = MPFR_RND)
  { set_prec(count_prec(x)); return r_set_2exp(x, e, rnd); }
  //@}

public:
  /// \name arithmetic functions -- addition (raw version)
  //@{
  /// addition for <tt>BigFloat+BigFloat</tt>  (raw version)
  int r_add(const BigFloat& x, const BigFloat& y, rnd_t rnd = MPFR_RND)
  { return mpfr_add(mp(), x.mp(), y.mp(), rnd); }
  /// addition for <tt>BigFloat+int</tt> (raw version)
  int r_add(const BigFloat& x, int y, rnd_t rnd = MPFR_RND)
  { return r_add(x, static_cast<long>(y), rnd); }
  /// addition for <tt>BigFloat+unsigned int</tt> (raw version)
  int r_add(const BigFloat& x, unsigned int y, rnd_t rnd = MPFR_RND)
  { return r_add(x, static_cast<unsigned long>(y), rnd); }
  /// addition for <tt>BigFloat+long</tt> (raw version)
  int r_add(const BigFloat& x, long y, rnd_t rnd = MPFR_RND)
  { return mpfr_add_si(mp(), x.mp(), y, rnd); }
  /// addition for <tt>BigFloat+unsigned long</tt> (raw version)
  int r_add(const BigFloat& x, unsigned long y, rnd_t rnd = MPFR_RND)
  { return mpfr_add_ui(mp(), x.mp(), y, rnd); }
  /// addition for <tt>BigFloat+double</tt> (raw version)
  int r_add(const BigFloat& x, double y, rnd_t rnd = MPFR_RND)
  { return r_add(x, BigFloat(y), rnd); }
  /// addition for <tt>BigFloat+BigInt</tt> (raw version)
  int r_add(const BigFloat& x, const BigInt& y, rnd_t rnd = MPFR_RND)
  { return mpfr_add_z(mp(), x.mp(), y.mp(), rnd); }
  /// addition for <tt>BigFloat+BigRat</tt> (raw version)
  int r_add(const BigFloat& x, const BigRat& y, rnd_t rnd = MPFR_RND)
  { return mpfr_add_q(mp(), x.mp(), y.mp(), rnd); }
  /// addition for <tt>T+BigFloat</tt> (raw version)
  template <typename T> 
  int r_add(const T& x, const BigFloat& y, rnd_t rnd = MPFR_RND)
  { return r_add(y, x, rnd); }
  //@}

  /// \name arithmetic functions -- addition (fixed version)
  //@{
  /// addition for <tt>BigFloat+BigFloat</tt> (fixed version)
  int add(const BigFloat& x, const BigFloat& y,prec_t prec,rnd_t rnd=MPFR_RND){
    if (&x == this || &y == this) { // if one of inputs are output
       if (prec > get_prec()) {
         prec_round (prec, rnd);
         return r_add(x,y,rnd);
       } else if (prec < get_prec()) {
         r_add(x,y,rnd);
         return prec_round (prec, rnd);
       } else {
         return r_add(x,y,rnd);
       }
/*
       BigFloat result(0, prec);
       int r = result.r_add(x, y, rnd);
       swap(result); return r;
*/
     } else {
       set_prec(prec); return r_add(x, y, rnd); 
     }
  }
  /// addition for <tt>BigFloat+T</tt> (fixed version)
  template <typename T> 
  int add(const BigFloat& x, const T& y, prec_t prec, rnd_t rnd = MPFR_RND) {
    if (&x == this) { // if x is same as output
//*
      if (prec > get_prec()) {
         prec_round (prec, rnd);
         return r_add(x,y,rnd);
       } else if (prec < get_prec()) {
         r_add(x,y,rnd);
         return prec_round (prec, rnd);
       } else {
         return r_add(x,y,rnd);
       }
/*/
       BigFloat result(0, prec);
       int r = result.r_add(x, y, rnd);
       swap(result); return r;
//*/
	} else {
       set_prec(prec); return r_add(x, y, rnd); 
     }
  }
  /// addition for <tt>T+BigFloat</tt> (fixed version)
  template <typename T> 
  int add(const T& x, const BigFloat& y, prec_t prec, rnd_t rnd = MPFR_RND)
  { return add(y, x, prec, rnd); }
  //@}

  /// \name arithmetic functions -- addition (auto version)
  //@{
  /// addition for <tt>BigFloat+BigFloat</tt> (auto version)
  int add(const BigFloat& x, const BigFloat& y, rnd_t rnd = MPFR_RND) {
    int r = add(x, y, add_prec(x, y), rnd);  // call fixed version
    remove_trailing_zeros(); return r;       // remove extra zeros
  }
  /// addition for <tt>BigFloat+T</tt> (auto version)
  template <typename T> 
  int add(const BigFloat& x, const T& y, rnd_t rnd = MPFR_RND) { 
    int r = add(x, y, add_prec(x, count_prec(y)), rnd);  // call fixed version
    remove_trailing_zeros(); return r;                   // remove extra zeros
  }
  /// addition for <tt>T+BigFloat</tt> (auto version)
  template <typename T> 
  int add(const T& x, const BigFloat& y, rnd_t rnd = MPFR_RND)
  { return add(y, x, rnd); }
  //@}

  /// \name arithmetic functions -- subtraction (raw version)
  //@{
  /// subtraction for <tt>BigFloat-BigFloat</tt> (raw version)
  int r_sub(const BigFloat& x, const BigFloat& y, rnd_t rnd = MPFR_RND)
  { return mpfr_sub(mp(), x.mp(), y.mp(), rnd); }
  /// subtraction for <tt>BigFloat-int</tt> (raw version)
  int r_sub(const BigFloat& x, int y, rnd_t rnd = MPFR_RND)
  { return r_sub(x, static_cast<long>(y), rnd); }
  /// subtraction for <tt>BigFloat-unsigned int</tt> (raw version)
  int r_sub(const BigFloat& x, unsigned int y, rnd_t rnd = MPFR_RND)
  { return r_sub(x, static_cast<unsigned long>(y), rnd); }
  /// subtraction for <tt>BigFloat-long</tt> (raw version)
  int r_sub(const BigFloat& x, long y, rnd_t rnd = MPFR_RND)
  { return mpfr_sub_si(mp(), x.mp(), y, rnd); }
  /// subtraction for <tt>BigFloat-unsigned long</tt> (raw version)
  int r_sub(const BigFloat& x, unsigned long y, rnd_t rnd = MPFR_RND)
  { return mpfr_sub_ui(mp(), x.mp(), y, rnd); }
  /// subtraction for <tt>BigFloat-double</tt> (raw version)
  int r_sub(const BigFloat& x, double y, rnd_t rnd = MPFR_RND)
  { return r_sub(x, BigFloat(y), rnd); }
  /// subtraction for <tt>BigFloat-BigInt</tt> (raw version)
  int r_sub(const BigFloat& x, const BigInt& y, rnd_t rnd = MPFR_RND)
  { return mpfr_sub_z(mp(), x.mp(), y.mp(), rnd); }
  /// subtraction for <tt>BigFloat-BigRat</tt> (raw version)
  int r_sub(const BigFloat& x, const BigRat& y, rnd_t rnd = MPFR_RND)
  { return mpfr_sub_q(mp(), x.mp(), y.mp(), rnd); }
  /// subtraction for <tt>int-BigFloat</tt> (raw version)
  int r_sub(int x, const BigFloat& y, rnd_t rnd = MPFR_RND)
  { return r_sub(static_cast<long>(x), y, rnd); }
  /// subtraction for <tt>unsigned int-BigFloat</tt> (raw version)
  int r_sub(unsigned int x, const BigFloat& y, rnd_t rnd = MPFR_RND)
  { return r_sub(static_cast<unsigned long>(x), y, rnd); }
  /// subtraction for <tt>long-BigFloat</tt> (raw version)
  int r_sub(long x, const BigFloat& y, rnd_t rnd = MPFR_RND)
  { return mpfr_si_sub(mp(), x, y.mp(), rnd); }
  /// subtraction for <tt>unsigned long-BigFloat</tt> (raw version)
  int r_sub(unsigned long x, const BigFloat& y, rnd_t rnd = MPFR_RND)
  { return mpfr_ui_sub(mp(), x, y.mp(), rnd); }
  /// subtraction for <tt>double-BigFloat</tt> (raw version)
  int r_sub(double x, const BigFloat& y, rnd_t rnd = MPFR_RND)
  { return r_sub(BigFloat(x), y, rnd); }
  /// subtraction for <tt>BigInt-BigFloat</tt> (raw version)
  int r_sub(const BigInt& x, const BigFloat& y, rnd_t rnd = MPFR_RND)
  { return mpfr_z_sub(mp(), x.mp(), y.mp(), rnd); }
  /// subtraction for <tt>BigRat-BigFloat</tt> (raw version)
  int r_sub(const BigRat& x, const BigFloat& y, rnd_t rnd = MPFR_RND)
  { return mpfr_q_sub(mp(), x.mp(), y.mp(), rnd); }
  //@}
  
  /// \name arithmetic functions -- subtraction (fixed version)
  //@{
  /// subtraction for <tt>BigFloat-BigFloat</tt> (fixed version)
  int sub(const BigFloat& x, const BigFloat& y,prec_t prec,rnd_t rnd=MPFR_RND){
    if (&x == this || &y == this) { // if one of inputs are output
       if (prec > get_prec()) {
         prec_round (prec, rnd);
         return r_sub(x,y,rnd);
       } else if (prec < get_prec()) {
         r_sub(x,y,rnd);
         return prec_round (prec, rnd);
       } else {
         return r_sub(x,y,rnd);
       }
/*
       BigFloat result(0, prec);
       int r = result.r_sub(x, y, rnd);
       swap(result); return r;
*/     } else {
       set_prec(prec); 
       int r= r_sub(x, y, rnd);
       return r;
     }
  }
  /// subtraction for <tt>BigFloat-T</tt> (fixed version)
  template <typename T> 
  int sub(const BigFloat& x, const T& y, prec_t prec, rnd_t rnd = MPFR_RND) {
    if (&x == this) { // if x is same as output
       if (prec > get_prec()) {
         prec_round (prec, rnd);
         return r_sub(x,y,rnd);
       } else if (prec < get_prec()) {
         r_sub(x,y,rnd);
         return prec_round (prec, rnd);
       } else {
         return r_sub(x,y,rnd);
       }
/*
       BigFloat result(0, prec);
       int r = result.r_sub(x, y, rnd); 
       swap(result); return r;
*/     } else {
       set_prec(prec); return r_sub(x, y, rnd);  
     }
  }
  /// subtraction for <tt>T-BigFloat</tt> (fixed version)
  template <typename T> 
  int sub(const T& x, const BigFloat& y, prec_t prec, rnd_t rnd = MPFR_RND) {
    if (&y == this) { // if y is same as output
       if (prec > get_prec()) {
         prec_round (prec, rnd);
         return r_sub(x,y,rnd);
       } else if (prec < get_prec()) {
         r_sub(x,y,rnd);
         return prec_round (prec, rnd);
       } else {
         return r_sub(x,y,rnd);
       }
/*
       BigFloat result(0, prec);
       int r = result.r_sub(x, y, rnd); 
       swap(result); return r;
*/     } else {
       set_prec(prec); return r_sub(x, y, rnd);  
     }
  }
  //@}

  /// \name arithmetic functions -- subtraction (auto version)
  //@{
  /// subtraction for <tt>BigFloat-BigFloat</tt> (auto version)
  int sub(const BigFloat& x, const BigFloat& y, rnd_t rnd = MPFR_RND) {
    int r = sub(x, y, add_prec(x, y), rnd);  // call fixed version
    remove_trailing_zeros(); return r;       // remove extra zeros
  } 
  /// subtraction for <tt>BigFloat-T</tt> (auto version)
  template <typename T>  
  int sub(const BigFloat& x, const T& y, rnd_t rnd = MPFR_RND) {
    int r = sub(x, y, add_prec(x, count_prec(y)), rnd);  // call fixed version
    remove_trailing_zeros(); return r;                   // remove extra zeros
  } 
  /// subtraction for <tt>T-BigFloat</tt> (auto version)
  template <typename T> 
  int sub(const T& x, const BigFloat& y, rnd_t rnd = MPFR_RND) {
    int r = sub(x, y, add_prec(y, count_prec(x)), rnd);  // call fixed version
    remove_trailing_zeros(); return r;                   // remove extra zeros
  } 
  //@}

  /// \name arithmetic functions -- multiplication (raw version)
  //@{
  /// multiplication for <tt>BigFloat*BigFloat</tt>
  int r_mul(const BigFloat& x, const BigFloat& y, rnd_t rnd = MPFR_RND)
  { return mpfr_mul(mp(), x.mp(), y.mp(), rnd); }
  /// multiplication for <tt>BigFloat*int</tt>
  int r_mul(const BigFloat& x, int y, rnd_t rnd = MPFR_RND)
  { return r_mul(x, static_cast<long>(y), rnd); }
  /// multiplication for <tt>BigFloat*unsigned int</tt>
  int r_mul(const BigFloat& x, unsigned int y, rnd_t rnd = MPFR_RND)
  { return r_mul(x, static_cast<unsigned long>(y), rnd); }
  /// multiplication for <tt>BigFloat*long</tt>
  int r_mul(const BigFloat& x, long y, rnd_t rnd = MPFR_RND)
  { return mpfr_mul_si(mp(), x.mp(), y, rnd); }
  /// multiplication for <tt>BigFloat*unsigned long</tt>
  int r_mul(const BigFloat& x, unsigned long y, rnd_t rnd = MPFR_RND)
  { return mpfr_mul_ui(mp(), x.mp(), y, rnd); }
  /// multiplication for <tt>BigFloat*double</tt>
  int r_mul(const BigFloat& x, double y, rnd_t rnd = MPFR_RND)
  { return r_mul(x, BigFloat(y), rnd); }
  /// multiplication for <tt>BigFloat*BigInt</tt>
  int r_mul(const BigFloat& x, const BigInt& y, rnd_t rnd = MPFR_RND)
  { return mpfr_mul_z(mp(), x.mp(), y.mp(), rnd); }
  /// multiplication for <tt>BigFloat*BigRat</tt>
  int r_mul(const BigFloat& x, const BigRat& y, rnd_t rnd = MPFR_RND)
  { return mpfr_mul_q(mp(), x.mp(), y.mp(), rnd); }
  template <typename T> 
  int r_mul(const T& x, const BigFloat& y, rnd_t rnd = MPFR_RND)
  { return r_mul(y, x, rnd); }
  //@}

  /// \name arithmetic functions -- multiplication (fixed version)
  //@{
  /// multiplication for <tt>BigFloat*BigFloat</tt> (fixed version)
  int mul(const BigFloat& x, const BigFloat& y,prec_t prec,rnd_t rnd=MPFR_RND){
    if (&x == this || &y == this) { // if one of inputs are output
//*
       if (prec > get_prec()) {
         prec_round (prec, rnd);
         return r_mul(x,y,rnd);
       } else if (prec < get_prec()) {
         r_mul(x,y,rnd);
         return prec_round (prec, rnd);
       } else {
         return r_mul(x,y,rnd);
       }
/*/    BigFloat result(0, prec);
       int r = result.r_mul(x, y, rnd);
       swap(result); return r;
//*/
     } else {
       set_prec(prec); return r_mul(x, y, rnd);
     }
  }
  /// multiplication for <tt>BigFloat*T</tt> (fixed version)
  template <typename T>
  int mul(const BigFloat& x, const T& y, prec_t prec, rnd_t rnd = MPFR_RND) {
    if (&x == this) { // if x is same as output
//*
      if (prec > get_prec()) {
        prec_round (prec, rnd);
        return r_mul(x,y,rnd);
      } else if (prec < get_prec()) {
        r_mul(x,y,rnd);
        return prec_round (prec, rnd);
      } else {
        return r_mul(x,y,rnd);
      }
/*/   BigFloat result(0, prec);
      int r = result.r_mul(x, y, rnd);
      swap(result); return r;
//*/
     } else {
       set_prec(prec); return r_mul(x, y, rnd);
     }
  }
  /// multiplication for <tt>T*BigFloat</tt> (fixed version)
  template <typename T>
  int mul(const T& x, const BigFloat& y, prec_t prec, rnd_t rnd = MPFR_RND)
  { return mul(y, x, prec, rnd); }
  //@}

  /// \name arithmetic functions -- multiplication (auto version)
  //@{
  /// multiplication for <tt>BigFloat*BigFloat</tt> (auto version)

  int mul(const BigFloat& x, const BigFloat& y, rnd_t rnd = MPFR_RND) {
    int r = mul(x, y, mul_prec(x, y), rnd);  // call fixed version
    remove_trailing_zeros(); return r;       // remove extra zeros
  }



  /// multiplication for <tt>BigFloat*T</tt> (auto version)
  template <typename T>
  int mul(const BigFloat& x, const T& y, rnd_t rnd = MPFR_RND) {
    int r = mul(x, y, mul_prec(x, count_prec(y)), rnd);  // call fixed version
    remove_trailing_zeros(); return r;                   // remove extra zeros
  }
  /// multiplication for <tt>T*BigFloat</tt> (auto version)
  template <typename T>
  int mul(const T& x, const BigFloat& y, rnd_t rnd = MPFR_RND)
  { return mul(y, x, rnd); }
  //@}


  /// \name arithmetic functions -- division (raw version)
  //@{
  /// division for <tt>BigFloat/BigFloat</tt>
  int r_div(const BigFloat& x, const BigFloat& y, rnd_t rnd = MPFR_RND)
  { return mpfr_div(mp(), x.mp(), y.mp(), rnd); }
  /// division for <tt>BigFloat/int</tt>
  int r_div(const BigFloat& x, int y, rnd_t rnd = MPFR_RND)
  { return r_div(x, static_cast<long>(y), rnd); }
  /// division for <tt>BigFloat/unsigned int</tt>
  int r_div(const BigFloat& x, unsigned int y, rnd_t rnd = MPFR_RND)
  { return r_div(x, static_cast<unsigned long>(y), rnd); }
  /// division for <tt>BigFloat/long</tt>
  int r_div(const BigFloat& x, long y, rnd_t rnd = MPFR_RND)
  { return mpfr_div_si(mp(), x.mp(), y, rnd); }
  /// division for <tt>BigFloat/unsigned long</tt>
  int r_div(const BigFloat& x, unsigned long y, rnd_t rnd = MPFR_RND)
  { return mpfr_div_ui(mp(), x.mp(), y, rnd); }
  /// division for <tt>BigFloat/double</tt>
  int r_div(const BigFloat& x, double y, rnd_t rnd = MPFR_RND)
  { return r_div(x, BigFloat(y), rnd); }
  /// division for <tt>BigFloat/BigInt</tt>
  int r_div(const BigFloat& x, const BigInt& y, rnd_t rnd = MPFR_RND)
  { return mpfr_div_z(mp(), x.mp(), y.mp(), rnd); }
  /// division for <tt>BigFloat/BigRat</tt>
  int r_div(const BigFloat& x, const BigRat& y, rnd_t rnd = MPFR_RND)
  { return mpfr_div_q(mp(), x.mp(), y.mp(), rnd); }
  /// division for <tt>int/BigFloat</tt>
  int r_div(int x, const BigFloat& y, rnd_t rnd = MPFR_RND)
  { return r_div(static_cast<long>(x), y, rnd); }
  /// division for <tt>unsigned int/BigFloat</tt>
  int r_div(unsigned int x, const BigFloat& y, rnd_t rnd = MPFR_RND)
  { return r_div(static_cast<unsigned long>(x), y, rnd); }
  /// division for <tt>long/BigFloat</tt>
  int r_div(long x, const BigFloat& y, rnd_t rnd = MPFR_RND)
  { return mpfr_si_div(mp(), x, y.mp(), rnd); }
  /// division for <tt>unsigned long/BigFloat</tt>
  int r_div(unsigned long x, const BigFloat& y, rnd_t rnd = MPFR_RND)
  { return mpfr_ui_div(mp(), x, y.mp(), rnd); }
  /// division for <tt>double/BigFloat</tt>
  int r_div(double x, const BigFloat& y, rnd_t rnd = MPFR_RND)
  { return r_div(BigFloat(x), y, rnd); }
  /// division for <tt>BigInt/BigFloat</tt>
  int r_div(const BigInt& x, const BigFloat& y, rnd_t rnd = MPFR_RND)
  { return r_div(BigFloat(x), y, rnd); }
  /// division for <tt>BigRat/BigFloat</tt> (BigRat will be converted)
  int r_div(const BigRat& x, const BigFloat& y, rnd_t rnd = MPFR_RND)
  { return r_div(BigFloat(x), y, rnd); }
  //@}

// TODO: We need "exact division" as in BigInt::exactdiv(..)
// in order to implement exact polynomial discriminant with NT=BigFloat
// -- Jul 2010 (Chee)
//
  /// \name arithmetic functions -- division (fixed version)
  //@{
  /// division for <tt>BigFloat/BigFloat</tt> (fixed version)
  int div(const BigFloat& x, const BigFloat& y, 
          prec_t prec = getDefaultBFdivPrec(), rnd_t rnd = MPFR_RND) {
    if (&x == this || &y == this) { // if one of inputs are output
       if (prec > get_prec()) {
         prec_round (prec, rnd);
         return r_div(x,y,rnd);
       } else if (prec < get_prec()) {
         r_div(x,y,rnd);
         return prec_round (prec, rnd);
       } else {
         return r_div(x,y,rnd);
       }
/*       BigFloat result(0, prec);
       int r = result.r_div(x, y, rnd);
       swap(result); return r; 
*/     } else {
       set_prec(prec); return r_div(x, y, rnd);
     } 
  }  
  /// division for <tt>BigFloat/T</tt> (fixed version)
  template <typename T> 
  int div(const BigFloat& x, const T& y, 
          prec_t prec = getDefaultBFdivPrec(), rnd_t rnd = MPFR_RND) {
    if (&x == this) { // if x is same as output 
       if (prec > get_prec()) {
         prec_round (prec, rnd);
         return r_div(x,y,rnd);
       } else if (prec < get_prec()) {
         r_div(x,y,rnd);
         return prec_round (prec, rnd);
       } else {
         return r_div(x,y,rnd);
       }
/*       BigFloat result(0, prec);
       int r = result.r_div(x, y, rnd);
       swap(result); return r; 
*/     } else {
       set_prec(prec); return r_div(x, y, rnd);
     } 
  }  
  /// division for <tt>T/BigFloat</tt> (fixed version)
  template <typename T>
  int div(const T& x, const BigFloat& y, 
          prec_t prec = getDefaultBFdivPrec(), rnd_t rnd = MPFR_RND) {
    if (&y == this) { // if y is same as output
       if (prec > get_prec()) {
         prec_round (prec, rnd);
         return r_div(x,y,rnd);
       } else if (prec < get_prec()) {
         r_div(x,y,rnd);
         return prec_round (prec, rnd);
       } else {
         return r_div(x,y,rnd);
       }
/*       BigFloat result(0, prec);
       int r = result.r_div(x, y, rnd);
       swap(result); return r;
*/     } else {
       set_prec(prec); return r_div(x, y, rnd);
     }
  }
   /// division for <tt>T/BigFloat</tt> (fixed version)
  template <typename T>
  int div(const T& x, const T& y, 
          prec_t prec = getDefaultBFdivPrec(), rnd_t rnd = MPFR_RND) {
    return div(BigFloat(x), BigFloat(y));
  }
  //@}

  /// \name Sine functions (raw version)
  //@{
  /// Sine function for <tt>BigFloat</tt> (raw version)
  int r_sin(const BigFloat& x, rnd_t rnd = MPFR_RND)
  { return mpfr_sin(mp(), x.mp(), rnd); }
  /// Sine function for <tt>T</tt> (raw version)
  template<typename T>
  int r_sin(T x, rnd_t rnd = MPFR_RND)
  { return r_sin(BigFloat(x), rnd); }

  /// Cosine function for <tt>BigFloat</tt> (raw version)
  int r_cos(const BigFloat& x, rnd_t rnd = MPFR_RND)
  { return mpfr_cos(mp(), x.mp(), rnd); }
  /// Cosine function for <tt>T</tt> (raw version)
  template<typename T>
  int r_cos(T x, rnd_t rnd = MPFR_RND)
  { return r_cos(BigFloat(x), rnd); }

  /// Tangent function for <tt>BigFloat</tt> (raw version)
  int r_tan(const BigFloat& x, rnd_t rnd = MPFR_RND)
  { return mpfr_tan(mp(), x.mp(), rnd); }
  /// Tangent function for <tt>T</tt> (raw version)
  template<typename T>
  int r_tan(T x, rnd_t rnd = MPFR_RND)
  { return r_tan(BigFloat(x), rnd); }

  /// Cotangent function for <tt>BigFloat</tt> (raw version)
   int r_cot(const BigFloat& x, rnd_t rnd = MPFR_RND)
  { return mpfr_cot(mp(), x.mp(), rnd); }
  /// Cotangent function for <tt>T</tt> (raw version)
  template<typename T>
  int r_cot(T x, rnd_t rnd = MPFR_RND)
  { return r_cot(BigFloat(x), rnd); }

  /// Arcsine function for <tt>BigFloat</tt> (raw version)
   int r_asin(const BigFloat& x, rnd_t rnd = MPFR_RND)
  { return mpfr_asin(mp(), x.mp(), rnd); }
  /// Arcsine function for <tt>T</tt> (raw version)
  template<typename T>
  int r_asin(T x, rnd_t rnd = MPFR_RND)
  { return r_asin(BigFloat(x), rnd); }

  /// Arccosine function for <tt>BigFloat</tt> (raw version)
   int r_acos(const BigFloat& x, rnd_t rnd = MPFR_RND)
  { return mpfr_acos(mp(), x.mp(), rnd); }
  /// Arccosine function for <tt>T</tt> (raw version)
  template<typename T>
  int r_acos(T x, rnd_t rnd = MPFR_RND)
  { return r_acos(BigFloat(x), rnd); }

  /// Arctangent function for <tt>BigFloat</tt> (raw version)
   int r_atan(const BigFloat& x, rnd_t rnd = MPFR_RND)
  { return mpfr_atan(mp(), x.mp(), rnd); }
  /// Arctangent function for <tt>T</tt> (raw version)
  template<typename T>
  int r_atan(T x, rnd_t rnd = MPFR_RND)
  { return r_atan(BigFloat(x), rnd); }

  /// log function for <tt>BigFloat</tt> (raw version)
   int r_log(const BigFloat& x, rnd_t rnd = MPFR_RND)
  { return mpfr_log(mp(), x.mp(), rnd); }
  /// log function for <tt>T</tt> (raw version)
  template<typename T>
  int r_log(T x, rnd_t rnd = MPFR_RND)
  { return r_log(BigFloat(x), rnd); }  
  
  /// log2 function for <tt>BigFloat</tt> (raw version)
   int r_log2(const BigFloat& x, rnd_t rnd = MPFR_RND)
  { return mpfr_log2(mp(), x.mp(), rnd); }
  /// log_2 function for <tt>T</tt> (raw version)
  template<typename T>
  int r_log2(T x, rnd_t rnd = MPFR_RND)
  { return r_log2(BigFloat(x), rnd); }

  /// log10 function for <tt>BigFloat</tt> (raw version)
   int r_log10(const BigFloat& x, rnd_t rnd = MPFR_RND)
  { return mpfr_log10(mp(), x.mp(), rnd); }
  /// log_10 function for <tt>T</tt> (raw version)
  template<typename T>
  int r_log10(T x, rnd_t rnd = MPFR_RND)
  { return r_log10(BigFloat(x), rnd); }

  /// exponent function for <tt>BigFloat</tt> (raw version)
   int r_exp(const BigFloat& x, rnd_t rnd = MPFR_RND)
  { return mpfr_exp(mp(), x.mp(), rnd); }
  /// exponent function for <tt>T</tt> (raw version)
  template<typename T>
  int r_exp(T x, rnd_t rnd = MPFR_RND)
  { return r_exp(BigFloat(x), rnd); }

  /// 2 power function for <tt>BigFloat</tt> (raw version)
   int r_exp2(const BigFloat& x, rnd_t rnd = MPFR_RND)
  { return mpfr_exp2(mp(), x.mp(), rnd); }
  /// exponent function for <tt>T</tt> (raw version)
  template<typename T>
  int r_exp2(T x, rnd_t rnd = MPFR_RND)
  { return r_exp2(BigFloat(x), rnd); }

  /// 10 power function for <tt>BigFloat</tt> (raw version)
   int r_exp10(const BigFloat& x, rnd_t rnd = MPFR_RND)
  { return mpfr_exp10(mp(), x.mp(), rnd); }
  /// exponent function for <tt>T</tt> (raw version)
  template<typename T>
  int r_exp10(T x, rnd_t rnd = MPFR_RND)
  { return r_exp10(BigFloat(x), rnd); }

  //@}

  /// \name Sine functions (fixed version)
  //@
  /// Sine function for <tt>BigFloat</tt> (fixed version)
  int sin(const BigFloat& x, 
           prec_t prec, rnd_t rnd = MPFR_RND) {
    assert(prec>=2);
    if (&x == this) { // if x is same as ouput
       if (prec > get_prec()) {
         prec_round (prec, rnd);
         return r_sin(x,rnd);
       } else if (prec < get_prec()) {
         r_sin(x,rnd);
         return prec_round (prec, rnd);
       } else {
         return r_sin(x,rnd);
       }
    } else {
      set_prec(prec); return r_sin(x, rnd);
    }
  }
  /// Sine function for <tt>T</tt> (fixed version)
  template <typename T> 
  int sin(const T& x, prec_t prec, rnd_t rnd = MPFR_RND)
  { set_prec(prec); return r_sin(x, rnd); }
  //@}

  /// \name Cosine functions (fixed version)
  //@
  /// Cosine function for <tt>BigFloat</tt> (fixed version)
  int cos(const BigFloat& x, 
           prec_t prec, rnd_t rnd = MPFR_RND) {
    assert(prec>=2);
    if (&x == this) { // if x is same as ouput
       if (prec > get_prec()) {
         prec_round (prec, rnd);
         return r_cos(x,rnd);
       } else if (prec < get_prec()) {
         r_cos(x,rnd);
         return prec_round (prec, rnd);
       } else {
         return r_cos(x,rnd);
       }
    } else {
      set_prec(prec); return r_cos(x, rnd);
    }
  }
  /// Cosine function for <tt>T</tt> (fixed version)
  template <typename T> 
  int cos(const T& x, prec_t prec, rnd_t rnd = MPFR_RND)
  { set_prec(prec); return r_cos(x, rnd); }
  //@}

  /// \name Tangent functions (fixed version)
  //@
  /// Tangent function for <tt>BigFloat</tt> (fixed version)
  int tan(const BigFloat& x, 
           prec_t prec, rnd_t rnd = MPFR_RND) {
    assert(prec>=2);
    if (&x == this) { // if x is same as ouput
       if (prec > get_prec()) {
         prec_round (prec, rnd);
         return r_tan(x,rnd);
       } else if (prec < get_prec()) {
         r_tan(x,rnd);
         return prec_round (prec, rnd);
       } else {
         return r_tan(x,rnd);
       }
    } else {
      set_prec(prec); return r_tan(x, rnd);
    }
  }
  /// Tangent function for <tt>T</tt> (fixed version)
  template <typename T> 
  int tan(const T& x, prec_t prec, rnd_t rnd = MPFR_RND)
  { set_prec(prec); return r_tan(x, rnd); }
  //@}

  /// \name Cotangent functions (fixed version)
  //@
  /// Cotangent function for <tt>BigFloat</tt> (fixed version)
  int cot(const BigFloat& x, 
           prec_t prec, rnd_t rnd = MPFR_RND) {
    assert(prec>=2);
    if (&x == this) { // if x is same as ouput
       if (prec > get_prec()) {
         prec_round (prec, rnd);
         return r_cot(x,rnd);
       } else if (prec < get_prec()) {
         r_cot(x,rnd);
         return prec_round (prec, rnd);
       } else {
         return r_cot(x,rnd);
       }
    } else {
      set_prec(prec); return r_cot(x, rnd);
    }
  }
  /// Cotangent function for <tt>T</tt> (fixed version)
  template <typename T> 
  int cot(const T& x, prec_t prec, rnd_t rnd = MPFR_RND)
  { set_prec(prec); return r_cot(x, rnd); }
  //@}

  /// \name Arcsine functions (fixed version)
  //@
  /// Arcsine function for <tt>BigFloat</tt> (fixed version)
  int asin(const BigFloat& x, 
           prec_t prec, rnd_t rnd = MPFR_RND) {
    assert(prec>=2);
    if (&x == this) { // if x is same as ouput
       if (prec > get_prec()) {
         prec_round (prec, rnd);
         return r_asin(x,rnd);
       } else if (prec < get_prec()) {
         r_asin(x,rnd);
         return prec_round (prec, rnd);
       } else {
         return r_asin(x,rnd);
       }
    } else {
      set_prec(prec); return r_asin(x, rnd);
    }
  }
  /// Arcsine function for <tt>T</tt> (fixed version)
  template <typename T> 
  int asin(const T& x, prec_t prec, rnd_t rnd = MPFR_RND)
  { set_prec(prec); return r_asin(x, rnd); }
  //@}

  /// \name Arccosine functions (fixed version)
  //@
  /// Arccosine function for <tt>BigFloat</tt> (fixed version)
  int acos(const BigFloat& x, 
           prec_t prec, rnd_t rnd = MPFR_RND) {
    assert(prec>=2);
    if (&x == this) { // if x is same as ouput
       if (prec > get_prec()) {
         prec_round (prec, rnd);
         return r_acos(x,rnd);
       } else if (prec < get_prec()) {
         r_acos(x,rnd);
         return prec_round (prec, rnd);
       } else {
         return r_acos(x,rnd);
       }
    } else {
      set_prec(prec); return r_acos(x, rnd);
    }
  }
  /// Arccosine function for <tt>T</tt> (fixed version)
  template <typename T> 
  int acos(const T& x, prec_t prec, rnd_t rnd = MPFR_RND)
  { set_prec(prec); return r_acos(x, rnd); }
  //@}

  /// \name Arctangent functions (fixed version)
  //@
  /// Arctangent function for <tt>BigFloat</tt> (fixed version)
  int atan(const BigFloat& x, 
           prec_t prec, rnd_t rnd = MPFR_RND) {
    assert(prec>=2);
    if (&x == this) { // if x is same as ouput
       if (prec > get_prec()) {
         prec_round (prec, rnd);
         return r_atan(x,rnd);
       } else if (prec < get_prec()) {
         r_atan(x,rnd);
         return prec_round (prec, rnd);
       } else {
         return r_atan(x,rnd);
       }
    } else {
      set_prec(prec); return r_atan(x, rnd);
    }
  }
  /// Arctangent function for <tt>T</tt> (fixed version)
  template <typename T> 
  int atan(const T& x, prec_t prec, rnd_t rnd = MPFR_RND)
  { set_prec(prec); return r_atan(x, rnd); }
  //@}
  /// \name log2 functions (fixed version)
  //@
  /// log2 function for <tt>BigFloat</tt> (fixed version)
  int log2(const BigFloat& x, 
           prec_t prec, rnd_t rnd = MPFR_RND) {
    assert(prec>=2);
    if (&x == this) { // if x is same as ouput
       if (prec > get_prec()) {
         prec_round (prec, rnd);
         return r_log2(x,rnd);
       } else if (prec < get_prec()) {
         r_log2(x,rnd);
         return prec_round (prec, rnd);
       } else {
         return r_log2(x,rnd);
       }
    } else {
      set_prec(prec); return r_log2(x, rnd);
    }
  }
  /// log2 function for <tt>T</tt> (fixed version)
  template <typename T> 
  int log2(const T& x, prec_t prec, rnd_t rnd = MPFR_RND)
  { set_prec(prec); return r_log2(x, rnd); }
  //@}
  /// \name log functions (fixed version)
  //@
  /// log function for <tt>BigFloat</tt> (fixed version)
  int log(const BigFloat& x, 
           prec_t prec, rnd_t rnd = MPFR_RND) {
    assert(prec>=2);
    if (&x == this) { // if x is same as ouput
       if (prec > get_prec()) {
         prec_round (prec, rnd);
         return r_log(x,rnd);
       } else if (prec < get_prec()) {
         r_log(x,rnd);
         return prec_round (prec, rnd);
       } else {
         return r_log(x,rnd);
       }
    } else {
      set_prec(prec); return r_log(x, rnd);
    }
  }
  /// log function for <tt>T</tt> (fixed version)
  template <typename T> 
  int log(const T& x, prec_t prec, rnd_t rnd = MPFR_RND)
  { set_prec(prec); return r_log(x, rnd); }
  //@}
  //
  /// \name log10 functions (fixed version)
  //@
  /// log10 function for <tt>BigFloat</tt> (fixed version)
  int log10(const BigFloat& x, 
           prec_t prec, rnd_t rnd = MPFR_RND) {
    assert(prec>=2);
    if (&x == this) { // if x is same as ouput
       if (prec > get_prec()) {
         prec_round (prec, rnd);
         return r_log10(x,rnd);
       } else if (prec < get_prec()) {
         r_log10(x,rnd);
         return prec_round (prec, rnd);
       } else {
         return r_log10(x,rnd);
       }
    } else {
      set_prec(prec); return r_log10(x, rnd);
    }
  }
  /// log10 function for <tt>T</tt> (fixed version)
  template <typename T> 
  int log10(const T& x, prec_t prec, rnd_t rnd = MPFR_RND)
  { set_prec(prec); return r_log10(x, rnd); }
  //@}

  /// \name exponent functions (fixed version)
  //@
  /// exponent function for <tt>BigFloat</tt> (fixed version)
  int exp(const BigFloat& x, 
           prec_t prec, rnd_t rnd = MPFR_RND) {
    assert(prec>=2);
    if (&x == this) { // if x is same as ouput
       if (prec > get_prec()) {
         prec_round (prec, rnd);
         return r_exp(x,rnd);
       } else if (prec < get_prec()) {
         r_exp(x,rnd);
         return prec_round (prec, rnd);
       } else {
         return r_exp(x,rnd);
       }
    } else {
      set_prec(prec); return r_exp(x, rnd);
    }
  }
  /// exponent function for <tt>T</tt> (fixed version)
  template <typename T> 
  int exp(const T& x, prec_t prec, rnd_t rnd = MPFR_RND)
  { set_prec(prec); return r_exp(x, rnd); }
  //@}

  /// \name 2 power functions (fixed version)
  //@
  /// 2 power function for <tt>BigFloat</tt> (fixed version)
  int exp2(const BigFloat& x, 
           prec_t prec, rnd_t rnd = MPFR_RND) {
    assert(prec>=2);
    if (&x == this) { // if x is same as ouput
       if (prec > get_prec()) {
         prec_round (prec, rnd);
         return r_exp2(x,rnd);
       } else if (prec < get_prec()) {
         r_exp2(x,rnd);
         return prec_round (prec, rnd);
       } else {
         return r_exp2(x,rnd);
       }
    } else {
      set_prec(prec); return r_exp2(x, rnd);
    }
  }
  /// 2 power function for <tt>T</tt> (fixed version)
  template <typename T> 
  int exp2(const T& x, prec_t prec, rnd_t rnd = MPFR_RND)
  { set_prec(prec); return r_exp2(x, rnd); }
  //@}

  /// \name 10 power functions (fixed version)
  //@
  /// exponent function for <tt>BigFloat</tt> (fixed version)
  int exp10(const BigFloat& x, 
           prec_t prec, rnd_t rnd = MPFR_RND) {
    assert(prec>=2);
    if (&x == this) { // if x is same as ouput
       if (prec > get_prec()) {
         prec_round (prec, rnd);
         return r_exp10(x,rnd);
       } else if (prec < get_prec()) {
         r_exp10(x,rnd);
         return prec_round (prec, rnd);
       } else {
         return r_exp10(x,rnd);
       }
    } else {
      set_prec(prec); return r_exp10(x, rnd);
    }
  }
  /// 10 power function for <tt>T</tt> (fixed version)
  template <typename T> 
  int exp10(const T& x, prec_t prec, rnd_t rnd = MPFR_RND)
  { set_prec(prec); return r_exp10(x, rnd); }
  //@}
  
  
  
  /// \name square root functions (raw version)
  //@{
  /// square root function for <tt>BigFloat</tt> (raw version)
  int r_sqrt(const BigFloat& x, rnd_t rnd = MPFR_RND)
  { return mpfr_sqrt(mp(), x.mp(), rnd); }
  /// square root function for <tt>int</tt> (raw version)
  int r_sqrt(int x, rnd_t rnd = MPFR_RND)
  { return r_sqrt(static_cast<long>(x), rnd); }
  /// square root function for <tt>unsigned int</tt> (raw version)
  int r_sqrt(unsigned int x, rnd_t rnd = MPFR_RND)
  { return r_sqrt(static_cast<unsigned long>(x), rnd); }
  /// square root function for <tt>long</tt> (raw version)
  int r_sqrt(long x, rnd_t rnd = MPFR_RND)
  { assert(x>=0); return mpfr_sqrt_ui(mp(), x, rnd); }
  /// square root function for <tt>unsigned long</tt> (raw version)
  int r_sqrt(unsigned long x, rnd_t rnd = MPFR_RND)
  { return mpfr_sqrt_ui(mp(), x, rnd); }
  /// square root function for <tt>double</tt> (raw version)
  int r_sqrt(double x, rnd_t rnd = MPFR_RND)
  { return r_sqrt(BigFloat(x), rnd); }
  //@}

  /// \name square root functions (fixed version)
  //@{
  /// square root function for <tt>BigFloat</tt> (fixed version)
  //  Returns 0 if the sqrt is an exact operation,
  //  Else it is a rounded operation
  int sqrt(const BigFloat& x, 
           prec_t prec = getDefaultBFradicalPrec(), rnd_t rnd = MPFR_RND) {
    assert(prec>=2);
    if (&x == this) { // if x is same as ouput
       if (prec > get_prec()) {
         prec_round (prec, rnd);
         return r_sqrt(x,rnd);
       } else if (prec < get_prec()) {
         r_sqrt(x,rnd);
         return prec_round (prec, rnd);
       } else {
         return r_sqrt(x,rnd);
       }
/*       BigFloat result(0, prec);
      int r = result.r_sqrt(x, rnd);
      swap(result); return r;
*/    } else {
      set_prec(prec); return r_sqrt(x, rnd);
    }
  }
  /// square root function for <tt>T</tt> (fixed version)
  template <typename T> 
  int sqrt(const T& x, prec_t prec = getDefaultBFradicalPrec(), rnd_t rnd = MPFR_RND)
  { set_prec(prec); return r_sqrt(x, rnd); }
  //@}

  // REMARK: unlike Core 1, there is no version of sqrt where we already
  // know a very good approximation, and we only want to increase
  // its precision.  May be we should improve this.

  /// \name power functions (raw version)
  //@{
  /// power function for <tt>BigFloat^BigFloat</tt> (raw version)
  int r_pow(const BigFloat& x, const BigFloat& y, rnd_t rnd = MPFR_RND)
  { return mpfr_pow(mp(), x.mp(), y.mp(), rnd); }
  /// power function for <tt>BigFloat^int</tt> (raw version)
  int r_pow(const BigFloat& x, int y, rnd_t rnd = MPFR_RND)
  { return r_pow(x, static_cast<long>(y), rnd); }
  /// power function for <tt>BigFloat^unsigned int</tt> (raw version)
  int r_pow(const BigFloat& x, unsigned int y, rnd_t rnd = MPFR_RND)
  { return r_pow(x, static_cast<unsigned long>(y), rnd); }
  /// power function for <tt>BigFloat^long</tt> (raw version)
  int r_pow(const BigFloat& x, long y, rnd_t rnd = MPFR_RND)
  { return mpfr_pow_si(mp(), x.mp(), y, rnd); }
  /// power function for <tt>BigFloat^unsigned long</tt> (raw version)
  int r_pow(const BigFloat& x, unsigned long y, rnd_t rnd = MPFR_RND)
  { return mpfr_pow_ui(mp(), x.mp(), y, rnd); }
  /// power function for <tt>BigFloat^BigInt</tt> (raw version)
  int r_pow(const BigFloat& x, const BigInt& y, rnd_t rnd = MPFR_RND)
  { return mpfr_pow_z(mp(), x.mp(), y.mp(), rnd); }
  /// power function for <tt>unsigned long^BigFloat</tt> (raw version)
  int r_pow(unsigned long x, const BigFloat& y, rnd_t rnd = MPFR_RND)
  { return mpfr_ui_pow(mp(), x, y.mp(), rnd); }
  /// power function for <tt>unsigned long^unsigned long</tt> (raw version)
  int r_pow(unsigned long x, unsigned long y, rnd_t rnd = MPFR_RND)
  { return mpfr_ui_pow_ui(mp(), x, y, rnd); }
  //@}

  /// \name power functions (fixed version)
  //@{
  /// power function for <tt>BigFloat^BigFloat</tt> (fixed version)
  int pow(const BigFloat& x, const BigFloat& y,prec_t prec,rnd_t rnd=MPFR_RND){
    assert(prec>=2);
    if (&x == this || &y == this) { // if one of inputs are output
       if (prec > get_prec()) {
         prec_round (prec, rnd);
         return r_pow(x,y,rnd);
       } else if (prec < get_prec()) {
         r_pow(x,y,rnd);
         return prec_round (prec, rnd);
       } else {
         return r_pow(x,y,rnd);
       }
/*       BigFloat result(0, prec);
       int r = result.r_pow(x, y, rnd);
       swap(result); return r;
*/     } else {
       set_prec(prec); return r_pow(x, y, rnd);
     }
  }
  /// power function for <tt>BigFloat^T</tt> (fixed version)
  template <typename T>
  int pow(const BigFloat& x, const T& y, prec_t prec, rnd_t rnd = MPFR_RND) {
    assert(prec>=2);
    if (&x == this) { // if x is same as output 
       if (prec > get_prec()) {
         prec_round (prec, rnd);
         return r_pow(x,y,rnd);
       } else if (prec < get_prec()) {
         r_pow(x,y,rnd);
         return prec_round (prec, rnd);
       } else {
         return r_pow(x,y,rnd);
       }
/*       BigFloat result(0, prec);
       int r = result.r_pow(x, y, rnd);
       swap(result); return r;
*/     } else {
       set_prec(prec); return r_pow(x, y, rnd);
     }
  }
  /// power function for <tt>T^BigFloat</tt> (fixed version)
  template <typename T>
  int pow(const T& x, const BigFloat& y, prec_t prec, rnd_t rnd = MPFR_RND) {
    assert(prec>=2);
    if (&y == this) { // if y is same as output
       if (prec > get_prec()) {
         prec_round (prec, rnd);
         return r_pow(x,y,rnd);
       } else if (prec < get_prec()) {
         r_pow(x,y,rnd);
         return prec_round (prec, rnd);
       } else {
         return r_pow(x,y,rnd);
       }
/*       BigFloat result(0, prec);
       int r = result.r_pow(x, y, rnd);
       swap(result); return r;
*/     } else {
       set_prec(prec); return r_pow(x, y, rnd);
     }
  }
  //@}

  /// \name other arithmetic functions (raw version)
  //@{
  /// square
  int r_sqr(const BigFloat& x, rnd_t rnd = MPFR_RND)
  { return mpfr_sqr(mp(), x.mp(), rnd); }
  /// cubic root
  int r_cbrt(const BigFloat& x, rnd_t rnd = MPFR_RND)
  { return mpfr_cbrt(mp(), x.mp(), rnd); }
  /// kth root
  int r_root(const BigFloat& x, unsigned long k, rnd_t rnd = MPFR_RND)
  { return mpfr_root(mp(), x.mp(), k, rnd); }
  /// negation
  int r_neg(const BigFloat& x, rnd_t rnd = MPFR_RND)
  { return mpfr_neg(mp(), x.mp(), rnd); }
  /// absolute value
  int r_abs(const BigFloat& x, rnd_t rnd = MPFR_RND)
  { return mpfr_abs(mp(), x.mp(), rnd); }
  //@}
  
  /// \name other arithmetic functions (fixed version)
  //@{
  /// square (fixed version)
  int sqr(const BigFloat& x, prec_t prec, rnd_t rnd = MPFR_RND) {
    assert(prec>=2);
    if (&x == this) { // if y is same as output
       if (prec > get_prec()) {
         prec_round (prec, rnd);
         return r_sqr(x,rnd);
       } else if (prec < get_prec()) {
         r_sqr(x,rnd);
         return prec_round (prec, rnd);
       } else {
         return r_sqr(x,rnd);
       }
/*       BigFloat result(0, prec);
       int r = result.r_sqr(x, rnd);
       swap(result); return r;
*/     } else {
       set_prec(prec); return r_sqr(x, rnd);
     }
  }
  /// cubic root (fixed version)
  int cbrt(const BigFloat& x, 
           prec_t prec = getDefaultBFradicalPrec(), rnd_t rnd = MPFR_RND) {
    assert(prec>=2);
    if (&x == this) { // if y is same as output
       if (prec > get_prec()) {
         prec_round (prec, rnd);
         return r_cbrt(x,rnd);
       } else if (prec < get_prec()) {
         r_cbrt(x,rnd);
         return prec_round (prec, rnd);
       } else {
         return r_cbrt(x,rnd);
       }
/*       BigFloat result(0, prec);
       int r = result.r_cbrt(x, rnd);
       swap(result); return r;
*/     } else {
       set_prec(prec); return r_cbrt(x, rnd);
     }
  }
  /// kth root (fixed version)
  int root(const BigFloat& x, unsigned long k,
           prec_t prec = getDefaultBFradicalPrec(), rnd_t rnd=MPFR_RND) {
    assert(prec>=2);
    if (&x == this) { // if y is same as output
       if (prec > get_prec()) {
         prec_round (prec, rnd);
         return r_root(x,rnd);
       } else if (prec < get_prec()) {
         r_root(x,rnd);
         return prec_round (prec, rnd);
       } else {
         return r_root(x,rnd);
       }
/*       BigFloat result(0, prec);
       int r = result.r_root(x, k, rnd);
       swap(result); return r;
*/    } else {
       set_prec(prec); return r_root(x, k, rnd);
     }
  }
  /// negation (fixed version)
  int neg(const BigFloat& x, prec_t prec, rnd_t rnd = MPFR_RND) {
    assert(prec>=2);
    if (&x == this) { // if y is same as output
       if (prec > get_prec()) {
         prec_round (prec, rnd);
         return r_neg(x,rnd);
       } else if (prec < get_prec()) {
         r_neg(x,rnd);
         return prec_round (prec, rnd);
       } else {
         return r_neg(x,rnd);
       }
/*       BigFloat result(0, prec);
       int r = result.neg(x, rnd);
       swap(result); return r;
*/     } else {
       set_prec(prec); return r_neg(x, rnd);
     }
  }
  /// absolute value (fixed version)
  int abs(const BigFloat& x, prec_t prec, rnd_t rnd = MPFR_RND) {
    assert(prec>=2);
    if (&x == this) { // if y is same as output
       if (prec > get_prec()) {
         prec_round (prec, rnd);
         return r_abs(x,rnd);
       } else if (prec < get_prec()) {
         r_abs(x,rnd);
         return prec_round (prec, rnd);
       } else {
         return r_abs(x,rnd);
       }
/*       BigFloat result(0, prec);
       int r = result.r_abs(x, rnd);
       swap(result); return r;
*/     } else {
       set_prec(prec); return r_abs(x, rnd);
     }
  }
  //@}

  /// \name other arithmetic functions (auto version)
  //@{
  /// negation (auto version)
  int neg(const BigFloat& x, rnd_t rnd = MPFR_RND)
  { if (&x != this) set_prec(x.get_prec()); return r_neg(x, rnd); }
  /// absolute value (auto version)
  int abs(const BigFloat& x, rnd_t rnd = MPFR_RND)
  { if (&x != this) set_prec(x.get_prec()); return r_abs(x, rnd); }
  //@}

  /// \name shift functions
  /// self modifying versions are implemented as member functions
  /// non-self modifying versions are implemented as friend functions in BigFloat.inl
  /// left shift for <tt>int</tt>
  BigFloat& mul_2exp(int y)
  { return mul_2exp(static_cast<long>(y)); }
  /// left shift for <tt>unsigned int</tt>
  BigFloat& mul_2exp(unsigned int y)
  { return mul_2exp(static_cast<unsigned long>(y)); }
  /// left shift for <tt>long</tt>
  BigFloat& mul_2exp(long y)
  { mpfr_mul_2si(mp(), mp(), y, MPFR_RND); return (*this); }
  /// left shift for <tt>unsigned long</tt>
  BigFloat& mul_2exp(unsigned long y)
  { mpfr_mul_2ui(mp(), mp(), y, MPFR_RND); return(*this); }
  /// right shift for <tt>int</tt>
  BigFloat& div_2exp(int y)
  { return div_2exp(static_cast<long>(y)); }
  /// right shift for <tt>unsigned int</tt>
  BigFloat& div_2exp(unsigned int y)
  { return div_2exp(static_cast<unsigned long>(y)); }
  /// right shift for <tt>long</tt>
  BigFloat& div_2exp(long y)
  { mpfr_div_2si(mp(), mp(), y, MPFR_RND); return (*this); }
  /// right shift for <tt>unsigned long</tt>
  BigFloat& div_2exp(unsigned long y)
  { mpfr_div_2ui(mp(), mp(), y, MPFR_RND); return (*this); }
  /// divide by 2
  BigFloat& div2()
  { return div_2exp(1); }
  //@}

  /// \name comparison functions
  //@{
  /// compare with <tt>BigFloat</tt>
  int cmp(const BigFloat& x) const
  { return mpfr_cmp(mp(), x.mp()); }
  /// compare with <tt>int</tt>
  int cmp(int x) const
  { return cmp(static_cast<long>(x)); }
  /// compare with <tt>unsigned int</tt>
  int cmp(unsigned int x) const
  { return cmp(static_cast<unsigned long>(x)); }
  /// compare with <tt>long</tt>
  int cmp(long x) const
  { return mpfr_cmp_si(mp(), x); }
  /// compare with <tt>unsigned long</tt>
  int cmp(unsigned long x) const
  { return mpfr_cmp_ui(mp(), x); }
  /// compare with <tt>double</tt>
  int cmp(double x) const
  { return mpfr_cmp_d(mp(), x); }
  /// compare with <tt>BigInt</tt>
  int cmp(const BigInt& x) const
  { return mpfr_cmp_z(mp(), x.mp()); }
  /// compare with <tt>BigRat</tt>
  int cmp(const BigRat& x) const
  { return mpfr_cmp_q(mp(), x.mp()); }
  /// compare with \f$x*2^e\f$ for <tt>int</tt>
  int cmp_2exp(int x, exp_t e) const
  { return mpfr_cmp_si_2exp(mp(), x, e); }
  /// compare with \f$x*2^e\f$ for <tt>unsigned int</tt>
  int cmp_2exp(unsigned int x, exp_t e) const
  { return mpfr_cmp_ui_2exp(mp(), x, e); }
  /// compare with \f$x*2^e\f$ for <tt>long</tt>
  int cmp_2exp(long x, exp_t e) const
  { return mpfr_cmp_si_2exp(mp(), x, e); }
  /// compare with \f$x*2^e\f$ for <tt>unsigned long</tt>
  int cmp_2exp(unsigned long x, exp_t e) const
  { return mpfr_cmp_ui_2exp(mp(), x, e); }
  /// compare (in absolute value) with <tt>BigFloat</tt>
  int cmpabs(const BigFloat& x) const
  { return mpfr_cmpabs(mp(), x.mp()); }
  //@}

  /// \name conversion functions
  //@{
  /// return double value
  double get_d(rnd_t rnd = MPFR_RND) const
  { return mpfr_get_d(mp(), rnd); }
  /// find d and exp s.t. \f$d*2^{exp}\f$ with \f$0.5\le|d|<1\f$
  double get_d_2exp(long* exp, rnd_t rnd = MPFR_RND) const
  { return mpfr_get_d_2exp(exp, mp(), rnd); }
  /// return long value
  long get_si(rnd_t rnd = MPFR_RND) const
  { return mpfr_get_si(mp(), rnd); }
  /// return unsigned long value
  unsigned long get_ui(rnd_t rnd = MPFR_RND) const
  { return mpfr_get_ui(mp(), rnd); }
  /// return BigInt value
  BigInt get_z(rnd_t rnd = MPFR_RND) const
  { BigInt r; mpfr_get_z(r.mp(), mp(), rnd); return r; }
  /// return z and exp s.t. it equals \f$x*2^{exp}\f$
  exp_t get_z_exp(BigInt& z) const
  { return mpfr_get_z_exp(z.mp(), mp()); }
  /// return BigRat value
  BigRat get_q() const {
    if (this->sgn() == 0) return 0;
    BigInt x; exp_t e = get_z_exp(x);
    if (e >= 0) { // convert to integer
      x.mul_2exp(x, e); return BigRat(x);
    } else { // convert to rational
      BigRat q; q.div_2exp(x, -e); return q;
    }
  }

  /// get_str(n, b, roundmode)
  ///    n=0 means we treat the bigFloat as EXACT and print the
  ///           complete decimal representation (this is implemented by MPFR)
  ///     n=1 is invalid for MPFR, so we change n=2 in this case.
  ///     n>1: we output a decimal string up to this length, again
  ///     depending on MPFR logic.
  std::string get_str(size_t n=0,int b=10, rnd_t rnd=MPFR_RND) const {
    if (n == 1) n = 2UL; 	// If n=0, then MPFR will return
    				// a string with the full precision.
				// However, n=1 is invalid for MPFR.
    return mpfr2str(mp(),
                    (std::min)(
                      (unsigned long)n,
	               bits2digits(get_prec()+1)),
		    b,
  	            get_output_fmt(),
		    rnd,
		    get_output_showpoint(),
                    get_output_showpos(),
		    get_output_uppercase());
  }
  //@}
  
  /// return 2^{e}
  static BigFloat exp2(int e)
  { return BigFloat(1, e, 2); }

  /// \name helper functions
  //@{
  bool is_nan() const
  { return mpfr_nan_p(mp()) != 0; }
  bool is_inf() const
  { return mpfr_inf_p(mp()) != 0; }
  bool is_number() const
  { return mpfr_number_p(mp()) != 0; }
  bool is_zero() const
  { return mpfr_zero_p(mp()) != 0; }
  bool greater(const BigFloat& x) const
  { return mpfr_greater_p(mp(), x.mp()) != 0; }
  bool greaterequal(const BigFloat& x) const
  { return mpfr_greaterequal_p(mp(), x.mp()) != 0; }
  bool less(const BigFloat& x) const
  { return mpfr_less_p(mp(), x.mp()) != 0; }
  bool lessequal(const BigFloat& x) const
  { return mpfr_lessequal_p(mp(), x.mp()) != 0; }
  bool lessgreater(const BigFloat& x) const
  { return mpfr_lessgreater_p(mp(), x.mp()) != 0; }
  bool equal(const BigFloat& x) const
  { return mpfr_equal_p(mp(), x.mp()) != 0; }
  bool unordered(const BigFloat& x) const
  { return mpfr_unordered_p(mp(), x.mp()) != 0; }
  bool is_integer() const
  { return mpfr_integer_p(mp()) != 0; }
  bool is_ulong(rnd_t rnd = MPFR_RND) const
  { return mpfr_fits_ulong_p(mp(), rnd) != 0; }
  bool is_slong(rnd_t rnd = MPFR_RND) const
  { return mpfr_fits_slong_p(mp(), rnd) != 0; }
  bool is_uint(rnd_t rnd = MPFR_RND) const
  { return mpfr_fits_uint_p(mp(), rnd) != 0; }
  bool is_sint(rnd_t rnd = MPFR_RND) const
  { return mpfr_fits_sint_p(mp(), rnd) != 0; }
  bool is_ushort(rnd_t rnd = MPFR_RND) const
  { return mpfr_fits_ushort_p(mp(), rnd) != 0; }
  bool is_sshort(rnd_t rnd = MPFR_RND) const
  { return mpfr_fits_sshort_p(mp(), rnd) != 0; }
  //@}

  /// demotion and promotion functions
  //@{
  /// demote from Expr
  
  static const prec_t UNUSED_PREC = 0;
  /// To handle calls to "approx" for Expressions.
  BigFloat approx(prec_t arg1=UNUSED_PREC, prec_t arg2=UNUSED_PREC) {
    return (*this);
  }

  //@}
  
  /// miscellaneous functions
  //@{
  /// set to +infty
  void set_pos_inf() 
  { mpfr_set_inf(mp(), 1); }
  /// set to -infty
  void set_neg_inf() 
  { mpfr_set_inf(mp(), -1); }
  /// set to NaN
  void set_nan() 
  { mpfr_set_nan(mp()); }
  /// swap Function
  void swap(BigFloat& other)
  { mpfr_swap(mp(), other.mp()); }
  /// return sign
  int sgn() const
  { return mpfr_sgn(mp()); }
  /// return upper bound of MSB
  long uMSB() const
  { if (sgn () == 0) return 0;
    BigInt x; exp_t e = get_z_exp(x); return x.ceillg() + e; }
  /// return lower bound of MSB
  long lMSB() const
  { if (sgn () == 0) return MSB_MIN;
    BigInt x; exp_t e = get_z_exp(x); return x.floorlg() + e; }
  /// remove trailing zeros
  void remove_trailing_zeros()
  { mpfr_remove_trailing_zeros(mp()); }

  /// special constant functions
  //@{
  /// pi
  void pi(prec_t prec, rnd_t rnd = MPFR_RND) {
    set_prec(prec);
    mpfr_const_pi(mp(), rnd);
  }
  /// e
  void e(prec_t prec, rnd_t rnd = MPFR_RND) {
    set_prec(prec);
    BigFloat bf(1);
    mpfr_exp(mp(), bf.mp(), rnd);
  }

  void rint(const BigFloat& x, rnd_t rnd = MPFR_RND)
  { mpfr_rint(mp(), x.mp(), rnd); }
  void ceil(const BigFloat& x)
  { mpfr_ceil(mp(), x.mp()); }
  void floor(const BigFloat& x)
  { mpfr_floor(mp(), x.mp()); }
  void round(const BigFloat& x)
  { mpfr_round(mp(), x.mp()); }
  void trunc(const BigFloat& x)
  { mpfr_trunc(mp(), x.mp()); }
  void rint_ceil(const BigFloat& x, rnd_t rnd = MPFR_RND)
  { mpfr_rint_ceil(mp(), x.mp(), rnd); }
  void rint_floor(const BigFloat& x, rnd_t rnd = MPFR_RND)
  { mpfr_rint_floor(mp(), x.mp(), rnd); }
  void rint_round(const BigFloat& x, rnd_t rnd = MPFR_RND)
  { mpfr_rint_round(mp(), x.mp(), rnd); }
  void rint_trunc(const BigFloat& x, rnd_t rnd = MPFR_RND)
  { mpfr_rint_trunc(mp(), x.mp(), rnd); }
  void frac(const BigFloat& x, rnd_t rnd = MPFR_RND)
  { mpfr_frac(mp(), x.mp(), rnd); }
  void nexttoward(const BigFloat& x) 
  { mpfr_nexttoward(mp(), x.mp()); }
  void nextabove() 
  { mpfr_nextabove(mp()); }
  void nextbelow() 
  { mpfr_nextbelow(mp()); }

#define CORE_PREVENT_MACRO_SUBSTITUTION

  void min CORE_PREVENT_MACRO_SUBSTITUTION (const BigFloat& x, const BigFloat& y, rnd_t rnd = MPFR_RND)
  { mpfr_min(mp(), x.mp(), y.mp(), rnd); }
  void max CORE_PREVENT_MACRO_SUBSTITUTION (const BigFloat& x, const BigFloat& y, rnd_t rnd = MPFR_RND)
  { mpfr_max(mp(), x.mp(), y.mp(), rnd); }
  //@}

  // count the precision of a int
  static prec_t count_prec(int)
  { return INT_PREC; }
  // count the precision of a unsigned int
  static prec_t count_prec(unsigned int)
  { return INT_PREC; }
  // count the precision of a long
  static prec_t count_prec(long)
  { return INT_PREC; }
  // count the precision of a unsigned long
  static prec_t count_prec(unsigned long)
  { return INT_PREC; }
  // count the precision of a double
  static prec_t count_prec(double)
  { return DOUBLE_PREC; }
  // count the precision of a BigInt
  static prec_t count_prec(const BigInt& z)
  { return (std::max)(z.sizeinbase(2), INT_PREC); }
  // count the precision of a BigInt
  static prec_t count_prec(const BigRat&)
  { return DOUBLE_PREC; }
  static prec_t count_prec(const BigFloat& z)
  { return z.get_prec(); }
  // count the precision in a string representation
  //   prec <= len*ilogb(base) <= len*(1+ilogb(base))
  static prec_t count_prec(const char* str, int base = 10)
  { return strlen(str)*(1+ilogb(base)); }	
  static prec_t count_prec(const std::string& str, int base = 10)
  { return str.length()*(1+ilogb(base)); }	

  // count how many precision needed for exact addition/subtraction
  // between one bigfloat and one integer
  static prec_t add_prec(const BigFloat& x, prec_t prec) {

    if (get_wbf_mode() == true)
      return get_wbf_prec();

    exp_t diff = x.get_exp()- x.get_prec();
    if (diff >= 0)
      return (std::max)((prec_t)2, 1+(std::max)(x.get_prec() + diff, prec));
    else
      return (std::max)((prec_t)2, 1+(std::max)(x.get_prec(), prec - diff));

//    if (get_wbf_mode() == true)
//      return (std::min)(get_wbf_prec(), ret);
//    else
//      return (std::max)(2UL, 1+(std::max)(x.get_prec(), prec - diff));

  }
  // count how many precision needed for exact addition/subtraction
  // of two bigfloats, see lemma in Zilin's thesis
  static prec_t add_prec(const BigFloat& x, const BigFloat& y) {

    if (get_wbf_mode() == true)
      return get_wbf_prec();

    int diff = (int)x.get_exp()-(int)x.get_prec()
                 -(int)y.get_exp()+(int)y.get_prec();
    if (diff >= 0)
      return (std::max)((prec_t)2, 1+(std::max)(x.get_prec() + diff, y.get_prec()));
    else
      return (std::max)((prec_t)2, 1+(std::max)(x.get_prec(), y.get_prec() - diff));

//    if (get_wbf_mode() == true)
//      return (std::min)(get_wbf_prec(), ret);
//    else
//      return ret;

  }
  // count how many precision needed for exact multiplication
  // between one bigfloat and one integer
  static prec_t mul_prec(const BigFloat& x, prec_t prec) { 
    prec_t ret;	  
    ret = x.get_prec() + prec;
    if (get_wbf_mode() == true)
      return (std::min)((prec_t)get_wbf_prec(), ret);
    else
      return ret;
  }
  // count how many precision needed for exact multiplication of two bigfloats
  // see lemma in Zilin's thesis
  static prec_t mul_prec(const BigFloat& x, const BigFloat& y) {
    if (get_wbf_mode() == true)
      return get_wbf_prec();

    return x.get_prec() + y.get_prec();

  }

public: // C++ operators
  /// \name unary, increment, decrement operators
  //@{
  /// unary plus operator
  BigFloat operator+() const
  { return BigFloat(*this); }
  /// unary negation operator
  BigFloat operator-() const
  { BigFloat r; r.neg(*this); return r; }
  /// prefix increment operator
  BigFloat& operator++()  
  { add(*this, 1); return *this; }
  /// postfix increment operator
  BigFloat operator++(int) 
  { BigFloat r(*this); ++(*this); return r; }
  /// prefix decrement operator
  BigFloat& operator--()
  { sub(*this, 1); return *this; }
  /// postfix decrement operator
  BigFloat operator--(int)
  { BigFloat r(*this); --(*this); return r; }
  //@}

  /// \name assignment and compound assignment operators
  //@{
  /// assignment operator for <tt>BigFloat</tt>
  BigFloat& operator=(const BigFloat& rhs)
  { base_cls::operator=(rhs); return *this; }
  /// assignment operator for <tt>int</tt>
  BigFloat& operator=(int rhs)
  { set(rhs); return *this; }
  /// assignment operator for <tt>unsigned int</tt>
  BigFloat& operator=(unsigned int rhs)
  { set(rhs); return *this; }
  /// assignment operator for <tt>long</tt>
  BigFloat& operator=(long rhs)
  { set(rhs); return *this; }
  /// assignment operator for <tt>unsigned long</tt>
  BigFloat& operator=(unsigned long rhs)
  { set(rhs); return *this; }
  /// assignment operator for <tt>double</tt>
  BigFloat& operator=(double rhs)
  { set(rhs); return *this; }
  /// assignment operator for <tt>char*</tt>
  BigFloat& operator=(const char* rhs)
  { set(rhs); return *this; }
  /// assignment operator for <tt>std::string</tt>
  BigFloat& operator=(const std::string& rhs)
  { set(rhs); return *this; }
  /// assignment operator for <tt>BigInt</tt>
  BigFloat& operator=(const BigInt& rhs)
  { set(rhs); return *this; }
  /// assignment operator for <tt>BigRat</tt>
  BigFloat& operator=(const BigRat& rhs)
  { set(rhs); return *this; }

  /// compound assignment operator <tt>+=</tt>
  BigFloat& operator+=(const BigFloat& rhs)
  { add(*this, rhs); return *this; }
  /// compound assignment operator <tt>+=</tt>
  BigFloat& operator+=(int rhs)
  { add(*this, rhs); return *this; }
  /// compound assignment operator <tt>+=</tt>
  BigFloat& operator+=(unsigned int rhs)
  { add(*this, rhs); return *this; }
  /// compound assignment operator <tt>+=</tt>
  BigFloat& operator+=(long rhs)
  { add(*this, rhs); return *this; }
  /// compound assignment operator <tt>+=</tt>
  BigFloat& operator+=(unsigned long rhs)
  { add(*this, rhs); return *this; }
  /// compound assignment operator <tt>+=</tt>
  BigFloat& operator+=(double rhs)
  { add(*this, rhs); return *this; }
  /// compound assignment operator <tt>+=</tt>
  BigFloat& operator+=(const BigInt& rhs)
  { add(*this, rhs); return *this; }
  /// compound assignment operator <tt>+=</tt>
  BigFloat& operator+=(const BigRat& rhs)
  { add(*this, rhs); return *this; }

  /// compound assignment operator <tt>-=</tt>
  BigFloat& operator-=(const BigFloat& rhs)
  { sub(*this, rhs); return *this; }
  /// compound assignment operator <tt>-=</tt>
  BigFloat& operator-=(int rhs)
  { sub(*this, rhs); return *this; }
  /// compound assignment operator <tt>-=</tt>
  BigFloat& operator-=(unsigned int rhs)
  { sub(*this, rhs); return *this; }
  /// compound assignment operator <tt>-=</tt>
  BigFloat& operator-=(long rhs)
  { sub(*this, rhs); return *this; }
  /// compound assignment operator <tt>-=</tt>
  BigFloat& operator-=(unsigned long rhs)
  { sub(*this, rhs); return *this; }
  /// compound assignment operator <tt>-=</tt>
  BigFloat& operator-=(double rhs)
  { sub(*this, rhs); return *this; }
  /// compound assignment operator <tt>-=</tt>
  BigFloat& operator-=(const BigInt& rhs)
  { sub(*this, rhs); return *this; }
  /// compound assignment operator <tt>-=</tt>
  BigFloat& operator-=(const BigRat& rhs)
  { sub(*this, rhs); return *this; }

  /// compound assignment operator <tt>*=</tt>
  BigFloat& operator*=(const BigFloat& rhs)
  { mul(*this, rhs); return *this; }
  /// compound assignment operator <tt>*=</tt>
  BigFloat& operator*=(int rhs)
  { mul(*this, rhs); return *this; }
  /// compound assignment operator <tt>*=</tt>
  BigFloat& operator*=(unsigned int rhs)
  { mul(*this, rhs); return *this; }
  /// compound assignment operator <tt>*=</tt>
  BigFloat& operator*=(long rhs)
  { mul(*this, rhs); return *this; }
  /// compound assignment operator <tt>*=</tt>
  BigFloat& operator*=(unsigned long rhs)
  { mul(*this, rhs); return *this; }
  /// compound assignment operator <tt>*=</tt>
  BigFloat& operator*=(double rhs)
  { mul(*this, rhs); return *this; }
  /// compound assignment operator <tt>*=</tt>
  BigFloat& operator*=(const BigInt& rhs)
  { mul(*this, rhs); return *this; }
  /// compound assignment operator <tt>*=</tt>
  BigFloat& operator*=(const BigRat& rhs)
  { mul(*this, rhs); return *this; }

  /// compound assignment operator <tt>/=</tt>
  BigFloat& operator/=(const BigFloat& rhs)
  { div(*this, rhs); return *this; }
  /// compound assignment operator <tt>/=</tt>
  BigFloat& operator/=(int rhs)
  { div(*this, rhs); return *this; }
  /// compound assignment operator <tt>/=</tt>
  BigFloat& operator/=(unsigned int rhs)
  { div(*this, rhs); return *this; }
  /// compound assignment operator <tt>/=</tt>
  BigFloat& operator/=(long rhs)
  { div(*this, rhs); return *this; }
  /// compound assignment operator <tt>/=</tt>
  BigFloat& operator/=(unsigned long rhs)
  { div(*this, rhs); return *this; }
  /// compound assignment operator <tt>/=</tt>
  BigFloat& operator/=(double rhs)
  { div(*this, rhs); return *this; }
  /// compound assignment operator <tt>/=</tt>
  BigFloat& operator/=(const BigInt& rhs)
  { div(*this, rhs); return *this; }
  /// compound assignment operator <tt>/=</tt>
  BigFloat& operator/=(const BigRat& rhs)
  { div(*this, rhs); return *this; }

  /// compound assignment operator <tt><<=</tt>
  BigFloat& operator<<=(int i)
  { return mul_2exp(i); }
  /// compound assignment operator <tt><<=</tt>
  BigFloat& operator<<=(unsigned int ui)
  { return mul_2exp(ui); }
  /// compound assignment operator <tt><<=</tt>
  BigFloat& operator<<=(long l)
  { return mul_2exp(l); }
  /// compound assignment operator <tt><<=</tt>
  BigFloat& operator<<=(unsigned long ul)
  { return mul_2exp(ul); }
  /// compound assignment operator <tt>>>=</tt>
  BigFloat& operator>>=(int i)
  { return div_2exp(i); }
  /// compound assignment operator <tt>>>=</tt>
  BigFloat& operator>>=(unsigned int ui)
  { return div_2exp(ui); }
  /// compound assignment operator <tt>>>=</tt>
  BigFloat& operator>>=(long l)
  { return div_2exp(l); }
  /// compound assignment operator <tt>>>=</tt>
  BigFloat& operator>>=(unsigned long ul)
  { return div_2exp(ul); }
  //@}

public:

  operator BigRat() const { return get_q(); }
#ifndef CORE_DISABLE_OLDNAMES
  /// \name back-compatiable functions
  //@{ 
  /// Has Exact Division
  static bool hasExactDivision() { return false; }
  /// set value from <tt>const char*</tt>
  int set_str(const char* s, int base = 0) { return set(s, base); }
  /// set value from <tt>const string</tt>
  int set_str(const std::string s, int base = 0) { return set(s, base); }
  /// intValue
  int intValue() const { return static_cast<int>(get_si()); }
  /// longValue
  long longValue() const { return get_si(); }
  /// ulongValue 
  unsigned long ulongValue() const { return get_ui(); }
  /// doubleValue
  double doubleValue() const { return get_d(); }
  /// BigIntValue
  BigInt BigIntValue() const { return get_z(); }
  /// BigRatValue
  BigRat BigRatValue() const { return get_q(); }
  //@}
#endif
};



 

 

#include <CORE/BigFloat.inl>

CORE_END_NAMESPACE

#endif /*__CORE_BIGFLOAT_H__*/
