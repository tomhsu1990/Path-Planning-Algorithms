/****************************************************************************
 * BigInt.h -- Big Integer number class based on mpz in GMP
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
 * $Id: BigInt.h,v 1.18 2010/11/23 17:58:36 exact Exp $
 ***************************************************************************/
#ifndef __CORE_BIGINT_H__
#define __CORE_BIGINT_H__

#include <CORE/Gmpz.h>
#include <string>
#include <iostream>
#include <limits>
#include <assert.h>

extern std::istream& extract(std::istream &i, mpz_ptr z);

CORE_BEGIN_NAMESPACE

/* _gmp_alloc_cstr */
struct _gmp_alloc_cstr {
  char *str;
  _gmp_alloc_cstr(int len) { str = new char[len]; }
  ~_gmp_alloc_cstr() { delete[] str; }
};

#ifndef CORE_DISABLE_REFCOUNTING
  typedef RcGmpz BigIntBase;
#else
  typedef Gmpz BigIntBase;
#endif

/// \class BigInt BigInt.h
/// \brief BigInt is a big integer number class based on <tt>mpz</tt> in GMP
class BigInt : public BigIntBase {
  typedef BigIntBase base_cls;
public:
  /// \name constructors and destructor
  //@{
  /// default constructor
  BigInt() {}
  /// constructor for <tt>int</tt>
  BigInt(int i) : base_cls(static_cast<long>(i)) {}
  /// constructor for <tt>unsigned int</tt>
  BigInt(unsigned int i) : base_cls(static_cast<unsigned long>(i)) {}
  /// constructor for <tt>long</tt>
  BigInt(long i) : base_cls(i) {}
  /// constructor for <tt>unsigned long</tt>
  BigInt(unsigned long i) : base_cls(i) {}
  /// constructor for <tt>double</tt>
  BigInt(double i) : base_cls(i) {}
  /// constructor for <tt>char*</tt> (no implicit conversion)
  explicit BigInt(const char* s,int base=0) : base_cls(s,base) {}
  /// constructor for <tt>std::string</tt> (no implicit conversion)
  explicit BigInt(const std::string& s,int base=0) : base_cls(s.c_str(),base){} 
  // internal used by BigRat
  explicit BigInt(mpz_srcptr x) : base_cls(x) {}
  //@}

public:
  /// \name assignment functions
  //@{
  /// assignment function for <tt>BigInt</tt>
  void set(const BigInt& rhs) { base_cls::operator=(rhs); }
  /// assignment function for <tt>int</tt>
  void set(int i) { set(static_cast<long>(i)); }
  /// assignment function for <tt>unsigned int</tt>
  void set(unsigned int i) { set(static_cast<unsigned long>(i)); }
  /// assignment function for <tt>long</tt>
  void set(long i) { mpz_set_si(mp(), i); }
  /// assignment function for <tt>unsigned long</tt>
  void set(unsigned long i) { mpz_set_ui(mp(), i); }
  /// assignment function for <tt>double</tt>
  void set(double i) { mpz_set_d(mp(), i); }
  /// assignment function for <tt>char*</tt>
  int set(const char* str, int base = 0) 
  { return mpz_set_str(mp(), str, base); }
  /// assignment function for <tt>std::string&</tt>
  int set(const std::string& str, int base = 0) 
  { return mpz_set_str(mp(), str.c_str(), base); }
  //@}

  /// \name other arithmetic functions
  //@{
  /// self-negation function
  void neg()
  { mpz_neg(mp(), mp()); }
  /// negation function
  void neg(const BigInt& x)
  { mpz_neg(mp(), x.mp()); }
  /// absolute value function (self)
  void abs()
  { mpz_abs(mp(), mp()); }
  /// absolute value function
  void abs(const BigInt& x)
  { mpz_abs(mp(), x.mp()); }
  //@}

  /// \name arithmetic functions -- addition
  //@{
  /// addition for <tt>BigInt + BigInt</tt>
  void add(const BigInt& x, const BigInt& y) 
  { mpz_add(mp(), x.mp(), y.mp()); }
  /// addition for <tt>BigInt + int</tt>
  void add(const BigInt& x, int y) 
  { add(x, static_cast<long>(y)); }
  /// addition for <tt>BigInt + unsigned int</tt>
  void add(const BigInt& x, unsigned int y) 
  { add(x, static_cast<unsigned long>(y)); }
  /// addition for <tt>BigInt + long</tt>
  void add(const BigInt& x, long y) 
  { if (y>=0) mpz_add_ui(mp(), x.mp(), y); else mpz_sub_ui(mp(), x.mp(), -y); }
  /// addition for <tt>BigInt + unsigned long</tt>
  void add(const BigInt& x, unsigned long y) 
  { mpz_add_ui(mp(), x.mp(), y); }
  /// addition for <tt>BigInt + double</tt>
  void add(const BigInt& x, double y) { add(x, BigInt(y)); }
  /// addition for <tt>int + BigInt</tt>
  void add(int x, const BigInt& y) { add(y, x); }
  /// addition for <tt>unsigned int + BigInt</tt>
  void add(unsigned int x, const BigInt& y) { add(y, x); }
  /// addition for <tt>long + BigInt</tt>
  void add(long x, const BigInt& y) { add(y, x); }
  /// addition for <tt>unsigned long + BigInt</tt>
  void add(unsigned long x, const BigInt& y) { add(y, x); }
  /// addition for <tt>double + BigInt</tt>
  void add(double x, const BigInt& y) { add(y, x); }
  //@}

  /// \name arithmetic functions -- subtraction
  //@{
  /// subtraction for <tt>BigInt - BigInt</tt>
  void sub(const BigInt& x, const BigInt& y) 
  { mpz_sub(mp(), x.mp(), y.mp()); }
  /// subtraction for <tt>BigInt - int</tt>
  void sub(const BigInt& x, int y) 
  { sub(x, static_cast<long>(y)); }
  /// subtraction for <tt>BigInt - unsigned int</tt>
  void sub(const BigInt& x, unsigned int y) 
  { sub(x, static_cast<unsigned long>(y)); }
  /// subtraction for <tt>BigInt - long</tt>
  void sub(const BigInt& x, long y) 
  { if (y>=0) mpz_sub_ui(mp(), x.mp(), y); else mpz_add_ui(mp(), x.mp(), -y); }
  /// subtraction for <tt>BigInt - unsigned long</tt>
  void sub(const BigInt& x, unsigned long y) 
  { mpz_sub_ui(mp(), x.mp(), y); }
  /// subtraction for <tt>BigInt - double</tt>
  void sub(const BigInt& x, double y) { sub(x, BigInt(y)); }
  /// subtraction for <tt>int - BigInt</tt>
  void sub(int x, const BigInt& y) 
  { sub(static_cast<long>(x), y); }
  /// subtraction for <tt>unsigned int - BigInt</tt>
  void sub(unsigned int x, const BigInt& y) 
  { sub(static_cast<unsigned long>(x), y); }
  /// subtraction for <tt>long - BigInt</tt>
  void sub(long x, const BigInt& y) 
  { if (x>=0) mpz_ui_sub(mp(), x, y.mp()); else {add(y, -x); neg();} }
  /// subtraction for <tt>unsigned long - BigInt</tt>
  void sub(unsigned long x, const BigInt& y) 
  { mpz_ui_sub(mp(), x, y.mp()); }
  /// subtraction for <tt>double - BigInt</tt>
  void sub(double x, const BigInt& y) { sub(BigInt(x), y); }
  //@}

  /// \name arithmetic functions -- multiplication
  //@{
  /// multiplication for <tt>BigInt * BigInt</tt>
  void mul(const BigInt& x, const BigInt& y)
  { mpz_mul(mp(), x.mp(), y.mp()); }
  /// multiplication for <tt>BigInt * int</tt>
  void mul(const BigInt& x, int y)
  { mul(x, static_cast<long>(y)); }
  /// multiplication for <tt>BigInt * unsigned int</tt>
  void mul(const BigInt& x, unsigned int y)
  { mul(x, static_cast<unsigned long>(y)); }
  /// multiplication for <tt>BigInt * long</tt>
  void mul(const BigInt& x, long y)
  { if (y>=0) mpz_mul_ui(mp(), x.mp(), y); else {mul(x, -y); neg();} }
  /// multiplication for <tt>BigInt * unsigned long</tt>
  void mul(const BigInt& x, unsigned long y)
  { mpz_mul_ui(mp(), x.mp(), y); }
  /// multiplication for <tt>BigInt * double</tt>
  void mul(const BigInt& x, double y) { mul(x, BigInt(y)); }
  /// multiplication for <tt>int * BigInt</tt>
  void mul(int x, const BigInt& y) { mul(y, x); }
  /// multiplication for <tt>unsigned int * BigInt</tt>
  void mul(unsigned int x, const BigInt& y) { mul(y, x); }
  /// multiplication for <tt>long * BigInt</tt>
  void mul(long x, const BigInt& y) { mul(y, x); }
  /// multiplication for <tt>unsigned long * BigInt</tt>
  void mul(unsigned long x, const BigInt& y) { mul(y, x); }
  /// multiplication for <tt>double * BigInt</tt>
  void mul(double x, const BigInt& y) { mul(BigInt(x), y); }
  //@}

  /// \name arithmetic functions -- division
  //@{
  /// division for <tt>BigInt / BigInt</tt>
  void div(const BigInt& x, const BigInt& y)
  { mpz_tdiv_q(mp(), x.mp(), y.mp()); }
  /// division for <tt>BigInt / int</tt>
  void div(const BigInt& x, int y)
  { div(x, static_cast<long>(y)); }
  /// division for <tt>BigInt / unsigned int</tt>
  void div(const BigInt& x, unsigned int y)
  { div(x, static_cast<unsigned long>(y)); }
  /// division for <tt>BigInt / long</tt>
  void div(const BigInt& x, long y)
  { if (y>=0) mpz_tdiv_q_ui(mp(), x.mp(), y); else {div(x, -y); neg();} }
  /// division for <tt>BigInt / unsigned long</tt>
  void div(const BigInt& x, unsigned long y)
  { mpz_tdiv_q_ui(mp(), x.mp(), y); }
  /// division for <tt>BigInt / double</tt>
  void div(const BigInt& x, double y) { div(x, BigInt(y)); }
  /// division for <tt>int / BigInt</tt>
  void div(int x, const BigInt& y) { div(BigInt(x), y); }
  /// division for <tt>unsigned int / BigInt</tt>
  void div(unsigned int x, const BigInt& y) { div(BigInt(x), y); }
  /// division for <tt>long / BigInt</tt>
  void div(long x, const BigInt& y) { div(BigInt(x), y); }
  /// division for <tt>unsigned long / BigInt</tt>
  void div(unsigned long x, const BigInt& y) { div(BigInt(x), y); }
  /// division for <tt>double / BigInt</tt>
  void div(double x, const BigInt& y) { div(BigInt(x), y); }
  //@}

  /// \name arithmetic functions -- modular
  //@{
  /// modular for <tt>BigInt % BigInt</tt>
  void mod(const BigInt& x, const BigInt& y)
  { mpz_tdiv_r(mp(), x.mp(), y.mp()); }
  /// modular for <tt>BigInt % int</tt>
  void mod(const BigInt& x, int y)
  { mod(x, static_cast<long>(y)); }
  /// modular for <tt>BigInt % unsigned int</tt>
  void mod(const BigInt& x, unsigned int y)
  { mod(x, static_cast<unsigned long>(y)); }
  /// modular for <tt>BigInt % long</tt>
  void mod(const BigInt& x, long y)
  { if (y>=0) mpz_tdiv_r_ui(mp(), x.mp(), y); else {mod(x, -y);} }
  /// modular for <tt>BigInt % unsigned long</tt>
  void mod(const BigInt& x, unsigned long y)
  { mpz_tdiv_r_ui(mp(), x.mp(), y); }
  /// modular for <tt>BigInt % double</tt>
  void mod(const BigInt& x, double y) { mod(x, BigInt(y)); }
  /// modular for <tt>int % BigInt</tt>
  void mod(int x, const BigInt& y) { mod(BigInt(x), y); }
  /// modular for <tt>unsigned int % BigInt</tt>
  void mod(unsigned int x, const BigInt& y) { mod(BigInt(x), y); }
  /// modular for <tt>long % BigInt</tt>
  void mod(long x, const BigInt& y) { mod(BigInt(x), y); }
  /// modular for <tt>unsigned long % BigInt</tt>
  void mod(unsigned long x, const BigInt& y) { mod(BigInt(x), y); }
  /// modular for <tt>double % BigInt</tt>
  void mod(double x, const BigInt& y) { mod(BigInt(x), y); }
  //@}

  /// \name arithmetic functions -- exact division
  //@{
  /// exact division for <tt>BigInt / BigInt</tt>
  void divexact(const BigInt& x, const BigInt& y)
  { mpz_divexact(mp(), x.mp(), y.mp()); }
  /// exact division for <tt>BigInt / int</tt>
  void divexact(const BigInt& x, int y)
  { divexact(x, static_cast<long>(y)); }
  /// exact division for <tt>BigInt / unsigned int</tt>
  void divexact(const BigInt& x, unsigned int y)
  { divexact(x, static_cast<unsigned long>(y)); }
  /// exact division for <tt>BigInt / long</tt>
  void divexact(const BigInt& x, long y)
  { if (y>=0) mpz_divexact_ui(mp(), x.mp(), y); else {divexact(x,-y); neg();} }
  /// exact division for <tt>BigInt / unsigned long</tt>
  void divexact(const BigInt& x, unsigned long y)
  { mpz_divexact_ui(mp(), x.mp(), y); }
  /// exact division for <tt>BigInt / double</tt>
  void divexact(const BigInt& x, double y) { divexact(x, BigInt(y)); }
  //@}

  /// \name arithmetic functions -- division with remainder
  //@{
  /// division with remainder for <tt>BigInt / BigInt</tt>
  void divrem(BigInt& r, const BigInt& x, const BigInt& y)
  { mpz_fdiv_qr(mp(), r.mp(), x.mp(), y.mp()); }
  /// division with remainder for <tt>BigInt / int</tt>
  void divrem(BigInt& r, const BigInt& x, int y)
  { divrem(r, x, static_cast<long>(y)); }
  /// division with remainder for <tt>BigInt / unsigned int</tt>
  void divrem(BigInt& r, const BigInt& x, unsigned int y)
  { divrem(r, x, static_cast<unsigned long>(y)); }
  /// division with remainder for <tt>BigInt / long</tt>
  void divrem(BigInt& r, const BigInt& x, long y)
  {if (y>=0) mpz_fdiv_qr_ui(mp(),r.mp(),x.mp(),y); else {divrem(r,x,-y);neg();}}
  /// division with remainder for <tt>BigInt / unsigned long</tt>
  void divrem(BigInt& r, const BigInt& x, unsigned long y)
  { mpz_fdiv_qr_ui(mp(), r.mp(), x.mp(), y); }
  /// division with remainder for <tt>BigInt / double</tt>
  void divrem(BigInt& r, const BigInt& x, double y) { divrem(r, x, BigInt(y)); }

  /// \name squart root function
  //@{
  /// square root for <tt>BigInt</tt>
  void sqrt(const BigInt& x)
  { mpz_sqrt(mp(), x.mp()); }
  //@}

  /// \name power functions
  //@{
  /// power function for <tt>BigInt</tt>
  void pow(const BigInt& x, unsigned long y)
  { mpz_pow_ui(mp(), x.mp(), y); }
  /// power function for <tt>int</tt>
  void pow(int x, unsigned long y)
  { pow(static_cast<long>(x), y); }
  /// power function for <tt>unsigned int</tt>
  void pow(unsigned int x, unsigned long y)
  { pow(static_cast<unsigned long>(x), y); }
  /// power function for <tt>long</tt>
  void pow(long x, unsigned long y)
  { mpz_ui_pow_ui(mp(), (x>=0?x:-x), y); if (x<0&&x%2!=0) neg(); }
  /// power function for <tt>unsigned long</tt>
  void pow(unsigned long x, unsigned long y)
  { mpz_ui_pow_ui(mp(), x, y); }
  //@}
  
  /// \name shift functions
  //@{
  /// left shift
  void mul_2exp(const BigInt& x, int y)
  { mul_2exp(x, static_cast<long>(y)); }
  /// left shift
  void mul_2exp(const BigInt& x, unsigned int y)
  { mul_2exp(x, static_cast<unsigned long>(y)); }
  /// left shift
  void mul_2exp(const BigInt& x, long y)
  { if (y>=0) mpz_mul_2exp(mp(),x.mp(),y); else mpz_div_2exp(mp(),x.mp(),-y); }
  /// left shift
  void mul_2exp(const BigInt& x, unsigned long y)
  { mpz_mul_2exp(mp(), x.mp(), y); }
  /// right shift
  void div_2exp(const BigInt& x, int y)
  { div_2exp(x, static_cast<long>(y)); }
  /// right shift
  void div_2exp(const BigInt& x, unsigned int y)
  { div_2exp(x, static_cast<unsigned long>(y)); }
  /// right shift
  void div_2exp(const BigInt& x, long y)
  { if (y>=0) mpz_div_2exp(mp(),x.mp(),y); else mpz_mul_2exp(mp(),x.mp(),-y);}
  /// right shift
  void div_2exp(const BigInt& x, unsigned long y)
  { mpz_div_2exp(mp(), x.mp(), y); }
  //@}

  /// \name comparison functions
  //@{
  /// compare with <tt>BigInt</tt>
  int cmp(const BigInt& x) const
  { return mpz_cmp(mp(), x.mp()); }
  /// compare with <tt>int</tt>
  int cmp(int x) const
  { return cmp(static_cast<long>(x)); }
  /// compare with <tt>unsigned int</tt>
  int cmp(unsigned int x) const
  { return cmp(static_cast<unsigned long>(x)); }
  /// compare with <tt>long</tt>
  int cmp(long x) const
  { return mpz_cmp_si(mp(), x); }
  /// compare with <tt>unsigned long</tt>
  int cmp(unsigned long x) const
  { return mpz_cmp_ui(mp(), x); }
  /// compare with <tt>double</tt>
  int cmp(double x) const
  { return mpz_cmp_d(mp(), x); }
  //@}

  /// \name comparison functions (in absolute value)
  //@{
  /// compare (in absolute value) with <tt>BigInt</tt>
  int cmpabs(const BigInt& x) const
  { return mpz_cmpabs(mp(), x.mp()); }
  /// compare (in absolute value) with <tt>int</tt>
  int cmpabs(int x) const
  { return cmpabs(static_cast<long>(x)); }
  /// compare (in absolute value) with <tt>unsigned int</tt>
  int cmpabs(unsigned int x) const
  { return cmpabs(static_cast<unsigned long>(x)); }
  /// compare (in absolute value) with <tt>long</tt>
  int cmpabs(long x) const
  { return mpz_cmpabs_ui(mp(), x>=0? x:-x); }
  /// compare (in absolute value) with <tt>unsigned long</tt>
  int cmpabs(unsigned long x) const
  { return mpz_cmpabs_ui(mp(), x); }
  /// compare (in absolute value) with <tt>double</tt>
  int cmpabs(double x) const
  { return mpz_cmp_d(mp(), x); }
  //@}

  /// \name logical and bit manipulation functions
  //@{
  /// logical and
  void logical_and(const BigInt& x, const BigInt& y)
  { mpz_and(mp(), x.mp(), y.mp()); }
  /// logical ior
  void logical_ior(const BigInt& x, const BigInt& y)
  { mpz_ior(mp(), x.mp(), y.mp()); }
  /// logical xor
  void logical_xor(const BigInt& x, const BigInt& y)
  { mpz_xor(mp(), x.mp(), y.mp()); }
  /// logical com
  void logical_com(const BigInt& x)
  { mpz_com(mp(), x.mp()); }
  //@}
  
  /// \name conversion functions
  //@{
  /// return double value
  double get_d() const
  { return mpz_get_d(mp()); }
  /// find d and exp s.t. \f$d*2^{exp}\f$ with \f$0.5\le|d|<1\f$
  double get_d_2exp(long* exp) const
  { return mpz_get_d_2exp(exp, mp()); }
  /// return long value
  long get_si() const
  { return mpz_get_si(mp()); }
  /// return unsigned long value
  unsigned long get_ui() const
  { return mpz_get_ui(mp()); }
  /// return the string representation
  std::string get_str(int base = 10) const {
    int len = mpz_sizeinbase(mp(), base) + 2;
    _gmp_alloc_cstr tmp(len);
    return std::string(mpz_get_str(tmp.str, base, mp()));
  }
  /// get exponent of power 2
  unsigned long get_2exp() const
  { return mpz_scan1(mp(), 0); }
  /// get exponent of power k
  unsigned long get_k_exp(BigInt& m, unsigned long k) const
  { return mpz_remove(m.mp(), mp(), BigInt(k).mp()); }
  //@}
  
  /// \name miscellaneous functions
  //@{
  /// swap function 
  void swap(BigInt& other)
  { mpz_swap(mp(), other.mp()); }
  /// gcd function
  void gcd(const BigInt& x, const BigInt& y)
  { mpz_gcd(mp(), x.mp(), y.mp()); }
  /// return size in base
  size_t sizeinbase(int base = 2) const
  { return mpz_sizeinbase(mp(), base); }
  /// return \f$\lceil\log|x|\rceil\f$
  unsigned long ceillg() const
  { unsigned long len=sizeinbase(); return (get_2exp()==len-1)?(len-1):len; }
  /// return \f$\lfloor\log|x|\lfloor\f$
  unsigned long floorlg() const
  { return sizeinbase() - 1; }
  /// return sign
  int sgn() const
  { return mpz_sgn(mp()); }
  /// return upper bound of MSB
  long uMSB() const
  { return ceillg(); } 
  /// return lower bound of MSB
  long lMSB() const
  { return floorlg(); }
  /// binomial coefficient
  void binomial(unsigned long n, unsigned long k)
  { mpz_bin_uiui(mp(), n, k); } 
  //@}

  /// \name helper functions
  //@{
  /// return true if it is divisible
  bool is_divisible(const BigInt& x) const
  { return mpz_divisible_p(mp(), x.mp()) != 0; }  
  /// return true if it is odd
  bool is_odd() const
  { return mpz_odd_p(mp()) != 0; }
  /// return true if it is even
  bool is_even() const
  { return mpz_even_p(mp()) != 0; }
  /// return true if it fits unsigned long
  bool is_ulong() const
  { return mpz_fits_ulong_p(mp()) != 0; }
  /// return true if it fits signed long
  bool is_slong() const
  { return mpz_fits_slong_p(mp()) != 0; }
  /// return true if it fits unsigned int
  bool is_uint() const
  { return mpz_fits_uint_p(mp()) != 0; }
  /// return true if it fits signed int
  bool is_sint() const
  { return mpz_fits_sint_p(mp()) != 0; }
  /// return true if it fits unsigned short
  bool is_ushort() const
  { return mpz_fits_ushort_p(mp()) != 0; }
  /// return true if it fits signed short
  bool is_sshort() const
  { return mpz_fits_sshort_p(mp()) != 0; }
  //@}

public: // C++ operators
  /// \name unary, increment, decrement operators
  //@{
  /// unary plus operator
  BigInt operator+() const
  { return BigInt(*this); }
  /// unary negation operator
  BigInt operator-() const
  { BigInt r; r.neg(*this); return r; }
  /// prefix increment operator
  BigInt& operator++()
  { add(*this, 1); return *this; }
  /// postfix increment operator
  BigInt operator++(int)
  { BigInt r(*this); ++(*this); return r; }
  /// prefix decrement operator
  BigInt& operator--()
  { sub(*this, 1); return *this; }
  /// postfix decrement operator
  BigInt operator--(int)
  { BigInt r(*this); --(*this); return r; }
  //@}

  /// \name assignment and compound assignment operators
  //@{
  /// assignment operator for <tt>BigInt</tt>
  BigInt& operator=(const BigInt& rhs)
  { base_cls::operator=(rhs); return *this; }
  /// assignment operator for <tt>int</tt>
  BigInt& operator=(int rhs)
  { set(rhs); return *this; }
  /// assignment operator for <tt>unsigned int</tt>
  BigInt& operator=(unsigned int rhs)
  { set(rhs); return *this; }
  /// assignment operator for <tt>long</tt>
  BigInt& operator=(long rhs)
  { set(rhs); return *this; }
  /// assignment operator for <tt>unsigned long</tt>
  BigInt& operator=(unsigned long rhs)
  { set(rhs); return *this; }
  /// assignment operator for <tt>double</tt>
  BigInt& operator=(double rhs)
  { set(rhs); return *this; }
  /// assignment operator for <tt>char*</tt>
  BigInt& operator=(const char* rhs)
  { set(rhs); return *this; }
  /// assignment operator for <tt>std::string</tt>
  BigInt& operator=(const std::string& rhs)
  { set(rhs); return *this; }

  /// compound assignment operator <tt>+=</tt>
  BigInt& operator+=(const BigInt& rhs)
  { add(*this, rhs); return *this; }
  /// compound assignment operator <tt>+=</tt>
  BigInt& operator+=(int rhs)
  { add(*this, rhs); return *this; }
  /// compound assignment operator <tt>+=</tt>
  BigInt& operator+=(unsigned int rhs)
  { add(*this, rhs); return *this; }
  /// compound assignment operator <tt>+=</tt>
  BigInt& operator+=(long rhs)
  { add(*this, rhs); return *this; }
  /// compound assignment operator <tt>+=</tt>
  BigInt& operator+=(unsigned long rhs)
  { add(*this, rhs); return *this; }
  /// compound assignment operator <tt>+=</tt>
  BigInt& operator+=(double rhs)
  { add(*this, rhs); return *this; }

  /// compound assignment operator <tt>-=</tt>
  BigInt& operator-=(const BigInt& rhs)
  { sub(*this, rhs); return *this; }
  /// compound assignment operator <tt>-=</tt>
  BigInt& operator-=(int rhs)
  { sub(*this, rhs); return *this; }
  /// compound assignment operator <tt>-=</tt>
  BigInt& operator-=(unsigned int rhs)
  { sub(*this, rhs); return *this; }
  /// compound assignment operator <tt>-=</tt>
  BigInt& operator-=(long rhs)
  { sub(*this, rhs); return *this; }
  /// compound assignment operator <tt>-=</tt>
  BigInt& operator-=(unsigned long rhs)
  { sub(*this, rhs); return *this; }
  /// compound assignment operator <tt>-=</tt>
  BigInt& operator-=(double rhs)
  { sub(*this, rhs); return *this; }

  /// compound assignment operator <tt>*=</tt>
  BigInt& operator*=(const BigInt& rhs)
  { mul(*this, rhs); return *this; }
  /// compound assignment operator <tt>*=</tt>
  BigInt& operator*=(int rhs)
  { mul(*this, rhs); return *this; }
  /// compound assignment operator <tt>*=</tt>
  BigInt& operator*=(unsigned int rhs)
  { mul(*this, rhs); return *this; }
  /// compound assignment operator <tt>*=</tt>
  BigInt& operator*=(long rhs)
  { mul(*this, rhs); return *this; }
  /// compound assignment operator <tt>*=</tt>
  BigInt& operator*=(unsigned long rhs)
  { mul(*this, rhs); return *this; }
  /// compound assignment operator <tt>*=</tt>
  BigInt& operator*=(double rhs)
  { mul(*this, rhs); return *this; }

  /// compound assignment operator <tt>/=</tt>
  BigInt& operator/=(const BigInt& rhs)
  { div(*this, rhs); return *this; }
  /// compound assignment operator <tt>/=</tt>
  BigInt& operator/=(int rhs)
  { div(*this, rhs); return *this; }
  /// compound assignment operator <tt>/=</tt>
  BigInt& operator/=(unsigned int rhs)
  { div(*this, rhs); return *this; }
  /// compound assignment operator <tt>/=</tt>
  BigInt& operator/=(long rhs)
  { div(*this, rhs); return *this; }
  /// compound assignment operator <tt>/=</tt>
  BigInt& operator/=(unsigned long rhs)
  { div(*this, rhs); return *this; }
  /// compound assignment operator <tt>/=</tt>
  BigInt& operator/=(double rhs)
  { div(*this, rhs); return *this; }

  /// compound assignment operator <tt>%=</tt>
  BigInt& operator%=(const BigInt& rhs)
  { mod(*this, rhs); return *this; }
  /// compound assignment operator <tt>%=</tt>
  BigInt& operator%=(int rhs)
  { mod(*this, rhs); return *this; }
  /// compound assignment operator <tt>%=</tt>
  BigInt& operator%=(unsigned int rhs)
  { mod(*this, rhs); return *this; }
  /// compound assignment operator <tt>%=</tt>
  BigInt& operator%=(long rhs)
  { mod(*this, rhs); return *this; }
  /// compound assignment operator <tt>%=</tt>
  BigInt& operator%=(unsigned long rhs)
  { mod(*this, rhs); return *this; }
  /// compound assignment operator <tt>%=</tt>
  BigInt& operator%=(double rhs)
  { mod(*this, rhs); return *this; }

  /// compound assignment operator <tt>&=</tt>
  BigInt& operator&=(const BigInt& rhs)
  { logical_and(*this, rhs); return *this; }
  /// compound assignment operator <tt>|=</tt>
  BigInt& operator|=(const BigInt& rhs)
  { logical_ior(*this, rhs); return *this; }
  /// compound assignment operator <tt>^=</tt>
  BigInt& operator^=(const BigInt& rhs)
  { logical_xor(*this, rhs); return *this; }
  /// compound assignment operator <tt><<=</tt>
  BigInt& operator<<=(int i)
  { mul_2exp(*this, i); return *this; }
  /// compound assignment operator <tt><<=</tt>
  BigInt& operator<<=(unsigned int ui)
  { mul_2exp(*this, ui); return *this; }
  /// compound assignment operator <tt><<=</tt>
  BigInt& operator<<=(long l)
  { mul_2exp(*this, l); return *this; }
  /// compound assignment operator <tt><<=</tt>
  BigInt& operator<<=(unsigned long ul)
  { mul_2exp(*this, ul); return *this; }
  /// compound assignment operator <tt>>>=</tt>
  BigInt& operator>>=(int i)
  { div_2exp(*this, i); return *this; }
  /// compound assignment operator <tt>>>=</tt>
  BigInt& operator>>=(unsigned int ui)
  { div_2exp(*this, ui); return *this; }
  /// compound assignment operator <tt>>>=</tt>
  BigInt& operator>>=(long l)
  { div_2exp(*this, l); return *this; }
  /// compound assignment operator <tt>>>=</tt>
  BigInt& operator>>=(unsigned long ul)
  { div_2exp(*this, ul); return *this; }
  //@}

#ifndef CORE_DISABLE_OLDNAMES
  /// \name back-compatiable functions
  //@{
  /// Has Exact Division
  static bool hasExactDivision() { return false; } 
  /// set value from <tt>const char*</tt>
  int set_str(const char* s, int base = 0) { return set(s, base); }
  /// intValue
  int intValue() const { return static_cast<int>(get_si()); }
  /// longValue
  long longValue() const { return get_si(); }
  /// ulongValue
  unsigned long ulongValue() const { return get_ui(); }
  /// doubleValue
  double doubleValue() const { return get_d(); }
  //@}
#endif
};

/// \addtogroup BigIntArithmeticOperators
//@{
/// BigInt + BigInt
inline BigInt operator+(const BigInt& x, const BigInt& y)
{ BigInt r; r.add(x, y); return r; }
/// BigInt + int
inline BigInt operator+(const BigInt& x, int y)
{ BigInt r; r.add(x, y); return r; }
/// int + BigInt
inline BigInt operator+(int x, const BigInt& y)
{ BigInt r; r.add(x, y); return r; }
/// BigInt + unsigned int
inline BigInt operator+(const BigInt& x, unsigned int y)
{ BigInt r; r.add(x, y); return r; }
/// unsigned int + BigInt
inline BigInt operator+(unsigned int x, const BigInt& y)
{ BigInt r; r.add(x, y); return r; }
/// BigInt + long
inline BigInt operator+(const BigInt& x, long y)
{ BigInt r; r.add(x, y); return r; }
/// long + BigInt
inline BigInt operator+(long x, const BigInt& y)
{ BigInt r; r.add(x, y); return r; }
/// BigInt + unsigned long
inline BigInt operator+(const BigInt& x, unsigned long y)
{ BigInt r; r.add(x, y); return r; }
/// unsigned long + BigInt
inline BigInt operator+(unsigned long x, const BigInt& y)
{ BigInt r; r.add(x, y); return r; }
/// BigInt + double
inline BigInt operator+(const BigInt& x, double y)
{ BigInt r; r.add(x, y); return r; }
/// double + BigInt
inline BigInt operator+(double x, const BigInt& y)
{ BigInt r; r.add(x, y); return r; }

/// BigInt - BigInt
inline BigInt operator-(const BigInt& x, const BigInt& y)
{ BigInt r; r.sub(x, y); return r; }
/// BigInt - int
inline BigInt operator-(const BigInt& x, int y)
{ BigInt r; r.sub(x, y); return r; }
/// int - BigInt
inline BigInt operator-(int x, const BigInt& y)
{ BigInt r; r.sub(x, y); return r; }
/// BigInt - unsigned int
inline BigInt operator-(const BigInt& x, unsigned int y)
{ BigInt r; r.sub(x, y); return r; }
/// unsigned int - BigInt
inline BigInt operator-(unsigned int x, const BigInt& y)
{ BigInt r; r.sub(x, y); return r; }
/// BigInt - long
inline BigInt operator-(const BigInt& x, long y)
{ BigInt r; r.sub(x, y); return r; }
/// long - BigInt
inline BigInt operator-(long x, const BigInt& y)
{ BigInt r; r.sub(x, y); return r; }
/// BigInt - unsigned long
inline BigInt operator-(const BigInt& x, unsigned long y)
{ BigInt r; r.sub(x, y); return r; }
/// unsigned long - BigInt
inline BigInt operator-(unsigned long x, const BigInt& y)
{ BigInt r; r.sub(x, y); return r; }
/// BigInt - double
inline BigInt operator-(const BigInt& x, double y)
{ BigInt r; r.sub(x, y); return r; }
/// double - BigInt
inline BigInt operator-(double x, const BigInt& y)
{ BigInt r; r.sub(x, y); return r; }

/// BigInt * BigInt
inline BigInt operator*(const BigInt& x, const BigInt& y)
{ BigInt r; r.mul(x, y); return r; }
/// BigInt * int
inline BigInt operator*(const BigInt& x, int y)
{ BigInt r; r.mul(x, y); return r; }
/// int * BigInt
inline BigInt operator*(int x, const BigInt& y)
{ BigInt r; r.mul(x, y); return r; }
/// BigInt * unsigned int
inline BigInt operator*(const BigInt& x, unsigned int y)
{ BigInt r; r.mul(x, y); return r; }
/// unsigned int * BigInt
inline BigInt operator*(unsigned int x, const BigInt& y)
{ BigInt r; r.mul(x, y); return r; }
/// BigInt * long
inline BigInt operator*(const BigInt& x, long y)
{ BigInt r; r.mul(x, y); return r; }
/// long * BigInt
inline BigInt operator*(long x, const BigInt& y)
{ BigInt r; r.mul(x, y); return r; }
/// BigInt * unsigned long
inline BigInt operator*(const BigInt& x, unsigned long y)
{ BigInt r; r.mul(x, y); return r; }
/// unsigned long * BigInt
inline BigInt operator*(unsigned long x, const BigInt& y)
{ BigInt r; r.mul(x, y); return r; }
/// BigInt * double
inline BigInt operator*(const BigInt& x, double y)
{ BigInt r; r.mul(x, y); return r; }
/// double * BigInt
inline BigInt operator*(double x, const BigInt& y)
{ BigInt r; r.mul(x, y); return r; }

/// BigInt / BigInt
inline BigInt operator/(const BigInt& x, const BigInt& y)
{ BigInt r; r.div(x, y); return r; }
/// BigInt / int
inline BigInt operator/(const BigInt& x, int y)
{ BigInt r; r.div(x, y); return r; }
/// int / BigInt
inline BigInt operator/(int x, const BigInt& y)
{ BigInt r; r.div(x, y); return r; }
/// BigInt / unsigned int
inline BigInt operator/(const BigInt& x, unsigned int y)
{ BigInt r; r.div(x, y); return r; }
/// unsigned int / BigInt
inline BigInt operator/(unsigned int x, const BigInt& y)
{ BigInt r; r.div(x, y); return r; }
/// BigInt / long
inline BigInt operator/(const BigInt& x, long y)
{ BigInt r; r.div(x, y); return r; }
/// long / BigInt
inline BigInt operator/(long x, const BigInt& y)
{ BigInt r; r.div(x, y); return r; }
/// BigInt / unsigned long
inline BigInt operator/(const BigInt& x, unsigned long y)
{ BigInt r; r.div(x, y); return r; }
/// unsigned long / BigInt
inline BigInt operator/(unsigned long x, const BigInt& y)
{ BigInt r; r.div(x, y); return r; }
/// BigInt / double
inline BigInt operator/(const BigInt& x, double y)
{ BigInt r; r.div(x, y); return r; }
/// double / BigInt
inline BigInt operator/(double x, const BigInt& y)
{ BigInt r; r.div(x, y); return r; }

/// BigInt % BigInt
inline BigInt operator%(const BigInt& x, const BigInt& y)
{ BigInt r; r.mod(x, y); return r; }
/// BigInt % int
inline BigInt operator%(const BigInt& x, int y)
{ BigInt r; r.mod(x, y); return r; }
/// int % BigInt
inline BigInt operator%(int x, const BigInt& y)
{ BigInt r; r.mod(x, y); return r; }
/// BigInt % unsigned int
inline BigInt operator%(const BigInt& x, unsigned int y)
{ BigInt r; r.mod(x, y); return r; }
/// unsigned int % BigInt
inline BigInt operator%(unsigned int x, const BigInt& y)
{ BigInt r; r.mod(x, y); return r; }
/// BigInt % long
inline BigInt operator%(const BigInt& x, long y)
{ BigInt r; r.mod(x, y); return r; }
/// long % BigInt
inline BigInt operator%(long x, const BigInt& y)
{ BigInt r; r.mod(x, y); return r; }
/// BigInt % unsigned long
inline BigInt operator%(const BigInt& x, unsigned long y)
{ BigInt r; r.mod(x, y); return r; }
/// unsigned long % BigInt
inline BigInt operator%(unsigned long x, const BigInt& y)
{ BigInt r; r.mod(x, y); return r; }
/// BigInt % double
inline BigInt operator%(const BigInt& x, double y)
{ BigInt r; r.mod(x, y); return r; }
/// double % BigInt
inline BigInt operator%(double x, const BigInt& y)
{ BigInt r; r.mod(x, y); return r; }

/// BigInt & BigInt
inline BigInt operator&(const BigInt& x, const BigInt& y)
{ BigInt r; r.logical_and(x, y); return r; }
/// BigInt | BigInt
inline BigInt operator|(const BigInt& x, const BigInt& y)
{ BigInt r; r.logical_ior(x, y); return r; }
/// BigInt ^ BigInt
inline BigInt operator^(const BigInt& x, const BigInt& y)
{ BigInt r; r.logical_xor(x, y); return r; }

/// BigInt << int
inline BigInt operator<<(const BigInt& x, int y)
{ BigInt r; r.mul_2exp(x, y); return r; }
/// BigInt << unsigned int
inline BigInt operator<<(const BigInt& x, unsigned int y)
{ BigInt r; r.mul_2exp(x, y); return r; }
/// BigInt << long
inline BigInt operator<<(const BigInt& x, long y)
{ BigInt r; r.mul_2exp(x, y); return r; }
/// BigInt << unsigned long
inline BigInt operator<<(const BigInt& x, unsigned long y)
{ BigInt r; r.mul_2exp(x, y); return r; }
/// BigInt >> int
inline BigInt operator>>(const BigInt& x, int y)
{ BigInt r; r.div_2exp(x, y); return r; }
/// BigInt >> unsigned int
inline BigInt operator>>(const BigInt& x, unsigned int y)
{ BigInt r; r.div_2exp(x, y); return r; }
/// BigInt >> long
inline BigInt operator>>(const BigInt& x, long y)
{ BigInt r; r.div_2exp(x, y); return r; }
/// BigInt >> unsigned long
inline BigInt operator>>(const BigInt& x, unsigned long y)
{ BigInt r; r.div_2exp(x, y); return r; }
//@}

/// \addtogroup BigIntComparisonOperators
//@{
/// BigInt == BigInt
inline bool operator==(const BigInt& x, const BigInt& y)
{ return x.cmp(y) == 0; }
/// BigInt == int
inline bool operator==(const BigInt& x, int y)
{ return x.cmp(y) == 0; }
/// int == BigInt
inline bool operator==(int x, const BigInt& y)
{ return y.cmp(x) == 0; }
/// BigInt == unsigned int
inline bool operator==(const BigInt& x, unsigned int y)
{ return x.cmp(y) == 0; }
/// unsigned int == BigInt
inline bool operator==(unsigned int x, const BigInt& y)
{ return y.cmp(x) == 0; }
/// BigInt == long
inline bool operator==(const BigInt& x, long y)
{ return x.cmp(y) == 0; }
/// long == BigInt
inline bool operator==(long x, const BigInt& y)
{ return y.cmp(x) == 0; }
/// BigInt == unsigned long
inline bool operator==(const BigInt& x, unsigned long y)
{ return x.cmp(y) == 0; }
/// unsigned long == BigInt
inline bool operator==(unsigned long x, const BigInt& y)
{ return y.cmp(x) == 0; }
/// BigInt == double
inline bool operator==(const BigInt& x, double y)
{ return x.cmp(y) == 0; }
/// double == BigInt
inline bool operator==(double x, const BigInt& y)
{ return y.cmp(x) == 0; }

/// BigInt != BigInt
inline bool operator!=(const BigInt& x, const BigInt& y)
{ return x.cmp(y) != 0; }
/// BigInt != int
inline bool operator!=(const BigInt& x, int y)
{ return x.cmp(y) != 0; }
/// int != BigInt
inline bool operator!=(int x, const BigInt& y)
{ return y.cmp(x) != 0; }
/// BigInt != unsigned int
inline bool operator!=(const BigInt& x, unsigned int y)
{ return x.cmp(y) != 0; }
/// unsigned int != BigInt
inline bool operator!=(unsigned int x, const BigInt& y)
{ return y.cmp(x) != 0; }
/// BigInt != long
inline bool operator!=(const BigInt& x, long y)
{ return x.cmp(y) != 0; }
/// long != BigInt
inline bool operator!=(long x, const BigInt& y)
{ return y.cmp(x) != 0; }
/// BigInt != unsigned long
inline bool operator!=(const BigInt& x, unsigned long y)
{ return x.cmp(y) != 0; }
/// unsigned long != BigInt
inline bool operator!=(unsigned long x, const BigInt& y)
{ return y.cmp(x) != 0; }
/// BigInt != double
inline bool operator!=(const BigInt& x, double y)
{ return x.cmp(y) != 0; }
/// double != BigInt
inline bool operator!=(double x, const BigInt& y)
{ return y.cmp(x) != 0; }

/// BigInt >= BigInt
inline bool operator>=(const BigInt& x, const BigInt& y)
{ return x.cmp(y) >= 0; }
/// BigInt >= int
inline bool operator>=(const BigInt& x, int y)
{ return x.cmp(y) >= 0; }
/// int >= BigInt
inline bool operator>=(int x, const BigInt& y)
{ return y.cmp(x) <= 0; }
/// BigInt >= unsigned int
inline bool operator>=(const BigInt& x, unsigned int y)
{ return x.cmp(y) >= 0; }
/// unsigned int >= BigInt
inline bool operator>=(unsigned int x, const BigInt& y)
{ return y.cmp(x) <= 0; }
/// BigInt >= long
inline bool operator>=(const BigInt& x, long y)
{ return x.cmp(y) >= 0; }
/// long >= BigInt
inline bool operator>=(long x, const BigInt& y)
{ return y.cmp(x) <= 0; }
/// BigInt >= unsigned long
inline bool operator>=(const BigInt& x, unsigned long y)
{ return x.cmp(y) >= 0; }
/// unsigned long >= BigInt
inline bool operator>=(unsigned long x, const BigInt& y)
{ return y.cmp(x) <= 0; }
/// BigInt >= double
inline bool operator>=(const BigInt& x, double y)
{ return x.cmp(y) >= 0; }
/// double >= BigInt
inline bool operator>=(double x, const BigInt& y)
{ return y.cmp(x) <= 0; }

/// BigInt <= BigInt
inline bool operator<=(const BigInt& x, const BigInt& y)
{ return x.cmp(y) <= 0; }
/// BigInt <= int
inline bool operator<=(const BigInt& x, int y)
{ return x.cmp(y) <= 0; }
/// int <= BigInt
inline bool operator<=(int x, const BigInt& y)
{ return y.cmp(x) >= 0; }
/// BigInt <= unsigned int
inline bool operator<=(const BigInt& x, unsigned int y)
{ return x.cmp(y) <= 0; }
/// unsigned int <= BigInt
inline bool operator<=(unsigned int x, const BigInt& y)
{ return y.cmp(x) >= 0; }
/// BigInt <= long
inline bool operator<=(const BigInt& x, long y)
{ return x.cmp(y) <= 0; }
/// long <= BigInt
inline bool operator<=(long x, const BigInt& y)
{ return y.cmp(x) >= 0; }
/// BigInt <= unsigned long
inline bool operator<=(const BigInt& x, unsigned long y)
{ return x.cmp(y) <= 0; }
/// unsigned long <= BigInt
inline bool operator<=(unsigned long x, const BigInt& y)
{ return y.cmp(x) >= 0; }
/// BigInt <= double
inline bool operator<=(const BigInt& x, double y)
{ return x.cmp(y) <= 0; }
/// double <= BigInt
inline bool operator<=(double x, const BigInt& y)
{ return y.cmp(x) >= 0; }

/// BigInt > BigInt
inline bool operator>(const BigInt& x, const BigInt& y)
{ return x.cmp(y) > 0; }
/// BigInt > int
inline bool operator>(const BigInt& x, int y)
{ return x.cmp(y) > 0; }
/// int > BigInt
inline bool operator>(int x, const BigInt& y)
{ return y.cmp(x) < 0; }
/// BigInt > unsigned int
inline bool operator>(const BigInt& x, unsigned int y)
{ return x.cmp(y) > 0; }
/// unsigned int > BigInt
inline bool operator>(unsigned int x, const BigInt& y)
{ return y.cmp(x) < 0; }
/// BigInt > long
inline bool operator>(const BigInt& x, long y)
{ return x.cmp(y) > 0; }
/// long > BigInt
inline bool operator>(long x, const BigInt& y)
{ return y.cmp(x) < 0; }
/// BigInt > unsigned long
inline bool operator>(const BigInt& x, unsigned long y)
{ return x.cmp(y) > 0; }
/// unsigned long > BigInt
inline bool operator>(unsigned long x, const BigInt& y)
{ return y.cmp(x) < 0; }
/// BigInt > double
inline bool operator>(const BigInt& x, double y)
{ return x.cmp(y) > 0; }
/// double > BigInt
inline bool operator>(double x, const BigInt& y)
{ return y.cmp(x) < 0; }

/// BigInt < BigInt
inline bool operator<(const BigInt& x, const BigInt& y)
{ return x.cmp(y) < 0; }
/// BigInt < int
inline bool operator<(const BigInt& x, int y)
{ return x.cmp(y) < 0; }
/// int < BigInt
inline bool operator<(int x, const BigInt& y)
{ return y.cmp(x) > 0; }
/// BigInt < unsigned int
inline bool operator<(const BigInt& x, unsigned int y)
{ return x.cmp(y) < 0; }
/// unsigned int < BigInt
inline bool operator<(unsigned int x, const BigInt& y)
{ return y.cmp(x) > 0; }
/// BigInt < long
inline bool operator<(const BigInt& x, long y)
{ return x.cmp(y) < 0; }
/// long < BigInt
inline bool operator<(long x, const BigInt& y)
{ return y.cmp(x) > 0; }
/// BigInt < unsigned long
inline bool operator<(const BigInt& x, unsigned long y)
{ return x.cmp(y) < 0; }
/// unsigned long < BigInt
inline bool operator<(unsigned long x, const BigInt& y)
{ return y.cmp(x) > 0; }
/// BigInt < double
inline bool operator<(const BigInt& x, double y)
{ return x.cmp(y) < 0; }
/// double < BigInt
inline bool operator<(double x, const BigInt& y)
{ return y.cmp(x) > 0; }
//@}

/// \addtogroup BigIntIostreamOperators
//@{
/// istream operator for <tt>BigInt</tt>
inline std::istream& operator>>(std::istream& is, BigInt& x)
{ return ::extract(is, x.mp()); }
/// ostream operator for <tt>BigInt</tt>
inline std::ostream& operator<<(std::ostream& os, const BigInt& x)
{ return os << x.get_str(); }
//@}

/// \addtogroup BigIntGlobalFunctions
//@{
/// read from file
inline void readFromFile(BigInt& z, std::istream& in) {
  std::string str; char c;
  in >> c;
  if (c != 'i')
    core_error("BigInt readFromFile wrong format", __FILE__, __LINE__, false);
  in >> str;
  z.set(str);
}

/// write to file
inline void writeToFile(const BigInt& z, std::ostream& out, int base=10) {
  std::string str("i");
  if (base == 2)
    str += "0b";
  else if (base == 16)
    str += "0x";
  else if (base == 8)
    str += '0';
  str += z.get_str(base);
  out << str;
}
  
/// return a gmp_randstate_t structure
inline gmp_randstate_t* getRandstate() {
  static gmp_randstate_t rstate;
  static bool initialized = false;
  if (!initialized) {
    gmp_randinit(rstate, GMP_RAND_ALG_DEFAULT, 32L);
    initialized = true;
  }
  return &rstate;
}
/// randomize function
inline BigInt randomize(const BigInt& a)
{ BigInt r; mpz_urandomm(r.mp(), *getRandstate(), a.mp()); return r; }
//@}

#ifndef CORE_DISABLE_OLDNAMES 
/// \addtogroup BigIntBackCompatiableFunctions
//@{
/// comparison
inline int cmp(const BigInt& x, const BigInt& y) { return x.cmp(y); }
/// sign 
inline int sign(const BigInt& a) { return a.sgn(); }
inline int sgn(const BigInt& a) { return a.sgn(); }
/// abs
inline BigInt abs(const BigInt& a) { BigInt r; r.abs(a); return r; }
/// neg
inline BigInt neg(const BigInt& a) { BigInt r; r.neg(a); return r; }
/// negate
inline void negate(BigInt& a) { a.neg(a); }
/// cmpabs
inline int cmpabs(const BigInt& a, const BigInt& b) { return a.cmpabs(b); }
/// longValue
inline long longValue(const BigInt& a) { return a.longValue(); }
/// ulongValue
inline unsigned long ulongValue(const BigInt& a) { return a.ulongValue(); }
/// doubleValue
inline double doubleValue(const BigInt& a) { return a.doubleValue(); }
// isEven
inline bool isEven(const BigInt& z) { return z.is_even(); }
/// isOdd
inline bool isOdd(const BigInt& z) { return z.is_odd(); }
/// get exponent of power 2
inline unsigned long getBinExpo(const BigInt& z) { return z.get_2exp(); }
/// get exponent of power k
inline void getKaryExpo(const BigInt& z, BigInt& m, int& e, unsigned long k) 
{ e = z.get_k_exp(m, k); }
/// divisible(x,y) = "x | y"
inline bool isDivisible(const BigInt& x, const BigInt& y)
{ return x.is_divisible(y) != 0; }
/// exact div
inline void divexact(BigInt& z, const BigInt& x, const BigInt& y)
{ z.divexact(x, y); }
/// exact div
inline BigInt div_exact(const BigInt& x, const BigInt& y)
{ BigInt z; divexact(z, x, y); return z; }
/// gcd
inline BigInt gcd(const BigInt& a, const BigInt& b)
{ BigInt r; r.gcd(a, b); return r; }

/// div_rem
inline void div_rem(BigInt& q, BigInt& r, const BigInt& a, const BigInt& b)
{ q.divrem(r, a, b); } 

/// void pow(Bigint& a, const Bigint& b, unsigned long c):
///   3 argument version of power, modifying the first BigInt parameter
// Chee: I propose to reverse the names:
//   --pow() should look more like MPFR's pow (taking 3 arguments)
//   --power() should look more like math functions (taking 2 arguments)
//   I will do this name switch.  June 20, 2010.
//	WARNING: This solution is not entirely consistent with the
//	pow(...) function in math.h.
inline void pow(BigInt& c, const BigInt& a, unsigned long ul) 
{ c.pow(a, ul); }

/// BigInt power(const BigInt& a, unsigned long b):
///    2 argument version of power, returning a BigInt 
inline BigInt power(const BigInt& a, unsigned long ui) 
{ BigInt r; r.pow(a, ui); return r; }

/// BigInt pow(const BigInt& a, const BigInt& b):
///	2 argument version, returning a BigInt 
//   Chee: we would like this second argument to be BigInt for 
//   compatibility even though in practice, this argument should
//    fit into a long:--
inline BigInt power(const BigInt& a, const BigInt& b) 
{ BigInt r;
  BigInt longMax = (std::numeric_limits<long>::max)();
  assert(b < longMax);	// it does not make sense to have larger exponent anyways...
  r.pow(a, b.ulongValue()); return r; }

///  BigInt power(int a, int b)
// Need this to disambiguate from other power(,) functions!
inline BigInt power(int a, const BigInt& b) 
{ return power(BigInt(a), b); }

// bit length
inline int bitLength(const BigInt& a) 
{ return a.sizeinbase(2); }
// binomial
inline BigInt binomial(unsigned long n, unsigned long k)
{ BigInt r; r.binomial(n, k); return r; }
/// floorLg -- floor of log_2(a)
/** Convention: a=0, floorLg(a) returns -1. (!!changed: return 0) 
 *  This makes sense for integer a.
 */
inline long floorLg(const BigInt& a) { return a.floorlg(); }
inline long floorlg(const BigInt& a) { return a.floorlg(); }
/// ceilLg -- ceiling of log_2(a) where a=BigInt, int or long
/** Convention: a=0, ceilLg(a) returns -1. (!!changed: return 0 now)
 *  This makes sense for integer a.
 */
inline long ceilLg(const BigInt& a) { return a.ceillg(); }
inline long ceillg(const BigInt& a) { return a.ceillg(); }
//@}
#endif

CORE_END_NAMESPACE

#endif /*__CORE_BIGINT_H__*/
