/****************************************************************************
 * BigFloat2.h -- Big Floating-point number class providing arbitrary precision
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
 * $Id: BigFloat2.h,v 1.44 2010/11/23 17:58:36 exact Exp $
 ***************************************************************************/
#ifndef __CORE_BIGFLOAT2_H__
#define __CORE_BIGFLOAT2_H__

#include <CORE/BigFloat.h>
#include <CORE/Policies.h>
#include <sstream>

CORE_BEGIN_NAMESPACE

// BUG:
// 
//  1. all fixed version and auto version have problems when one of input 
//     are output since the set_prec function will reset the input value.
//	(Problem has been fixed, Jihun)

/// \class BigFloat2
/// \brief BigFloat2 is a floating-point interval class 
class BigFloat2 {
public: // public typedefs
  typedef BigInt   ZT;
  typedef BigRat   QT;
  typedef BigFloat FT;

public: //private:
  FT   m_l;      ///<- lower bound
  FT   m_r;      ///<- upper bound
  bool m_exact;  ///<- exact flag (when it is true, m_l is the exact value
                 ///                and m_r is undefined)

public:
  /// \name constructors (auto version)
  //@{
  /// default constructor
  BigFloat2() : m_l(0), m_exact(true)
  {}
  /// copy constructor
  BigFloat2(const BigFloat2& r) : m_l(r.m_l), m_r(r.m_r), m_exact(r.m_exact)
  {}
  explicit BigFloat2(const double l, const double r) : m_l(l), m_r(r), m_exact(false) {}
  /// constructor from two BigFloat
  explicit BigFloat2(const BigFloat& l, const BigFloat& r) : m_l(l),m_r(r),m_exact(false) {}
  /// constructor for <tt>int</tt> 
  explicit BigFloat2(const int v) : m_l(v), m_exact(true) {}
  /// constructor for <tt>unsigned int</tt> 
  explicit BigFloat2(const unsigned int v) : m_l(v), m_exact(true) {}
  /// constructor for <tt>long</tt> 
  explicit BigFloat2(const long v) : m_l(v), m_exact(true) {}
  /// constructor for <tt>unsigned long</tt> 
  explicit BigFloat2(const unsigned long v) : m_l(v), m_exact(true) {}
  /// constructor for <tt>double</tt> 
  explicit BigFloat2(const double v) : m_l(v), m_exact(true) {}
  /// constructor for <tt>BigInt</tt> 
  explicit BigFloat2(const BigInt& v) : m_l(v), m_exact(true) {}
  /// constructor for <tt>BigFloat</tt> 
  explicit BigFloat2(const BigFloat& v) : m_l(v), m_exact(true) {}
#ifndef CORE_LEVEL_1_NO_WRAPPERS
  /// constructor for <tt>DoulbeWrapper</tt> 
  explicit BigFloat2(const DoubleWrapper& v) : m_l(v.doubleValue()), m_exact(true) {}
#endif
  /// generic constructor for <tt>T</TT>
//  template <typename T> BigFloat2(const T& v) : m_l(v), m_exact(true) {}
  /// construct for <tt>QT = BigRat</tt>
  explicit BigFloat2(const QT& q) 
  { set(q, DOUBLE_PREC); }
  /// construct for <tt>char*</tt> (no implicit conversion)
  explicit BigFloat2(const char* str, int base = 10) 
  { set(str, base); }
  /// construct for <tt>std::string</tt> (no implicit conversion)
  explicit BigFloat2(const std::string& str, int base = 10) 
  { set(str, base); }
  //@}

  /// \name constructors (fixed version)
  //@{
  /// generic constructor with specified precision (reduce to call set())
//  template <typename T> BigFloat2(const T& v, prec_t prec) 
//  { set(v, prec); }
  /// constructor for <tt>char*</tt>
  BigFloat2(const char* str, int base, prec_t prec) 
  { set(str, base, prec); }
  /// constructor for <tt>std::string</tt>
  BigFloat2(const std::string& str, int base, prec_t prec) 
  { set(str, base, prec); }
  /// constructor with value \f$i*2^e\f$ for <tt>char</tt>
  BigFloat2(char i, exp_t e, prec_t prec) 
  { set_2exp(i, e, prec); }
  /// constructor with value \f$i*2^e\f$ for <tt>unsigned char</tt>
  BigFloat2(unsigned char i, exp_t e, prec_t prec) 
  { set_2exp(i, e, prec); }
  /// constructor with value \f$i*2^e\f$ for <tt>short</tt>
  BigFloat2(short i, exp_t e, prec_t prec) 
  { set_2exp(i, e, prec); }
  /// constructor with value \f$i*2^e\f$ for <tt>unsigned short</tt>
  BigFloat2(unsigned short i, exp_t e, prec_t prec) 
  { set_2exp(i, e, prec); }
  /// constructor with value \f$i*2^e\f$ for <tt>int</tt>
  BigFloat2(int i, exp_t e, prec_t prec) 
  { set_2exp(i, e, prec); }
  /// constructor with value \f$i*2^e\f$ for <tt>unsigned int</tt>
  BigFloat2(unsigned int i, exp_t e, prec_t prec) 
  { set_2exp(i, e, prec); }
  /// constructor with value \f$i*2^e\f$ for <tt>long</tt>
  BigFloat2(long i, exp_t e, prec_t prec) 
  { set_2exp(i, e, prec); }
  /// constructor with value \f$i*2^e\f$ for <tt>unsigned long</tt>
  BigFloat2(unsigned long i, exp_t e, prec_t prec) 
  { set_2exp(i, e, prec); }
  //@}

public: 
  /// \name precision accessors
  //@{
  /// return current precision
  prec_t get_prec() const 
  { return is_exact()?m_l.get_prec():(std::max)(m_l.get_prec(),m_r.get_prec()); }
  /// set current precision
  void set_prec(prec_t prec)
  { m_l.set_prec(prec); if (!is_exact()) m_r.set_prec(prec); }
  //@}

  /// \name exactness accessors
  //@{
  /// return exactness
  bool is_exact() const
  { return m_exact; };
  /// set exactness
  void set_exact(bool b)
  { m_exact = b; }

  /// if it's exact, return true
  /// if not marked exact, but exact, set flag and return true;
  /// else return false;
  bool check_exactness() 
  { 
    if ( m_exact ) return true;
    if ( m_l == m_r ) { set_exact( true ); return true; }    
    return false;
  }
  //@}

  /// check if m_r >= m_l
  /// if not, swap bounds
  void validate() 
  { 
    if ( m_l > m_r ) { m_l.swap(m_r); }    
  }
  //@}


  /// \name assignment functions (raw version)
  //@{
  /// assignment function for <tt>BigFloat2</tt>
  bool r_set(const BigFloat2& x)
  { return _set_f<RawArithmeticPolicy>(x); }
  /// assignment function for <tt>char*</tt>
  bool r_set(const char* x, int base = 10)
  { return _set_str<RawArithmeticPolicy>(x, base); }
  /// assignment function for <tt>std::string</tt>
  bool r_set(const std::string& x, int base = 10)
  { return _set_str<RawArithmeticPolicy>(x.c_str(), base); }
  /// generic assignment function for <tt>T</tt>
  template <typename T> bool r_set(const T& x)
  { return _set<RawArithmeticPolicy, T>(x); }
  /// set value to be \f$i*2^e\f$ for <tt>int</tt>
  bool r_set_2exp(int i, exp_t e)
  { return _set_2exp_si<RawArithmeticPolicy>(long(i), e); }
  /// set value to be \f$i*2^e\f$ for <tt>unsigned int</tt>
  bool r_set_2exp(unsigned int i, exp_t e)
  { return _set_2exp_ui<RawArithmeticPolicy>((unsigned long)(i), e); }
  /// set value to be \f$i*2^e\f$ for <tt>long</tt>
  bool r_set_2exp(long i, exp_t e)
  { return _set_2exp_si<RawArithmeticPolicy>(i, e); }
  /// set value to be \f$i*2^e\f$ for <tt>unsigned long</tt>
  bool r_set_2exp(unsigned long i, exp_t e)
  { return _set_2exp_ui<RawArithmeticPolicy>(i, e); }
  //@}

  /// \name assignment functions (fixed version)
  //@{
  /// assignment function for <tt>BigFloat2</tt>
  bool set(const BigFloat2& x, prec_t prec)
  { return _set_f<FixedArithmeticPolicy>(x, prec); }
  /// assignment function for <tt>char*</tt>
  bool set(const char* x, int base, prec_t prec)
  { return _set_str<FixedArithmeticPolicy>(x, base, prec); }
  /// assignment function for <tt>std::string</tt>
  bool set(const std::string& x, int base, prec_t prec)
  { return _set_str<FixedArithmeticPolicy>(x.c_str(), base, prec); }
  /// assignment function for <tt>T</tt>
  template <typename T> bool set(const T& x, prec_t prec)
  { return _set<FixedArithmeticPolicy, T>(x, prec); }

  /// set value to be \f$i*2^e\f$ for <tt>int</tt>
  bool set_2exp(int i, exp_t e, prec_t prec)
  { return _set_2exp_si<FixedArithmeticPolicy>(long(i), e, prec); }
  /// set value to be \f$i*2^e\f$ for <tt>unsigned int</tt>
  bool set_2exp(unsigned int i, exp_t e, prec_t prec)
  { return _set_2exp_ui<FixedArithmeticPolicy>((unsigned long)(i), e, prec); }
  /// set value to be \f$i*2^e\f$ for <tt>long</tt>
  bool set_2exp(long i, exp_t e, prec_t prec)
  { return _set_2exp_si<FixedArithmeticPolicy>(i, e, prec); }
  /// set value to be \f$i*2^e\f$ for <tt>unsigned long</tt>
  bool set_2exp(unsigned long i, exp_t e, prec_t prec)
  { return _set_2exp_ui<FixedArithmeticPolicy>(i, e, prec); }
  //@}

  /// \name assignment functions (auto version)
  //@{
  /// assignment function for <tt>BigFloat2</tt>
  bool set(const BigFloat2& x)
  { return _set_f<AutoArithmeticPolicy>(x); }
  /// assignment function for <tt>char*</tt>
  bool set(const char* x, int base = 10)
  { return _set_str<AutoArithmeticPolicy>(x, base); }
  /// assignment function for <tt>std::string</tt>
  bool set(const std::string& x, int base = 10)
  { return _set_str<AutoArithmeticPolicy>(x.c_str(), base); }
  /// assignment function for <tt>T</tt>
  template <typename T> bool set(const T& x)
  { return _set<AutoArithmeticPolicy, T>(x); }

  /// set value to be \f$i*2^e\f$ for <tt>int</tt>
  bool set_2exp(int i, exp_t e)
  { return _set_2exp_si<AutoArithmeticPolicy>(long(i), e); }
  /// set value to be \f$i*2^e\f$ for <tt>unsigned int</tt>
  bool set_2exp(unsigned int i, exp_t e)
  { return _set_2exp_ui<AutoArithmeticPolicy>((unsigned long)(i), e); }
  /// set value to be \f$i*2^e\f$ for <tt>long</tt>
  bool set_2exp(long i, exp_t e)
  { return _set_2exp_si<AutoArithmeticPolicy>(i, e); }
  /// set value to be \f$i*2^e\f$ for <tt>unsigned long</tt>
  bool set_2exp(unsigned long i, exp_t e)
  { return _set_2exp_ui<AutoArithmeticPolicy>(i, e); }
  //@}

  /// \name arithmetic functions (raw version)
  //@{
  /// negation for <tt>BigFloat2</tt>
  bool r_neg(const BigFloat2& x)
  { return _neg_f<RawArithmeticPolicy>(x); }
  /// negation for <tt>T</tt>
  template <typename T> bool r_neg(const T& x)
  { return _neg<RawArithmeticPolicy, T>(x); }

  /// square root for <tt>BigFloat2</tt>
  bool r_sqrt(const BigFloat2& x)
  { return _sqrt_f<RawArithmeticPolicy>(x); }
  /// square root for <tt>T</tt>
  template <typename T> bool r_sqrt(const T& x)
  { return _sqrt<RawArithmeticPolicy, T>(x); }

  /// cubic root for <tt>BigFloat2</tt>
  bool r_cbrt(const BigFloat2& x)
  { return _cbrt_f<RawArithmeticPolicy>(x); }
  // /// cubic root for <tt>T</tt>
  // template <typename T> bool r_cbrt(const T& x)
  // { return _cbrt<RawArithmeticPolicy, T>(x); }

  /// k-th root for <tt>BigFloat2</tt>
  bool r_root(const BigFloat2& x, unsigned long k)
  { return _root_f<RawArithmeticPolicy>(x, k); }
  // /// k-th root for <tt>T</tt>
  // template <typename T> bool r_root(const T& x, unsigned long k)
  // { return _root<RawArithmeticPolicy, T>(x, k); }

  /// sine for <tt>BigFloat2</tt>
  bool r_sin(const BigFloat2& x)
  { return _sin_f<RawArithmeticPolicy>(x); }
  /// sine for <tt>T</tt>
  template <typename T> bool r_sin(const T& x)
  { return _sin<RawArithmeticPolicy, T>(x); }

  /// cosine for <tt>BigFloat2</tt>
  bool r_cos(const BigFloat2& x)
  { return _cos_f<RawArithmeticPolicy>(x); }
  /// cosine for <tt>T</tt>
  template <typename T> bool r_cos(const T& x)
  { return _cos<RawArithmeticPolicy, T>(x); }

  /// tangent for <tt>BigFloat2</tt>
  bool r_tan(const BigFloat2& x)
  { return _tan_f<RawArithmeticPolicy>(x); }
  /// tangent for <tt>T</tt>
  template <typename T> bool r_tan(const T& x)
  { return _tan<RawArithmeticPolicy, T>(x); }

  /// cotangent for <tt>BigFloat2</tt>
  bool r_cot(const BigFloat2& x)
  { return _cot_f<RawArithmeticPolicy>(x); }
  /// cotangent for <tt>T</tt>
  template <typename T> bool r_cot(const T& x)
  { return _cot<RawArithmeticPolicy, T>(x); }

  /// arcsine for <tt>BigFloat2</tt>
  bool r_asin(const BigFloat2& x)
  { return _asin_f<RawArithmeticPolicy>(x); }
  /// arcsine for <tt>T</tt>
  template <typename T> bool r_asin(const T& x)
  { return _asin<RawArithmeticPolicy, T>(x); }

  /// arccosine for <tt>BigFloat2</tt>
  bool r_acos(const BigFloat2& x)
  { return _acos_f<RawArithmeticPolicy>(x); }
  /// sin for <tt>T</tt>
  template <typename T> bool r_acos(const T& x)
  { return _acos<RawArithmeticPolicy, T>(x); }

  /// arctangent for <tt>BigFloat2</tt>
  bool r_atan(const BigFloat2& x)
  { return _atan_f<RawArithmeticPolicy>(x); }
  /// arctangent for <tt>T</tt>
  template <typename T> bool r_atan(const T& x)
  { return _atan<RawArithmeticPolicy, T>(x); }

  /// log for <tt>BigFloat2</tt>
  bool r_log(const BigFloat2& x)
  { return _log_f<RawArithmeticPolicy>(x); }
  /// sin for <tt>T</tt>
  template <typename T> bool r_log(const T& x)
  { return _log<RawArithmeticPolicy, T>(x); }  
  
  /// log2 for <tt>BigFloat2</tt>
  bool r_log2(const BigFloat2& x)
  { return _log2_f<RawArithmeticPolicy>(x); }
  /// sin for <tt>T</tt>
  template <typename T> bool r_log2(const T& x)
  { return _log2<RawArithmeticPolicy, T>(x); }

  /// log10 for <tt>BigFloat2</tt>
  bool r_log10(const BigFloat2& x)
  { return _log10_f<RawArithmeticPolicy>(x); }
  /// sin for <tt>T</tt>
  template <typename T> bool r_log10(const T& x)
  { return _log10<RawArithmeticPolicy, T>(x); }

  /// exponent of e for <tt>BigFloat2</tt>
  bool r_exp(const BigFloat2& x)
  { return _exp_f<RawArithmeticPolicy>(x); }
  /// exponent for <tt>T</tt>
  template <typename T> bool r_exp(const T& x)
  { return _exp<RawArithmeticPolicy, T>(x); }

  /// exponent of 2 for <tt>BigFloat2</tt>
  bool r_exp2(const BigFloat2& x)
  { return _exp2_f<RawArithmeticPolicy>(x); }
  /// exponent for <tt>T</tt>
  template <typename T> bool r_exp2(const T& x)
  { return _exp2<RawArithmeticPolicy, T>(x); }

  /// exponent of 10 for <tt>BigFloat2</tt>
  bool r_exp10(const BigFloat2& x)
  { return _exp10_f<RawArithmeticPolicy>(x); }
  /// exponent for <tt>T</tt>
  template <typename T> bool r_exp10(const T& x)
  { return _exp10<RawArithmeticPolicy, T>(x); }

  /// addition/subtraction for <tt>BigFloat2</tt>
  bool r_addsub(const BigFloat2& x, const BigFloat2& y, bool is_add)
  { return is_add ? r_add(x, y) : r_sub(x, y); }
  /// addition/subtraction for <tt>BigFloat2, T</tt>
  template <typename T>
  bool r_addsub(const BigFloat2& x, const T& y, bool is_add)
  { return is_add ? r_add(x, y) : r_sub(x, y); }
  /// addition/subtraction for <tt>T, BigFloat2</tt>
  template <typename T>
  bool r_addsub(const T& x, const BigFloat2& y, bool is_add)
  { return is_add ? r_add(x, y) : r_sub(x, y); }

  /// addition for <tt>BigFloat2+BigFloat2</tt>
  bool r_add(const BigFloat2& x, const BigFloat2& y)
  { return _add_f<RawArithmeticPolicy>(x, y); }
  /// addition for <tt>BigFloat2+T</tt>
  template <typename T> bool r_add(const BigFloat2& x, const T& y)
  { return _add<RawArithmeticPolicy, T>(x, y); }
  /// addition for <tt>T+BigFloat2</tt>
  template <typename T> bool r_add(const T& x, const BigFloat2& y)
  { return _add<RawArithmeticPolicy, T>(x, y); }

  /// subtraction for <tt>BigFloat2-BigFloat2</tt>
  bool r_sub(const BigFloat2& x, const BigFloat2& y)
  { return _sub_f<RawArithmeticPolicy>(x, y); }
  /// subtraction for <tt>BigFloat2-T</tt>
  template <typename T> bool r_sub(const BigFloat2& x, const T& y)
  { return _sub<RawArithmeticPolicy, T>(x, y); }
  /// subtraction for <tt>T-BigFloat2</tt>
  template <typename T> bool r_sub(const T& x, const BigFloat2& y)
  { return _sub<RawArithmeticPolicy, T>(x, y); }

  /// multiplication for <tt>BigFloat2*BigFloat2</tt>
  bool r_mul(const BigFloat2& x, const BigFloat2& y)
  { return _mul_f<RawArithmeticPolicy>(x, y); }
  /// multiplication for <tt>BigFloat2*T</tt>
  template <typename T> bool r_mul(const BigFloat2& x, const T& y)
  { return _mul<RawArithmeticPolicy, T>(x, y); }
  /// multiplication for <tt>T*BigFloat2</tt>
  template <typename T> bool r_mul(const T& x, const BigFloat2& y)
  { return _mul<RawArithmeticPolicy, T>(x, y); }

  /// division for <tt>BigFloat2/BigFloat2</tt>
  bool r_div(const BigFloat2& x, const BigFloat2& y)
  { return _div_f<RawArithmeticPolicy>(x, y); }
  /// division for <tt>BigFloat2/T</tt>
  template <typename T> bool r_div(const BigFloat2& x, const T& y)
  { return _div<RawArithmeticPolicy, T>(x, y); }
  /// division for <tt>T/BigFloat2</tt>
  template <typename T> bool r_div(const T& x, const BigFloat2& y)
  { return _div<RawArithmeticPolicy, T>(x, y); }
  //@}

  /// \name arithmetic functions (fixed version)
  //@{
  /// negation for <tt>BigFloat2</tt>
  bool neg(const BigFloat2& x, prec_t prec)
  { return _neg_f<FixedArithmeticPolicy>(x, prec); }
  /// negation for <tt>T</tt>
  template <typename T> bool neg(const T& x, prec_t prec)
  { return _neg<FixedArithmeticPolicy, T>(x, prec); }

  /// square root for <tt>BigFloat2</tt>
  bool sqrt(const BigFloat2& x, prec_t prec = getDefaultBFradicalPrec())
  { return _sqrt_f<FixedArithmeticPolicy>(x, prec); }
  /// square root for <tt>T</tt>
  template <typename T> bool sqrt(const T& x, prec_t prec = getDefaultBFradicalPrec())
  { return _sqrt<FixedArithmeticPolicy, T>(x, prec); }

  /// cubic root for <tt>BigFloat2</tt>
  bool cbrt(const BigFloat2& x, prec_t prec = getDefaultBFradicalPrec())
  { return _cbrt_f<FixedArithmeticPolicy>(x, prec); }
  // /// cubic root for <tt>T</tt>
  // template <typename T> bool cbrt(const T& x, prec_t prec = getDefaultBFradicalPrec())
  // { return _cbrt<FixedArithmeticPolicy, T>(x, prec); }

  /// k-th root for <tt>BigFloat2</tt>
  bool root(const BigFloat2& x, unsigned long k, prec_t prec = getDefaultBFradicalPrec())
  { return _root_f<FixedArithmeticPolicy>(x, k, prec); }
  // /// k-th root for <tt>T</tt>
  // template <typename T> 
  // bool root(const T& x, unsigned long k, prec_t prec = getDefaultBFradicalPrec())
  // { return _root<FixedArithmeticPolicy, T>(x, k, prec); }

  /// sine for <tt>BigFloat2</tt>
  bool sin(const BigFloat2& x, prec_t prec)
  { return _sin_f<FixedArithmeticPolicy>(x, prec); }
  /// sine for <tt>T</tt>
  template <typename T> bool sin(const T& x, prec_t prec)
  { return _sin<FixedArithmeticPolicy, T>(x, prec); }

  /// cosine for <tt>BigFloat2</tt>
  bool cos(const BigFloat2& x, prec_t prec)
  { return _cos_f<FixedArithmeticPolicy>(x, prec); }
  /// cose for <tt>T</tt>
  template <typename T> bool cos(const T& x, prec_t prec)
  { return _cos<FixedArithmeticPolicy, T>(x, prec); }

  /// tane for <tt>BigFloat2</tt>
  bool tan(const BigFloat2& x, prec_t prec)
  { return _tan_f<FixedArithmeticPolicy>(x, prec); }
  /// tane for <tt>T</tt>
  template <typename T> bool tan(const T& x, prec_t prec)
  { return _tan<FixedArithmeticPolicy, T>(x, prec); }

  /// coe for <tt>BigFloat2</tt>
  bool cot(const BigFloat2& x, prec_t prec)
  { return _cot_f<FixedArithmeticPolicy>(x, prec); }
  /// cote for <tt>T</tt>
  template <typename T> bool cot(const T& x, prec_t prec)
  { return _cot<FixedArithmeticPolicy, T>(x, prec); }

  /// arcsine for <tt>BigFloat2</tt>
  bool asin(const BigFloat2& x, prec_t prec)
  { return _asin_f<FixedArithmeticPolicy>(x, prec); }
  /// arcsine for <tt>T</tt>
  template <typename T> bool asin(const T& x, prec_t prec)
  { return _asin<FixedArithmeticPolicy, T>(x, prec); }

  /// arcosine for <tt>BigFloat2</tt>
  bool acos(const BigFloat2& x, prec_t prec)
  { return _acos_f<FixedArithmeticPolicy>(x, prec); }
  /// arccosine for <tt>T</tt>
  template <typename T> bool acos(const T& x, prec_t prec)
  { return _acos<FixedArithmeticPolicy, T>(x, prec); }

  /// arctangent for <tt>BigFloat2</tt>
  bool atan(const BigFloat2& x, prec_t prec)
  { return _atan_f<FixedArithmeticPolicy>(x, prec); }
  /// arctangent for <tt>T</tt>
  template <typename T> bool atan(const T& x, prec_t prec)
  { return _atan<FixedArithmeticPolicy, T>(x, prec); }

  /// log for <tt>BigFloat2</tt>
  bool log(const BigFloat2& x, prec_t prec)
  { return _log_f<FixedArithmeticPolicy>(x, prec); }
  /// log2 for <tt>T</tt>
  template <typename T> bool log(const T& x, prec_t prec)
  { return _log<FixedArithmeticPolicy, T>(x, prec); }  
  
  /// log2 for <tt>BigFloat2</tt>
  bool log2(const BigFloat2& x, prec_t prec)
  { return _log2_f<FixedArithmeticPolicy>(x, prec); }
  /// log2 for <tt>T</tt>
  template <typename T> bool log2(const T& x, prec_t prec)
  { return _log2<FixedArithmeticPolicy, T>(x, prec); }

  /// log10 for <tt>BigFloat2</tt>
  bool log10(const BigFloat2& x, prec_t prec)
  { return _log10_f<FixedArithmeticPolicy>(x, prec); }
  /// log2 for <tt>T</tt>
  template <typename T> bool log10(const T& x, prec_t prec)
  { return _log10<FixedArithmeticPolicy, T>(x, prec); }  
  
  /// exponent for <tt>BigFloat2</tt>
  bool exp(const BigFloat2& x, prec_t prec)
  { return _exp_f<FixedArithmeticPolicy>(x, prec); }
  /// exponent for <tt>T</tt>
  template <typename T> bool exp(const T& x, prec_t prec)
  { return _exp<FixedArithmeticPolicy, T>(x, prec); }

  /// exponent of 2 for <tt>BigFloat2</tt>
  bool exp2(const BigFloat2& x, prec_t prec)
  { return _exp2_f<FixedArithmeticPolicy>(x, prec); }
  /// exponent for <tt>T</tt>
  template <typename T> bool exp2(const T& x, prec_t prec)
  { return _exp2<FixedArithmeticPolicy, T>(x, prec); }  
  
  /// exponent of 10 for <tt>BigFloat2</tt>
  bool exp10(const BigFloat2& x, prec_t prec)
  { return _exp10_f<FixedArithmeticPolicy>(x, prec); }
  /// exponent for <tt>T</tt>
  template <typename T> bool exp10(const T& x, prec_t prec)
  { return _exp10<FixedArithmeticPolicy, T>(x, prec); }  

  /// addition/subtraction for <tt>BigFloat2</tt>
  bool addsub(const BigFloat2& x, const BigFloat2& y, prec_t prec, bool is_add)
  { return is_add ? add(x, y, prec) : sub(x, y, prec); }
  /// addition/subtraction for <tt>BigFloat2, T</tt>
  template <typename T>
  bool addsub(const BigFloat2& x, const T& y, prec_t prec, bool is_add)
  { return is_add ? add(x, y, prec) : sub(x, y, prec); }
  /// addition/subtraction for <tt>T, BigFloat2</tt>
  template <typename T>
  bool addsub(const T& x, const BigFloat2& y, prec_t prec, bool is_add)
  { return is_add ? add(x, y, prec) : sub(x, y, prec); }

  /// addition for <tt>BigFloat2+BigFloat2</tt>
  bool add(const BigFloat2& x, const BigFloat2& y, prec_t prec)
  { return _add_f<FixedArithmeticPolicy>(x, y, prec); }
  /// addition for <tt>BigFloat2+T</tt>
  template <typename T> bool add(const BigFloat2& x, const T& y, prec_t prec)
  { return _add<FixedArithmeticPolicy, T>(x, y, prec); }
  /// addition for <tt>T+BigFloat2</tt>
  template <typename T> bool add(const T& x, const BigFloat2& y, prec_t prec)
  { return _add<FixedArithmeticPolicy, T>(x, y, prec); }

  /// subtraction for <tt>BigFloat2-BigFloat2</tt>
  bool sub(const BigFloat2& x, const BigFloat2& y, prec_t prec)
  { return _sub_f<FixedArithmeticPolicy>(x, y, prec); }
  /// subtraction for <tt>BigFloat2-T</tt>
  template <typename T> bool sub(const BigFloat2& x, const T& y, prec_t prec)
  { return _sub<FixedArithmeticPolicy, T>(x, y, prec); }
  /// subtraction for <tt>T-BigFloat2</tt>
  template <typename T> bool sub(const T& x, const BigFloat2& y, prec_t prec)
  { return _sub<FixedArithmeticPolicy, T>(x, y, prec); }

  /// multiplication for <tt>BigFloat2*BigFloat2</tt>
  bool mul(const BigFloat2& x, const BigFloat2& y, prec_t prec)
  { return _mul_f<FixedArithmeticPolicy>(x, y, prec); }
  /// multiplication for <tt>BigFloat2*T</tt>
  template <typename T> bool mul(const BigFloat2& x, const T& y, prec_t prec)
  { return _mul<FixedArithmeticPolicy, T>(x, y, prec); }
  /// multiplication for <tt>T*BigFloat2</tt>
  template <typename T> bool mul(const T& x, const BigFloat2& y, prec_t prec)
  { return _mul<FixedArithmeticPolicy, T>(x, y, prec); }

  /// division for <tt>BigFloat2/BigFloat2</tt>
  bool div(const BigFloat2& x, const BigFloat2& y, prec_t prec = getDefaultBFdivPrec())
  { return _div_f<FixedArithmeticPolicy>(x, y, prec); }
  /// division for <tt>BigFloat2/T</tt>
  template <typename T> 
  bool div(const BigFloat2& x, const T& y, prec_t prec = getDefaultBFdivPrec())
  { return _div<FixedArithmeticPolicy, T>(x, y, prec); }
  /// division for <tt>T/BigFloat2</tt>
  template <typename T> 
  bool div(const T& x, const BigFloat2& y, prec_t prec = getDefaultBFdivPrec())
  { return _div<FixedArithmeticPolicy, T>(x, y, prec); }
  //@}

  /// \name arithmetic functions (auto version)
  /// NOTE: no auto version for sqrt, cbrt, root, div
  //@{
  /// negation for <tt>BigFloat2</tt>
  bool neg(const BigFloat2& x)
  { return _neg_f<AutoArithmeticPolicy>(x); }
  /// negation for <tt>T</tt>
  template <typename T> bool neg(const T& x)
  { return _neg<AutoArithmeticPolicy, T>(x); }

  /// addition/subtraction for <tt>BigFloat2</tt>
  bool addsub(const BigFloat2& x, const BigFloat2& y, bool is_add)
  { return is_add ? add(x, y) : sub(x, y); }
  /// addition/subtraction for <tt>BigFloat2, T</tt>
  template <typename T>
  bool addsub(const BigFloat2& x, const T& y, bool is_add)
  { return is_add ? add(x, y) : sub(x, y); }
  /// addition/subtraction for <tt>T, BigFloat2</tt>
  template <typename T>
  bool addsub(const T& x, const BigFloat2& y, bool is_add)
  { return is_add ? add(x, y) : sub(x, y); }

  /// addition for <tt>BigFloat2+BigFloat2</tt>
  bool add(const BigFloat2& x, const BigFloat2& y)
  { return _add_f<AutoArithmeticPolicy>(x, y); }
  /// addition for <tt>BigFloat2+T</tt>
  template <typename T> bool add(const BigFloat2& x, const T& y)
  { return _add<AutoArithmeticPolicy, T>(x, y); }
  /// addition for <tt>T+BigFloat2</tt>
  template <typename T> bool add(const T& x, const BigFloat2& y)
  { return _add<AutoArithmeticPolicy, T>(x, y); }

  /// subtraction for <tt>BigFloat2-BigFloat2</tt>
  bool sub(const BigFloat2& x, const BigFloat2& y)
  { return _sub_f<AutoArithmeticPolicy>(x, y); }
  /// subtraction for <tt>BigFloat2-T</tt>
  template <typename T> bool sub(const BigFloat2& x, const T& y)
  { return _sub<AutoArithmeticPolicy, T>(x, y); }
  /// subtraction for <tt>T-BigFloat2</tt>
  template <typename T> bool sub(const T& x, const BigFloat2& y)
  { return _sub<AutoArithmeticPolicy, T>(x, y); }

  /// multiplication for <tt>BigFloat2*BigFloat2</tt>
  bool mul(const BigFloat2& x, const BigFloat2& y)
  { return _mul_f<AutoArithmeticPolicy>(x, y); }
  /// multiplication for <tt>BigFloat2*T</tt>
  template <typename T> bool mul(const BigFloat2& x, const T& y)
  { return _mul<AutoArithmeticPolicy, T>(x, y); }
  /// multiplication for <tt>T*BigFloat2</tt>
  template <typename T> bool mul(const T& x, const BigFloat2& y)
  { return _mul<AutoArithmeticPolicy, T>(x, y); }
  //@}

  /// \special functions
  //@{
  /// pi function
  // May 2010 (Jihun/Chee): MPFR bug???
  // In this code -- the rounding up and down always gives the same value, but
  // this should NEVER be true!   So we will need to fudge around MPFR,
  // to get a correct non-zero interval.
  bool pi(prec_t prec) {
  /* The following should be the correct code, if MPFR were correctly rounded:
    m_l.pi(prec, GMP_RNDD);
    m_r.pi(prec, GMP_RNDU);
  */
  // Here is the fudged code, and it only depends on MPFR giving us
  // a value that is within the specified precision, regardless of rounding!
	  m_l.pi(prec+1, GMP_RNDD);
	  m_l.nextbelow(); // use nice builtin function to get previous BigFloat value
	  m_r.pi(prec+1, GMP_RNDU);
	  m_r.nextabove(); // use nice builtin function to get next BigFloat value
	  set_exact(false);
    return false;
  }
  /// e function
  // May 2010 (Jihun/Chee): same problem from MPFR as in pi() above
  bool e(prec_t prec) {
  /* this should be the correct code, if MPFR were correctly rounded:
     m_l.e(prec, GMP_RNDD);  
     m_l.e(prec, GMP_RNDU);
  */
  // Here is the fudged code to get around MPFR:
	m_l.e(prec+1, GMP_RNDD);
	m_l.nextbelow(); 
	m_r.e(prec+1, GMP_RNDU);
	m_r.nextabove(); 
	set_exact(false);
    return false;
  }
  //@}

public:
  /// \name conversion functions
  //@{
  /// return double value
  double get_d() const {
    return getCenter().get_d();
  }
  /// return the string representation
  ///    get_str( n, b) returns a string of length n in base b.
  //      
  std::string get_str(size_t digits = 0, int base = 10) const {
    if(is_exact())
      return get_f().get_str(digits,base);
    else {
      long valprec;
      FT f = get_f(valprec);
      if(valprec > 0) 
        return f.get_str(bits2digits(f.get_prec()+1),base);
      else {
        std::ostringstream oss;
        oss << "0+/-10e" << bits2digits(-valprec);
        return oss.str();
      }  
    }
/*
    BigFloat bf = getCenter();
    BigFloat radi = abs_radius();

std::cerr<< "diam and radius = "
	<< abs_diam()
        << std::endl
        << div2(abs_diam())
        << std::endl;

    BigFloat logValue = -log_10(radi,100,GMP_RNDU);
    int errDigits = logValue.intValue();

std::cerr << "bf=" << bf
	<< std::endl
	<< "radi=" << radi
	<< std::endl
	<< "logValue=" << logValue
	<< std::endl
	<< "errDigits=" << errDigits
	<< std::endl;

    if (errDigits > log_10(bf,100,GMP_RNDU).intValue()) { // the error is larger than the magnitude of center value
      std::ostringstream oss;
      oss << "0+/-10e" << -errDigits;
      return oss.str();
    } else
      return get_f().get_str(errDigits,10);
*/
  }
  /// return <tt>BigInt</tt> value
  ZT get_z() const
  { return m_l.get_z(); }
  /// return <tt>BigRat</tt> value
  QT get_q() const {
    return m_l.get_q();
  }
  /// return <tt>BigFloat</tt> value
  FT get_f() const {
    if (is_exact()) 
      return m_l;
    else {
      // get error precision
      long bits = abs_diam().get_exp();
      
      FT result(getCenter());

      // round up to validate precision
      long valprec = result.get_exp() - bits;

      // jihun 2006: in case of m_r.get_exp() > m_l.get_exp(),
      // e.g. m_l = 0.1 * 2^-426, m_r = 0.1 * 2^-423, bits = -423
      // I dont know how to deal with this.
      // We need to check if m_l and m_r has the same exponent value
      if (valprec > 0)
        result.prec_round((std::max)(valprec, (long)2));
      else if (valprec <= 0) {
	core_error("This BigFloat2 has no implicit error representation",
	 __FILE__, __LINE__, false); // a warning, not necessarily an error!
      }

        // COMMENT: A BigFloat say 0.1234 has an implicit error of 0.0001.
	//          So we say 0.1234 is an "implicit error representation"
	//          of the BigFloat2 [0.1233,0.1235].
	//          However, this concept is an I/O issue and should not concern us here.
      return result;
    }
  }
   /// return <tt>BigFloat</tt> value with a valprec
  FT get_f(long& valprec) const {
    if (is_exact()) 
      return m_l;
    else {
      // get error precision
      long bits = abs_diam().get_exp();
      
      FT result(getCenter());

      // round up to validate precision
      valprec = result.get_exp() - bits;

      // jihun 2006: in case of m_r.get_exp() > m_l.get_exp(),
      // e.g. m_l = 0.1 * 2^-426, m_r = 0.1 * 2^-423, bits = -423
      // is this a problem?
      // Answer: it is not a problem.
      if (valprec > 0)
        result.prec_round((std::max)(valprec, (long)2));
      else if (valprec < 0) {
	core_error("This BigFloat2 has no implicit error representation",
	 __FILE__, __LINE__, false); // a warning, not necessarily an error!
      }
        // COMMENT: A BigFloat say 0.1234 has an implicit error of 0.0001.
	//          So we say 0.1234 is an "implicit error representation"
	//          of the BigFloat2 [0.1233,0.1235].  However, this concept
	//          is an I/O issue and should not concern us here.
      return result;
    }
  }
  //@}

  /// \name miscellaneous functions
  //@{
  /// swap function
  void swap(BigFloat2 &rhs)
  { m_l.swap(rhs.m_l); m_r.swap(rhs.m_r); std::swap(m_exact, rhs.m_exact); }

  /// return absolute diameter
  FT abs_diam() const
  { FT diam; abs_diam(diam); return diam; }

  /// get absolute diameter, return 0 if exact
  int abs_diam(FT& diam) const {
    if (is_exact()) {
      return diam.set(0);
    } else {
      diam.set_prec(get_prec()); return diam.sub(m_r, m_l); 
    }
  }
  FT abs_radius() const {
    return div2(abs_diam());
  }
  void makeCeilExact() { 
    if (is_exact()) return;  
    m_l = m_r; m_r = 0; m_exact = true;
  }
  void makeFloorExact() { 
    if (is_exact()) return;  
    m_r = 0; m_exact = true;
  }
  void makeExact() { 
    set_exact(true);
  }
  FT get_max() const {
    if (is_exact()) return m_l;
    if (has_sign())
      if (sgn() > 0) return m_r;
      else return m_l;
    else
      if (abs(m_l) > abs(m_r)) return m_l;
      else return m_r;
  }
  FT get_min() const {
    if (is_exact()) return m_l;
    if (has_sign())
      if (sgn() > 0) return m_l;
      else return m_r;
    else
      if (abs(m_l) > abs(m_r)) return m_r;
      else return m_l;
  }
  BigFloat getLeft() const {
    return m_l;
  }
  BigFloat getRight() const {
    if (is_exact())
      return m_l;
    else
      return m_r;
  }
  BigFloat getCenter() const {
    if (is_exact()) return m_l;
    BigFloat ret = div2(m_l + m_r);
    return ret;
  }
  BigFloat centerize() {
    if (is_exact()) return m_l;  
    m_l.add(m_l, m_r);
    m_l.div_2exp(1);
    m_r = 0; m_exact = true;
    return m_l;
  }
  /// check whether contains zero
  bool has_zero() const
  { return is_exact() ? (m_l.sgn()==0) : (m_l.sgn()<=0 && m_r.sgn()>=0); } 
  bool isZeroIn() const
  { return has_zero(); }
  bool has_sign() const
  { return !has_zero(); }
  /// return sign
  int sgn() const {
    if (is_exact()) {
      return m_l.sgn();
    } else if (m_l.sgn() > 0) {
      return 1;
    } else if (m_r.sgn() < 0) {
      return -1;
    } else {
#ifndef NDEBUG
      std::cerr << "BigFloat2 Warning: cannot get correct sign!" << std::endl;
      std::cerr << "m_l:m_r=" << m_l << ":" << m_r << std::endl;
#endif
      return 0;
    }
  }
  /// return upper bound of MSB
  long uMSB() const
  { return is_exact() ? m_l.uMSB() : m_r.uMSB(); }
  /// return lower bound of MSB
  long lMSB() const
  { return m_l.lMSB(); }
  //@}

  /// \name helper functions
  //@{
  /// set to be zero
  void set_zero()
  { m_l.set(0); m_exact = true; }
  /// set to be \f$[-\infty, +\infty]\f$
  void set_inf()
  { m_l.set_neg_inf(); m_r.set_pos_inf(); m_exact = false; }
  /// set to be NaN
  void set_nan()
  { m_l.set_nan(); m_exact = true; }
  bool is_integer() const
  { return is_exact() && m_l.is_integer(); }
  //@}

public: // C++ style operators
  /// \name unary, increment, decrement operators
  //@{ 
  /// unary plus operator
  BigFloat2 operator+() const
  { return BigFloat2(*this); }
  /// unary negation operator
  BigFloat2 operator-() const
  { BigFloat2 r; r.neg(*this); return r; }
  /// prefix increment operator
  BigFloat2& operator++()
  { *this += 1; return *this; }
  /// postfix increment operator
  BigFloat2 operator++(int)
  { BigFloat2 r(*this); ++(*this); return r; }
  /// prefix decrement operator
  BigFloat2& operator--()
  { *this -= 1; return *this; }
  /// postfix decrement operator
  BigFloat2 operator--(int)
  { BigFloat2 r(*this); --(*this); return r; }
  //@}

  /// \name assignment and compound assignment operators (call auto version)
  //@{
  /// assignment operator for <tt>BigFloat2</tt>
  BigFloat2& operator=(const BigFloat2& rhs)
  { if (&rhs != this) set(rhs); return *this; }
  /// generic assignment operator for <tt>T</tt>
  template <typename T> BigFloat2& operator=(const T& rhs)
  { set(rhs); return *this; }
  /// generic compound assignment operator <tt>+=</tt>
  template <typename T> BigFloat2& operator+=(const T& rhs) 
  { BigFloat2 t; t.add(*this, rhs); t.swap(*this); return *this; }
  /// generic compound assignment operator <tt>-=</tt>
  template <typename T> BigFloat2& operator-=(const T& rhs) 
  { BigFloat2 t; t.sub(*this, rhs); t.swap(*this); return *this; }
  /// generic compound assignment operator <tt>*=</tt>
  template <typename T> BigFloat2& operator*=(const T& rhs) 
  { BigFloat2 t; t.mul(*this, rhs); t.swap(*this); return *this; }
  /// generic compound assignment operator <tt>/=</tt>
  template <typename T> BigFloat2& operator/=(const T& rhs) 
  { BigFloat2 t; t.div(*this, rhs); t.swap(*this); return *this; }
  //@}

#ifdef CORE_OLDNAMES
  /// \name back-compatiable functions
  //@{
  /// Has Exact Division
  static bool hasExactDivision() { return false; }
  /// sign function
  /** \note This is only the sign of the mantissa, it can be taken to be
      the sign of the BigFloat2 only if !(isZeroIn()). */
  int sign() const { return sgn(); }
  /// check whether contains zero
  /** \return true if contains zero, otherwise false */
  bool isZeroIn() const { return has_zero(); }
  /// set value from <tt>const char*</tt> (base = 10)
  void fromString(const char* s) { set(s, 10); }
  /// convert to <tt>std::string</tt> (base = 10)
  std::string toString() const { return get_str(); }
  /// convert to <tt>std::string</tt> (base = 10)
  std::string str() const { return toString(); }
  /// return int value
  int intValue() const { return static_cast<int>(get_d()); }
  /// return long value
  int longValue() const { return static_cast<long>(get_d()); }
  /// return float value
  float floatValue() const { return static_cast<float>(get_d()); }
  /// return double value
  double doubleValue() const { return get_d(); }
  /// return BigInt value
  BigInt BigIntValue() const { return get_z(); }
  /// return BigRat value
  BigRat BigRatValue() const { return get_q(); }
  //@}
#endif
private:
  // assignment
  template <template <typename, typename, typename> class Policy>
  bool _set_f(const BigFloat2& x, prec_t prec = 0);
  template <template <typename, typename, typename> class Policy>
  bool _set_str(const char* x, int base, prec_t prec = 0);
  template <template <typename, typename, typename> class Policy, typename T>
  bool _set(const T& x, prec_t prec = 0);
  template <template <typename, typename, typename> class Policy>
  bool _set_2exp_si(long, exp_t e, prec_t prec = 0);
  template <template <typename, typename, typename> class Policy>
  bool _set_2exp_ui(unsigned long, exp_t e, prec_t prec = 0);

  // negation
  template <template <typename, typename, typename> class Policy>
  bool _neg_f(const BigFloat2& x, prec_t prec = 0);
  template <template <typename, typename, typename> class Policy, typename T>
  bool _neg(const T& x, prec_t prec = 0);
  
  // square root
  template <template <typename, typename, typename> class Policy>
  bool _sqrt_f(const BigFloat2& x, prec_t prec = 0);
  template <template <typename, typename, typename> class Policy, typename T>
  bool _sqrt(const T& x, prec_t prec = 0);
  
  // cubic root
  template <template <typename, typename, typename> class Policy>
  bool _cbrt_f(const BigFloat2& x, prec_t prec = 0);
  template <template <typename, typename, typename> class Policy, typename T>
  bool _cbrt(const T& x, prec_t prec = 0);
  
  // k-th root
  template <template <typename, typename, typename> class Policy>
  bool _root_f(const BigFloat2& x, unsigned long k, prec_t prec = 0);
  template <template <typename, typename, typename> class Policy, typename T>
  bool _root(const T& x, unsigned long k, prec_t prec = 0);
  
  // sine
  template <template <typename, typename, typename> class Policy>
  bool _sin_f(const BigFloat2& x, prec_t prec = 0);
  template <template <typename, typename, typename> class Policy, typename T>
  bool _sin(const T& x, prec_t prec = 0);
 
  // cosine
  template <template <typename, typename, typename> class Policy>
  bool _cos_f(const BigFloat2& x, prec_t prec = 0);
  template <template <typename, typename, typename> class Policy, typename T>
  bool _cos(const T& x, prec_t prec = 0);
 
  // tangent
  template <template <typename, typename, typename> class Policy>
  bool _tan_f(const BigFloat2& x, prec_t prec = 0);
  template <template <typename, typename, typename> class Policy, typename T>
  bool _tan(const T& x, prec_t prec = 0);
 
  // cotangent
  template <template <typename, typename, typename> class Policy>
  bool _cot_f(const BigFloat2& x, prec_t prec = 0);
  template <template <typename, typename, typename> class Policy, typename T>
  bool _cot(const T& x, prec_t prec = 0);
 
  // arcsine
  template <template <typename, typename, typename> class Policy>
  bool _asin_f(const BigFloat2& x, prec_t prec = 0);
  template <template <typename, typename, typename> class Policy, typename T>
  bool _asin(const T& x, prec_t prec = 0);
 
  // arccosine
  template <template <typename, typename, typename> class Policy>
  bool _acos_f(const BigFloat2& x, prec_t prec = 0);
  template <template <typename, typename, typename> class Policy, typename T>
  bool _acos(const T& x, prec_t prec = 0);
 
  // arctan
  template <template <typename, typename, typename> class Policy>
  bool _atan_f(const BigFloat2& x, prec_t prec = 0);
  template <template <typename, typename, typename> class Policy, typename T>
  bool _atan(const T& x, prec_t prec = 0);

  // log
  template <template <typename, typename, typename> class Policy>
  bool _log_f(const BigFloat2& x, prec_t prec = 0);
  template <template <typename, typename, typename> class Policy, typename T>
  bool _log(const T& x, prec_t prec = 0);
  
  // log2
  template <template <typename, typename, typename> class Policy>
  bool _log2_f(const BigFloat2& x, prec_t prec = 0);
  template <template <typename, typename, typename> class Policy, typename T>
  bool _log2(const T& x, prec_t prec = 0);
 
  // log10
  template <template <typename, typename, typename> class Policy>
  bool _log10_f(const BigFloat2& x, prec_t prec = 0);
  template <template <typename, typename, typename> class Policy, typename T>
  bool _log10(const T& x, prec_t prec = 0);

  // exponent
  template <template <typename, typename, typename> class Policy>
  bool _exp_f(const BigFloat2& x, prec_t prec = 0);
  template <template <typename, typename, typename> class Policy, typename T>
  bool _exp(const T& x, prec_t prec = 0);
 
  // exponent of 2
  template <template <typename, typename, typename> class Policy>
  bool _exp2_f(const BigFloat2& x, prec_t prec = 0);
  template <template <typename, typename, typename> class Policy, typename T>
  bool _exp2(const T& x, prec_t prec = 0);

  // exponent of 10
  template <template <typename, typename, typename> class Policy>
  bool _exp10_f(const BigFloat2& x, prec_t prec = 0);
  template <template <typename, typename, typename> class Policy, typename T>
  bool _exp10(const T& x, prec_t prec = 0);

  // addition
  template <template <typename, typename, typename> class Policy>
  bool _add_f(const BigFloat2& x, const BigFloat2& y, prec_t prec = 0);
  template <template <typename, typename, typename> class Policy, typename T>
  bool _add(const BigFloat2& x, const T& y, prec_t prec = 0);
  template <template <typename, typename, typename> class Policy, typename T>
  bool _add(const T& x, const BigFloat2& y, prec_t prec = 0);

  // subtraction 
  template <template <typename, typename, typename> class Policy>
  bool _sub_f(const BigFloat2& x, const BigFloat2& y, prec_t prec = 0);
  template <template <typename, typename, typename> class Policy, typename T>
  bool _sub(const BigFloat2& x, const T& y, prec_t prec = 0);
  template <template <typename, typename, typename> class Policy, typename T>
  bool _sub(const T& x, const BigFloat2& y, prec_t prec = 0);

  // multiplication
  template <template <typename, typename, typename> class Policy>
  bool _mul_f(const BigFloat2& x, const BigFloat2& y, prec_t prec = 0);
  template <template <typename, typename, typename> class Policy, typename T>
  bool _mul(const BigFloat2& x, const T& y, prec_t prec = 0);
  template <template <typename, typename, typename> class Policy, typename T>
  bool _mul(const T& x, const BigFloat2& y, prec_t prec = 0);

  // division
  template <template <typename, typename, typename> class Policy>
  bool _div_f(const BigFloat2& x, const BigFloat2& y, prec_t prec = 0);
  template <template <typename, typename, typename> class Policy, typename T>
  bool _div(const BigFloat2& x, const T& y, prec_t prec = 0);
  template <template <typename, typename, typename> class Policy, typename T>
  bool _div(const T& x, const BigFloat2& y, prec_t prec = 0);
};

/// \addtogroup BigFloat2ArithmeticOperators
//@{
/// BigFloat2 + BigFloat2
inline BigFloat2 operator+(const BigFloat2& x, const BigFloat2& y)
{ BigFloat2 r; r.add(x, y); return r; }
/// BigFloat2 + T
//template <typename T> inline BigFloat2 operator+(const BigFloat2& x, const T& y)
//{ BigFloat2 r; r.add(x, y); return r; }
/// T + BigFloat2
//template <typename T> inline BigFloat2 operator+(const T& x, const BigFloat2& y)
//{ BigFloat2 r; r.add(x, y); return r; }

/// BigFloat2 - BigFloat2
inline BigFloat2 operator-(const BigFloat2& x, const BigFloat2& y)
{ BigFloat2 r; r.sub(x, y); return r; }

#ifndef AF_DONT_DEFINE_MINUS_FOR_BIGFLOAT_T // AF: the following templates lead to mismatches
/// BigFloat2 - T
template <typename T> inline BigFloat2 operator-(const BigFloat2& x, const T& y)
{ BigFloat2 r; r.sub(x, y); return r; }
/// T - BigFloat2
template <typename T> inline BigFloat2 operator-(const T& x, const BigFloat2& y)
{ BigFloat2 r; r.sub(x, y); return r; }
#endif

/// BigFloat2 * BigFloat2
inline BigFloat2 operator*(const BigFloat2& x, const BigFloat2& y)
{ BigFloat2 r; r.mul(x, y); return r; }
/// BigFloat2 * T
template <typename T> inline BigFloat2 operator*(const BigFloat2& x, const T& y)
{ BigFloat2 r; r.mul(x, y); return r; }
/// T * BigFloat2
template <typename T> inline BigFloat2 operator*(const T& x, const BigFloat2& y)
{ BigFloat2 r; r.mul(x, y); return r; }

/// BigFloat2 / BigFloat2 (w/ default precision getDefaultBFdivPrec())
inline BigFloat2 operator/(const BigFloat2& x, const BigFloat2& y)
{ BigFloat2 r; r.div(x, y); return r; }
/// BigFloat2 / T (w/ default precision getDefaultBFdivPrec())
template <typename T> inline BigFloat2 operator/(const BigFloat2& x, const T& y)
{ BigFloat2 r; r.div(x, y); return r; }
/// T / BigFloat2 (w/ default precision getDefaultBFdivPrec())
template <typename T> inline BigFloat2 operator/(const T& x, const BigFloat2& y)
{ BigFloat2 r; r.div(x, y); return r; }
//@}

/// \addtogroup BigFloat2IostreamOperators
//@{
/// istream operator for <tt>BigFloat2</tt>
inline std::istream& operator>>(std::istream& is, BigFloat2& x)
{ BigFloat2::FT tmp; is >> tmp; x.set(tmp); return is; }
/// ostream operator for <tt>BigFloat2</tt>
inline std::ostream& operator<<(std::ostream& os, const BigFloat2& x)
{
  return os << x.get_str();
}
//@}

/// \addtogroup BigFloat2GlobalFunctions
//@{
/// square root
inline BigFloat2 sqrt(const BigFloat2& x, prec_t prec = getDefaultBFradicalPrec())
{ BigFloat2 r; r.sqrt(x, prec); return r; }
/// cubic root
inline BigFloat2 cbrt(const BigFloat2& x, prec_t prec = getDefaultBFradicalPrec())
{ BigFloat2 r; r.cbrt(x, prec); return r; }
/// k-th root
inline BigFloat2 root(const BigFloat2& x, unsigned long k, prec_t prec = getDefaultBFradicalPrec())
{ BigFloat2 r; r.root(x, k, prec); return r; }
//@}

// include inline functions (private)
#include <CORE/BigFloat2.inl>
inline int sign(const BigFloat2& x) 
{ return x.sgn(); }

#ifdef CORE_OLDNAMES 
/// \addtogroup BigFloat2BackCompatiableFunctions
//@{
//@}
#endif

CORE_END_NAMESPACE

#endif /*__CORE_BIGFLOAT2_H__*/
