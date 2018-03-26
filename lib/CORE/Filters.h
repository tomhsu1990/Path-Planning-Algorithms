/****************************************************************************
 * Filters.h -- Floating-point filters
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
 * $Id: Filters.h,v 1.11 2010/11/23 17:58:37 exact Exp $
 ***************************************************************************/
#ifndef __CORE_FILTERS_H__
#define __CORE_FILTERS_H__

#include <iostream>

#include <cmath>
#include <cfloat>
#include <CORE/BigFloat2.h>

#if defined (_MSC_VER) || defined (__MINGW32__) // add support for MinGW
  #define finite(x)	_finite(x)
  #define ilogb(x)	(int)_logb(x)
#endif

#if defined(sun) || defined(__sun)
  #include <ieeefp.h>
#endif

CORE_BEGIN_NAMESPACE

/// \class DummyFilter
/// \brief a dummy filter, do nothing
class DummyFilter {
  typedef DummyFilter thisClass;
public:
#ifdef CORE_DEBUG_FILTER
  void dump() const {}
#endif
  bool is_ok() const { return false; }
  int sign() const { assert(0); return 0; }
  long uMSB() const { assert(0); return 0; }
  long lMSB() const { assert(0); return 0; }
  template <typename T> void set(const T&) {}
  void neg(const thisClass&) {}
  void sqrt(const thisClass&) {}
  void cbrt(const thisClass&) {}
  void root(const thisClass&, unsigned long) {}
  void addsub(const thisClass&, const thisClass&, bool) {}
  void mul(const thisClass&, const thisClass&) {}
  void div(const thisClass&, const thisClass&) {}
};

extern bool fpFilterFlag;

// Chee: this is needed to avoid the isfinite() error below:
//extern bool isfinite (float x);
//extern bool isfinite (double x);
//extern bool isfinite (long double x);

/// turn floating-point filter on/off
inline bool setFpFilterFlag(bool f) {
  bool oldf = fpFilterFlag;
  fpFilterFlag = f;
  return oldf;
}

// constants
const int IEEE_DOUBLE_PREC = 52;
const double DBL_INFTY = ::ldexp(DBL_MAX, 1);
const double CORE_EPS = ::ldexp(1.0, -IEEE_DOUBLE_PREC);

// k-th root for double (using BigFloat for now)
inline double root(double x, unsigned long k) 
{ BigFloat r(x); r.root(x, k); return r.get_d(); }

/// \class filteredFp Filter.h
/// \brief filteredFp is a simple filtered floating point number.
///        It is based on the Burnikel-Funke-Schirra (BFS) filter scheme.
///        We do not use IEEE exception mechanism here.
template <typename Kernel = BigFloat>
class BfsFilter {
  double fpVal;         // approximate double value for some "real value"
  double maxAbs;        // if (|fpVal| > maxAbs * ind * 2^{-52}) then
  int ind;              // sign of value is sign(fpVal).  Else, don't know.
  bool _isok;
  int _sign;
  long _uMSB;
  long _lMSB;
  prec_t _r_prec;
  prec_t _a_prec;
  // REFERENCE: Burnikel, Funke, Schirra (BFS) filter
  // Feb'07, Chee: in isOK(), you used the test "|fpVal| >= maxAbs * ind * 2^{-52}" 
  // which seems to be correct (i.e., not |fpVal| > maxAbs * ind * 2^{-52})
  // Apr'07, Jihun : Above is wrong. The correct test is |fpVal| > maxAbs * ind * 2^{-52}
  //  This is because when value is very close to zero, fpVal = maxAbs = 0
  typedef BfsFilter thisClass;
  typedef typename Kernel::ZT ZT;
  typedef typename Kernel::QT QT;
  typedef typename Kernel::FT FT;
private:
  void compute_cache () {
    double Val = maxAbs*ind*CORE_EPS;
  //  Sep'14, Chee: finite(fpVal) is deprecated on MacOS!  so must use isfinite(fpVal)
  //  Furthermore, there is a error of this nature:
  //	"there are no arguments to 'isfinite' that depend on a template parameter"
  //	See:
  //http://stackoverflow.com/questions/9941987/there-are-no-arguments-that-depend-on-a-template-parameter
  //  Solution is to declare the "extern bool isfinite()" above.
    _isok = isfinite(fpVal)&&(::fabs(fpVal)>Val);
    if (!_isok) return;
    _sign = (fpVal == 0.0) ? 0 : (fpVal > 0.0 ? 1: -1);
    _lMSB = (sign() == 0) ? MSB_MIN : long(ilogb(::fabs(fpVal) - Val));
    _uMSB = (sign() == 0) ? 0 : long(ilogb(::fabs(fpVal) + Val)+1);
    //_a_prec = (std::min)((std::max)(msb_t(-floorlg(Val)),msb_t(2)), 1024);
    //_r_prec = (std::min)((std::max)(long(get_a_prec()) + uMSB(), msb_t(2)), msb_t(IEEE_DOUBLE_PREC));
  }
public:
  BfsFilter()
  { _isok = false; }
#ifdef CORE_DEBUG_FILTER
  void dump() const 
  { std::cerr<<"[fpVal,maxAbs,ind]="<<fpVal<<","<<maxAbs<<","<<ind<<std::endl; }
#endif
  bool is_ok() const 
  { return fpFilterFlag&&_isok; }
  int sign() const 
  { assert(_isok); return _sign; }
  long lMSB() const 
  { assert(_isok); return _lMSB; }
  long uMSB() const 
  { assert(_isok); return _uMSB; }
  double get_value() const 
  { assert(_isok); return fpVal; }
  prec_t get_r_prec() const
  { assert(_isok); return _r_prec; }
  prec_t get_a_prec() const
  { assert(_isok); return _a_prec; }
  double r_approx(int prec) const
  { assert(_isok); return fpVal; }
  double a_approx(int prec) const
  { assert(_isok); return fpVal; }

  /// from "Exact Geometric Predicates using Cascaded Computation" P176
  ///
  ///   An input value $x$ exactly representable by a double has the 
  ///   floating-point approximation $\tidle{x}=x$, the supremum 
  ///   $\tidle{x_{sup}}=|x|$ and the index 0. Otherwise, 
  ///   $\tidle{x}=round(x)$, the supremum 
  ///   $\tidle{x_{sup}}=|\tidle{x}|=|round(x)|$ and the index 1.
  void set(long value) { 
    fpVal = value; maxAbs = value > 0 ? value : (-value); 
    ind = (sizeof(long) > 4 && ceillg(value) >= 53) ? 1 : 0;
    compute_cache();
  }
  void set(unsigned long value) {
    fpVal = value; maxAbs = value; 
    ind = (sizeof(unsigned long) > 4 && ceillg(value) >= 53) ? 1 : 0;
    compute_cache();
  }
  void set(double value) {
    fpVal = value; maxAbs = ::fabs(value); ind = 0;
    compute_cache();
  }
  void set(const ZT& value) { 
    fpVal = value.get_d(); maxAbs = ::fabs(fpVal); 
    ind = value.uMSB() >= 53 ? 1 : 0; 
    compute_cache();
  }
  void set(const QT& value) { 
    fpVal = value.get_d(); maxAbs = ::fabs(fpVal); 
    ind = 1; //value.uMSB() >= 53 ? 1 : 0; // ??? denonimator has to be power of 2
    compute_cache();
  }
  void set(const FT& value) {
    fpVal = value.get_d(); maxAbs = ::fabs(fpVal); 
    ind = value.get_prec() >= 53 ? 1 : 0;
    compute_cache();
  }
  void set(const Kernel& value) {
    fpVal = value.get_d(); maxAbs = ::fabs(fpVal); 
    ind = value.get_prec() >= 53 ? 1 : 0;
    compute_cache();
  }

  // negation
  void neg(const thisClass& child) {
    fpVal = -child.fpVal; maxAbs = child.maxAbs; ind = child.ind;
    compute_cache();
  }

  // square root
  void sqrt(const thisClass& child) {
    if (child.fpVal > 0.0) {
      fpVal = ::sqrt(child.fpVal); maxAbs = child.maxAbs / child.fpVal * fpVal;
    } else {
      fpVal = 0.0; maxAbs = ::ldexp(::sqrt(child.maxAbs), 26);
    }
    ind = 1 + child.ind;
    compute_cache();
  }

  // cubic root
  void cbrt(const thisClass& child) {
    if (child.fpVal > 0.0) {
      fpVal = ::cbrt(child.fpVal); maxAbs = child.maxAbs / child.fpVal * fpVal;
    } else {
      fpVal = 0.0; maxAbs = ::ldexp(::cbrt(child.maxAbs), 18);
    }
    ind = 1 + child.ind;
    compute_cache();
  }

  // k-th root
  // TODO: root() function for double, using newton ??
  void root(const thisClass& child, unsigned long k) {
    if (child.fpVal > 0.0) {
      fpVal = CORE_NS::root(child.fpVal, k); 
      maxAbs = child.maxAbs / child.fpVal * fpVal;
    } else {
      fpVal = 0.0; 
      maxAbs = ::ldexp(CORE_NS::root(child.maxAbs, k), (IEEE_DOUBLE_PREC+k-1)/k);
    }
    ind = 1 + child.ind;
    compute_cache();
  }

  // addition/subtraction
  void addsub(const thisClass& f, const thisClass& s, bool is_add) {
    fpVal = is_add ? f.fpVal + s.fpVal : f.fpVal - s.fpVal;
    maxAbs = f.maxAbs + s.maxAbs;
    ind = 1 + (f.ind > s.ind ? f.ind : s.ind);
    compute_cache();
  }

  // multiplication
  void mul(const thisClass& f, const thisClass& s) {
    fpVal = f.fpVal * s.fpVal;
    maxAbs = f.maxAbs * s.maxAbs + DBL_MIN;
    ind = 1 + f.ind + s.ind; 
    compute_cache();
  }

  // division
  void div(const thisClass& f, const thisClass& s) {
    double xxx = ::fabs(s.fpVal) / s.maxAbs - (s.ind+1)*CORE_EPS;
    if (xxx > 0) {
      fpVal = f.fpVal / s.fpVal;
      maxAbs = (::fabs(f.fpVal)/::fabs(s.fpVal)+f.maxAbs/s.maxAbs)/xxx+DBL_MIN;
      ind = 1 + (f.ind > s.ind+1 ? f.ind : s.ind+1);
    } else {
      fpVal = DBL_INFTY;
      maxAbs = 0.0;
      ind = 0;
    }
    compute_cache();
  }
};

CORE_END_NAMESPACE

#endif /*__CORE_FILTERS_H__*/
