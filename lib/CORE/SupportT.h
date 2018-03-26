/******************************************************************
 * File: SupportT.h
 * Synopsis:
 *      Support functions to provide a uniform API over trignometric
 *      functions so they can be used  by the complex (i.e., complex numbers)
 *      classes. Note that
 *      at level 2 the precision is determined by DEFAULT_PREC.
 *
 * Written by
 *       Narayan Kamath (kamath.narayan@gmail.com) 2010
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *****************************************************************/


#ifndef SUPPORT_T_H_
#define SUPPORT_T_H_

// NOTE(narayan) :
// This file is a giant hack. This is to gloss over
// discrepancies in API between various operation levels.
// This namespace exposes a uniform API :
//
// pi<NT>()
// sinT<NT>(const T &)
// cosT<NT>(const T &)
// atanT<NT>(const T &)
//
// The main purpose of this class is to convert complex
// numbers in and out of polar form. Note that the trignometric
// operations have inherent exactness issues (except at Level 3)
// so the conversion will not be exact. The level of precision
// can be controlled by setting complex_support::DEFAULT_PREC
// at the start of execution.
namespace complex_support {

// Change this if you need a higher level of precision for
// sin and cos.
static prec_t DEFAULT_PREC = DOUBLE_PREC;

// NOTE(narayan) :
//
// If these functions are moved after the point of instantiation,
// the compiler does not seem to find the specializations, and this
// results in weird linker errors.
//
// This function should be moved to core at some point.
// For now i'll leave it here to avoid too much trouble.
template <typename NT> inline NT piT() {
  NT type;
  type.pi(DEFAULT_PREC, MPFR_RND);
  return type;
}
template < > inline machine_double piT() {
  DoubleWrapper p;
  p.pi();
  return p.doubleValue();
}

// arctan
// ------
template <typename NT> inline NT atanT(const NT &in) {
  return atan(in);
}

template <> inline machine_double atanT(const machine_double &in) {
  return ::atan(in);
}

// sin
// ---
template <typename NT> inline NT sinT(const NT &in) {
  return sin(in);
}
template <> inline machine_double sinT(const machine_double &in) {
  return ::sin(in);
}
// cos
// ---
template <typename NT> inline NT cosT(const NT &in) {
  return cos(in);
}

template <> inline machine_double cosT(const machine_double &in) {
  return ::cos(in);
}

#if CORE_LEVEL >= 2 && CORE_LEVEL != 3
template <> inline BigFloat atanT(const BigFloat &in) {
  BigFloat type;
  type.atan(in, DEFAULT_PREC, MPFR_RND);
  return type;
}
template <> inline BigFloat cosT(const BigFloat &in) {
  BigFloat type;
  type.cos(in, DEFAULT_PREC, MPFR_RND);
  return type;
}
template <> inline BigFloat sinT(const BigFloat &in) {
  BigFloat type;
  type.sin(in, DEFAULT_PREC, MPFR_RND);
  return type;
}
#endif


#if CORE_LEVEL >= 3
// Define a specialization for Expr.
template < > inline Expr piT() {
  // This is defined in Expr.h
  return CORE::pi();
}
#endif

}  // namespace complex_suppot.


#endif /* SUPPORTT_H_ */
