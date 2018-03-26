/****************************************************************************
 * Policies.h -- Computation policies for BigFloat2
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
 * $Id: Policies.h,v 1.10 2010/05/06 20:49:30 exact Exp $
 ***************************************************************************/
#ifndef __CORE_POLICIES_H__
#define __CORE_POLICIES_H__

CORE_BEGIN_NAMESPACE

// raw version (user need to set precision for result before call)
template <typename T, typename TL = T, typename TR = T>
struct RawArithmeticPolicy {
  static bool set(T& z, const TL& x, prec_t, rnd_t rnd = MPFR_RND)
  { return z.r_set(x, rnd) == 0; }
  static bool set(T& z, const char* x, int b, prec_t, rnd_t rnd = MPFR_RND)
  { return z.r_set(x, b, rnd) == 0; }
  static bool set_2exp(T& z, long x, exp_t e, prec_t, rnd_t rnd)
  { return z.r_set_2exp(x, e, rnd) == 0; }
  static bool set_2exp(T& z, unsigned long x, exp_t e, prec_t, rnd_t rnd)
  { return z.r_set_2exp(x, e, rnd) == 0; }

  static bool neg(T& z, const TL& x, prec_t, rnd_t rnd)
  { return z.r_neg(x, rnd) == 0; }
  static bool sqrt(T& z, const TL& x, prec_t, rnd_t rnd)
  { return z.r_sqrt(x, rnd) == 0; }
  static bool cbrt(T& z, const TL& x, prec_t, rnd_t rnd)
  { return z.r_cbrt(x, rnd) == 0; }
  static bool root(T& z, const TL& x, unsigned long k, prec_t, rnd_t rnd)
  { return z.r_root(x, k, rnd) == 0; }

  static bool sin(T& z, const TL& x, prec_t, rnd_t rnd)
  { return z.r_sin(x, rnd) == 0; }
  static bool cos(T& z, const TL& x, prec_t, rnd_t rnd)
  { return z.r_cos(x, rnd) == 0; }
  static bool tan(T& z, const TL& x, prec_t, rnd_t rnd)
  { return z.r_tan(x, rnd) == 0; }
  static bool cot(T& z, const TL& x, prec_t, rnd_t rnd)
  { return z.r_cot(x, rnd) == 0; }
  static bool asin(T& z, const TL& x, prec_t, rnd_t rnd)
  { return z.r_asin(x, rnd) == 0; }
  static bool acos(T& z, const TL& x, prec_t, rnd_t rnd)
  { return z.r_acos(x, rnd) == 0; }
  static bool atan(T& z, const TL& x, prec_t, rnd_t rnd)
  { return z.r_atan(x, rnd) == 0; }
  static bool log(T& z, const TL& x, prec_t, rnd_t rnd)
  { return z.r_log(x, rnd) == 0; }
  static bool log2(T& z, const TL& x, prec_t, rnd_t rnd)
  { return z.r_log2(x, rnd) == 0; }
  static bool log10(T& z, const TL& x, prec_t, rnd_t rnd)
  { return z.r_log10(x, rnd) == 0; }
  static bool exp(T& z, const TL& x, prec_t, rnd_t rnd)
  { return z.r_exp(x, rnd) == 0; }
  static bool exp2(T& z, const TL& x, prec_t, rnd_t rnd)
  { return z.r_exp2(x, rnd) == 0; }
  static bool exp10(T& z, const TL& x, prec_t, rnd_t rnd)
  { return z.r_exp10(x, rnd) == 0; }

  static bool add(T& z, const TL& x, const TR& y, prec_t, rnd_t rnd)
  { return z.r_add(x, y, rnd) == 0; }
  static bool sub(T& z, const TL& x, const TR& y, prec_t, rnd_t rnd)
  { return z.r_sub(x, y, rnd) == 0; }
  static bool mul(T& z, const TL& x, const TR& y, prec_t, rnd_t rnd)
  { return z.r_mul(x, y, rnd) == 0; }
  static bool div(T& z, const TL& x, const TR& y, prec_t, rnd_t rnd)
  { return z.r_div(x, y, rnd) == 0; }
};

// fix version (call with specified precision)
template <typename T, typename TL = T, typename TR = T>
struct FixedArithmeticPolicy {
  static bool set(T& z, const TL& x, prec_t prec, rnd_t rnd = MPFR_RND)
  { return z.set(x, prec, rnd) == 0; }
  static bool set(T& z, const char* x, int b, prec_t prec, rnd_t rnd = MPFR_RND)
  { return z.set(x, b, prec, rnd) == 0; }
  static bool set_2exp(T& z, long x, exp_t e, prec_t prec, rnd_t rnd)
  { return z.set_2exp(x, e, prec, rnd) == 0; }
  static bool set_2exp(T& z, unsigned long x, exp_t e, prec_t prec, rnd_t rnd)
  { return z.set_2exp(x, e, prec, rnd) == 0; }

  static bool neg(T& z, const TL& x, prec_t prec, rnd_t rnd)
  { return z.neg(x, prec, rnd) == 0; }
  static bool sqrt(T& z, const TL& x, prec_t prec, rnd_t rnd)
  { return z.sqrt(x, prec, rnd) == 0; }
  static bool cbrt(T& z, const TL& x, prec_t prec, rnd_t rnd)
  { return z.cbrt(x, prec, rnd) == 0; }
  static bool root(T& z, const TL& x, unsigned long k, prec_t prec, rnd_t rnd)
  { return z.root(x, k, prec, rnd) == 0; }

  static bool sin(T& z, const TL& x, prec_t prec, rnd_t rnd)
  { return z.sin(x, prec, rnd) == 0; }
  static bool cos(T& z, const TL& x, prec_t prec, rnd_t rnd)
  { return z.cos(x, prec, rnd) == 0; }
  static bool tan(T& z, const TL& x, prec_t prec, rnd_t rnd)
  { return z.tan(x, prec, rnd) == 0; }
  static bool cot(T& z, const TL& x, prec_t prec, rnd_t rnd)
  { return z.cot(x, prec, rnd) == 0; }
  static bool asin(T& z, const TL& x, prec_t prec, rnd_t rnd)
  { return z.asin(x, prec, rnd) == 0; }
  static bool acos(T& z, const TL& x, prec_t prec, rnd_t rnd)
  { return z.acos(x, prec, rnd) == 0; }
  static bool atan(T& z, const TL& x, prec_t prec, rnd_t rnd)
  { return z.atan(x, prec, rnd) == 0; }
  static bool log(T& z, const TL& x, prec_t prec, rnd_t rnd)
  { return z.log(x, prec, rnd) == 0; }
  static bool log2(T& z, const TL& x, prec_t prec, rnd_t rnd)
  { return z.log2(x, prec, rnd) == 0; }
  static bool log10(T& z, const TL& x, prec_t prec, rnd_t rnd)
  { return z.log10(x, prec, rnd) == 0; }
  static bool exp(T& z, const TL& x, prec_t prec, rnd_t rnd)
  { return z.exp(x, prec, rnd) == 0; }
  static bool exp2(T& z, const TL& x, prec_t prec, rnd_t rnd)
  { return z.exp2(x, prec, rnd) == 0; }
  static bool exp10(T& z, const TL& x, prec_t prec, rnd_t rnd)
  { return z.exp10(x, prec, rnd) == 0; }

  static bool add(T& z, const TL& x, const TR& y, prec_t prec, rnd_t rnd)
  { return z.add(x, y, prec, rnd) == 0; }
  static bool sub(T& z, const TL& x, const TR& y, prec_t prec, rnd_t rnd)
  { return z.sub(x, y, prec, rnd) == 0; }
  static bool mul(T& z, const TL& x, const TR& y, prec_t prec, rnd_t rnd)
  { return z.mul(x, y, prec, rnd) == 0; }
  static bool div(T& z, const TL& x, const TR& y, prec_t prec, rnd_t rnd)
  { return z.div(x, y, prec, rnd) == 0; }
};

// auto version (precision will be determinated automatically)
template <typename T, typename TL = T, typename TR = T>
struct AutoArithmeticPolicy {
  static bool set(T& z, const TL& x, prec_t, rnd_t rnd = MPFR_RND)
  { return z.set(x, rnd) == 0; }
  static bool set(T& z, const char* x, int b, prec_t, rnd_t rnd = MPFR_RND)
  { return z.set(x, b, rnd) == 0; }
  static bool set_2exp(T& z, long x, exp_t e, prec_t, rnd_t rnd)
  { return z.set_2exp(x, e, rnd) == 0; }
  static bool set_2exp(T& z, unsigned long x, exp_t e, prec_t, rnd_t rnd)
  { return z.set_2exp(x, e, rnd) == 0; }

  static bool neg(T& z, const TL& x, prec_t, rnd_t rnd)
  { return z.neg(x, rnd) == 0; }

  static bool add(T& z, const TL& x, const TR& y, prec_t, rnd_t rnd)
  { return z.add(x, y, rnd) == 0; }
  static bool sub(T& z, const TL& x, const TR& y, prec_t, rnd_t rnd)
  { return z.sub(x, y, rnd) == 0; }
  static bool mul(T& z, const TL& x, const TR& y, prec_t, rnd_t rnd)
  { return z.mul(x, y, rnd) == 0; }
};

CORE_END_NAMESPACE

#endif /*__CORE_POLICIES_H__*/
