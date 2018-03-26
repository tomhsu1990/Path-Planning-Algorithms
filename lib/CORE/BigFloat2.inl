/****************************************************************************
 * BigFloat2.inl -- Inline functions for BigFloat2
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
 * $Id: BigFloat2.inl,v 1.12 2010/07/13 14:52:38 exact Exp $
 ***************************************************************************/
#define BF_RNDD GMP_RNDD
#define BF_RNDU GMP_RNDU

////////////////////////////////////////////////////////////////////////////////
/// assignment -- set(BigFloat2)
template <template <typename, typename, typename> class Policy>
inline bool BigFloat2::_set_f(const BigFloat2& x, prec_t prec) {
  if (x.is_exact())
    return this->_set<Policy, FT>(x.m_l, prec);
  else {
    Policy<FT, FT, FT>::set(m_l, x.m_l, prec, BF_RNDD);
    Policy<FT, FT, FT>::set(m_r, x.m_r, prec, BF_RNDU);
    set_exact(false);
  }
  return is_exact();
}
/// assignment -- set(const char*, int base)
template <template <typename, typename, typename> class Policy>
inline bool BigFloat2::_set_str(const char* x, int base, prec_t prec) {
  // since MPFR set_str() function cannot tell the exactness of result,
  // we need compare the two results after conversion
  Policy<FT, FT, FT>::set(m_l, x, base, prec, BF_RNDD);
  Policy<FT, FT, FT>::set(m_r, x, base, prec, BF_RNDU);
  set_exact(m_l.cmp(m_r));
  return is_exact();
}
/// assignment -- set(T)
template <template <typename, typename, typename> class Policy, typename T>
inline bool BigFloat2::_set(const T& x, prec_t prec) {
  set_exact(Policy<FT, T, FT>::set(m_l, x, prec, BF_RNDD));
  if (!is_exact()) Policy<FT, T, FT>::set(m_r, x, prec, BF_RNDU);
  return is_exact();
}
/// assignment -- set_2exp(long)
template <template <typename, typename, typename> class Policy>
inline bool BigFloat2::_set_2exp_si(long x, exp_t e, prec_t prec) {
  set_exact(Policy<FT, FT, FT>::set_2exp(m_l, x, e, prec, BF_RNDD));
  if (!is_exact()) Policy<FT, FT, FT>::set_2exp(m_l, x, e, prec, BF_RNDU);
  return is_exact();
}
/// assignment -- set_2exp(long)
template <template <typename, typename, typename> class Policy>
inline bool BigFloat2::_set_2exp_ui(unsigned long x, exp_t e, prec_t prec) {
  set_exact(Policy<FT, FT, FT>::set_2exp(m_l, x, e, prec, BF_RNDD));
  if (!is_exact()) Policy<FT, FT, FT>::set_2exp(m_l, x, e, prec, BF_RNDU);
  return is_exact();
}

////////////////////////////////////////////////////////////////////////////////
/// negation -- neg(BigFloat2)
template <template <typename, typename, typename> class Policy>
inline bool BigFloat2::_neg_f(const BigFloat2& x, prec_t prec) {
  if (x.is_exact())
    return this->_neg<Policy, FT>(x.m_l, prec);
  else {
    Policy<FT, FT, FT>::neg(m_l, x.m_r, prec, BF_RNDD);
    Policy<FT, FT, FT>::neg(m_r, x.m_l, prec, BF_RNDU);
    set_exact(false);
  }
  return is_exact();
}
/// negation -- neg(T)
template <template <typename, typename, typename> class Policy, typename T>
inline bool BigFloat2::_neg(const T& x, prec_t prec) {
  set_exact(Policy<FT, T, FT>::neg(m_l, x, prec, BF_RNDD));
  if (!is_exact()) Policy<FT, T, FT>::neg(m_r, x, prec, BF_RNDU);
  return is_exact();
}

////////////////////////////////////////////////////////////////////////////////
/// square root -- sqrt(BigFloat2)
template <template <typename, typename, typename> class Policy>
inline bool BigFloat2::_sqrt_f(const BigFloat2& x, prec_t prec) {
  if (x.is_exact())
    return this->_sqrt<Policy, FT>(x.m_l, prec);
  else {
    Policy<FT, FT, FT>::sqrt(m_l, x.m_l, prec, BF_RNDD);
    Policy<FT, FT, FT>::sqrt(m_r, x.m_r, prec, BF_RNDU);
    set_exact(false);
  }
  return is_exact();
}
/// square root -- sqrt(T)
template <template <typename, typename, typename> class Policy, typename T>
inline bool BigFloat2::_sqrt(const T& x, prec_t prec) {
  set_exact(Policy<FT, T, FT>::sqrt(m_l, x, prec, BF_RNDD));
  if (!is_exact()) Policy<FT, T, FT>::sqrt(m_r, x, prec, BF_RNDU);
  return is_exact();
}

////////////////////////////////////////////////////////////////////////////////
/// cubic root -- cbrt(BigFloat2)
template <template <typename, typename, typename> class Policy>
inline bool BigFloat2::_cbrt_f(const BigFloat2& x, prec_t prec) {
  if (x.is_exact())
    return this->_cbrt<Policy, FT>(x.m_l, prec);
  else {
    Policy<FT, FT, FT>::cbrt(m_l, x.m_l, prec, BF_RNDD);
    Policy<FT, FT, FT>::cbrt(m_r, x.m_r, prec, BF_RNDU);
    set_exact(false);
  }
  return is_exact();
}
/// cubic root -- cbrt(T)
template <template <typename, typename, typename> class Policy, typename T>
inline bool BigFloat2::_cbrt(const T& x, prec_t prec) {
  set_exact(Policy<FT, T, FT>::cbrt(m_l, x, prec, BF_RNDD));
  if (!is_exact()) Policy<FT, T, FT>::cbrt(m_r, x, prec, BF_RNDU);
  return is_exact();
}

////////////////////////////////////////////////////////////////////////////////
/// k-th root -- root(BigFloat2)
template <template <typename, typename, typename> class Policy>
inline bool BigFloat2::_root_f(const BigFloat2& x, unsigned long k, prec_t prec) {
  if (x.is_exact())
    return this->_root<Policy, FT>(x.m_l, k, prec);
  else {
    Policy<FT, FT, FT>::root(m_l, x.m_l, k, prec, BF_RNDD);
    Policy<FT, FT, FT>::root(m_r, x.m_r, k, prec, BF_RNDU);
    set_exact(false);
  }
  return is_exact();
}
/// k-th root -- root(T)
template <template <typename, typename, typename> class Policy, typename T>
inline bool BigFloat2::_root(const T& x, unsigned long k, prec_t prec) {
  set_exact(Policy<FT, T, FT>::root(m_l, x, k, prec, BF_RNDD));
  if (!is_exact()) Policy<FT, T, FT>::root(m_r, x, k, prec, BF_RNDU);
  return is_exact();
}

////////////////////////////////////////////////////////////////////////////////
/// sine -- sin(BigFloat2)
template <template <typename, typename, typename> class Policy>
inline bool BigFloat2::_sin_f(const BigFloat2& x, prec_t prec) {
  if (x.is_exact())
    return this->_sin<Policy, FT>(x.m_l, prec);
  else {
    Policy<FT, FT, FT>::sin(m_l, x.m_l, prec, BF_RNDD);
    Policy<FT, FT, FT>::sin(m_r, x.m_r, prec, BF_RNDU);
    set_exact(false);
  }
  return is_exact();
}
/// sine -- sin(T)
template <template <typename, typename, typename> class Policy, typename T>
inline bool BigFloat2::_sin(const T& x, prec_t prec) {
  set_exact(Policy<FT, T, FT>::sin(m_l, x, prec, BF_RNDD));
  if (!is_exact()) Policy<FT, T, FT>::sin(m_r, x, prec, BF_RNDU);
  return is_exact();
}

////////////////////////////////////////////////////////////////////////////////
/// cosine -- cos(BigFloat2)
template <template <typename, typename, typename> class Policy>
inline bool BigFloat2::_cos_f(const BigFloat2& x, prec_t prec) {
  if (x.is_exact()) {
    return this->_cos<Policy, FT>(x.m_l, prec);
  } else {
    Policy<FT, FT, FT>::cos(m_l, x.m_l, prec, BF_RNDU);
    Policy<FT, FT, FT>::cos(m_r, x.m_r, prec, BF_RNDD);

    set_exact(false);
  }
  return is_exact();
}
/// cosine -- cos(T)
template <template <typename, typename, typename> class Policy, typename T>
inline bool BigFloat2::_cos(const T& x, prec_t prec) {
  set_exact(Policy<FT, T, FT>::cos(m_l, x, prec, BF_RNDU));
  if (!is_exact()) Policy<FT, T, FT>::cos(m_r, x, prec, BF_RNDD);
  return is_exact();
}

////////////////////////////////////////////////////////////////////////////////
/// tangent -- tan(BigFloat2)
template <template <typename, typename, typename> class Policy>
inline bool BigFloat2::_tan_f(const BigFloat2& x, prec_t prec) {
  if (x.is_exact())
    return this->_tan<Policy, FT>(x.m_l, prec);
  else {
    Policy<FT, FT, FT>::tan(m_l, x.m_l, prec, BF_RNDD);
    Policy<FT, FT, FT>::tan(m_r, x.m_r, prec, BF_RNDU);
    set_exact(false);
  }
  return is_exact();
}
/// tangent -- tan(T)
template <template <typename, typename, typename> class Policy, typename T>
inline bool BigFloat2::_tan(const T& x, prec_t prec) {
  set_exact(Policy<FT, T, FT>::tan(m_l, x, prec, BF_RNDD));
  if (!is_exact()) Policy<FT, T, FT>::tan(m_r, x, prec, BF_RNDU);
  return is_exact();
}

////////////////////////////////////////////////////////////////////////////////
/// cotangent -- cot(BigFloat2)
template <template <typename, typename, typename> class Policy>
inline bool BigFloat2::_cot_f(const BigFloat2& x, prec_t prec) {
  if (x.is_exact())
    return this->_cot<Policy, FT>(x.m_l, prec);
  else {
    Policy<FT, FT, FT>::cot(m_l, x.m_l, prec, BF_RNDD);
    Policy<FT, FT, FT>::cot(m_r, x.m_r, prec, BF_RNDU);
    set_exact(false);
  }
  return is_exact();
}
/// cotangent -- cot(T)
template <template <typename, typename, typename> class Policy, typename T>
inline bool BigFloat2::_cot(const T& x, prec_t prec) {
  set_exact(Policy<FT, T, FT>::cot(m_l, x, prec, BF_RNDD));
  if (!is_exact()) Policy<FT, T, FT>::cot(m_r, x, prec, BF_RNDU);
  return is_exact();
}

////////////////////////////////////////////////////////////////////////////////
/// arcsine -- asin(BigFloat2)
template <template <typename, typename, typename> class Policy>
inline bool BigFloat2::_asin_f(const BigFloat2& x, prec_t prec) {
  if (x.is_exact())
    return this->_asin<Policy, FT>(x.m_l, prec);
  else {
    Policy<FT, FT, FT>::asin(m_l, x.m_l, prec, BF_RNDD);
    Policy<FT, FT, FT>::asin(m_r, x.m_r, prec, BF_RNDU);
    set_exact(false);
  }
  return is_exact();
}
/// arcsine -- asin(T)
template <template <typename, typename, typename> class Policy, typename T>
inline bool BigFloat2::_asin(const T& x, prec_t prec) {
  set_exact(Policy<FT, T, FT>::asin(m_l, x, prec, BF_RNDD));
  if (!is_exact()) Policy<FT, T, FT>::asin(m_r, x, prec, BF_RNDU);
  return is_exact();
}

////////////////////////////////////////////////////////////////////////////////
/// arccosine -- acos(BigFloat2)
template <template <typename, typename, typename> class Policy>
inline bool BigFloat2::_acos_f(const BigFloat2& x, prec_t prec) {
  if (x.is_exact())
    return this->_acos<Policy, FT>(x.m_l, prec);
  else {
    Policy<FT, FT, FT>::acos(m_l, x.m_l, prec, BF_RNDD);
    Policy<FT, FT, FT>::acos(m_r, x.m_r, prec, BF_RNDU);
    set_exact(false);
  }
  return is_exact();
}
/// arccosine -- acos(T)
template <template <typename, typename, typename> class Policy, typename T>
inline bool BigFloat2::_acos(const T& x, prec_t prec) {
  set_exact(Policy<FT, T, FT>::acos(m_l, x, prec, BF_RNDD));
  if (!is_exact()) Policy<FT, T, FT>::acos(m_r, x, prec, BF_RNDU);
  return is_exact();
}

////////////////////////////////////////////////////////////////////////////////
/// arctangent -- atan(BigFloat2)
template <template <typename, typename, typename> class Policy>
inline bool BigFloat2::_atan_f(const BigFloat2& x, prec_t prec) {
  if (x.is_exact())
    return this->_atan<Policy, FT>(x.m_l, prec);
  else {
    Policy<FT, FT, FT>::atan(m_l, x.m_l, prec, BF_RNDD);
    Policy<FT, FT, FT>::atan(m_r, x.m_r, prec, BF_RNDU);
    set_exact(false);
  }
  return is_exact();
}
/// arctangent -- atan(T)
template <template <typename, typename, typename> class Policy, typename T>
inline bool BigFloat2::_atan(const T& x, prec_t prec) {
  set_exact(Policy<FT, T, FT>::atan(m_l, x, prec, BF_RNDD));
  if (!is_exact()) Policy<FT, T, FT>::atan(m_r, x, prec, BF_RNDU);
  return is_exact();
}
////////////////////////////////////////////////////////////////////////////////
/// log base e -- log(BigFloat2)
template <template <typename, typename, typename> class Policy>
inline bool BigFloat2::_log_f(const BigFloat2& x, prec_t prec) {
  if (x.is_exact())
    return this->_log<Policy, FT>(x.m_l, prec);
  else {
    Policy<FT, FT, FT>::log(m_l, x.m_l, prec, BF_RNDD);
    Policy<FT, FT, FT>::log(m_r, x.m_r, prec, BF_RNDU);
    set_exact(false);
  }
  return is_exact();
}
/// log base e -- log(T)
template <template <typename, typename, typename> class Policy, typename T>
inline bool BigFloat2::_log(const T& x, prec_t prec) {
  set_exact(Policy<FT, T, FT>::log(m_l, x, prec, BF_RNDD));
  if (!is_exact()) Policy<FT, T, FT>::log(m_r, x, prec, BF_RNDU);
  return is_exact();
}

////////////////////////////////////////////////////////////////////////////////
/// log base 2 -- log2(BigFloat2)
template <template <typename, typename, typename> class Policy>
inline bool BigFloat2::_log2_f(const BigFloat2& x, prec_t prec) {
  if (x.is_exact())
    return this->_log2<Policy, FT>(x.m_l, prec);
  else {
    Policy<FT, FT, FT>::log2(m_l, x.m_l, prec, BF_RNDD);
    Policy<FT, FT, FT>::log2(m_r, x.m_r, prec, BF_RNDU);
    set_exact(false);
  }
  return is_exact();
}
/// log base 2 -- log2(T)
template <template <typename, typename, typename> class Policy, typename T>
inline bool BigFloat2::_log2(const T& x, prec_t prec) {
  set_exact(Policy<FT, T, FT>::log2(m_l, x, prec, BF_RNDD));
  if (!is_exact()) Policy<FT, T, FT>::log2(m_r, x, prec, BF_RNDU);
  return is_exact();
}

////////////////////////////////////////////////////////////////////////////////
/// log base 10 -- log2(BigFloat2)
template <template <typename, typename, typename> class Policy>
inline bool BigFloat2::_log10_f(const BigFloat2& x, prec_t prec) {
  if (x.is_exact())
    return this->_log10<Policy, FT>(x.m_l, prec);
  else {
    Policy<FT, FT, FT>::log10(m_l, x.m_l, prec, BF_RNDD);
    Policy<FT, FT, FT>::log10(m_r, x.m_r, prec, BF_RNDU);
    set_exact(false);
  }
  return is_exact();
}
/// log base 2 -- log2(T)
template <template <typename, typename, typename> class Policy, typename T>
inline bool BigFloat2::_log10(const T& x, prec_t prec) {
  set_exact(Policy<FT, T, FT>::log10(m_l, x, prec, BF_RNDD));
  if (!is_exact()) Policy<FT, T, FT>::log10(m_r, x, prec, BF_RNDU);
  return is_exact();
}

////////////////////////////////////////////////////////////////////////////////
/// exponent -- exp(BigFloat2)
template <template <typename, typename, typename> class Policy>
inline bool BigFloat2::_exp_f(const BigFloat2& x, prec_t prec) {
  if (x.is_exact())
    return this->_exp<Policy, FT>(x.m_l, prec);
  else {
    Policy<FT, FT, FT>::exp(m_l, x.m_l, prec, BF_RNDD);
    Policy<FT, FT, FT>::exp(m_r, x.m_r, prec, BF_RNDU);
    set_exact(false);
  }
  return is_exact();
}
/// exponent -- exp(T)
template <template <typename, typename, typename> class Policy, typename T>
inline bool BigFloat2::_exp(const T& x, prec_t prec) {
  set_exact(Policy<FT, T, FT>::exp(m_l, x, prec, BF_RNDD));
  if (!is_exact()) Policy<FT, T, FT>::exp(m_r, x, prec, BF_RNDU);
  return is_exact();
}

////////////////////////////////////////////////////////////////////////////////
/// exponent of 2-- exp2(BigFloat2)
template <template <typename, typename, typename> class Policy>
inline bool BigFloat2::_exp2_f(const BigFloat2& x, prec_t prec) {
  if (x.is_exact())
    return this->_exp2<Policy, FT>(x.m_l, prec);
  else {
    Policy<FT, FT, FT>::exp2(m_l, x.m_l, prec, BF_RNDD);
    Policy<FT, FT, FT>::exp2(m_r, x.m_r, prec, BF_RNDU);
    set_exact(false);
  }
  return is_exact();
}
/// exponent of 2-- exp(T)
template <template <typename, typename, typename> class Policy, typename T>
inline bool BigFloat2::_exp2(const T& x, prec_t prec) {
  set_exact(Policy<FT, T, FT>::exp2(m_l, x, prec, BF_RNDD));
  if (!is_exact()) Policy<FT, T, FT>::exp2(m_r, x, prec, BF_RNDU);
  return is_exact();
}

////////////////////////////////////////////////////////////////////////////////
/// exponent of 10-- exp10(BigFloat2)
template <template <typename, typename, typename> class Policy>
inline bool BigFloat2::_exp10_f(const BigFloat2& x, prec_t prec) {
  if (x.is_exact())
    return this->_exp10<Policy, FT>(x.m_l, prec);
  else {
    Policy<FT, FT, FT>::exp10(m_l, x.m_l, prec, BF_RNDD);
    Policy<FT, FT, FT>::exp10(m_r, x.m_r, prec, BF_RNDU);
    set_exact(false);
  }
  return is_exact();
}
/// exponent of 2-- exp(T)
template <template <typename, typename, typename> class Policy, typename T>
inline bool BigFloat2::_exp10(const T& x, prec_t prec) {
  set_exact(Policy<FT, T, FT>::exp10(m_l, x, prec, BF_RNDD));
  if (!is_exact()) Policy<FT, T, FT>::exp10(m_r, x, prec, BF_RNDU);
  return is_exact();
}

////////////////////////////////////////////////////////////////////////////////
/// addition -- (BigFloat2 + BigFloat2)
template <template <typename, typename, typename> class Policy>
inline bool BigFloat2::_add_f(const BigFloat2& x, const BigFloat2& y, prec_t prec) {
  if (x.is_exact())
    return this->_add<Policy, FT>(x.m_l, y, prec);
  else if (y.is_exact())
    return this->_add<Policy, FT>(x, y.m_l, prec);
  else {
    Policy<FT, FT, FT>::add(m_l, x.m_l, y.m_l, prec, BF_RNDD);
    Policy<FT, FT, FT>::add(m_r, x.m_r, y.m_r, prec, BF_RNDU);
    assert(m_r!=m_l);
    set_exact(false);
    return is_exact();
  }
}
/// addition -- (BigFloat2 + T)
template <template <typename, typename, typename> class Policy, typename T>
inline bool BigFloat2::_add(const BigFloat2& x, const T& y, prec_t prec) {
  set_exact(Policy<FT, FT, T>::add(m_l, x.m_l, y, prec, BF_RNDD));
  if (!is_exact() || !x.is_exact()) {
    Policy<FT, FT, T>::add(m_r, x.m_r, y, prec, BF_RNDU);
    assert(m_r!=m_l);
    set_exact(false);
  }
  return is_exact();
}
/// addition -- (T + BigFloat2)
template <template <typename, typename, typename> class Policy, typename T>
inline bool BigFloat2::_add(const T& x, const BigFloat2& y, prec_t prec) {
  set_exact(Policy<FT, T, FT>::add(m_l, x, y.m_l, prec, BF_RNDD));
  if (!is_exact() || !y.is_exact()) {
    Policy<FT, T, FT>::add(m_r, x, y.getRight(), prec, BF_RNDU);
    assert(m_r!=m_l);
    set_exact(false);
  } 
  return is_exact();
}

////////////////////////////////////////////////////////////////////////////////
/// subtraction -- (BigFloat2 - BigFloat2)
template <template <typename, typename, typename> class Policy>
inline bool BigFloat2::_sub_f(const BigFloat2& x, const BigFloat2& y, prec_t prec) {
  if (x.is_exact())
    return this->_sub<Policy, FT>(x.m_l, y, prec);
  else if (y.is_exact())
    return this->_sub<Policy, FT>(x, y.m_l, prec);
  else {
    Policy<FT, FT, FT>::sub(m_l, x.m_l, y.m_r, prec, BF_RNDD);
    Policy<FT, FT, FT>::sub(m_r, x.m_r, y.m_l, prec, BF_RNDU);
    assert(m_r!=m_l);
    set_exact(false);
    return is_exact();
  }
}
/// subtraction -- (BigFloat2 - T)
template <template <typename, typename, typename> class Policy, typename T>
inline bool BigFloat2::_sub(const BigFloat2& x, const T& y, prec_t prec) {
  set_exact(Policy<FT, FT, T>::sub(m_l, x.m_l, y, prec, BF_RNDD));
  if (!is_exact() || !x.is_exact()) {
    Policy<FT, FT, T>::sub(m_r, x.getRight(), y, prec, BF_RNDU);
    assert(m_r!=m_l);
    set_exact(false);
  }
  return is_exact();
}
/// subtraction -- (T - BigFloat2)
template <template <typename, typename, typename> class Policy, typename T>
inline bool BigFloat2::_sub(const T& x, const BigFloat2& y, prec_t prec) {
  set_exact(Policy<FT, T, FT>::sub(m_l, x, y.getRight(), prec, BF_RNDD));
  if (!is_exact() || !y.is_exact()) {
    Policy<FT, T, FT>::sub(m_r, x, y.m_l, prec, BF_RNDU);
    assert(m_r!=m_l);
    set_exact(false);
  } 
  return is_exact();
}

////////////////////////////////////////////////////////////////////////////////
/// multiplication -- (BigFloat2 * BigFloat2)
template <template <typename, typename, typename> class Policy>
inline bool BigFloat2::_mul_f(const BigFloat2& x, const BigFloat2& y, prec_t prec) {
  if (x.is_exact())
    return this->_mul<Policy, FT>(x.m_l, y, prec);
  else if (y.is_exact())
    return this->_mul<Policy, FT>(x, y.m_l, prec);
  else {
    typedef Policy<FT, FT, FT> P;
    if (x.m_l.sgn() >= 0) {
      if (y.m_l.sgn() >= 0) {
        P::mul(m_l, x.m_l, y.m_l, prec, BF_RNDD);
        P::mul(m_r, x.m_r, y.m_r, prec, BF_RNDU);
      } else if (y.m_r.sgn() <= 0) {
        P::mul(m_l, x.m_r, y.m_l, prec, BF_RNDD);
        P::mul(m_r, x.m_l, y.m_r, prec, BF_RNDU);
      } else {
        P::mul(m_l, x.m_r, y.m_l, prec, BF_RNDD);
        P::mul(m_r, x.m_r, y.m_r, prec, BF_RNDU);
      }
    } else if (x.m_r.sgn() <= 0) {
      if (y.m_l.sgn() >= 0) {
        P::mul(m_l, x.m_l, y.m_r, prec, BF_RNDD);
        P::mul(m_r, x.m_r, y.m_l, prec, BF_RNDU);
      } else if (y.m_r.sgn() <= 0) {
        P::mul(m_l, x.m_r, y.m_r, prec, BF_RNDD);
        P::mul(m_r, x.m_l, y.m_l, prec, BF_RNDU);
      } else {
        P::mul(m_l, x.m_l, y.m_r, prec, BF_RNDD);
        P::mul(m_r, x.m_l, y.m_l, prec, BF_RNDU);
      }
    } else {
      if (y.m_l.sgn() >= 0) {
        P::mul(m_l, x.m_l, y.m_r, prec, BF_RNDD);
        P::mul(m_r, x.m_r, y.m_r, prec, BF_RNDU);
      } else if (y.m_r.sgn() <= 0) {
        P::mul(m_l, x.m_r, y.m_l, prec, BF_RNDD);
        P::mul(m_r, x.m_l, y.m_l, prec, BF_RNDU);
      } else {
        FT tmp;
        // compute min{x.m_l*y.m_r, x.m_r*y.m_l}
        P::mul(m_l, x.m_l, y.m_r, prec, BF_RNDD);
        tmp.set_prec(m_l.get_prec());
        P::mul(tmp, x.m_r, y.m_l, prec, BF_RNDD);
        if (m_l.cmp(tmp) > 0) m_l.swap(tmp);
        // compute max{x.m_r*y.m_r, x.m_l*y.m_l}
        P::mul(m_r, x.m_r, y.m_r, prec, BF_RNDU);
        tmp.set_prec(m_r.get_prec());
        P::mul(tmp, x.m_l, y.m_l, prec, BF_RNDU);
        if (m_r.cmp(tmp) < 0) m_r.swap(tmp);
      }
    }
    assert(m_r!=m_l);
    set_exact(false);
    return is_exact();
  }
}
/// multiplication -- (BigFloat2 * T)
template <template <typename, typename, typename> class Policy, typename T>
inline bool BigFloat2::_mul(const BigFloat2& x, const T& y, prec_t prec) {
  typedef Policy<FT, FT, T> P;
  if (x.is_exact()) {
    set_exact(P::mul(m_l, x.m_l, y, prec, BF_RNDD));
    if (!is_exact()) P::mul(m_r, x.m_l, y, prec, BF_RNDU);
  } else {
    P::mul(m_l, (y>0?x.m_l:x.m_r), y, prec, BF_RNDD);
    P::mul(m_r, (y>0?x.m_r:x.m_l), y, prec, BF_RNDU);
    set_exact(m_l==m_r);
  } 
  return is_exact();
}
/// multiplication -- (T * BigFloat2)
template <template <typename, typename, typename> class Policy, typename T>
inline bool BigFloat2::_mul(const T& x, const BigFloat2& y, prec_t prec) {
  typedef Policy<FT, T, FT> P;
  if (y.is_exact()) {
    set_exact(P::mul(m_l, x, y.m_l, prec, BF_RNDD));
    if (!is_exact()) P::mul(m_r, x, y.m_l, prec, BF_RNDU);
  } else {
    P::mul(m_l, x, (x>0?y.m_l:y.m_r), prec, BF_RNDD);
    P::mul(m_r, x, (x>0?y.m_r:y.m_l), prec, BF_RNDU);
    set_exact(m_l==m_r);
  } 
  return is_exact();
}

////////////////////////////////////////////////////////////////////////////////
/// division -- (BigFloat2 / BigFloat2)
template <template <typename, typename, typename> class Policy>
inline bool BigFloat2::_div_f(const BigFloat2& x, const BigFloat2& y, prec_t prec) {
  if (x.is_exact())
    return this->_div<Policy, FT>(x.m_l, y, prec);
  else if (y.is_exact()) {
	  return this->_div<Policy, FT>(x, y.m_l, prec);; }
  else if (y.has_zero()) {
    set_inf();
  } else {
   // Jihun, July 2010
   // Heck, y is corrupted sometime
   // y.m_r < y.m_l
   // need to correct this
   BigFloat2 z(y); z.validate();
   typedef Policy<FT, FT, FT> P;
    if (x.m_l.sgn() >= 0) {
      if (z.m_l.sgn() >= 0) {
        P::div(m_l, x.m_l, z.m_r, prec, BF_RNDD);
        P::div(m_r, x.m_r, z.m_l, prec, BF_RNDU);
      } else if (z.m_r.sgn() <= 0) {
        P::div(m_l, x.m_r, z.m_r, prec, BF_RNDD);
        P::div(m_r, x.m_l, z.m_l, prec, BF_RNDU);
      }
    } else if (x.m_r.sgn() <= 0) {
      if (z.m_l.sgn() >= 0) {
        P::div(m_l, x.m_l, z.m_l, prec, BF_RNDD);
        P::div(m_r, x.m_r, z.m_r, prec, BF_RNDU);
      } else if (z.m_r.sgn() <= 0) {
        P::div(m_l, x.m_r, z.m_l, prec, BF_RNDD);
        P::div(m_r, x.m_l, z.m_r, prec, BF_RNDU);
      }
    } else {
      if (z.m_l.sgn() > 0) {
        P::div(m_l, x.m_l, z.m_l, prec, BF_RNDD);
        P::div(m_r, x.m_r, z.m_l, prec, BF_RNDU);
      } else if (z.m_r.sgn() < 0) {
        P::div(m_l, x.m_r, z.m_r, prec, BF_RNDD);
        P::div(m_r, x.m_l, z.m_r, prec, BF_RNDU);
      }
    }
    assert(m_r!=m_l);
    set_exact(false);
  }
  return is_exact();
}
/// division -- (BigFloat2 / T)
template <template <typename, typename, typename> class Policy, typename T>
inline bool BigFloat2::_div(const BigFloat2& x, const T& y, prec_t prec) {
  typedef Policy<FT, FT, T> P;
  if (x.is_exact()) {
    set_exact(P::div(m_l, x.m_l, y, prec, BF_RNDD));
    if (!is_exact()) P::div(m_r, x.m_l, y, prec, BF_RNDU);
  } else {
    P::div(m_l, (y>0?x.m_l:x.m_r), y, prec, BF_RNDD);
    P::div(m_r, (y>0?x.m_r:x.m_l), y, prec, BF_RNDU);
    assert(m_r!=m_l);
    set_exact(false);
  }
  return is_exact();
}
/// division -- (T / BigFloat2)
template <template <typename, typename, typename> class Policy, typename T>
inline bool BigFloat2::_div(const T& x, const BigFloat2& y, prec_t prec) {
  typedef Policy<FT, T, FT> P;
  if (y.is_exact()) {
    set_exact(P::div(m_l, x, y.m_l, prec, BF_RNDD));
    if (!is_exact()) P::div(m_r, x, y.m_l, prec, BF_RNDU);
  } else {
    P::div(m_l, x, (x>0?y.m_r:y.m_l), prec, BF_RNDD);
    P::div(m_r, x, (x>0?y.m_l:y.m_r), prec, BF_RNDU);
    set_exact(m_l==m_r); // for the case of x=0;
  } 
  return is_exact();
}
