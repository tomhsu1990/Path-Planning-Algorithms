/****************************************************************************
 * MpfrIO.cpp -- I/O functions for MPFR
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
 * $Id: MpfrIO.cpp,v 1.12 2009/02/23 19:34:10 exact Exp $
 ***************************************************************************/
#include <mpfr.h>
#include <string>
#include <cstring>
#include <iostream>
/* Define BITS_PER_MP_LIMB
   Can't use sizeof(mp_limb_t) since it should be a preprocessor constant */
#if defined(GMP_NUMB_BITS) /* GMP 4.1.2 or above */
# define BITS_PER_MP_LIMB  (GMP_NUMB_BITS+GMP_NAIL_BITS)
#elif defined (__GMP_BITS_PER_MP_LIMB) /* Older versions 4.x.x */
# define BITS_PER_MP_LIMB  __GMP_BITS_PER_MP_LIMB
# define GMP_NUMB_BITS BITS_PER_MP_LIMB
# ifndef GMP_NAIL_BITS
#  define GMP_NAIL_BITS 0
# endif
#else
# error "Could not detect BITS_PER_MP_LIMB. Try with gmp internal files."
#endif
  
#define MPFR_MANT(x)      ((x)->_mpfr_d)
#define MPFR_PREC(x)      ((x)->_mpfr_prec)
#define MPFR_LIMB_SIZE(x) ((MPFR_PREC((x))-1)/BITS_PER_MP_LIMB+1)

// remove trailing zeros (only by limbs)
void mpfr_remove_trailing_zeros(mpfr_t x) {
  unsigned int xn = MPFR_LIMB_SIZE(x); 
  mp_limb_t* xp = MPFR_MANT(x);
  unsigned int i = 0;
  while (i < xn && xp[i] == 0) i++;
  if (i > 0 && i < xn) mpfr_round_prec(x, GMP_RNDN, (xn-i)*BITS_PER_MP_LIMB); 
}

// Jun. 1 2017: extract has a bug cannot be built
/*
// use mpf to read mpfr from istream
extern std::istream& extract(std::istream& i, mpf_t x);
std::istream& operator>>(std::istream& is, mpfr_ptr f) {
  mpf_t tmp; 
  mpf_init2(tmp, mpfr_get_prec(f));
  ::extract(is,tmp);
  mpfr_set_f(f, tmp, __gmp_default_rounding_mode);
  mpf_clear(tmp); 
  return is;
}
*/

// convert integer to string
static char* itoa(int val, int base) {
  static char buf[32] = {0};
  int i = 30;
  for (; val && i; --i, val /= base)
    buf[i] = "0123456789abcdef"[val % base];
  return &buf[i+1];
}

/* _mpfr_alloc_cstr */
struct _mpfr_alloc_cstr {
  char *str;
  _mpfr_alloc_cstr() : str(0) { }
  void alloc(int len) { delete[] str; str = new char[len]; }
  ~_mpfr_alloc_cstr() { delete[] str; }
};

// IO format:
//
// 1. we use exact precision bits returned from cout.precision(), while
//    C++ IO use it as bits after decimal point.
// 2. we use exact digits to output exponents, while C++ IO use 2 or more 
//    digits.

// convert mpfr to string
std::string mpfr2str(
  mpfr_srcptr mp,   // mpfr
  size_t ndigits,   // number of digits
  int base,         // base
  int fmt,          // format
  mp_rnd_t rnd,     // round mode
  bool showpoint,   // show decimal point
  bool showpos,     // show '+' sign
  bool uppercase    // show in uppercase
) {
using namespace std;
  // check zero
  bool is_zero = mpfr_zero_p(mp) != 0;
  // if zero and fmt = 0
  if (is_zero && fmt == 0) return std::string("0");

  mp_exp_t exp;
  char *str;
  _mpfr_alloc_cstr tmp;
  unsigned long len;

  // get the digits 
  if (ndigits == 0) { // impossible to predetermine the size of string
    str = mpfr_get_str(0, &exp, base, ndigits, mp, rnd);
  } else {
    len = (std::max)(ndigits + 2UL, 7UL);
    tmp.alloc(len);
    str = mpfr_get_str(tmp.str, &exp, base, (std::max)(ndigits, size_t(2)), mp, rnd);
  }
  
  // if failed, return empty string
  if (str == 0) return std::string();

  // find the start position implicit radix point
  std::string::size_type first = (str[0] == '-') ? 1 : 0;
  
  // choose format by value
  len = strlen(str) - first;
  if (fmt == 0) {
    unsigned long e = (exp >= 0) ? exp : (-exp);
    if (e > len) // too small/big
      fmt = 2;
    else
      fmt = 1; 
  }
  
  // now formatting
  std::string result;

  if (fmt == 1) { // fixed format
    // counting trailing zeros
    while (len > 1 && str[len+first-1] == '0')
      len --;

    result.assign(str, len+first);
    if (exp > 0) {  
      // BUG fixed: if ((unsigned long)exp > len-first) (12/14/06)
      //   result.append(exp - len+first, '0');
      if ((unsigned long)exp > len) // integer need padding 0
        result.append(exp - len, '0');
      else if ((unsigned long)exp < len) // float point value
        result.insert(first + exp, ".");
    } else if (exp < 0) {
      result.insert(first, -exp, '0');
      result.insert(first, ".");
    } else
      result.insert(first, "0.");
  } else { // scientific format
    // counting trailing zeros
    while (len > 1 && str[len+first-1] == '0')
      len --;

    result.assign(str, len+first);
    // put the decimal before the second digit
    result.insert(first+1, ".");
    // decrease the exponent by 1
    if (!is_zero) --exp;
    // output exponent
    if (base <= 10)
      result.append(uppercase ? "E" : "e");
    else
      result.append("@");
    if (exp == 0) {
      result.append("+0");
    } else if (exp > 0) {
      result.append("+"); result.append(itoa(exp, base));
    } else {
      result.append("-"); result.append(itoa(-exp, base));
    }
  }
  // insert "+" if show pos
  if (first == 0 && showpos) result.insert(0, "+");

  // free memory allocated by mpfr_get_str
  if (ndigits == 0) mpfr_free_str(str);

  return result;
}

