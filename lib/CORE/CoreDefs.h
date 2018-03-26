/****************************************************************************
 * CoreDefs.h -- Global parameters for the Core Library
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
 * $Id: CoreDefs.h,v 1.29 2010/11/23 17:58:36 exact Exp $
 ***************************************************************************/
#ifndef __CORE_COREDEFS_H__
#define __CORE_COREDEFS_H__
#include <climits>
#include <string>
#include <stack>

CORE_BEGIN_NAMESPACE

#ifndef sign_t
typedef int sign_t;
#endif

#ifndef msb_t
typedef long msb_t;
#endif
  
#define MSB_MIN LONG_MIN
#define MSB_MAX LONG_MAX
  
#define LOG2_10 3.3219280948873626 /* upper approximation of log(10)/log(2) */
#define LOG10_2 0.30102999566398120 /* upper approximation of log(2)/log(10) */
  
#define CORE_INFTY LONG_MAX
#define CORE_posInfty LONG_MAX
#define CORE_negInfty LONG_MIN

extern long defAbsPrec;
extern long defRelPrec;
extern long defInputDigits;
extern long defOutputDigits;
extern long defBFdivRelPrec;
extern long defBFradicalRelPrec;   // this replaces defBFsqrtAbsPrec in Core1
                                   // Should also be used for any rootOf()
extern long defWBFPrec;		   // weak bigfloat precision
extern bool WBFenabled;
extern unsigned long cutOffBound;  // arbitary cutoff for ABSOLUTE precision
extern unsigned long escapeBound;  // this is to "escape" in transcendental
                                   //  evaluation when we have no root bounds
extern bool coretest_error;	   // error flag for unit test
extern bool coretest_verbose;	   // verbose flag for unit test
extern int  coreErrorFlag;		// This is 0 if no error,
								// core_error() will increment this whenever it
								// called with err=true.
								// Users can reset this and test for this value.
extern std::string coreErrorMsg; // Message that will be printed by core_error()
								// User can set this message before a suspected error.

extern bool rationalReduceFlag;  // We bring this back from Core1. We make rational nodes as leaf nodes by computing expression


inline unsigned long digits2bits(unsigned long digits)
{
  // Jan 2010 (Jihun) bug fix: 
  // Added test for CORE_INFTY because digits = CORE_INFTY is
  //     interpreted to mean that we want a rational representation!
  //     This digits is also converted to bits (=prec) which is
  //     also checked for CORE_INFTY when converting from string to 
  //     expression value (i.e., BigFloat if prec<CORE_INFTY, or BigRat if
  //     prec=CORE_INFTY).
  //
  // We do the following min(x,y) on a pair of double values,
  //    so that if if digits=CORE_INFTY, we will return CORE_INFTY
  if (digits == CORE_INFTY) return CORE_INFTY;
  return (unsigned long)(std::min)(digits*LOG2_10, (double)CORE_INFTY);
}

inline unsigned long bits2digits(unsigned long bits)
{
  if(bits == CORE_INFTY) return CORE_INFTY; // added: Jan 2010
  return (unsigned long)(bits*LOG10_2);
}


extern void core_error(std::string msg, std::string file, int lineno, bool err);

/// This sets the global variable defRelPrec and defAbsPrec.
//  PROBLEM IS this is "composite precision" which we don't really
//  support in Core2.  So, one of these two must be CORE_INFTY.
//  But in Expr.h, we see that approx() is implemented so that
//     If defAbsPrec<CORE_INFTY, then we compute to defAbsPrec (absolute prec)
//     Else defAbsPrec=CORE_INFTY, we compute to defRelPrec (relative prec).
inline void setDefaultPrecision(long r, long a)
{ defRelPrec = r; defAbsPrec = a; }
inline void setDefaultRelPrecision(long r) {
  defRelPrec = r; 
  defAbsPrec = CORE_INFTY;
  if (defRelPrec == CORE_INFTY)
    core_error("Relative Prec and Absolute Prec are both CORE_INFTY", __FILE__, __LINE__, false);
}
inline void setDefaultAbsPrecision(long a) {
  defAbsPrec = a; 
  defRelPrec = CORE_INFTY;
  if (defAbsPrec == CORE_INFTY)
    core_error("Relative Prec and Absolute Prec are both CORE_INFTY", __FILE__, __LINE__, false);
}
inline void setDefaultComPrecision(long r, long a) {
  defAbsPrec = a;
  defRelPrec = r;
  if (defAbsPrec == CORE_INFTY && defRelPrec == CORE_INFTY)
    core_error("Relative Prec and Absolute Prec are both CORE_INFTY", __FILE__, __LINE__, false);
}

// THIS controls the precision for converting input values into Core internal representation
//    -- if this value is CORE_INFTY, then we must not have any error in the
//       the internal representation.  This may require a BigRat value
inline long getDefaultInputDigits()
{ return defInputDigits; }
inline void setDefaultInputDigits(long digits)
{ defInputDigits = digits; }

inline long get_def_input_digits()
{ return defInputDigits; }
inline bool is_infty(long l)
{ return l == CORE_INFTY; }

inline long getDefaultBFdivPrec()
{ return defBFdivRelPrec; }
inline void setDefaultBFdivPrec(long p)
{ defBFdivRelPrec= p; }

inline long getDefaultBFradicalPrec()
{ return defBFradicalRelPrec; }
inline void setDefaultBFradicalPrec(long p)
{ defBFradicalRelPrec= p; }

/// get bound (in bits) for global absolute error iteration
/// this bound is used in refine()
inline unsigned long get_cut_off_bound()
{ return cutOffBound; }
inline unsigned long set_cut_off_bound(unsigned long p) {
  unsigned long ret = cutOffBound;
  cutOffBound = p; 
  return ret;
}

inline bool get_rational_reduce_flag() {
  return rationalReduceFlag;
}

inline bool set_rational_reduce_flag(const bool b) {
  bool ret = rationalReduceFlag;
  rationalReduceFlag = b;
  return ret;
}

/// escape bound is used as root bound for transcendental nodes
inline unsigned long get_escape_bound()
{ return escapeBound; }
inline unsigned long set_escape_bound(unsigned long p) {
  unsigned long ret = escapeBound;
  escapeBound = p; 
  return ret;
}

/// COMMENTS NEEDED...
// Jihun (Oct 2007)
inline void set_wbf_mode(bool bEnable) {
  WBFenabled = bEnable;
}

/// COMMENTS NEEDED...
// Jihun (Oct 2007)
inline bool get_wbf_mode() {
  return WBFenabled;
}

// Temporary test. Weak BigFloat precision control functions.
inline void set_wbf_prec(unsigned long p) {
  setDefaultBFdivPrec(p);
  defWBFPrec = p;
}

inline unsigned long get_wbf_prec() {
  return defWBFPrec;
}

// An enum type for "Sign Conditions" -- we should use
//  this more systematically throughout Core.
//  For instance, many of our functions that return a Boolean value
//  should actually return a SignCondition.   E.g., zero() should
//  return _ZERO or _NONZERO (or perhaps _UNKNOWN).   --Chee 9/13/2010

enum SignCondition { _POS, _NEG, _ZERO, _NONNEG, _NONPOS, _NONZERO, _UNKNOWN };

// code for managing core Environment. Haag 10/08/2008

enum GlobalParamType { _WBFenabled, _defWBFPrec};

class coreEnv {

  bool setEnv_WBFenabled;
  bool Env_WBFenabled;
  bool prevEnv_WBFenabled;

  bool setEnv_defWBFPrec;
  long Env_defWBFPrec;
  long prevEnv_defWBFPrec;

public:
  coreEnv();
  void Env_set_wbf_mode( bool v) { setEnv_WBFenabled = true; Env_WBFenabled = v;}
  bool Env_get_wbf_mode() { return  Env_WBFenabled;}
  void Env_unset_wbf_mode() { setEnv_WBFenabled = false;}
  bool Env_test_wbf_mode() { return setEnv_WBFenabled;}
  bool prevEnv_get_wbf_mode() { return  prevEnv_WBFenabled;}

  void Env_set_wbf_prec( long v) { setEnv_defWBFPrec = true; Env_defWBFPrec = v;}
  long Env_get_wbf_prec() { return Env_defWBFPrec; }
  void Env_unset_wbf_prec() { setEnv_defWBFPrec = false;}
  bool Env_test_wbf_prec() { return setEnv_defWBFPrec;}
  long prevEnv_get_wbf_prec() { return prevEnv_defWBFPrec; }

  long apply();
};


long restoreEnv(); 

// End core Environment code



CORE_END_NAMESPACE

#endif /*__CORE_COREDEFS_H__*/
