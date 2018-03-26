/****************************************************************************
 * Promote.h -- Templates providing type promotion
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
 * $Id: Promote.h,v 1.8 2010/09/23 20:06:22 exact Exp $
 ***************************************************************************/
#ifndef __CORE_PROMOTE_H__
#define __CORE_PROMOTE_H__

#include <CORE/Config.h>

CORE_BEGIN_NAMESPACE

/// hasExactDivision()
///   CHECKING if NT has exact division
///   NOTE: it is important that the compiler does not try to
///   prove theorems about arithmetic identities like "x*(y/x) == y"
///   USAGE:  If you want to check if a number type NT has exact division, 
///           do for example,
///            if (hasExactDivision< NT >::check()) ...
///   		We use this in Polynomial<NT> class.
template < class NT >
struct hasExactDivision {
// This default function is supposed to work for NT other than BigRat or Expr
  static bool check() {	
     return false;
  }
};

template<> struct hasExactDivision<Expr> {
  static bool check() {
     return true;
  }
};
template<> struct hasExactDivision<BigRat> {
  static bool check() {
     return false;
  }
};

/// hasStringConstructor()
///   CHECKING if NT has string constructor
///   NOTE: it is important that the compiler does not try to
///   convert string to NT
///   USAGE:  If you want to check if a number type NT has string converter, 
///           do for example,
///            if (hasStringConstructor< NT >::check()) ...
///   		We use this in Curves.tcc class.
template <class NT>
struct hasStringConstructor {
  static bool check() {
    return true;
  }
  static NT convert(const char* str) {
    return NT(str);
  }
};

template <> struct hasStringConstructor<int> {
  static bool check() {
    return false;
  }
  static int convert(const char* str) {
    return atoi(str);
  } 
};

template <> struct hasStringConstructor<long> {
  static bool check() {
    return false;
  }
  static long convert(const char* str) {
    return atol(str);
  } 
};

template <> struct hasStringConstructor<double> {
  static bool check() {
    return false;
  } 
   static double convert(const char* str) {
    return atof(str);
  }
};



template<typename T1, typename T2>
class Promotion;

template<typename T>
class Promotion<T, T> {
  public:
    typedef T ResultT;
};

#define MAX_TYPE(T1, T2)		\
  typename Promotion<T1, T2>::ResultT

#define DEFINE_MAX_TYPE(T1, T2, Tr)	\
  template<> class Promotion<T1, T2> {	\
    public:				\
      typedef Tr ResultT;		\
  };					\
  template<> class Promotion<T2, T1> {	\
    public:				\
      typedef Tr ResultT;		\
  };

/*
 * For example:
 *
 * DEFINE_MAX_TYPE(BigInt, BigRat, BigRat)   	// define the promotion
 *
 * template<typename T1, typename T2> 		// define function f with type templates
 *   MAX_TYPE(T1, T2) f(T1& , T2& );
 *
 * or
 *
 * template<typename T1, typename T2> 		// define function f with type templates
 *   const MAX_TYPE(T1, T2)& f(T1& , T2& );
 *
 * BigInt  a  =  1;
 * BigRat  b  = "1/3";
 * BigRat  c  =  f(a, b);			// or, typename Promotion<BigInt, BigRat>::ResultT c = f(a,b);
 *
 * REMARK: this mechanism is used by the eval function for polynomial evaluation (see Poly.tcc)
 * where the two types are NT (type of coefficients) and N (type of evaluation point).
 */

/* 
 * primary types: (11)
 *
 * 	bool, 
 *	char, unsigned char, 
 *	short, unsigned short,
 * 	int, unsigned int, 
 *	long, unsigned long, 
 *	float, double
 *
 * CORE types: (5)
 *
 * 	BigInt < BigFloat < BigRat < Real < Expr
 *
 *      (NOTE: BigFloat here must be error-free)
 *
 */

class BigInt;
class BigFloat;
class BigRat;

DEFINE_MAX_TYPE(long, BigInt, BigInt)
DEFINE_MAX_TYPE(long, BigFloat, BigFloat)
DEFINE_MAX_TYPE(long, BigRat, BigRat)
DEFINE_MAX_TYPE(long, Expr, Expr)

DEFINE_MAX_TYPE(int, BigInt, BigInt)
DEFINE_MAX_TYPE(int, BigFloat, BigFloat)
DEFINE_MAX_TYPE(int, BigRat, BigRat)
DEFINE_MAX_TYPE(int, Expr, Expr)

DEFINE_MAX_TYPE(BigInt, BigFloat, BigFloat)
DEFINE_MAX_TYPE(BigInt, BigRat, BigRat)
DEFINE_MAX_TYPE(BigInt, Expr, Expr)

DEFINE_MAX_TYPE(BigFloat, BigRat, BigRat)
DEFINE_MAX_TYPE(BigFloat, Expr, Expr)

DEFINE_MAX_TYPE(BigRat, Expr, Expr)

CORE_END_NAMESPACE

#endif /*__CORE_PROMOTE_H__*/
