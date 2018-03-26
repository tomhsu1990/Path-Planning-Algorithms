/****************************************************************************
 * CORE.h -- The main inclusion file for the Core Library
 *
 * 	$Id: CORE.h,v 1.25 2010/07/21 16:18:08 exact Exp $
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
 ***************************************************************************/
#ifndef __CORE_CORE_H__
#define __CORE_CORE_H__

#ifndef CORE_LEVEL
  #define CORE_LEVEL 3
#endif

#include <CORE/Config.h>
#include <CORE/CoreDefs.h>

// level 1 with Wrappers. 
// If Wrappers are to be disabled the
// CORE_LEVEL_1_NO_WRAPPERS must be defined.
#ifndef CORE_LEVEL_1_NO_WRAPPERS
#include <CORE/Wrappers.h>
#endif

// REMARK: why do we include all these files 
// 	if CORE_LEVEL=1?
//
// level 1
#include <CORE/BigInt.h>
#include <CORE/BigRat.h>
#include <CORE/BigFloat.h>
#include <CORE/BFInterval.h>
// level 2
#include <CORE/BigFloat2.h>
// level 3
#include <CORE/Expr.h>

// include inline functions for ExprRep
#include <CORE/ExprRep.inl>

// timer
#include <CORE/Timer.h>

typedef double machine_double;
typedef long machine_long;

#if CORE_LEVEL == 1
  // LEVEL_1 with wrappers provides wrappers over fundamental
  // data types, and exposes the commonly used part of the  API
  // that the BigInt / BigFloat and BigRat classes expose. Define
  // this quantity if for some reason you would not like to use
  // wrappers at Level 1.
  #ifndef CORE_LEVEL_1_NO_WRAPPERS
    // LongWrappers are disabled for now because the implementation
    // is incomplete. So BigInts become longs.
    //
    // NOTE (Narayan): On many platforms, sizeof(long) == sizeof(int)
    // so we need "long long", available as int64_t on most unix
    // variants.
    //
    // #define BigInt LongWrapper
    // #define long LongWrapper

    #define BigInt long
    #define BigFloat DoubleWrapper
    #define BigRat DoubleWrapper
    #define Expr DoubleWrapper
    #define CORE::double DoubleWrapper
  #else
    // Use plain simple datatypes.
    #define BigInt long
    #define BigFloat double
    #define BigRat double
    #define Expr double
  #endif
#elif CORE_LEVEL == 2
  #undef CORE::double
  #define CORE::double BigFloat
  #undef CORE::long
  #define CORE::long BigInt
  #define Expr BigFloat
#elif CORE_LEVEL == 3
  #undef CORE::double
  #define CORE::double Expr
  #undef CORE::long
  #define CORE::long BigInt
#elif CORE_LEVEL == 5	// this is a hack to allow users (Willi Mann)
			// to get "Level 3 doubles" (i.e., Expr)
			// and at the same time use "unsigned long" with no error.
  #undef CORE::double
  #define CORE::double Expr
#endif

//#ifndef CORE_NO_AUTOMATIC_NAMESPACE
//  using namespace CORE;
//#endif

#endif /*__CORE_CORE_H__*/
