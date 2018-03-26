/****************************************************************************
 * VExpr.h -- EGC number class providing guarranteed precision
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
 * $Id: VExpr.h,v 1.1 2010/01/14 13:53:50 exact Exp $
 ***************************************************************************/
#ifndef __CORE_EXPR_H__
#define __CORE_EXPR_H__

#include <CORE/VExprRep.h>

CORE_BEGIN_NAMESPACE

/// \class ExprT
/// Kernel -- internal representation 
template <typename RootBd, typename Filter, typename Kernel>
class VExprT {
public: // public typedefs
  typedef RootBd     RootBdT;
  typedef Filter     FilterT;
  typedef Kernel     KernelT;

  typedef typename Kernel::ZT ZT;
  typedef typename Kernel::QT QT;
  typedef typename Kernel::FT FT;
  typedef Kernel              KT;

private: // private typedefs
  typedef InterfaceRepT<RootBd, Filter, Kernel> InterfaceRep;

public:
 /// helper function for constructing Polynomial node witb BFInterval
  template <class NT>
  friend VExprT rootOf(const BiPoly<NT>& f, const BiPoly<NT>& g) {
    return VExprT(new InterfaceRepT<RootBd, Filter, Kernel, NT>(f, g));
  }
  /// return internal rep
  VExprRep* rep() const
  { return m_rep; }

  VExprRep* m_rep; ///<- internal representation
}; // end if ExprT

CORE_END_NAMESPACE

///////////////////////////////////////////////////////////////////////////
// Definition of VExpr
///////////////////////////////////////////////////////////////////////////

#include <CORE/RootBounds.h>
#include <CORE/Filters.h>

CORE_BEGIN_NAMESPACE

// BFMSS-Kary root bound + BFS filter + BigFloat2
typedef VExprT<BfmsskRootBd_BigFloat<BigFloat2>, BfsFilter<BigFloat2>, BigFloat2> VExpr;
// BFMSS-Kary root bound + BFS filter + BigFloat2
//typedef ExprT<BfmsskRootBd<BigFloat2>, BfsFilter<BigFloat2>, BigFloat2> VExpr;
// BFMSS root bound + BFS filter + BigFloat2
typedef VExprT<BfmssRootBd<BigFloat2>, BfsFilter<BigFloat2>, BigFloat2> VExpr_old;
// BFMSS root bound + Dummy filter + BigFloat2
//typedef ExprT<BfmssRootBd<BigFloat2>, DummyFilter, BigFloat2> Expr;

// Dummy root bound + Dummy filter + BigFloat2
//typedef ExprT<BfmssRootBd<BigFloat2>, DummyRootBd<10>, BigFloat> Expr;


CORE_END_NAMESPACE
#endif // __CORE_EXPR_H__
