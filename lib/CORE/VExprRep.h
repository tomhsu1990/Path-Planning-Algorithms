/****************************************************************************
 * VExprRep.h -- Internal Representation for Expr
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
 * $Id: VExprRep.h,v 1.1 2010/01/14 13:53:50 exact Exp $
 ***************************************************************************/
#ifndef __CORE_EXPRREP_H__
#define __CORE_EXPRREP_H__

#include <CORE/CoreAux.h>
#include <CORE/poly/Curves.h>
#include <bitset>
#include <iostream>
#include <sstream>

CORE_BEGIN_NAMESPACE

/// \class ExprRepT
/// \brief represent an expression node in DAG
template <typename RootBd, typename Filter, typename Kernel>
class VExprRepT {
  typedef VExprRepT<RootBd, Filter, Kernel> thisClass;
protected:
  VExprRepT() : m_ref_counter(1)
  {} 
  virtual ~VExprRepT()
  {}
public: // reference counting
  void inc_ref() 
  { ++m_ref_counter; }
  void dec_ref()
  { if (--m_ref_counter == 0) delete this; }
  const int get_ref() const
  { return m_ref_counter; }
private:
  int m_ref_counter;
};

/// \class InterfaceRepT 
/// \brief Interface node
template <typename RootBd, typename Filter, typename Kernel,typename NT>
class InterfaceRepT : public VExprRepT<RootBd, Filter, Kernel> {
  typedef VExprRepT<RootBd, Filter, Kernel> VExprRep;
protected:
  Kernel x,y;
public:
  InterfaceRepT(const BiPoly<NT>& f, const BiPoly<NT>& g) {
    x = 1; y = 0;
  }

  virtual ~InterfaceRepT() 
  {}

protected:
};

CORE_END_NAMESPACE

#endif /*__CORE_EXPRREP_H__*/
