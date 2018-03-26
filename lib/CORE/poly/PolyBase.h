/****************************************************************************
 * PolyBase.h -- Polynomial Base Class
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
 * $Id: PolyBase.h,v 1.9 2010/11/23 18:00:52 exact Exp $
 ***************************************************************************/
#ifndef __CORE_POLYBASE_H__
#define __CORE_POLYBASE_H__

#include <CORE/Config.h>
#include <vector>

#ifndef CORE_DISABLE_REFCOUNTING
  #include <CORE/RefCounting.h>
#endif

CORE_BEGIN_NAMESPACE

/// \class PolyBase PolyBase.h
/// \brief PolyBase is the base class of polynomial class
template <typename NT>
class PolyBase
#ifndef CORE_DISABLE_REFCOUNTING
  : public RcRepImpl<PolyBase<NT> >
#endif
{
public:
  typedef std::vector<NT> VecNT;
  /// \name constructors, destructor and assignment operator
  //@{
  /// default constructor
  PolyBase() : _deg(-1), _coeff(0) {}
  /// copy constructor
  PolyBase(const PolyBase& rhs) : _deg(rhs._deg), _coeff(0) {
    if (_deg >= 0) {
      _coeff = new NT[_deg+1];
      for (int i=0; i<=_deg; ++i) _coeff[i] = rhs._coeff[i];
    }
  }
  /// copy constructor
  template <class T>
  PolyBase(const PolyBase<T>& rhs) : _deg(rhs.degree()), _coeff(0) {
    if (_deg >= 0) {
      _coeff = new NT[_deg+1];
      for (int i=0; i<=_deg; ++i) _coeff[i] = rhs.coeff()[i];
    }
  }
  /// destructor
  ~PolyBase()
  { if (_deg >= 0) delete[] _coeff; }

  /// constructor for unit polynomial with nominal degree n
  PolyBase(int n) : _deg(n), _coeff(0) {
    if (_deg >= 0) {
      _coeff = new NT[_deg+1]; _coeff[0] = 1;
      for (int i=1; i<=_deg; ++i) _coeff[i] = 0;
    }
  }
  /// constructor with coeff array
  PolyBase(int n, NT* coef) : _deg(n), _coeff(0) {
    assert(coef != 0);
    if (_deg >= 0) {
      _coeff = new NT[_deg+1];
      for (int i=0; i<=_deg; ++i) _coeff[i] = coef[i];
    }
  }
  /// constructor with coeff vector
  PolyBase(int n, const VecNT& coef) : _deg(coef.size()-1), _coeff(0) {
    if (_deg >= 0) {
      _coeff = new NT[_deg+1];
      for (int i=0; i<=_deg; ++i) _coeff[i] = coef[i];
    }
  }
  PolyBase(int n, const char* s[]) : _deg(n), _coeff(0) {
    if (_deg >= 0) {
      _coeff = new NT[_deg+1];
      for (int i=0; i<=_deg; ++i) _coeff[i] = s[i];
    }
  }
  /// constructor from <tt>char*</tt> (no implicit conversion)
  PolyBase(const char* s, char c) { std::cout << "error";}

  /// assignment operator for <tt>PolyBase</tt>
  PolyBase& operator=(const PolyBase& rhs) {
    if (this != &rhs) {
      _deg = rhs._deg;
      if (_deg >= 0) {
 	 delete[] _coeff;
        _coeff = new NT[_deg+1];
        for (int i=0; i<=_deg; ++i) _coeff[i] = rhs._coeff[i];
      }
    }
    return *this;
  }
  /// assignment function for coeff array
  // REMARK: this function should probably be called "setCoeffs"
  // (contrast to "setCoeff" which sets one coefficient).
  void set(int n, NT* coef) { 
    assert(coef != 0);
    _deg = n;
    /*
     * Bug: if n=-1, the following does not do deletion.
     *   Case n=-1 is useful to contract a zero polynomial
     *   whose nominal degree might be zero.  Chee(11/6/2010)
      if (_deg >= 0) {
        delete[] _coeff;
        _coeff = new NT[_deg+1];
        for (int i=0; i<=_deg; ++i) _coeff[i] = coef[i];
      }
     * */
    if (_deg >= 0) {
      delete[] _coeff;
      _coeff = new NT[_deg+1];
      for (int i=0; i<=_deg; ++i) _coeff[i] = coef[i];
    }
  }
  /// assignment function for coeff vector
  void set(int n, const VecNT& coef) {
    assert(coef != 0);
    _deg = n;
    delete[] _coeff;
    if (_deg >= 0) {
      _coeff = new NT[_deg+1];
      for (int i=0; i<=_deg; ++i) _coeff[i] = coef[i];
    }
  }
  //@}
public:
  /// return the coeff (const)
  const NT* coeff() const { return _coeff; }
  /// return the coeff
  NT* coeff() { return _coeff; }
  /// return the degree
  const int& degree() const { return _deg; }
private:
  int _deg;
  NT* _coeff;
};

#ifndef CORE_DISABLE_REFCOUNTING
/// \class RcPolyBase PolyBase.h
/// \brief RcPolyBase is the base class of polynomial class 
///        with reference counting
template <typename NT>
class RcPolyBase : public RcImpl<PolyBase<NT> > {
  typedef RcImpl<PolyBase<NT> > base_cls;
public:
  typedef std::vector<NT> VecNT;
  /// \name constructors, destructor and assignment operator
  //@{
  /// default constructor
  RcPolyBase() : base_cls(new PolyBase<NT>()) {}
  /// copy constructor
  RcPolyBase(const RcPolyBase& rhs) : base_cls(rhs._rep) 
  { this->_rep->inc_rc(); }
  /// destructor
  ~RcPolyBase() { this->_rep->dec_rc(); }

  /// constructor for unit polynomial with nominal degree n
  RcPolyBase(int n) : base_cls(new PolyBase<NT>(n)) {}
  /// constructor with coeff array
  RcPolyBase(int n, NT* coef) : base_cls(new PolyBase<NT>(n, coef)) {}
  /// constructor with coeff vector
  RcPolyBase(int n, const VecNT& coef): base_cls(new PolyBase<NT>(n, coef)) {}
  /// constructor from <tt>char*</tt> (no implicit conversion)
  RcPolyBase(const char* s, char c) : base_cls(new PolyBase<NT>(s, c)) {}
  /// assignment operator for <tt>PolyBase</tt>
  RcPolyBase& operator=(const RcPolyBase& rhs) {
    if (this != &rhs) { 
      this->_rep->dec_rc(); this->_rep = rhs._rep; this->_rep->inc_rc(); 
    } 
    return *this;
  }
  /// assignment function for coeff array
  void set(int n, NT* coef) { base_cls::set(n, coef); }
  /// assignment function for coeff vector
  void set(int n, const VecNT& coef) { base_cls::set(n, coef); }
  //@}
protected:
  /// return the coeff (const)
  const NT* coeff() const { return ((const PolyBase<NT>*)this->_rep)->_coeff; }
  /// return the coeff
  NT* coeff() { this->make_copy(); return this->_rep->_coeff; }
  /// return the degree
  const int& degree() const { return ((const PolyBase<NT>*)this->_rep)->_deg; }
};
#endif

CORE_END_NAMESPACE

#endif /*__CORE_POLYBASE_H__*/
