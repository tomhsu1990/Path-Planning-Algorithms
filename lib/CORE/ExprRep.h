/****************************************************************************
 * ExprRep.h -- Internal Representation for Expr
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
 * $Id: ExprRep.h,v 1.52 2010/11/23 17:58:37 exact Exp $
 ***************************************************************************/
#ifndef __CORE_EXPRREP_H__
#define __CORE_EXPRREP_H__

#include <CORE/CoreAux.h>
#include <CORE/poly/Sturm.h>
#include <CORE/poly/Descartes.h>
#include <bitset>
#include <iostream>
#include <sstream>

CORE_BEGIN_NAMESPACE

enum NODE_NUMTYPE {
  NODE_NT_INTEGER,
  NODE_NT_DYADIC,
  NODE_NT_RATIONAL,
  NODE_NT_ALGEBRAIC,
  NODE_NT_TRANSCENDENTAL,
};

enum NODE_OPTYPE {
  NODE_OP_CONST,
  NODE_OP_ADD,
  NODE_OP_SUB,
  NODE_OP_MUL,
  NODE_OP_DIV,
  NODE_OP_NEG,

  NODE_OP_RATIONAL,
  
  NODE_OP_ANARY,
  NODE_OP_ALGEBRAIC,
  NODE_OP_TRANSCENDENTAL,
};

#define PREC_MIN        MPFR_PREC_MIN
#define DEF_INIT_PREC   53UL

#define    LIST_MODE 0
#define    TREE_MODE 1

#define    SIMPLE_LEVEL 0
#define    DETAIL_LEVEL 1

#define    OPERATOR_ONLY 0
#define    VALUE_ONLY 1
#define    OPERATOR_VALUE 2
#define    FULL_DUMP 3

/// \Pre-defined macros for user own operations
/// \Unary node
#define BEGIN_DEFINE_UNARY_NODE(cls_name)                     \
template <typename RootBd, typename Filter, typename Kernel>  \
class cls_name : public UnaryOpRepT<RootBd, Filter, Kernel> { \
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;           \
  typedef UnaryOpRepT<RootBd, Filter, Kernel> UnaryOpRep;     \
  using UnaryOpRep::child;                                    \
  using ExprRep::filter;                                      \
  using ExprRep::rootBd;                                      \
  using ExprRep::sign;                                        \
  using ExprRep::uMSB;                                        \
  using ExprRep::lMSB;                                        \
  using ExprRep::appValue;                                    \
  using ExprRep::abs2rel;                                     \
  using ExprRep::numType;                                     \
public:                                                       \
  cls_name(ExprRep* c) : UnaryOpRep(c)                        \
  { compute_filter(); compute_numtype(); }                    \
protected:

/// \Binary node
#define BEGIN_DEFINE_BINARY_NODE(cls_name)                    \
template <typename RootBd, typename Filter, typename Kernel>  \
class cls_name : public BinaryOpRepT<RootBd, Filter, Kernel> {\
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;           \
  typedef BinaryOpRepT<RootBd, Filter, Kernel> BinaryOpRep;   \
  using BinaryOpRep::first;                                   \
  using BinaryOpRep::second;                                  \
  using ExprRep::filter;                                      \
  using ExprRep::rootBd;                                      \
  using ExprRep::sign;                                        \
  using ExprRep::uMSB;                                        \
  using ExprRep::lMSB;                                        \
  using ExprRep::appValue;                                    \
  using ExprRep::abs2rel;                                     \
  using ExprRep::numType;                                     \
public:                                                       \
  cls_name(ExprRep* f, ExprRep* s, bool b = false)            \
   : BinaryOpRep(f, s, b)                                     \
  {compute_filter(); compute_numtype(); }                     \
protected:

/// \Anary node
#define BEGIN_DEFINE_KNARY_NODE(cls_name)                     \
template <typename RootBd, typename Filter, typename Kernel>  \
class cls_name : public AnaryOpRepT<RootBd, Filter, Kernel> { \
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;           \
  typedef AnaryOpRepT<RootBd, Filter, Kernel> AnaryOpRep;     \
  using AnaryOpRep::children;                                 \
  using ExprRep::filter;                                      \
  using ExprRep::rootBd;                                      \
  using ExprRep::sign;                                        \
  using ExprRep::uMSB;                                        \
  using ExprRep::lMSB;                                        \
  using ExprRep::appValue;                                    \
  using ExprRep::abs2rel;                                     \
  using ExprRep::numType;                                     \
public:                                                       \
  cls_name(const std::vector<ExprRep*>& c =                   \
           std::vector<ExprRep*>())                           \
   : AnaryOpRep(c)                                            \
  { compute_filter(); compute_numtype(); }                    \
protected:


#define END_DEFINE_UNARY_NODE };

#define BEGIN_DEFINE_RULE(fun_name)                           \
	virtual void fun_name() {

#define BEGIN_DEFINE_RULE_INSERT                              \
public:                                                       \
	void insert(ExprRep* c) {                             

#define BEGIN_DEFINE_RULE_FILTER                              \
	void compute_filter() {
#define BEGIN_DEFINE_RULE_NUMTYPE                             \
	void compute_numtype() {

#define BEGIN_DEFINE_RULE_SIGN                                \
	virtual bool compute_sign() {
#define BEGIN_DEFINE_RULE_UMSB                                \
	virtual bool compute_uMSB() {
#define BEGIN_DEFINE_RULE_LMSB                                \
	virtual bool compute_lMSB() {

#define BEGIN_DEFINE_RULE_R_APPROX                            \
	virtual bool compute_r_approx(prec_t prec) {
#define BEGIN_DEFINE_RULE_A_APPROX                            \
	virtual bool compute_a_approx(prec_t prec) {

#define BEGIN_DEFINE_RULE_ROOTBD                              \
	virtual void compute_rootbd() {

#define END_DEFINE_RULE }

#define DEFINE_UNARY_FUNCTION(fun_name, cls_name)             \
  template<typename T>                                        \
  Expr<T> fun_name(const Expr<T>& e)                          \
  { return new cls_name<T>(e.rep()); }				

#define DEFINE_BINARY_FUNCTION(fun_name, cls_name)            \
  template<typename T>                                        \
  Expr<T> fun_name(const Expr<T>& f, const Expr<T>& s)        \
  { return new cls_name<T>(f.rep(), s.rep()); }				

#define DEFINE_KNARY_FUNCTION(fun_name, cls_name)             \
  template<typename T>                                        \
  Expr<T> fun_name(const std::vector<Expr<T>>& e)             \
  { return new cls_name<T>(e);: }

/// \class ExprRepT
/// \brief represent an expression node in DAG
template <typename RootBd, typename Filter, typename Kernel>
class ExprRepT {
  typedef ExprRepT<RootBd, Filter, Kernel> thisClass;
protected:
  ExprRepT() : m_nodeinfo(0), m_ref_counter(1)
  {} 
  virtual ~ExprRepT()
  { if (m_nodeinfo) delete m_nodeinfo; }

public: // public methods
  /// return the current precision (relative)
  prec_t get_prec() const
  { return appValue().get_prec(); }

  /// return sign
  sign_t get_sign() {
    // if filter works
    // filter is turned off for transcendental nodes
    // it is temporary. eventually we will implement filter for such nodes.
    if (filter().is_ok() && !is_transcendental()) {
      return filter().sign();
    }
    
    // initialize nodeinfo if necessary
    if (!m_nodeinfo)
      init_nodeinfo();
    else { 
      // if cache has sign
      if (flags().test(fSign)) {
        return sign();
      }

      // if kernel has sign
      // Jihun, July 2010
      // Added is_transcendental check
      // This is heck. There's problem of returning a wrong sign
      // from Log node. It is gone when the node is approximated with enough precision
      if (flags().test(fInit) && appValue().has_sign() && !is_transcendental()) {
        return appValue().sgn();
      }
    }
    // do exact evaluation
    if (!compute_sign()) 
      refine();
	
    flags().set(fSign);
    return sign(); 
  }//get_sign()

  /// return uMSB
  msb_t get_uMSB() { 
    // if filter works
    if (filter().is_ok() && !is_transcendental())
      return filter().uMSB();
    // initialize nodeinfo if necessary
    if (!m_nodeinfo)
      init_nodeinfo();
    else {
      // if cache has uMSB
      if (flags().test(fuMSB))
        return uMSB();
      // if kernel is initialized
      if (flags().test(fInit)&&!is_transcendental())
        return appValue().uMSB();
    }
    // do exact evaluation
    if (!compute_uMSB())
      refine();   // Chee: This is an overkill!!  
    
    flags().set(fuMSB);
    return uMSB(); 
  }// get_uMSB


  /// return lMSB
  /** lMSB of an expression x satisfies the relation $|x|\ge 2^{lMSB}$.
   * since MSB(x) is defined to be $\floor{\lg |x|}$
   */
  msb_t get_lMSB() {
    // if filter works
    if (filter().is_ok() && !is_transcendental())
      return filter().lMSB();
    // initialize nodeinfo if necessary
    if (!m_nodeinfo)
      init_nodeinfo();
    else {
      // if cache has lMSB
      if (flags().test(flMSB))
        return lMSB();
      // if kernel is initialized
      if (flags().test(fInit)&&!is_transcendental())
        return appValue().lMSB();
    }
    // do exact evaluation
    if (!compute_lMSB())
      refine();
    
    flags().set(flMSB);
    return lMSB(); 
  }// get_lMSB

  /// return root bound
  RootBd& get_rootBd() { 
   if (!m_nodeinfo) 
     init_nodeinfo();

   if (!flags().test(fRootBd)) {
      compute_rootBd();
      flags().set(fRootBd);
    }
    
    return rootBd(); 
  }
  /// return approximated value w/ relative precision
  /// -- it calls compute_r_approx()
  Kernel& r_approx(prec_t prec) {
    // if filter works and has enough precision
    //if (filter().is_ok() && filter().get_r_prec() > prec)
    //{ return filter().get_r_prec(); }
    // initialize nodeinfo if necessary
    if (!m_nodeinfo) init_nodeinfo();
    // do exact evaluation
    if (is_approx_needed(prec)) {
      if (compute_r_approx(prec)) flags().set(fExact);
      flags().set(fInit);
    }
    return appValue();
  }
  /// return approximated value w/ absolute precision
  Kernel& a_approx(prec_t prec) {
    // initialize nodeinfo if necessary
    if (!m_nodeinfo) init_nodeinfo();
    // if filter works and has enough precision
    //if (filter().is_ok() && filter().get_a_prec() > prec)
    //{ return filter().get_a_prec(); }
    // do exact evaluation
    if (is_approx_needed(abs2rel(prec))) {
      if (compute_a_approx(prec)) flags().set(fExact);
      flags().set(fInit);
    }
    return appValue();
  }
  /// check whether compute_a(r)_approx is needed
  bool is_approx_needed(prec_t prec)
  { return !flags().test(fInit) || (!is_exact() && get_prec() < prec); }

  /// check whether appValue() is exact
  bool is_exact() const
  { return flags().test(fExact); }

  void rootBd_init() {
    rootBd().set_visit(false);
    for (size_t i=0; i < get_children_size(); i++) {
      if (get_child(i)->rootBd().get_visit())
        get_child(i)->rootBd_init();
    }
  }

  void compute_Deg(unsigned long& deg) {
    deg *= rootBd().get_deg();
    rootBd().set_visit(true);
    for (size_t i=0; i < get_children_size(); i++) {
      if (get_child(i)->rootBd().get_visit() == false)
        get_child(i)->compute_Deg(deg);
    }
  }

  // Jihun (Jan, 2010) : We need to compute the degree everytime in order to
  // obtain optimal degree bound.

  unsigned long get_Deg() {
    unsigned long deg = 1; 
    compute_Deg(deg);
    rootBd_init();

    return deg;
  } 

  void dump(int level) {
    if (!m_nodeinfo) {
      std::cout << "[uninitialized node]";
      return;
    }
    if (level == OPERATOR_ONLY) {
#ifdef CORE_DEBUG
      std::cout << op();
#endif
    } else if (level == VALUE_ONLY) {
      std::cout << appValue();
    } else if (level == OPERATOR_VALUE) {
#ifdef CORE_DEBUG
      std::cout << op() << "[val: " << appValue() << "]";
#endif
    } else if (level == FULL_DUMP) {
      std::cout 
#ifdef CORE_DEBUG
      << op()
#endif
      << "[val: "  << appValue() << "; "
      << "lMSB: " << get_lMSB() << "; "
      << "uMSB: " << get_uMSB() << "; "
      << "sign: " << get_sign() << "; "
      << "numType: " << numType() << "; "
      << "rootBd: " << get_rootBd().dump() << "; "
      << "]";
    }
    // note that str() return an array not properly terminated!
  }

  int getDepth() {
    int maxChildDepth(0);

    for(size_t i=0;i<get_children_size();++i)
      maxChildDepth = (std::max)(maxChildDepth, get_child(i)->getDepth());

    return maxChildDepth+1;
  }
  
  void debugList(int level, int depthLimit) {
    if (depthLimit <= 0)
      return;
    if (level == SIMPLE_LEVEL) {
      std::cout << "("; dump(OPERATOR_VALUE);
    } else if (level == DETAIL_LEVEL) {
      std::cout << "("; dump(FULL_DUMP);
    }  
    for (size_t i = 0; i < get_children_size(); i++)
      get_child(i)->debugList(level, depthLimit - 1);
    std::cout << ")";
  }

  void debugTree(int level, int indent, int depthLimit) {
    if (depthLimit <= 0)
      return;
    for (int i = 0; i<indent; i++ )
      std::cout << "  ";
    std::cout << "|_";
    if (level == SIMPLE_LEVEL)
      dump(OPERATOR_VALUE);
    else if (level == DETAIL_LEVEL)
      dump(FULL_DUMP);
    std::cout << std::endl;
    for (size_t i = 0; i < get_children_size(); i++)
      get_child(i)->debugTree(level, indent + 2, depthLimit - 1);
  }

  prec_t get_constructive_bound() {
    RootBd& rootBd = get_rootBd();
    return rootBd.get_bound(get_Deg());
    // Jan, 2010 : We have to unroll the following to avoid the m_nodeinfo
    // initialization error. Probably a compiler bug.
    //return get_rootBd().get_bound(get_Deg());
  } 				      
  ////////////////////////////////////////////////// 
  /// refine() is called by compute_sign and compute_lMSB
  /// 	-- it computes the root bound and then improves the
  ///   current approximation until the root bound.
  ///   -- bound_type=0 if cut_off_bound
  /// 		     =1 if escape_bound
  ///                =2 if constructive_root_bound
  //
  // Aug 2006: Jihun pulled this refine() function out of AddSubRep.
  // Oct 2006: Jihun/Chee implemented new refine() from Zilin's thesis (p.41)
  //////////////////////////////////////////////////

  void refine() {
#ifdef CORE_DEBUG_ROOTBOUND
      std::cerr << std::endl << "In refine," << std::endl;
#endif
  
    //Step 1
    unsigned long prec = 53, bound = 53;
    int bound_type;

    if (!kernel_initialized())
      initialize_kernel();
    else
      prec = get_prec()*2;
    //Step 2
    if (is_transcendental()) {
      bound = get_escape_bound();
      bound_type = 1;
    } else {
      if (get_rootBd().is_constructive())
			bound = get_constructive_bound();
      bound_type = 2;
    }
    // Step 3
    if (bound > get_cut_off_bound()) {
#ifdef CORE_DEBUG_ROOTBOUND
      std::cerr << "root bound is bigger than cut off bound." << std::endl;
      std::cerr << "root bound=" << bound << std::endl;
      std::cerr << "cut off bound=" << get_cut_off_bound() << std::endl;
#endif
      bound = get_cut_off_bound(); bound_type = 0;
    }
#ifdef CORE_DEBUG_ROOTBOUND
      std::cerr << "refine bound type=" << bound_type << std::endl;
      std::cerr << "refine bound=" << bound << std::endl;
#endif
    // Step 4
	bound = (std::max)(bound,53UL);
    do {
	  unsigned long realPrec = (std::min)(prec,bound);
//	  unsigned long realPrec = prec;
#ifdef CORE_DEBUG_ROOTBOUND
      std::cerr << "refine prec=" << realPrec << std::endl;
#endif
      a_approx(realPrec);
      if (appValue().has_sign()) {
        set_flags(appValue().sgn(), appValue().uMSB(), appValue().lMSB());
#ifdef CORE_DEBUG_ROOTBOUND
        std::cerr << "found sign =" << appValue().sgn() << std::endl;
#endif
        return;
      }
      prec <<= 1;
    } while (prec < 2*bound);
#ifdef CORE_DEBUG_ROOTBOUND
    std::cerr << "root bound=" << get_constructive_bound() << std::endl;
    std::cerr << "root type=" << bound_type << std::endl;
    rootBd().dump();
#endif

    // Step 5
    if (bound_type < 2) {
      if (bound_type==0) {
        core_error("Hit Cut-off Bound", __FILE__, __LINE__, false);
        core_error("ZERO ASSERTION: correctnes conditioned on value = 0", __FILE__, __LINE__, false);
      } else { // bound_type=1
        core_error("Hit Escape Bound", __FILE__, __LINE__, false);
        core_error("ZERO ASSERTION: correctnes conditioned on value = 0", __FILE__, __LINE__, false);
      }
    }
    appValue().set(0);
    set_flags(0, 0, MSB_MIN);
  }
/*
  void refine() {
    prec_t rootbd;
    if(get_rootBd().is_constructive())
      rootbd = get_rootBd().get_bound(get_Deg());
#ifdef CORE_DEBUG_ROOTBOUND
    std::cout << "\nrootbd degree upperbound = " << get_Deg() << std::endl;
    std::cout << "rootbd bit = " << rootbd << std::endl;
#endif
    
    for (prec_t prec=(std::min)(rootbd, DEF_INIT_PREC); prec<2*rootbd; prec<<=1) {
#ifdef CORE_DEBUG_ROOTBOUND
      std::cerr << "rootbd prec=" << prec << std::endl;
#endif
      a_approx(prec);
      if (!appValue().has_zero()) {
        set_flags(appValue().sgn(), appValue().uMSB(), appValue().lMSB());
#ifdef CORE_DEBUG_ROOTBOUND
        std::cerr << "found sign =" << appValue().sgn() << std::endl;
#endif
        return;
      }
    }
#ifdef CORE_DEBUG_ROOTBOUND
    std::cerr << "root bound=" << rootBd().get_bound(get_Deg()) << std::endl;
    rootBd().dump();
#endif
    set_flags(0, 0, MSB_MIN);
  }//refine
*/

public: 
  /// convert absolute precision to relative precision
  /// WARNING: we currently do not check for overflows
  // TODO: have a "strict mode" where we do check overflow!
  // Note this output is at least PREC_MIN (=2).
  prec_t abs2rel(prec_t prec) {
    //core_assert(prec >= MPFR_PREC_MIN && prec <= MPFR_PREC_MAX)
    prec_t ret= (std::max)(long(prec) + get_uMSB(), long(PREC_MIN)); 
    if (ret < MPFR_PREC_MIN || ret > MPFR_PREC_MAX)
      std::cout << "abs2rel error" << std::endl;
    return (std::max)(long(prec) + get_uMSB(), long(PREC_MIN)); 
  }
  /// convert relative precision to absolute precision
  /// WARNING: we currently do not check for overflows
  prec_t rel2abs(prec_t prec) {
    //core_assert(prec >= MPFR_PREC_MIN && prec <= MPFR_PREC_MAX)
    prec_t ret = (std::max)(long(prec) - get_lMSB(), long(PREC_MIN)); 
    if (ret < MPFR_PREC_MIN || ret > MPFR_PREC_MAX)
      std::cout << "rel2abs error" << std::endl;
    return (std::max)(long(prec) - get_lMSB(), long(PREC_MIN)); 
  }

protected: 
  /// set flags
  void set_flags(const sign_t& s, const msb_t& u, const msb_t& l) {
    assert(m_nodeinfo);
    sign() = s; flags().set(fSign);
    uMSB() = u; flags().set(fuMSB);
    lMSB() = l; flags().set(flMSB);
  }

  template <typename V>
  void init_value(const V& value) {
    if (!m_nodeinfo)
      new_nodeinfo();
    // initialize approximated value
    if (this->appValue().set(value)) this->flags().set(fExact);
    // set flag
    flags().set(fInit);
  }
protected: // overridable methods
  /// initialize nodeinfo
  virtual void init_nodeinfo()
  { new_nodeinfo(); }
  /// compute sign
  virtual bool compute_sign()
  { return false;}
  /// compute uMSB
  virtual bool compute_uMSB()
  { return false;}
  /// compute lMSB
  virtual bool compute_lMSB()
  { return false;}
  /// compute rootBd
  virtual void compute_rootBd()
  {}
#ifdef CORE_DEBUG
  virtual std::string op()
  { return std::string("undefined"); }
#endif 

  /// compute relative approximation (default)
  /// return true if the approximation can be represented exactly by MPFR) 
  virtual bool compute_r_approx(prec_t prec)
  { return compute_a_approx(rel2abs(prec)); } 

  /// compute absolute approximation (default)
  /// return true if the approximation can be represented exactly by MPFR) 
  virtual bool compute_a_approx(prec_t prec)
  { return compute_r_approx(abs2rel(prec)); }
public: 
  virtual size_t get_children_size() const
  { return 0; }
  virtual thisClass* get_child(size_t t) const
  { assert(0); return 0; }

protected: 
  bool is_transcendental() const 
  { return numType() == NODE_NT_TRANSCENDENTAL; }
  
  enum CachedFlags {
    fSign, fuMSB, flMSB, fRootBd, 
    fInit, fExact, numFlags
  };
  typedef std::bitset<numFlags> flag_t;

  const bool kernel_initialized() const {
    return m_nodeinfo != 0;
  }
  
  void initialize_kernel() {
    new_nodeinfo();
  }

  const flag_t& flags() const { 
    assert(m_nodeinfo);
    return m_nodeinfo->m_flags;
  }
  flag_t& flags()  { 
    assert(m_nodeinfo);
    return m_nodeinfo->m_flags;
  }

  const sign_t& sign() const { 
    assert(m_nodeinfo);
   return m_nodeinfo->m_sign;
  }
  sign_t& sign()  { 
    assert(m_nodeinfo);
    return m_nodeinfo->m_sign;
  }

  const msb_t& uMSB() const  { 
    assert(m_nodeinfo);
    return m_nodeinfo->m_uMSB;
  }
  msb_t& uMSB()  { 
    assert(m_nodeinfo);
    return m_nodeinfo->m_uMSB;
  }

  const msb_t& lMSB() const {
    assert(m_nodeinfo);
    return m_nodeinfo->m_lMSB;
  }
  msb_t& lMSB()  { 
    assert(m_nodeinfo);
    return m_nodeinfo->m_lMSB;
  }

  const RootBd& rootBd() const  { 
    assert(m_nodeinfo);
    return m_nodeinfo->m_rootBd;
  }
  RootBd& rootBd()  { 
    assert(m_nodeinfo);
    return m_nodeinfo->m_rootBd;
  }

  void new_nodeinfo() { m_nodeinfo = new NodeInfo(); }
public:
  Kernel& appValue()  { 
    assert(m_nodeinfo);
    return m_nodeinfo->m_appValue;
  }
  const Kernel& appValue() const  { 
    assert(m_nodeinfo);
    return m_nodeinfo->m_appValue;
  }
	
  Filter& filter() { return m_filter; }
  const Filter& filter() const { return m_filter; }

  NODE_NUMTYPE& numType() { return _numType; }
  const NODE_NUMTYPE& numType() const { return _numType; }

  NODE_OPTYPE& opType() { return _opType; }
  const NODE_OPTYPE& opType() const { return _opType; }

  virtual BigInt getZTVal() { return 0; }
  virtual BigFloat getFTVal() { return 0; }
  virtual BigRat getQTVal() { return 0; }
public:
  /// node information
  struct NodeInfo {
    flag_t  m_flags;    ///<- flags
    sign_t  m_sign;     ///<- sign
    msb_t   m_uMSB;     ///<- upper bound of Most Significant Bit
    msb_t   m_lMSB;     ///<- low bound of Most Significant Bit
    RootBd  m_rootBd;   ///<- root bound
    Kernel  m_appValue; ///<- apprixmated value
  };

  NodeInfo* m_nodeinfo; ///<- node information
  Filter    m_filter;   ///<- filter
  NODE_NUMTYPE _numType; ///<- number type
  NODE_OPTYPE  _opType; ///<- operator type

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

/// These 12 functions were introduced to support getZTVal(), getFTVal(), getQTVal()
///  which are in turn used in getExactVal(), for use in ReduceToRational.
inline BigInt getBigIntVal(const BigRat& v) { return 0; }
inline BigInt getBigIntVal(const BigFloat& v) { return 0; }
inline BigInt getBigIntVal(const BigInt& v) { return v; }
inline BigInt getBigIntVal(const long& v) { return (BigInt)v; }

inline BigFloat getBigFloatVal(const BigRat& v) { return 0; }
inline BigFloat getBigFloatVal(const BigFloat& v) { return v; }
inline BigFloat getBigFloatVal(const BigInt& v) { return (BigFloat)v; }
inline BigFloat getBigFloatVal(const long& v) { return (BigFloat)v; }

inline BigRat getBigRatVal(const BigRat& v) { return v; }
inline BigRat getBigRatVal(const BigFloat& v) { return (BigRat)v; }
inline BigRat getBigRatVal(const BigInt& v) { return (BigRat)v; }
inline BigRat getBigRatVal(const long& v) { return (BigRat)v; }

/// \class ConstRepT 
/// \brief constant node
template <typename RootBd, typename Filter, typename Kernel, typename T>
class ConstRepT : public ExprRepT<RootBd, Filter, Kernel> {
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;
  typedef RootBd* id_rootbd_t;
  using ExprRep::filter;
  using ExprRep::sign;
  using ExprRep::uMSB;
  using ExprRep::lMSB;
  using ExprRep::appValue;
  using ExprRep::rootBd;
  using ExprRep::set_flags;
  using ExprRep::new_nodeinfo;
  using ExprRep::init_value;
  using ExprRep::init_nodeinfo;
  using ExprRep::numType;
  using ExprRep::opType;
protected:
  T value;
public:
  template <typename V>
  ConstRepT(const V& v, NODE_NUMTYPE t) : value(v) {
    filter().set(v);
    numType() = t;
    opType() = NODE_OP_CONST;
    init_nodeinfo();
  }  
  virtual ~ConstRepT() 
  {}

protected:
  virtual void init_nodeinfo() {
    init_value(value);
  }

  virtual bool compute_sign()
  { sign() = sgn(value); return true; }
  virtual bool compute_uMSB()
  { uMSB() = ceillg(value); return true; }
  virtual bool compute_lMSB()
  { lMSB() = floorlg(value); return true; }
  virtual void compute_rootBd()
  { rootBd().set(value); }
  virtual bool compute_r_approx(prec_t prec)
  { return appValue().set(value, prec); }
public:
  T getValue() const { return value; }
  /// In the following getBigXXXValue(), if the return type is smaller than
  /// type T, then it should be an error, but we return 0.  
  BigInt getZTVal() {return getBigIntVal(value); }
  BigFloat getFTVal() {return getBigFloatVal(value); }
  BigRat getQTVal() { return getBigRatVal(value); }

#ifdef CORE_DEBUG
  virtual std::string op()
  { return std::string("Const"); }
#endif 

private:
  // generic version for BigInt, BigRat, Mpfr, BigFloat
  template <typename V> void init_flags(const V& v) 
  { set_flags(v.sgn(), v.uMSB(), v.lMSB()); }
  // specialized version for long
  void init_flags(long v) 
  { set_flags(sgn(v), ceillg(v), floorlg(v)); }
  // specialized version for unsigned long
  void init_flags(unsigned long v) 
  { set_flags(sgn(v), ceillg(v), floorlg(v)); }
  // specialized version for double
  void init_flags(double v) 
  { set_flags(sgn(v), ceillg(v), floorlg(v)); }
};

/// \class ConstPolyRepT 
/// \brief Const Poly node
template <typename RootBd, typename Filter, typename Kernel,typename NT>
class ConstPolyRepT : public ExprRepT<RootBd, Filter, Kernel> {
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;
  using ExprRep::filter;
  using ExprRep::sign;
  using ExprRep::uMSB;
  using ExprRep::lMSB;
  using ExprRep::appValue;
  using ExprRep::rootBd;
  using ExprRep::numType;
  using ExprRep::opType;
protected:
  Kernel value;
  Descartes<NT> ss; ///< internal Sturm sequences
  BFInterval I; ///< current interval contains the real value
public:
  ConstPolyRepT(const Polynomial<NT>& p, int n) : ss(p) {
    I = ss.isolateRoot(n);
    // check whether n-th root exists
    if (I.first == 1 && I.second == 0) {
      core_error("root index out of bound", __FILE__, __LINE__, true);
      assert(0);
    }
    // refine initial interval to absolte error of 2^53
    I = ss.newtonRefine(I, 54);
    // we get an exact root
    if (I.first == I.second) value = I.first;
    else value.set(BigFloat2(I.first, I.second), 53);
    filter().set(value);
    numType() = NODE_NT_ALGEBRAIC;
    opType() = NODE_OP_CONST;
  }

  ConstPolyRepT(const Polynomial<NT>& p, const BFInterval& II) : ss(p), I(II) {
    BFVecInterval v;
    ss.isolateRoots(I.first, I.second, v);
    if (v.size() != 1) {
      core_error("root interval is not isolating", __FILE__, __LINE__, true);
      assert(0);
    }
    I = v.front();
    // refine initial interval to absolute error of 2^53
    I = ss.newtonRefine(I, 54);
    // we get an exact root
    if (I.first == I.second) value.set(I.first, 53);
    else value.set(BigFloat2(I.first, I.second), 53);
    filter().set(value);
    numType() = NODE_NT_ALGEBRAIC;
  }

  virtual ~ConstPolyRepT() 
  {}

protected:
  virtual bool compute_sign()
  { sign() = value.sgn(); return true; }
  virtual bool compute_uMSB()
  { uMSB() = value.uMSB(); return true; }
  virtual bool compute_lMSB()
  { lMSB() = value.lMSB(); return true; }
  virtual void compute_rootBd()
  { rootBd().set(value); }
  virtual bool compute_a_approx(prec_t prec) {
    I = ss.newtonRefine(I, prec+1);
    assert(I.first <= I.second && I.first*I.second >= 0);
    if (I.first == I.second) value.set(I.first, prec);
    else value.set(BigFloat2(I.first, I.second), prec);
    return appValue().set(value);
  }
#ifdef CORE_DEBUG
  virtual std::string op()
  { return std::string("rootOf"); }
#endif 
};

/// \class PiRepT
template <typename RootBd, typename Filter, typename Kernel>
class PiRepT : public ExprRepT<RootBd, Filter, Kernel> {
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;
  using ExprRep::filter;
  using ExprRep::sign;
  using ExprRep::uMSB;
  using ExprRep::lMSB;
  using ExprRep::appValue;
  using ExprRep::rootBd;
  using ExprRep::numType;
  using ExprRep::opType;
public:
  PiRepT() {
    compute_filter();
    compute_numtype();
    opType() = NODE_OP_CONST;
  }
  virtual ~PiRepT() 
  {}
protected:
  void compute_filter() {
    filter().set(3.1415926535897932384626433832795028F);// Chee: this is strange,
    // I expect filter values to have error bounds...
  }
  void compute_numtype() {
    numType() = NODE_NT_TRANSCENDENTAL;
  }
  virtual bool compute_sign() 
  { sign() = 1; return true;}
  virtual bool compute_uMSB() 
  { uMSB() = 2; return true;}
  virtual bool compute_lMSB() 
  { lMSB() = 1; return true;}
  virtual void compute_rootBd()
  { rootBd().set(3.14); }
  virtual bool compute_r_approx(prec_t prec)
  { return appValue().pi(prec);}
#ifdef CORE_DEBUG
  virtual std::string op()
  { return std::string("Pi"); }
#endif  
};

/// \class ERepT
template <typename RootBd, typename Filter, typename Kernel>
class ERepT : public ExprRepT<RootBd, Filter, Kernel> {
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;
  using ExprRep::filter;
  using ExprRep::sign;
  using ExprRep::uMSB;
  using ExprRep::lMSB;
  using ExprRep::appValue;
  using ExprRep::rootBd;
  using ExprRep::numType;
  using ExprRep::opType;
public:
  ERepT() {
    compute_filter();
    compute_numtype();
    opType() = NODE_OP_CONST;
  }
  virtual ~ERepT() 
  {}
protected:
  void compute_filter() {
    filter().set(2.72);// Chee: this is strange -- filters should be an interval
    // need error bound?
  }
  void compute_numtype() {
    numType() = NODE_NT_TRANSCENDENTAL;
  }
  virtual bool compute_sign() 
  { sign() = 1; return true;}
  virtual bool compute_uMSB() 
  { uMSB() = 2; return true;}
  virtual bool compute_lMSB() 
  { lMSB() = 1; return true;}
  virtual void compute_rootBd()
  { rootBd().set(2.72); }
  virtual bool compute_r_approx(prec_t prec)
  { return appValue().e(prec+1);}
};

/// \class UnaryOpRepT
/// \brief unary operation expression node
template <typename RootBd, typename Filter, typename Kernel>
class UnaryOpRepT : public ExprRepT<RootBd, Filter, Kernel> {
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;
public:
  UnaryOpRepT(ExprRep* c) : child(c)
  { child->inc_ref(); }
  virtual ~UnaryOpRepT()
  { child->dec_ref(); }
  virtual size_t get_children_size() const
  { return 1; }
  virtual ExprRep* get_child(size_t t) const
  { return child; }

protected:
  /// check whether the operation and child node are both exact
  bool check_exact(bool ret)
  { return ret && child->is_exact(); }
  ExprRepT<RootBd, Filter, Kernel>* child; /// <- pointer to the child node
};

/// \class NegRepT
/// \brief negation node
template <typename RootBd, typename Filter, typename Kernel>
class NegRepT : public UnaryOpRepT<RootBd, Filter, Kernel> {
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;
  typedef UnaryOpRepT<RootBd, Filter, Kernel> UnaryOpRep;
  typedef RootBd* id_rootbd_t;
  using UnaryOpRep::child;
  using UnaryOpRep::check_exact;
  using ExprRep::filter;
  using ExprRep::sign;
  using ExprRep::uMSB;
  using ExprRep::lMSB;
  using ExprRep::appValue;
  using ExprRep::rootBd;
  using ExprRep::abs2rel;
  using ExprRep::numType;
  using ExprRep::opType;
public:
  NegRepT(ExprRep* c) : UnaryOpRep(c)  {
    filter().neg(child->filter());
    numType() = child->numType();
    opType() = NODE_OP_NEG;
  }  
  virtual ~NegRepT() 
  {}
protected:
  /* notice here we use lazy evaluation */
  virtual bool compute_sign() 
  { sign() = -child->get_sign(); return true;}
  virtual bool compute_uMSB() 
  { uMSB() = child->get_uMSB(); return true;}
  virtual bool compute_lMSB() 
  { lMSB() = child->get_lMSB(); return true;}
  virtual void compute_rootBd()
  { rootBd().neg(child->get_rootBd()); }
  virtual bool compute_r_approx(prec_t prec) {
    return check_exact(appValue().neg(child->r_approx(prec), prec));
  } 
  virtual bool compute_a_approx(prec_t prec) {
    return this->check_exact(appValue().neg(child->a_approx(prec),
			    abs2rel(prec)));
  }
  BigInt getZTVal() { return -child->getZTVal(); }
  BigFloat getFTVal() { return -child->getFTVal(); }
  BigRat getQTVal() { return -child->getQTVal(); }
#ifdef CORE_DEBUG
  virtual std::string op()
  { return std::string("Neg"); }
#endif 
};

/// \class SqrtRepT
/// \brief square root node
template <typename RootBd, typename Filter, typename Kernel>
class SqrtRepT : public UnaryOpRepT<RootBd, Filter, Kernel> {
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;
  typedef UnaryOpRepT<RootBd, Filter, Kernel> UnaryOpRep;
  typedef RootBd* id_rootbd_t;
  using UnaryOpRep::child;
  using UnaryOpRep::check_exact;
  using ExprRep::filter;
  using ExprRep::sign;
  using ExprRep::uMSB;
  using ExprRep::lMSB;
  using ExprRep::appValue;
  using ExprRep::rootBd;
  using ExprRep::abs2rel;
  using ExprRep::init_value;  
  using ExprRep::numType;  
  using ExprRep::opType;
public:
  SqrtRepT(ExprRep* c) : UnaryOpRep(c) 
  { filter().sqrt(child->filter()); 
    numType() = (std::max)(NODE_NT_ALGEBRAIC, child->numType());
    opType() = NODE_OP_ALGEBRAIC;
    if (child->get_sign()<0) 
      core_error("sqrt of negative value", __FILE__, __LINE__, true);
    if (child->get_sign()==0)
      init_value(0);
  }
  virtual ~SqrtRepT() 
  {}
protected:
  virtual bool compute_sign() 
  { sign() = child->get_sign(); return true;}
  virtual bool compute_uMSB() 
  { uMSB() = (child->get_uMSB()+1) >> 1; return true;}
  virtual bool compute_lMSB() 
  { lMSB() = child->get_lMSB() >> 1; return true;}
  virtual void compute_rootBd()
  { rootBd().root(child->get_rootBd(), 2); }
  virtual bool compute_r_approx(prec_t prec) {
    // compute_r_approx is called by r_approx only --
    // if child's sign is zero, compute_r_approx is not called	  
    return check_exact(appValue().sqrt(child->r_approx(prec*2), prec));
  }
#ifdef CORE_DEBUG
  virtual std::string op()
  { return std::string("Sqrt"); }
#endif 
};

/// \class CbrtRepT
/// \brief cubic root node
template <typename RootBd, typename Filter, typename Kernel>
class CbrtRepT : public UnaryOpRepT<RootBd, Filter, Kernel> {
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;
  typedef UnaryOpRepT<RootBd, Filter, Kernel> UnaryOpRep;
  typedef RootBd* id_rootbd_t;
  using UnaryOpRep::child;
  using UnaryOpRep::check_exact;
  using ExprRep::filter;
  using ExprRep::sign;
  using ExprRep::uMSB;
  using ExprRep::lMSB;
  using ExprRep::appValue;
  using ExprRep::rootBd;
  using ExprRep::abs2rel;
  using ExprRep::init_value;
  using ExprRep::numType;
  using ExprRep::opType;
public:
  CbrtRepT(ExprRep* c) : UnaryOpRep(c) 
  { filter().cbrt(child->filter());
    numType() = (std::max)(NODE_NT_ALGEBRAIC, child->numType());
    opType() = NODE_OP_ALGEBRAIC;
    if (child->get_sign()==0)
      init_value(0);
  }
  virtual ~CbrtRepT() 
  {}
protected:
  virtual bool compute_sign() 
  { sign() = child->get_sign(); return true;}
  virtual bool compute_uMSB() 
  { uMSB() = (child->get_uMSB()+2)/3; return true;}
  virtual bool compute_lMSB() 
  { lMSB() = child->get_lMSB()/3; return true;}
  virtual void compute_rootBd()
  { rootBd().root(child->get_rootBd(), 3); }
  virtual bool compute_r_approx(prec_t prec) {
    return check_exact(appValue().cbrt(child->r_approx(prec*3), prec));
  }
#ifdef CORE_DEBUG
  virtual std::string op()
  { return std::string("Cbrt"); }
#endif 
};

/// \class Sin
/// \brief geometric sine node
template <typename RootBd, typename Filter, typename Kernel>
class SinRepT : public UnaryOpRepT<RootBd, Filter, Kernel> {
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;
  typedef UnaryOpRepT<RootBd, Filter, Kernel> UnaryOpRep;
  typedef RootBd* id_rootbd_t;
  using UnaryOpRep::child;
  using UnaryOpRep::check_exact;
  using ExprRep::filter;
  using ExprRep::sign;
  using ExprRep::uMSB;
  using ExprRep::lMSB;
  using ExprRep::appValue;
  using ExprRep::rootBd;
  using ExprRep::abs2rel;
  using ExprRep::init_value;
  using ExprRep::numType;
  using ExprRep::opType;
public:
  SinRepT(ExprRep* c) : UnaryOpRep(c) {
    numType() = (std::max)(NODE_NT_TRANSCENDENTAL, child->numType());
    opType() = NODE_OP_TRANSCENDENTAL;
    if (child->get_sign()==0)
      init_value(0);
  }
  virtual ~SinRepT() 
  {}
protected:
  virtual bool compute_sign() 
  { return false;}
  virtual bool compute_uMSB() 
  { uMSB() = 0; return true;}
  virtual bool compute_lMSB() 
  { return false;}
  virtual bool compute_a_approx(prec_t prec) {
    return check_exact(appValue().sin(child->a_approx(prec+1), abs2rel(prec+1)));
  }
#ifdef CORE_DEBUG
  virtual std::string op()
  { return std::string("Sin"); }
#endif 
};
/// \class Cos
/// \brief geometric cosine node
template <typename RootBd, typename Filter, typename Kernel>
class CosRepT : public UnaryOpRepT<RootBd, Filter, Kernel> {
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;
  typedef UnaryOpRepT<RootBd, Filter, Kernel> UnaryOpRep;
  typedef RootBd* id_rootbd_t;
  using UnaryOpRep::child;
  using UnaryOpRep::check_exact;
  using ExprRep::filter;
  using ExprRep::sign;
  using ExprRep::uMSB;
  using ExprRep::lMSB;
  using ExprRep::appValue;
  using ExprRep::rootBd;
  using ExprRep::abs2rel;
  using ExprRep::init_value;
  using ExprRep::numType;
  using ExprRep::opType;
public:
  CosRepT(ExprRep* c) : UnaryOpRep(c) {
    numType() = (std::max)(NODE_NT_TRANSCENDENTAL, child->numType());
    opType() = NODE_OP_TRANSCENDENTAL;
    if (child->get_sign()==0)
      init_value(1);
  }
  virtual ~CosRepT() 
  {}
protected:
  virtual bool compute_sign() 
  { return false;}
  virtual bool compute_uMSB() 
  { uMSB() = 0; return true;}
  virtual bool compute_lMSB() 
  { return false;}
  virtual bool compute_a_approx(prec_t prec) {
    return check_exact(appValue().cos(child->a_approx(prec+1), abs2rel(prec+1)));
  }
#ifdef CORE_DEBUG
  virtual std::string op()
  { return std::string("Cos"); }
#endif 
};
/// \class Tan
/// \brief geometric tangent node
template <typename RootBd, typename Filter, typename Kernel>
class TanRepT : public UnaryOpRepT<RootBd, Filter, Kernel> {
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;
  typedef UnaryOpRepT<RootBd, Filter, Kernel> UnaryOpRep;
  typedef RootBd* id_rootbd_t;
  using UnaryOpRep::child;
  using UnaryOpRep::check_exact;
  using ExprRep::filter;
  using ExprRep::sign;
  using ExprRep::uMSB;
  using ExprRep::lMSB;
  using ExprRep::appValue;
  using ExprRep::rootBd;
  using ExprRep::abs2rel;
  using ExprRep::init_value;
  using ExprRep::numType;
  using ExprRep::opType;
public:
  TanRepT(ExprRep* c) : UnaryOpRep(c) {
    numType() = (std::max)(NODE_NT_TRANSCENDENTAL, child->numType());
    opType() = NODE_OP_TRANSCENDENTAL;
    if (child->get_sign()==0)
      init_value(0);
  }
  virtual ~TanRepT() 
  {}
protected:
  virtual bool compute_sign() 
  { return false;}
  virtual bool compute_uMSB() 
  { return false;}
  virtual bool compute_lMSB() 
  { return false;}
  virtual bool compute_a_approx(prec_t prec) {
    return check_exact(appValue().tan(child->a_approx(prec+3), abs2rel(prec+1)));
  }
#ifdef CORE_DEBUG
  virtual std::string op()
  { return std::string("Tan"); }
#endif 
};
/// \class Cot
/// \brief geometric cotangent node
template <typename RootBd, typename Filter, typename Kernel>
class CotRepT : public UnaryOpRepT<RootBd, Filter, Kernel> {
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;
  typedef UnaryOpRepT<RootBd, Filter, Kernel> UnaryOpRep;
  typedef RootBd* id_rootbd_t;
  using UnaryOpRep::child;
  using UnaryOpRep::check_exact;
  using ExprRep::filter;
  using ExprRep::sign;
  using ExprRep::uMSB;
  using ExprRep::lMSB;
  using ExprRep::appValue;
  using ExprRep::rootBd;
  using ExprRep::abs2rel;
  using ExprRep::init_value;
  using ExprRep::numType;
  using ExprRep::opType;
public:
  CotRepT(ExprRep* c) : UnaryOpRep(c) {
    numType() = (std::max)(NODE_NT_TRANSCENDENTAL, child->numType());
    opType() = NODE_OP_TRANSCENDENTAL;
  }
  virtual ~CotRepT() 
  {}
protected:
  virtual bool compute_sign() 
  { return false;}
  virtual bool compute_uMSB() 
  { return false;}
  virtual bool compute_lMSB() 
  { return false;}
  virtual bool compute_a_approx(prec_t prec) {
    return check_exact(appValue().cot(child->a_approx(prec+2), abs2rel(prec+1)));
  }
#ifdef CORE_DEBUG
  virtual std::string op()
  { return std::string("Cot"); }
#endif 
};
/// \class ASin
/// \brief geometric arcsine node
template <typename RootBd, typename Filter, typename Kernel>
class ASinRepT : public UnaryOpRepT<RootBd, Filter, Kernel> {
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;
  typedef UnaryOpRepT<RootBd, Filter, Kernel> UnaryOpRep;
  typedef RootBd* id_rootbd_t;
  using UnaryOpRep::child;
  using UnaryOpRep::check_exact;
  using ExprRep::filter;
  using ExprRep::sign;
  using ExprRep::uMSB;
  using ExprRep::lMSB;
  using ExprRep::appValue;
  using ExprRep::rootBd;
  using ExprRep::abs2rel;
  using ExprRep::init_value;
  using ExprRep::numType;
  using ExprRep::opType;
public:
  ASinRepT(ExprRep* c) : UnaryOpRep(c) {
    numType() = (std::max)(NODE_NT_TRANSCENDENTAL, child->numType());
    opType() = NODE_OP_TRANSCENDENTAL;
  }
  virtual ~ASinRepT() 
  {}
protected:
  virtual bool compute_sign() 
  { sign() = child->get_sign(); return true; }
  virtual bool compute_uMSB() 
  { uMSB() = 1; return true;}
  virtual bool compute_lMSB() 
  { lMSB() = floorlg(abs(child->a_approx(2).get_min())); return true;}
  virtual bool compute_a_approx(prec_t prec) {
    return check_exact(appValue().asin(child->a_approx(prec+2), abs2rel(prec+1)));
  }
#ifdef CORE_DEBUG
  virtual std::string op()
  { return std::string("ASin"); }
#endif 
};
/// \class ACos
/// \brief geometric arccosine node
template <typename RootBd, typename Filter, typename Kernel>
class ACosRepT : public UnaryOpRepT<RootBd, Filter, Kernel> {
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;
  typedef UnaryOpRepT<RootBd, Filter, Kernel> UnaryOpRep;
  typedef RootBd* id_rootbd_t;
  using UnaryOpRep::child;
  using UnaryOpRep::check_exact;
  using ExprRep::filter;
  using ExprRep::sign;
  using ExprRep::uMSB;
  using ExprRep::lMSB;
  using ExprRep::appValue;
  using ExprRep::rootBd;
  using ExprRep::abs2rel;
  using ExprRep::init_value;
  using ExprRep::numType;
  using ExprRep::opType;
public:
  ACosRepT(ExprRep* c) : UnaryOpRep(c) {
    numType() = (std::max)(NODE_NT_TRANSCENDENTAL, child->numType());
    opType() = NODE_OP_TRANSCENDENTAL;
  }
  virtual ~ACosRepT() 
  {}
protected:
  virtual bool compute_sign() 
  { sign() = 1; return true;}
  virtual bool compute_uMSB() 
  { uMSB() = 2; return true;}
  virtual bool compute_lMSB() 
  { return false;}
  virtual bool compute_a_approx(prec_t prec) {
    return check_exact(appValue().acos(child->a_approx(prec+2), abs2rel(prec+1)));
  }
#ifdef CORE_DEBUG
  virtual std::string op()
  { return std::string("ACos"); }
#endif 
};
/// \class ATan
/// \brief geometric arctangent node
template <typename RootBd, typename Filter, typename Kernel>
class ATanRepT : public UnaryOpRepT<RootBd, Filter, Kernel> {
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;
  typedef UnaryOpRepT<RootBd, Filter, Kernel> UnaryOpRep;
  typedef RootBd* id_rootbd_t;
  using UnaryOpRep::child;
  using UnaryOpRep::check_exact;
  using ExprRep::filter;
  using ExprRep::sign;
  using ExprRep::uMSB;
  using ExprRep::lMSB;
  using ExprRep::appValue;
  using ExprRep::rootBd;
  using ExprRep::abs2rel;
  using ExprRep::init_value;
  using ExprRep::numType;
  using ExprRep::opType;
public:
  ATanRepT(ExprRep* c) : UnaryOpRep(c) {
    numType() = (std::max)(NODE_NT_TRANSCENDENTAL, child->numType());
    opType() = NODE_OP_TRANSCENDENTAL;
  }
  virtual ~ATanRepT() 
  {}
protected:
  virtual bool compute_sign() 
  { return false;}
  virtual bool compute_uMSB() 
  { uMSB() = 1; return true;}
  virtual bool compute_lMSB() 
  { return false;}
  virtual bool compute_a_approx(prec_t prec) {
    return check_exact(appValue().atan(child->a_approx(prec+2), abs2rel(prec+2)));
  }
#ifdef CORE_DEBUG
  virtual std::string op()
  { return std::string("ATan"); }
#endif 
};
/// \class Exp
/// \brief geometric power of e node
template <typename RootBd, typename Filter, typename Kernel>
class ExpRepT : public UnaryOpRepT<RootBd, Filter, Kernel> {
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;
  typedef UnaryOpRepT<RootBd, Filter, Kernel> UnaryOpRep;
  typedef RootBd* id_rootbd_t;
  using UnaryOpRep::child;
  using UnaryOpRep::check_exact;
  using ExprRep::filter;
  using ExprRep::sign;
  using ExprRep::uMSB;
  using ExprRep::lMSB;
  using ExprRep::appValue;
  using ExprRep::rootBd;
  using ExprRep::abs2rel;
  using ExprRep::init_value;
  using ExprRep::numType;
  using ExprRep::opType;
public:
  ExpRepT(ExprRep* c) : UnaryOpRep(c) {
    numType() = (std::max)(NODE_NT_TRANSCENDENTAL, child->numType());
    opType() = NODE_OP_TRANSCENDENTAL;
    if (child->get_sign()==0)
      init_value(1);
    child->a_approx(2);
  }
  virtual ~ExpRepT() 
  {}
protected:
  virtual bool compute_sign() 
  { sign() = 1; return true;}
  virtual bool compute_uMSB() {
    if (child->get_sign())
      uMSB() = child->appValue().getRight().get_si() * 2;
    else if(child->get_sign() < 0)
      uMSB() = child->appValue().getLeft().get_si();
    else
      uMSB() = 0;
    return 1;
  }
  virtual bool compute_lMSB() {
    if (child->get_sign() < 0)
      lMSB() = child->appValue().getRight().get_si() * 2;
    else if(child->get_sign() > 0)
      lMSB() = child->appValue().getLeft().get_si();
    else
      lMSB() = 0;
    return 1;
  }
  virtual bool compute_a_approx(prec_t prec) {
    prec_t L = prec + 2 * (abs(child->appValue().get_max()).get_ui() + 1);
    return check_exact(appValue().exp(child->a_approx(L), abs2rel(prec+1)));
  }
#ifdef CORE_DEBUG
  virtual std::string op()
  { return std::string("Exp"); }
#endif 
};
/// \class Exp2
/// \brief geometric power of 2 node
template <typename RootBd, typename Filter, typename Kernel>
class Exp2RepT : public UnaryOpRepT<RootBd, Filter, Kernel> {
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;
  typedef UnaryOpRepT<RootBd, Filter, Kernel> UnaryOpRep;
  typedef RootBd* id_rootbd_t;
  using UnaryOpRep::child;
  using UnaryOpRep::check_exact;
  using ExprRep::filter;
  using ExprRep::sign;
  using ExprRep::uMSB;
  using ExprRep::lMSB;
  using ExprRep::appValue;
  using ExprRep::rootBd;
  using ExprRep::abs2rel;
  using ExprRep::init_value;
  using ExprRep::numType;
  using ExprRep::opType;
public:
  Exp2RepT(ExprRep* c) : UnaryOpRep(c) {
    numType() = (std::max)(NODE_NT_TRANSCENDENTAL, child->numType());
    opType() = NODE_OP_TRANSCENDENTAL;
    if (child->get_sign()==0)
      init_value(1);
    child->a_approx(2);
  }
  virtual ~Exp2RepT() 
  {}
protected:
  virtual bool compute_sign() 
  { sign() = 1; return true;}
  virtual bool compute_uMSB() {
    if (child->get_sign())
      uMSB() = child->appValue().getRight().get_si() * 2;
    else if(child->get_sign() < 0)
      uMSB() = child->appValue().getLeft().get_si();
    else
      uMSB() = 0;
    return 1;
  }
  virtual bool compute_lMSB() {
    if (child->get_sign() < 0)
      lMSB() = child->appValue().getRight().get_si() * 2;
    else if(child->get_sign() > 0)
      lMSB() = child->appValue().getLeft().get_si();
    else
      lMSB() = 0;
    return 1;
  }
  virtual bool compute_a_approx(prec_t prec) {
    prec_t L = prec + 2 * (abs(child->appValue().get_max()).get_ui() + 1);
    return check_exact(appValue().exp2(child->a_approx(L), abs2rel(prec+1)));
  }
#ifdef CORE_DEBUG
  virtual std::string op()
  { return std::string("Exp2"); }
#endif 
};
/// \class Exp10
/// \brief geometric power of 10 node
template <typename RootBd, typename Filter, typename Kernel>
class Exp10RepT : public UnaryOpRepT<RootBd, Filter, Kernel> {
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;
  typedef UnaryOpRepT<RootBd, Filter, Kernel> UnaryOpRep;
  typedef RootBd* id_rootbd_t;
  using UnaryOpRep::child;
  using UnaryOpRep::check_exact;
  using ExprRep::filter;
  using ExprRep::sign;
  using ExprRep::uMSB;
  using ExprRep::lMSB;
  using ExprRep::appValue;
  using ExprRep::rootBd;
  using ExprRep::abs2rel;
  using ExprRep::init_value;
  using ExprRep::numType;
  using ExprRep::opType;
public:
  Exp10RepT(ExprRep* c) : UnaryOpRep(c) {
    numType() = (std::max)(NODE_NT_TRANSCENDENTAL, child->numType());
    opType() = NODE_OP_TRANSCENDENTAL;
    if (child->get_sign()==0)
      init_value(1);
    child->a_approx(2);
  }
  virtual ~Exp10RepT() 
  {}
protected:
  virtual bool compute_sign() 
  { sign() = 1; return true;}
  virtual bool compute_uMSB() {
    if (child->get_sign())
      uMSB() = child->appValue().getRight().get_si() * 2;
    else if(child->get_sign() < 0)
      uMSB() = child->appValue().getLeft().get_si();
    else
      uMSB() = 0;
    return 1;
  }
  virtual bool compute_lMSB() {
    if (child->get_sign() < 0)
      lMSB() = child->appValue().getRight().get_si() * 2;
    else if(child->get_sign() > 0)
      lMSB() = child->appValue().getLeft().get_si();
    else
      lMSB() = 0;
    return 1;
  }
  virtual bool compute_a_approx(prec_t prec) {
    prec_t L = prec + 2 * (abs(child->appValue().get_max()).get_ui() + 1);
    return check_exact(appValue().exp10(child->a_approx(L), abs2rel(prec+1)));
  }
#ifdef CORE_DEBUG
  virtual std::string op()
  { return std::string("Exp10"); }
#endif 
};
/// \class Log2
/// \brief geometric logarithm base 2
template <typename RootBd, typename Filter, typename Kernel>
class Log2RepT : public UnaryOpRepT<RootBd, Filter, Kernel> {
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;
  typedef UnaryOpRepT<RootBd, Filter, Kernel> UnaryOpRep;
  typedef RootBd* id_rootbd_t;
  using UnaryOpRep::child;
  using UnaryOpRep::check_exact;
  using ExprRep::filter;
  using ExprRep::sign;
  using ExprRep::uMSB;
  using ExprRep::lMSB;
  using ExprRep::appValue;
  using ExprRep::rootBd;
  using ExprRep::abs2rel;
  using ExprRep::init_value;
  using ExprRep::numType;
  using ExprRep::opType;
public:
  Log2RepT(ExprRep* c) : UnaryOpRep(c) {
    numType() = (std::max)(NODE_NT_TRANSCENDENTAL, child->numType());
    opType() = NODE_OP_TRANSCENDENTAL;
    if (child->get_sign()<=0)
      core_error("log2 of negative value", __FILE__, __LINE__, true);
    child->a_approx(2);
  }
  virtual ~Log2RepT() 
  {}
protected:
  virtual bool compute_sign() {
/*    if (child->appValue().get_min() > 1)
    { sign() = 1; return true; }
    else if (child->appValue().get_max() < 1)
    { sign() = -1; return true; }
    else if (child->appValue().is_exact() && child->appValue().get_f() == 1)
    { sign() = 0; return true; }
    else
*/      return false;
  }
  virtual bool compute_uMSB() 
  { uMSB() = ceillg(child->get_uMSB()); return true; }
  virtual bool compute_lMSB() 
  { lMSB() = floorlg(child->get_lMSB()); return true; }
  virtual bool compute_a_approx(prec_t prec) {
    return check_exact(appValue().log2(child->a_approx(prec+1-floorlg(child->a_approx(2).get_min())), abs2rel(prec+1)));
  }
#ifdef CORE_DEBUG
  virtual std::string op()
  { return std::string("log2"); }
#endif 
};

/// \class Log
/// \brief geometric logarithm base e
template <typename RootBd, typename Filter, typename Kernel>
class LogRepT : public UnaryOpRepT<RootBd, Filter, Kernel> {
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;
  typedef UnaryOpRepT<RootBd, Filter, Kernel> UnaryOpRep;
  typedef RootBd* id_rootbd_t;
  using UnaryOpRep::child;
  using UnaryOpRep::check_exact;
  using ExprRep::filter;
  using ExprRep::sign;
  using ExprRep::uMSB;
  using ExprRep::lMSB;
  using ExprRep::appValue;
  using ExprRep::rootBd;
  using ExprRep::abs2rel;
  using ExprRep::init_value;
  using ExprRep::numType;
  using ExprRep::opType;
public:
  LogRepT(ExprRep* c) : UnaryOpRep(c) {
    numType() = (std::max)(NODE_NT_TRANSCENDENTAL, child->numType());
    opType() = NODE_OP_TRANSCENDENTAL;
    if (child->get_sign()<=0)
      core_error("log2 of negative value", __FILE__, __LINE__, true);
    child->a_approx(2);
  }
  virtual ~LogRepT() 
  {}
protected:
  virtual bool compute_sign() {
/*    if (child->get_sign() > 1)
    { sign() = 1; return true; }
    else if (child->get_sign() < 1)
    { sign() = -1; return true; }
    else if (child->appValue().is_exact() && child->appValue().get_f() == 1)
    { sign() = 0; return true; }
    else
*/      return false;
  }
  virtual bool compute_uMSB() 
  { uMSB() = ceillg(child->get_uMSB()); return true; }
  virtual bool compute_lMSB() 
  { lMSB() = floorlg(child->get_lMSB()); return true; }
  virtual bool compute_a_approx(prec_t prec) {
    return check_exact(appValue().log(child->a_approx(prec+1-floorlg(child->a_approx(2).get_min())), abs2rel(prec+1)));
  }
#ifdef CORE_DEBUG
  virtual std::string op()
  { return std::string("log"); }
#endif 
};

/// \class Log10
/// \brief geometric logarithm base 10
template <typename RootBd, typename Filter, typename Kernel>
class Log10RepT : public UnaryOpRepT<RootBd, Filter, Kernel> {
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;
  typedef UnaryOpRepT<RootBd, Filter, Kernel> UnaryOpRep;
  typedef RootBd* id_rootbd_t;
  using UnaryOpRep::child;
  using UnaryOpRep::check_exact;
  using ExprRep::filter;
  using ExprRep::sign;
  using ExprRep::uMSB;
  using ExprRep::lMSB;
  using ExprRep::appValue;
  using ExprRep::rootBd;
  using ExprRep::abs2rel;
  using ExprRep::init_value;
  using ExprRep::numType;
  using ExprRep::opType;
public:
  Log10RepT(ExprRep* c) : UnaryOpRep(c) {
    numType() = (std::max)(NODE_NT_TRANSCENDENTAL, child->numType());
    opType() = NODE_OP_TRANSCENDENTAL;
    if (child->get_sign()<=0)
      core_error("log2 of negative value", __FILE__, __LINE__, true);
    child->a_approx(2);
  }
  virtual ~Log10RepT() 
  {}
protected:
  virtual bool compute_sign() {
/*    if (child->appValue().get_min() > 1)
    { sign() = 1; return true; }
    else if (child->appValue().get_max() < 1)
    { sign() = -1; return true; }
    else if (child->appValue().is_exact() && child->appValue().get_f() == 1)
    { sign() = 0; return true; }
    else
*/      return false;
  }
  virtual bool compute_uMSB() 
  { uMSB() = ceillg(child->get_uMSB()); return true; }
  virtual bool compute_lMSB() 
  { lMSB() = floorlg(child->get_lMSB()); return true; }
  virtual bool compute_a_approx(prec_t prec) {
    return check_exact(appValue().log10(child->a_approx(prec+1-floorlg(child->a_approx(2).get_min())), abs2rel(prec+1)));
  }
#ifdef CORE_DEBUG
  virtual std::string op()
  { return std::string("log10"); }
#endif 
};

/// \class RootRepT
/// \brief k-th root node
template <typename RootBd, typename Filter, typename Kernel>
class RootRepT : public UnaryOpRepT<RootBd, Filter, Kernel> {
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;
  typedef UnaryOpRepT<RootBd, Filter, Kernel> UnaryOpRep;
  typedef RootBd* id_rootbd_t;
  using UnaryOpRep::child;
  using UnaryOpRep::check_exact;
  using ExprRep::filter;
  using ExprRep::sign;
  using ExprRep::uMSB;
  using ExprRep::lMSB;
  using ExprRep::appValue;
  using ExprRep::rootBd;
  using ExprRep::abs2rel;
  using ExprRep::init_value;
  using ExprRep::numType;
  using ExprRep::opType;
public:
  RootRepT(ExprRep* c, unsigned long k) : UnaryOpRep(c), m_k(k) 
  { if (k == 0)
      core_error("0-th root is undefined", __FILE__, __LINE__, true);
    if (k%2 == 0 && c->get_sign() < 0)
      core_error("even root of negative value", __FILE__, __LINE__, true);
    filter().root(child->filter(), k);
    if (child->get_sign()==0)
      init_value(0);
    numType() = (std::max)(NODE_NT_ALGEBRAIC, child->numType());
    opType() = NODE_OP_ALGEBRAIC;
  }
  virtual ~RootRepT() 
  {}
protected:
  virtual bool compute_sign() 
  { sign() = child->get_sign(); return true;}
  virtual bool compute_uMSB() 
  { uMSB() = (child->get_uMSB()+m_k-1)/m_k; return true;}
  virtual bool compute_lMSB() 
  { lMSB() = child->get_lMSB()/m_k; return true;}
  virtual void compute_rootBd()
  { rootBd().root(child->get_rootBd(), m_k); }
  virtual bool compute_r_approx(prec_t prec) {
    return check_exact(appValue().root(child->r_approx(prec*m_k),
			    m_k, prec));
  }
#ifdef CORE_DEBUG
  virtual std::string op()
  { return std::string("Radical"); }
#endif 
private:
  unsigned long m_k;
};

/// \class BinaryOpRepT
/// \brief binary operation expression node
template <typename RootBd, typename Filter, typename Kernel>
class BinaryOpRepT : public ExprRepT<RootBd, Filter, Kernel> {
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;
public:
  BinaryOpRepT(ExprRep* f, ExprRep* s, bool is_self) : first(f), second(s)
  { if (!is_self) first->inc_ref(); second->inc_ref(); }
  virtual ~BinaryOpRepT()
  { first->dec_ref(); second->dec_ref(); }
  virtual size_t get_children_size() const
  { return 2; }
  virtual ExprRep* get_child(size_t i) const
  { return i==0 ? first : second; } 
  /// check whether the operation and child nodes are both exact
  bool check_exact(bool ret)
  { return ret && first->is_exact() && second->is_exact(); }
protected:
  ExprRep* first;  /// <- pointer to the first child node
  ExprRep* second; /// <- pointer to the second child node
};// BinaryOpRepT

/// \class AddSubRepT
/// \brief add/sub node -- this is the most critical type of node, and
/// combines addition and subtraction.
template <typename RootBd, typename Filter, typename Kernel, bool is_add>
class AddSubRepT : public BinaryOpRepT<RootBd, Filter, Kernel> {
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;
  typedef BinaryOpRepT<RootBd, Filter, Kernel> BinaryOpRep;
  typedef RootBd* id_rootbd_t;

  using BinaryOpRep::check_exact;	// This line added Apr24'2013 (Chee)
  using BinaryOpRep::first; 
  using BinaryOpRep::second; 
  using ExprRep::filter;
  using ExprRep::rootBd;
  using ExprRep::sign;
  using ExprRep::uMSB;
  using ExprRep::lMSB;
  using ExprRep::appValue;
  using ExprRep::abs2rel;
  using ExprRep::set_flags;
  using ExprRep::numType;
  using ExprRep::opType;
  
//  static int debugVariable;
public:
  AddSubRepT(ExprRep* f, ExprRep* s, bool b = false) : BinaryOpRep(f, s, b) {
    filter().addsub(first->filter(),second->filter(), is_add);
    numType() = (std::max)(first->numType(), second->numType());
    opType() = is_add ? NODE_OP_ADD : NODE_OP_SUB;
  }
  virtual ~AddSubRepT() 
  {}
protected:
  // precondition : this is only called after filter fails, approx value fails, and no cache.
  // The goal of compute_sign() is to force the sign computation of the children, and then if
  //    this does not determine the current sign, we try to use the
  //    children's lMSB and uMSB, without doing refine().   So compute_sign() COULD fail and return false.
  //    If it succeeds, it will store the sign in cache.
  virtual bool compute_sign() {
    sign_t sf = first->get_sign();
    sign_t ss = second->get_sign();
    if (!is_add) ss = -ss;
    if (sf == 0) // first operand is zero
      sign() = ss;
    else if (ss == 0) // second operand is zero
      sign() = sf;
    else if (sf == ss) // same sign
      sign() = sf;
    else { // different sign
     if (first->get_lMSB() > second->get_uMSB())
        sign() = sf;
      else if (first->get_uMSB() < second->get_lMSB())
        sign() = ss ;
      else // unknown sign!
        return false;
    }
    return true;
  }
  virtual bool compute_uMSB() {
 // sign is too expensive and should be avoided at all cost:
  /*
    sign_t sf = first->get_sign();
    sign_t ss = second->get_sign();
    if (!is_add) ss = -ss;
    if (sf == 0) // first operand is zero
      uMSB() = second->get_uMSB();
    else if (ss == 0) // second operand is zero
      uMSB() = first->get_uMSB();
    else {
      msb_t uf = first->get_uMSB();
      msb_t us = second->get_uMSB();
      uMSB() = (std::max)(uf, us) + 1;
      if (sf == ss) uMSB() += 1;
    }
  /*/
    msb_t uf = first->get_uMSB();
    msb_t us = second->get_uMSB();
    uMSB() = (std::max)(uf, us) + 1;
  //*/
    return true;
  }
  virtual bool compute_lMSB() {
// Sep 8, 2006: Jihun/Chee rewrote to avoid computing signs as long as possible
// Jan 18, 2010 : Jihun/Chee we roll back to the old implementation, because 
// testIdent fails when we turn off the filter and compute_lMSB goes wrong.
    /*
    msb_t lf = first->get_lMSB();
    msb_t ls = second->get_lMSB();
    if (lf > second->get_uMSB()+1) { // note that we do need the "+1" to ensure a factor of 2 gap
      lMSB() = lf - 1; return true;
    } 
    if (ls > first->get_uMSB()+1) {
      lMSB() = ls - 1; return true;
    }
    sign_t sf = first->get_sign();
    if (sf == 0) {// first operand is zero
      lMSB() = second->get_lMSB(); return true;
    }
    sign_t ss = second->get_sign();
	// In core prelease, this next line was missing. 
	if(!is_add) ss = -ss;
    if (ss == 0) {// second operand is zero
      lMSB() = first->get_lMSB(); return true;
    }
   if (sf == ss) {// same sign
     if (lf == ls)
       lMSB() = lf + 1;
     else
       lMSB() = (std::max)(lf, ls);
   } else { // different sign
       return false; // failed to compute lMSB if different sign and we cannot
       			// decide if first->MSB is different from second-MSB
   }
   /*/
    // The following code is right, but inefficient. Keep for debugging:
    sign_t sf = first->get_sign();
    sign_t ss = second->get_sign();
    if (!is_add) ss = -ss;
    if (sf == 0) // first operand is zero
      lMSB() = second->get_lMSB();
    else if (ss == 0) // second operand is zero
      lMSB() = first->get_lMSB();
    else {
      msb_t lf = first->get_lMSB();
      msb_t ls = second->get_lMSB();
      if (sf == ss) {// same sign
        if (lf == ls)
          lMSB() = lf + 1;
        else
          lMSB() = (std::max)(lf, ls);
      } else { // different sign
        if (lf > second->get_uMSB())
          lMSB() = lf - 1;
        else if (ls > first->get_uMSB())
          lMSB() = ls - 1;
        else // unknown lMSB()!
          return false;
      }
    }
    //*/
    return true;
  } //compute_lMSB

  virtual void compute_rootBd()
  {rootBd().addsub(first->get_rootBd(),second->get_rootBd());}
  virtual bool compute_a_approx(prec_t prec) {
    return this->check_exact(appValue().addsub(first->a_approx(prec+2),
	  second->a_approx(prec+2), abs2rel(prec+1), is_add)); 
  }
  BigInt getZTVal() { return first->getZTVal() + second->getZTVal() * (is_add ? 1 : -1); }
  BigFloat getFTVal() { return first->getFTVal() + second->getFTVal() * (is_add ? 1 : -1); }
  BigRat getQTVal() { return first->getQTVal() + second->getQTVal() * (is_add ? 1 : -1); }
#ifdef CORE_DEBUG
  virtual std::string op()
  { return std::string("AddSub"); }
#endif 
}; //AddSubRepT

/// \class MulRepT
/// \brief multiplication node
template <typename RootBd, typename Filter, typename Kernel>
class MulRepT : public BinaryOpRepT<RootBd, Filter, Kernel> {
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;
  typedef BinaryOpRepT<RootBd, Filter, Kernel> BinaryOpRep;
  typedef RootBd* id_rootbd_t;
  //using BinaryOpRep::check_exact;	// This line added Apr24'2013 (Chee)
  using BinaryOpRep::check_exact;
  using BinaryOpRep::first; 
  using BinaryOpRep::second; 
  using ExprRep::filter;
  using ExprRep::rootBd;
  using ExprRep::sign;
  using ExprRep::uMSB;
  using ExprRep::lMSB;
  using ExprRep::appValue;
  using ExprRep::numType;
  using ExprRep::opType;
public:
  MulRepT(ExprRep* f, ExprRep* s, bool b = false) : BinaryOpRep(f, s, b) {
    filter().mul(first->filter(), second->filter()); 
    numType() = (std::max)(first->numType(), second->numType());
    opType() = NODE_OP_MUL;
  }  
  virtual ~MulRepT() 
  {}
protected:
  virtual bool compute_sign() 
  { sign() = first->get_sign() * second->get_sign(); return true;}
  virtual bool compute_uMSB() 
  { uMSB() = first->get_uMSB() + second->get_uMSB(); return true;}
  virtual bool compute_lMSB() 
  { lMSB() = first->get_lMSB() + second->get_lMSB(); return true;}
  virtual void compute_rootBd()
  { rootBd().mul(first->get_rootBd(), second->get_rootBd()); }
  virtual bool compute_r_approx(prec_t prec) {
   return check_exact(appValue().mul(first->r_approx(prec+2),
			   second->r_approx(prec+2), prec+1));
  }
  BigInt getZTVal() { return first->getZTVal() * second->getZTVal(); }
  BigFloat getFTVal() { return first->getFTVal() * second->getFTVal(); }
  BigRat getQTVal() {
	  std::cerr << "QTVal=" << (first->getQTVal() * second->getQTVal())
		  	<< std::endl;
	  std::cerr << "first=" << first->getQTVal()
		  	<< std::endl;
	  std::cerr << "second=" << second->getQTVal()
		  	<< std::endl;
	  return first->getQTVal() * second->getQTVal(); }
#ifdef CORE_DEBUG
  virtual std::string op()
  { return std::string("Mul"); }
#endif 
};

/// \class DivRepT
/// \brief division node
template <typename RootBd, typename Filter, typename Kernel>
class DivRepT : public BinaryOpRepT<RootBd, Filter, Kernel> {
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;
  typedef BinaryOpRepT<RootBd, Filter, Kernel> BinaryOpRep;
  typedef RootBd* id_rootbd_t;
  using BinaryOpRep::check_exact;	// This line added Apr24'2013 (Chee)
  //using BinaryOpRep::check_exact;	// 
  using BinaryOpRep::first; 
  using BinaryOpRep::second; 
  using ExprRep::filter;
  using ExprRep::rootBd;
  using ExprRep::sign;
  using ExprRep::uMSB;
  using ExprRep::lMSB;
  using ExprRep::appValue;
  using ExprRep::numType;
  using ExprRep::opType;
public:
  DivRepT(ExprRep* f, ExprRep* s, bool b = false) : BinaryOpRep(f, s, b)
  { filter().div(first->filter(), second->filter());
    numType() = (std::max)(NODE_NT_RATIONAL, (std::max)(first->numType(), second->numType()));
    opType() = NODE_OP_DIV;
    if (second->get_sign() == 0)
      core_error("divide by zero", __FILE__, __LINE__, true);
  }
  virtual ~DivRepT() 
  {}
protected:
  virtual bool compute_sign() 
  { sign() = first->get_sign() * second->get_sign(); return true;}
  virtual bool compute_uMSB() 
  { uMSB() = first->get_uMSB() - second->get_lMSB(); return true;}
  virtual bool compute_lMSB() 
  { lMSB() = first->get_lMSB() - second->get_uMSB(); return true;}
  virtual void compute_rootBd()
  { rootBd().div(first->get_rootBd(), second->get_rootBd()); }
  virtual bool compute_r_approx(prec_t prec) {
    return this->check_exact(appValue().div(first->r_approx(prec+2),
			    second->r_approx(prec+2), prec+1));
    // Chee: April 2013: added "this->" to check_exact to try to remove warning.
  }
  BigInt getZTVal() { return first->getZTVal() / second->getZTVal(); }
  BigFloat getFTVal() { return first->getFTVal() / second->getFTVal(); }
  BigRat getQTVal() { return first->getQTVal() / second->getQTVal(); }
#ifdef CORE_DEBUG
  virtual std::string op()
  { return std::string("Div"); }
#endif 
};

/// \class AnaryOpRepT
/// \brief k-nary operation expression node
template <typename RootBd, typename Filter, typename Kernel>
class AnaryOpRepT : public ExprRepT<RootBd, Filter, Kernel> {
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;
  using ExprRep::opType;
public:
  AnaryOpRepT() {}
  AnaryOpRepT(const std::vector<ExprRep*>& c, bool is_self = false) : children(c)
  {
    opType() = NODE_OP_ANARY;
    for (size_t i=0; i<children.size(); i++) {
      children[i]->inc_ref();
	}
	if (children.size()>0 && is_self) children[0]->dec_ref();
  }
  virtual ~AnaryOpRepT() {
    for (size_t i=0; i<children.size(); i++) 
      children[i]->dec_ref();
  }
  virtual size_t get_children_size() const
  { return children.size(); }
  virtual ExprRep* get_child(size_t i) const
  { return children[i]; } 
protected:
  std::vector<ExprRep*> children;  /// <- vector of pointers to children nodes
};

/// \class SumRepT
/// \brief summation node
template <typename RootBd, typename Filter, typename Kernel>
class SumOpRepT : public AnaryOpRepT<RootBd, Filter, Kernel> {
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;
  typedef SumOpRepT<RootBd, Filter, Kernel> SumOpRep;
  typedef AddSubRepT<RootBd, Filter, Kernel, true> AddOpRep;
  typedef AnaryOpRepT<RootBd, Filter, Kernel> AnaryOpRep;
  typedef RootBd* id_rootbd_t;
  using AnaryOpRep::children;
  using ExprRep::filter;
  using ExprRep::rootBd;
  using ExprRep::sign;
  using ExprRep::uMSB;
  using ExprRep::lMSB;
  using ExprRep::appValue;
  using ExprRep::abs2rel;
  using ExprRep::numType;
public:
  SumOpRepT(ExprRep* first, ExprRep* second, bool is_self = false) {
    insert(first, is_self);
    insert(second);
  }    
    
  SumOpRepT(const std::vector<ExprRep*>& c = std::vector<ExprRep*>(), bool is_self = false)
  : AnaryOpRep(c, is_self)  {
    compute_filter();
    compute_numtype();
  }
  virtual ~SumOpRepT() {}
  
  void insert(ExprRep* c, bool is_self = false) {
    children.push_back(c);
    
    if (!is_self) c->inc_ref();
    
    filter().addsub(filter(), c->filter(), true);
    numType() = (std::max)(numType(), c->numType());
  }
protected:
  void compute_filter () {
    filter().set(0L);
    for(size_t i=0; i<children.size(); i++)
      filter().addsub(filter(), children[i]->filter(), true);
  }
  void compute_numtype () {
    numType() = NODE_NT_INTEGER;
    for(size_t i=0; i<children.size(); i++)
      numType() = (std::max)(numType(), children[i]->numType());
  }

  virtual bool compute_sign() {
    return false;
  }
  virtual bool compute_uMSB() {
    uMSB() = children[0]->get_uMSB();
    for(size_t i = 1; i < children.size(); i++) {
      uMSB() = (std::max)(uMSB(), children[i]->get_uMSB()) + 1;
    }
    return true;
  }
  virtual bool compute_lMSB() {
    return false;
  }
  virtual void compute_rootBd() {
    rootBd().addsub(children[0]->get_rootBd(), children[1]->get_rootBd());
    for (size_t i = 2; i < children.size(); i++) {
      rootBd().addsub(rootBd(), children[i]->get_rootBd());
    }
  }
  virtual bool compute_a_approx(prec_t prec) {
    int n = ceillg(children.size()) + 1;
    Kernel Value;
    appValue().set(0);

    bool exact = true;
    for (size_t i=0; i<children.size(); i++) {
       bool r = Value.add(appValue(), children[i]->a_approx(prec+n), abs2rel(prec+n));
       if (!children[i]->is_exact() || !r) exact = false;
       Value.swap(appValue());
    }
    return exact;
  }
#ifdef CORE_DEBUG
  virtual std::string op()
  { return std::string("Sum"); }
#endif 
public:
  bool is_self_mergable(prec_t absprec = DEF_INIT_PREC) {
    if (this->get_ref() > 1) return false;
    if (!is_approx_needed(abs2rel(absprec))) return false;

    return true;
  }
}; 

/// \class Product
/// \brief Product node
template <typename RootBd, typename Filter, typename Kernel>
class ProdOpRepT : public AnaryOpRepT<RootBd, Filter, Kernel> {
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;
  typedef AnaryOpRepT<RootBd, Filter, Kernel> AnaryOpRep;
  typedef RootBd* id_rootbd_t;
  using AnaryOpRep::children;
  using ExprRep::filter;
  using ExprRep::sign;
  using ExprRep::uMSB;
  using ExprRep::lMSB;
  using ExprRep::appValue;
  using ExprRep::rootBd;
  using ExprRep::numType;
public:
  ProdOpRepT(const std::vector<ExprRep*>& c = std::vector<ExprRep*>(), bool is_self = false)
  : AnaryOpRep(c, is_self)  {
    compute_filter();
    compute_numtype();  
  }
  virtual ~ProdOpRepT() {}
  
  void insert(ExprRep* c, bool is_self = false) {
    children.push_back(c);
    if (!is_self) c->inc_ref();
    filter().mul(filter(), c->filter());
  }
protected:
  void compute_filter () {
    filter().set(1L);
    for(size_t i=0; i<children.size(); i++)
      filter().mul(filter(), children[i]->filter());
  }
  void compute_numtype () {
    numType() = NODE_NT_INTEGER;
    for(size_t i=0; i<children.size(); i++)
      numType() = (std::max)(numType(), children[i]->numType());
  }

  virtual bool compute_sign() {
    sign_t tmpsgn = children[0]->get_sign();
    for(size_t i=1; i<children.size(); i++) {
      tmpsgn *= children[i]->get_sign();
    }
    sign() = tmpsgn;
    return true;
  }
  virtual bool compute_uMSB() {
    msb_t tmpuMSB = children[0]->get_uMSB();
    for(size_t i=1; i<children.size(); i++) {
      tmpuMSB += children[i]->get_uMSB();
    }
    uMSB() = tmpuMSB;
    return true;
  }
  virtual bool compute_lMSB() {
    msb_t tmplMSB = children[0]->get_lMSB();
    for(size_t i=1; i<children.size(); i++) {
      tmplMSB += children[i]->get_lMSB();
    }
    lMSB() = tmplMSB;
    return true;
  }
  virtual void compute_rootBd() {
    rootBd().mul(children[0]->get_rootBd(), children[1]->get_rootBd());
    for (size_t i=2; i<children.size(); i++) {
      rootBd().mul(rootBd(), children[i]->get_rootBd());
    }
  }
  virtual bool compute_r_approx(prec_t prec) {
    int n = ceillg(children.size()) + 1;
    Kernel Value;
    appValue().set(1);
    
    bool exact = true;
    for (size_t i=0; i<children.size(); i++) {
       bool r = Value.mul(appValue(), children[i]->r_approx(prec+n), prec+n);
       if (!children[i]->is_exact() || !r) exact = false;
       Value.swap(appValue());
    }
    return exact;
  }
#ifdef CORE_DEBUG
  virtual std::string op()
  { return std::string("Prod"); }
#endif 
}; 
/*
/// \class NegRepT macro version
/// \brief negation node
BEGIN_DEFINE_UNARY_NODE(NegRepT)
  BEGIN_DEFINE_RULE_FILTER
    filter().neg(child->m_filter);
  END_DEFINE_RULE

  BEGIN_DEFINE_RULE_NUMTYPE
    this->numType() = this->child->numType();
  END_DEFINE_RULE

  BEGIN_DEFINE_RULE_SIGN
    sign() = -child->get_sign(); return true;
  END_DEFINE_RULE

  BEGIN_DEFINE_RULE_UMSB
    uMSB() = child->get_uMSB(); return true;
  END_DEFINE_RULE

  BEGIN_DEFINE_RULE_LMSB
    lMSB() = child->get_lMSB(); return true;
  END_DEFINE_RULE

  BEGIN_DEFINE_RULE_R_APPROX
    return check_exact(appValue().neg(child->r_approx(prec), prec));
  END_DEFINE_RULE
    
  BEGIN_DEFINE_RULE_A_APPROX
    return check_exact(appValue().neg(child->a_approx(prec), abs2rel(prec)));
  END_DEFINE_RULE

  BEGIN_DEFINE_RULE_ROOTBD
    rootBd().neg(child->get_rootBd());
  END_DEFINE_RULE

END_DEFINE_UNARY_NODE

/// \class DivRepT macro version
/// \brief division node
BEGIN_DEFINE_BINARY_NODE(DivRepT)
  BEGIN_DEFINE_RULE_FILTER
    filter().div(first->filter(), second->filter());
  END_DEFINE_RULE

  BEGIN_DEFINE_RULE_NUMTYPE
    this->numType() = this->child->numType();
  END_DEFINE_RULE

  BEGIN_DEFINE_RULE_SIGN
    sign() = first->get_sign() * second->get_sign(); return true;
  END_DEFINE_RULE

  BEGIN_DEFINE_RULE_UMSB
    uMSB() = first->get_uMSB() - second->get_lMSB(); return true;
  END_DEFINE_RULE

  BEGIN_DEFINE_RULE_LMSB
    lMSB() = first->get_lMSB() - second->get_uMSB(); return true;
  END_DEFINE_RULE

  BEGIN_DEFINE_RULE_ROOTBD
    rootBd().div(first->get_rootBd(), second->get_rootBd());
  END_DEFINE_RULE

  BEGIN_DEFINE_RULE_R_APPROX
    return check_exact(appValue().div(first->r_approx(prec+2), second->r_approx(prec+2), prec+1));
  END_DEFINE_RULE
END_DEFINE_UNARY_NODE
*/

CORE_END_NAMESPACE

#endif /*__CORE_EXPRREP_H__*/
