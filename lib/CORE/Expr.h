/****************************************************************************
 * Expr.h -- EGC number class providing guarranteed precision
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
 * $Id: Expr.h,v 1.63 2010/11/23 17:58:36 exact Exp $
 ***************************************************************************/

/* **************************************************
 * COMMENTS on Expr.h:
 *
 * 	The main class here is a templated class called ExprT.
 * 	It takes 3 classes as template arguments:
 * 		ExprT< RootBound, Filter, Kernel >
 * 	where
 * 		RootBound computes root bounds for expressions
 * 		Filter provides a cheap floating point filter
 * 		Kernel is a number type for approximations
 *
 * 	We define Expr as a typedef (see line 550 below).
 * 	E.g.,
 *          typedef ExprT <   BfmsskRootBd_double<BigFloat2>,
 *                            BfsFilter<BigFloat2>,
 *                            BigFloat2
 *                        >   Expr;
 *
 *      Note that the RootBound and Filter classes are
 *      themselves templated, and must have the Kernel class as
 *      template argument.
 *
 * In the Core 2.0,
 *	--the only Kernel available is BigFloat2
 * 	--we have two Filters: BfsFilter or the DummyFilter (i.e., no filter).
 *	--we have two main variants of filters, either BFMSS or k-BFMSS
 *		(the latter is the k-ary version of BFMSS)
 *
 *  	 In 2010, we debugged the BFMSS and kary-BFMSS bounds because the
 *       the root bounds were originally stored as machine longs, and they
 *       could overflow and lead to error.   We experimented with two
 *       alternatives for machine longs: either machine doubles or BigFloats.
 *       This gave rise to variants of the BFMSS and k-BFMSS bounds,
 *	 are qualified by "_double" or "_BigFloat".  Although _double
 *	 is not full-proof from overflow, it seems sufficient to overcome the
 *	 errors encountered by long.  Although the tightest bounds are
 *	 provided by _BigFloat, we deem _double a good compromise.
 *
 * -- Jihun/Chee ( June 2010)
 *
 * TODO: Unfortunately, other root bounds which were previously
 *       implemented in Core 1.7 (Li-Yap Bound, Mahler Measure)
 *       have not been implemented.
 * **************************************************/

#ifndef __CORE_EXPR_H__
#define __CORE_EXPR_H__

#include <CORE/ExprRep.h>

CORE_BEGIN_NAMESPACE

/// \class ExprT
/// Kernel -- internal representation 
template <typename RootBd, typename Filter, typename Kernel>
class ExprT {
public: // public typedefs
  typedef RootBd     RootBdT;
  typedef Filter     FilterT;
  typedef Kernel     KernelT;

  typedef typename Kernel::ZT ZT;
  typedef typename Kernel::QT QT;
  typedef typename Kernel::FT FT;
  typedef Kernel              KT;

private: // private typedefs
  typedef ExprT<RootBd, Filter, Kernel> thisClass;
  typedef ExprRepT<RootBd, Filter, Kernel> ExprRep;
  // zero ary operators
  typedef ConstRepT<RootBd, Filter, Kernel, long> ConstLongRep;
  typedef ConstRepT<RootBd, Filter, Kernel, unsigned long> ConstULongRep;
  typedef ConstRepT<RootBd, Filter, Kernel, double> ConstDoubleRep;
  typedef ConstRepT<RootBd, Filter, Kernel, ZT> ConstZTRep;
  typedef ConstRepT<RootBd, Filter, Kernel, QT> ConstQTRep;
  typedef ConstRepT<RootBd, Filter, Kernel, FT> ConstFTRep;
  typedef PiRepT<RootBd, Filter, Kernel> PiRep;
  typedef ERepT<RootBd, Filter, Kernel> ERep;
  //typedef ConstRepT<RootBd, Filter, Kernel, KT> ConstKTRep;
  // binary operators
  typedef AddSubRepT<RootBd, Filter, Kernel, true> AddRep;
  typedef AddSubRepT<RootBd, Filter, Kernel, false> SubRep;
  typedef MulRepT<RootBd, Filter, Kernel> MulRep;
  typedef DivRepT<RootBd, Filter, Kernel> DivRep;
  // unary operators
  typedef NegRepT<RootBd, Filter, Kernel> NegRep;
  typedef SqrtRepT<RootBd, Filter, Kernel> SqrtRep;
  typedef CbrtRepT<RootBd, Filter, Kernel> CbrtRep;
  typedef RootRepT<RootBd, Filter, Kernel> RootRep;
  // transcendental unary operators
  typedef SinRepT<RootBd, Filter, Kernel> SinRep;
  typedef CosRepT<RootBd, Filter, Kernel> CosRep;
  typedef TanRepT<RootBd, Filter, Kernel> TanRep;
  typedef CotRepT<RootBd, Filter, Kernel> CotRep;
  typedef ASinRepT<RootBd, Filter, Kernel> ASinRep;
  typedef ACosRepT<RootBd, Filter, Kernel> ACosRep;
  typedef ATanRepT<RootBd, Filter, Kernel> ATanRep;
  typedef LogRepT<RootBd, Filter, Kernel> LogRep;
  typedef Log2RepT<RootBd, Filter, Kernel> Log2Rep;
  typedef Log10RepT<RootBd, Filter, Kernel> Log10Rep;
  typedef ExpRepT<RootBd, Filter, Kernel> ExpRep;
  typedef Exp2RepT<RootBd, Filter, Kernel> Exp2Rep;
  typedef Exp10RepT<RootBd, Filter, Kernel> Exp10Rep;
  // anary operators
  typedef SumOpRepT<RootBd, Filter, Kernel> SumOpRep;
  typedef ProdOpRepT<RootBd, Filter, Kernel> ProdOpRep;
  typedef SumOpRepT<RootBd, Filter, Kernel> SumRep;
  typedef ProdOpRepT<RootBd, Filter, Kernel> ProdRep;
  // TO DO:
  // DiamondRep (generalization of RootOfRep?)
  // DeterminantRep

protected:
  ExprRep* GetZeroConst() {
    static ConstLongRep zero(0L, NODE_NT_INTEGER);
    zero.inc_ref();
    return &zero;
  }

public:
  ExprT() : m_rep(GetZeroConst()) { }
  ExprT(char c) : m_rep(new ConstLongRep(long(c), NODE_NT_INTEGER)) {}
  ExprT(unsigned char uc) : m_rep(new ConstULongRep((unsigned long)(uc), NODE_NT_INTEGER)) {}
  ExprT(short s) : m_rep(new ConstLongRep(long(s), NODE_NT_INTEGER)) {}
  ExprT(unsigned short us) : m_rep(new ConstULongRep((unsigned long)(us), NODE_NT_INTEGER)) {}
  ExprT(int i) : m_rep(new ConstLongRep(long(i), NODE_NT_INTEGER)) {}
  ExprT(unsigned int ui) : m_rep(new ConstULongRep((unsigned long)(ui), NODE_NT_INTEGER)) {}
  ExprT(long l) : m_rep(new ConstLongRep(l, NODE_NT_INTEGER)) {}
  ExprT(unsigned long ul) : m_rep(new ConstULongRep(ul, NODE_NT_INTEGER)) {}
  ExprT(float f) : m_rep(new ConstDoubleRep(double(f), NODE_NT_DYADIC)) {}
  ExprT(double d) : m_rep(new ConstDoubleRep(d, NODE_NT_DYADIC)) {}
  ExprT(const ZT& z) : m_rep(new ConstZTRep(z, NODE_NT_INTEGER)) {}
  ExprT(const FT& f) : m_rep(new ConstFTRep(f, NODE_NT_DYADIC)) {}
  ExprT(const QT& q) { 
    FT f;
    (f.set(q)==0)?
    (m_rep=new ConstFTRep(f, NODE_NT_DYADIC)):
    (m_rep=new ConstQTRep(q, NODE_NT_RATIONAL));
  }
  ExprT(const KT& k) : m_rep(new ConstFTRep(k.get_f(), NODE_NT_DYADIC)) {}
  ExprT(const char* s, int base = 10, prec_t prec = digits2bits(getDefaultInputDigits())) 
  { construct_from_string(s, base, prec); }
  ExprT(const std::string& s, int base = 10, prec_t prec = digits2bits(getDefaultInputDigits())) 
  { construct_from_string(s.c_str(), base, prec); }
public:
  ExprT(ExprRep* r) : m_rep(r) {}
  ExprT(const ExprT& r) : m_rep(r.m_rep) { m_rep->inc_ref(); }
  ~ExprT() { m_rep->dec_ref(); }
  ExprT& operator=(const ExprT& r) {
    if (&r != this) { m_rep->dec_ref(); m_rep = r.m_rep; m_rep->inc_ref(); }
    return *this;
  }
public:
  ExprT& operator+=(const ExprT& e)
  { m_rep = new AddRep(m_rep, e.m_rep, true); return *this; }
  ExprT& operator-=(const ExprT& e)
  { m_rep = new SubRep(m_rep, e.m_rep, true); return *this; }
  ExprT& operator*=(const ExprT& e)
  { m_rep = new MulRep(m_rep, e.m_rep, true); return *this; }
  ExprT& operator/=(const ExprT& e)
  { m_rep = new DivRep(m_rep, e.m_rep, true); return *this; }

  ExprT operator+() const 
  { return ExprT(*this); }
  ExprT operator-() const 
  { if(get_rational_reduce_flag()) return getExactVal(- *this);
    return ExprT(new NegRep(m_rep)); }

  ExprT& operator++()
  { *this += 1; return *this; }
  ExprT operator++(int)
  { ExprT r(*this); ++(*this); return r; }
  ExprT& operator--()
  { *this -= 1; return *this; }
  ExprT operator--(int)
  { ExprT r(*this); --(*this); return r; }

  /// addition
  friend ExprT operator+(const ExprT& e1, const ExprT& e2) 
//  { return ExprT(new AddRep(e1.rep(), e2.rep())); }
  { 
    ExprT e(new AddRep(e1.rep(), e2.rep())); 
    if(get_rational_reduce_flag()) return getExactVal(e);
    return e;
  }

  template <typename T>
  friend ExprT operator+(const ExprT& e1, const T& v)
  { return ExprT(new AddRep(e1.rep(), ExprT(v).rep())); }
  template <typename T>
  friend ExprT operator+(const T& v, const ExprT& e2) 
  { return ExprT(new AddRep(ExprT(v).rep(), e2.rep())); }
  
  /// subtraction
  friend ExprT operator-(const ExprT& e1, const ExprT& e2)
  { 
    ExprT e(new SubRep(e1.rep(), e2.rep())); 
    if(get_rational_reduce_flag()) return getExactVal(e);
    return e;
  }
  template <typename T>
  friend ExprT operator-(const ExprT& e1, const T& v)
  { return ExprT(new SubRep(e1.rep(), ExprT(v).rep())); }
  template <typename T>
  friend ExprT operator-(const T& v, const ExprT& e2)
  { return ExprT(new SubRep(ExprT(v).rep(), e2.rep())); }

  /// multiplication
  friend ExprT operator*(const ExprT& e1, const ExprT& e2)
  { 
    ExprT e(new MulRep(e1.rep(), e2.rep())); 
    if(get_rational_reduce_flag()) return getExactVal(e);
    return e;
  }
//{ return ExprT(new MulRep(e1.rep(), e2.rep())); }
  template <typename T>
  friend ExprT operator*(const ExprT& e1, const T& v)
  { return ExprT(new MulRep(e1.rep(), ExprT(v).rep())); }
  template <typename T>
  friend ExprT operator*(const T& v, const ExprT& e2)
  { return ExprT(new MulRep(ExprT(v).rep(), e2.rep())); }

 /// division
  friend ExprT operator/(const ExprT& e1, const ExprT& e2)
  { 
    ExprT e(new DivRep(e1.rep(), e2.rep())); 
    if(get_rational_reduce_flag()) return getExactVal(e);
    return e;
  }
  template <typename T>
  friend ExprT operator/(const ExprT& e1, const T& v)
  { return ExprT(new DivRep(e1.rep(), ExprT(v).rep())); }
  template <typename T>
  friend ExprT operator/(const T& v, const ExprT& e2)
  { return ExprT(new DivRep(ExprT(v).rep(), e2.rep())); }

  /// square root
  friend ExprT sqrt(const ExprT& e)
  { return ExprT(new SqrtRep(e.rep())); }
      // TO DO:
	  // Mar 9, 2009 (Chee/Jihun): in Core2, Level 3, "long" becomes "BigInt".
	  // This caused ambiguity for sqrt(BigInt) in graham.cpp.  The
	  // following templated definition for sqrt(const T& v) removed the
	  // compilation error, but the result is wrong.  Still to be fixed.
  template <typename T>
  friend ExprT sqrt(const T& v)
  { return ExprT(new SqrtRep(ExprT(v).rep())); }

  /// cube root
  friend ExprT cbrt(const ExprT& e)
  { return ExprT(new CbrtRep(e.rep())); }
  friend ExprT dbrt(const ExprT& e)
  { return ExprT(new CbrtRep(e.rep())); }

  /// kth-root
  friend ExprT root(const ExprT& e, unsigned long k)
  { return ExprT(new RootRep(e.rep(), k)); }

  /// pi
  static void pi(ExprT& e) { // Jihun's solution to the friend injection problem for pi
    e = ExprT(new PiRep());
  }
//  friend ExprT pi(const ExprT& e=)
//  { return ExprT(new PiRep()); }

  /// e
  static void e(ExprT& e) // Jihun's solution to the friend injection problem for e()
  { e = ExprT(new ERep()); }

  /// sine
  friend ExprT sin(const ExprT& e)
  { return ExprT(new SinRep(e.rep())); }

  /// power of e
  friend ExprT exp(const ExprT& e)
  { return ExprT(new ExpRep(e.rep())); }

  /// power of 2
  friend ExprT exp2(const ExprT& e)
  { return ExprT(new Exp2Rep(e.rep())); }

  /// power of 10
  friend ExprT exp10(const ExprT& e)
  { return ExprT(new Exp10Rep(e.rep())); }

  /// log
  friend ExprT log(const ExprT& e)
  { return ExprT(new LogRep(e.rep())); }

  /// log2
  friend ExprT log2(const ExprT& e)
  { return ExprT(new Log2Rep(e.rep())); }

  /// log10
  friend ExprT log10(const ExprT& e)
  { return ExprT(new Log10Rep(e.rep())); }

  /// cosine
  friend ExprT cos(const ExprT& e)
  { return ExprT(new CosRep(e.rep())); }

  /// tangent
  friend ExprT tan(const ExprT& e) {
      return sin(e)/cos(e);
  }

  /// cotangent
  friend ExprT cot(const ExprT& e) {
      return cos(e)/sin(e);
  }

  /// arcsine
  friend ExprT asin(const ExprT& e) {
    if (e.abs() > 1) {
      core_error("asin out of range", __FILE__, __LINE__, true);
      return ExprT(0);
    }
    else if (e.abs() >= 1 / sqrt(ExprT(2)))
      return 2 * e.sign() * asin(sqrt((1 - sqrt(1 - e*e)) / 2));
    else if (e.abs() >= 0.5)
      return atan(e / sqrt(1-e*e));
    else
      return ExprT(new ASinRep(e.rep()));
  }

  /// arcsocine
  friend ExprT acos(const ExprT& e) {
    if (e.abs() > 1) {
      core_error("acos out of range", __FILE__, __LINE__, true);
      return ExprT(0);
    }
    else if (e.abs() >= 1 / sqrt(ExprT(2)))
      return 2 * e.sign() * acos(sqrt((1 + sqrt(1 - e*e)) / 2));
    else if (e.abs() >= 0.5)
      return ExprT(new PiRep()) / 2 - asin(e);
    else
      return ExprT(new ACosRep(e.rep()));
  }

  /// arctangent
  friend ExprT atan(const ExprT& e) {
    if (e.abs() >= 1)
      return asin(e / (sqrt(pow(e,2) + 1)));
    else
      return ExprT(new ATanRep(e.rep()));
  }

  // radical -- alternative name for root(n,k)
  template<class NT>
  friend ExprT radical(const NT& n, unsigned long k) {
    assert(n>=0 && k>=1);
    /*  This code is slower because root calls MPFR
     *  while this code uses our own Newton iteration.
    if (n==0 || n == 1 || k ==1) return n;
    Polynomial<NT> Q(k);
    Q.setCoeff(0, -n);
    Q.setCoeff(k, 1);
    return rootOf(Q);
    */
    return root(ExprT(n),k);
  }

  /// helper function for constructing Polynomial node (n-th node)
  template <class NT>
  friend ExprT rootOf(const Polynomial<NT>& p, int n = 0) {
    return ExprT(new ConstPolyRepT<RootBd, Filter, Kernel, NT>(p, n));
  }
  /// helper function for constructing Polynomial node witb BFInterval
  template <class NT>
  friend ExprT rootOf(const Polynomial<NT>& p, const BFInterval& I) {
    return ExprT(new ConstPolyRepT<RootBd, Filter, Kernel, NT>(p, I));
  }
  /// helper function for constructing Polynomial node with pair of BigFloats
  template <class NT, class T>
  friend ExprT rootOf(const Polynomial<NT>& p, const T& x, const T& y) {
    return ExprT(new ConstPolyRepT<RootBd, Filter, Kernel, NT>(p, BFInterval(x, y)));
  }

  friend ExprT power(const ExprT& e, long k) {
    if (k==0)  return 1;
    else if (k==1) return e;
    else {
      bool sign = true;
      if (k < 0) 
      { sign = false; k = -k; }
      std::vector<ExprRep*> c;
      ProdRep* newRep = new ProdRep(c);
      for (long i=0; i < k; i++)
        newRep->insert (e.rep());
      if (sign)
        return ExprT(newRep);
      else
        return ExprT(1) / ExprT(newRep);
    }
  }
  friend ExprT pow(const ExprT& e, long k) {
    return power(e,k);
  }
  
  // NOTE: this summation operator is an anary operator
  // in the sense the number of arguments to be summed ranges
  // from fun(n) for n=start to n=end
  //
  // In the future, we should automate this so that 
  // Core can automatically restructure any subexpression that
  // consists entirely of the addition operator.
  //
  template <typename Function, typename T>
  friend ExprT summation(Function fun, T start, T end) {
    SumOpRep* rep = new SumOpRep;
    for (T it = start; it <= end; ++it)
      rep->insert(fun(it).rep());

    return ExprT(rep);
  }

  // NOTE: this product operator is an anary operator
  // which is the multiplication analog of the summation operator
  //
  template <typename Function, typename T>
  friend ExprT product(Function fun, T start, T end) {
    ProdOpRep* rep = new ProdOpRep;
    for (T it = start; it <= end; ++it)
      rep->insert(fun(it).rep());

    return ExprT(rep);
  }

  /// GetExactVal(e) returns an expression e' whose value equals Value of e
  /// but in case the NUMTYPE of e is <= RATIONAL, then e' is a leaf node.
  /// This function is called when we perform reduceToRational
  friend ExprT getExactVal(const ExprT& e) {
    switch(e.rep()->numType()) {
      case NODE_NT_INTEGER:
	      std::cout << "getExactVal INT Type" << std::endl;
	return ExprT(e.rep()->getZTVal());
        break;
      case NODE_NT_DYADIC:
	      std::cout << "getExactVal DYADIC Type" << std::endl;
	return ExprT(e.rep()->getFTVal());
        break;
      case NODE_NT_RATIONAL:
	      std::cout << "getExactVal RAT Type" << std::endl;
	return ExprT(e.rep()->getQTVal());
        break;
      default:
	return e;
    }
  }
 
  /// compare function
  int cmp(const ExprT& e) const
  { return m_rep == e.m_rep ? 0 : SubRep(m_rep, e.m_rep).get_sign(); }

  friend bool operator==(const ExprT& x, const ExprT& y)
  { return x.cmp(y) == 0; }
  friend bool operator!=(const ExprT& x, const ExprT& y)
  { return x.cmp(y) != 0; }
  friend bool operator>=(const ExprT& x, const ExprT& y)
  { return x.cmp(y) >= 0; }
  friend bool operator<=(const ExprT& x, const ExprT& y)
  { return x.cmp(y) <= 0; }
  friend bool operator<(const ExprT& x, const ExprT& y)
  { return x.cmp(y) < 0; }
  friend bool operator>(const ExprT& x, const ExprT& y)
  { return x.cmp(y) > 0; }

  friend std::istream& operator>>(std::istream& is, ExprT& x)
  { 
    std::string str;
    is >> str;
    x.construct_from_string(str.c_str(), 10, digits2bits(getDefaultInputDigits()));
    return is;	// we do not check if input string is valid...
  }
  
  friend std::ostream& operator<<(std::ostream& os, const ExprT& x) {
    ExprT* p = const_cast<ExprT*>(&x);
    if (p->sign()) os << p->approx(defRelPrec,defAbsPrec);
    else os << "0"; return os;
  }
  std::string toString() {
    ExprT* p = const_cast<ExprT*>(this);
    if (p->sign()) return p->approx2(defRelPrec,defAbsPrec).get_str();
    else return std::string("0");
}
//@}

public: // public methods
  /// return relative approximation
  KT& r_approx(prec_t prec)
  { return m_rep->r_approx(prec); }
  /// return absolute approximation
  KT& a_approx(prec_t prec)
  { return m_rep->a_approx(prec); }
  
  /// return a BigFloat approximation \f$[r, \infty]\f$ or \f$[\infty, a]\f$
  /// Note: Compare to approx2(..) that returns a BigFloat2
  FT approx(prec_t r_prec = defRelPrec, prec_t a_prec = defAbsPrec) {
    if (a_prec == CORE_INFTY)
      return r_approx(r_prec).get_f();
    else if (r_prec == CORE_INFTY)
      return a_approx(a_prec).get_f();
    else 
      return a_approx( (std::min)(a_prec, m_rep->rel2abs(r_prec)) ).get_f();
  }
  /// return a BigFloat2 approximation \f$[r, \infty]\f$ or \f$[\infty, a]\f$
  /// Note: Compare to approx(..) that returns a BigFloat only
  KT& approx2(prec_t r_prec = defRelPrec, prec_t a_prec = defAbsPrec) {
    if (a_prec == CORE_INFTY)
      return r_approx(r_prec);
    else if (r_prec == CORE_INFTY)
      return a_approx(a_prec);
    else 
      return a_approx( (std::min)(a_prec, m_rep->rel2abs(r_prec)) );
  }


  /// return integer value 
  int intValue() const {
    return (int)longValue();
  }
  /// return long value 
  long longValue() const {
    return (long)const_cast<ExprT*>(this)->a_approx(2).get_d();
  }
  /// return float value 
  float floatValue() const {
    return (float)doubleValue();
  }
  /// return double value 
  double doubleValue() const {
    return const_cast<ExprT*>(this)->r_approx(52).get_d();
  }
  /// return BigInt value 
  BigInt BigIntValue() const {
    return const_cast<ExprT*>(this)->a_approx(2).get_z(); 
  }
  /// return BigRat value 
  BigRat BigRatValue(const prec_t r=defRelPrec, const prec_t a=defAbsPrec) const {
    return const_cast<ExprT*>(this)->approx2(r,a).get_q();
  }
  /// return BigFloatValue
  FT BigFloatValue(const prec_t r=defRelPrec, const prec_t a=defAbsPrec) const {
    return const_cast<ExprT*>(this)->approx(r,a);
  }
  /// return BigFloat2Value
  KT BigFloat2Value(const prec_t r=defRelPrec, const prec_t a=defAbsPrec) const {
    return const_cast<ExprT*>(this)->approx2(r,a);
  }
  /// double interval
  void doubleInterval(double& lb, double& ub) {
    ExprT* p = const_cast<ExprT*>(this);
    KT interval = p->r_approx(52);
    lb = interval.getLeft().get_d(BF_RNDD);
    ub = interval.getRight().get_d(BF_RNDU);
  }

  /// return sign (dirty cast)
  sign_t sign() const
  { return const_cast<ExprT*>(this)->m_rep->get_sign(); }
  /// return upper bound of MSB (dirty cast)
  msb_t uMSB() const
  { return const_cast<ExprT*>(this)->m_rep->get_uMSB(); }
  /// return lower bound of MSB (dirty cast)
  msb_t lMSB() const
  { return const_cast<ExprT*>(this)->m_rep->get_lMSB(); }
  /// absolute value
  ExprT abs() const 
  { return sign() >= 0 ? +(*this) : -(*this); }
  

  /// return internal rep
  ExprRep* rep() const
  { return m_rep; }

  /// Debug Help Functions
  void debug(int mode, int level, int depthLimit) {
    std::cout << "-------- Expr debug() -----------" << std::endl;
    std::cout << "rep = " << rep() << std::endl;
    if (mode == LIST_MODE)
      rep()->debugList(level, depthLimit);
    else if (mode == TREE_MODE)
      rep()->debugTree(level, 0, depthLimit);
    else
      core_error("unknown debugging mode", __FILE__, __LINE__, false);
    std::cout << "---- End Expr debug(): " << std::endl;
  }
  
private:
  void construct_from_string(const char* str, int base, prec_t prec) {
    if (strchr(str, '/') != 0)
      m_rep = new ConstQTRep(QT(str), NODE_NT_RATIONAL);
    else if (strchr(str, '.') != 0 && is_infty(prec)) {
      m_rep = new ConstQTRep(stringToBigRat(str), NODE_NT_RATIONAL);
    } else if (is_infty(prec))
      m_rep = new ConstZTRep(ZT(str), NODE_NT_INTEGER);
    else
      m_rep = new ConstFTRep(FT(str, base, prec), NODE_NT_DYADIC);
  }
private:
  ExprRep* m_rep; ///<- internal representation
}; // end if ExprT

CORE_END_NAMESPACE

///////////////////////////////////////////////////////////////////////////
// Definition of Expr
///////////////////////////////////////////////////////////////////////////

#include <CORE/RootBounds.h>
#include <CORE/Filters.h>

CORE_BEGIN_NAMESPACE

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// VERY IMPORTANT: This is where we typedef "Expr" by choosing the
//       RootBound, Filter, and Kernel Modules
// for the templated class ExprT.  
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// YOU CAN CHOOSE the version of Expr you like by uncommenting one line below:
// 	We recommend using CASE 2-a (k-ary BFMSS Bound based on double).
// 	
//=======================================================
// CASE 1-a: BFMSS root bound (based on long) + dummy filter
//
//    typedef ExprT<BfmssRootBd<BigFloat2>, DummyFilter, BigFloat2> Expr;
//
// CASE 1-b: BFMSS root bound (based on long) + BFS filter
//
//    typedef ExprT<BfmssRootBd<BigFloat2>, BfsFilter<BigFloat2>, BigFloat2> Expr;
//
// CASE 2-a: kary-BFMSS root bound (based on double) + BFS filter
//
 typedef ExprT<BfmsskRootBd_double<BigFloat2>, BfsFilter<BigFloat2>, BigFloat2> Expr;
//
// CASE 2-b: kary-BFMSS root bound (based on BigFloat) + BFS filter
//
//   typedef ExprT<BfmsskRootBd_BigFloat<BigFloat2>, BfsFilter<BigFloat2>, BigFloat2> Expr;
//
// CASE 2-c: kary-BFMSS root bound (based on long) + BFS filter
//
//    typedef ExprT<BfmsskRootBd<BigFloat2>, BfsFilter<BigFloat2>, BigFloat2> Expr;
//
// CASE 2-d: kary-BFMSS root bound (based on double) + Dummy filter
//
//    typedef ExprT<BfmsskRootBd_double<BigFloat2>, DummyFilter, BigFloat2> Expr;
//
// CASE 3: Dummy root bound + Dummy filter (NOTE: Apparently not implemented!)
//
//   typedef ExprT<DummyRootBd<10>, DummyFilter, BigFloat2> Expr;
//
//=======================================================
//

/// pi()
inline Expr pi() {
  Expr x; Expr::pi(x); return x; 
}

/// e()
inline Expr e() {
  Expr x; Expr::e(x); return x;
}

/// absolute value
inline Expr abs(const Expr& x) {
  return x.abs();
}
/// absolute value (same as abs)
inline Expr fabs(const Expr& x) {
  return abs(x);
}

/// sign of Expr
inline sign_t sign(const Expr& x) {
  return x.sign();
}

/// Convert to double value
inline double Todouble(const Expr& e, prec_t r = defRelPrec, prec_t a=defAbsPrec) {
  Expr* p = const_cast<Expr*>(&e);
  if (p->sign()) {
    p->approx(r,a);
    return e.doubleValue();
  } else
    return 0;
}

/// convert Expr to BigFloat2
inline BigFloat2 ToBigFloat2(const Expr& e, prec_t r = defRelPrec, prec_t a=defAbsPrec) {
  Expr* p = const_cast<Expr*>(&e);
  if (p->sign()) {
    p->approx(r,a);
    return e.BigFloat2Value();
  } else
    return BigFloat2(0);
}

inline BigInt ToBigInt(const Expr& e, prec_t r = defRelPrec, prec_t a=defAbsPrec) {
  Expr* p = const_cast<Expr*>(&e);
  if (p->sign()) {
    p->approx(r,a);
    return e.BigIntValue();
  } else
    return 0;
}

inline BigRat ToBigRat(const Expr& e, prec_t r = defRelPrec, prec_t a=defAbsPrec) {
  Expr* p = const_cast<Expr*>(&e);
  if (p->sign()) {
    p->approx(r,a);
    return e.BigRatValue();
  } else
    return 0;
}

inline bool isDivisible(const Expr& x, const Expr& y) {
  Expr e = x/y;
  return ((e - ToBigInt(e, CORE_INFTY, 2)) == 0);
}

inline Expr gcd(const Expr& x, const Expr& y) {
  return 1;
}

inline Expr div_exact(const Expr& x, const Expr& y) {
  return (x/y).approx(CORE_INFTY, 2);
}

inline BigInt floor(const Expr& x) {
  BigInt r = ToBigInt(x, CORE_INFTY, 2);
  if (x - r >= 0)
    return r;
  else
    return --r;
}

inline BigInt ceil(const Expr& x) {
  return -floor(-x);
}

inline long floorLg(const Expr& x) {
  if (x < 1)
    return -ceilLg(ceil(1/x));
  else
    return floorLg(floor(x));
}

inline long ceilLg(const Expr& x) {
  if (x < 1)
    return -floorLg(floor(1/x));
  else
    return ceilLg(ceil(x));
}


CORE_END_NAMESPACE
#endif // __CORE_EXPR_H__
