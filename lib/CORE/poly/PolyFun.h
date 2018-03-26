#ifndef __CORE_POLYFUN_H__
#define __CORE_POLYFUN_H__

//#include <CORE/poly/Poly.h>
//#include <CORE/BigFloat.h>
//#include <CORE/Expr.h>
#include <CORE/Promote.h>
#include <CORE/CoreAux.h>

CORE_BEGIN_NAMESPACE
typedef long extLong;

///Various forms of evaluation. Some specific to the type of the point of evaluation.
//{@

/// Evaluation Function (generic version, always returns the exact value).
///
///  This evaluation function is easy to use, but may not be efficient
///  when you have BigRat or Expr values.
///
/// User must be aware that the return type of eval is Max of Types NT and T.
///
/// E.g., If NT is BigRat, and T is Expr then Max(NT,T)=Expr. 
/// 	
/// REMARK: If NT is BigFloat, it is assumed that the BigFloat is error-free.  
template <typename NT, typename T>
MAX_TYPE(NT, T) eval(Polynomial<NT> &p,const T& f) {	// evaluation
  typedef MAX_TYPE(NT, T) ResultT;
  int deg = p.getTrueDegree();
  if (deg == -1)
    return ResultT(0);
  if (deg == 0 || f == 0)
    return ResultT(p.coeff()[0]);
  ResultT val(0);
  ResultT ff(f);
  for (int i=deg; i>=0; i--) {
    val *= ff;
    val += ResultT(p.coeff()[i]);	
  }
  return val;
}//eval

/// Approximate Evaluation of Polynomials
/// 	the coefficients of the polynomial are approximated to some
///	specified composite precision (r,a).
/// @param  f evaluation point 
/// @param  r relative precision to which the coefficients are evaluated
/// @param  a absolute precision to which the coefficients are evaluated
/// @return a BigFloat with error containing value of the polynomial.
///     If zero is in this BigFloat interval, then we don't know the sign.
//
// 	ASSERT: NT = BigRat or Expr
//
template <class NT>
BigFloat2 evalApprox(Polynomial <NT> &p, const BigFloat& f, 
	const extLong& r=defRelPrec, const extLong& a=defAbsPrec) {// evaluation
  int deg = p.getTrueDegree();
  
  if (deg == -1)
    return BigFloat2(0);
  if (deg == 0 || f == 0)
    return ToBigFloat2(p.coeff()[0], r);
    

  BigFloat2 val(0), c, ff(f);
  for (int i=deg; i>=0; i--) {
    c = ToBigFloat2(p.coeff()[i], r);	
    val *= ff; 
    val += c;
  }
  return val;
}//evalApprox

// This BigInt version of evalApprox should never be called...
inline 
BigFloat2 evalApprox( Polynomial<BigInt> &p, const BigFloat& f,
	const extLong& r, const extLong& a) {// evaluation
  assert(0);
  return BigFloat2(0);
}



/**
 * Evaluation at a BigFloat value
 * using "filter" only when NT is BigRat or Expr.
 * Using filter means we approximate the polynomials
 * coefficients using BigFloats.  If this does not give us
 * the correct sign, we will resort to an "exact" evaluation
 * using Expr.
 *
 * If NT <= BigFloat, we just call eval().
 *
   We use the following heuristic estimates of precision for coefficients:

      r = 1 + lg(|P|_\infty) + lg(d+1)  		if f <= 1
      r = 1 + lg(|P|_\infty) + lg(d+1) + d*lg|f| 	if f > 1
      
   if the filter fails, then we use Expr to do evaluation.

   This function is mainly called by Newton iteration (which
   has some estimate for oldMSB from previous iteration).

   @param p polynomial to be evaluated
   @param val the evaluation point
   @param oldMSB an rough estimate of the lg|p(val)|
   @return bigFloat interval contain p(val), with the correct sign

 ***************************************************/
template <class NT>
BigFloat2 evalExactSign(Polynomial<NT> &p, const BigFloat& val,
	 const extLong& oldMSB = 54) {
  if (p.getTrueDegree() == -1)
    return BigFloat2(0);
  
  if (hasExactDivision<NT>::check()) { // currently, only to detect NT=Expr and NT=BigRat
    extLong r;
    r = 1 + height(p).uMSB() + ceilLg(long(p.getTrueDegree()+1));
    if (val > 1)
      r += p.getTrueDegree() * val.uMSB();
    r += (std::max)(extLong(0), oldMSB);

    BigFloat2 rVal = evalApprox(p, val, r);
    if (rVal.isZeroIn()) {
      Expr eVal = eval(p, Expr(val));	// eval gives exact value
      return ToBigFloat2(eVal, 54, CORE_INFTY);
    } else 
      return rVal;
  } else {
    return ToBigFloat2(eval(p, val));
  }
  
  assert(0);
  return BigFloat2(0);
}//evalExactSign
  
/// Cauchy Upper Bound on Roots.
// -- ASSERTION: NT is an integer type
template < class NT >
BigFloat CauchyUpperBound(const Polynomial<NT> &p) {
  if (zeroP(p))
    return BigFloat(0);
  NT mx = 0;
  int deg = p.getTrueDegree();
  for (int i = 0; i < deg; ++i) {
    mx = (std::max)(mx, abs(p.coeff()[i]));
  }
  Expr e = mx;
  e /= Expr(abs(p.coeff()[deg]));
  e.approx(CORE_INFTY, 2);
  // get an absolute approximate value with error < 1/4
  return (e.BigFloatValue() + 2);
}

/// Moebius Transform
template<class NT, class T>
inline Polynomial<MAX_TYPE(NT,T)>
moebiusTransform (const Polynomial<NT>& _poly,
                  const T& a, const T& b,
                  const T& c, const T& d) {
  typedef MAX_TYPE(NT, T) ResultT;
  ResultT coeff[2];
  Polynomial<ResultT> _maxpoly(_poly);

  coeff[0] = a; coeff[1] = b;
  Polynomial<ResultT> bXplusa(1, coeff);

  coeff[0] = c; coeff[1] = d;
  Polynomial<ResultT> dXplusc(1, coeff); 

  //return composeHornerBinary(_poly, bXplusa, dXplusc);
  return composeBinary(_maxpoly, bXplusa, dXplusc);
}

/// Moebius Transform
inline Polynomial<BigFloat>
moebiusTransform (const Polynomial<BigFloat>& _poly,
                  const BigFloat& a, const BigFloat& b,
                  const BigFloat& c, const BigFloat& d) {
  BigFloat coeff[2];

  coeff[0] = a; coeff[1] = b;
  Polynomial<BigFloat> bXplusa(1, coeff);

  coeff[0] = c; coeff[1] = d;
  Polynomial<BigFloat> dXplusc(1, coeff); 

  return composeHornerBinary(_poly, bXplusa, dXplusc);
  //return composeBinary(_poly, bXplusa, dXplusc);
}

//@}
CORE_END_NAMESPACE

#endif /*__CORE_POLYFUN_H__*/
