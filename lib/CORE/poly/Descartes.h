#ifndef CORE_DESCARTES_H
#define CORE_DESCARTES_H

#include <CORE/poly/Poly.h>

CORE_BEGIN_NAMESPACE

/***************************************************
 * Descartes Class:
 *   this is modeled after the Sturm class.
 *
 * Author: Vikram Sharma, Jihun Yu,  (Dec 2006)
 ***************************************************/

template <typename T>
inline void shift(T* coeff, int deg, T* shifted){
  //This is the ascending coefficient method suggested by
  //Krandick in Isolierung reeller Nullstellen von Polynomen
  //( English version is called Isolation of Polynomial Real Roots)
  //to compute the Taylor shift of a polynomial by one.

  for(int i=0; i<= deg ; i++)
    shifted[i] = coeff[i];
  
  for(int i=0; i<= deg-1;i++)
    for(int j=deg-1; j>=i; j--)
      shifted[j]+=shifted[j+1];

}

template <typename T>
inline void half(T* coeff, int deg, T* halved){
  for(int i=0; i<= deg; i++)
    halved[i] = coeff[i]*(BigInt(1)<<(deg-i));
}

//Given the coefficient sequence coeff of some polynomial P
//this computes the coefficient sequence contracted of the polynomial
//P(\lambda X). Since lambda can be of a different type than the
//coefficients of P we need introduce another typename to resolve this.
//We assume that T2 is more general than T1 so that the conversion from
// the latter to the former can take place without error.
template <typename T1, typename T2, typename T3>
inline void contract(T1* coeff, int deg, T2 lambda, T3* contracted){
  T2 pow=1;
  for(int i=0; i <= deg ; i++){
    contracted[i] = coeff[i]*T1(pow);
    pow*=lambda;
  }
}

//void shift(BigFloat* coeff, int deg, BigFloat lambda, BigFloat* shifted){}

//Computes the Taylor shift by a constant lambda. Confer the comments
//for contract above.
template <typename T1, typename T2, typename T3>
inline void shift(T1* coeff, int deg, T2 lambda, T3* shifted){

  for(int i=0; i<= deg ; i++)
    shifted[i] = coeff[i];
  
  if(lambda != 0){
    for(int i=0; i<= deg-1;i++)
      for(int j=deg-1; j>=i; j--)
	shifted[j]+= T3(lambda) * shifted[j+1];
  }
}

template <typename T>
inline int shiftAndSigncount(T* coeff, int deg){
  //This is the ascending coefficient method suggested by
  //Krandick in Isolierung reeller Nullstellen von Polynomen
  //( English version is called Isolation of Polynomial Real Roots)
  //to compute the Taylor shift of a polynomial by one.
  //The advantage of this method is that it computes the coefficient of x^i
  //in n-i steps. Thus we can check for sign variation as we compute
  //the coefficients.

  //First reverse the polynomial
  T* temp = new T[deg + 1];
  for(int i=0; i<= deg ; i++)
    temp[i] = coeff[deg-i];

  
  int num=0, i;
  int lastsign =0;//The last non-zero sign
  int currsign;//The sign of the current coefficient

  //Compute the sequence of coefficients and simultaneously get
  //the number of sign variations. If the sign variations are greater
  //than one then break out of the loop, i.e. return num.
  for(i=0; i<= deg-1;i++){
    for(int j=deg-1; j>=i; j--)
      temp[j]+=temp[j+1];

    currsign = sign(temp[i]);
    if(currsign !=0){
      if (lastsign * currsign < 0) num++;
      lastsign = currsign;//lastsign is always non-zero except for the starting.
    }
     if(num > 1) return -1;
  }
  //To account the sign variation between the previous non-zero coefficient
  // and deg coeffecient.
  //This is done if the number of sign variations were less than or
  //equal to one from the above loop.
  if(lastsign*sign(temp[deg]) < 0) num++;
  //  if(sign(temp[deg-1])*sign(temp[deg]) < 0) num++;

  delete[] temp;

  if(num ==2 ) return -1;
  return num;
}
//Isolates real roots of P in the interval (0,1).
//The roots of P are in bijective correspondence with the roots of
//the original input polynomial P_{in}. More precisely,
//P has a root in (0,1) iff P_{in} has a root in I.
//We also ensure that the list of isolating intervals in v
//is sorted. The intervals in v are of two types:
//1) Both the end points of the interval are the same. This happens
//only if the end point is a root of P.
//2) Otherwise the interval represents an open interval which contains a root.
// Note, the end points may be roots but they will be distinct in this case.
//Advantage of this representation is that we don't have to compute the separation
//bound.
template <typename T>
inline void isolateRoots_unit(Polynomial<T>& P, const BFInterval I, int deg,
                    BFVecInterval &v) {
  
 int num = shiftAndSigncount(P.coeff(), deg);
  //std::cout << "sign variations after shift " << num << std::endl;
  //std::cout <<"Interval is ["<<I.first << ":"<< I.second << "]"<<std::endl;
  //std::cout << "Polynomial is "; P.dump() ; std::cout <<std::endl;
  if(num == 0) return;
  else if(num == 1)
    if (I.first > I.second)
      v.push_back(BFInterval(I.second, I.first));
    else
      v.push_back(I);
  else{
    BigFloat m = div2(I.second + I.first);
    T* temp1 = new T[deg +1];
    T* temp2 = new T[deg +1];

    half(P.coeff(),deg, temp1);
    Polynomial<T> Q(deg, temp1);
    //    std::cout <<"After halving polynomial is "; Q.dump(); std::cout<<std::endl;
    
    shift(temp1, deg, temp2);
    Polynomial<T> R(deg, temp2);
    //    std::cout<<"Inside isolateRoots: second polynomial "<< R << std::endl;
    
    BFInterval J = std::make_pair(I.first, m);
    BFInterval JJ = std::make_pair(m, I.second);
    isolateRoots_unit(Q, J, deg, v);
    if(R.coeff()[0] == 0) 
      v.push_back(std::make_pair(m,m));
    isolateRoots_unit(R, JJ, deg, v);
  }
}

inline void isolateRoots(Polynomial<BigInt>& P, BFInterval I, BFVecInterval& v){
  int deg = P.getTrueDegree();
  if(deg == 0)
    std::cout<< "Polynomial is a constant" << std::endl;
  
  BigFloat a = I.first, b=I.second;
  
  int n = P.getTrueDegree();

  BigFloat* temp1 = new BigFloat[n+1];
  BigFloat* temp2 = new BigFloat[n+1];

  shift(P.coeff(), n, a, temp1);
  contract(temp1, n, b-a, temp2);

  if(P.evalSign(a) == 0)
    v.push_back(std::make_pair(a, a));

  Polynomial<BigFloat> R(n, temp2);
  //std::cout <<"Corresponding polynomial with roots in unit interval ";R.dump();
  //std::cout <<std::endl;

  delete[] temp1;
  delete[] temp2;
  isolateRoots_unit(R, I, n, v);

 if(P.evalSign(b) == 0)
    v.push_back(std::make_pair(b, b));
}

template <typename T>
inline void isolateRoots(Polynomial<T>& P, BFInterval I, BFVecInterval& v){
  int deg = P.getTrueDegree();
  if(deg == 0)
    std::cout<< "Polynomial is a constant" << std::endl;
  
  BigFloat a = I.first, b=I.second;

  int n = P.getTrueDegree();
  T temp1[n + 1], temp2[n+1];

  shift(P.coeff(), n, a, temp1);
  contract(temp1, n, b-a, temp2);

  if(evalExactSign(P, a).sgn() == 0)
    v.push_back(std::make_pair(a, a));

  Polynomial<T> R(n, temp2);

  isolateRoots_unit(R, I, n, v);

  if(evalExactSign(P, b).sgn() == 0)
    v.push_back(std::make_pair(b, b));
}

inline void isolateRoots(Polynomial<int>& P, BFInterval I, BFVecInterval& v)
{ Polynomial<BigInt> poly(P); isolateRoots(poly, I, v); }

inline void isolateRoots(Polynomial<long>& P, BFInterval I, BFVecInterval& v)
{ Polynomial<BigInt> poly(P); isolateRoots(poly, I, v); }

//An alternative approach to isolating all roots. Here we separate the positive
//and negative roots separately. This is slightly more efficient than the method
//above.
template <typename T>
inline void isolateRoots(Polynomial<T>& P, BFVecInterval& v)
{
  int deg = P.getTrueDegree();
  if(deg == 0)
    std::cout<< "Polynomial is a constant" << std::endl;

  //Compute an upper bound on the positive roots of P
  BigInt B = CauchyBound(P);

  bool isZeroRoot = false;
  
  if(P.coeff()[0] == 0)
    isZeroRoot = true;

  T* temp1 = new T[deg+1];

  //Construct the polynomial whose roots in the unit interval correspond
  //with the roots of P in (0, B)
  contract(P.coeff(), deg, B, temp1);
  Polynomial<T> Q(deg, temp1);
  
  BFVecInterval vPos;//Stores the intervals for the positive roots of P
  BFInterval I = BFInterval(0,B);
  isolateRoots_unit(Q, I, deg, vPos);

  //Flip the signs of the coefficients of Q to construct a polynomial whose
  //roots in the unit interval correspond with the roots of P in (-B, 0). 
  for(int i=1; i<= deg; i++){
    if(i % 2 != 0)
      Q.coeff()[i] *=-1;
  }
  
  I = BFInterval(0, -B);
  isolateRoots_unit(Q, I, deg, v);
  //This ensures that the interval corresponding to the roots are in sorted order 
  if(isZeroRoot)
    v.push_back(std::make_pair(0,0));
  
  for (BFVecInterval::const_iterator it = vPos.begin(); it != vPos.end(); ++it)
    v.push_back(*it);

    delete[] temp1;
}

inline void isolateRoots(Polynomial<long>& P, BFVecInterval& v)
{ Polynomial<BigInt> poly(P); isolateRoots(poly, v); }

inline void isolateRoots(Polynomial<int>& P, BFVecInterval& v)
{ Polynomial<BigInt> poly(P); isolateRoots(poly, v); }

template <class NT>
class Descartes {
private:
  int len;      // len is 1 less than number of non-zero entries in array seq.
  		//     I.e., len + 1 = length of the Sturm Sequence
                // N.B. When len = -1 or len = 0 are special,
                //     the array seq is not used!
                //     Hence, one must test these special cases
  static int N_STOP_ITER;    // Stop IterE after this many iterations. This
                             // is initialized below, outside the Newton class
  Polynomial<NT> _poly;
  Polynomial<NT> _poly_derivative;
  bool NEWTON_DIV_BY_ZERO;  // this is set to true when divide by 0 in Newton
public:
  Descartes() : NEWTON_DIV_BY_ZERO(false) {}
  Descartes(Polynomial<NT> p) : _poly(p), NEWTON_DIV_BY_ZERO(false){
    len = p.getTrueDegree();
    if (len < 0) return;
    
    _poly.sqFreePart();
    _poly_derivative = differentiate(_poly);  
  }
  Descartes(Descartes& d) : _poly(d._poly), NEWTON_DIV_BY_ZERO(false),
  	len(d.len){}

  void setPoly(Polynomial<NT> p) {
    len = p.getTrueDegree();
    if (len < 0) return;

    _poly = p;
    _poly.sqFreePart();
    _poly_derivative = differentiate(_poly);  
  }

  //Isolates roots of P in the CLSOED interval I; assumes P is square-free.
  //This is achieved by computing a polynomial Q whose roots in the
  //unit interval are in bijective correspondence with the roots of
  //P in I. More precisely, Q(X) = P((b-a)X + a) where I=(a,b).
  //This is obtained by first computing the Taylor shift by a of the polynomial
  //and then doing a contraction by b-a.
  //Again v contains the sorted list of intervals.
  void isolateRoots(const BFInterval &I, BFVecInterval& v)
  { CORE_NS::isolateRoots(_poly, I, v); }

  void isolateRoots(const BigFloat &x, const BigFloat &y,
                    BFVecInterval &v) {
    isolateRoots(BFInterval(x, y), v);
   }

  // isolateRoots(v)
  ///   isolates all roots and returns them in v
  /**   v is a vector of isolated intervals
   */
  void isolateRoots(BFVecInterval& v) {
    CORE_NS::isolateRoots(_poly, v);
  }
  // isolateRoot(i)
  ///   Isolates the i-th smallest root 
  ///         If i<0, isolate the (-i)-th largest root
  ///   Defaults to i=0 (i.e., the smallest positive root a.k.a. main root)
  BFInterval isolateRoot(int i = 0) {
    if (len <= 0) 
       return BFInterval(1,0);   // ERROR CONDITION
    if (i == 0)
      return mainRoot();
    BigFloat bd = CauchyUpperBound(_poly);
    return isolateRoot(i, -bd, bd);
  }
  // isolateRoot(i, x, y)
  ///   isolates the i-th smallest root in [x,y]
  /**   If i is negative, then we want the i-th largest root in [x,y]
   *    We assume i is not zero.
   */
  BFInterval isolateRoot(int i, BFInterval &I)
  { return isolateRoot(i, I.first, I.second); }

  BFInterval isolateRoot(int i, const BigFloat& x, const BigFloat& y) {
    BFVecInterval v;

    // We isolate the positive and negative roots separately, thus ensuring
    // that zero is not contained within an interval.
    if(sign(x) == sign(y))
      isolateRoots(BFInterval(x,y), v);
    else{
      isolateRoots(BFInterval(x, 0), v);
      if(_poly.coeff()[0] == NT(0)) // zero is a root of P
        v.erase(v.end()-BFVecInterval::iterator::difference_type(1), v.end()); // erase the entry corresponding to zero in
                                    // v since the next call we add it again
      isolateRoots(BFInterval(0, y), v);
    }

  
    int n= v.size(); //the precise number of real roots in I
    if (i < 0) {//then we want the n-i+1 root
      i += n+1;
      if (i <= 0)
        return BFInterval(1,0); // ERROR CONDITION
    }
    if(i > n)
      return BFInterval(1,0);  // ERROR CONDITION INDICATED

    //Now 0< i <= n
    if (n == 1) {// Thus there is only one root in (x,y). Moreover, the way
                 // we isolated the root we are sure that the interval in v
                 // does not contain any zero. 
      return *(v.begin());
    }

    //Otherwise traverse v and return the i'th interval
    if(i == n)//slight optimization
      return (*(v.end()-BFVecInterval::iterator::difference_type(1)));

    //Now 1 <= i < n
    int count =1;
    for (BFVecInterval::const_iterator it = v.begin(); ; ++it) {
      if(count == i)
        return (*it);
      else
        count ++;
    }
  }

  // same as isolateRoot(i).
  BFInterval diamond(int i) {
    return isolateRoot(i);
  }

  // First root above
  BFInterval firstRootAbove(const BigFloat &e) {
    if (len <= 0)
       return BFInterval(1,0);   // ERROR CONDITION
    return isolateRoot(1, e, CauchyUpperBound(_poly));
  }

  // Main root (i.e., first root above 0)
  BFInterval mainRoot() {
    if (len <= 0)
       return BFInterval(1,0);   // ERROR CONDITION
    return isolateRoot(1, 0, CauchyUpperBound(_poly));
  }

  // First root below
  BFInterval firstRootBelow(const BigFloat &e) {
    if (len <= 0)
       return BFInterval(1,0);   // ERROR CONDITION
    BigFloat bd = CauchyUpperBound(_poly); // bd is exact
    int n = numberOfRoots(-bd, e);
    if (n <= 0)
      return BFInterval(1,0);
    //BigFloat bdBF = BigFloat(ceil(bd));
    BigFloat bdBF;
    bdBF.ceil(bd);
    if (n == 1)
      return BFInterval(-bdBF, e);
    return isolateRoot(n, -bdBF, e);
  }
  
  int signVar(BFInterval& I)
  { return signVar(I.first, I.second); }
  
  int signVar(const BigFloat &x, const BigFloat &y) {
    return signVariationofCoeff(
           moebiusTransform(_poly, x, y, BigFloat(1), BigFloat(1))
	   );
  }
  int numberOfRoots() {
    BFVecInterval v;
    isolateRoots(v);
    return v.size();
  }

  int numberOfRoots(BFInterval& I)
  { return numberOfRoots(I.first, I.second); }

  int numberOfRoots(const BigFloat &x, const BigFloat &y) {
    BFVecInterval v;
    isolateRoots(x, y, v);
    return v.size();
  }

  BFInterval refine(const BFInterval& I, int aprec) {
    return refine(I.first, I.second, aprec);
  }

  BFInterval refine(const BigFloat &x, const BigFloat &y, int aprec) {
    assert(x<=y);
    BFInterval retI(std::make_pair(x, y));
    BigFloat eps = BigFloat::exp2(-aprec);
    BigFloat mid;

    sign_t x_sign = sign(evalExactSign(_poly, retI.first));

    while (retI.second - retI.first > eps) {
      mid = div2(retI.second + retI.first);
      sign_t mid_sign = sign(evalExactSign(_poly,mid));
      if (mid_sign == 0) {
        retI.first = retI.second = mid;
        return retI;
      }
      if (x_sign * mid_sign < 0) {
        retI.second = mid;
      } else {
        retI.first = mid;
        x_sign = mid_sign;
      }
    }
    return retI;
  }

  void refineAllRoots(BFVecInterval &v, int aprec) {

    BFVecInterval v1;
    BFInterval J;

    if (v.empty())
      isolateRoots(v);

    for (BFVecInterval::iterator it = v.begin();
         it != v.end(); ++it) {        // Iterate through all the intervals
      //refine them to the given precision aprec
      J = refine(*it, aprec);
      if (NEWTON_DIV_BY_ZERO) {
        J.first = 1;
        J.second = 0;   // indicating divide by zero
      }
      v1.push_back(std::make_pair(J.first, J.second));
    }
    v.swap(v1);
  }//refineAllRoots

  // This is the new version of "refineAllRoots"
  //    	based on Newton iteration
  // It should be used instead of refineAllRoots!
  void newtonRefineAllRoots( BFVecInterval &v, int aprec) {

    BFVecInterval v1;
    BFInterval  J;

    if (v.empty())
      isolateRoots(v);

    for (BFVecInterval::iterator it = v.begin();
         it != v.end(); ++it) {        // Iterate through all the intervals
      //refine them to the given precision aprec
      J = newtonRefine(*it, aprec);
      if (NEWTON_DIV_BY_ZERO) {
        J.first = 1;
        J.second = 0;   // indicating divide by zero
      }
      v1.push_back(std::make_pair(J.first, J.second));
    }
    v.swap(v1);
  }//End of newtonRefineAllRoots

  /** val = newtonIterN(n, bf, del, err, fuMSB, ffuMSB)
   * 
   *    val is the root after n iterations of Newton
   *       starting from initial value of bf and is exact.
   *    fuMSB and ffuMSB are precision parameters for the approximating
   *		the coefficients of the underlyinbg polynomial, f(x).
   *    	THEY are used ONLY if the coefficients of the polynomial
   *		comes from a field (in particular, Expr or BigRat).
   *		We initially approximate the coefficients of f(x) to fuMSB 
   *		relative bits, and f'(x) to ffuMSB relative bits.
   *		The returned values of fuMSB and ffuMSB are the final
   *		precision used by the polynomial evaluation algorithm.
   *    Return by reference, "del" (difference between returned val and value
   *       in the previous Newton iteration)
   *
   *    Also, "err" is returned by reference and bounds the error in "del".
   *
   *    IMPORTANT: we assume that when x is an exact BigFloat,
   *    then Polynomial<NT>::eval(x) will be exact!
   *    But current implementation of eval() requires NT <= BigFloat.
   * ****************************************************/    

  BigFloat newtonIterN(long n, const BigFloat& bf, BigFloat& del,
	BigFloat& err, extLong& fuMSB, extLong& ffuMSB) {
    if (len <= 0) return bf;   // Nothing to do!  User must
                               // check this possibility!
    BigFloat val = bf;  

    // newton iteration
    for (int i=0; i<n; i++) {
      ////////////////////////////////////////////////////
      // Filtered Eval
      ////////////////////////////////////////////////////
      BigFloat2 ff = evalExactSign(_poly_derivative,val, 3*ffuMSB); //3 is a slight hack
      ffuMSB = ff.uMSB();
      //ff is guaranteed to have the correct sign as the exact evaluation.
      ////////////////////////////////////////////////////

      if (ff.sgn() == 0) {
        NEWTON_DIV_BY_ZERO = true;
        del = 0;
        core_error("Zero divisor in Newton Iteration",
                __FILE__, __LINE__, false);
        return 0;
      }

      ////////////////////////////////////////////////////
      // Filtered Eval
      ////////////////////////////////////////////////////
      BigFloat2 f= evalExactSign(_poly, val, 3*fuMSB); //3 is a slight hack
      fuMSB = f.uMSB();
      ////////////////////////////////////////////////////

      if (f.sgn() == 0) {
        NEWTON_DIV_BY_ZERO = false;
        del = 0;    // Indicates that we have reached the exact root
		    //    This is because eval(val) is exact!!!
        return val; // val is the exact root, before the last iteration
      }
      del = (f/ff).getLeft(); // But the accuracy of "f/ff" must be controllable
		    // by the caller...
      err = BigFloat(del,getDefaultBFdivPrec());
      val -= del;
    }
    return val;
  }//newtonIterN

  //Another version of newtonIterN which does not return the error 
  //and passing the uMSB as arguments; it is easier for the user to call
  //this.
  BigFloat newtonIterN(long n, const BigFloat& bf, BigFloat& del){
    BigFloat err;
    extLong fuMSB=0, ffuMSB=0;
    return newtonIterN(n, bf, del, err, fuMSB, ffuMSB);
  }

  // v = newtonIterE(prec, bf, del, fuMSB, ffuMSB)
  //
  //    return the value v which is obtained by Newton iteration
  //    until del.uMSB < -prec, starting from initial value of bf.
  //    Returned value is an exact BigFloat.
  //    We guarantee at least one Newton step (so del is defined).
  //
  //	   The parameters fuMSB and ffuMSB are precision parameters for
  //	   evaluating coefficients of f(x) and f'(x), used similarly
  //	   as described above for newtonIterN(....)
  //
  //    Return by reference "del" (difference between returned val and value
  //       in the previous Newton iteration).  This "del" is an upper bound
  //       on the last (f/f')-value in Newton iteration.
  //
  //    IN particular, if v is in the Newton zone of a root z^*, then z^* is
  //       guaranteed to lie inside [v-del, v+del].
  //
  //    Note that this is dangerous unless you know that bf is already
  //       in the Newton zone.  So we use the global N_STOP_ITER to
  //       prevent infinite loop.

  BigFloat newtonIterE(int prec, const BigFloat& bf, BigFloat& del, 
	extLong& fuMSB, extLong& ffuMSB) {
    // usually, prec is positive
    int count = N_STOP_ITER; // upper bound on number of iterations
    int stepsize = 1;
    BigFloat val = bf;
    BigFloat err = 0;

    do {
      val = newtonIterN(stepsize, val, del, err, fuMSB, ffuMSB);
      count -= stepsize;
      stepsize++; // heuristic
    } while ((del != 0) && ((del.uMSB() >= -prec) && (count >0))) ;

    if (count == 0) core_error("newtonIterE: reached count=0",
		    	__FILE__, __LINE__, true);
    //del = BigFloat(core_abs(del.m()), err, del.exp() );
    //del.makeCeilExact();
    del += err;
    return val;
  }

  //Another version of newtonIterE which avoids passing the uMSB's.
  BigFloat newtonIterE(int prec, const BigFloat& bf, BigFloat& del){
    extLong fuMSB=0, ffuMSB=0;
    return newtonIterE(prec, bf, del, fuMSB, ffuMSB);
  }


  //newtonRefine(J, a) 
  //
  //    ASSERT(J is an isolating interval for some root x^*)
  //
  //    ASSERT(J.first and J.second are exact BigFloats)
  //
  //    Otherwise, the boundaries of the interval are not well defined.
  //    We will return a refined interval with exact endpoints,
  //    still called J, containing x^* and
  //
  // 			|J| < 2^{-a}.
  //
  // 	TO DO: write a version of newtonRefine(J, a, sign) where
  // 	sign=J.first.sign(), since you may already know the sign
  // 	of J.first.  This will skip the preliminary stuff in the
  // 	current version.
  //
  BFInterval newtonRefine(BFInterval &J, int aprec) {

#ifdef CORE_DEBUG
    if (sign(evalExactSign(_poly,J.first)) * sign(evalExactSign(_poly,J.second)) > 0){
      std::cout <<" ERROR! Root is not in the Input Interval " << std::endl;
      std::cout <<" Polynomial is " << _poly << std::endl;
      std::cout <<" Interval is [" << J.first << ", " << J.second << "]" << std::endl;
      std::cout <<" Sign evaluation of a lower bound on poly [" << evalExactSign(_poly,J.first).getLeft()
                << ", " << evalExactSign(_poly,J.first).getRight() << "]" << std::endl;
      std::cout <<" Sign evaluation of a upper bound on poly [" << evalExactSign(_poly,J.second).getLeft()
                << ", " << evalExactSign(_poly,J.second).getRight() << "]" << std::endl;
      return BFInterval(1,0);
    }
#endif
 
#ifdef CORE_DEBUG_NEWTON
    std::cout << "In newtonRefine, input J=[" << J.first
	<< ", " << J.second << "] precision = " << aprec << std::endl;
#endif

    if (len <= 0) return J;   // Nothing to do!  User must
                               // check this possibility!
      

    if((J.second - J.first).uMSB() < -aprec){
      return (J);
    }
    int xSign, leftSign, rightSign;

    leftSign = sign(evalExactSign(_poly, J.first));
    if (leftSign == 0) {
      J.second = J.first;
      return J;
    }

    rightSign = sign(evalExactSign(_poly, J.second));
    if (rightSign == 0) {
      J.first = J.second;
      return J;
    }

    assert( leftSign * rightSign < 0 );

    //N is number of times Newton is called without checking
    // whether the result is still in the interval or not
    #define NO_STEPS 2
    // REMARK: NO_STEPS=1 is incorrect, as it may lead to
    //      linear convergence (it is somewhat similar to Dekker-Brent's
    //      idea of guaranteeing that bisection does not
    //	    destroy the superlinear convergence of Newton.
    int N = NO_STEPS;

    BigFloat x, del, olddel, temp;
    BigFloat err;
    BigFloat yap = yapsBound(_poly);

    BigFloat old_width = J.second - J.first;
    x = div2(J.second + J.first);

    // initial estimate for the evaluation of filter to floating point precision
    extLong fuMSB=54, ffuMSB=54;

    //MAIN WHILE LOOP. We ensure that J always contains the root

    while ( !smaleBoundTest(_poly, _poly_derivative, x) && 
	    (J.second - J.first) > yap &&
	   (J.second - J.first).uMSB() >= -aprec) {

     x = newtonIterN(N, x, del, err, fuMSB, ffuMSB);

      if ((del == 0)&&(NEWTON_DIV_BY_ZERO == false)) {  // reached exact root!
        J.first = J.second = x;
        return J;
      }

      BigFloat left(x), right(x);
      if (del>0) {
      	left -= del; right += del;
      } else {
      	left += del; right -= del;
      }

      // update interval
      if ((left > J.first)&&(left <J.second)) {
	  int lSign = sign(evalExactSign(_poly, left));
          if (lSign == leftSign)  // leftSign=sign of J.first
            J.first = left;
	  else if (lSign == 0) {
            J.first = J.second = left;
            return J;
          } else {
	    J.second = left;
          }	
      }
      if ((right < J.second)&&(right >J.first)) {
	  int rSign = sign(evalExactSign(_poly, right));
          if (rSign == rightSign)
            J.second = right;
	  else if (rSign == 0) {
            J.first = J.second = right;
            return J;
          } else {
            J.first = right;
          }
      }
      BigFloat width = J.second - J.first;

      //left and right are exact, since x is exact.
      if (width*2 <= old_width && !NEWTON_DIV_BY_ZERO) {
                                  // we can get a better root:

	// No, it is not necessary to update x to
	// the midpoint of the new interval J.
	// REASON: basically, it is hard to be smarter than Newton's method!
	// Newton might bring x very close to one endpoint, but it can be
	// because the root is near there!  In any case,
	// by setting x to the center of J, you only gain at most
	// one bit of accuracy, but you stand to loose an
	// arbitrary amount of bits of accuracy if you are unlucky!
	// So I will comment out the next line.  --Chee (Aug 9, 2004).
	// 
	// x = (J.second + J.first).div2();
	if (J.first > x || J.second < x)
	  x = div2(J.second + J.first);

	old_width = width; // update width

        N ++;      // be more confident or aggressive
	           //  (perhaps we should double N)
		   //
      } else {// Either NEWTON_DIV_BY_ZERO=true
	      // Or width has not decreased sufficiently
	x = div2(J.second + J.first);//Reset x to midpoint since it was the
	                                //value from a failed Newton step
	xSign = sign(evalExactSign(_poly, x));
	if (xSign == rightSign) {
	  J.second = x;
	} else if (xSign == leftSign) {
	  J.first = x;
	} else { // xSign must be 0
	  J.first = J.second = x; return J;
	}
	x = div2(J.second + J.first);

	old_width.div2(); // update width
	
	// reduce value of N:
	N = (std::max)(N-1, NO_STEPS);   // N must be at least NO_STEPS
      }
    }//MAIN WHILE LOOP

   if((J.second - J.first).uMSB() >= -aprec){ // The interval J
	            //still hasn't reached the required precision.
	            //But we know the current value of x (call it x_0)
		    //is in the strong Newton basin of the
		    //root x^* (because it passes Smale's bound)
      //////////////////////////////////////////////////////////////////
      //Both x_0 and the root x^* are in the interval J.
      //Let NB(x^*) be the strong Newton basin of x^*.  By definition,
      //x_0 is in NB(x^*) means that:
      //
      //    x_0 is in NB(x^*) iff |x_n-x^*| \le 2^{-2^{n}+1} |x_0-x^*|
      //    
      // where x_n is the n-th iterate of Newton.  
      //    
      //  LEMMA 1: if x_0  \in NB(x^*) then 
      //               |x_0 - x^*| <= 2|del|      (*)
      //  and
      //               |x_1 - x^*| <= |del|       (**)
      //
      //  where del = -f(x_0)/f'(x_0) and x_1 = x_0 + del
      //Proof:
      //Since x_0 is in the strong Newton basin, we have
      //         |x_1-x^*| <= |x_0-x^*|/2.        (***)
      //The bound (*) is equivalent to
      //         |x_0-x^*|/2 <= |del|.
      //This is equivalent to
      //         |x_0-x^*| - |del| <= |x_0-x^*|/2,
      //which follows from
      //         |x_0-x^* + del| <= |x_0-x^*|/2,
      //which is equivalent to (***).  
      //The bound (**) follows from (*) and (***).
      //QED
      //
      //  COMMENT: the above derivation refers to the exact values,
      //  but what we compute is X_1 where X_1 is an approximation to
      //  x_1.  However, let us write X_1 = x_0 - DEL, where DEL is
      //  an approximation to del.  
      //
      //  LEMMA 2:  If |DEL| >= |del|,
      //  then (**) holds with X_1 and DEL in place of x_1 and del.
      //  
      //  NOTE: We implemented this DEL in newtonIterE.   

#ifdef CORE_DEBUG
      std::cout << "Inside Newton Refine: Refining Part " << std::endl;

      if((J.second - J.first) > yap)
	std::cout << "Smales Bound satisfied " << std::endl;
      else
	std::cout << "Chees Bound satisfied " << std::endl;
#endif
      xSign = sign(evalExactSign(_poly, x));
      if(xSign == 0){
	J.first = J.second = x; 
	return J; // missing before!
      }

      //int k = clLg((-(J.second - J.first).lMSB() + aprec).asLong());
      x = newtonIterE(aprec, x, del, fuMSB, ffuMSB);
      xSign = sign(evalExactSign(_poly, x));

      if(xSign == leftSign){//Root is greater than x
	J.first = x;
	J.second = x + abs(del);  // justified by Lemma 2 above
      }else if(xSign == rightSign){//Root is less than x
	J.first = x - abs(del);   // justified by Lemma 2 above
	J.second = x;
      }else{//x is the root
	J.first = J.second = x;
      }
    }

#ifdef CORE_DEBUG
    if (sign(evalExactSign(_poly,J.first)) * sign(evalExactSign(_poly,J.second)) > 0){
      std::cout <<" ERROR! Root is not in the Output Interval " << std::endl;
      std::cout <<" Polynomial is " << _poly << std::endl;
      std::cout <<" Interval is [" << J.first << ", " << J.second << "]" << std::endl;
      std::cout <<" Sign evaluation of a lower bound on poly [" << evalExactSign(_poly,J.first).getLeft()
                << ", " << evalExactSign(_poly,J.first).getRight() << "]" << std::endl;
      std::cout <<" Sign evaluation of a upper bound on poly [" << evalExactSign(_poly,J.second).getLeft()
                << ", " << evalExactSign(_poly,J.second).getRight() << "]" << std::endl;
      return BFInterval(1,0);
    }
    if(J.second - J.first >  BigFloat::exp2(-aprec)) {
      std::cout << "ERROR! Newton Refine failed to achieve desired precision" << std::endl;
      return BFInterval(1,0);
    }
#endif

      return(J);
 }//End of newton refine
}; //class Descartes

template <class NT>
int Descartes<NT>::N_STOP_ITER = 10000;


// testNewtonDescartes( Poly, aprec, n)
//   will run the Newton-Descartes refinement to isolate the roots of Poly
//         until absolute precision aprec.
//   n is the predicated number of roots
//      (will print an error message if n is wrong)
template<class NT>
inline void testNewtonDescartes(const Polynomial<NT>&P, int prec, int n = -1) {
  Descartes<NT> D (P);
  BFVecInterval v;
  D.newtonRefineAllRoots(v, prec);
  std::cout << "   Number of roots is " << v.size();
  if ((n >= 0) & (v.size() == (unsigned)n))
    std::cout << " (CORRECT!)" << std::endl;
  else
    std::cout << " (ERROR!) " << std::endl;

  int i = 0;
  for (BFVecInterval::iterator it = v.begin();
       it != v.end(); ++it) {
    std::cout << ++i << "th Root is in ["
    << it->first << " ; " << it->second << "]" << std::endl;
    if(it->second - it->first <= (BigFloat(1)/power(BigFloat(2), prec)))
      std::cout << " (CORRECT!) Precision attained" << std::endl;
    else
      std::cout << " (ERROR!) Precision not attained" << std::endl;
  }
}// testNewtonDescartes

template<class NT>
inline int signVar(Polynomial<NT> poly, BFInterval I) {
  return signVariationofCoeff(
         moebiusTransform(poly, I.first, I.second, BigFloat(1), BigFloat(1))
  );
}

CORE_END_NAMESPACE

#endif
