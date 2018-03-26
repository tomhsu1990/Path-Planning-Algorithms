#ifndef __CORE_POLY_H__
#define __CORE_POLY_H__

#include <CORE/Expr.h>
#include <CORE/poly/PolyBase.h>
#include <string>
#include <iostream>

// NOTE: This IntervalT.h is  used in evalH and evalCF method(see below)
//	Vikram, Sept. 2012
#include "CORE/IntervalT.h"

CORE_BEGIN_NAMESPACE

/// \class Polynomial Poly.h
/// \brief A template polynomial class
/**
*      REPRESENTATION:
*      --Each polynomial has a nominal "degree" (this
*              is an upper bound on the true degree, which
*              is determined by the first non-zero coefficient).
*      --coefficients are parametrized by some number type "NT".
*      --coefficients are stored in the "coeff" array of
*              length "degree + 1".  
*              CONVENTION: coeff[i] is the coefficient of X^i.  So, a
*                          coefficient list begins with the constant term.
*      --IMPORTANT CONVENTION:
*              the zero polynomial has degree -1
*              while nonzero constant polynomials have degree 0.
* 
*      FUNCTIONALITY:
*      --Polynomial Ring Operations (+,-,*)
*      --Power
*      --Evaluation
*      --Differentiation
*      --Remainder, Quotient 
*      --GCD
*      --Polynomial Composition (done!  Philip Davidson, Dec 2006)
*      --Resultant, Discriminant (planned)
*      --file I/O (planned)
*/

template <class NT>
inline void strToNT(const char* input, NT& output) {
  output = NT(input);
}

template<>
inline void strToNT(const char* input, int& output) {
  output = atoi(input);
}

template<>
inline void strToNT(const char* input, long& output) {
  output = atol(input);
}

template<>
inline void strToNT(const char* input, float& output) {
  output = (float)atof(input);
}

template<>
inline void strToNT(const char* input, double& output) {
  output = atof(input);
}

template <typename NT>
class Polynomial : public 
#ifndef CORE_DISABLE_REFCOUNTING 
PolyBase<NT> 
#else
RcPolyBase<NT> 		// reference count version
#endif
{
#ifndef CORE_DISABLE_REFCOUNTING 
typedef PolyBase<NT> base_cls;
#else
typedef RcPolyBase<NT> base_cls; // reference count version
#endif
private:

/// METHODS USED BY STRING CONSTRUCTOR
//{@

//Sets the input Polynomial to X^n

void constructX(int n, Polynomial<NT>& P){
	Polynomial<NT> q(n);//Nominal degree n
	q.setCoeff(n,NT(1));
	if (n>0) q.setCoeff(0,NT(0));
	P = q;
	}

/* REMARKS:
 * (1) Strictly speaking, the parsing for "getnumber" ought to be
 * 	done by the number type NT. Future work.
 * (2) We generalize the input number to allow at most one "." or
 * 	at most one "/".  So you can input "12.34" or "12/34" 
 * 	but not "1.2/34" or "1.2/3.4".  In any case, we let the
 * 	number type NT do the parsing (including detecting error).
 * */
///Returns in P the coefficient starting from start
int getnumber(const char* c, int start, unsigned int len,
		  Polynomial<NT> & P){
	  int j=0;
	  char *temp = new char[len];
	  while(is_float_or_rational(c[j+start])){
	    temp[j]=c[j+start];j++;
	  }
	  temp[j] = '\0';
	  // Jihun, 2010. This does not work for NT to be machine numbertype such as INT, LONG
	//  NT cf = NT(temp);
	  NT cf; strToNT<NT>(temp,cf);
	
	  Polynomial<NT> q(0);
	  q.setCoeff(0, cf);
	  P = q;
	  delete[] temp;
	  return (j-1+start);//Exactly the length of the number
	}

///checks whether a character can be in a float value or rational value.
bool is_float_or_rational(char c){
	  if (c == '0' || c == '1' || c == '2' || c == '3' || c == '4' ||
		c == '5' || c == '6' || c == '7' || c == '8' || c == '9')
		return true;
	  if (c == '.' || c == '/') return true;
	  return false;
	}

///Returns as integer the number starting from start in c

int getint(const char* c, int start, unsigned int len, int & n){
	int j=0;
	char *temp = new char[len];
	while(is_float_or_rational(c[j+start])){
	temp[j]=c[j+start];j++;
	}
	temp[j] = '\n';
	n = atoi(temp);
	delete[] temp;
	return (j-1+start);//Exactly the length of the number
	}

///Given a string starting with an open parentheses returns the place
/// which marks the end of the corresponding closing parentheses.
//Strings of the form (A).
int matchparen(const char* cstr, int start){
	int count = 0;
	int j=start;
	
	do{
	if(cstr[j] == '('){
	count++;
	}
	if(cstr[j] == ')'){
	count--;
	}
	j++;      
	}while(count != 0 );//j is one more than the matching ')'
	
	return j-1;
	}

int getbasicterm(std::string & s, Polynomial<NT> & P){
	const char * cstr = s.c_str();
	unsigned int len = s.length();
	int i=0;
	//Polynomial<NT> * temp = new Polynomial<NT>();
	
	if(is_float_or_rational(cstr[i])){
	i = getnumber(cstr, i, len, P);
	}else if(cstr[i] == 'x'||cstr[i] == 'X'){
	constructX(1, P);
	}else if(cstr[i] =='('){
	int oldi = i;
	i = matchparen(cstr, i);
	std::string t = s.substr(oldi+1, i -oldi -1);
	P = getpoly(t);
	}else{
	std::cout <<"ERROR IN PARSING BASIC TERM" << std::endl;
	}
	//i+1 points to the beginning of next syntactic object in the string.
	if(cstr[i+1] == '^'){
	int n;
	i = getint(cstr, i+2, len, n);
	P.power(n);
	}
	return i;
	}

int getterm(std::string & s, Polynomial<NT> & P){
	unsigned int len = s.length();
	if(len == 0){// Zero Polynomial
	P=Polynomial<NT>();
	return 0;
	}
	unsigned int ind, oind;
	const char* cstr =s.c_str();
	std::string t;
	//P will be used to accumulate the product of basic terms.
	ind = getbasicterm(s, P);
	while(ind != len-1 && cstr[ind + 1]!='+' && cstr[ind + 1]!='-' ){
	//Recursively get the basic terms till we reach the end or see
	// a '+' or '-' sign.
	if(cstr[ind + 1] == '*'){
	t = s.substr(ind + 2, len - ind -2);
	oind = ind + 2;
	}else{
	t = s.substr(ind + 1, len -ind -1);
	oind = ind + 1;
	}
	
	Polynomial<NT> R;
	ind = oind + getbasicterm(t, R);//Because the second term is the offset in
	//t
	P *= R;
	}
	
	return ind;
	}

Polynomial<NT> getpoly(std::string & s){
	//Remove white spaces from the string
	size_t cnt=s.find(' ',0);
	while(cnt != std::string::npos){
		s.erase(cnt, 1);
		cnt = s.find(' ', cnt);
	}

	size_t len = s.length();
	if(len <= 0){//Zero Polynomial
		return Polynomial<NT>();
	}

	//To handle the case when there is one '=' sign
	//Suppose s is of the form s1 = s2. Then we assign s to
	//s1 + (-1)(s2) and reset len
	size_t loc;
	if((loc=s.find('=',0)) != std::string::npos){
		s.replace(loc,1,1,'+');
		std::string s3 = "(-1)(";
		s.insert(loc+1, s3);
		len = s.length();
		s.insert(len, 1, ')');
	}
	len = s.length();

	const char *cstr = s.c_str();
	std::string t;
	Polynomial<NT> P;
	// P will be the polynomial in which we accumulate the
	// sum and difference of the different terms.
	unsigned int ind;
	if(cstr[0] == '-'){
		t = s.substr(1, len);
		ind = getterm(t,P) + 1;
		P.negate();
	}else{
		ind = getterm(s, P);
	}
	unsigned int oind =0;//the string between oind and ind is a term
	while(ind != len -1){
		Polynomial<NT> R;
		t = s.substr(ind + 2, len -ind -2);
		oind = ind;
		ind = oind + 2 + getterm(t, R);
		if(cstr[oind + 1] == '+')
			P += R;
		else if(cstr[oind + 1] == '-')
			P -= R;
		else
			std::cerr << "ERROR IN PARSING POLY! " << std::endl;
	}
	
	    return (P);
	  }

//This is the main function to call to construct
//	a polynomial from a string:
//
void constructFromString(std::string & s, char myX) {
	    if(myX != 'x' || myX != 'X'){
	      //Replace myX with 'x'.
	      size_t loc = s.find(myX, 0);
	      while(loc != std::string::npos){
		s.replace(loc,1,1,'x');
		loc = s.find(myX, loc+1);
	      }
	    }
	
	    //coeff = NULL;//Did this to ape the constructor from polynomial above
	    *this = getpoly(s);
	  }

  //@}
  
public:
  typedef std::vector<NT> VecNT;
  /// \name Constructors and Destructor
  //@{
  /// default constructor (zero polynomial)
  Polynomial() {}
  /// copy constructor
  Polynomial(const Polynomial &rhs) : base_cls(rhs) {}
  /// copy constructor (from possibly different NT polynomial type)
  template <class T>
    Polynomial(const Polynomial<T> &rhs) : base_cls(rhs) {}
  /// constructor of the Unit Polynomial of nominal deg n>=0
  Polynomial(int n) : base_cls(n) {}
  /// constructor with coeff array
  Polynomial(int n, NT* coef) : base_cls(n, coef) {}
  /// constructor with coeff vector
  Polynomial(const VecNT & coef) : base_cls(coef.size(), coef) {}
  Polynomial(int n, const char* s[]) : base_cls(n, s) {}
  /// constructor from char*.
  ///The BNF syntax is the following:-
  ///    [poly] -> [term]| [term] '+/-' [poly] |
  ///    		'-' [term] | '-' [term] '+/-' [poly]
  ///    [term] -> [basic term] | [basic term] [term] | [basic term]*[term]
  ///    [basic term] -> [number] | 'x' | [basic term] '^' [number]
  ///                    | '(' [poly] ')' 
  ///COMMENT: 
  ///  [number] is assumed to be a BigInt; in the future, we probably
  ///  want to generalize this to BigFloat, etc.
  ///
  Polynomial(const char* s, char myX='x') {
    std::string ss(s);
	constructFromString(ss, myX);
  }
  /// constructor from std::string
  Polynomial(const std::string& s, char myX='x') {
    std::string ss(s);
	constructFromString(ss, myX);
  }
  //@}

  /// \name help functions
  //@{
  /// return coeff (const)
  const NT* coeff() const { return base_cls::coeff(); }
  /// return coeff (non-const)
  NT* coeff() { return base_cls::coeff(); }
  /// return the degree
  int degree() const { return base_cls::degree(); }
  /// return the degree (same as degree())
  int getDegree() const { return base_cls::degree(); }
  /// return the true degree
  int getTrueDegree() const {
    int i = getDegree();
	NT c;
    while (i>=0 && sign(coeff()[i]) == 0) --i;// We are assuming that i>=0 is tested first
	/*  if(i >= 0) // The commented line above is correct. We need to check before the while loop whether i >= 0. Vikram Jan 2013
	  {
		c = coeff()[i];
    	while (i>=0 && sign(c) == 0) {
	    --i;
	    c = coeff()[i];
		}
	  }*/
    return i;
  }

  /// return the coeff
  const NT* getCoeffs() const
  { return coeff(); }
  /// return the coeff
  NT* getCoeffs() 
  { return coeff(); }

  /// set the entire coefficent array
  /// -- n should be the degree and cc's length should n+1
  /// -- SPECIAL CASE: if n=-1, then cc is ignored (so could be null)
  void setCoeffs(int n, const NT& cc)
  { return set(n, cc); }

  /// return the i-th coeff
  const NT& getCoeff(int i) const
  { assert(i>=0 && i<=getDegree()); return coeff()[i]; }
  /// set the i-th coeff
  bool setCoeff(int i, const NT& cc)
  { if (i<0 || i>getDegree()) return false; coeff()[i] = cc; return true; }
  
  /// return the leading coeff
  const NT& getLeadCoeff() const 
  { return getCoeff(getTrueDegree()); }
  /// return the last non-zero coeff
  const NT& getTailCoeff() const {
    int i = 0;
    while (i<=getTrueDegree() && getCoeff(i) == 0)
      ++i;
    return getCoeff(i); 
  }
  /// Expands the nominal degree to n;
  ///      Returns n if nominal degree is changed to n
  ///      Else returns -2
  int expand(int n) {
    if ((n <= getDegree())||(n < 0))
      return -2;
    int i;
    NT * c = new NT[n+1];
    for (i = 0; i<= getDegree(); i++)
      c[i] = coeff()[i];
    for (i = getDegree()+1; i<=n; i++)
      c[i] = 0;
    this->set(n, c);

    // This fixes a memory leak, note that set( ) makes
    // a copy of C and does not use C by itself.
    delete[] c;
    return n;
  }
  //@}

  /// contract() gets rid of leading zero coefficients
  ///	and returns the new (true) degree;
  ///	It returns -2 if this is a no-op

  int contract() {
    int d = getTrueDegree();
    if (d == degree())
      return (-2);  // nothing to do

    NT * c = new NT[d+1];
    for (int i = 0; i<= d; i++)
      c[i] = coeff()[i];
    this->set(d, c);

    delete[] c;
    return d;
  }

  /// \name assignment and compound assignment operators (Polynomial arithmetic)
  //@{
  /// assignment operator for <tt>Polynomial</tt>
  Polynomial& operator=(const Polynomial& rhs)
  { base_cls::operator=(rhs); return *this; }
  /// compound assignment operator <tt>+=</tt>
  Polynomial& operator+=(const Polynomial& rhs) {
    int d = rhs.getDegree();
    if (d > getDegree()) expand(d);
    NT c;
    for (int i=0; i<=d; ++i){
      // coeff()[i] += rhs.coeff()[i];		//const error...
	c = rhs.coeff()[i];			// fix...
	coeff()[i] += c;
    }
    return *this;
  }
  /// compound assignment operator <tt>-=</tt>
  Polynomial& operator-=(const Polynomial& rhs) {
    int d = rhs.getDegree();
    if (d > getDegree()) expand(d);
    NT c;
    for (int i=0; i<=d; ++i){
      // coeff()[i] -= rhs.coeff()[i];		// const error..
	c = rhs.coeff()[i];			// fix...
	coeff()[i] -= c;
    }
    return *this;
  }
  /// compound assignment operator <tt>*=</tt>
  Polynomial& operator*=(const Polynomial& rhs);
  
  static  int COEFF_PER_LINE;           // pretty print parameters
  static const char * INDENT_SPACE;             // pretty print parameters

  /// Four special polynomials:
  /// Zero() = 0
  static const Polynomial<NT> & Zero();
  /// Unity() = 1
  static const Polynomial<NT> & Unity();
  /// NegUnity() = -1
  static const Polynomial<NT> & NegUnity();
  /// Identity = x
  static const Polynomial<NT> & Identity();

  ///Multiply by a scalar
  Polynomial<NT> & mulScalar( const NT & c) {
	if (c==NT(0)) {
		*this = Polynomial<NT>();
		return *this;
	}

    if(c!= NT(1)){
		NT cc = c;
	
		for (int i = 0; i<=degree() ; i++){
			//	coeff()[i] *= c;		// const error..
			coeff()[i] *= cc;				// fix...
		}
	}
	  return *this;
  }

  ///exactDivScalar is to divide exactly by a scalar
  ///  It is assumed that NT has an exact division operation
  // NOTE: this function is not fully tested or used currently
  Polynomial<NT> & exactDivScalar( const NT & c) {
    for (int i = 0; i<=degree() ; i++)
      coeff()[i] /= c;
    return *this;
  }

  /// mulXpower: Multiply by X^i (COULD be a divide if i<0)
  /// if s >= 0, then this is equivalent to
  /// multiplying by X^s;  if s < 0, to dividing by X^s
  Polynomial<NT> & mulXpower(int s) {

  if (s==0)
    return *this;
  int d = s+getTrueDegree();
  if (d < 0) {
    *this = Polynomial<NT>();
    return *this;
  }
  NT * c = new NT[d+1];
  if (s>0)
    for (int j=0;  j <= d; j++) {
      if (j <= degree())
        c[d-j] = coeff()[d-s-j];
      else
        c[d-j] = 0;
    }
  if (s<0) {
    for (int j=0; j <= d; j++)
      c[d-j] = coeff()[d-s-j];  // since s<0, (d-s-j) > (d-j) !!
  }
  this->set(d,c);
  
  delete[] c;
  return *this;
  }//mulXpower


/// REDUCE STEP (helper for PSEUDO-REMAINDER function)
/// Let THIS=(*this) be the current polynomial, and P be the input
///	argument for reduceStep.  Let R be returned polynomial.
///	R has the special form as a binomial,
///		R = C + X*M
///	where C is a constant and M a monomial (= coeff * some power of X).
///	Moreover, THIS is transformed to a new polynomial, THAT, which
///	is given by
/// 		(C * THIS) = M * P  + THAT
///	MOREOVER: deg(THAT) < deg(THIS) unless deg(P)>deg(Q).
///	Basically, C is a power of the leading coefficient of P.
///	REMARK: R is NOT defined as C+M, because in case M is
///	a constant, then we cannot separate C from M.
///	Furthermore, R.mulScalar(-1) gives us M.
Polynomial<NT> reduceStep (
  const Polynomial<NT>& p) {
  // 	Chee: Omit the next 2 contractions as unnecessary
  // 	since reduceStep() is only called by pseudoRemainder().
  // 	Also, reduceStep() already does a contraction before returning.
  // p.contract();	
  // contract();	// first contract both polynomials
  Polynomial<NT> q(p);		// q is initially a copy of p
  //	but is eventually M*P
  int pDeg  = q.degree();
  int myDeg = degree();
  if (pDeg == -1)
    return *(new Polynomial());  // Zero Polynomial
  // NOTE: pDeg=-1 (i.e., p=0) is really an error condition!
  if (myDeg < pDeg)
    return *(new Polynomial(0));  // Unity Polynomial
  // i.e., C=1, M=0.
  // Now (myDeg >= pDeg).  Start to form the Return Polynomial R=C+X*M
  Polynomial<NT> R(myDeg - pDeg + 1);  // deg(M)= myDeg - pDeg
  q.mulXpower(myDeg - pDeg);  	 // q is going to become M*P

  NT myLC = coeff()[myDeg];	  // LC means "leading coefficient"
  NT qLC = q.coeff()[myDeg];  // p also has degree "myDeg" (qLC non-zero)
  NT LC;

  //  NT must support
  //  isDivisible(x,y), gcd(x,y), div_exact(x,y) in the following:
  //  ============================================================
  if (isDivisible(myLC, qLC)) { // myLC is divisible by qLC
    LC = div_exact(myLC, qLC);	 
    R.setCoeff(0, 1);  		 //  C = 1,

    R.setCoeff(R.degree(), LC); //  M = LC * X^(myDeg-pDeg)
    q.mulScalar(LC); 	  //  q = M*P. 
  }
  else if (isDivisible(qLC, myLC)) { // qLC is divisible by myLC
    LC = div_exact(qLC, myLC);	 //
    if ((LC != 1) && (LC != -1)) { // IMPORTANT: unlike the previous
      // case, we need an additional condition
      // that LC != -1.  THis is because
      // if (LC = -1), then we have qLC and
      // myLC are mutually divisible, and
      // we would be updating R twice!
      R.setCoeff(0, LC); 	   // C = LC, 
      R.setCoeff(R.degree(), 1);     // M = X^(myDeg-pDeg)(THIS WAS NOT DONE)
      mulScalar(LC); 	   	   // THIS => THIS * LC

    }
  } else {  			// myLC and qLC are not mutually divisible
    NT g = gcd(qLC, myLC); 	// This ASSUMES gcd is defined in NT !!
    //NT g = 1;  			// In case no gcd is available
    if (g == 1) {
      R.setCoeff(0, qLC);	  	// C = qLC
      R.setCoeff(R.degree(), myLC);	 // M = (myLC) * X^{myDeg-pDeg}
      mulScalar(qLC);	 		// forming  C * THIS
      q.mulScalar(myLC);		// forming  M * P
    } else {
      NT qLC1= div_exact(qLC,g);
      NT myLC1= div_exact(myLC,g);
      R.setCoeff(0, qLC1);	  	// C = qLC/g
      R.setCoeff(R.degree(), myLC1);	// M = (myLC/g) * X^{myDeg-pDeg}
      mulScalar(qLC1);	 	// forming  C * THIS
      q.mulScalar(myLC1);		// forming  M * P
    }
  }
  (*this) -= q;		// THAT = (C * THIS) - (M * P)

  contract();

  return R;		// Returns R = C + X*M
}// reduceStep

/// Poly prem(Poly B) 
/// Poly pseudoRemainder(Poly B) 
/// Poly pseudoRemainder(Poly B, NT C)
///	-- 3 versions of pseudo remainder (prem) and pseudo quotient (pquo)
///     -- it modifies *this into prem(*this,B) and returns pquo.
/// BUT "pseudoRemainder" is not exactly the "prem" of theory.
///     Basically, "A.prem(B)" :=
///		A *= (lc(B))^{deg(A)-deg(B)+1};
///		A.pseudoRemainder(B);
//
/// THEORY: PSEUDO-REMAINDER and PSEUDO-QUOTIENT.
/// Let A = (*this) and B be the argument polynomial.
/// Let Quo be the returned polynomial, 
/// and let the final value of (*this) be Rem.
/// Also, C is the constant that we maintain.
/// We are computing A divided by B.  The relation we guarantee is 
/// 		(C * A) = (Quo * B)  + Rem
/// where deg(Rem) < deg(B).  So Rem is the Pseudo-Remainder
/// and Quo is the Pseudo-Quotient.
/// Moreover, C is uniquely determined (we won't spell it out)
/// except to note that
///	C divides D = (LC)^{deg(A)-deg(B)+1}
///	where LC is the leading coefficient of B.
/// NOTE: 1. Normally, Pseudo-Remainder is defined so that C = D
/// 	 So be careful when using our algorithm.
/// 	 2. We provide a version of pseudoRemainder which does not
/// 	 require C as argument.  [For efficiency, we should provide this
///	 version directly, instead of calling the other version!]

Polynomial<NT> pseudoRemainder (
  const Polynomial<NT>& B) {
	NT temp;	// dummy argument to be discarded
	return pseudoRemainder(B, temp);
}//pseudoRemainder

Polynomial<NT> prem (
  const Polynomial<NT>& B) {
  int m=getTrueDegree();
  int n=B.getTrueDegree();
  if (m<n) {//switch roles of A and B
    Polynomial<NT> BB=*this;
    (*this) = B;
    B = BB;
    int mm=n;
    n = m;  m = mm;
  }  
  (*this) *= power(B.coeff()[n],m-n+1);
  return pseudoRemainder(B);
}//prem

Polynomial<NT> pseudoRemainder (
  const Polynomial<NT>& B, NT & C) { 
  contract();         // Let A = (*this).  Contract A.
  Polynomial<NT> tmpB(B);
  tmpB.contract();    // local copy of B
  C = *(new NT(1));  // Initialized to C=1.
  int Bdeg = B.degree();
  
  if (Bdeg == -1)  {
    std::cout << "ERROR in Polynomial<NT>::pseudoRemainder :\n" <<
    "    -- divide by zero polynomial" << std::endl;
    return Polynomial(0);  // Unit Polynomial (arbitrary!)
  }
  int deg = degree();
  
  if (Bdeg > deg) {
    return Polynomial<NT>(); // Zero Polynomial
    // CHECK: 1*THIS = 0*B + THAT,  deg(THAT) < deg(B)
  }

  Polynomial<NT> Quo;  // accumulate the return polynomial, Quo
  Polynomial<NT> tmpQuo;
  while (degree() >= Bdeg) {  // INVARIANT: C*A = B*Quo + (*this)
    tmpQuo = reduceStep(tmpB);  // Let (*this )be (*oldthis), which
			        // is transformed into (*newthis). Then,
                                //     c*(*oldthis) = B*m + (*newthis)
                                // where tmpQuo = c + X*m
    // Hence,   C*A =   B*Quo +   (*oldthis)      -- the old invariant
    //        c*C*A = c*B*Quo + c*(*oldthis)
    //              = c*B*Quo + (B*m + (*newthis))
    //              = B*(c*Quo + m)  + (*newthis)
    // i.e, to update invariant, we do C->c*C,  Quo --> c*Quo + m.
    C *= tmpQuo.coeff()[0];	    // C = c*C
    Quo.mulScalar(tmpQuo.coeff()[0]); // Quo -> c*Quo
    tmpQuo.mulXpower(-1);           // tmpQuo is now equal to m
    Quo += tmpQuo;                  // Quo -> Quo + m 
  }

  return Quo;	// Quo is the pseudo-quotient
}//pseudoRemainder

/// Returns the negative of pseudoRemainder
/// 	(self-modification)
Polynomial<NT> & negPseudoRemainder (
  const Polynomial<NT>& B) {
	NT temp;	// dummy argument to be discarded
	pseudoRemainder(B, temp);
	if (temp < 0) return (*this);
	return negate();
}

///Unary minus.
//Used in construction from strings.
Polynomial<NT> & operator-() {	// unary minus
  for (int i=0; i<=degree(); i++)
    coeff()[i] *= -1;
  return *this;
}
///Self-power.
//Used in construction from strings.
Polynomial<NT> & power(unsigned int n) {	// self-power
  if (n == 0) {
    NT* c = new NT[1];
    c[0] = 1;
    this->set(0, c);
    delete[] c;
  } else {
    Polynomial<NT> p = *this;
    for (unsigned int i=1; i<n; i++) {
      *this *= p;		// Would a binary power algorithm be better?
    }
  }
  return *this;
}

/// differentiation
Polynomial<NT> & differentiate() {	// self-differentiation
  contract();
  int deg = degree();
  if (deg >= 0) {
    NT * c = new NT[deg];
    for (int i=1; i<=deg; i++)
      c[i-1] = coeff()[i] * i;
    deg--;
    this->set(deg,c);

    // This fixes a memory leak, note that set( ) makes a copy
    // of c and does not use c by itself.
    delete[] c;
  }

  return *this;
}// differentiation

/// multi-differentiate
Polynomial<NT> & differentiate(int n) {
  assert(n >= 0);
  for (int i=1; i<=n; i++)
    this->differentiate();
  return *this;
} // multi-differentiate

  // Reductions of polynomials (NT must have gcd function)
  Polynomial sqFreePart(); // Square free part of P is P/gcd(P,P'). Return gcd
  Polynomial & primPart();   // Primitive Part of *this (which is changed)


  //@}

/// Useful member functions.
//{@

/// 	Reverses the list of coefficients
void reverse() {
  NT tmp;
  for (int i=0; i<= degree()/2; i++) {
    tmp = coeff()[i];
    coeff()[i] =   coeff()[degree()-i];
    coeff()[degree()-i] = tmp;
  }
}//reverse

/// negate: multiplies the polynomial by -1.
// 	Chee: 4/29/04 -- added negate() to support negPseudoRemainder(B)
Polynomial<NT> & negate() {
  for (int i=0; i<= degree(); i++) 
    coeff()[i] *= -1;  	// every NT must be able to construct from -1
  return *this;
}//negate

/// makeTailCoeffNonzero
/// 	Divide (*this) by X^k, so that the tail coeff becomes non-zero.
///	The value k is returned.  In case (*this) is 0, return k=-1.
///	Otherwise, if (*this) is unchanged, return k=0.
int makeTailCoeffNonzero() {
  int k=-1;
  int deg = degree();
  for (int i=0; i <= deg; i++) {
    if (coeff()[i] != 0) {
      k=i;
      break;
    }
  }
  if (k <= 0)
    return k;	// return either 0 or -1
  deg -=k;		// new (lowered) degree
  NT * c = new NT[1+deg];
  for (int i=0; i<=deg; i++)
    c[i] = coeff()[i+k];

  this->set(deg, c);

  delete[] c;
  return k;
}//

/// filedump(string msg, ostream os, string com, string com2):
///      Dumps polynomial to output stream os
///      msg is any message.
//      NOTE: Default is com="#", which is placed at start of each 
//            output line. 
void filedump(std::ostream & os,
        std::string msg,
	std::string commentString,
        std::string commentString2) const {
	  int d= getTrueDegree();
	  if (msg != "") os << commentString << msg << std::endl;
	  int i=0;
	  if (d == -1) { // if zero polynomial
	    os << commentString << "0";
	    return;
	  }
	  for (; i<=d;  ++i)	// get to the first non-zero coeff
	    if (coeff()[i] != 0)
	      break;
	  int termsInLine = 1;
	
	  // OUTPUT the first nonzero term
	  os << commentString;
	  if (coeff()[i] == 1) {	// special cases when coeff[i] is
	    if (i>1) os << "x^" << i;	// either 1 or -1 
	    else if (i==1) os << "x" ;
	    else os << "1";
	  } else if (coeff()[i] == -1) {
	    if (i>1) os << "-x^" << i;
	    else if (i==1) os << "-x" ;
	    else os << "-1";
	  } else {				// non-zero coeff
	    os << coeff()[i];
	    if (i>1) os << "*x^" << i;
	    else if (i==1) os << "x" ;
	  } 
	  // OUTPUT the remaining nonzero terms
	  for (i++ ; i<= getTrueDegree(); ++i) {
	    if (coeff()[i] == 0) 
	      continue;
	    termsInLine++;
	    if (termsInLine % Polynomial<NT>::COEFF_PER_LINE == 0) {
	      os << std::endl;
	      os << commentString2;
	    }
	    if (coeff()[i] == 1) {	// special when coeff[i] = 1
	      if (i==1) os << " + x";
	      else os << " + x^" << i;
	    } else if (coeff()[i] == -1) { // special when coeff[i] = -1
	      if (i==1) os << " - x";
	      else os << " -x^" << i;
	    } else {			// general case
	      if(coeff()[i] > 0){	
		os << " + ";
	        os << coeff()[i];
	      }else
	        os << coeff()[i];
	
	      if (i==1) os << "*x";
	      else os << "*x^" << i;
	    } 
	  }
	}//filedump

/// dump(ofstream, message, com1, com2) -- dump to output stream
/// NOTE: see also toString(c) for a more estatically pleasing output!
void dump(std::ofstream & ofs,
	std::string msg = "",
	std::string commentString = "# ",
        std::string commentString2 = "# ") const {
  filedump(ofs, msg, commentString, commentString2);
}

/// dump(message) 	-- to std output
/// NOTE: see also toString(c) for a more estatically pleasing output!
void dump(std::string msg = "", std::string com = "# ",
	std::string com2 = "# ") const {
  filedump(std::cout, msg, com, com2);
}

/// Dump of Maple Code for Polynomial
void mapleDump() const {
	  if (zeroP(*this)) {
	    std::cout << 0 << std::endl;
	    return;
	  }
	  std::cout << coeff()[0];
	  for (int i = 1; i<= getTrueDegree(); ++i) {
	    std::cout << " + (" << coeff()[i] << ")";
	    std::cout << "*x^" << i;
	  }
	  std::cout << std::endl;
	}//mapleDump

/// toString(c) 
///    c is the optional char for variable name (default to 'x').
///    returns a string representation of polynomial, starting with
///    the leading coefficient.  See also toStringRev(c)
///    This is the inverse of Polynomail<NT>::getpoly(string).
std::string toString(char var='x'){
  int d= getTrueDegree();
  //
  if (d == -1) return std::string("0"); 
  // FOR EXTRA ROBUSTNESS,
  // 	we check the possibility that
  //	the zero polynomial is regarded as having degree 0:
  if (d == 0 && getCoeff(0) == 0) return std::string("0"); 
  int i=d;
  std::ostringstream oss;
  std::string s;
  NT c;
  c = getCoeff(i);	// first term
  if (d==0) {// we now know that c != 0
	oss << c;
	s = oss.str();
	return s;}
  if (d==1) {
	if (c == 1) oss << var;
	else if (c== -1) oss << "-" << var;
  	else oss << c << var;
  } else 	// d>1
	if (c == 1) oss << var << "^" << d;
	else if (c == -1) oss << "-" << var << "^" << i;
  	else oss << c << var << "^" << i;
  i=d; i--;
  for (; i>= 2; i--) {	// general term
	c = getCoeff(i);
	if (c != 0) {
	  if (c == 1)
  		oss << " + " << var << "^" << i;
	  else if (c == -1)
  		oss << " - " << var << "^" << i;
	  else if (c > 0)
  		oss << " + " << c << var << "^" << i;
	  else {
		c *= NT(-1);
  		oss << " - " << c << var << "^" << i;
	  }}} //for i>=2
  if (d > 1) { // if d=1, we must not output anything (else duplicate)
	  c = getCoeff(1); // linear term
	  if (c != 0) {
		if (c == 1)
	  		oss << " + " << var;
		else if (c == -1)
	  		oss << " - " << var;
		else if (c > 0)
	  		oss << " + " << c << var;
		else {
			c *= NT(-1);
	  		oss << " - " << c << var;
		}}} // if d>=1
  if (d > 0) { // if d=0, must not output here (else duplicate)
	   c = getCoeff(0); // constant term
	   if (c != 0) {
		if (c>0)
	  		oss << " + " << c;
		else {
			c *= NT(-1);
	  		oss << " - " << c;
		}}} // if d>=0
  s=oss.str();
  return s;
}//toString()

/// toStringRev(c)
///    This is the version of toString() where we begin with the
///    tail coefficient and end up in the leading coeffient.
std::string toStringRev(char var='x'){
  int d= getTrueDegree();
  //
  if (d == -1) return std::string("0"); 
  std::ostringstream oss;
  bool seenFirstTerm = false;
  NT c;
  if (d >= 0){ 
  	c = getCoeff(0);
	if (c != 0) {
	  seenFirstTerm = true;
	  oss << c;
	}}
  if (d >= 1){
	c = getCoeff(1);
	if (c != 0) {
	  if (seenFirstTerm) {
		  if (c>0) oss << " + ";
		  else oss << " - ";
	  } else
		  seenFirstTerm = true;
	  if (c == 1 || c == -1) oss << var;
	  else {
		  if (c<0) c *= NT(-1);
		  oss << c << var;
	  }}}
  int i=2; 
  for (; i<= d; i++) {
	c = getCoeff(i);
	if (c != 0) {
	  if (seenFirstTerm) {
		if (c>0) oss << " + ";
		else oss << " - ";
	  } else
		  seenFirstTerm = true;
	  if (c < 0) c *= NT(-1);
	  if (c == 1) oss << var << "^" << i;
	  else oss << c << var << "^" << i;
  	}}
  std::string s=oss.str();
  return s;
}//toStringRev


//@}

///Evaluation: Assumes that the point of evaluation has the same
///type as the coefficients.
template<typename T>
T eval(const T& x) const {	// evaluation
  int deg = degree();
  if (deg == -1)
    return T(0);
  if (deg == 0)
    return T(coeff()[0]);

  T val(coeff()[deg]);
  for (int i = deg-1; i>=0; i--) {
    val *= x;
    val += T(coeff()[i]);	
  }
  return val;
}//eval

///Evaluation: Assumes that the interval type is the same
///type as the coefficients.
// Horner's rule based interval evaluation scheme. 
// NOT REQUIRED. Invoke eval<Interval<T> > instead.
/*template<typename T>
IntervalT<T> evalH(const IntervalT<T> & x) const {	// interval evaluation
	int deg = degree();
	if (deg == -1)
		return IntervalT<T>(0);
	if (deg == 0)
		return IntervalT<T>(coeff()[0]);

	IntervalT<T> val(coeff()[deg]);
	for (int i = deg-1; i>=0; i--) {
		val *= x;
		val += IntervalT<T>(coeff()[i]);	
	}
	return val;
		
}//evalH
*/
	
///Evaluation: Assumes that the interval type is the same
///type as the coefficients.
// Slope form based interval evaluation scheme.
// Drawback: This is not inclusion monotone, i.e., if interval I is a included
// in another interval J, then evalSlopeForm(I) is not necessarily included in evalSlopeForm(J). See Stahl's thesis.
// REMARK: this is a misnomer -- it should be called the "Slope Form"
//          See Stahl's thesis, p.72.
//          The Centered Form is called Taylor Form in Stahl.
//          Should also implement what Stahl calls "Bicentered Mean Value Form".
template<typename T>
IntervalT<T> evalSlopeForm(const IntervalT<T> & x) const {	// interval evaluation
	int deg = degree();
	if (deg == -1)
		return IntervalT<T>(0);
	if (deg == 0)
		return IntervalT<T>(coeff()[0]);
	T xmid = x.mid();
	
	// tmp is the quotient obtained by dividing the polynomial (*this) by (x-xmid).
	// (*this) = (*this).eval(xmid) + (x-xmid)*tmp.
	Polynomial<T> tmp(deg-1);
	tmp.coeff()[deg-1] = (*this).coeff()[deg];
	for (int i=deg-1; i>0; i--) {
		tmp.coeff()[i-1] = (*this).coeff()[i]+xmid*tmp.coeff()[i];
	}

	T fmid = (*this).coeff()[0] + xmid*tmp.coeff()[0];// val is (*this) evaluated at xmid
	
	IntervalT<T> val =
	    IntervalT<T>(fmid)+(x-IntervalT<T>(xmid))*tmp.template eval<IntervalT<T> >(x);
	
	return val;
		
	}//eval
	

	///Evaluation: Assumes that the interval type is the same
	///type as the coefficients.
	// Centered form evaluation. This approach does a Taylor shift
	// at the midpoint of the interval rather than differentiating and
	// computing all the derivatives.
	template<typename T>
	IntervalT<T> evalCenteredForm(const IntervalT<T> & x) const {	// interval evaluation
		int deg = degree();
		if (deg == -1)
			return IntervalT<T>(0);
		if (deg == 0)
			return IntervalT<T>(coeff()[0]);
		T xmid = x.mid();

		// Compute the Taylor coefficients at xmid
		T shifted[deg+1];
		for(int i=0; i<= deg ; i++)
			shifted[i] = coeff()[i];
		
		if(xmid != 0){
			for(int i=0; i<= deg-1;i++)
				for(int j=deg-1; j>=i; j--)
					shifted[j]+= (xmid) * shifted[j+1];
		}
				
		T w = x.halfwidth(); // Half the width of interval
		
		// Compute the error term in centered form.
		// The absolute value of the Taylor coefficients is used.
		T val = abs(shifted[deg]);
		for (int i = deg-1; i>0; i--) {
			val *= w;
			val += abs(shifted[i]); 	
		}
		val*=w;
		return IntervalT<T>(shifted[0]-val, shifted[0]+val);		
	}		
	
template<typename T>
int evalSign(const T& f) {
  int deg = getTrueDegree();
  if (deg == -1)
    return 0;
  if (deg == 0 || f == 0)
    return sign(coeff()[0]);
  T val(coeff()[deg]);
  for (int i=deg-1; i>=0; i--) {
    val *= f;
    val += T(coeff()[i]);	
  }
  return sign(val);
}//evalSign 

template<typename T>
  T operator() ( const T& x ) const {  return eval( x ) ; }

}; // class Polynomial 

/// Polynomial + Polynomial
template <typename NT>
inline 
Polynomial<NT> operator+(const Polynomial<NT>& x, const Polynomial<NT>& y) {
  return Polynomial<NT>(x) += y;
}
/// Polynomial - Polynomial
template <typename NT>
inline 
Polynomial<NT> operator-(const Polynomial<NT>& x, const Polynomial<NT>& y) {
  return Polynomial<NT>(x) -= y;
}
/// Polynomial * Polynomial
template <typename NT>
inline 
Polynomial<NT> operator*(const Polynomial<NT>& x, const Polynomial<NT>& y) {
  if ( x.getDegree() < 0 || y.getDegree() < 0 ) 
    return Polynomial<NT>(std::vector<NT>(1, NT(0)));

  int d = x.getDegree() + y.getDegree();
  NT* c = new NT[d+1];
  for(int i=0; i<d+1; ++i)
    c[i] = NT(0);

  for (int i=0; i<=x.getDegree(); ++i)
    for (int j=0; j<=y.getDegree(); ++j)
      c[i+j] += x.coeff()[i] * y.coeff()[j];
  Polynomial<NT> result(d, c);
  delete[] c;
  return result;
}// operator*

template <typename NT>
inline
Polynomial<NT>& Polynomial<NT>::operator*=(const Polynomial<NT>& rhs) {
  *this = *this * rhs; return *this;
}

template < class NT >
inline Polynomial<NT> power(const Polynomial<NT>& p, int n) { // power
  return Polynomial<NT>(p).power(n);
}

template <class NT>
int Polynomial<NT>::COEFF_PER_LINE  = 4;           // pretty print parameters
template <class NT>
const char* Polynomial<NT>::INDENT_SPACE ="   ";  // pretty print parameters

// Four Special (static) Polynomials: 0, 1, -1, x:
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//Zero() = 0
template < class NT >
inline const Polynomial<NT> & Polynomial<NT>::Zero() {
  static Polynomial<NT> pZ;	// Should not use "zeroP" which is the predicate
  				// to test if a polynomial is zero!
  return pZ;
}//Zero

//Unity() = 1
template < class NT >
inline const Polynomial<NT> & Polynomial<NT>::Unity() {
  static NT c[] = {1};
  static Polynomial<NT> one(0, c);
  return one;
}//Unity

//NegUnity() = -1
template < class NT >
inline const Polynomial<NT> & Polynomial<NT>::NegUnity() {
  static NT c[] = {-1};
  static Polynomial<NT> negOne(0, c);
  return negOne;
}//Unity

//Identity() = x
template < class NT >
inline const Polynomial<NT> & Polynomial<NT>::Identity() {
  static NT c[] = {0, 1};
  static Polynomial<NT> ident(1, c);
  return ident;
}//Identity

// inline functions
/// zeroP(p) is true if p is identically zero.
template <typename NT>
inline bool zeroP(const Polynomial<NT>& p)
{ return p.getTrueDegree() == -1;
}// zeroP

/// unitP(p) is true if p is the unit element 1
template <typename NT>
inline bool unitP(const Polynomial<NT>& p)
{ return p.getTrueDegree() == 0 && p.coeff()[0] == 1;
}// unitP

// stream i/o
template <typename NT>
std::ostream& operator<<(std::ostream& o, const Polynomial<NT>& p) {
  o <<   "Polynomial<NT> ( deg = " << p.getDegree() ;
  if (p.getDegree() >= 0) {
    o << "," << std::endl;
    o << ">  coeff c0,c1,... = " << p.coeff()[0];
    for (int i=1; i<= p.getDegree(); i++)
      o << ", " <<  p.coeff()[i] ;
  }
  o << ")" << std::endl;
  return o;
}//operator<<

template < class NT >
std::istream& operator>>(std::istream& is, Polynomial<NT>& p) {
  // read degree and coeff
  int degree;
  is >> degree;
  NT* coeff = new NT[degree+1];
  for (int i=0; i<= degree; i++)
    is >> coeff[i];
  // set coefficients for poly
  p.set(degree, coeff);
  delete coeff;
  return is;
}//operator>>

// For internal use only:
// Checks that c*A = B*m + AA 
// 	where A=(*oldthis) and AA=(*newthis)
template <class NT>
inline
Polynomial<NT> testReduceStep(const Polynomial<NT>& A, 
	const Polynomial<NT>& B) {
std::cout << "+++++++++++++++++++++TEST REDUCE STEP+++++++++++++++++++++\n";
  Polynomial<NT> cA(A);
  Polynomial<NT> AA(A);
  Polynomial<NT> quo;
  quo = AA.reduceStep(B);	        // quo = c + X*m  (m is monomial, c const)
                                // where c*A = B*m + (*newthis)
std::cout << "A = " << A << std::endl;
std::cout << "B = " << B << std::endl;
  cA.mulScalar(quo.coeff[0]);    // A -> c*A
  Polynomial<NT> m(quo);
  m.mulXpower(-1);            // m's value is now m
std::cout << "c = " << quo.coeff[0] << std::endl;
std::cout << "c + xm = " << quo << std::endl;
std::cout << "c*A = " << cA << std::endl;
std::cout << "AA = " << AA << std::endl;
std::cout << "B*m = " << B*m << std::endl;
std::cout << "B*m + AA = " << B*m + AA << std::endl;
  if (cA == (B*m + AA))
	  std::cout << "CORRECT inside testReduceStep" << std::endl;
  else
	  std::cout << "ERROR inside testReduceStep" << std::endl;
std::cout << "+++++++++++++++++END TEST REDUCE STEP+++++++++++++++++++++\n";
  return quo;
}// testReduceStep

/// GCD, content, primitive and square-free parts
//{@

/// Divisibility predicate for polynomials.
/// isDivisible(P,Q) returns true iff Q divides P
// To FIX: the predicate name is consistent with Expr.h but not with BigInt.h
template <class NT>
inline bool isDivisible(Polynomial<NT> p, Polynomial<NT> q) {
  if(zeroP(p))
    return true;
  if(zeroP(q))
    return false;  // We should really return error!
  if(p.getTrueDegree() < q.getTrueDegree())
    return false;
  p.pseudoRemainder(q);
  if(zeroP(p))
    return true;
  else
    return false;
}//isDivisible

///Content of a polynomial P.
///      -- content(P) is just the gcd of all the coefficients
///      -- REMARK: by definition, content(P) is non-negative
///                 We rely on the fact that NT::gcd always
///                 return a non-negative value!!!
template <class NT>
inline NT content(const Polynomial<NT>& p) {
  if(zeroP(p))
    return 0;
  int d = p.getTrueDegree();
  if(d == 0){
    if(p.coeff()[0] > 0)
      return p.coeff()[0];
    else
      return -p.coeff()[0];
  }

  NT content = p.coeff()[d];
  for (int i=d-1; i>=0; i--) {
    content = gcd(content, p.coeff()[i]);
    if(content == 1) break;   // remark: gcd is non-negative, by definition
  }
  //if (p.coeff[0] < 0) return -content;(BUG!)
  return content;
}//content

/// Primitive Part-- primPart(P) is just P/content(P)
//  (*this) is transformed to primPart and returned
//      -- Should we return content(P) instead? [SHOULD IMPLEMENT THIS]
// IMPORTANT: we require that content(P)>0, hence
//      the coefficients of primPart(P) does 
//      not change sign; this is vital for use in Sturm sequences
template <class NT>
Polynomial<NT> & Polynomial<NT>::primPart() {
  // ASSERT: GCD must be provided by NT
  int d = getTrueDegree();
  assert (d >= 0);
  if (d == 0) {
    if (coeff()[0] > 0) coeff()[0] = 1;
    else coeff()[0] = -1;
    return *this;
  }

  NT g = content(*this);
  if (g == 1 && coeff()[d] > 0)
     return (*this);
  for (int i=0; i<=d; i++) {
    coeff()[i] =  div_exact(coeff()[i], g);
  }
  return *this;
}// primPart

///GCD of two polynomials.
///   --Assumes that the coeffecient ring NT has a gcd function
///   --Returns the gcd with a positive leading coefficient (*)
///     otherwise division by gcd causes a change of sign affecting
///     Sturm sequences.
//   --To Check: would a non-recursive version be much better?
template <class NT>
inline Polynomial<NT> gcd(const Polynomial<NT>& p, const Polynomial<NT>& q) {
  // If the first polynomial has a smaller degree then the second,
  // then change the order of calling
  if(p.getTrueDegree() < q.getTrueDegree())
    return gcd(q,p);
  // If any polynomial is zero then the gcd is the other polynomial
  if(zeroP(q)){
    if(zeroP(p))
       return p;
    else{
       if(p.getCoeff(p.getTrueDegree()) < 0){
         return Polynomial<NT>(p).negate();
       }else
         return p;	// If q<>0, then we know p<>0
   }
  }
  Polynomial<NT> temp0(p);
  Polynomial<NT> temp1(q);

  // We want to compute:
  // gcd(p,q) = gcd(content(p),content(q)) * gcd(primPart(p), primPart(q))
  NT cont0 = content(p);	// why is this temporary needed?
  NT cont1 = content(q);
  NT cont = gcd(cont0,cont1);
  temp0.primPart();
  temp1.primPart();
  temp0.pseudoRemainder(temp1);
  return (gcd(temp1, temp0).mulScalar(cont));
}//gcd

// sqFreePart()
//      -- this is a self-modifying operator!
//      -- Let P =(*this) and Q=square-free part of P.
//      -- (*this) is transformed into P, and gcd(P,P') is returned
// NOTE: The square-free part of P is defined to be P/gcd(P, P')
template <class NT>
Polynomial<NT>  Polynomial<NT>::sqFreePart() {
  
  int d = getTrueDegree();
  if(d <= 1) // linear polynomials and constants are square-free
    return *this;
  Polynomial<NT> temp(*this);
  Polynomial<NT> R = gcd(*this, temp.differentiate()); // R = gcd(P, P')
  // If P and P' have a constant gcd, then P is squarefree
  if(R.getTrueDegree() == 0)
    return (Polynomial<NT>(0)); // returns the unit polynomial as gcd
  
  (*this)=pseudoRemainder(R); // (*this) is transformed to P/R, the sqfree part
  //Note: This is up to multiplication by units
  return (R); // return the gcd
}//sqFreePart()

//@}

/// Useful friend functions for Polynomial<NT> class
//{@

/// differentiation
template <class NT>
inline Polynomial<NT> differentiate(const Polynomial<NT> & p) {	  // differentiate
  Polynomial<NT> q(p);
  return q.differentiate();
}

/// multi-differentiation
template <class NT>
inline 
Polynomial<NT> differentiate(const Polynomial<NT> & p, int n) {//multi-differentiate
  Polynomial<NT> q(p);
  assert(n >= 0);
  for (int i=1; i<=n; i++)
    q.differentiate();
  return q;
}

/// comparison (equality of two polynomials) 
template <class NT>
inline 
bool operator==(const Polynomial<NT>& p, const Polynomial<NT>& q) {	// ==

  int d = p.getTrueDegree();
	
  if (d != q.getTrueDegree()) return false;	// not equal degree
	
  if (d == -1) return true;			// zero polynomial
  // ASSERT(d >= 0)
  for (int i = 0; i <= d; i++)
    if (p.coeff()[i] != q.coeff()[i])
      return false;	// unequal 
  return true; 		// equal
}// operator==

/// non-equality comparison
template <class NT>
inline bool operator!=(const Polynomial<NT>& p, const Polynomial<NT>& q) {	// !=
  return (!(p == q));
}// operator !=


/// res(p, q)
///   -- returns the resultant of polynomial p and q.
///   Theory: Since we use "pseudoRemainder" (not true "prem")
///     we must modify the original algorithm.
///	If C * P = Q * R + S, where C is a constant and S = prem(P,Q), m=deg(P),
///	      n=deg(Q), ll = deg(S), and b=lc(Q),
///	then res(P,Q) = (-1)^(mn) b^(m-ll) res(Q,S)/C^(n).
///
///	The base case being res(P,Q) = Q^(deg(P)) if Q is a constant, zero otherwise
// Jul'10: the signs of resultant are now correctly computed
// 	(previously we didn't care about sign!
template <class NT>
inline NT res(Polynomial<NT> p, Polynomial<NT> q) {

  int m, n;
  m = p.getTrueDegree();
  n = q.getTrueDegree();

  if(m == -1 || n == -1) return 0;  // We define res(p,q)=0   if pq=0.

  if(m == 0) { 			    // by rule (i) in Chap VI.4 (Yap's bk)
    if ((n%2) == 0) return p.coeff()[0];
    else return (-1)* p.coeff()[0];
  }
  if(n == 0) { 			    // by rules (i) and (ii) in Chap VI.4
    if ((m%2) == 0) return q.coeff()[0];
    else return (-1)* q.coeff()[0];
  }

  if(n > m) {
     if ((m*n % 2)==0) return (res(q,p)); // bug fix (Narayan/Chee, Jul 7'10)
     else return (-1)*res(q,p);
  }

  NT b = q.coeff()[n];		// Leading coefficient of Q
  NT C;
  //NT lc = p.coeff()[m], C;	// Leading coefficient of P, unused!

  p.pseudoRemainder(q, C);  // p is now equal to pseudoRemainder(p,q)!

  int ll = p.getTrueDegree();

  if (((m*n)%2)==0) return power(b,m-ll)*res(q,p)/power(C,n);
  else return (-1)*(power(b,m-ll)*res(q,p)/power(C,n));
}//res

/// disc(P,Q)
///   returns the discriminant of polynomials P,Q.
///   WARNING: this is a cheap discriminant ignoring
///   sign and dividing by the leading coefficient of P.
// NOTE: to divide by lc(P), we need to implement exactdiv in all
// the NT's.  Currently, only BigInt has exactdiv.  Trivial for BigRat,
// but need it for BigFloat, DoubleWrapper, LongWrapper, IntWrapper.
template <class NT>
inline NT disc(Polynomial<NT> p) {
  Polynomial<NT> q=p; 
  return res(p, q.differentiate());  // ought to divide by lc(p), and by (-1)^s.
}//disc

///i^th Principal Subresultant Coefficient (psc) of two polynomials.
template <class NT>
inline NT psc(int i,Polynomial<NT> p, Polynomial<NT> q) {

  assert(i >= 0);
  if(i == 0)
     return res(p,q);

  int m = p.getTrueDegree();
  int n = q.getTrueDegree();

  if(m == -1 || n == -1) return 0;
  if(m < n) return psc(i, q, p);

  if ( i == n) //This case occurs when i > degree of gcd
    return power(q.coeff()[n], m - n);

  if(n < i && i <= m)
     return 0;

  NT b = q.coeff()[n];  	//Leading coefficient of Q
  NT lc = p.coeff()[m], C;	// lc unused...


  p.pseudoRemainder(q, C);

  if(zeroP(p) && i < n)//This case occurs when i < deg(gcd)
     return 0;

  if(zeroP(p) && i == n)//This case occurs when i=deg(gcd) might be constant
    return power(q.coeff()[n], m - n);

  int l = p.getTrueDegree();
  return power(NT(-1),(m-i)*(n-i))*power(b,m-l)*psc(i,q,p)/power(C, n-i);

}//psc

///factorI(p,m)
///    computes the polynomial q which containing all roots
///    of multiplicity m of a given polynomial P
///    Used to determine the nature of intersection of two algebraic curves
///    The algorithm is given in Wolperts Thesis
template <class NT>
inline Polynomial<NT> factorI(Polynomial<NT> p, int m){
  int d=p.getTrueDegree();
  Polynomial<NT> *w = new Polynomial<NT>[d+1];
  w[0] = p;
  Polynomial<NT> temp;

  for(int i = 1; i <=d ; i ++){
     temp = w[i-1];
     w[i] = gcd(w[i-1],temp.differentiate());
  }

  Polynomial<NT> *u = new Polynomial<NT>[d+1];
  u[d] = w[d-1];
  for(int i = d-1; i >=m; i--){
	temp = power(u[i+1],2);
     for(int j=i+2; j <=d; j++){
        temp *= power(u[j],j-i+1);
     }
     u[i] = w[i-1].pseudoRemainder(temp);//i-1 since the array starts at 0
  } 
	
	delete[] w;
  return u[m];
}//factorI

//////////////////////////////////////////////////
// FUNCTIONS for polynomial composition (by Philip Davidson)
//////////////////////////////////////////////////

// Univariate = Univariate( Univariate ) 

// standard composition, for P = c_i X^i

template < class NT >
Polynomial<NT> composeNaive ( const Polynomial<NT>& P, const Polynomial<NT>& X ) { 

  Polynomial<NT> P_X;
  int i, deg_P;
  NT ci;
 
  deg_P = P.getTrueDegree();
  
  for ( i = 0 ; i <= deg_P ; i++ ) { 

    ci = P.getCoeff(i);      
    if ( ci == 0 ) continue;        // null term ! requires == 0 defined for NT .. ?isZero?
    // maybe move pKi outside and reuse, if assignment clears properly

    Polynomial<NT> pKi = X;         // initialize
    pKi.power(i);                   // generate term i
    P_X += pKi.mulScalar(ci);       // accumulate

  }
  return P_X;
}


// Univariate = Univariate( Univariate ) 

// standard composition, for P = c_i X^i

template < class NT >
Polynomial<NT> composeHorner ( const Polynomial<NT>& P, const Polynomial<NT>& X ) { 

  int i, deg_P;
  NT ci;
 
  deg_P = P.getTrueDegree();
  if ( deg_P < 0 ) return Polynomial<NT>::Zero();

  ci = P.getCoeff( deg_P );
  Polynomial<NT> P_X = Polynomial<NT>( 0, &ci );
  
  for ( i = deg_P-1 ; i >= 0 ; --i ) { 

    P_X *= X;

    ci = P.getCoeff(i);
    P_X += Polynomial<NT>( 0, &ci );

  }
  return P_X;
}

// which should equal 

template < class NT >
Polynomial<NT> compose ( const Polynomial<NT>& P, const Polynomial<NT>& X ) { 
  return P(X);
}


template < class NT >
Polynomial<NT> composeHornerBinary ( const Polynomial<NT>& P, const Polynomial<NT>& A, const Polynomial<NT>& B ) { 
  
  int i, deg_P;
  NT ci;
 
  deg_P = P.getTrueDegree();

  // leading coeff of P
  ci = P.getCoeff(deg_P);
  Polynomial<NT> P_AB = Polynomial<NT>(0, &ci);
  // temporary B for powering
  Polynomial<NT> pB = B;

  for ( i = deg_P - 1 ; i >= 0 ; --i ) { 
    P_AB *= A;
    ci = P.getCoeff(i);
    P_AB += pB * Polynomial<NT>(0, &ci);
    pB *= B;
  }
  return P_AB;
}


template < class NT >
Polynomial<NT> composeBinary ( const Polynomial<NT>& P, const Polynomial<NT>& A, const Polynomial<NT>& B ) { 
  
  Polynomial<NT> P_AB;
  int i, deg_P;
  NT ci;
 
  deg_P = P.getTrueDegree();

  for ( i = 0 ; i <= deg_P ; i++ ) { 

    ci = P.getCoeff(i);      
    if ( ci == 0 ) continue;

    Polynomial<NT> pKA = A;    //initialize
    Polynomial<NT> pKB = B;
    pKA.power(i);              //generate term i
    pKB.power(deg_P-i);
    pKA *= pKB;                
    P_AB += pKA.mulScalar(ci); //accumulate
  }
  return P_AB;
}

/// Bounds
//{@

//NOTE: In most of the bounds below we assume that there is a constructor
//of the form BigFloat2(NT) for all NT's.

/// An iterative version of computing Cauchy Bound from Erich Kaltofen.
// See the writeup under collab/sep/.
template < class NT >
BigInt CauchyBound(const Polynomial<NT> &p) {
  int deg = p.getTrueDegree();
  BigInt B(1);
  BigFloat lhs(0), rhs(1);
  while (true) {
    /* compute \sum_{i=0}{deg-1}{|a_i|B^i} */
    lhs = 0;
    for (int i=deg-1; i>=0; i--) {
      lhs *= B;
      lhs += ToBigFloat2(abs(p.coeff()[i])).getLeft();//Need a lower bound on 
                                                    //the absolute value
    }
    //lhs /= abs(p.coeff()[deg]);
    //lhs.makeFloorExact();
    /* compute B^{deg} */
    if (rhs * ToBigFloat2(abs(p.coeff()[deg])).getRight() <= lhs) {
      B <<= 1;
      rhs *= (BigInt(1)<<deg);
    } else
      break;
  }
  return B;
}



///Another iterative bound which is at least as good as the above bound
///by Erich Kaltofen.
template < class NT >
BigInt UpperBound(Polynomial<NT> &p) {
  int deg = p.getTrueDegree();

  BigInt B(1);
  BigFloat lhsPos(0), lhsNeg(0), rhs(1);
  while (true) {
    /* compute \sum_{i=0}{deg-1}{|a_i|B^i} */
    lhsPos = lhsNeg = 0;
    for (int i=deg-1; i>=0; i--) {
      if (p.coeff()[i]>0) {
      	lhsPos = lhsPos * B + ToBigFloat2(p.coeff()[i]).getLeft();
      	lhsNeg = lhsNeg * B;
      } else {
      	lhsNeg = lhsNeg * B - ToBigFloat2(p.coeff()[i]).getLeft();
      	lhsPos = lhsPos * B;
      } 
    }

    /*lhsNeg /= abs(p.coeff()[deg]);
    lhsPos /= abs(p.coeff()[deg]);
    lhsPos.makeCeilExact();
    lhsNeg.makeCeilExact();*/
    //We can avoid the above steps by multiplying rhs by the leading coefficient
    //and then compare the result. Though we have to take get a BigFloat2 
    //approximation that is an upper bound on the leading coefficient.
    //
    /* compute B^{deg} */
    if (rhs * ToBigFloat2(abs(p.coeff()[deg])).getRight() <= (std::max)(lhsPos,lhsNeg)) {
      B <<= 1;
      rhs *= (BigInt(1)<<deg);
    } else
      break;
  }
  return B;
}

// Cauchy Lower Bound on Roots
// -- ASSERTION: NT is an integer type
template < class NT >
BigFloat CauchyLowerBound(const Polynomial<NT> &p) {
  if ((zeroP(p)) || p.coeff()[0] == 0)
    return BigFloat(0);
  NT mx = 0;
  int deg = p.getTrueDegree();
  for (int i = 1; i <= deg; ++i) {
    mx = core_max(mx, abs(p.coeff()[i]));
  }
  BigFloat2 e = ToBigFloat2(abs(p.coeff()[0]))/ ToBigFloat2(abs(p.coeff()[0]) + mx);
  return (div2(e.getLeft()));
}

// Separation bound for polynomials that may have multiple roots.
// We use the Rump-Schwartz bound.
//
//    ASSERT(the return value is an exact BigFloat and a Lower Bound)
//
template < class NT >
BigFloat sepBound(Polynomial<NT> &p) {
  BigInt d;
  BigFloat e;
  int deg = p.getTrueDegree();

  CORE::pow(d, BigInt(deg), ((deg)+4)/2);
  e = CORE::power(height(p).getRight()+1, deg);

  return (BigFloat2(1)/(e*2*d)).getLeft();
        // BUG fix: ``return 1/e*2*d'' was wrong
        // NOTE: the relative error in this division (1/(e*2*d))
        //   is defBFdivRelPrec (a global variable), but
        //   since this is always positive, we are OK.
        //   But to ensure that defBFdivRelPrec is used,
        //   we must make sure that e and d are exact.
        // Finally, using "makeFloorExact()" is OK because
        //   the mantissa minus error (i.e., m-err) will remain positive
        //   as long as the relative error (defBFdivRelPrec) is >1.
}

//@}

  // A Smale bound which is an \'a posteriori condition. Applying 
  // Newton iteration to any point z satisfying this condition we are 
  // sure to converge to the nearest root in a certain interval of z.
  // The condition is for all k >= 2,
  //    | \frac{p^(k)(z)}{k!p'(z)} |^{1\(k-1)} < 1/8 * |\frac{p'(z)}{p(z)}|
  // Note: code below has been streamlined (Chee)
  /*
    bool smaleBound(const Polynomial<NT> * p, BigFloat z){
    int deg = p[0].getTrueDegree();
    BigFloat max, temp, temp1, temp2;
    temp2 = p[1].eval(z);
    temp = core_abs(temp2/p[0].eval(z))/8;
    BigInt fact_k = 2;
    for(int k = 2; k <= deg; k++){
      temp1 = core_abs(p[k].eval(z)/(fact_k*temp2)); 
      if(k-1 == 2)
	temp1 = sqrt(temp1);
      else
	temp1 = root(temp1, k-1);
      if(temp1 >= temp) return false; 
    }
    return true;
    }
   */

  //An easily computable Smale's point estimate for Newton as compared to the
  //one above. The criterion is
  //
  // ||f||_{\infty} * \frac{|f(z)|}{|f'(z)|^2} 
  //                * \frac{\phi'(|z|)^2}{\phi(|z|)}  < 0.03
  // where
  //           \phi(r) = \sum_{i=0}{m}r^i,
  //           m = deg(f)
  //
  //It is given as Theorem B in [Smale86].
  //Reference:- Chapter 8 in Complexity and Real Computation, by
  //            Blum, Cucker, Shub and Smale
  //
  //For our implementation we calculate an upper bound on
  //the second fraction in the inequality above.  For r>0,
  //
  //    \phi'(r)^2     m^2 (r^m + 1)^2
  //     ---------  <  -------------------          (1)
  //    \phi(r)        (r-1) (r^{m+1} - 1)
  //
  // Alternatively, we have
  // 
  //    \phi'(r)^2     (mr^{m+1} + 1)^2
  //     ---------  <  -------------------          (2)
  //    \phi(r)        (r-1)^3 (r^{m+1} - 1)
  //
  // The first bound is better when r > 1.
  // The second bound is better when r << 1.
  // Both bounds (1) and (2) assumes r is not equal to 1.
  // When r=1, the exact value is
  //
  //    \phi'(r)^2     m^2 (m + 1)
  //     ---------  =  -----------                  (3)
  //    \phi(r)            4
  // 
  // REMARK: smaleBoundTest(z) actually computes an upper bound
  // 	on alpha(f,z), and compares it to 0.02 (then our theory
  // 	says that z is a robust approximate zero).
  //
  template <class NT>
  bool smaleBoundTest(Polynomial<NT> _poly, Polynomial<NT> _poly_derivative, const BigFloat& z){
    //assert(z.isExact());   // the bound only makes sense for exact z

#ifdef CORE_DEBUG_ROOTBOUND
    std::cout <<"Computing Smale's bound = " <<  std::endl;
#endif
    if(evalExactSign(_poly, z).sgn() == 0)// Reached the exact root.
      return true;

    BigFloat2 signEval = evalExactSign(_poly_derivative,z);
    BigFloat fprime = signEval.getLeft();
    if (signEval.sgn() < 0) fprime *= -1;
    if (fprime == 0) return false;  // z is a critical value!
    signEval = evalExactSign(_poly,z);
    BigFloat temp = signEval.getRight();
    if (signEval.sgn() < 0) temp *= -1;
    temp.div(temp, power(fprime, 2), getDefaultBFdivPrec(), BF_RNDU);
    temp = temp*height(_poly).getLeft();  // remains exact
    //Thus, temp >=  ||f||_{\infty} |\frac{f(z)}{f'(z)^2}|

    int m = _poly.getTrueDegree();    
    BigFloat x = core_abs(z);
    if (x==1)   {// special case, using (3)
            temp *= m*m*(m+1);
            temp.div_2exp(2);
	    return (temp < 0.02);
    }

    BigFloat temp1;
    if (x>1) { // use formula (1)
      //temp1 = power(m* (power(x, m)+1), 2);          // m^2*(x^m + 1)^2
      //temp1 /= ((x - 1)*(power(x, m+1) - 1));        // formula (1)
      temp1.div(power(m* (power(x, m)+1), 2), ((x - 1)*(power(x, m+1) - 1)), 
        getDefaultBFdivPrec(), BF_RNDU);        // formula (1)
    } else {  // use formula (2)
      //temp1 = power(m*(power(x, m+1) +1), 2);        // (m*x^{m+1} + 1)^2
      //temp1 /= (power(x - 1,3)*(power(x, m+1) -1));  // formula (2)
      temp1.div(power(m*(power(x, m+1) +1), 2), 
        (power(x-1,3)*(power(x, m+1)-1)), 
        getDefaultBFdivPrec(), BF_RNDU); // formula (2)
    }


    if(temp * temp1 < 0.03)          // make temp1 exact!
      return true;
    else
      return false;
  }//smaleBoundTest


  // yapsBound(p)
  // 	returns a bound on size of isolating interval of polynomial p
  // 	which is guaranteed to be in the Newton Zone.
  //    N.B. p MUST be square-free
  //
  //   Reference: Theorem 6.37, p.184 of Yap's book
  //   	   [Fundamental Problems of Algorithmic Algebra]

  template<class NT>
  BigFloat yapsBound(const Polynomial<NT> & p) {
    unsigned int deg = p.getTrueDegree();
    return  BigFloat(1)/(1 + power(BigFloat(deg), 3*deg+9)
               *power(2+height(p).getRight(),6*deg));
  }

///Norms on Polynomials
//{@
/// height function
/// @return a BigFloat with error
template < class NT >
BigFloat2 height(const Polynomial<NT> &p) {
  if (zeroP(p))
    return BigFloat2(0);
  int deg = p.getTrueDegree();
  NT ht = abs(p.coeff()[0]);

  for (int i = 1; i<= deg; i++) {
    NT absVal = abs(p.coeff()[i]);
    if (ht < absVal)
      ht = absVal;
  }

  return ToBigFloat2(ht);
}


/// length function
/// @return a BigFloat with error
template < class NT >
BigFloat2 length(Polynomial<NT> &p) {
  if (zeroP(p))
    return BigFloat2(0); //Sep'2014, Chee: fixed fatal error: "return BigFloat(0);" 
  int deg = p.getTrueDegree();
  NT length = 0;
  for (int i = 0; i< deg; i++)
    length += abs(p.coeff()[i]*p.coeff()[i]);
  return sqrt(ToBigFloat2(length));
}

template<class NT>
int signVariationofCoeff(const Polynomial<NT> p) {
  int deg = p.getTrueDegree();
  int count = 0;
  int last_sign = sign(p.coeff()[deg]);
  for (int i=deg-1; i>=0; i--) {
    if (sign(p.coeff()[i]) * last_sign < 0) {
      count++;
      last_sign *= -1;
    }
  }
  return count;
}

//@}
CORE_END_NAMESPACE

#endif /*__CORE_POLY_H__*/
