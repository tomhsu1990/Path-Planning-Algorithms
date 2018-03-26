/****************************************************************************
 * CoreAux.h -- Auxilliary functions for the Core Library
 *
 * WHAT IS IN THIS FILE?
 * 	-- low level bit manipulation routines (like ceiling log_2, etc)
 * 	-- unit test routines
 * 	-- standard error routine
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
 * $Id: CoreAux.h,v 1.31 2010/11/23 17:58:36 exact Exp $
 ***************************************************************************/
#ifndef __CORE_COREAUX_H__
#define __CORE_COREAUX_H__

#include <cmath>
#include <iostream>
#include <cstdlib>

#include <CORE/BigRat.h>

#ifdef _MSC_VER
#include <float.h>
#define ilogb(x) (int)_logb(x)
#endif
CORE_BEGIN_NAMESPACE

/// Writes out an error or warning message in the local file CORE_DIAGFILE
/** If last argument (err) is TRUE, then this is considered an error
 *  (not just warning).  In this case, the message is also printed in
 *  std::cerr, using std::perror().  We also print coreErrorMsg.
 *  */
void core_error(std::string msg, std::string file, int lineno, bool err = true);

/// Sets the message to be printed by core_error().
std::string set_error_message(std::string msg);

/// This is for debugging messages
inline void core_debug(std::string msg){
  std::cout << __FILE__ << "::" << __LINE__ << ": " << msg
            << std::endl;
}
template <class T1, class T2>
void core_test(T1 ans, T2 unknown, std::string msg1 = "", std::string msg2 = "") {
  if (ans != unknown) {
    std::cout << "ERROR!!! " << msg1 << std::endl;
    std::cout << "  answer : " << ans << std::endl;
    std::cout << "  unknown : " << unknown << std::endl;
    coretest_error = true;
  }
  else if (coretest_verbose) {
    std::cout << "CORRECT!!! " << msg2 << std::endl;
    std::cout << "  answer : " << ans << std::endl;
    std::cout << "  unknown : " << unknown << std::endl;
  }
}

// help inline functions for size_t
// 2/24/09: Willi Mann suggested replacing size_t by unsigned int since size_t is defined
// 	as unsigend int on some platform and as unsigned long on others.
//msb_t ceillg(size_t v);
msb_t ceillg(unsigned int v);

// help inline functions for long
inline sign_t sgn(long v)
{ return v==0 ? 0 : (v>0 ? 1 : -1); }
inline bool isDivisible(long x, long y)
{ return x % y == 0; }
inline long sign(long a)
{ return a==0 ? 0 : a > 0 ? 1 : -1; }
inline long abs(long x)
{ return (x>=0) ? x : (-x); }
inline long div_exact(long x, long y)
{ return x/y; }
msb_t ceillg(long v);
msb_t floorlg(long v);
inline msb_t ceilLg(long v)
{ return ceillg(v); }
inline msb_t floorLg(long v)
{ return floorlg(v); }
long gcd(long x, long y);

// help inline functions for int
inline sign_t sgn(int v)
{ return v==0 ? 0 : (v>0 ? 1 : -1); }
inline bool isDivisible(int x, int y)
{ return x % y == 0; }
inline int sign(int v)
{ return v==0 ? 0 : (v>0 ? 1 : -1); }
inline int abs(int x)
{ return (x>=0) ? x : (-x); }
inline int div_exact(int x, int y)
{ return x/y; }
inline msb_t ceillg(int v)
{ return ceillg(long(v)); }
inline msb_t floorlg(int v)
{ return floorlg(long(v)); }
inline msb_t ceilLg(int v)
{ return ceillg(v); }
inline msb_t floorLg(int v)
{ return floorlg(v); }
int gcd(int x, int y);
inline int (max)(int x, int y)
{ return x>=y ? x : y; }
inline int power(int x, int y)
{ return (int)::pow((double)x,(double)y); }

// help inline functions for unsigned long
inline sign_t sgn(unsigned long v)
{ return v==0 ? 0 : 1; }
unsigned long gcd(unsigned long x, unsigned long y);
/* size_t and unsigned long are the same in some OS and different in others.
 * When they are different, the compiler does not complain.
 * When they are the same, the compiler will complain about
 * the following functions (ceillg and floorlg) because of duplicated definitions
 * So we decide to comment this out:
 * */
// 2/24/09: remove the #if following Willi Mann's suggestion above, to remove size_t.
// #if defined (gnu) || (cyg)
msb_t ceillg(unsigned long v);
msb_t floorlg(unsigned long v);
// #endif


// help inline functions for double
inline sign_t sgn(double v)
{ return v==0 ? 0 : (v>0 ? 1 : -1); }
inline sign_t sign(double v)
{ return sgn(v); }
inline msb_t ceillg(double v)
{
 return ilogb(v)+1; 
}

inline msb_t floorlg(double v)
{ return ilogb(v); }
inline double ToDouble(double v)
{ return v; }

/// template function returns the absolute value
template <class T>
inline const T core_abs(const T& a) {
  return ((a < T(0)) ? -a : a);
}

template <class T>
inline const T core_max(const T& a, const T& b) {
  return ((a < b) ? b : a);
}

template <class T>
inline  void core_swap(T& a, T& b) {
  T tmp;
  tmp = a;
  a = b;
  b = tmp;
}

////////////////////////////////////////////////////
// Functions to check if two decimal strings are
// compatible.   Here is the definition:
//     For any decimal string SX, it has a nominal value X
//     and also an uncertainty EX.
//     E.g., SX= 1.234 then X=1.234 and EX=0.001
//     E.g., SX= 1200e-2 then X=12.00 and EX=0.01
//     E.g., SX= 12 then X=12 and EX=1
//
// The string SX then defines an interval [X-EX, X+EX].
// 
//     When we say two strings, SX and SY are compatible, we mean
//     that the corresponding intervals overlap.
//
//     To decide compatibility of SX and SY, we first compute
//     X, EX, Y, EY.  Wlog, let X <= Y.  Then
//     SX and SY are compatible iff
//           X + EX >= Y - EY.
//     REMARK: we could try to work with the stricter notion of X+EX > Y-EY
//     but that might cause trouble with some of our other implementations.
//
//
//     E.g., 1.200 is compatible with 1.19999 and 1.201
//     E.g., 1.234 is compatible with 1.235 and 1.233
//
////////////////////////////////////////////////////

/// Returns the uncertainty log_10(SX) in an output decimal string SX
inline int getUncertainty( std::string& strin) {
  int u = 0;  		// eventually we want to return u
  bool dot=false;	// have we seen decimal point yet?
  int j=0;		// position of most significant digit
  if ((strin[0] == '+')|| (strin[0] == '-')) j++;	// Takes care of sign
  while (strin[j] == '0') j++;   // Takes care of initial zeros: 00123.456

  for (size_t i = j; i < strin.size(); i++) {	
    if (strin[i] == 'e' || strin[i] == 'E') {  
      u += atoi(strin.substr(i+1,strin.size()-i).c_str());
      break;
    } 
    if (strin[i] == '.')
     dot = true;   	// found dot
    if (strin[i] >= '0' && strin[i] <= '9') {
      if (dot) u--;
    }// else error!
  }//for
  return u;
}  

// helper function for isCompatible
// this extract mantissa and exponent (and sign) from decimal string
// i.e., if value is m 10^e, then return (m, e).
//   Note: mantissa may contain a tail string of zeros:
//     e.g., 0.000 will have mantisa = 000
//           12.34000 will have mantissa of 1234000
//     In short, the length of the mantissa is the significance of the input
//     number.

inline void getDigits(std::string& strin, std::string& strout,
	int& exponent, int& sign) {
	// strin = string to be analyzed,
	// strout = string representing an integer (removing preceding 0's,
	// 	dot and exponents)
  exponent = 0;			// value of exponent
  bool dot = false;		// if decimal point is found
  int j = 0;			// position of most significant digit
  bool isZero = true;           // to detect the zero value

  sign = 1; 	//default sign
  if ((strin[0] == '+')|| (strin[0] == '-')) {
    if (strin[0] == '-') {
       sign = -1;
    }
    j++;			// Takes care of sign
  }
  while (strin[j] == '0') j++;    // Takes care of 00012.345 

  for (size_t i = j; i < strin.size(); i++) {	
    if (strin[i] == 'e' || strin[i] == 'E') {  
      exponent += atoi(strin.substr(i+1,strin.size()-i).c_str());
      break;
    } 
    if (strin[i] == '.') {
     dot = true;		// found dot
     i++;
     while (strout.size()==0 && strin[i] == '0') { // Can only remove 0's if strout is non-empty
	     i++;	// Takes care of 0.0009
	     exponent--;
     }
    }
    if (strin[i] >= '0' && strin[i] <= '9') {
      strout += strin[i]; 
      if (strin[i] != '0') isZero = false;
      if (dot) exponent--;	// takes account of decimal point 
    }// else error!
  }//for
  if (strout.size() == 0) {
     strout="0";
  }
  if (isZero) sign=0;
}  

/// Function to convert a decimal string into a BigRat
inline BigRat stringToBigRat(const char* pStrIn) {
  std::string digitIn;
  int expIn;
  int sgn;
 
  std::string strIn(pStrIn);
  getDigits(strIn, digitIn, expIn, sgn);
  BigRat XIn(digitIn);
  BigInt powIn;
  powIn.pow(10, abs(expIn));
  if (expIn > 0)
    XIn *= powIn;
  else
    XIn /= powIn;
  return(sgn * XIn);
}

inline BigRat stringToBigRat(std::string & strIn) {
  return stringToBigRat(strIn.c_str());
}

////////////////////////////////////////////////////
// Jan 2010: decided to write a general version
// of "compatible(s,t)" where we return the index i of
// the first position where strings s and t are no longer
// compatible.
//  	To a first approximation, we can now define "isCompatible(s,t)"
// to be true iff compatible(s,t) = min(length(s),length(t)).
// But this definition function does not work if either
// s or t is the value 0.  In this case, we return -1.
// Otherwise, "compatible(s,t)" returns a non-negative value.
//
// The main function for users to call is therefore
//       isCompatible(str1, str2)
// defined below.
////////////////////////////////////////////////////

inline int compatible(std::string s_1, std::string s_2) {
  std::string mantissa_1,mantissa_2;
  int exp_1,exp_2;
  int msb_1,msb_2;
  int sgn_1,sgn_2;

  getDigits(s_1,mantissa_1,exp_1, sgn_1);
  getDigits(s_2,mantissa_2,exp_2, sgn_2);

  if (sgn_1 * sgn_2 == 0) return -1;

  msb_1 = mantissa_1.length()+exp_1;
  msb_2 = mantissa_2.length()+exp_2;

  // Must not return if either input value is 0:
  if ((sgn_1 * sgn_2) < 0) return 0;
  else if (msb_1 != msb_2)  return 0;

  int minLength = (std::min)(mantissa_1.length(), mantissa_2.length());

  int i=0;
  while(mantissa_1[i]==mantissa_2[i]) {
    i++; if(i>=minLength) return i; 
  }

  switch(mantissa_1[i]-mantissa_2[i]) {
    case 1:
      i++; if(i>=minLength) return i; 
      while(mantissa_1[i]=='0' && mantissa_2[i]=='9') {
        i++; if(i>=minLength) return i; 
      }
      break;
    case -1:
      i++; if(i>=minLength) return i; 
       while(mantissa_2[i]=='0' && mantissa_1[i]=='9') {
        i++; if(i>=minLength) return i; 
      }
      break;
  }
  return i;
}

/// Helper function to check if two decimal string values are compatible
/// Every decimal string has an implicit uncertainty, and 2 strings are
/// compatible if they could be equal within this uncertainty.

inline bool isCompatibleEfficient(std::string& strIn, std::string& strAns) {
  int uIn=getUncertainty(strIn);
  int uAns=getUncertainty(strAns);
  
  BigInt temp;
  temp.pow(10,abs(uIn));
  BigRat uncertaintyIn(temp);
  temp.pow(10,abs(uAns));
  BigRat uncertaintyAns(temp);
  if (uIn < 0)
    uncertaintyIn.inv();  
  if (uAns < 0)
    uncertaintyAns.inv();  

  BigRat XIn=stringToBigRat(strIn);
  BigRat XAns=stringToBigRat(strAns);

  if (XIn >= XAns) 
	  return (XIn - uncertaintyIn < XAns + uncertaintyAns);
  else
	  return (XIn + uncertaintyIn > XAns - uncertaintyAns);
}//isCompatible



/// isCompatible(s1,s2) is the MAIN function for users to call.
///   It calls
///      compatible(s1,s2) and  isCompatibleEfficient(s1,s2)
///   as helper functions.
///
inline bool isCompatible(std::string s_1,std::string s_2) {
  std::string mantissa_1,mantissa_2;
  int exp_1,exp_2;
  int sgn_1,sgn_2;

  getDigits(s_1,mantissa_1,exp_1,sgn_1);
  getDigits(s_2,mantissa_2,exp_2,sgn_2);

  int minLength = (std::min)(mantissa_1.length(), mantissa_2.length());

  if (sgn_1*sgn_2 == 0) return isCompatibleEfficient(s_1, s_2);
	//:in case sgn_1*sgn_2=0, we cannot use "compatible" information...
  if (sgn_1*sgn_2 < 0) return(false);	// different signs are not compatible
  return (compatible(s_1,s_2) == minLength); 
}


#ifdef CORE_OLDNAMES 
/// \addtogroup GlobalBackCompatiableFunctions
//@{
inline long ceilLg(long a) { return ceilLg(BigInt(a)); }
inline long ceilLg(int a) { return ceilLg(BigInt(a)); }
//@}
#endif

CORE_END_NAMESPACE

#endif /*__CORE_COREAUX_H__*/
