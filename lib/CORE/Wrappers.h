/*
 * Wrappers.h 
 * 	Wrappers for fundamental machine types that implement
 *	the same interface as provided by BigFloat, BigInt, BigRat
 *	and Expr classes.  This allows a Level 2 or Level 3 program
 *	to be compilable in Level 1.
 *
 * ***WARNING*** : This implementation is incomplete, and is still 
 * experimental. The section that is complete is correct, and has been
 * built with the explicit purpose of making the CXY (progs/mesh/) code
 * work correctly. The missing section of the BigFloat / BigInt API can
 * be added to these classes 
 *
 * Author: Narayan Kamath
 * Since Core 2.0
 */
#ifndef __CORE_WRAPPERS_H__
#define __CORE_WRAPPERS_H__

#include <CORE/Mpfr.h>  // for the definition of prec_t
#include <cmath>
#include <stdlib.h>
#include <sstream>

CORE_BEGIN_NAMESPACE

/* ***************************************************
 * Wrapper class for machine long type.
 * It is currently incomplete
 * apart from an implementation of pow, and is disabled
 * by default. Will be implemented if a real need
 * arises.
 *
 * NOTE: The compiler generated copy constructor and default constructor
 * have not been suppressed.
 *************************************************** */
class LongWrapper {
 public:
  long l_val;
  operator long() { return l_val; }
};

/// power(l, p)
///   --power function for LongWrapper
//
// NOTE (Narayan): This function is a magnet for unwanted automatic
// type conversions.
//
// NOTE (Narayan): If this function was not defined "inline" then many compilers
// would break. This is because CORE.h will be compiled into every
// translation unit (normally a C++ file that gets converted to a .o file)
// and each of those files will contain a definition of this function,
// and when multiple such files are linked together, there will be
// an ambiguous definition of the said function as the exact same
// symbol occurs in multiple object files.
//
// This is a problem because the "inline" keyword is a request
// to the compiler, and not an order. Future compiler implementations have
// the potential to break the library.
//
inline LongWrapper power(const LongWrapper &l, unsigned long p) {
  LongWrapper ret;
  ret.l_val = (long)::pow((double)l.l_val, (double)p);  // calls pow of std math library
  return ret;
}


/* ***************************************************
 * END of LongWrapper Definitions
 *************************************************** */


/* ***************************************************
 * Wrapper class for machine int type.
 * It is currently incomplete
 * apart from an implementation of pow. 
 *
 * NOTE: The compiler generated copy constructor and default constructor
 * have not been suppressed.
 *************************************************** */

class IntWrapper {
 public:
  int i_val;
  operator int() { return i_val; }

  // Chee (July 2012): basic datatype conversions,
  // 		copied from DoubleWrapper (needed by tBiPoly.cpp):
  //		==================================================
  // Conversion from fundamental data types, the default C++
  // promotions / demotions occur.
  IntWrapper() : i_val(0) { }
  IntWrapper(const double &rhs) : i_val(rhs) { }
  IntWrapper(const int rhs) : i_val(rhs) { }
  IntWrapper(const long rhs) : i_val(rhs) { }
  IntWrapper(const LongWrapper &w) : i_val(w.l_val) { }
  // 		IntWrapper(const DoubleWrapper &d) : i_val(d.d_val) { }
  // Duplication of Promote.h: (?)
  IntWrapper(const char* str) : i_val(atoi(str)) {}

  // Chee (July 2012): *= operator,
  // 		copied from DoubleWrapper (needed by tBiPoly.cpp):
  //		==================================================
  /// compound assignment operator <tt>*=</tt>
  //		IntWrapper& operator*=(const IntWrapper& rhs)
  //		{  i_val *= rhs.d_val; return *this; }
  /// compound assignment operator <tt>*=</tt>
  IntWrapper& operator*=(const int rhs)
  {  i_val *= rhs; return *this; }
  /// compound assignment operator <tt>*=</tt>
  IntWrapper& operator*=(const unsigned int rhs)
  {  i_val *= rhs; return *this;  }
  /// compound assignment operator <tt>*=</tt>
  IntWrapper& operator*=(const long rhs)
  {  i_val *= rhs; return *this;  }
  /// compound assignment operator <tt>*=</tt>
  IntWrapper& operator*=(const unsigned long rhs)
  {  i_val *= rhs; return *this;  }
  /// compound assignment operator <tt>*=</tt>
  IntWrapper& operator*=(const double rhs)
  {  i_val *= rhs; return *this;  }
  
  // Chee (July 2012): += operator,
  // 		modeled after *= operator of IntWrapper above (needed by Poly.h):
  //		==================================================
  /// compound assignment operator <tt>+=</tt>
  //		IntWrapper& operator+=(const IntWrapper& rhs)
  //		{  i_val += rhs.d_val; return *this; }
  /// compound assignment operator <tt>+=</tt>
  IntWrapper& operator+=(const int rhs)
  {  i_val += rhs; return *this; }
  /// compound assignment operator <tt>+=</tt>
  IntWrapper& operator+=(const unsigned int rhs)
  {  i_val += rhs; return *this;  }
  /// compound assignment operator <tt>+=</tt>
  IntWrapper& operator+=(const long rhs)
  {  i_val += rhs; return *this;  }
  /// compound assignment operator <tt>+=</tt>
  IntWrapper& operator+=(const unsigned long rhs)
  {  i_val += rhs; return *this;  }
  /// compound assignment operator <tt>+=</tt>
  IntWrapper& operator+=(const double rhs)
  {  i_val += rhs; return *this;  }

  // Chee (July 2012): -= operator,
  // 		modeled after += operator of IntWrapper above (needed by Poly.h):
  //		==================================================
  /// compound assignment operator <tt>-=</tt>
  //		IntWrapper& operator-=(const IntWrapper& rhs)
  //		{  i_val -= rhs.d_val; return *this; }
  /// compound assignment operator <tt>-=</tt>
  IntWrapper& operator-=(const int rhs)
  {  i_val -= rhs; return *this; }
  /// compound assignment operator <tt>-=</tt>
  IntWrapper& operator-=(const unsigned int rhs)
  {  i_val -= rhs; return *this;  }
  /// compound assignment operator <tt>-=</tt>
  IntWrapper& operator-=(const long rhs)
  {  i_val -= rhs; return *this;  }
  /// compound assignment operator <tt>-=</tt>
  IntWrapper& operator-=(const unsigned long rhs)
  {  i_val -= rhs; return *this;  }
  /// compound assignment operator <tt>-=</tt>
  IntWrapper& operator-=(const double rhs)
  {  i_val -= rhs; return *this;  }

};//IntWrapper Class

/// power(i,p)
///    --power function for IntWrapper
inline IntWrapper power(const IntWrapper &i, unsigned long p) {
  IntWrapper ret;
  ret.i_val = (int)::pow((double)i.i_val, (double)p);  // calls pow of std math library
  return ret;
}

  // Chee (July 2012): != operator,
  // 		modeled after == operator for IntWrapper (needed by inc/CORE/poly/Poly.h):
  //		==================================================
  /// \addtogroup IntWrapperComparisonOperators
  //@{
  /// IntWrapper  != IntWrapper
  inline bool operator!=(const IntWrapper & x, const IntWrapper & y)
  { return x.i_val != y.i_val; }
  /// IntWrapper  != int
  inline bool operator!=(const IntWrapper & x, int y)
  { return x.i_val != y; }
  /// int != IntWrapper
  inline bool operator!=(int x, const IntWrapper & y)
  { return x != y.i_val; }
  /// IntWrapper  != unsigned int	(WARNING...)
  //		inline bool operator!=(const IntWrapper & x, unsigned int y)
  //		{ return x.i_val != y; }
  /// unsigned int != IntWrapper	(WARNING...)
  //		inline bool operator!=(unsigned int x, const IntWrapper & y)
  //		{ return x != y.i_val;; }
  /// IntWrapper  != long
  inline bool operator!=(const IntWrapper & x, long y)
  { return x.i_val != y; }
  /// long != IntWrapper
  inline bool operator!=(long x, const IntWrapper & y)
  { return x != y.i_val; }
  /// IntWrapper  != unsigned long	(WARNING...)
  //		inline bool operator!=(const IntWrapper & x, unsigned long y)
  //		{ return x.i_val != y; }
  /// unsigned long != IntWrapper	(WARNING...)
  //		inline bool operator!=(unsigned long x, const IntWrapper & y)
  //		{ return x != y.i_val;; }
  /// DoubleWrapper  != double
  //		inline bool operator!=(const DoubleWrapper & x, double y)
  //		{ return x.d_val != y; }
  /// double != DoubleWrapper
  //		inline bool operator!=(double x, const DoubleWrapper & y)
  //		{ return x != y.d_val; }
  /// IntWrapper  != LongWrapper
  inline bool operator!=(const IntWrapper & x, const LongWrapper& y)
  { return x.i_val != y.l_val; }
  /// LongWrapper != IntWrapper
  inline bool operator!=(const LongWrapper& x, const IntWrapper & y)
  { return x.l_val != y.i_val; }

  //@}
 
  // Chee (July 2012): == operator,
  // 		copied from DoubleWrapper (needed by tBiPoly.cpp):
  //		==================================================
  /// \addtogroup IntWrapperComparisonOperators
  //@{
  /// IntWrapper  == IntWrapper
  inline bool operator==(const IntWrapper & x, const IntWrapper & y)
  { return x.i_val == y.i_val; }
  /// IntWrapper  == int
  inline bool operator==(const IntWrapper & x, int y)
  { return x.i_val == y; }
  /// int == IntWrapper
  inline bool operator==(int x, const IntWrapper & y)
  { return x == y.i_val; }
  /// IntWrapper  == unsigned int	(WARNING...)
  //		inline bool operator==(const IntWrapper & x, unsigned int y)
  //		{ return x.i_val == y; }
  /// unsigned int == IntWrapper	(WARNING...)
  //		inline bool operator==(unsigned int x, const IntWrapper & y)
  //		{ return x == y.i_val;; }
  /// IntWrapper  == long
  inline bool operator==(const IntWrapper & x, long y)
  { return x.i_val == y; }
  /// long == IntWrapper
  inline bool operator==(long x, const IntWrapper & y)
  { return x == y.i_val; }
  /// IntWrapper  == unsigned long	(WARNING...)
  //		inline bool operator==(const IntWrapper & x, unsigned long y)
  //		{ return x.i_val == y; }
  /// unsigned long == IntWrapper	(WARNING...)
  //		inline bool operator==(unsigned long x, const IntWrapper & y)
  //		{ return x == y.i_val;; }
  /// DoubleWrapper  == double
  //		inline bool operator==(const DoubleWrapper & x, double y)
  //		{ return x.d_val == y; }
  /// double == DoubleWrapper
  //		inline bool operator==(double x, const DoubleWrapper & y)
  //		{ return x == y.d_val; }
  /// IntWrapper  == LongWrapper
  inline bool operator==(const IntWrapper & x, const LongWrapper& y)
  { return x.i_val == y.l_val; }
  /// LongWrapper == IntWrapper
  inline bool operator==(const LongWrapper& x, const IntWrapper & y)
  { return x.l_val == y.i_val; }

  //@}
  
/* */
// Chee (July 2012): modeling after DoubleWrapper ArithmeticOperators
//
/// \addtogroup IntWrapper ArithmeticOperators
//@{
/// IntWrapper * IntWrapper

/// IntWrapper * IntWrapper
inline IntWrapper operator*(const IntWrapper& x, const IntWrapper& y)
{ IntWrapper r; r.i_val = x.i_val * y.i_val;  return r; }
/// IntWrapper * int
inline IntWrapper operator*(const IntWrapper& x, const int y)
{ IntWrapper r; r.i_val = x.i_val * y;  return r; }
/// int * IntWrapper
inline IntWrapper operator*(const int x, const IntWrapper& y)
{ IntWrapper r; r.i_val = x * y.i_val;  return r; }
/// IntWrapper * unsigned int
inline IntWrapper operator*(const IntWrapper& x, const unsigned int y)
{ IntWrapper r; r.i_val = x.i_val * y;  return r; }
/// unsigned int * IntWrapper
inline IntWrapper operator*(const unsigned int x, const IntWrapper& y)
{ IntWrapper r; r.i_val = x * y.i_val; return r; }
/// IntWrapper * long
inline IntWrapper operator*(const IntWrapper& x, const long y)
{ IntWrapper r;  r.i_val = x.i_val * y; return r; }
/// long * IntWrapper
inline IntWrapper operator*(const long x, const IntWrapper& y)
{ IntWrapper r; r.i_val = x * y.i_val;  return r; }
/// IntWrapper * unsigned long
inline IntWrapper operator*(const IntWrapper& x, const unsigned long y)
{ IntWrapper r; r.i_val = x.i_val * y;  return r; }
/// unsigned long * IntWrapper
inline IntWrapper operator*(const unsigned long x, const IntWrapper& y)
{ IntWrapper r; r.i_val = x * y.i_val;  return r; }
/// IntWrapper * double
inline IntWrapper operator*(const IntWrapper& x, const double y)
{ IntWrapper r; r.i_val = x.i_val * y;  return r; }
/// double * IntWrapper
inline IntWrapper operator*(const double x, const IntWrapper& y)
{ IntWrapper r; r.i_val = x * y.i_val;  return r; }
/// IntWrapper * LongWrapper
inline IntWrapper operator*(const IntWrapper& x, const LongWrapper& y)
{ IntWrapper r; r.i_val = x.i_val * y.l_val;  return r; }
/// LongWrapper * IntWrapper
inline IntWrapper operator*(const LongWrapper& x, const IntWrapper& y)
{ IntWrapper r; r.i_val = x.l_val * y.i_val; return r; }
 
//@}

/* */

/* ***************************************************
 * END of IntWrapper Definitions
 *************************************************** */

/* ***************************************************
 * Wrapper class for machine double type.
 * It is currently incomplete, but is expected to be filled as needed.
 * But it is rather more complete
 * than the LongWrapper or IntWrapper classes.
 *
 * NOTE: The compiler generated copy constructor and default constructor
 * have not been suppressed.
 *************************************************** */
class DoubleWrapper {
 public:
  // The only data member of this class, deliberately defined
  // as public, as there is no real advantage to defining it as
  // private, as there are public functions that offer unrestricted
  // access to it anyway.
  double d_val;

  // Conversion from the corresponding string representations.
  DoubleWrapper(const char *string_rep) {
    std::istringstream ss(string_rep);
    ss >> d_val;
  }
  DoubleWrapper(const std::string &string_rep) {
    std::istringstream ss(string_rep);
    ss >> d_val;
  }

  // Conversion from fundamental data types, the default C++
  // promotions / demotions occur.
  DoubleWrapper() : d_val(0.0f) { }
  DoubleWrapper(const double &rhs) : d_val(rhs) { }
  DoubleWrapper(const int rhs) : d_val(rhs) { }
  DoubleWrapper(const long rhs) : d_val(rhs) { }
  DoubleWrapper(const LongWrapper &w) : d_val(w.l_val) { }
  DoubleWrapper(const DoubleWrapper &d) : d_val(d.d_val) { }

  // The following two functions belong to the old
  // "API but turn out to be called most frequently.
  double doubleValue() const {
    return d_val;
  }
  int intValue() const {
    return static_cast<int>(d_val);
  }

  // TODO : fix this implementation.
  long lMSB() const {
    return 0L;
  }

  static const prec_t UNUSED_PREC = 0;
  static const rnd_t UNUSED_RND = (rnd_t)0;
  /// To handle calls to "approx" for Expressions.
  DoubleWrapper approx(prec_t arg1=UNUSED_PREC, prec_t arg2=UNUSED_PREC) {
    return (*this);
  }

  void pi(prec_t prec = UNUSED_PREC, rnd_t rnd = UNUSED_RND) {
    // This might seem like a bit of a hack, this is taken from the GNU
    // C library header <cmath> . But its not a part of the standard, so
    // we cannot expect it to be present on all systems (and indeed it is
    // not on Mac OS X) .
    d_val = 3.14159265358979323846;
  }

  /// unary plus operator
  DoubleWrapper operator+() const
  { return *this; }
  /// unary negation operator
  DoubleWrapper operator-() const
  { DoubleWrapper r(-(this->d_val)); return r; }
  /// prefix increment operator
  DoubleWrapper& operator++()
  { d_val++; return *this; }
  /// prefix decrement operator
  DoubleWrapper& operator--()
  { d_val--; return *this; }

  /// \name assignment and compound assignment operators
  //@{
  /// assignment operator for <tt>BigFloat</tt>
  DoubleWrapper& operator=(const DoubleWrapper& rhs)
  {  d_val = rhs.d_val; return *this; }
  /// assignment operator for <tt>int</tt>
  DoubleWrapper& operator=(int rhs)
  {  d_val = rhs; return *this; }
  /// assignment operator for <tt>unsigned int</tt>
  DoubleWrapper& operator=(unsigned int rhs)
  {  d_val = rhs; return *this; }
  /// assignment operator for <tt>long</tt>
  DoubleWrapper& operator=(long rhs)
  {  d_val = rhs; return *this; }
  /// assignment operator for <tt>unsigned long</tt>
  DoubleWrapper& operator=(unsigned long rhs)
  {  d_val = rhs; return *this; }
  /// assignment operator for <tt>double</tt>
  DoubleWrapper& operator=(double rhs)
  {  d_val = rhs; return *this; }
  /// assignment operator for <tt>char*</tt>
  DoubleWrapper& operator=(const char* rhs)
  {  std::istringstream ss(rhs); ss >> d_val; return *this; }
  /// assignment operator for <tt>std::string</tt>
  DoubleWrapper& operator=(const std::string& rhs)
  {  std::istringstream ss(rhs); ss >> d_val; return *this; }

  /// compound assignment operator <tt>+=</tt>
  DoubleWrapper& operator+=(const DoubleWrapper& rhs)
  {  d_val += rhs.d_val; return *this; }
  /// compound assignment operator <tt>+=</tt>
  DoubleWrapper& operator+=(int rhs)
  {  d_val += rhs; return *this; }
  /// compound assignment operator <tt>+=</tt>
  DoubleWrapper& operator+=(unsigned int rhs)
  {  d_val += rhs; return *this; }
  /// compound assignment operator <tt>+=</tt>
  DoubleWrapper& operator+=(long rhs)
  {  d_val += rhs; return *this; }
  /// compound assignment operator <tt>+=</tt>
  DoubleWrapper& operator+=(unsigned long rhs)
  {  d_val += rhs; return *this; }
  /// compound assignment operator <tt>+=</tt>
  DoubleWrapper& operator+=(double rhs)
  {  d_val += rhs; return *this; }

  /// compound assignment operator <tt>-=</tt>
  DoubleWrapper& operator-=(const DoubleWrapper& rhs)
  {  d_val -= rhs.d_val; return *this; }
  /// compound assignment operator <tt>-=</tt>
  DoubleWrapper& operator-=(const int rhs)
  {  d_val -= rhs; return *this; }
  /// compound assignment operator <tt>-=</tt>
  DoubleWrapper& operator-=(const unsigned int rhs)
  {  d_val -= rhs; return *this; }
  /// compound assignment operator <tt>-=</tt>
  DoubleWrapper& operator-=(const long rhs)
  {  d_val -= rhs; return *this; }
  /// compound assignment operator <tt>-=</tt>
  DoubleWrapper& operator-=(const unsigned long rhs)
  {  d_val -= rhs; return *this; }
  /// compound assignment operator <tt>-=</tt>
  DoubleWrapper& operator-=(const double rhs)
  {  d_val -= rhs; return *this; }

 /// compound assignment operator <tt>*=</tt>
 DoubleWrapper& operator*=(const DoubleWrapper& rhs)
 {  d_val *= rhs.d_val; return *this; }
 /// compound assignment operator <tt>*=</tt>
 DoubleWrapper& operator*=(const int rhs)
 {  d_val *= rhs; return *this; }
 /// compound assignment operator <tt>*=</tt>
 DoubleWrapper& operator*=(const unsigned int rhs)
 {  d_val *= rhs; return *this;  }
 /// compound assignment operator <tt>*=</tt>
 DoubleWrapper& operator*=(const long rhs)
 {  d_val *= rhs; return *this;  }
 /// compound assignment operator <tt>*=</tt>
 DoubleWrapper& operator*=(const unsigned long rhs)
 {  d_val *= rhs; return *this;  }
 /// compound assignment operator <tt>*=</tt>
 DoubleWrapper& operator*=(const double rhs)
 {  d_val *= rhs; return *this;  }

 /// compound assignment operator <tt>/=</tt>
 DoubleWrapper& operator/=(const DoubleWrapper& rhs)
 {  d_val /= rhs.d_val; return *this; }
 /// compound assignment operator <tt>/=</tt>
 DoubleWrapper& operator/=(const int rhs)
 {  d_val /= rhs; return *this; }
 /// compound assignment operator <tt>/=</tt>
 DoubleWrapper& operator/=(const unsigned int rhs)
 {  d_val /= rhs; return *this; }
 /// compound assignment operator <tt>/=</tt>
 DoubleWrapper& operator/=(const long rhs)
 {  d_val /= rhs; return *this; }
 /// compound assignment operator <tt>/=</tt>
 DoubleWrapper& operator/=(const unsigned long rhs)
 {  d_val /= rhs; return *this; }
 /// compound assignment operator <tt>/=</tt>
 DoubleWrapper& operator/=(const double rhs)
 {  d_val /= rhs; return *this; }

 //@}
};//DoubleWrapper Class


/// \addtogroup BigFloatArithmeticOperators
//@{
/// BigFloat + BigFloat

inline DoubleWrapper operator+(const DoubleWrapper& x, const DoubleWrapper& y)
{ DoubleWrapper r; r.d_val = x.d_val + y.d_val; return r; }
/// DoubleWrapper + int
inline DoubleWrapper operator+(const DoubleWrapper& x, const  int y)
{ DoubleWrapper r; r.d_val = x.d_val + y; return r; }
/// int + DoubleWrapper
inline DoubleWrapper operator+(const int x, const DoubleWrapper& y)
{ DoubleWrapper r; r.d_val = x + y.d_val; return r; }
/// DoubleWrapper + unsigned int
inline DoubleWrapper operator+(const DoubleWrapper& x, const unsigned int y)
{ DoubleWrapper r; r.d_val = x.d_val + y; return r; }
/// unsigned int + DoubleWrapper
inline DoubleWrapper operator+(const unsigned int x, const DoubleWrapper& y)
{ DoubleWrapper r; r.d_val = x + y.d_val;  return r; }
/// DoubleWrapper + long
inline DoubleWrapper operator+(const DoubleWrapper& x, const long y)
{ DoubleWrapper r; r.d_val = x.d_val + y;  return r; }
/// long + DoubleWrapper
inline DoubleWrapper operator+(const long x, const DoubleWrapper& y)
{ DoubleWrapper r; r.d_val = x + y.d_val;  return r; }
/// DoubleWrapper + unsigned long
inline DoubleWrapper operator+(const DoubleWrapper& x, const unsigned long y)
{ DoubleWrapper r; r.d_val = x.d_val + y;  return r; }
/// unsigned long + DoubleWrapper
inline DoubleWrapper operator+(const unsigned long x, const DoubleWrapper& y)
{ DoubleWrapper r; r.d_val = x + y.d_val;  return r; }
/// DoubleWrapper + double
inline DoubleWrapper operator+(const DoubleWrapper& x, const double y)
{ DoubleWrapper r; r.d_val = x.d_val + y;  return r; }
/// double + DoubleWrapper
inline DoubleWrapper operator+(const double x, const DoubleWrapper& y)
{ DoubleWrapper r; r.d_val = x + y.d_val;  return r; }
/// DoubleWrapper + LongWrapper
inline DoubleWrapper operator+(const DoubleWrapper& x, const LongWrapper& y)
{ DoubleWrapper r; r.d_val = x.d_val + y.l_val;  return r; }
/// LongWrapper + DoubleWrapper
inline DoubleWrapper operator+(const LongWrapper& x, const DoubleWrapper& y)
{ DoubleWrapper r; r.d_val = x.l_val + y.d_val;  return r; }

/// DoubleWrapper - DoubleWrapper
inline DoubleWrapper operator-(const DoubleWrapper& x, const DoubleWrapper& y)
{ DoubleWrapper r; r.d_val = x.d_val - y.d_val; return r; }
/// DoubleWrapper - int
inline DoubleWrapper operator-(const DoubleWrapper& x, const int y)
{ DoubleWrapper r; r.d_val = x.d_val - y;  return r; }
/// int - DoubleWrapper
inline DoubleWrapper operator-(const int x, const DoubleWrapper& y)
{ DoubleWrapper r; r.d_val = x - y.d_val;  return r; }
/// DoubleWrapper - unsigned int
inline DoubleWrapper operator-(const DoubleWrapper& x, const unsigned int y)
{ DoubleWrapper r; r.d_val = x.d_val - y;  return r; }
/// unsigned int - DoubleWrapper
inline DoubleWrapper operator-(const unsigned int x, const DoubleWrapper& y)
{ DoubleWrapper r; r.d_val = x - y.d_val; return r; }
/// DoubleWrapper - long
inline DoubleWrapper operator-(const DoubleWrapper& x, const long y)
{ DoubleWrapper r; r.d_val = x.d_val - y;  return r; }
/// long - DoubleWrapper
inline DoubleWrapper operator-(const long x, const DoubleWrapper& y)
{ DoubleWrapper r; r.d_val = x - y.d_val; return r; }
/// DoubleWrapper - unsigned long
inline DoubleWrapper operator-(const DoubleWrapper& x, const unsigned long y)
{ DoubleWrapper r; r.d_val = x.d_val - y; return r; }
/// unsigned long - DoubleWrapper
inline DoubleWrapper operator-(const unsigned long x, const DoubleWrapper& y)
{ DoubleWrapper r; r.d_val = x - y.d_val; return r; }
/// DoubleWrapper - double
inline DoubleWrapper operator-(const DoubleWrapper& x, const double y)
{ DoubleWrapper r; r.d_val = x.d_val - y; return r; }
/// double - DoubleWrapper
inline DoubleWrapper operator-(const double x, const DoubleWrapper& y)
{ DoubleWrapper r; r.d_val = x - y.d_val; return r; }
/// DoubleWrapper - LongWrapper
inline DoubleWrapper operator-(const DoubleWrapper& x, const LongWrapper& y)
{ DoubleWrapper r; r.d_val = x.d_val - y.l_val; return r; }
/// LongWrapper - DoubleWrapper
inline DoubleWrapper operator-(const LongWrapper& x, const DoubleWrapper& y)
{ DoubleWrapper r; r.d_val = x.l_val - y.d_val;  return r; }

/// DoubleWrapper * DoubleWrapper
inline DoubleWrapper operator*(const DoubleWrapper& x, const DoubleWrapper& y)
{ DoubleWrapper r; r.d_val = x.d_val * y.d_val;  return r; }
/// DoubleWrapper * int
inline DoubleWrapper operator*(const DoubleWrapper& x, const int y)
{ DoubleWrapper r; r.d_val = x.d_val * y;  return r; }
/// int * DoubleWrapper
inline DoubleWrapper operator*(const int x, const DoubleWrapper& y)
{ DoubleWrapper r; r.d_val = x * y.d_val;  return r; }
/// DoubleWrapper * unsigned int
inline DoubleWrapper operator*(const DoubleWrapper& x, const unsigned int y)
{ DoubleWrapper r; r.d_val = x.d_val * y;  return r; }
/// unsigned int * DoubleWrapper
inline DoubleWrapper operator*(const unsigned int x, const DoubleWrapper& y)
{ DoubleWrapper r; r.d_val = x * y.d_val; return r; }
/// DoubleWrapper * long
inline DoubleWrapper operator*(const DoubleWrapper& x, const long y)
{ DoubleWrapper r;  r.d_val = x.d_val * y; return r; }
/// long * DoubleWrapper
inline DoubleWrapper operator*(const long x, const DoubleWrapper& y)
{ DoubleWrapper r; r.d_val = x * y.d_val;  return r; }
/// DoubleWrapper * unsigned long
inline DoubleWrapper operator*(const DoubleWrapper& x, const unsigned long y)
{ DoubleWrapper r; r.d_val = x.d_val * y;  return r; }
/// unsigned long * DoubleWrapper
inline DoubleWrapper operator*(const unsigned long x, const DoubleWrapper& y)
{ DoubleWrapper r; r.d_val = x * y.d_val;  return r; }
/// DoubleWrapper * double
inline DoubleWrapper operator*(const DoubleWrapper& x, const double y)
{ DoubleWrapper r; r.d_val = x.d_val * y;  return r; }
/// double * DoubleWrapper
inline DoubleWrapper operator*(const double x, const DoubleWrapper& y)
{ DoubleWrapper r; r.d_val = x * y.d_val;  return r; }
/// DoubleWrapper * LongWrapper
inline DoubleWrapper operator*(const DoubleWrapper& x, const LongWrapper& y)
{ DoubleWrapper r; r.d_val = x.d_val * y.l_val;  return r; }
/// LongWrapper * DoubleWrapper
inline DoubleWrapper operator*(const LongWrapper& x, const DoubleWrapper& y)
{ DoubleWrapper r; r.d_val = x.l_val * y.d_val; return r; }

/// DoubleWrapper / DoubleWrapper
inline DoubleWrapper operator/(const DoubleWrapper& x, const DoubleWrapper& y)
{ DoubleWrapper r; r.d_val = x.d_val / y.d_val;  return r; }
/// DoubleWrapper / int
inline DoubleWrapper operator/(const DoubleWrapper& x, int y)
{ DoubleWrapper r; r.d_val = x.d_val / y; return r; }
/// int / DoubleWrapper
inline DoubleWrapper operator/(int x, const DoubleWrapper& y)
{ DoubleWrapper r; r.d_val = x / y.d_val; return r; }
/// DoubleWrapper / unsigned int
inline DoubleWrapper operator/(const DoubleWrapper& x, unsigned int y)
{ DoubleWrapper r; r.d_val = x.d_val / y;  return r; }
/// unsigned int / DoubleWrapper
inline DoubleWrapper operator/(unsigned int x, const DoubleWrapper& y)
{ DoubleWrapper r; r.d_val = x / y.d_val;  return r; }
/// DoubleWrapper / long
inline DoubleWrapper operator/(const DoubleWrapper& x, long y)
{ DoubleWrapper r; r.d_val = x.d_val / y;  return r; }
/// long / DoubleWrapper
inline DoubleWrapper operator/(long x, const DoubleWrapper& y)
{ DoubleWrapper r; r.d_val = x / y.d_val;  return r; }
/// DoubleWrapper / unsigned long
inline DoubleWrapper operator/(const DoubleWrapper& x, unsigned long y)
{ DoubleWrapper r; r.d_val = x.d_val / y;  return r; }
/// unsigned long / DoubleWrapper
inline DoubleWrapper operator/(unsigned long x, const DoubleWrapper& y)
{ DoubleWrapper r; r.d_val = x / y.d_val;  return r; }
/// DoubleWrapper / double
inline DoubleWrapper operator/(const DoubleWrapper& x, double y)
{ DoubleWrapper r; r.d_val = x.d_val / y;  return r; }
/// double / DoubleWrapper
inline DoubleWrapper operator/(double x, const DoubleWrapper& y)
{ DoubleWrapper r; r.d_val = x / y.d_val;  return r; }
/// DoubleWrapper / LongWrapper
inline DoubleWrapper operator/(const DoubleWrapper& x, const LongWrapper& y)
{ DoubleWrapper r; r.d_val = x.d_val / y.l_val;  return r; }
/// LongWrapper / DoubleWrapper
inline DoubleWrapper operator/(const LongWrapper& x, const DoubleWrapper& y)
{ DoubleWrapper r; r.d_val = x.l_val / y.d_val;  return r; }
//@}

/// \addtogroup DoubleWrapperComparisonOperators
//@{
/// DoubleWrapper  == DoubleWrapper
inline bool operator==(const DoubleWrapper & x, const DoubleWrapper & y)
{ return x.d_val == y.d_val; }
/// DoubleWrapper  == int
inline bool operator==(const DoubleWrapper & x, int y)
{ return x.d_val == y; }
/// int == DoubleWrapper
inline bool operator==(int x, const DoubleWrapper & y)
{ return x == y.d_val; }
/// DoubleWrapper  == unsigned int
inline bool operator==(const DoubleWrapper & x, unsigned int y)
{ return x.d_val == y; }
/// unsigned int == DoubleWrapper
inline bool operator==(unsigned int x, const DoubleWrapper & y)
{ return x == y.d_val;; }
/// DoubleWrapper  == long
inline bool operator==(const DoubleWrapper & x, long y)
{ return x.d_val == y; }
/// long == DoubleWrapper
inline bool operator==(long x, const DoubleWrapper & y)
{ return x == y.d_val;; }
/// DoubleWrapper  == unsigned long
inline bool operator==(const DoubleWrapper & x, unsigned long y)
{ return x.d_val == y; }
/// unsigned long == DoubleWrapper
inline bool operator==(unsigned long x, const DoubleWrapper & y)
{ return x == y.d_val;; }
/// DoubleWrapper  == double
inline bool operator==(const DoubleWrapper & x, double y)
{ return x.d_val == y; }
/// double == DoubleWrapper
inline bool operator==(double x, const DoubleWrapper & y)
{ return x == y.d_val; }
/// DoubleWrapper  == LongWrapper
inline bool operator==(const DoubleWrapper & x, const LongWrapper& y)
{ return x.d_val == y.l_val; }
/// LongWrapper == DoubleWrapper
inline bool operator==(const LongWrapper& x, const DoubleWrapper & y)
{ return x.l_val == y.d_val; }

/// DoubleWrapper  != DoubleWrapper
inline bool operator!=(const DoubleWrapper & x, const DoubleWrapper & y)
{ return !(x == y); }
/// DoubleWrapper  != int
inline bool operator!=(const DoubleWrapper & x, int y)
{ return !(x == y); }
/// int != DoubleWrapper
inline bool operator!=(int x, const DoubleWrapper & y)
{ return !(x == y); }
/// DoubleWrapper  != unsigned int
inline bool operator!=(const DoubleWrapper & x, unsigned int y)
{ return !(x == y); }
/// unsigned int != DoubleWrapper
inline bool operator!=(unsigned int x, const DoubleWrapper & y)
{ return !(x == y); }
/// DoubleWrapper  != long
inline bool operator!=(const DoubleWrapper & x, long y)
{ return !(x == y); }
/// long != DoubleWrapper
inline bool operator!=(long x, const DoubleWrapper & y)
{ return !(x == y); }
/// DoubleWrapper  != unsigned long
inline bool operator!=(const DoubleWrapper & x, unsigned long y)
{ return !(x == y); }
/// unsigned long != DoubleWrapper
inline bool operator!=(unsigned long x, const DoubleWrapper & y)
{ return !(x == y); }
/// DoubleWrapper  != double
inline bool operator!=(const DoubleWrapper & x, double y)
{ return !(x == y); }
/// double != DoubleWrapper
inline bool operator!=(double x, const DoubleWrapper & y)
{ return !(x == y); }
/// DoubleWrapper  != LongWrapper
inline bool operator!=(const DoubleWrapper & x, const LongWrapper& y)
{ return !(x == y); }
/// LongWrapper != DoubleWrapper
inline bool operator!=(const LongWrapper& x, const DoubleWrapper & y)
{ return !(x == y); }

/// DoubleWrapper  >= DoubleWrapper
inline bool operator>=(const DoubleWrapper & x, const DoubleWrapper & y)
{ return x.d_val >= y.d_val; }
/// DoubleWrapper  >= int
inline bool operator>=(const DoubleWrapper & x, int y)
{ return x.d_val >= y; }
/// int >= DoubleWrapper
inline bool operator>=(int x, const DoubleWrapper & y)
{ return x >= y.d_val; }
/// DoubleWrapper  >= unsigned int
inline bool operator>=(const DoubleWrapper & x, unsigned int y)
{ return x.d_val >= y; }
/// unsigned int >= DoubleWrapper
inline bool operator>=(unsigned int x, const DoubleWrapper & y)
{ return x >= y.d_val; }
/// DoubleWrapper  >= long
inline bool operator>=(const DoubleWrapper & x, long y)
{ return x.d_val >= y; }
/// long >= DoubleWrapper
inline bool operator>=(long x, const DoubleWrapper & y)
{ return x >= y.d_val; }
/// DoubleWrapper  >= unsigned long
inline bool operator>=(const DoubleWrapper & x, unsigned long y)
{ return x.d_val >= y; }
/// unsigned long >= DoubleWrapper
inline bool operator>=(unsigned long x, const DoubleWrapper & y)
{ return x >= y.d_val; }
/// DoubleWrapper  >= double
inline bool operator>=(const DoubleWrapper & x, double y)
{ return x.d_val >= y; }
/// double >= DoubleWrapper
inline bool operator>=(double x, const DoubleWrapper & y)
{ return x >= y.d_val; }
/// DoubleWrapper  >= LongWrapper
inline bool operator>=(const DoubleWrapper & x, const LongWrapper& y)
{ return x.d_val >= y.l_val; }
/// LongWrapper >= DoubleWrapper
inline bool operator>=(const LongWrapper& x, const DoubleWrapper & y)
{ return x.l_val >= y.d_val; }

/// DoubleWrapper  <= DoubleWrapper
inline bool operator<=(const DoubleWrapper & x, const DoubleWrapper & y)
{ return x.d_val <= y.d_val; }
/// DoubleWrapper  <= int
inline bool operator<=(const DoubleWrapper & x, int y)
{ return x.d_val <= y; }
/// int <= DoubleWrapper
inline bool operator<=(int x, const DoubleWrapper & y)
{ return x <= y.d_val; }
/// DoubleWrapper  <= unsigned int
inline bool operator<=(const DoubleWrapper & x, unsigned int y)
{ return x.d_val <= y; }
/// unsigned int <= DoubleWrapper
inline bool operator<=(unsigned int x, const DoubleWrapper & y)
{ return x <= y.d_val; }
/// DoubleWrapper  <= long
inline bool operator<=(const DoubleWrapper & x, long y)
{ return x.d_val <= y; }
/// long <= DoubleWrapper
inline bool operator<=(long x, const DoubleWrapper & y)
{ return x <= y.d_val; }
/// DoubleWrapper  <= unsigned long
inline bool operator<=(const DoubleWrapper & x, unsigned long y)
{ return x.d_val <= y; }
/// unsigned long <= DoubleWrapper
inline bool operator<=(unsigned long x, const DoubleWrapper & y)
{ return x <= y.d_val; }
/// DoubleWrapper  <= double
inline bool operator<=(const DoubleWrapper & x, double y)
{ return x.d_val <= y; }
/// double <= DoubleWrapper
inline bool operator<=(double x, const DoubleWrapper & y)
{ return x <= y.d_val; }
/// DoubleWrapper  <= LongWrapper
inline bool operator<=(const DoubleWrapper & x, const LongWrapper& y)
{ return x.d_val <= y.l_val; }
/// LongWrapper <= DoubleWrapper
inline bool operator<=(const LongWrapper& x, const DoubleWrapper & y)
{ return x.l_val <= y.d_val; }

/// DoubleWrapper  > DoubleWrapper
inline bool operator>(const DoubleWrapper & x, const DoubleWrapper & y)
{ return x.d_val > y.d_val; }
/// DoubleWrapper  > int
inline bool operator>(const DoubleWrapper & x, int y)
{ return x.d_val > y; }
/// int > DoubleWrapper
inline bool operator>(int x, const DoubleWrapper & y)
{ return x > y.d_val; }
/// DoubleWrapper  > unsigned int
inline bool operator>(const DoubleWrapper & x, unsigned int y)
{ return x.d_val > y; }
/// unsigned int > DoubleWrapper
inline bool operator>(unsigned int x, const DoubleWrapper & y)
{ return x > y.d_val; }
/// DoubleWrapper  > long
inline bool operator>(const DoubleWrapper & x, long y)
{ return x.d_val > y; }
/// long > DoubleWrapper
inline bool operator>(long x, const DoubleWrapper & y)
{ return x > y.d_val; }
/// DoubleWrapper  > unsigned long
inline bool operator>(const DoubleWrapper & x, unsigned long y)
{ return x.d_val > y; }
/// unsigned long > DoubleWrapper
inline bool operator>(unsigned long x, const DoubleWrapper & y)
{ return x > y.d_val; }
/// DoubleWrapper  > double
inline bool operator>(const DoubleWrapper & x, double y)
{ return x.d_val > y; }
/// double > DoubleWrapper
inline bool operator>(double x, const DoubleWrapper & y)
{ return x > y.d_val; }
/// DoubleWrapper  > LongWrapper
inline bool operator>(const DoubleWrapper & x, const LongWrapper& y)
{ return x.d_val > y.l_val; }
/// LongWrapper > DoubleWrapper
inline bool operator>(const LongWrapper& x, const DoubleWrapper & y)
{ return x.l_val > y.d_val; }

/// DoubleWrapper  < DoubleWrapper
inline bool operator<(const DoubleWrapper & x, const DoubleWrapper & y)
{ return x.d_val < y.d_val; }
/// DoubleWrapper  < int
inline bool operator<(const DoubleWrapper & x, int y)
{ return x.d_val < y; }
/// int < DoubleWrapper
inline bool operator<(int x, const DoubleWrapper & y)
{ return x < y.d_val; }
/// DoubleWrapper  < unsigned int
inline bool operator<(const DoubleWrapper & x, unsigned int y)
{ return x.d_val < y; }
/// unsigned int < DoubleWrapper
inline bool operator<(unsigned int x, const DoubleWrapper & y)
{ return x < y.d_val; }
/// DoubleWrapper  < long
inline bool operator<(const DoubleWrapper & x, long y)
{ return x.d_val < y; }
/// long < DoubleWrapper
inline bool operator<(long x, const DoubleWrapper & y)
{ return x < y.d_val; }
/// DoubleWrapper  < unsigned long
inline bool operator<(const DoubleWrapper & x, unsigned long y)
{ return x.d_val < y; }
/// unsigned long < DoubleWrapper
inline bool operator<(unsigned long x, const DoubleWrapper & y)
{ return x < y.d_val; }
/// DoubleWrapper  < double
inline bool operator<(const DoubleWrapper & x, double y)
{ return x.d_val < y; }
/// double < DoubleWrapper
inline bool operator<(double x, const DoubleWrapper & y)
{ return x < y.d_val; }
/// DoubleWrapper  < LongWrapper
inline bool operator<(const DoubleWrapper & x, const LongWrapper& y)
{ return x.d_val < y.l_val; }
/// LongWrapper < DoubleWrapper
inline bool operator<(const LongWrapper& x, const DoubleWrapper & y)
{ return x.l_val < y.d_val; }
//@}

/// \addtogroup DoubleWrapperIostreamOperators
//@{
/// istream operator for <tt>DoubleWrapper</tt>
inline std::istream& operator>>(std::istream& is, DoubleWrapper& x)
{
  return (is >> x.d_val);
}
/// ostream operator for <tt>DoubleWrapper</tt>
inline std::ostream& operator<<(std::ostream& os, const DoubleWrapper& x) {
  return os << x.d_val;
} //@}

inline int sign(const DoubleWrapper& a) {
  if (a.d_val > 0) {
	return 1;
  } else if (a.d_val < 0) {
	return -1;
  }

  return 0;
}

// Chee: added this...
inline DoubleWrapper power(const DoubleWrapper &x, unsigned long p) {
  DoubleWrapper r;
  r.d_val = ::pow(x.d_val, p);
  return r;
}
inline DoubleWrapper pow(const DoubleWrapper &x, unsigned long p) {
  return power(x,p);
}
inline DoubleWrapper sqrt(const DoubleWrapper &x)
{ DoubleWrapper r(::sqrt(x.d_val)); return r; }

inline DoubleWrapper fabs(const DoubleWrapper &x) 
{ DoubleWrapper r(::fabs(x.d_val)); return r; }

inline DoubleWrapper atan(const DoubleWrapper &x)
{ DoubleWrapper r(::atan(x.d_val)); return r; }

inline DoubleWrapper sin(const DoubleWrapper &x)
{ DoubleWrapper r(::sin(x.d_val)); return r; }

inline DoubleWrapper cos(const DoubleWrapper &x)
{ DoubleWrapper r(::cos(x.d_val)); return r; }

// NOTE(narayan): This is needed by people who rely on the contract
// that the std math library provides.
//
// inline DoubleWrapper pow(const DoubleWrapper &x, unsigned long p) {
//   DoubleWrapper ret;
//   ret.d_val = ::pow(x.d_val, p);
//   return ret;
// }

CORE_END_NAMESPACE

/* ***************************************************
 * END of DoubleWrapper Definitions
 *************************************************** */

#endif /*__CORE_WRAPPERS_H__*/
