/****************************************************************************
 * BigFloat.inl -- Inline functions for BigFloat
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
 * $Id: BigFloat.inl,v 1.32 2010/11/23 17:58:36 exact Exp $
 ***************************************************************************/

/// \addtogroup BigFloatArithmeticOperators
//@{
/// BigFloat + BigFloat
inline BigFloat operator+(const BigFloat& x, const BigFloat& y)
{ BigFloat r; r.add(x, y); return r; }
/// BigFloat + int
inline BigFloat operator+(const BigFloat& x, int y)
{ BigFloat r; r.add(x, y); return r; }
/// int + BigFloat
inline BigFloat operator+(int x, const BigFloat& y)
{ BigFloat r; r.add(x, y); return r; }
/// BigFloat + unsigned int
inline BigFloat operator+(const BigFloat& x, unsigned int y)
{ BigFloat r; r.add(x, y); return r; }
/// unsigned int + BigFloat
inline BigFloat operator+(unsigned int x, const BigFloat& y)
{ BigFloat r; r.add(x, y); return r; }
/// BigFloat + long
inline BigFloat operator+(const BigFloat& x, long y)
{ BigFloat r; r.add(x, y); return r; }
/// long + BigFloat
inline BigFloat operator+(long x, const BigFloat& y)
{ BigFloat r; r.add(x, y); return r; }
/// BigFloat + unsigned long
inline BigFloat operator+(const BigFloat& x, unsigned long y)
{ BigFloat r; r.add(x, y); return r; }
/// unsigned long + BigFloat
inline BigFloat operator+(unsigned long x, const BigFloat& y)
{ BigFloat r; r.add(x, y); return r; }
/// BigFloat + double
inline BigFloat operator+(const BigFloat& x, double y)
{ BigFloat r; r.add(x, y); return r; }
/// double + BigFloat
inline BigFloat operator+(double x, const BigFloat& y)
{ BigFloat r; r.add(x, y); return r; }
/// BigFloat + BigInt
inline BigFloat operator+(const BigFloat& x, const BigInt& y)
{ BigFloat r; r.add(x, y); return r; }
/// BigInt + BigFloat
inline BigFloat operator+(const BigInt& x, const BigFloat& y)
{ BigFloat r; r.add(x, y); return r; }
/// BigFloat + BigRat
inline BigFloat operator+(const BigFloat& x, const BigRat& y)
{ BigFloat r; r.add(x, y); return r; }
/// BigRat + BigFloat
inline BigFloat operator+(const BigRat& x, const BigFloat& y)
{ BigFloat r; r.add(x, y); return r; }

/// BigFloat - BigFloat
inline BigFloat operator-(const BigFloat& x, const BigFloat& y)
{ BigFloat r; r.sub(x, y); return r; }
/// BigFloat - int
inline BigFloat operator-(const BigFloat& x, int y)
{ BigFloat r; r.sub(x, y); return r; }
/// int - BigFloat
inline BigFloat operator-(int x, const BigFloat& y)
{ BigFloat r; r.sub(x, y); return r; }
/// BigFloat - unsigned int
inline BigFloat operator-(const BigFloat& x, unsigned int y)
{ BigFloat r; r.sub(x, y); return r; }
/// unsigned int - BigFloat
inline BigFloat operator-(unsigned int x, const BigFloat& y)
{ BigFloat r; r.sub(x, y); return r; }
/// BigFloat - long
inline BigFloat operator-(const BigFloat& x, long y)
{ BigFloat r; r.sub(x, y); return r; }
/// long - BigFloat
inline BigFloat operator-(long x, const BigFloat& y)
{ BigFloat r; r.sub(x, y); return r; }
/// BigFloat - unsigned long
inline BigFloat operator-(const BigFloat& x, unsigned long y)
{ BigFloat r; r.sub(x, y); return r; }
/// unsigned long - BigFloat
inline BigFloat operator-(unsigned long x, const BigFloat& y)
{ BigFloat r; r.sub(x, y); return r; }
/// BigFloat - double
inline BigFloat operator-(const BigFloat& x, double y)
{ BigFloat r; r.sub(x, y); return r; }
/// double - BigFloat
inline BigFloat operator-(double x, const BigFloat& y)
{ BigFloat r; r.sub(x, y); return r; }
/// BigFloat - BigInt
inline BigFloat operator-(const BigFloat& x, const BigInt& y)
{ BigFloat r; r.sub(x, y); return r; }
/// BigInt - BigFloat
inline BigFloat operator-(const BigInt& x, const BigFloat& y)
{ BigFloat r; r.sub(x, y); return r; }
/// BigFloat - BigRat
inline BigFloat operator-(const BigFloat& x, const BigRat& y)
{ BigFloat r; r.sub(x, y); return r; }
/// BigRat - BigFloat
inline BigFloat operator-(const BigRat& x, const BigFloat& y)
{ BigFloat r; r.sub(x, y); return r; }

/// BigFloat * BigFloat
inline BigFloat operator*(const BigFloat& x, const BigFloat& y)
{ BigFloat r; r.mul(x, y); return r; }
/// BigFloat * int
inline BigFloat operator*(const BigFloat& x, int y)
{ BigFloat r; r.mul(x, y); return r; }
/// int * BigFloat
inline BigFloat operator*(int x, const BigFloat& y)
{ BigFloat r; r.mul(x, y); return r; }
/// BigFloat * unsigned int
inline BigFloat operator*(const BigFloat& x, unsigned int y)
{ BigFloat r; r.mul(x, y); return r; }
/// unsigned int * BigFloat
inline BigFloat operator*(unsigned int x, const BigFloat& y)
{ BigFloat r; r.mul(x, y); return r; }
/// BigFloat * long
inline BigFloat operator*(const BigFloat& x, long y)
{ BigFloat r; r.mul(x, y); return r; }
/// long * BigFloat
inline BigFloat operator*(long x, const BigFloat& y)
{ BigFloat r; r.mul(x, y); return r; }
/// BigFloat * unsigned long
inline BigFloat operator*(const BigFloat& x, unsigned long y)
{ BigFloat r; r.mul(x, y); return r; }
/// unsigned long * BigFloat
inline BigFloat operator*(unsigned long x, const BigFloat& y)
{ BigFloat r; r.mul(x, y); return r; }
/// BigFloat * double
inline BigFloat operator*(const BigFloat& x, double y)
{ BigFloat r; r.mul(x, y); return r; }
/// double * BigFloat
inline BigFloat operator*(double x, const BigFloat& y)
{ BigFloat r; r.mul(x, y); return r; }
/// BigFloat * BigInt
inline BigFloat operator*(const BigFloat& x, const BigInt& y)
{ BigFloat r; r.mul(x, y); return r; }
/// BigInt * BigFloat
inline BigFloat operator*(const BigInt& x, const BigFloat& y)
{ BigFloat r; r.mul(x, y); return r; }
/// BigFloat * BigRat
inline BigFloat operator*(const BigFloat& x, const BigRat& y)
{ BigFloat r; r.mul(x, y); return r; }
/// BigRat * BigFloat
inline BigFloat operator*(const BigRat& x, const BigFloat& y)
{ BigFloat r; r.mul(x, y); return r; }

/// BigFloat / BigFloat
inline BigFloat operator/(const BigFloat& x, const BigFloat& y)
{ BigFloat r; r.div(x, y); return r; }
/// BigFloat / int
inline BigFloat operator/(const BigFloat& x, int y)
{ BigFloat r; r.div(x, y); return r; }
/// int / BigFloat
inline BigFloat operator/(int x, const BigFloat& y)
{ BigFloat r; r.div(x, y); return r; }
/// BigFloat / unsigned int
inline BigFloat operator/(const BigFloat& x, unsigned int y)
{ BigFloat r; r.div(x, y); return r; }
/// unsigned int / BigFloat
inline BigFloat operator/(unsigned int x, const BigFloat& y)
{ BigFloat r; r.div(x, y); return r; }
/// BigFloat / long
inline BigFloat operator/(const BigFloat& x, long y)
{ BigFloat r; r.div(x, y); return r; }
/// long / BigFloat
inline BigFloat operator/(long x, const BigFloat& y)
{ BigFloat r; r.div(x, y); return r; }
/// BigFloat / unsigned long
inline BigFloat operator/(const BigFloat& x, unsigned long y)
{ BigFloat r; r.div(x, y); return r; }
/// unsigned long / BigFloat 
inline BigFloat operator/(unsigned long x, const BigFloat& y)
{ BigFloat r; r.div(x, y); return r; }
/// BigFloat / double
inline BigFloat operator/(const BigFloat& x, double y)
{ BigFloat r; r.div(x, y); return r; }
/// double / BigFloat 
inline BigFloat operator/(double x, const BigFloat& y)
{ BigFloat r; r.div(x, y); return r; }
/// BigFloat / BigInt
inline BigFloat operator/(const BigFloat& x, const BigInt& y)
{ BigFloat r; r.div(x, y); return r; }
/// BigInt / BigFloat
inline BigFloat operator/(const BigInt& x, const BigFloat& y)
{ BigFloat r; r.div(x, y); return r; }
/// BigFloat / BigRat
inline BigFloat operator/(const BigFloat& x, const BigRat& y)
{ BigFloat r; r.div(x, y); return r; }
/// BigRat / BigFloat
inline BigFloat operator/(const BigRat& x, const BigFloat& y)
{ BigFloat r; r.div(x, y); return r; }
//@}

/// BigFloat  << int
inline BigFloat operator<<(const BigFloat & x, int y)
{ BigFloat r(x); r.mul_2exp(y); return r; }
/// BigFloat  << unsigned int
inline BigFloat operator<<(const BigFloat & x, unsigned int y)
{ BigFloat r(x); r.mul_2exp(y); return r; }
/// BigFloat  << long
inline BigFloat operator<<(const BigFloat & x, long y)
{ BigFloat r(x); r.mul_2exp(y); return r; }
/// BigFloat  << unsigned long
inline BigFloat operator<<(const BigFloat & x, unsigned long y)
{ BigFloat r(x); r.mul_2exp(y); return r; }
/// BigFloat  >> int
inline BigFloat operator>>(const BigFloat & x, int y)
{ BigFloat r(x); r.div_2exp(y); return r; }
/// BigFloat  >> unsigned int
inline BigFloat operator>>(const BigFloat & x, unsigned int y)
{ BigFloat r(x); r.div_2exp(y); return r; }
/// BigFloat  >> long
inline BigFloat operator>>(const BigFloat & x, long y)
{ BigFloat r(x); r.div_2exp(y); return r; }
/// BigFloat  >> unsigned long
inline BigFloat operator>>(const BigFloat & x, unsigned long y)
{ BigFloat r(x); r.div_2exp(y); return r; }
//@}

/// \addtogroup BigFloatComparisonOperators
//@{
/// BigFloat  == BigFloat 
inline bool operator==(const BigFloat & x, const BigFloat & y)
{ return x.cmp(y) == 0; }
/// BigFloat  == int
inline bool operator==(const BigFloat & x, int y)
{ return x.cmp(y) == 0; }
/// int == BigFloat 
inline bool operator==(int x, const BigFloat & y)
{ return y.cmp(x) == 0; }
/// BigFloat  == unsigned int
inline bool operator==(const BigFloat & x, unsigned int y)
{ return x.cmp(y) == 0; }
/// unsigned int == BigFloat 
inline bool operator==(unsigned int x, const BigFloat & y)
{ return y.cmp(x) == 0; }
/// BigFloat  == long
inline bool operator==(const BigFloat & x, long y)
{ return x.cmp(y) == 0; }
/// long == BigFloat 
inline bool operator==(long x, const BigFloat & y)
{ return y.cmp(x) == 0; }
/// BigFloat  == unsigned long
inline bool operator==(const BigFloat & x, unsigned long y)
{ return x.cmp(y) == 0; }
/// unsigned long == BigFloat 
inline bool operator==(unsigned long x, const BigFloat & y)
{ return y.cmp(x) == 0; }
/// BigFloat  == double
inline bool operator==(const BigFloat & x, double y)
{ return x.cmp(y) == 0; }
/// double == BigFloat 
inline bool operator==(double x, const BigFloat & y)
{ return y.cmp(x) == 0; }
/// BigFloat  == BigInt
inline bool operator==(const BigFloat & x, const BigInt& y)
{ return x.cmp(y) == 0; }
/// BigInt == BigFloat 
inline bool operator==(const BigInt& x, const BigFloat & y)
{ return y.cmp(x) == 0; }
/// BigFloat  == BigRat
inline bool operator==(const BigFloat & x, const BigRat& y)
{ return x.cmp(y) == 0; }
/// BigRat == BigFloat 
inline bool operator==(const BigRat& x, const BigFloat & y)
{ return y.cmp(x) == 0; }

/// BigFloat  != BigFloat 
inline bool operator!=(const BigFloat & x, const BigFloat & y)
{ return x.cmp(y) != 0; }
/// BigFloat  != int
inline bool operator!=(const BigFloat & x, int y)
{ return x.cmp(y) != 0; }
/// int != BigFloat 
inline bool operator!=(int x, const BigFloat & y)
{ return y.cmp(x) != 0; }
/// BigFloat  != unsigned int
inline bool operator!=(const BigFloat & x, unsigned int y)
{ return x.cmp(y) != 0; }
/// unsigned int != BigFloat 
inline bool operator!=(unsigned int x, const BigFloat & y)
{ return y.cmp(x) != 0; }
/// BigFloat  != long
inline bool operator!=(const BigFloat & x, long y)
{ return x.cmp(y) != 0; }
/// long != BigFloat 
inline bool operator!=(long x, const BigFloat & y)
{ return y.cmp(x) != 0; }
/// BigFloat  != unsigned long
inline bool operator!=(const BigFloat & x, unsigned long y)
{ return x.cmp(y) != 0; }
/// unsigned long != BigFloat 
inline bool operator!=(unsigned long x, const BigFloat & y)
{ return y.cmp(x) != 0; }
/// BigFloat  != double
inline bool operator!=(const BigFloat & x, double y)
{ return x.cmp(y) != 0; }
/// double != BigFloat 
inline bool operator!=(double x, const BigFloat & y)
{ return y.cmp(x) != 0; }
/// BigFloat  != BigInt
inline bool operator!=(const BigFloat & x, const BigInt& y)
{ return x.cmp(y) != 0; }
/// BigInt != BigFloat 
inline bool operator!=(const BigInt& x, const BigFloat & y)
{ return y.cmp(x) != 0; }
/// BigFloat  != BigRat
inline bool operator!=(const BigFloat & x, const BigRat& y)
{ return x.cmp(y) != 0; }
/// BigRat != BigFloat 
inline bool operator!=(const BigRat& x, const BigFloat & y)
{ return y.cmp(x) != 0; }

/// BigFloat  >= BigFloat 
inline bool operator>=(const BigFloat & x, const BigFloat & y)
{ return x.cmp(y) >= 0; }
/// BigFloat  >= int
inline bool operator>=(const BigFloat & x, int y)
{ return x.cmp(y) >= 0; }
/// int >= BigFloat 
inline bool operator>=(int x, const BigFloat & y)
{ return y.cmp(x) <= 0; }
/// BigFloat  >= unsigned int
inline bool operator>=(const BigFloat & x, unsigned int y)
{ return x.cmp(y) >= 0; }
/// unsigned int >= BigFloat 
inline bool operator>=(unsigned int x, const BigFloat & y)
{ return y.cmp(x) <= 0; }
/// BigFloat  >= long
inline bool operator>=(const BigFloat & x, long y)
{ return x.cmp(y) >= 0; }
/// long >= BigFloat 
inline bool operator>=(long x, const BigFloat & y)
{ return y.cmp(x) <= 0; }
/// BigFloat  >= unsigned long
inline bool operator>=(const BigFloat & x, unsigned long y)
{ return x.cmp(y) >= 0; }
/// unsigned long >= BigFloat 
inline bool operator>=(unsigned long x, const BigFloat & y)
{ return y.cmp(x) <= 0; }
/// BigFloat  >= double
inline bool operator>=(const BigFloat & x, double y)
{ return x.cmp(y) >= 0; }
/// double >= BigFloat 
inline bool operator>=(double x, const BigFloat & y)
{ return y.cmp(x) <= 0; }
/// BigFloat  >= BigInt
inline bool operator>=(const BigFloat & x, const BigInt& y)
{ return x.cmp(y) >= 0; }
/// BigInt >= BigFloat 
inline bool operator>=(const BigInt& x, const BigFloat & y)
{ return y.cmp(x) <= 0; }
/// BigFloat  >= BigRat
inline bool operator>=(const BigFloat & x, const BigRat& y)
{ return x.cmp(y) >= 0; }
/// BigRat >= BigFloat 
inline bool operator>=(const BigRat& x, const BigFloat & y)
{ return y.cmp(x) <= 0; }

/// BigFloat  <= BigFloat 
inline bool operator<=(const BigFloat & x, const BigFloat & y)
{ return x.cmp(y) <= 0; }
/// BigFloat  <= int
inline bool operator<=(const BigFloat & x, int y)
{ return x.cmp(y) <= 0; }
/// int <= BigFloat 
inline bool operator<=(int x, const BigFloat & y)
{ return y.cmp(x) >= 0; }
/// BigFloat  <= unsigned int
inline bool operator<=(const BigFloat & x, unsigned int y)
{ return x.cmp(y) <= 0; }
/// unsigned int <= BigFloat 
inline bool operator<=(unsigned int x, const BigFloat & y)
{ return y.cmp(x) >= 0; }
/// BigFloat  <= long
inline bool operator<=(const BigFloat & x, long y)
{ return x.cmp(y) <= 0; }
/// long <= BigFloat 
inline bool operator<=(long x, const BigFloat & y)
{ return y.cmp(x) >= 0; }
/// BigFloat  <= unsigned long
inline bool operator<=(const BigFloat & x, unsigned long y)
{ return x.cmp(y) <= 0; }
/// unsigned long <= BigFloat 
inline bool operator<=(unsigned long x, const BigFloat & y)
{ return y.cmp(x) >= 0; }
/// BigFloat  <= double
inline bool operator<=(const BigFloat & x, double y)
{ return x.cmp(y) <= 0; }
/// double <= BigFloat 
inline bool operator<=(double x, const BigFloat & y)
{ return y.cmp(x) >= 0; }
/// BigFloat  <= BigInt
inline bool operator<=(const BigFloat & x, const BigInt& y)
{ return x.cmp(y) <= 0; }
/// BigInt <= BigFloat 
inline bool operator<=(const BigInt& x, const BigFloat & y)
{ return y.cmp(x) >= 0; }
/// BigFloat  <= BigRat
inline bool operator<=(const BigFloat & x, const BigRat& y)
{ return x.cmp(y) <= 0; }
/// BigRat <= BigFloat 
inline bool operator<=(const BigRat& x, const BigFloat & y)
{ return y.cmp(x) >= 0; }

/// BigFloat  > BigFloat 
inline bool operator>(const BigFloat & x, const BigFloat & y)
{ return x.cmp(y) > 0; }
/// BigFloat  > int
inline bool operator>(const BigFloat & x, int y)
{ return x.cmp(y) > 0; }
/// int > BigFloat 
inline bool operator>(int x, const BigFloat & y)
{ return y.cmp(x) < 0; }
/// BigFloat  > unsigned int
inline bool operator>(const BigFloat & x, unsigned int y)
{ return x.cmp(y) > 0; }
/// unsigned int > BigFloat 
inline bool operator>(unsigned int x, const BigFloat & y)
{ return y.cmp(x) < 0; }
/// BigFloat  > long
inline bool operator>(const BigFloat & x, long y)
{ return x.cmp(y) > 0; }
/// long > BigFloat 
inline bool operator>(long x, const BigFloat & y)
{ return y.cmp(x) < 0; }
/// BigFloat  > unsigned long
inline bool operator>(const BigFloat & x, unsigned long y)
{ return x.cmp(y) > 0; }
/// unsigned long > BigFloat 
inline bool operator>(unsigned long x, const BigFloat & y)
{ return y.cmp(x) < 0; }
/// BigFloat  > double
inline bool operator>(const BigFloat & x, double y)
{ return x.cmp(y) > 0; }
/// double > BigFloat 
inline bool operator>(double x, const BigFloat & y)
{ return y.cmp(x) < 0; }
/// BigFloat  > BigInt
inline bool operator>(const BigFloat & x, const BigInt& y)
{ return x.cmp(y) > 0; }
/// BigInt > BigFloat 
inline bool operator>(const BigInt& x, const BigFloat & y)
{ return y.cmp(x) < 0; }
/// BigFloat  > BigRat
inline bool operator>(const BigFloat & x, const BigRat& y)
{ return x.cmp(y) > 0; }
/// BigRat > BigFloat 
inline bool operator>(const BigRat& x, const BigFloat & y)
{ return y.cmp(x) < 0; }

/// BigFloat  < BigFloat 
inline bool operator<(const BigFloat & x, const BigFloat & y)
{ return x.cmp(y) < 0; }
/// BigFloat  < int
inline bool operator<(const BigFloat & x, int y)
{ return x.cmp(y) < 0; }
/// int < BigFloat 
inline bool operator<(int x, const BigFloat & y)
{ return y.cmp(x) > 0; }
/// BigFloat  < unsigned int
inline bool operator<(const BigFloat & x, unsigned int y)
{ return x.cmp(y) < 0; }
/// unsigned int < BigFloat 
inline bool operator<(unsigned int x, const BigFloat & y)
{ return y.cmp(x) > 0; }
/// BigFloat  < long
inline bool operator<(const BigFloat & x, long y)
{ return x.cmp(y) < 0; }
/// long < BigFloat 
inline bool operator<(long x, const BigFloat & y)
{ return y.cmp(x) > 0; }
/// BigFloat  < unsigned long
inline bool operator<(const BigFloat & x, unsigned long y)
{ return x.cmp(y) < 0; }
/// unsigned long < BigFloat 
inline bool operator<(unsigned long x, const BigFloat & y)
{ return y.cmp(x) > 0; }
/// BigFloat  < double
inline bool operator<(const BigFloat & x, double y)
{ return x.cmp(y) < 0; }
/// double < BigFloat 
inline bool operator<(double x, const BigFloat & y)
{ return y.cmp(x) > 0; }
/// BigFloat  < BigInt
inline bool operator<(const BigFloat & x, const BigInt& y)
{ return x.cmp(y) < 0; }
/// BigInt < BigFloat 
inline bool operator<(const BigInt& x, const BigFloat & y)
{ return y.cmp(x) > 0; }
/// BigFloat  < BigRat
inline bool operator<(const BigFloat & x, const BigRat& y)
{ return x.cmp(y) < 0; }
/// BigRat < BigFloat 
inline bool operator<(const BigRat& x, const BigFloat & y)
{ return y.cmp(x) > 0; }
//@}

/// \addtogroup BigFloatIostreamOperators
//@{
/// istream operator for <tt>BigFloat</tt>
inline std::istream& operator>>(std::istream& is, BigFloat& x)
{
  // Jihun Nov,2010
  // this is hack.
  // Users should not call this operator with infinity input digits
  // Bigrat number type must be used instead.
  if(getDefaultInputDigits() == CORE_INFTY) {
    x.set_prec(52); // default double precision
  } else {
    x.set_prec(digits2bits(getDefaultInputDigits()));
  }
  return is >> x.mp();
}
/// ostream operator for <tt>BigFloat</tt>
inline std::ostream& operator<<(std::ostream& os, const BigFloat& x) {
  return os <<
  mpfr2str(x.mp(),
      (std::min)(
      (unsigned long)get_output_precision(os),
      bits2digits(x.get_prec()+1)),
    get_output_base(os),
    get_output_fmt(os),
    get_output_rounding_mode(),
    get_output_showpoint(os),
    get_output_showpos(os),
    get_output_uppercase(os));
} //@}

/// \addtogroup BigFloatGlobalFunctions
//@{
/// square root
inline BigFloat sqrt(const BigFloat& x, prec_t prec = getDefaultBFradicalPrec())
{ BigFloat r(0, prec); r.sqrt(x); return r; }
/// cubic root
inline BigFloat cbrt(const BigFloat& x, prec_t prec = getDefaultBFradicalPrec())
{ BigFloat r(0, prec); r.cbrt(x); return r; }
/// k-th root
inline BigFloat root(const BigFloat& x, unsigned long k, prec_t prec = getDefaultBFradicalPrec())
{ BigFloat r(0, prec); r.root(x, k); return r; }
inline BigFloat div(const BigFloat& x, const BigFloat& y, prec_t prec = getDefaultBFradicalPrec())
{ BigFloat r(0, prec); r.div(x, y); return r; }
/// left shift for <tt>unsigned int</tt>
inline BigFloat mul_2exp(const BigFloat& x, int y)
{ BigFloat r(x); return r.mul_2exp(y); }
/// left shift for <tt>unsigned int</tt>
inline BigFloat mul_2exp(const BigFloat& x, unsigned int y)
{ BigFloat r(x); return r.mul_2exp(y); }
/// left shift for <tt>long</tt>
inline BigFloat mul_2exp(const BigFloat& x, long y)
{ BigFloat r(x); return r.mul_2exp(y); }
/// left shift for <tt>unsigned long</tt>
inline BigFloat mul_2exp(const BigFloat& x, unsigned long y)
{ BigFloat r(x); return r.mul_2exp(y); }
/// right shift for <tt>int</tt>
inline BigFloat div_2exp(const BigFloat& x, int y)
{ BigFloat r(x); return r.div_2exp(y); }
/// right shift for <tt>unsigned int</tt>
inline BigFloat div_2exp(const BigFloat& x, unsigned int y)
{ BigFloat r(x); return r.div_2exp(y); }
/// right shift for <tt>long</tt>
inline BigFloat div_2exp(const BigFloat& x, long y)
{ BigFloat r(x); return r.div_2exp(y); }
/// right shift for <tt>unsigned long</tt>
inline BigFloat div_2exp(const BigFloat& x, unsigned long y)
{ BigFloat r(x); return r.div_2exp(y); }
/// divide by 2
inline BigFloat div2(const BigFloat& x)
{ BigFloat r(x); return r.div2(); }

/// minStar(m,n) returns the min-star of m and n
inline long minStar(long m, long n) {
  if (m*n <= 0) return 0;
  if (m>0)
    return (std::min)(m, n);
  else
    return (std::max)(m, n);
}
/// read from file
inline void readFromFile(BigFloat& z, std::istream& in) {
  char c; char dummy[256];
  int base;
  std::string str;
  do {
    in >> c;
    if (c != 'f')
      in.getline(dummy, 256);
  } while (c != 'f');
  in >> base;
  in >> str; str.erase(0, 1);
  z.set(str, base);
}
/// write to file
inline void writeToFile(const BigFloat& z, std::ostream& out, int base=10) {
  out << 'f';
  out << base;
  out << '|' << z.get_str(0, base) << std::endl; 
}

/// isDivisible(a,b) = "is a divisible by b"
/**     Assuming that a and  b are in canonized forms.
        Defined to be true if mantissa(b) | mantissa(a) && 
        exp(b) = min*(exp(b), exp(a)).
 *      This concepts assume a and b are exact BigFloat.
 */
inline bool isDivisible(const BigFloat& a, const BigFloat& b) {
  // assert: a and b are exact BigFloats.
  BigInt m_a, m_b;
  exp_t e_a = a.get_z_exp(m_a);
  exp_t e_b = b.get_z_exp(m_b);
  if (sign(m_a) == 0) return true;
  if (sign(m_b) == 0) return false;

  unsigned long bin_a = getBinExpo(m_a);
  unsigned long bin_b = getBinExpo(m_b);
  
  m_a >>= bin_a;
  m_b >>= bin_b;
  e_a += bin_a;
  e_b += bin_b;

  long dx = minStar(e_a, e_b);
  return isDivisible(m_a, m_b) && (dx == e_b); 
}

inline bool isDivisible(double x, double y) {
  //Are these exact?
  return isDivisible(BigFloat(x), BigFloat(y));
}

/// div_exact(x,y) returns the BigFloat quotient of x divided by y
/**     This is defined only if isDivisible(x,y).
 */
// Chee (8/1/2004)   The definition of div_exact(x,y) 
//   ensures that Polynomials<NT> works with NT=BigFloat and NT=double:
// Jihun (7/3/2006)  This was not exactly working when NT=BigFloat
//   z needs more precision than default.  Problem fixed.
inline BigFloat div_exact(const BigFloat& x, const BigFloat& y) {
  BigFloat z;

  assert (isDivisible(x,y));
 
  BigInt m_x, m_y;
  x.get_z_exp(m_x);
  y.get_z_exp(m_y);

  unsigned long bin_x = getBinExpo(m_x);
  unsigned long bin_y = getBinExpo(m_y);

  if ((prec_t)bin_x > x.get_prec()) bin_x = x.get_prec();
  if ((prec_t)bin_y > y.get_prec()) bin_y = y.get_prec(); 

  long valprec = (x.get_prec()-bin_x)-(y.get_prec()-bin_y)+1;

  if (valprec < 0) valprec = -valprec;
  if (valprec < 2) valprec = 2;
  
  z.div(x,y,valprec);
  return z;
}

inline BigFloat div_exact(double x, double y) {
  return div_exact(BigFloat(x), BigFloat(y));
}

// Remark: there is another notion of "exact division" for BigFloats,
//      and that is to make the division return an "exact" BigFloat
//      i.e., err()=0.  

/// gcd(a,b) =  BigFloat(gcd(a.mantissa,b.matissa), min(a.exp(), b.exp()) )
inline BigFloat gcd(const BigFloat& a, const BigFloat& b) {
  BigInt m_a, m_b;
  exp_t e_a = a.get_z_exp(m_a);
  exp_t e_b = b.get_z_exp(m_b);
  if (sign(m_a) == 0) return core_abs(b);
  if (sign(m_b) == 0) return core_abs(a);

  unsigned long bin_a = getBinExpo(m_a);
  unsigned long bin_b = getBinExpo(m_b);

  m_a >>= bin_a;
  m_b >>= bin_b;
  e_a += bin_a;
  e_b += bin_b;

  BigInt r = gcd(m_a, m_b);
  long dx = minStar(e_a, e_b);

  // return x*2^{dx}
  BigFloat x(r);
  x.mul_2exp(dx);
  return x;
}

inline BigFloat fabs(const BigFloat& a) {
  return abs(a);
}

inline double Todouble(const BigFloat& a) {
  return a.doubleValue();
}

inline BigFloat log10(const BigFloat& a, prec_t prec, rnd_t rnd = MPFR_RND) {
  BigFloat r; r.log10(a,prec,rnd); return r;
}

//@}

#ifndef CORE_DISABLE_OLDNAMES 
/// \addtogroup BigFloatBackCompatiableFunctions
//@{
/// comparison
inline int cmp(const BigFloat& x, const BigFloat& y) { return x.cmp(y); }
/// sign 
inline int sign(const BigFloat& a) { return a.sgn(); }
inline int sgn(const BigFloat& a) { return a.sgn(); }
/// abs
inline BigFloat abs(const BigFloat& a) { BigFloat r; r.abs(a); return r; }
/// neg
inline BigFloat neg(const BigFloat& a) { BigFloat r; r.neg(a); return r; }
/// pow 
inline BigFloat power(const BigFloat& a, unsigned long p) 
{ BigFloat r; r.pow(a, p, (std::max)(a.get_prec()*p, 2UL)); return r; }
/// floorlg
inline long floorLg(const BigFloat& a) { return a.lMSB(); }
inline long floorlg(const BigFloat& a) { return a.lMSB(); }
/// ceillg
inline long ceilLg(const BigFloat& a) { return a.uMSB(); }
inline long ceillg(const BigFloat& a) { return a.uMSB(); }
//@}
#endif

