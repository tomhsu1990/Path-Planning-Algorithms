/****************************************************************************
 * RootBounds.h -- Constructive root bounds
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
 * $Id: RootBounds.h,v 1.18 2010/11/23 17:58:37 exact Exp $
 ***************************************************************************/
#ifndef __CORE_ROOTBOUNDS_H__
#define __CORE_ROOTBOUNDS_H__

#include <iostream>
#include <sstream>
#include <string>
#include <limits>

CORE_BEGIN_NAMESPACE

/// Dummy Root Bound (minimal root bound class)
template <int bound = 1000>
class DummyRootBd {
  typedef DummyRootBd thisClass;
public:
#ifdef CORE_DEBUG_ROOTBOUND
  void dump() const {}
#endif
  bool is_constructive() const { return false; }
  unsigned long get_bound() const { return bound; }
  template <typename T> void set(const T& value) {}
  void neg(const thisClass& child) {}
  void root(const thisClass& child, unsigned long k) {}
  void addsub(const thisClass& f, const thisClass& s) {}
  void mul(const thisClass& f, const thisClass& s) {}
  void div(const thisClass& f, const thisClass& s) {}
};

const double log_5 = std::log(double(5)) / std::log(double(2));
//const unsigned long Inf = std::numeric_limits<unsigned long>::infinity();

template <typename Kernel = BigFloat2>
class BfmsskRootBd_BigFloat {
  typedef BfmsskRootBd_BigFloat thisClass;
  typedef typename Kernel::ZT ZT;
  typedef typename Kernel::QT QT;
  typedef typename Kernel::FT FT;

  bool m_visit;
  unsigned long d;
  // u25 and l25 are maintained as log base 2 values.
  FT u25;
  FT l25;
  unsigned long v2p;
  unsigned long v2m;
  unsigned long v5p;
  unsigned long v5m;
  const prec_t defaultPrec;
  
  FT lg5; // const value of log 5 base 2.

FT value; // for debugging
bool debug;

  void cancelPowers() {
	subtractMin(v2p,v2m); subtractMin(v5p,v5m);
  }
  void subtractMin(unsigned long& x, unsigned long& y) {
	if(x > y) { x -= y; y = 0; } else { y -= x; x = 0; }
  }
public:
  BfmsskRootBd_BigFloat() : m_visit(false), d(1), u25(0), l25(0), v2p(0), v2m(0), v5p(0), v5m(0), defaultPrec(53) {
    lg5.log2(5,defaultPrec,BF_RNDU);
	debug = false;
  }
  void dump()  
  { std::cerr<<"[d,u25,l25]="<<d<<","<<u25<<","<<l25<<std::endl; 
    std::cerr<<"[v2p,v2m,v5p,v5m]="<<v2p<<","<<v2m<<","<<v5p<<","<<v5m<<std::endl;
    std::cerr<<"lg5(v5p-v5m)="<<lg5*v5p-lg5*v5m<<std::endl;
if(debug)
std::cerr<<"value="<<value<<std::endl;
  }
  unsigned long get_deg()
  { return d; }
  void set_visit(bool visit)
  { m_visit = visit; }
  bool get_visit()
  { return m_visit; }
  bool is_constructive() const
  { return true; }
  unsigned long get_bound(unsigned long deg) {
    FT bound = l25+u25*(deg-1UL)+v2m+lg5*v5m-v2p-lg5*v5p;
	if(bound<=0) return 1UL;

    // bound cannot be zero, because of mpfr requirement.
        return (std::max)(bound.get_ui(BF_RNDU), 1UL);
  }
  void set(long value) 			{ set(ZT(value)); }
  void set(unsigned long value) { set(ZT(value)); }
  void set(double value) 		{ set(FT(value)); }
  void set(const ZT& value) {
if(debug)
this->value = value;
    if (value == 0) return;
	
	ZT temp(value);

    v2p = temp.get_k_exp(temp, 2UL);
    v5p = temp.get_k_exp(temp, 5UL);

    v2m = v5m = 0;
    
    u25.log2(abs(temp),defaultPrec,BF_RNDU);
    l25 = 0;	
  }
  
  void set(const QT& value) {
if(debug)
this->value = value;  
    if (value == 0) return;

    ZT num = value.num();
    ZT den = value.den();

    v2p = num.get_k_exp(num, 2UL);
    v2m = den.get_k_exp(den, 2UL);

    v5p = num.get_k_exp(num, 5UL);
    v5m = den.get_k_exp(den, 5UL);

	cancelPowers();
	
    u25.log2(abs(num),defaultPrec,BF_RNDU);
    l25.log2(abs(den),defaultPrec,BF_RNDU);
  }
  void set(const FT& value) {
// TODO : 
// deal with the case when the value is integer
if(debug)
this->value = value;
    if(value == 0) return; 
	FT temp(value); temp.remove_trailing_zeros();	
    ZT x; exp_t e = temp.get_z_exp(x);

	// x has powers of 2,which is from the precision of f.
	// extract powers of 2 and 5(e2 and e5)
    // then compute the true exponent e' = e2 + e.
/*
	// Old routine
    if (e >= 0)
      v2p = e;
    else
      v2m = -e;
	
    v5p = 0;
	v5m = 0;

*/
	e += x.get_k_exp(x,2UL);	
    if (e >= 0)
      v2p = e;
    else
      v2m = -e;
	
    v5p = x.get_k_exp(x,5UL);
	v5m = 0;

	cancelPowers();
	
    u25.log2(abs(x),defaultPrec,BF_RNDU);
    l25 = 0;
  }
  void set(const Kernel& value) {
    if (value.sgn() >= 0)
      set(value.getLeft());
    else
      set(value.getRight());
  }
  void neg(const thisClass& child) {
    u25 = child.u25; l25 = child.l25;
    v2p = child.v2p; v2m = child.v2m;
    v5p = child.v5p; v5m = child.v5m;
  }

  void root(thisClass& child, unsigned long k) {
    d = k;

    if (child.v2p + lg5*child.v5p + child.u25 >=
        child.v2m + lg5*child.v5m + child.l25) {
      unsigned long vtilda2 = child.v2p + (k - 1UL) * child.v2m;
      v2p = vtilda2 / k;
      v2m = child.v2m;
      unsigned long vmod2;
      vmod2 = vtilda2 - k * v2p; // == vtilda2 % k
      unsigned long vtilda5 = child.v5p + (k - 1UL) * child.v5m;
      v5p = vtilda5 / k;
      v5m = child.v5m;
      unsigned long vmod5;
      vmod5 = vtilda5 - k * v5p; // == vtilda5 % k

      u25.div(child.u25 + (k - 1UL) * child.l25 + vmod2 + lg5*(vmod5) + 1UL,k,defaultPrec,BF_RNDU);
      l25 = child.l25;
    } else {
      unsigned long vtilda2 = (k - 1UL) * child.v2p + child.v2m;
      v2p = child.v2p;
      v2m = vtilda2 / k;
      unsigned long vmod2;
      vmod2 = vtilda2 - k * v2m; // == vtilda2 % k
      unsigned long vtilda5 = (k - 1UL) * child.v5p + child.v5m;
      v5p = child.v5p;
      v5m = vtilda5 / k;
      unsigned long vmod5;
      vmod5 = vtilda5 - k * v5m; // == vtilda5 % k

      u25 = child.u25;
      l25.div((k - 1UL) * child.u25 + child.l25 + vmod2 + lg5*(vmod5) + 1UL,k,defaultPrec,BF_RNDU);
    }

	roundprec();
  }
  void addsub(thisClass& f, thisClass& s) {
    v2p = (std::min)(f.v2p + s.v2m, f.v2m + s.v2p);
    v2m = f.v2m + s.v2m;
    v5p = (std::min)(f.v5p + s.v5m, f.v5m + s.v5p);
    v5m = f.v5m + s.v5m;
	
	FT A = f.v2p + s.v2m - v2p + lg5*(f.v5p + s.v5m - v5p) + f.u25 + s.l25;
	FT B = f.v2m + s.v2p - v2p + lg5*(f.v5m + s.v5p - v5p) + f.l25 + s.u25;
if(debug){
std::cerr << "addsub node" << std::endl;
	
std::cerr << std::endl << "dump of f" << std::endl; f.dump();
std::cerr << std::endl << "dump of s" << std::endl; s.dump();
std::cerr << std::endl;
std::cerr << "B="<<B<<std::endl;
std::cerr << "A="<<A<<std::endl;	
std::cerr << "B-A="<<B-A<<std::endl;
}

	FT exp_2; exp_2.exp2(B-A,defaultPrec,BF_RNDU);
	FT B_; B_.log2(1+exp_2,defaultPrec,BF_RNDU);

    u25 = A+B_;
    l25 = f.l25 + s.l25;
if(debug){
std::cerr << "dump of self" << std::endl;
dump();
}

	roundprec();
  }
  void mul(thisClass& f,thisClass& s) {
    v2p = f.v2p + s.v2p;
    v2m = f.v2m + s.v2m;
    v5p = f.v5p + s.v5p;
    v5m = f.v5m + s.v5m;
    u25 = f.u25 + s.u25;
    l25 = f.l25 + s.l25;
if(debug){
std::cerr << "mull node" << std::endl;	
std::cerr << std::endl << "dump of f" << std::endl; f.dump();
std::cerr << std::endl << "dump of s" << std::endl; s.dump();
std::cerr << std::endl;
}
	roundprec();
  }
  void div(thisClass& f, thisClass& s) { 
    v2p = f.v2p + s.v2m;
    v2m = f.v2m + s.v2p;
    v5p = f.v5p + s.v5m;
    v5m = f.v5m + s.v5p;
    u25 = f.u25 + s.l25;
    l25 = f.l25 + s.u25;

if(debug){
std::cerr << "div node" << std::endl;	
std::cerr << std::endl << "dump of f" << std::endl; f.dump();
std::cerr << std::endl << "dump of s" << std::endl; s.dump();
std::cerr << std::endl;
}
	roundprec();
  }
private:
  void roundprec() {
//    u25.remove_trailing_zeros(); 
//    l25.remove_trailing_zeros(); 
    u25.prec_round(defaultPrec,BF_RNDU); 
	l25.prec_round(defaultPrec,BF_RNDU);
    return;
  }
};

template <typename Kernel = BigFloat2>
class BfmsskRootBd_double {
public:
  typedef BfmsskRootBd_double thisClass;
  typedef typename Kernel::ZT ZT;
  typedef typename Kernel::QT QT;
  typedef typename Kernel::FT FT;
  bool			m_visit;
  unsigned long d;
  
  // u25 and l25 are maintained as upper bounds of the lg values.
  double u25, l25;
  unsigned long v2p, v2m, v5p, v5m;

  void cancelPowers() {
	subtractMin(v2p,v2m); subtractMin(v5p,v5m);
  }
  void subtractMin(unsigned long& x, unsigned long& y) {
	if(x > y) { x -= y; y = 0; } else { y -= x; x = 0; }
  }
  
  // Quick implementation of rounding up operators of machine double type.
  // rounding up for double type.
  double add(const double x, const double y) {
	assert(x>=0&&y>=0);
	const static double eps = std::numeric_limits<double>::epsilon();
	const static double factor = (1+4*eps);
	return (x+y)*factor;
  }  
  double mul(const double x, const double y) {
	assert(x>=0&&y>=0);
	const static double eps = std::numeric_limits<double>::epsilon();
	const static double factor = (1+4*eps);
	return (x*y)*factor;
  }
  double div(const double x, const double y) {
	assert(x>=0&&y>=0);
	const static double eps = std::numeric_limits<double>::epsilon();
	const static double factor = (1+4*eps);
	return (x/y)*factor;
  }
  
public:
  BfmsskRootBd_double() : m_visit(false), d(1), u25(0), l25(0), v2p(0), v2m(0), v5p(0), v5m(0) {}
  void dump() const 
  { std::cerr<<"[d,u25,l25]="<<d<<","<<u25<<","<<l25<<std::endl; 
    std::cerr<<"[v2p,v2m,v5p,v5m]="<<v2p<<","<<v2m<<","<<v5p<<","<<v5m<<std::endl;
    std::cerr<<"Lg5(v5p-v5m)="<<log_5*v5p-log_5*v5m<<std::endl;
  }
  unsigned long get_deg() const
  { return d; }
  void set_visit(bool visit)
  { m_visit = visit; }
  bool get_visit() const 
  { return m_visit; }
  bool is_constructive() const
  { return true; }
  unsigned long get_bound(unsigned long deg) {
    unsigned long bound;

	//bound := l25 + u25 * (deg - 1UL) + v2m - v2p - lg5*(v5p - v5m);
//	unsigned long boundPos = static_cast<unsigned long>(std::ceil(l25 + u25 * (deg - 1UL) + v2m + log_5*v5m));
	unsigned long boundPos = static_cast<unsigned long>(std::ceil(add(add(add(l25,mul(u25,(deg - 1UL))),mul(log_5,v5m)),v2m)));
//	unsigned long boundNeg = static_cast<unsigned long>(std::floor(v2p + log_5*v5p));
	unsigned long boundNeg = static_cast<unsigned long>(std::floor(add(v2p,mul(log_5,v5p))));

	bound = boundPos>boundNeg ? boundPos-boundNeg : 0; 

    // bound cannot be zero, because of mpfr requirement.
        return (std::max)(bound, 1UL);
  }
  
  void set(long value) 			{ set(ZT(value)); }
  void set(unsigned long value) { set(ZT(value)); }
  void set(double value) 		{ set(FT(value)); }
  void set(const ZT& value) {
    if (value == 0) return;

    ZT temp(value);
	  
    v2p = temp.get_k_exp(temp, 2UL);
    v5p = temp.get_k_exp(temp, 5UL);

    v2m = v5m = 0;
    
	FT u25temp; u25temp.log2(abs(temp),54,BF_RNDU);
    u25 = u25temp.get_d(BF_RNDU);
//	u25 = temp.ceillg();
    l25 = 0;
  }
  void set(const QT& value) {
    if (value == 0) return;

    ZT num = value.num();
    ZT den = value.den();

    v2p = num.get_k_exp(num, 2UL);
    v2m = den.get_k_exp(den, 2UL);

    v5p = num.get_k_exp(num, 5UL);
    v5m = den.get_k_exp(den, 5UL);
	
	cancelPowers();

	FT u25temp; u25temp.log2(abs(num),54,BF_RNDU);
    u25 = u25temp.get_d(BF_RNDU);
	FT l25temp; l25temp.log2(abs(den),54,BF_RNDU);
    l25 = l25temp.get_d(BF_RNDU);
//    u25 = num.ceillg();
//    l25 = den.ceillg();
  }
  void set(const FT& value) {
    if (value == 0) return;

    ZT x; exp_t e = value.get_z_exp(x);
    if (e >= 0)
      v2p = e;
    else
      v2m = -e;
	
	v2p += x.get_k_exp(x,2UL);
    v5p = x.get_k_exp(x,5UL);
	v5m = 0;
	
	cancelPowers();

 	FT u25temp; u25temp.log2(abs(x),54,BF_RNDU);
    u25 = u25temp.get_d(BF_RNDU);
//   u25 = x.ceillg();
    l25 = 0;
  }
  void set(const Kernel& value) {
    if (value.sgn() >= 0)
      set(value.getLeft());
    else
      set(value.getRight());
  }
  void neg(const thisClass& child) {
    u25 = child.u25; l25 = child.l25;
    v2p = child.v2p; v2m = child.v2m;
    v5p = child.v5p; v5m = child.v5m;
  }
  void root(const thisClass& child, const unsigned long k) {
    d = k;

//    if (child.v2p + log_5*child.v5p + child.u25 >=
//        child.v2m + log_5*child.v5m + child.l25) {
    if (add(add(child.v2p,mul(log_5,child.v5p)),child.u25) >=
        add(add(child.v2m,mul(log_5,child.v5m)),child.l25)) {
      unsigned long vtilda2 = child.v2p + (k - 1UL) * child.v2m;
      v2p = vtilda2 / k;
      v2m = child.v2m;
      unsigned long vmod2;
      vmod2 = vtilda2 - k * v2p; // == vtilda2 % k
      unsigned long vtilda5 = child.v5p + (k - 1UL) * child.v5m;
      v5p = vtilda5 / k;
      v5m = child.v5m;
      unsigned long vmod5;
	  vmod5 = vtilda5 - k * v5p; // == vtilda5 % k

//      u25 = (child.u25 + (k - 1UL) * child.l25 + vmod2 + log_5*(vmod5) + 1UL) / k;
      u25 = div(add(add(add(child.u25,mul(k - 1UL,child.l25)),add(vmod2,mul(log_5,vmod5))),1),k);
      l25 = child.l25;
    } else {
      unsigned long vtilda2 = (k - 1UL) * child.v2p + child.v2m;
      v2p = child.v2p;
      v2m = vtilda2 / k;
      unsigned long vmod2;
      vmod2 = vtilda2 - k * v2m; // == vtilda2 % k
      unsigned long vtilda5 = (k - 1UL) * child.v5p + child.v5m;
      v5p = child.v5p;
      v5m = vtilda5 / k;
      unsigned long vmod5;
      vmod5 = vtilda5 - k * v5m; // == vtilda5 % k

      u25 = child.u25;
//      l25 = ((k - 1UL) * child.u25 + child.l25 + vmod2 + log_5*(vmod5) + 1UL) / k;
      l25 = div(add(add(mul(k - 1UL,child.u25),child.l25),add(add(vmod2,mul(log_5,vmod5)),1)), k);
    }
	
	//if(IsCancellingPower) cancelPowers();
  }
  void addsub(const thisClass& f, const thisClass& s) {
    v2p = (std::min)(f.v2p + s.v2m, f.v2m + s.v2p);
    v2m = f.v2m + s.v2m;
    v5p = (std::min)(f.v5p + s.v5m, f.v5m + s.v5p);
    v5m = f.v5m + s.v5m;

	// Do not need to check the signs of parameters of ceilLg5
	// They are always nonnegative!!!
	// u25 := 1 + max(f.v2p + s.v2m - v2p + lg5*(f.v5p + s.v5m - v5p) + f.u25 + s.l25,
	//				  f.v2m + s.v2p - v2p + lg5*(f.v5m + s.v5p - v5p) + s.u25 + f.l25)
	
//    u25 = 1 + (std::max)(
//  	    f.v2p + s.v2m - v2p + log_5*(f.v5p + s.v5m-v5p) + f.u25 + s.l25,
//	    f.v2m + s.v2p - v2p + log_5*(f.v5m + s.v5p-v5p) + f.l25 + s.u25); 
    u25 = 1 + (std::max)(
  	    add(add(f.v2p + s.v2m - v2p,mul(log_5,f.v5p + s.v5m-v5p)),add(f.u25,s.l25)),
	    add(add(f.v2m + s.v2p - v2p,mul(log_5,f.v5m + s.v5p-v5p)),add(f.l25,s.u25))
		); 
    l25 = f.l25 + s.l25;
/*
std::cerr << "addsub node" << std::endl;
	
std::cerr << std::endl << "dump of f" << std::endl; f.dump();
std::cerr << std::endl << "dump of s" << std::endl; s.dump();
std::cerr << std::endl << "dump of self" << std::endl; dump();
std::cerr << std::endl;
*/
	//if(IsCancellingPower) cancelPowers();
  }
  void mul(thisClass& f,thisClass& s) {
    v2p = f.v2p + s.v2p;
    v2m = f.v2m + s.v2m;
    v5p = f.v5p + s.v5p;
    v5m = f.v5m + s.v5m;
    u25 = f.u25 + s.u25;
    l25 = f.l25 + s.l25;
/*
std::cerr << "mul node" << std::endl;
	
std::cerr << std::endl << "dump of f" << std::endl; f.dump();
std::cerr << std::endl << "dump of s" << std::endl; s.dump();
std::cerr << std::endl << "dump of self" << std::endl; dump();
std::cerr << std::endl;
*/
	//if(IsCancellingPower) cancelPowers();
  }
  void div(thisClass& f, thisClass& s) { 
    v2p = f.v2p + s.v2m;
    v2m = f.v2m + s.v2p;
    v5p = f.v5p + s.v5m;
    v5m = f.v5m + s.v5p;
    u25 = f.u25 + s.l25;
    l25 = f.l25 + s.u25;
/*
std::cerr << "div node" << std::endl;
	
std::cerr << std::endl << "dump of f" << std::endl; f.dump();
std::cerr << std::endl << "dump of s" << std::endl; s.dump();
std::cerr << std::endl << "dump of self" << std::endl; dump();
std::cerr << std::endl;
*/
	//if(IsCancellingPower) cancelPowers();
  }
};

template <typename Kernel = BigFloat2>
class BfmsskRootBd {
public:
  typedef BfmsskRootBd thisClass;
  typedef typename Kernel::ZT ZT;
  typedef typename Kernel::QT QT;
  typedef typename Kernel::FT FT;
  bool			m_visit;
  unsigned long d;
  
  // u25 and l25 are maintained as upper bounds of the lg values.
  unsigned long	u25, l25;
  unsigned long v2p, v2m, v5p, v5m;

  // Utility functions    
  //return the ceil of log2(5^x)
  //x should be non-negative
  //Still not exact. do not take account the rounding effects.
  //For log_5, we need to precompute upper and lower bounds in double type.
  //Rounding up and down needs to be added.
  unsigned long ceilLg5(const unsigned long x) {
	return static_cast<unsigned long>(std::ceil(log_5*x));
//	static unsigned long ceilLg5 = ceilLg(5);
//  return ceilLg5*x;
  }
  
  //return the floor of log2(5^x)
  //x should be non-negative
  unsigned long floorLg5(const unsigned long x) {
	return static_cast<unsigned long>(std::floor(log_5*x));
//	static unsigned long floorLg5 = floorLg(5);
//  return floorLg5*x;
  }
  void cancelPowers() {
	subtractMin(v2p,v2m); subtractMin(v5p,v5m);
  }
  void subtractMin(unsigned long& x, unsigned long& y) {
	if(x > y) { x -= y; y = 0; } else { y -= x; x = 0; }
  }
public:
  BfmsskRootBd() : m_visit(false), d(1), u25(0), l25(0), v2p(0), v2m(0), v5p(0), v5m(0) {}
  void dump()  
  { std::cerr<<"[d,u25,l25]="<<d<<","<<u25<<","<<l25<<std::endl; 
    std::cerr<<"[v2p,v2m,v5p,v5m]="<<v2p<<","<<v2m<<","<<v5p<<","<<v5m<<std::endl;
    std::cerr<<"Lg5(v5p-v5m)="<<ceilLg5(v5p)-floorLg5(v5m)<<std::endl;
  }
  unsigned long get_deg() const
  { return d; }
  void set_visit(bool visit)
  { m_visit = visit; }
  bool get_visit() const 
  { return m_visit; }
  bool is_constructive() const
  { return true; }
  unsigned long get_bound(unsigned long deg) {
    unsigned long bound;

	//bound := l25 + u25 * (deg - 1UL) + v2m - v2p - lg5*(v5p - v5m);
	unsigned long boundPos = l25 + u25 * (deg - 1UL) + v2m + ceilLg5(v5m);
	unsigned long boundNeg = v2p + floorLg5(v5p);

	bound = boundPos>boundNeg ? boundPos-boundNeg : 0; 

    // bound cannot be zero, because of mpfr requirement.
        return (std::max)(bound, 1UL);
  }
  
  void set(long value) 			{ set(ZT(value)); }
  void set(unsigned long value) { set(ZT(value)); }
  void set(double value) 		{ set(FT(value)); }
  void set(const ZT& value) {
    if (value == 0) return;

    ZT temp(value);
	  
    v2p = temp.get_k_exp(temp, 2UL);
    v5p = temp.get_k_exp(temp, 5UL);

    v2m = v5m = 0;
    
    u25 = temp.ceillg();
    l25 = 0;
  }
  void set(const QT& value) {
    if (value == 0) return;

    ZT num = value.num();
    ZT den = value.den();

    v2p = num.get_k_exp(num, 2UL);
    v2m = den.get_k_exp(den, 2UL);

    v5p = num.get_k_exp(num, 5UL);
    v5m = den.get_k_exp(den, 5UL);
	
	cancelPowers();

    u25 = num.ceillg();
    l25 = den.ceillg();
  }
  void set(const FT& value) {
    if (value == 0) return;

    ZT x; exp_t e = value.get_z_exp(x);
    if (e >= 0)
      v2p = e;
    else
      v2m = -e;
	
	v2p += x.get_k_exp(x,2UL);
    v5p = x.get_k_exp(x,5UL);
	v5m = 0;
	
	cancelPowers();

    u25 = x.ceillg();
    l25 = 0;
  }
  void set(const Kernel& value) {
    if (value.sgn() >= 0)
      set(value.getLeft());
    else
      set(value.getRight());
  }
  void neg(const thisClass& child) {
    u25 = child.u25; l25 = child.l25;
    v2p = child.v2p; v2m = child.v2m;
    v5p = child.v5p; v5m = child.v5m;
  }
  void root(const thisClass& child, const unsigned long k) {
    d = k;

    if (child.v2p + ceilLg5(child.v5p) + child.u25 >=
        child.v2m + ceilLg5(child.v5m) + child.l25) {
      unsigned long vtilda2 = child.v2p + (k - 1UL) * child.v2m;
      v2p = vtilda2 / k;
      v2m = child.v2m;
      unsigned long vmod2;
      vmod2 = vtilda2 - k * v2p; // == vtilda2 % k
      unsigned long vtilda5 = child.v5p + (k - 1UL) * child.v5m;
      v5p = vtilda5 / k;
      v5m = child.v5m;
      unsigned long vmod5;
	  vmod5 = vtilda5 - k * v5p; // == vtilda5 % k

      u25 = (child.u25 + (k - 1UL) * child.l25 + vmod2 + ceilLg5(vmod5) + 1UL) / k + 1;
      l25 = child.l25;
    } else {
      unsigned long vtilda2 = (k - 1UL) * child.v2p + child.v2m;
      v2p = child.v2p;
      v2m = vtilda2 / k;
      unsigned long vmod2;
      vmod2 = vtilda2 - k * v2m; // == vtilda2 % k
      unsigned long vtilda5 = (k - 1UL) * child.v5p + child.v5m;
      v5p = child.v5p;
      v5m = vtilda5 / k;
      unsigned long vmod5;
      vmod5 = vtilda5 - k * v5m; // == vtilda5 % k

      u25 = child.u25;
      l25 = ((k - 1UL) * child.u25 + child.l25 + vmod2 + ceilLg5(vmod5) + 1UL) / k + 1;
    }
	
	//if(IsCancellingPower) cancelPowers();
  }
  void addsub(const thisClass& f, const thisClass& s) {
    v2p = (std::min)(f.v2p + s.v2m, f.v2m + s.v2p);
    v2m = f.v2m + s.v2m;
    v5p = (std::min)(f.v5p + s.v5m, f.v5m + s.v5p);
    v5m = f.v5m + s.v5m;

	// Do not need to check the signs of parameters of ceilLg5
	// They are always nonnegative!!!
	// u25 := 1 + max(f.v2p + s.v2m - v2p + lg5*(f.v5p + s.v5m - v5p) + f.u25 + s.l25,
	//				  f.v2m + s.v2p - v2p + lg5*(f.v5m + s.v5p - v5p) + s.u25 + f.l25)
	
    u25 = 1UL + (std::max)(
  	    f.v2p + s.v2m - v2p + ceilLg5(f.v5p + s.v5m-v5p) + f.u25 + s.l25,
	    f.v2m + s.v2p - v2p + ceilLg5(f.v5m + s.v5p-v5p) + f.l25 + s.u25); 
    l25 = f.l25 + s.l25;

	//if(IsCancellingPower) cancelPowers();
  }
  void mul(thisClass& f,thisClass& s) {
    v2p = f.v2p + s.v2p;
    v2m = f.v2m + s.v2m;
    v5p = f.v5p + s.v5p;
    v5m = f.v5m + s.v5m;
    u25 = f.u25 + s.u25;
    l25 = f.l25 + s.l25;

	//if(IsCancellingPower) cancelPowers();
  }
  void div(thisClass& f, thisClass& s) { 
    v2p = f.v2p + s.v2m;
    v2m = f.v2m + s.v2p;
    v5p = f.v5p + s.v5m;
    v5m = f.v5m + s.v5p;
    u25 = f.u25 + s.l25;
    l25 = f.l25 + s.u25;

	//if(IsCancellingPower) cancelPowers();
  }
};


/// BFMSS Root Bound
template <typename Kernel = BigFloat2>
class BfmssRootBd {
  typedef BfmssRootBd thisClass;
  typedef thisClass* id_rootbd_t;
  typedef typename Kernel::ZT ZT;
  typedef typename Kernel::QT QT;
  typedef typename Kernel::FT FT;
  bool    m_visit;
  unsigned long d_e;
  unsigned long u_e;
  unsigned long l_e;
public:
  BfmssRootBd() : m_visit(false), d_e(1) {}
  void dump() const 
  { std::cout<<"[d_e,u_e,l_e]="<<d_e<<","<<u_e<<","<<l_e<<std::endl; }
public:
  unsigned long get_deg()
  { return d_e; }
  void set_visit(bool visit)
  { m_visit = visit; }
  bool get_visit()
  { return m_visit; }
  bool is_constructive() const
  { return true; }
  /// BFMSS root bound could be zero, so set it to be 1 in that case
  unsigned long get_bound(unsigned long deg) const
  { return (std::max)(l_e + (deg - 1) * u_e, 1UL); }
  void set(long value)
  { u_e = ceillg(value); l_e = 0; }
  void set(unsigned long value)
  { u_e = ceillg(value); l_e = 0; }
  void set(double value)
  { set(FT(value)); }
  void set(const ZT& value)
  { u_e = value.ceillg(); l_e = 0; }
  void set(const QT& value)
  { u_e = value.num().ceillg(); l_e = value.den().ceillg(); }
  void set(const FT& value) {
    ZT x; exp_t e = value.get_z_exp(x);
    if (e >= 0) { // convert to integer
      x.mul_2exp(x, e); set(x);
    } else { // convert to rational
      QT q; q.div_2exp(x, -e); set(q);
    }
  }
  void set(const Kernel& value) {
    if (value.sgn() >= 0)
      set(value.getLeft());
    else
      set(value.getRight());
  }
  void neg(const thisClass& child) {
    u_e = child.u_e; l_e = child.l_e;
  }
  void root(thisClass& child, unsigned long k) {
    if (child.u_e >= child.l_e) {
      u_e = (child.u_e + (k-1)*child.l_e + (k-1)) / k + 1;
      l_e = child.l_e;
    } else {
      u_e = child.u_e;
      l_e = ((k-1)*child.u_e + child.l_e + (k-1)) / k + 1;
    }
    d_e = k;
  }
  void addsub(thisClass& f, thisClass& s) {
    u_e = (std::max)(f.u_e + s.l_e, f.l_e + s.u_e) + 1; 
    l_e = f.l_e + s.l_e;
  }
  void mul(thisClass& f,thisClass& s) {
    u_e = f.u_e + s.u_e;
    l_e = f.l_e + s.l_e;
  }
  void div(thisClass& f, thisClass& s) { 
    u_e = f.u_e + s.l_e;
    l_e = f.l_e + s.u_e;
  }
};

/// Minimum Root Bound (root bound class which taking minimum of two root bounds)
template <typename RootBd1, typename RootBd2>
class MinRootBd{
  typedef MinRootBd thisClass;
  RootBd1 m_rootBd1;
  RootBd2 m_rootBd2;
public:
  void dump() const 
  { m_rootBd1.dump(); m_rootBd2.dump(); }
  bool is_constructive() const 
  { return m_rootBd1.is_constructive() || m_rootBd2.is_constructive(); }
  unsigned long get_bound() const 
  { return (std::min)(m_rootBd1.get_bound(), m_rootBd2.get_bound()); }
  template <typename T> void set(const T& value) 
  { m_rootBd1.set(value), m_rootBd2.set(value); }
  void neg(const thisClass& child) 
  { m_rootBd1.neg(child), m_rootBd2.neg(child); }
  void root(const thisClass& child, unsigned long k)
  { m_rootBd1.root(child, k), m_rootBd2.root(child, k); }
  void addsub(const thisClass& f, const thisClass& s)
  { m_rootBd1.addsub(f, s), m_rootBd2.addsub(f, s); }
  void mul(const thisClass& f, const thisClass& s)
  { m_rootBd1.mul(f, s), m_rootBd2.mul(f, s); }
  void div(const thisClass& f, const thisClass& s)
  { m_rootBd1.div(f, s), m_rootBd2.div(f, s); }
};

CORE_END_NAMESPACE

#endif /*__CORE_ROOTBOUNDS_H__*/
