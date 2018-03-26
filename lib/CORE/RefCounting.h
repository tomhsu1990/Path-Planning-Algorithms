#ifndef __CORE_REFCOUNTING_H__
#define __CORE_REFCOUNTING_H__

/*
 *     1. This file defines two templated classes:
 *               RCRepImpl<class N>
 *     to create Reps of the class N.  The basic functions provided by
 *     this class is reference counting.   The other class is
 *               RCImpl<class T>
 *     for implementing the envelop-letter paradigm for a class whose Rep
 *     is the class T.  So, T is the "letter", and RCImpl<T> the "envelop".
 *
 *     2. All Rep classes (BigIntRep, BigFloatRep, BigRatRep, ExprRep, etc)
 *     are derived from RCRepImpl<N>.  E.g.,
 *
 *         class BigRatRep : public RCRepImp<BigRatRep> {
 *         ...
 *         }
 *     (Note the recursive use of "BigRatRep").
 *
 *     3. All Number classes (BigInt, BigFloat, BigRat, Expr, etc)
 *     are derived from RCImpl<T>.  E.g.
 *
 *         typedef RCImpl<BigRatRep> RCBigRat;
 *         class BigRat : public RCBigRat {
 *         ...
 *         }
 */

CORE_BEGIN_NAMESPACE

template <class T>
class RcRepImpl {
public:
  RcRepImpl() : _ref_counter(1) {}
  int get_rc() const { return _ref_counter; }
  void inc_rc() { ++_ref_counter; }
  void dec_rc() { if (--_ref_counter == 0) delete static_cast<T*>(this); }
private:
  int _ref_counter;
};

template <class T>
class RcImpl {
public:
  RcImpl(T* p) : _rep(p) {}
  int get_rc() const { return _rep->get_rc(); }
  void make_copy() {
    if (get_rc() > 1) { // if the object is shared, then clone it!
      _rep->dec_rc();
      _rep = _rep ? new T(*_rep) : 0;
    }
  }
protected:
  T* _rep;
};

CORE_END_NAMESPACE

#endif /*__CORE_REFCOUNTING_H__*/
