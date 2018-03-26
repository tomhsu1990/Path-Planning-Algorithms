/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CORE (http://cs.nyu.edu/exact/core/); you may
 * redistribute it under the terms of the Q Public License version 1.0.
 * See the file LICENSE.QPL distributed with CORE.
 *
 * Licensees holding a valid commercial license may use this file in
 * accordance with the commercial license agreement provided with the
 * software.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * File: CoreIO.h
 *
 * 	Most of these code are to convert from mpfr's number IO.
 *
 * Written by 
 *       Zilin Du <zilin@cs.nyu.edu>
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $Source: /home/exact/cvsroot/exact/corelib2/src/CoreIo.cpp,v $
 * $Revision: 1.1 $ $Date: 2006/11/10 21:09:28 $
 ***************************************************************************/


#ifndef __CORE_COREIO_H__
#define __CORE_COREIO_H__

#include <CORE/Config.h>
#include <CORE/CoreDefs.h>
#include <CORE/CoreAux.h>
#include <iostream>
#include <iomanip>

CORE_BEGIN_NAMESPACE

extern rnd_t def_output_rounding_mode;
extern int def_output_base;

////////////////////////////////////////
// output getters

// get output rounding mode
inline rnd_t get_output_rounding_mode()
{ return def_output_rounding_mode; }

// get output precision
inline int get_output_precision(const std::ostream& os = std::cout)
{ return os.precision(); }

// get output base (use def_output_base first)
inline int get_output_base(const std::ostream& os = std::cout) {
  if (def_output_base != 0) return def_output_base;
  else if (os.flags() & std::ios::oct) return 8;
  else if (os.flags() & std::ios::hex) return 16;
  else return 10;
}

// get output fmt 
inline int get_output_fmt(const std::ostream& os = std::cout) {
  if (os.flags() & std::ios::fixed) return 1;
  else if (os.flags() & std::ios::scientific) return 2;
  else return 0;
}

// get output showpoint (decimal point)
inline bool get_output_showpoint(const std::ostream& os = std::cout)
{ return os.flags() & std::ios::showpoint; }

// get output showpos (positive sign)
inline bool get_output_showpos(const std::ostream& os = std::cout)
{ return os.flags() & std::ios::showpos; }

// get output uppercase 
inline bool get_output_uppercase(const std::ostream& os = std::cout)
{ return os.flags() & std::ios::uppercase; }

////////////////////////////////////////
// output setters

// set default output rounding mode
inline void set_output_rounding_mode(rnd_t rnd)
{ def_output_rounding_mode = rnd; }

// set output precision
inline void set_output_precision(int prec, std::ostream& os = std::cout)
{ os.precision(prec); }

// set dec output
inline void set_dec_output(std::ostream& os = std::cout)
{ os.setf(std::ios::dec, std::ios::basefield); }
// set hex output
inline void set_hex_output(std::ostream& os = std::cout)
{ os.setf(std::ios::hex, std::ios::basefield); }
// set oct output
inline void set_oct_output(std::ostream& os = std::cout)
{ os.setf(std::ios::oct, std::ios::basefield); }
// set output base
inline void set_output_base(const std::ios::fmtflags& f, std::ostream& os = std::cout)
{ os.setf(f, std::ios::basefield); }
// set output base
inline void set_output_base(int base, std::ostream& os = std::cout) {
  if (base == 8) { set_oct_output(os); def_output_base = 0; }
  else if (base == 10) { set_dec_output(os); def_output_base = 0; }
  else if (base == 16) { set_hex_output(os); def_output_base = 0; }
  else { set_dec_output(os); def_output_base = base; }
}


// set default output fmt
inline void set_default_output(std::ostream& os = std::cout)
{ os.unsetf(std::ios::fixed); os.unsetf(std::ios::scientific);}
// set fixed output fmt
inline void set_fixed_output(std::ostream& os = std::cout)
{ os.setf(std::ios::fixed, std::ios::floatfield);}
// set scientific output fmt
inline void set_scientific_output(std::ostream& os = std::cout)
{ os.setf(std::ios::scientific, std::ios::floatfield);}
// set output fmt
inline void set_output_fmt(const std::ios::fmtflags& f, std::ostream& os = std::cout)
{ os.setf(f, std::ios::floatfield); }
// set output fmt
inline void set_output_fmt(int fmt, std::ostream& os = std::cout) {
  if (fmt == 1) set_fixed_output(os);
  else if (fmt == 2) set_scientific_output(os);
  else set_default_output(os);
}

// show decimal point
inline void show_dec_point(std::ostream& os = std::cout)
{ os.setf(std::ios::showpoint); }
// hide decimal point
inline void hide_dec_point(std::ostream& os = std::cout)
{ os.unsetf(std::ios::showpoint); }
// set output decimal point
inline void set_output_showpoint(bool showpoint, std::ostream& os = std::cout) 
{ if (showpoint) show_dec_point(os); else hide_dec_point(os); }

// show positive sign
inline void show_pos(std::ostream& os = std::cout)
{ os.setf(std::ios::showpos); }
// hide positive sign
inline void hide_pos(std::ostream& os = std::cout)
{ os.unsetf(std::ios::showpos); }
// set output positve sign
inline void set_output_showpos(bool showpos, std::ostream& os = std::cout) 
{ if (showpos) show_pos(os); else hide_pos(os); }

#ifndef CORE_DISABLE_OLDNAMES
inline void setPositionalFormat(std::ostream& os = std::cout)
{ set_fixed_output(os); }
inline void setScientificFormat(std::ostream& os = std::cout)
{ set_scientific_output(os); }
inline void setDefaultOutputDigits(long p, std::ostream& o = std::cout)
{ set_output_precision(p, o); }
/// This is the simple interface for setting both the output precision
/// and the precision for internal computation to the same precision in digits.
inline void set_def_output_and_computing_digits(long p, std::ostream& o = std::cout)
{
  if (p == CORE_INFTY) { // this really should not happen...
    core_error("DefOutputDigits cannot be CORE_INFTY", __FILE__, __LINE__, false);
    p = bits2digits(defRelPrec);
  }
  setDefaultOutputDigits(p, o);
  setDefaultRelPrecision(digits2bits(p));// be careful!  p can be infinite for input prec
}
inline long getDefaultOutputDigits()
{ return get_output_precision(); }
#endif

CORE_END_NAMESPACE

#endif /*__CORE_COREIO_H__*/
