/****************************************************************************
 * Timer.h -- A C++ class providing timing
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
 * $Id: Timer.h,v 1.7 2010/11/23 17:58:37 exact Exp $
 ***************************************************************************/
#ifndef __CORE_TIMER_H__
#define __CORE_TIMER_H__

#include <CORE/Config.h>
#include <ctime>
#ifndef _MSC_VER
#include <sys/resource.h>
#endif

CORE_BEGIN_NAMESPACE

/// \class Timer Timer.h
/// \brief a timer using clock() function, return in seconds, less precise
class Timer {
public:
  Timer() : startClock(0), clocks(0) {}
  /// start timer
  void start() { startClock = clock(); }
  /// stop timer
  void stop() { clocks = clock() - startClock; }
  /// return in clocks
  long getClocks() { return clocks; }
  /// return in seconds
  float getSeconds() { return (float)clocks/CLOCKS_PER_SEC; }
private:
  long startClock;
  long clocks;
};

/// \class Timer2 Timer.h
/// \brief a timer using getrusage() function, return in mseconds, more precise
class Timer2 {
public:
  Timer2() : startClock(0), clocks(0) {}
  /// start timer
  void start() { startClock = cputime(); }
  /// stop timer
  void stop() { clocks = cputime() - startClock; }
  /// return in microseconds
  long get_mseconds() { return clocks; }
  /// return in seconds
  float get_seconds() { return (float)clocks/1000; }
  /// return in seconds
  float getSeconds() { return (float)clocks/1000; }
private:
  long startClock;
  long clocks;
  static long cputime () {
#ifndef _MSC_VER
    struct rusage rus; getrusage (0, &rus);
    return rus.ru_utime.tv_sec * 1000 + rus.ru_utime.tv_usec / 1000;
#else 
    return 0;
#endif
  }
};

CORE_END_NAMESPACE

#endif /*__CORE_TIMER_H__*/
