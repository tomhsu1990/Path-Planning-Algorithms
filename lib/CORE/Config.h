/****************************************************************************
 * Config.h -- Configuration Macros
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
 * $Id: Config.h,v 1.3 2006/03/03 17:19:58 exact Exp $
 ***************************************************************************/
#ifndef __CORE_CONFIG_H__
#define __CORE_CONFIG_H__

// define version number
#define CORE_VERSION            2
#define CORE_VERSION_MINOR      0
#define CORE_VERSION_PATCHLEVEL 0

// macros for defining namespace
#define CORE_BEGIN_NAMESPACE    namespace CORE {
#define CORE_END_NAMESPACE      };
#define CORE_NS                 CORE

// diable old names
//#define CORE_DISABLE_OLDNAMES

// disable reference counting
//#define CORE_DISABLE_REFCOUNTING

// disable memory pool
//#define CORE_DISABLE_MEMPOOL

// disable debug
//#define CORE_DISABLE_DEBUG

// debug filter
//#define CORE_DEBUG_FILTER 1

// debug root bound
//#define CORE_DEBUG_ROOTBOUND 1

// disable automatic namespace
//#define CORE_NO_AUTOMATIC_NAMESPACE

#endif /*__CORE_CONFIG_H__*/
