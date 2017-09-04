// Copyright (C) 2003-2007, INRIA
// Author(s): Vivien Mallet
//
// This file is part of SeldonData library, used for data processing.
//
// SeldonData is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// SeldonData is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
// more details.
//
// For more information, visit the SeldonData home page:
//      http://vivienmallet.net/lib/seldondata/

#ifndef FILE_SELDONDATA_SELDONDATA_HXX

// Blitz++.
#include <blitz/array.h>
using namespace blitz;

// Talos.
#include <Talos.hxx>
using namespace Talos;

#include <string>
#include <sstream>
#include <exception>
#include <stdexcept>

//! SeldonData namespace.
namespace SeldonData
{
}


//////////////////
// DEBUG LEVELS //
//////////////////

#ifdef SELDONDATA_DEBUG_LEVEL_4
#ifndef SELDONDATA_DEBUG_LEVEL_3
#define SELDONDATA_DEBUG_LEVEL_3
#endif
#endif

#ifdef SELDONDATA_DEBUG_LEVEL_3
#ifndef SELDONDATA_DEBUG_CHECK_INDICES
#define SELDONDATA_DEBUG_CHECK_INDICES
#endif
#ifndef SELDONDATA_DEBUG_LEVEL_2
#define SELDONDATA_DEBUG_LEVEL_2
#endif
#endif

#ifdef SELDONDATA_DEBUG_LEVEL_2
#ifndef SELDONDATA_DEBUG_CHECK_DIMENSIONS
#define SELDONDATA_DEBUG_CHECK_DIMENSIONS
#endif
#ifndef SELDONDATA_DEBUG_LEVEL_1
#define SELDONDATA_DEBUG_LEVEL_1
#endif
#endif

#ifdef SELDONDATA_DEBUG_LEVEL_1
#ifndef SELDONDATA_DEBUG_CHECK_MEMORY
#define SELDONDATA_DEBUG_CHECK_MEMORY
#endif
#ifndef SELDONDATA_DEBUG_CHECK_IO
#define SELDONDATA_DEBUG_CHECK_IO
#endif
#ifndef SELDONDATA_DEBUG_LEVEL_0
#define SELDONDATA_DEBUG_LEVEL_0
#endif
#endif


#include "Errors.hxx"


////////////
// MACROS //
////////////

//! Convenient structure to catch exceptions.
/*!
  Use TRY and END to catch exceptions thrown by SeldonData:

  [...]
  TRY;

  Code that could throw an exception.

  END;

  [...]

  If an exception is caught, explanations are displayed to identify
  the problem (name of the involved function and comments).

*/

#ifdef TRY
#undef TRY
#endif
#define TRY try					\
    {

#ifdef END
#undef END
#endif
#define END							\
  }								\
    catch(SeldonData::Error& Err)				\
      {								\
	Err.What();						\
	return 1;						\
      }								\
    catch (std::exception& Err)					\
      {								\
	cout << "C++ exception: " << Err.what() << endl;	\
	return 1;						\
      }								\
    catch (std::string& str)					\
      {								\
	cout << str << endl;					\
	return 1;						\
      }								\
    catch (const char* str)					\
      {								\
	cout << str << endl;					\
	return 1;						\
      }								\
    catch(...)							\
      {								\
	cout << "Unknown error..." <<endl;			\
	return 1;						\
      }

// To get 'min' and 'max' functions.
#include <algorithm>

#include "Function.hxx"
#include "Grid.cxx"
#include "Data.cxx"
#include "Format.cxx"
#include "Functions.cxx"

#define FILE_SELDONDATA_SELDONDATA_HXX
#endif
