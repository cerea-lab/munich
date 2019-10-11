// Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
// Author(s): Vivien Mallet
//
// This file is part of AtmoData library, a tool for data processing in
// atmospheric sciences.
//
// AtmoData is developed in the INRIA - ENPC joint project-team CLIME and in
// the ENPC - EDF R&D joint laboratory CEREA.
//
// AtmoData is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// AtmoData is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// For more information, visit the AtmoData home page:
//      http://cerea.enpc.fr/polyphemus/atmodata.html


#ifndef ATMODATA_FILE_ATMODATA_HXX

#include "SeldonData.hxx"
using namespace SeldonData;

//! AtmoData namespace.
namespace AtmoData
{
}

#include "Common.hxx"

#include "Meteorology.cxx"
#include "Kz.cxx"
#include "TimeDiagnosis.cxx"
#include "Photolysis.cxx"
#include "Emissions.cxx"
#include "Deposition.cxx"
#include "CoordTransform.cxx"
#include "Transform.cxx"
#include "Polair.cxx"
#include "Aerosol.hxx"
#include "MEGAN.cxx"

#include "Format.cxx"

#include "Errors.cxx"

#define ATMODATA_FILE_ATMODATA_HXX
#endif
