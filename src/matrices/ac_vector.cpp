//LIC// ====================================================================
//LIC// This file forms part of SciCell++, an object-oriented, 
//LIC// framework for the the simulation of biological and physical
//LIC// phenomena modelled as continuous or discrete processes.
//LIC// 
//LIC// You can find a copy at https://github.com/tachidok/scicellxx
//LIC// 
//LIC//    Version 0.6.0
//LIC//
//LIC// 31/10/2022
//LIC// 
//LIC// SciCell++ Copyright (C) 2016-2022 Julio César Pérez Sansalvador
//LIC// 
//LIC// This framework is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU GENERAL PUBLIC LICENSE
//LIC// published by the Free Software Foundation; either version 3 of
//LIC// the License, or (at your option) any later version.
//LIC// 
//LIC// This framework is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//LIC// GNU GENERAL PUBLIC LICENSE for more details.
//LIC// 
//LIC// You should have received a copy of the GNU GENERAL PUBLIC LICENSE
//LIC// along with this framework; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The author may be contacted at jcp.sansalvador@inaoep.mx
//LIC// 
//LIC// ====================================================================
// IN THIS FILE: Implementation of an abstract class to represent
// vectors

#include "ac_vector.h"

namespace scicellxx
{

 // ===================================================================
 // Empty constructor
 // ===================================================================
 ACVector::ACVector() 
  : NValues(0), Is_own_memory_allocated(false), Delete_vector(true), Is_column_vector(true)
 { }
 
 // ===================================================================
 // Constructor to create an n size zero vector
 // ===================================================================
 ACVector::ACVector(const unsigned long n, bool is_column_vector)
  : NValues(n), Is_own_memory_allocated(false), Delete_vector(true),
    Is_column_vector(is_column_vector)
 { }
 
 // ===================================================================
 // Destructor
 // ===================================================================
ACVector::~ACVector()
{ }
 
}
