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
#include "cc_inverse_multiquadric_poisson.h"

namespace scicellxx
{
 
 // ===================================================================
 // Constructor
 // ===================================================================
 CCInverseMultiquadricPoisson::CCInverseMultiquadricPoisson()
 {
  
 }
 
 // ===================================================================
 // Destructor
 // ===================================================================
 CCInverseMultiquadricPoisson::~CCInverseMultiquadricPoisson()
 {
  
 }

// ===================================================================
 /// Radial function
 // ===================================================================
 const Real CCInverseMultiquadricPoisson::Lpsi(const Real r)
 {
  const Real er2 = (epsilon*r)*(epsilon*r);
  return 4*epsilon*epsilon*std::exp(-er2)*(er2-1);
 }
 
}
