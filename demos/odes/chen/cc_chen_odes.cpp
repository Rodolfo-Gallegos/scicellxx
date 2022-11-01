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
/** \file This file implements the CCChenODEs class
 */
#include "cc_chen_odes.h"

namespace scicellxx
{
 // ===================================================================
 // Constructor, sets the number of odes.
 // ===================================================================
 CCChenODEs::CCChenODEs(Real _a, Real _b, Real _c)
  : ACODEs(3)
 {  
  // Configure parameters
  a = _a;
  b = _b;
  c = _c;
 }
 
 // ===================================================================
 // Empty destructor
 // ===================================================================
 CCChenODEs::~CCChenODEs()
 { }
 
 // ===================================================================
 /// Evaluates the system of odes at time "t".
 // ===================================================================
 void CCChenODEs::evaluate_time_derivatives(const Real t,
                                            CCData &u,
                                            CCData &dudt,
                                            const unsigned k)
 {
  // -----------------
  /// \frac{du_{1}}{dt} = a*(u_{2} - u_{1})
  /// \frac{du_{2}}{dt} = (c-a)*u_{1} - u_{1}*u_{3} + c*u_{2}
  /// \frac{du_{3}}{dt} = u_{1}*u_{2} - b*u_{3}
  
  dudt(0) = a*(u(1,k) - u(0,k));
  dudt(1) = (c-a)*u(0,k) - (u(0,k)*u(2,k)) + (c*u(1,k));
  dudt(2) = (u(0,k)*u(1,k)) - (b*u(2,k));
  
 }
 
}
