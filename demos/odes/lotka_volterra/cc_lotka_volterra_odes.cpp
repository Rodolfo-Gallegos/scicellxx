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
/** \file This file implements the CCLotkaVolterraODEs class
 */
#include "cc_lotka_volterra_odes.h"

namespace scicellxx
{
 // ===================================================================
 // Constructor, sets the number of odes.
 // ===================================================================
 CCLotkaVolterraODEs::CCLotkaVolterraODEs(Real _a, Real _b,
                                          Real _c, Real _d)
  : ACODEs(2)
 {  
  // Configure parameters
  a = _a;
  b = _b;
  c = _c;
  d = _d;
 }
 
 // ===================================================================
 // Empty destructor
 // ===================================================================
 CCLotkaVolterraODEs::~CCLotkaVolterraODEs()
 { }
 
 // ===================================================================
 // Evaluates the system of odes at time "t", using the history values
 // of u at index k
 // ===================================================================
 void CCLotkaVolterraODEs::evaluate_time_derivatives(const Real t,
                                                     CCData &u,
                                                     CCData &dudt,
                                                     const unsigned k)
 {
  // -----------------
  // u(0,k) Number of prey at history index k
  // u(1,k) Number of predators at history index k
  // -----------------
  // dudt(0) Rate of change or prey with respect to time
  // dudt(1) Rate of change of predators with respecto to time
  
  /// \frac{du_{1}}{dt} = a*u_{1} - b*u_{1}*u_{2}
  /// \frac{du_{2}}{dt} = -c*u_{2} + d*u_{1}*u_{2}
  dudt(0) = a*u(0,k) - b*u(0,k)*u(1,k);
  dudt(1) = -c*u(1,k) + d*u(0,k)*u(1,k);
  
 }
 
}
