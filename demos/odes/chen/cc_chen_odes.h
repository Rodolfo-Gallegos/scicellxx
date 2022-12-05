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
#ifndef CCCHENODES_H
#define CCCHENODES_H

// Include SciCell++ libraries
#include "../../../src/scicellxx.h"

namespace scicellxx
{

 /// \class CCChenODEs cc_chen_odes.h
    
 /// This class implements the Chen ODEs
 ///
 /// \frac{du_{1}}{dt} = a*(u_{2} - u_{1})
 /// \frac{du_{2}}{dt} = (c-a)*u_{1} - u_{1}*u_{3} + c*u_{2}
 /// \frac{du_{3}}{dt} = u_{1}*u_{2} - b*u_{3}
 class CCChenODEs : public virtual ACODEs
 {
 
 public:
  
  /// Constructor
  CCChenODEs(Real _a, Real _b, Real _c);
  
  /// Empty destructor
  virtual ~CCChenODEs();
  
  /// Evaluates the system of odes at time 't', using the history
  /// values of u at index k
  void evaluate_time_derivatives(const Real t, CCData &u, CCData &dudt, const unsigned k = 0);
  
 protected:
  
  // Specific ODEs parameters
  // ------------------------------------
  Real a;
  Real b;
  Real c;
  
 };
 
}

#endif // #ifndef CCCHENODES_H
