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
#ifndef CCODESLOTKAVOLTERRAODES_H
#define CCODESLOTKAVOLTERRAODES_H

// Include general/common includes, utilities and initialisation
#include "../../../src/general/common_includes.h"
#include "../../../src/general/utilities.h"
#include "../../../src/general/initialise.h"
// The class used to store the values of u and dudt
#include "../../../src/data_structures/cc_data.h"
// The class implementing the interfaces for the ODEs
#include "../../../src/data_structures/ac_odes.h"

namespace scicellxx
{

 /// \class CCLotkaVolterraODEs cc_lotka_volterra_odes.h
    
 /// This class implements the simplest version of the Lotka-Volterra
 /// equations
 ///
 /// \frac{du_{1}}{dt} = a*u_{1} - b*u_{1}*u_{2}
 /// \frac{du_{2}}{dt} = -c*u_{2} + d*u_{1}*u_{2}
 class CCLotkaVolterraODEs : public virtual ACODEs
 {
 
 public:
  
  /// Constructor
  CCLotkaVolterraODEs(Real _a, Real _b, Real _c, Real _d);
  
  /// Empty destructor
  virtual ~CCLotkaVolterraODEs();
  
  /// Evaluates the system of odes at time 't', using the history
  /// values of u at index k
  void evaluate_derivatives(const Real t, CCData &u, CCData &dudt, const unsigned k = 0);
  
 protected:
  
  /// Copy constructor (we do not want this class to be
  /// copiable). Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
 CCLotkaVolterraODEs(const CCLotkaVolterraODEs &copy)
  : ACODEs(copy)
  {
   BrokenCopy::broken_copy("CCLotkaVolterraODEs");
  }
  
  /// Assignment operator (we do not want this class to be
  /// copiable. Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  void operator=(const CCLotkaVolterraODEs &copy)
   {
    BrokenCopy::broken_assign("CCLotkaVolterraODEs");
   }
  
  // Specific ODEs parameters
  // ------------------------------------
  // The prey grow rate
  Real a;
  // The predator-prey interaction rate on the prey death
  Real b;
  // The predator death rate
  Real c;
  // The predator-prey interaction rate on the predator grow
  Real d;
  
 };
 
}

#endif // #ifndef CCODESLOTKAVOLTERRAODES_H
