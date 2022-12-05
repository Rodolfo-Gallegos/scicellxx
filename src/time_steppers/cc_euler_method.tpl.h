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
#ifndef CCEULERMETHOD_TPL_H
#define CCEULERMETHOD_TPL_H

#include "ac_time_stepper.h"

namespace scicellxx
{

 /// @class CCEulerMethod cc_euler_method.tpl.h
 /// This class implements Euler's method to integrate ODE's
 template<class EQUATIONS_TYPE>
 class CCEulerMethod : public virtual ACTimeStepper<EQUATIONS_TYPE>
 {
 
 public:

  /// Constructor
  CCEulerMethod();
  
  /// Empty destructor
  virtual ~CCEulerMethod();
  
  /// Applies Eulers method to the given odes from the current time
  /// "t" to the time "t+h". The values of u at time t+h will be
  /// stored at index k (default k = 0).
  void time_step(EQUATIONS_TYPE &odes, const Real h, const Real t,
                 CCData &u, const unsigned k = 0);
  
 protected:
  
  /// Copy constructor (we do not want this class to be
  /// copiable). Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  CCEulerMethod(const CCEulerMethod &copy)
   : ACTimeStepper<EQUATIONS_TYPE>()
   {
    BrokenCopy::broken_copy("CCEulerMethod");
   }
  
  /// Assignment operator (we do not want this class to be
  /// copiable. Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  void operator=(const CCEulerMethod &copy)
   {
    BrokenCopy::broken_assign("CCEulerMethod");
   }

 };

}
 
#endif // #ifndef CCEULERMETHOD_TPL_H
