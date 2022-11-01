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
/// IN THIS FILE: The definition of the concrete class
/// CCNewtonsMethodForAdamsMoulton2 to apply Newton's methods for
/// Adams-Moulton 2 or Trapezoidal Rule method

/// Check whether the class has been already defined
#ifndef CCNEWTONSMETHODFORADAMSMOULTON2_TPL_H
#define CCNEWTONSMETHODFORADAMSMOULTON2_TPL_H

#include "../general/general.h"

#include "cc_jacobian_and_residual_for_adams_moulton_2.h"
#include "ac_newtons_method_for_implicit_time_stepper.h"

namespace scicellxx
{
 
 /// A concrete class that implements Newton's method for time stepping
 /// methods
 template<class EQUATIONS_TYPE>
 class CCNewtonsMethodForAdamsMoulton2 : public virtual ACNewtonsMethodForImplicitTimeStepper<EQUATIONS_TYPE>
 {
  
 public:
   
  /// Constructor
  CCNewtonsMethodForAdamsMoulton2();
   
  /// Empty destructor
  ~CCNewtonsMethodForAdamsMoulton2();
   
 protected:
  
  /// Performs actions before initial converngence check
  void actions_before_initial_convergence_check();
   
 private:
  
  /// Copy constructor (we do not want this class to be copiable because
  /// it contains dynamically allocated variables, A in this
  /// case). Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  CCNewtonsMethodForAdamsMoulton2(const CCNewtonsMethodForAdamsMoulton2 &copy)
   : ACNewtonsMethodForImplicitTimeStepper<EQUATIONS_TYPE>()
   {
    BrokenCopy::broken_copy("CCNewtonsMethodForAdamsMoulton2");
   }
   
  /// Copy constructor (we do not want this class to be copiable because
  /// it contains dynamically allocated variables, A in this
  /// case). Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  void operator=(const CCNewtonsMethodForAdamsMoulton2 &copy)
   {
    BrokenCopy::broken_assign("CCNewtonsMethodForAdamsMoulton2");
   }
   
  /// The Jacobian and Residual computation strategy that implements
  /// Adams-Moulton or Trapezoidal rule method
  CCJacobianAndResidualForAdamsMoulton2<EQUATIONS_TYPE> Jacobian_and_residual_for_adams_moulton_2;
   
 };
 
}

#endif /// #ifndef CCNEWTONSMETHODFORADAMSMOULTON2_TPL_H

