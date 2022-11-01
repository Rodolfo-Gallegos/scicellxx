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
#ifndef CCBACKWARDEULERMETHOD_TPL_H
#define CCBACKWARDEULERMETHOD_TPL_H

// Include matrices and vectors
#include "../matrices/matrices.h"

#include "ac_time_stepper.h"

/// Newton's method
#include "cc_newtons_method_for_backward_euler.h"

/// Time stepper to compute the initial guess for Newton's method
///#include "cc_euler_method.h"
#include "cc_runge_kutta_4_method.h"

/// To allow the setting of the strategy for the computation of the
/// Jacobian of the ODEs
#include "ac_jacobian_and_residual_for_implicit_time_stepper.h"

namespace scicellxx
{
 
 /// @class CCBackwardEulerMethod cc_backward_euler_method.h This
 /// class implements Backward Euler's method to integrate ODE's
 template<class EQUATIONS_TYPE>
 class CCBackwardEulerMethod : public virtual ACTimeStepper<EQUATIONS_TYPE>
 {

 public:
  
  /// Constructor
  CCBackwardEulerMethod();
  
  /// Empty destructor
  virtual ~CCBackwardEulerMethod();
  
  /// Applies Backward Euler method to the given odes from the current
  /// time "t" to the time "t+h". The values of u at time t+h will be
  /// stored at index k (default k = 0).
  void time_step(EQUATIONS_TYPE &odes, const Real h, const Real t,
                 CCData &u, const unsigned k = 0);
  
  /// Set the strategy for the computation of the Jacobian of the ODEs (if known)
  inline void set_strategy_for_odes_jacobian(ACJacobianAndResidualForImplicitTimeStepper<EQUATIONS_TYPE> *jacobian_strategy_for_odes_pt)
  {Newtons_method.set_strategy_for_odes_jacobian(jacobian_strategy_for_odes_pt);}
  
 protected:
  
  /// Copy constructor (we do not want this class to be
  /// copiable). Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
 CCBackwardEulerMethod(const CCBackwardEulerMethod &copy)
  : ACTimeStepper<EQUATIONS_TYPE>()
   {
    BrokenCopy::broken_copy("CCBackwardEulerMethod");
   }
  
  /// Assignment operator (we do not want this class to be
  /// copiable. Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  void operator=(const CCBackwardEulerMethod &copy)
   {
    BrokenCopy::broken_assign("CCBackwardEulerMethod");
   }
  
  // Newton's method for backward Euler
  CCNewtonsMethodForBackwardEuler<EQUATIONS_TYPE> Newtons_method;
  
  // The time stepper used to compute the initial guess
  //CCEulerMethod Time_stepper_initial_guess;

  // NOTE: We decided to use a RK4 method as the initial guess method
  // to reduce accuracy errors given by Euler's methods
  
  // The time stepper used to compute the initial guess
  CCRK4Method<EQUATIONS_TYPE> Time_stepper_initial_guess;
  
 };
 
}

#endif // #ifndef CCBACKWARDEULERMETHOD_TPL_H
