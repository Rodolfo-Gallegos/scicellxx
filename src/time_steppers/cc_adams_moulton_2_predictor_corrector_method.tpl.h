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
#ifndef CCADAMSMOULTON2PREDICTORCORRECTORMETHOD_TPL_H
#define CCADAMSMOULTON2PREDICTORCORRECTORMETHOD_TPL_H

#include "ac_predictor_corrector_time_stepper.h"

// Include matrices and vectors
#include "../matrices/matrices.h"

#if 0
// Time stepper to compute the initial guess for Newton's method
#include "cc_runge_kutta_4_method.h"
#endif // #if 0

namespace scicellxx
{

 /// @class CCAdamsMoulton2PCMethod
 /// cc_adams_moulton_2_predictor_corrector_method.h This class
 /// implements Adams-Moulton 2 method or Trapezoidal rule as a
 /// predictor corrector to integrate ODE's.
 template<class EQUATIONS_TYPE>
 class CCAdamsMoulton2PCMethod : public virtual ACPredictorCorrectorTimeStepper<EQUATIONS_TYPE>
 {
 
 public:

  /// Constructor
  CCAdamsMoulton2PCMethod();
  
  /// Empty destructor
  virtual ~CCAdamsMoulton2PCMethod();
  
  /// Applies Adams-Moulton 2 method implemented as
  /// Predictor-Corrector to the given odes from the current time "t"
  /// to the time "t+h".  The values of u at time t+h will be stored
  /// at index k (default k = 0).
  void time_step(EQUATIONS_TYPE &odes, const Real h, const Real t,
                 CCData &u, const unsigned k = 0);
  
 protected:
 
  /// Copy constructor (we do not want this class to be
  /// copiable). Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  CCAdamsMoulton2PCMethod(const CCAdamsMoulton2PCMethod &copy)
     : ACPredictorCorrectorTimeStepper<EQUATIONS_TYPE>()
   //   : ACTimeStepperForODEs() -- try this HERE
  //   : ACTimeStepper()
   {
    BrokenCopy::broken_copy("CCAdamsMoulton2PCMethod");
   }
 
  /// Assignment operator (we do not want this class to be
  /// copiable. Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  void operator=(const CCAdamsMoulton2PCMethod &copy)
   {
    BrokenCopy::broken_assign("CCAdamsMoulton2PCMethod");
   }
  
#if 0
  // The time stepper used to compute the initial guess
  CCRK4Method Time_stepper_initial_guess;
#endif // #if 0
  
 };

}
 
#endif // #ifndef CCADAMSMOULTON2PREDICTORCORRECTORMETHOD_TPL_H
