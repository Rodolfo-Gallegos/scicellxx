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
#ifndef CCFACTORYTIMESTEPPER_TPL_H
#define CCFACTORYTIMESTEPPER_TPL_H

// Include the integration methods for odes (time steppers)
#include "ac_time_stepper.h"
#include "cc_euler_method.h"
#include "cc_runge_kutta_4_method.h"
#include "cc_backward_euler_predictor_corrector_method.h"
#include "cc_adams_moulton_2_predictor_corrector_method.h"
#include "cc_backward_euler_method.h"
#include "cc_adams_moulton_2_method.h"
#include "cc_bdf_2_method.h"
#include "cc_adaptive_runge_kutta_45F_method.h"
#include "cc_adaptive_runge_kutta_45DP_method.h"

namespace scicellxx
{

 /// @class CCFactoryTimeStepper cc_factory_time_stepper.h

 /// This class implements a factory for the integration methods for
 /// ODEs (time steppers)
 template<class EQUATIONS_TYPE>
 class CCFactoryTimeStepper
 {
  
 public:
  
  /// Empty constructor
  CCFactoryTimeStepper();
  
  /// Empty destructor
  virtual ~CCFactoryTimeStepper();
  
  /// Returns the specified time stepper (integration method)
  ACTimeStepper<EQUATIONS_TYPE>* create_time_stepper(std::string time_stepper_name);
  
 protected:
  
  /// Copy constructor (we do not want this class to be
  /// copiable). Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  CCFactoryTimeStepper(const CCFactoryTimeStepper &copy)
   {
    BrokenCopy::broken_copy("CCFactoryTimeStepper");
   }
  
  /// Assignment operator (we do not want this class to be
  /// copiable. Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  void operator=(const CCFactoryTimeStepper &copy)
   {
    BrokenCopy::broken_assign("CCFactoryTimeStepper");
   }
 
 };

}
 
#endif // #ifndef CCFACTORYTIMESTEPPER_TPL_H

