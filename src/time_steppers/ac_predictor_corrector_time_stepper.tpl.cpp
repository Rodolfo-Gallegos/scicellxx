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
#include "ac_predictor_corrector_time_stepper.tpl.h"

namespace scicellxx
{

 // ===================================================================
 // Empty constructor
 // ===================================================================
 template<class EQUATIONS_TYPE>
 ACPredictorCorrectorTimeStepper<EQUATIONS_TYPE>::ACPredictorCorrectorTimeStepper()
  : ACTimeStepper<EQUATIONS_TYPE>(),
    Maximum_iterations(DEFAULT_PREDICTOR_CORRECTOR_TIME_STEPPER_MAXIMUM_ITERATIONS),
    Maximum_tolerance(DEFAULT_PREDICTOR_CORRECTOR_TIME_STEPPER_MAXIMUM_TOLERANCE),
    Minimum_tolerance(DEFAULT_PREDICTOR_CORRECTOR_TIME_STEPPER_MINIMUM_TOLERANCE),
    Output_messages(false),
    Fixed_number_of_iterations(false),
    Perform_final_evaluation(false)
 {
  // Sets the default new time step strategy
  set_default_configuration();
 }
 
 // ===================================================================
 // Empty destructor
 // ===================================================================
 template<class EQUATIONS_TYPE>
 ACPredictorCorrectorTimeStepper<EQUATIONS_TYPE>::~ACPredictorCorrectorTimeStepper()
 { 

 }

 // ===================================================================
 // Set the default configuration
 // ===================================================================
 template<class EQUATIONS_TYPE>
 void ACPredictorCorrectorTimeStepper<EQUATIONS_TYPE>::set_default_configuration()
 {
  set_default_maximum_iterations();
  set_default_maximum_tolerance(); 
  set_default_minimum_tolerance();
  
  // By default output messages are disabled
  disable_output_messages();
  // By default disable fixed number of operations
  disable_fixed_number_of_corrections();
  // By default disable the last evaluation
  disable_final_evaluation();
 }
 
}
