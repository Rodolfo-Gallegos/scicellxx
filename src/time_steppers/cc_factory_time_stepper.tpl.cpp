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
#include "cc_factory_time_stepper.tpl.h"

namespace scicellxx
{

 // ===================================================================
 /// Empty constructor
 // ===================================================================
 template<class EQUATIONS_TYPE>
 CCFactoryTimeStepper<EQUATIONS_TYPE>::CCFactoryTimeStepper()
 { 

 }

 // ===================================================================
 /// Empty destructor
 // ===================================================================
 template<class EQUATIONS_TYPE>
 CCFactoryTimeStepper<EQUATIONS_TYPE>::~CCFactoryTimeStepper()
 { 

 }

 // ===================================================================
 /// Returns the specified time stepper (integration method)
 // ===================================================================
 template<class EQUATIONS_TYPE>
 ACTimeStepper<EQUATIONS_TYPE>* CCFactoryTimeStepper<EQUATIONS_TYPE>::create_time_stepper(std::string time_stepper_name)
 {
  // Get the string and change it to lower case 
  std::transform(time_stepper_name.begin(), time_stepper_name.end(),
                 time_stepper_name.begin(), ::tolower);
  
  // ------------------------------------------------------
  // Check what time stepper method we need to create
  // ------------------------------------------------------
  // Euler method
  if (time_stepper_name.compare("euler")==0)
   {
    return new CCEulerMethod<EQUATIONS_TYPE>();
   }
  // Runge-Kutta 4 method
  else if (time_stepper_name.compare("rk4")==0)
   {
    return new CCRK4Method<EQUATIONS_TYPE>();
   }
  // Backward-Euler as Predictor-Corrector method
  else if (time_stepper_name.compare("bepc")==0)
   {
    return new CCBackwardEulerPCMethod<EQUATIONS_TYPE>();
   }
  // Adams-Moulton 2 as Predictor-Corrector method
  else if (time_stepper_name.compare("am2pc")==0)
   {
    return new CCAdamsMoulton2PCMethod<EQUATIONS_TYPE>();
   }
  // Backward Euler method
  else if (time_stepper_name.compare("bdf1")==0)
   {
    return new CCBackwardEulerMethod<EQUATIONS_TYPE>();
   }
  // Adams-Moulton 2 or Trapezoidal Rule method
  else if (time_stepper_name.compare("am2")==0)
   {
    return new CCAdamsMoulton2Method<EQUATIONS_TYPE>();
   }
  // BDF 2 method
  else if (time_stepper_name.compare("bdf2")==0)
   {
    return new CCBDF2Method<EQUATIONS_TYPE>();
   }
  // Runge-Kutta 4(5) Fehlberg method
  else if (time_stepper_name.compare("rk45f")==0)
   {
    return new CCAdaptiveRK45FMethod<EQUATIONS_TYPE>();
   }
  // Runge-Kutta 4(5) Dormand-Prince method
  else if (time_stepper_name.compare("rk45dp")==0)
   {
    return new CCAdaptiveRK45DPMethod<EQUATIONS_TYPE>();
   }
  else
   {
    std::ostringstream error_message;
    error_message << "The time stepper (integration method) you want to use is not "
                  << "implemented yet. Please implement it yourself or select"
                  << "another one\n\n"
                  << "Availables ones\n"
                  << "- Euler (euler)\n"
                  << "- Runge-Kutta 4 (rk4)\n"
                  << "- Backward Euler - Predictor-Corrector (bepc)\n"
                  << "- Adams-Moulton 2 - Predictor-Corrector (am2pc)\n"
                  << "- Backward Euler - Fully Implicit (bdf1)\n"
                  << "- Adams-Moulton 2 - Fully Implicit (am2)\n"
                  << "- Backward Differentiation Formula 2 - Fully Implicit (bdf2)\n"
                  << "- Adaptive Runge-Kutta 4(5) Fehlberg (rk45f)\n"
                  << "- Adaptive Runge-Kutta 4(5) Dormand-Prince (rk45dp)\n"
                  << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
 }

}

