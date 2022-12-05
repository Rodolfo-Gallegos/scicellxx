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
/// IN THIS FILE: Implementation of the abstract class
/// ACNewtonMethodForImplicitTimeStepper

#include "ac_newtons_method_for_implicit_time_stepper.tpl.h"

namespace scicellxx
{
 
 // ===================================================================
 /// Empty constructor
 // ===================================================================
 template<class EQUATIONS_TYPE>
 ACNewtonsMethodForImplicitTimeStepper<EQUATIONS_TYPE>::ACNewtonsMethodForImplicitTimeStepper()
  : CCNewtonsMethod(),
    ODEs_pt(NULL),
    U_pt(NULL),
    Data_for_jacobian_and_residual_has_been_set(false)
 {
  
 }
 
 // ===================================================================
 /// Empty destructor
 // ===================================================================
 template<class EQUATIONS_TYPE>
 ACNewtonsMethodForImplicitTimeStepper<EQUATIONS_TYPE>::~ACNewtonsMethodForImplicitTimeStepper()
 {
  
 }
 
 // ===================================================================
 /// Performs actions before Newton's method step
 // ===================================================================
 template<class EQUATIONS_TYPE>
 void ACNewtonsMethodForImplicitTimeStepper<EQUATIONS_TYPE>::actions_before_newton_step()
 {
  
 }
 
 // ===================================================================
 /// Performs actions after Newton's method step
 // ===================================================================
 template<class EQUATIONS_TYPE>
 void ACNewtonsMethodForImplicitTimeStepper<EQUATIONS_TYPE>::actions_after_newton_step()
 {
  // Update U, store the new values at index 'History_index'
  const unsigned k = History_index;
  const unsigned long n_data = U_pt->n_values();
  for (unsigned long i = 0; i < n_data; i++)
   {
    U_pt->value(i,k)=this->x_pt()->value(i);
   }
  
 }
 
 // ===================================================================
 /// Set data for Jacobian and residual computation. The odes, the time
 /// step 'h', the current time 't', the values of 'u' and the index
 /// where the values of 'u' at time 't+h' will be stored
 // ===================================================================
 template<class EQUATIONS_TYPE>
 void ACNewtonsMethodForImplicitTimeStepper<EQUATIONS_TYPE>::
 set_data_for_jacobian_and_residual(EQUATIONS_TYPE *odes_pt, const Real h, const Real t,
                                    CCData *u_pt, const unsigned k)
 {
  // Set the odes
  ODEs_pt = odes_pt;
  
  // Set the time step 
  Time_step = h;
  
  // Set the constant time
  Current_time = t;
  
  // Set the storage of the data
  U_pt = u_pt;
  
  // Set the index of where the values of u at time 't+h' should be
  // stored
  History_index = k;
  
  // Change the flag indicating that the data for the computation of
  // the Jacobian and the residual has been set
  Data_for_jacobian_and_residual_has_been_set = true;
  
 }
 
 // ===================================================================
 /// Set the strategy to compute the ODE's Jacobian
 // ===================================================================
 template<class EQUATIONS_TYPE>
 void ACNewtonsMethodForImplicitTimeStepper<EQUATIONS_TYPE>::set_strategy_for_odes_jacobian(ACJacobianAndResidualForImplicitTimeStepper<EQUATIONS_TYPE> *jacobian_strategy_for_odes_pt)
 {
  ACJacobianAndResidualForImplicitTimeStepper<EQUATIONS_TYPE> *cache_jacobian_strategy_pt = dynamic_cast<ACJacobianAndResidualForImplicitTimeStepper<EQUATIONS_TYPE> *>(this->jacobian_and_residual_strategy_pt());
  if (cache_jacobian_strategy_pt != NULL)
   {
    cache_jacobian_strategy_pt->set_strategy_for_odes_jacobian(jacobian_strategy_for_odes_pt);
   }
  else
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The dynamic cast was not succesful\n"
                  << "From [ACJacobianAndResidual<MAT_TYPE,VEC_TYPE> *]\n"
                  << "To [ACJacobianAndResidualForImplicitTimeStepper<MAT_TYPE, VEC_TYPE> *]\n" 
                  << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
 }
 
}

