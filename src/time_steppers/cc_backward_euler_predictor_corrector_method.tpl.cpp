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
#include "cc_backward_euler_predictor_corrector_method.tpl.h"

namespace scicellxx
{

 // ===================================================================
 // Constructor
 // ===================================================================
 template<class EQUATIONS_TYPE>
 CCBackwardEulerPCMethod<EQUATIONS_TYPE>::CCBackwardEulerPCMethod()
  : ACPredictorCorrectorTimeStepper<EQUATIONS_TYPE>()
 {
  
  // Sets the number of history values
  this->N_history_values = 2;
  
  // Enable output messages
  //enable_output_messages();
  
 }
 
 // ===================================================================
 // Empty destructor
 // ===================================================================
 template<class EQUATIONS_TYPE>
 CCBackwardEulerPCMethod<EQUATIONS_TYPE>::~CCBackwardEulerPCMethod()
 {
 
 }
 
 // ===================================================================
 // Applies Backward Eulers method implemented as Predictor-Corrector
 // to the given odes from the current time "t" to the time "t+h".
 // The values of u at time t+h will be stored at index k (default k =
 // 0).
 // ===================================================================
 template<class EQUATIONS_TYPE>
 void CCBackwardEulerPCMethod<EQUATIONS_TYPE>::time_step(EQUATIONS_TYPE &odes,
                                                         const Real h,
                                                         const Real t,
                                                         CCData &u,
                                                         const unsigned k)
 {
#ifdef SCICELLXX_PANIC_MODE
  // Check if the ode has the correct number of history values to
  // apply Backwards-Eulers method
  const unsigned n_history_values = u.n_history_values();
  if (n_history_values < this->N_history_values)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The number of history values is less than\n"
                  << "the required by Backward Euler's method\n"
                  << "Required number of history values: "
                  << this->N_history_values << "\n"
                  << "Number of history values: "
                  << n_history_values << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
#endif // #ifdef SCICELLXX_PANIC_MODE
  
  // The method is implemented following an P(EC)^k E with the final
  // evaluation step as optionalw
  
  unsigned n_iterations = 0;
  
  // Get the number of odes
  const unsigned n_equations = odes.n_equations();
  
  // The residual vector
#ifdef SCICELLXX_USES_ARMADILLO
  CCVectorArmadillo local_error_vector(n_equations);
#else
  CCVector local_error_vector(n_equations);
#endif // #ifdef SCICELLXX_USES_ARMADILLO
  
  // Initialise local error with 0
  Real local_error = 0;
  
  // -----------------------------------------------------------------
  // -- Prediction phase --
  // -----------------------------------------------------------------
  // Temporary vector to store the evaluation of the odes.
  CCData dudt(n_equations);
  // Evaluate the ODE at time "t" using the current values of "u"
  // stored in index k
  odes.evaluate_time_derivatives(t, u, dudt, k);
  
  // Store the PREDICTED value by the external time stepper. Copy the
  // initial values from u
  CCData u_p(u);
  
  // Prediction step (Forward Euler)
  for (unsigned i = 0; i < n_equations; i++)
   {
    u_p(i,k) = u(i,k) + (h * (dudt(i)));
   }
  
  // -----------------------------------------------------------------
  // -- Evaluation phase
  // -----------------------------------------------------------------
  // -- Temporary vector to store the evaluation of the odes with the
  // -- predicted values.
  CCData dudt_p(n_equations);

  // Cache values
  const Real max_tol = this->maximum_tolerance();
  const Real min_tol = this->minimum_tolerance();
  const unsigned max_iter = this->maximum_iterations();
  bool out_msg = this->output_messages();
  bool fixed_iter = this->fixed_number_of_iterations();
  bool do_final_eval = this->perform_final_evaluation();
  
  
  do {
   // Evaluate the ODE at time "t+h" using the predicted values of
   // "u_p" stored at index k=0 because u_p has not history values
   odes.evaluate_time_derivatives(t+h, u_p, dudt_p, 0);
   
   // -----------------------------------------------------------------
   // -- Correction phase
   // -----------------------------------------------------------------
   for (unsigned i = 0; i < n_equations; i++)
   {
    u_p(i,k) = u(i,k) + (h * dudt(i));
   }
   
   // Compute error
   for (unsigned i = 0; i < n_equations; i++)
    {
     local_error_vector(i) = (u_p(i,k) - u(i,k)) / u_p(i,k);
    }
   
   // Get the maximum norm
   local_error = local_error_vector.norm_inf();
   // Is local error smaller than allowed tolerance
   if (local_error < min_tol)
    {
     if (out_msg)
      {
       scicellxx_output << "Local error is smaller than minimum tolerance value ["
                       << local_error << "] < [" << min_tol << "]" << std::endl;
      }
    }
   
   // Increase the number of iterations
   n_iterations++;

   if (n_iterations >= max_iter)
    {
     if (out_msg)
      {
       scicellxx_output << "Maximum number of iterations reached ["<< max_iter
                       <<"], local error [" << local_error << "], maximum_tolerance ["
                       << max_tol << "]\n"
                       << "You can change the maximum number of iterations by calling the method\n"
                       << "set_maximum_iterations()\n"
                       << "You can change the maximum tolerance by calling the method\n"
                       << "set_new_maximum_tolerance()" << std::endl;
      }
     
    }
   
   // Check whether a fixed number of iterations is enabled
   if (fixed_iter)
    {
     // Force local error to be greater than maximum tolerance
     local_error = max_tol + 1.0;
    }
   
   // Check whether reaching maximum number of iteratios or error in
   // tolerance ranges
  }while(local_error > max_tol && n_iterations < max_iter);
  
  // Perform a last evaluation such that the strategy becomes in a
  // E(PC)^k E
  if (do_final_eval)
   {
    // Evaluate the ODE at time "t" using the current values of "u"
    // stored in index k
    odes.evaluate_time_derivatives(t, u_p, dudt, k);
   }
  
  // Shift values to the right to provide storage for the new values
  u.shift_history_values();
  
  // Copy the values to the original vector
  for (unsigned i = 0; i < n_equations; i++)
   {
    u(i,k) = u_p(i,k);
   }
  
 }
 
}
