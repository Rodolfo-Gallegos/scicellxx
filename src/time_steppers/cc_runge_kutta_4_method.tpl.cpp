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
#include "cc_runge_kutta_4_method.tpl.h"

namespace scicellxx
{

 // ===================================================================
 // Constructor
 // ===================================================================
 template<class EQUATIONS_TYPE>
 CCRK4Method<EQUATIONS_TYPE>::CCRK4Method()
  : ACTimeStepper<EQUATIONS_TYPE>()
 {
 
  // Sets the number of history values
  this->N_history_values = 2;
 
 }
 
 // ===================================================================
 // Empty destructor
 // ===================================================================
 template<class EQUATIONS_TYPE>
 CCRK4Method<EQUATIONS_TYPE>::~CCRK4Method()
 {
 
 }
 
 // ===================================================================
 // Applies Runge-Kutta 4 method to the given odes from the current
 // time "t" to the time "t+h". The values of u at time t+h will be
 // stored at index k (default k = 0).
 // ===================================================================
 template<class EQUATIONS_TYPE>
 void CCRK4Method<EQUATIONS_TYPE>::time_step(EQUATIONS_TYPE &odes,
                                             const Real h,
                                             const Real t,
                                             CCData &u,
                                             const unsigned k)
 {
#ifdef SCICELLXX_PANIC_MODE
  // Check if the ode has the correct number of history values to
  // apply Runge-Kutta 4 method
  const unsigned n_history_values = u.n_history_values();
  if (n_history_values < this->N_history_values)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The number of history values is less than\n"
                  << "the required by RK4 method" << std::endl
                  << "Required number of history values: "
                  << this->N_history_values << std::endl
                  << "Number of history values: "
                  << n_history_values << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
#endif // #ifdef SCICELLXX_PANIC_MODE
  
  // Get the number of odes
  const unsigned n_equations = odes.n_equations();
  
  // Temporary vector to store the evaluation of the odes.
  CCData dudt(n_equations);
  
  // Evaluate the ODE at time "t" using the current values of "u"
  odes.evaluate_time_derivatives(t, u, dudt, k);
  
  // Temporary vector to store the K_i evaluations proper of
  // Runge-Kutta methods
  CCData K1(n_equations);
  CCData K2(n_equations);
  CCData K3(n_equations);
  CCData K4(n_equations);
  
  // Create a copy of the u vector
  CCData u_copy(u);

  // --------------------------------------------------------------------
  // Runge-Kutta 4 method
  // --------------------------------------------------------------------
  // K1 = y(t, u)
  // K2 = y(t + \frac{1}{2}h, u + h(\frac{1}{2}K1))
  // K3 = y(t + \frac{1}{2}h, u + h(\frac{1}{2}K2))
  // K4 = y(t + h, u + hK3)
  // u(t+h) = u(t) + h(\frac{1}{6} K1 + \frac{1}{3} K2 + \frac{1}{3} K3 + \frac{1}{6} K4)
  // -------------------------------------------------------------------- 
  // Butcher tableau
  // --------------------------------------------------------------------
  // 0           |
  // \frac{1}{2} | \frac{1}{2}
  // \frac{1}{2} | 0           \frac{1}{2}
  // 1           | 0           0           1
  // --------------------------------------------------------------
  //             | \frac{1}{6} \frac{1}{3} \frac{1}{3} \frac{1}{6}
  
  // --------------------------------------------------------------------
  // Evaluate the ODE at time "t" using the current values of "u"
  odes.evaluate_time_derivatives(t, u, K1, k);
  // --------------------------------------------------------------------
  // Evaluate the ODE at time "t+(h/2)" and with "u+(h/2)K1"
  const Real h_half = h*0.5;
  for (unsigned i = 0; i < n_equations; i++)
   {
    u_copy(i,k) = u(i,k)+h_half*K1(i);
   }
  odes.evaluate_time_derivatives(t+h_half, u_copy, K2, k);
  
  // --------------------------------------------------------------------
  // Evaluate the ODE at time "t+(h/2)" and with "u+(h/2)K2"
  for (unsigned i = 0; i < n_equations; i++)
   {
    u_copy(i,k) = u(i,k)+h_half*K2(i);
   }
  odes.evaluate_time_derivatives(t+h_half, u_copy, K3, k);
  
  // -------------------------------------------------------------------- 
  // Evaluate the ODE at time "t+h" and with "u+hK3"
  for (unsigned i = 0; i < n_equations; i++)
   {
    u_copy(i,k) = u(i,k)+h*K3(i);
   }
  odes.evaluate_time_derivatives(t+h, u_copy, K4);
  
  // Shift values to the right to provide storage for the new values
  u.shift_history_values();
  
  // Once the derivatives have been obtained compute the new "u" as
  // the weighted sum of the K's
  for (unsigned i = 0; i < n_equations; i++)
   {
    u(i,k) = u(i,k+1) + (h/6.0)*(K1(i) + 2.0*K2(i) + 2.0*K3(i) + K4(i));
   }
 
 }

}
