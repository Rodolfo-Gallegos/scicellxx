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
/// IN THIS FILE: Implementation of the concrete class
/// CCNewtonMethodForBDF2

#include "cc_newtons_method_for_bdf_2.tpl.h"

namespace scicellxx
{

 // ===================================================================
 /// Constructor
 // ===================================================================
 template<class EQUATIONS_TYPE>
 CCNewtonsMethodForBDF2<EQUATIONS_TYPE>::CCNewtonsMethodForBDF2()
  : ACNewtonsMethodForImplicitTimeStepper<EQUATIONS_TYPE>()
 {
  // Set the Jacobian and residual strategy for Newton's method (used
  // for parent class)
  this->set_jacobian_and_residual_strategy(&Jacobian_and_residual_for_bdf_2);
 }
 
 // ===================================================================
 /// Empty destructor
 // ===================================================================
 template<class EQUATIONS_TYPE>
 CCNewtonsMethodForBDF2<EQUATIONS_TYPE>::~CCNewtonsMethodForBDF2()
 {
  
 }
 
 // ===================================================================
 /// Performs actions before initial convergence check
 // ===================================================================
 template<class EQUATIONS_TYPE>
 void CCNewtonsMethodForBDF2<EQUATIONS_TYPE>::actions_before_initial_convergence_check()
 {
  // Get the odes
  EQUATIONS_TYPE *odes_pt = this->odes_pt();
  // Get the time step
  const Real h = this->time_step();
  // Get the current time
  const Real t = this->current_time();
  // Get the u values
  CCData *u_pt = this->u_pt();
  // Get the index of the history values at time 't+h'
  const unsigned k = this->history_index();
  // Set the data for the computation of the jacobian and the residual
  Jacobian_and_residual_for_bdf_2.set_data_for_jacobian_and_residual(odes_pt, h, t, u_pt, k);
 }
 
}

