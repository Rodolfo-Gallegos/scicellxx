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
#ifndef SCICELLXX_TIMESTEPPERS_H
#define SCICELLXX_TIMESTEPPERS_H

// Include the integration methods (time steppers)
#include "ac_time_stepper.h"

#include "ac_adaptive_new_step_size_strategy.h"
#include "ac_adaptive_time_stepper.h"
#include "ac_jacobian_and_residual_for_implicit_time_stepper.h"
#include "ac_newtons_method_for_implicit_time_stepper.h"
#include "ac_predictor_corrector_time_stepper.h"

#include "cc_adams_moulton_2_method.h"
#include "cc_adams_moulton_2_predictor_corrector_method.h"
#include "cc_adaptive_new_step_size_half_double.h"
#include "cc_adaptive_runge_kutta_45DP_method.h"
#include "cc_adaptive_runge_kutta_45F_method.h"
#include "cc_backward_euler_method.h"
#include "cc_backward_euler_predictor_corrector_method.h"
#include "cc_bdf_2_method.h"
#include "cc_euler_method.h"
#include "cc_jacobian_and_residual_for_adams_moulton_2.h"
#include "cc_jacobian_and_residual_for_backward_euler.h"
#include "cc_jacobian_and_residual_for_bdf_2.h"
#include "cc_jacobian_by_fd_and_residual_from_odes.h"
//#include "cc_jacobian_by_fd_from_odes_with_shampine_h_computation.h"
#include "cc_newtons_method_for_adams_moulton_2.h"
#include "cc_newtons_method_for_backward_euler.h"
#include "cc_newtons_method_for_bdf_2.h"
#include "cc_runge_kutta_4_method.h"

// Include the factory class
#include "cc_factory_time_stepper.h"
 
#endif // #ifndef SCICELLXX_TIMESTEPPERS_H

