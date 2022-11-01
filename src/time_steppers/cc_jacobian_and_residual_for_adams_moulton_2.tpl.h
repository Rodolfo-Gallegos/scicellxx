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
#ifndef CCJACOBIANANDRESIDUALFORADAMSMOULTON2_TPL_H
#define CCJACOBIANANDRESIDUALFORADAMSMOULTON2_TPL_H

#include "ac_jacobian_and_residual_for_implicit_time_stepper.h"

// One of the possible Jacobian strategies to compute FY
#include "cc_jacobian_by_fd_and_residual_from_odes.h"

namespace scicellxx
{
 
 /// A concrete class to compute the Jacobian matrix and the residual
 /// vector for Adams-Moulton 2 or Trapezoidal Rule used for Newton's
 /// method
 template<class EQUATIONS_TYPE>
 class CCJacobianAndResidualForAdamsMoulton2 : virtual public ACJacobianAndResidualForImplicitTimeStepper<EQUATIONS_TYPE>
 {
  
 public:
  
  /// Constructor
  CCJacobianAndResidualForAdamsMoulton2();
  
  /// Empty destructor
  ~CCJacobianAndResidualForAdamsMoulton2();
  
  /// In charge of computing the Jacobian (virtual function
  /// implementation)
  void compute_jacobian();
  
  /// In charge of computing the residual
  void compute_residual();
   
  /// Enables the evaluation of the odes using the values of u at time
  /// 't'
  inline void enable_evaluation_of_odes_with_u_old_values() {Evaluate_odes_with_u_old_values=true;}
  
 private:
   
  /// Copy constructor (we do not want this class to be copiable because
  /// it contains dynamically allocated variables, A in this
  /// case). Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  CCJacobianAndResidualForAdamsMoulton2(const CCJacobianAndResidualForAdamsMoulton2 &copy)
   : ACJacobianAndResidualForImplicitTimeStepper<EQUATIONS_TYPE>()
   {
    BrokenCopy::broken_copy("CCJacobianAndResidualForAdamsMoulton2");
   }
   
  /// Copy constructor (we do not want this class to be copiable because
  /// it contains dynamically allocated variables, A in this
  /// case). Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  void operator=(const CCJacobianAndResidualForAdamsMoulton2 &copy)
   {
    BrokenCopy::broken_assign("CCJacobianAndResidualForAdamsMoulton2");
   }
   
  /// Temporary vector to store the evaluation of the odes with the
  /// values of u at time 't'. These are kept constant during Newton's
  /// iterations
  CCData *DUDT_old;
   
  /// A flag to indicate whether to evaluate the odes with the values
  /// of u at time 't'. These values are kept constant during Newton's
  /// iterations
  bool Evaluate_odes_with_u_old_values;
  
  /// A strategy to compute the Jacobian $\frac{\partial
  /// \mathbf{F}(\mathbf{Y})}{\partial \mathbf{Y}}$, where $\mathbf{Y}
  /// = \mathbf{y}_{k+1}$.
  CCJacobianByFDAndResidualFromODEs<EQUATIONS_TYPE> Jacobian_by_FD_strategy;
  
 };
 
}

#endif // #ifndef CCJACOBIANANDRESIDUALFORADAMSMOULTON2_TPL_H

