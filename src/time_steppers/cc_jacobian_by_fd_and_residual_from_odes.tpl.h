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
#ifndef CCJACOBIANBYFDANDRESIDUALFROMODES_TPL_H
#define CCJACOBIANBYFDANDRESIDUALFROMODES_TPL_H

#include "../data_structures/data_structures.h"

#include "ac_jacobian_and_residual_for_implicit_time_stepper.h"

namespace scicellxx
{

 /// A concrete class to compute the Jacobian matrix using Finite
 /// Differences from a set of ODES
 template<class EQUATIONS_TYPE>
 class CCJacobianByFDAndResidualFromODEs : virtual public ACJacobianAndResidualForImplicitTimeStepper<EQUATIONS_TYPE>
 {
  
 public:
  
  /// Empty constructor
  CCJacobianByFDAndResidualFromODEs();
   
  /// Empty destructor
  ~CCJacobianByFDAndResidualFromODEs();
   
  /// In charge of computing the Jacobian using Finite Differences
  /// (virtual function implementation)
  void compute_jacobian();
   
  /// In charge of computing the residual
  void compute_residual();
   
 private:
   
  /// Copy constructor (we do not want this class to be copiable
  /// because it contains dynamically allocated variables, A in this
  /// case). Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  CCJacobianByFDAndResidualFromODEs(const CCJacobianByFDAndResidualFromODEs &copy)
   : ACJacobianAndResidualForImplicitTimeStepper<EQUATIONS_TYPE>()
   {
    BrokenCopy::broken_copy("CCJacobianByFDAndResidualFromODEs");
   }
   
  /// Copy constructor (we do not want this class to be copiable
  /// because it contains dynamically allocated variables, A in this
  /// case). Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  void operator=(const CCJacobianByFDAndResidualFromODEs &copy)
   {
    BrokenCopy::broken_assign("CCJacobianByFDAndResidualFromODEs");
   }
  
 };
 
}

#endif // #ifndef CCJACOBIANBYFDANDRESIDUALFROMODES_TPL_H

