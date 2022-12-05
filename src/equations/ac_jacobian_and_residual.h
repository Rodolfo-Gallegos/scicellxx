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

#ifndef ACJACOBIANANDRESIDUAL_H
#define ACJACOBIANANDRESIDUAL_H

#include "../general/general.h"
#include "../data_structures/data_structures.h"
#include "../matrices/matrices.h"

namespace scicellxx
{
 /// An abstract class used as a template for the algorithms to
 /// compute the Jacobian matrix
 class ACJacobianAndResidual
 {
  
 public:
  
  /// Constructor
  ACJacobianAndResidual();
  
  /// Destructor
  virtual ~ACJacobianAndResidual();
  
  /// In charge of computing the Jacobian based on the particular
  /// strategy implemented in the derived class
  virtual void compute_jacobian() = 0;
  
  /// Get access to the Jacobian
  inline ACMatrix *jacobian_pt() {return Jacobian_pt;}
  
  /// In charge of computing the residual vector based on the
  /// particular strategy implemented in the derived class
  virtual void compute_residual() = 0;
   
  /// Get access to the residual
  inline ACVector *residual_pt() {return Residual_pt;}
   
 protected:
  
  /// Storage for the Jacobian matrix
  ACMatrix *Jacobian_pt;
  
  /// Storage for the residual vector
  ACVector *Residual_pt;
  
 private:
   
  /// Copy constructor (we do not want this class to be
  /// copiable. Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  ACJacobianAndResidual(const ACJacobianAndResidual &copy)
   {
    BrokenCopy::broken_copy("ACJacobianAndResidual");
   }
 
  /// Copy constructor (we do not want this class to be
  /// copiable. Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  void operator=(const ACJacobianAndResidual &copy)
   {
    BrokenCopy::broken_assign("ACJacobianAndResidual");
   }
  
 };
 
}

#endif // #ifndef ACJACOBIANANDRESIDUAL_H

