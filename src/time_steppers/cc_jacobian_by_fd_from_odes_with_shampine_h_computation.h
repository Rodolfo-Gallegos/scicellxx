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
#ifndef CCJACOBIANBYFDFROMODES_H
#define CCJACOBIANBYFDFROMODES_H

#include "ac_jacobian.h"
#include "../data_structures/ac_odes.h"

namespace scicellxx
{

 // A concrete class to compute the Jacobian matrix using Finite
 // Differences from a set of ODES
 class CCJacobianByFDFromODEs : virtual public ACJacobian
 {
  
 public:
  
  // Constructor
  CCJacobianByFDFromODEs();
  
  // Destructor
  ~CCJacobianByFDFromODEs();
   
  // In charge of computing the Jacobian using Finite Differences
  void compute_jacobian(const double time);
   
  // Set the ODEs to compute the Jacobian
  void set_ODEs(const ACODEs &odes);

 protected:

  // Stores the set of ODEs used to compute the Jacobian
  ACODEs ODEs;

  // A flag to indicate whether the ODEs have been set or not
  bool Set_of_ODEs_have_been_set;
   
 private:
  
  // Copy constructor (we do not want this class to be copiable because
  // it contains dynamically allocated variables, A in this
  // case). Check
  // http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  CCJacobianByFDFromODEs(const CCJacobianByFDFromODEs &copy)
   {
    BrokenCopy::broken_copy("CCJacobianByFDFromODEs");
   }
 
  // Copy constructor (we do not want this class to be copiable because
  // it contains dynamically allocated variables, A in this
  // case). Check
  // http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  void operator=(const CCJacobianByFDFromODEs &copy)
   {
    BrokenCopy::broken_assign("CCJacobianByFDFromODEs");
   }
   
 };
 
}

#endif // #ifndef CCJACOBIANBYFDFROMODES_H

