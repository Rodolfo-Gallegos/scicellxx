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
#ifndef CCNODEIMQPOISSON_H
#define CCNODEIMQPOISSON_H

#include "../../general/general.h"

#include "../rbf/cc_inverse_multiquadric.h"

namespace scicellxx
{
 /// A class than implements an Inverse Muliquadric RBF Node for the
 /// solution of the Poisson equation
 template<class NODE_TYPE>
  class CCNodeIMQPoisson : public virtual NODE_TYPE
 {


  ////////////////// HERE // We do not always need boundary nodes, we
  //////////////////only need boundary nodes when they are at the
  //////////////////boundary, we expect that we can create a class of
  //////////////////CCNodeIMQPoisson with the boundary features of the
  //////////////////nodes, and other class CCNodeIMQPoisson without
  //////////////////the boundary features of the node

  
 public:
  
  /// Constructor
  CCBoundaryNode(const unsigned boundary, const unsigned dimension, const unsigned n_variables, const unsigned n_history_values=1);
  
  /// Empty destructor
  virtual ~CCBoundaryNode();

  ///////////////////////*-*-*-*-*-*-*-*-*-*
  
 public:
  
  // Constructor
  CCNodeIMQPoisson();
  
  // Destructor
  virtual ~CCNodeIMQPoisson();
  
  // Laplacian psi
  const Real Lpsi(const Real r);
  
 private:
  
  // Copy constructor (we do not want this class to be
  // copiable. Check
  // http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  CCNodeIMQPoisson(const CCNodeIMQPoisson &copy)
   {
    BrokenCopy::broken_copy("CCNodeIMQPoisson");
   }
  
  // Copy constructor (we do not want this class to be
  // copiable. Check
  // http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  void operator=(const CCNodeIMQPoisson &copy)
   {
    BrokenCopy::broken_assign("CCNodeIMQPoisson");
   }
  
  /// The RBF
  CCInverseMultiquadric 
  
  
 };
 
}

#endif // #ifndef CCNODEIMQPOISSON_H
