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
/// IN THIS FILE: Implementation of an abstract class to solve linear
/// systems of equations, this class is inhereted by any concrete
/// implementations of linear solvers

#include "ac_linear_solver.h"

namespace scicellxx
{

 // ===================================================================
 /// Constructor
 // ===================================================================
 ACLinearSolver::ACLinearSolver() 
  : Matrix_A_has_been_set(false)
 { }
 
 // ===================================================================
 /// Constructor where we specify the matrix A of size m X n
 // ===================================================================
 ACLinearSolver::ACLinearSolver(ACMatrix *const matrix_pt)
 {
  set_matrix_A(matrix_pt);
  
  // Set the flag to indicate that the matrix A has been set
  Matrix_A_has_been_set = true;
 }
 
 // ===================================================================
 /// Empty destructor
 // ===================================================================
 ACLinearSolver::~ACLinearSolver()
 {
  // Deallocate memory
  clean_up();
 }

 // ===================================================================
 /// Set the matrix A
 // ===================================================================
 void ACLinearSolver::set_matrix_A(ACMatrix *const matrix_pt)
 {
  // First clean any other previously stored matrix
  clean_up();
  
  // Set matrix A
  A_pt = matrix_pt;
  
  // Set the flag to indicate that the matrix A has been set
  Matrix_A_has_been_set = true;
 }

 // ===================================================================
 /// Clean up for any dynamically stored data
 // ===================================================================
 void ACLinearSolver::clean_up()
 {
  // Check whether the matrix has been set
  if (Matrix_A_has_been_set)
   {
    // Delete the content of the matrix
    //A_pt->clean_up();
    
    // We do not have to delete the values of the matrix A, we only
    // unset the pointer
    A_pt = 0;
    
    // Mark the matrix as not been set
    Matrix_A_has_been_set = false;
   }
 
 }

}
