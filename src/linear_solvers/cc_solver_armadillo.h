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
/// IN THIS FILE: The definition of the concrete class
/// CCSolverArmadillo to solve systems of equations. This class calls
/// the methods solve() or spsolve() from Armadillo to perform the
/// solution of the system of equations.

/// Check whether the class has been already defined
#ifndef CCSOLVERARMADILLO_H
#define CCSOLVERARMADILLO_H

// Include the header from inherited class
#include "ac_linear_solver.h"

namespace scicellxx
{
 
 /// A concrete class for solving a linear system of equations. This
 /// class uses the methods solve() or spsolve() from Armadillo to
 /// perform the solution of the system of equations.
 class CCSolverArmadillo : public virtual ACLinearSolver
 {
  
 public:
  
  /// Empty constructor
  CCSolverArmadillo();
  
  /// Constructor where we specify the matrix A
  CCSolverArmadillo(ACMatrix *const A_mat_pt);
  
  /// Empty destructor
  ~CCSolverArmadillo();
  
  /// Solves a system of equations with input A_mat. We specify the
  /// right-hand side B and the X matrices where the results are
  /// returned. We assume that the input/output matrices have the
  /// correct dimensions: A_mat.n_columns() x A_mat.n_rows() for B, and
  /// A_mat.n_rows() x A_mat.n_columns() for X.
  void solve(ACMatrix *const A_mat_pt, const ACMatrix *const B_pt, ACMatrix *const X_pt);
   
  /// Solves a system of equations with input A_mat. We specify the
  /// right-hand side b and the x vector where the result is
  /// returned. We assume that the input/output vectors have the
  /// correct dimensions: A_mat.n_columns() for b, and A_mat.n_rows()
  /// for x.
  void solve(ACMatrix *const A_mat_pt, const ACVector *const b_pt, ACVector *const x_pt);
   
  /// Solve a system of equations with the already stored matrix A. We
  /// specify the right-hand side B and the X matrices where the
  /// results are returned, B and X may be 1-column matrices
  /// (vectors). We assume that the input/output matrices have the
  /// correct dimensions.
  void solve(const ACMatrix *const B_pt, ACMatrix *const X_pt);
   
  /// Solve a system of equations with the already stored matrix A. We
  /// specify the right-hand side b and the x vectors where the result
  /// is returned. We assume that the input/output vectors have the
  /// correct dimensions: A.ncolumns() for b, and A.nrows() for x.
  void solve(const ACVector *const b_pt, ACVector *const x_pt);
   
 private:
 
  /// Copy constructor (we do not want this class to be copiable
  /// because it contains dynamically allocated variables, A in this
  /// case). Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
 CCSolverArmadillo(const CCSolverArmadillo &copy)
  : ACLinearSolver()
   {
    BrokenCopy::broken_copy("CCSolverArmadillo");
   }
  
  /// Copy constructor (we do not want this class to be copiable because
  /// it contains dynamically allocated variables, A in this
  /// case). Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  void operator=(const CCSolverArmadillo &copy)
   {
    BrokenCopy::broken_assign("CCSolverArmadillo");
   }
   
 };
 
}

#endif // #ifndef CCSOLVERARMADILLO_H
