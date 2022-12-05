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
/// IN THIS FILE: The definition of an abstract class to solve linear
/// systems of equations. Concrete or specific solvers MUST inherent
/// from this abstract class and implement the methods solve() and
/// resolve()

/// IN THIS FILE: The definition of the base class for linear solvers
#ifndef ACLINEARSOLVER_H
#define ACLINEARSOLVER_H

#include "../general/general.h"

#include "../matrices/matrices.h"

namespace scicellxx
{
 
 /// Abstract class to solve linear systems of equations, this class is
 /// inhereted by any concrete implementations of linear solvers.
 class ACLinearSolver
 {
  
 public:
  
  /// Constructor
  ACLinearSolver();
  
  /// Constructor where we specify the matrix A
  ACLinearSolver(ACMatrix * const matrix_pt);
  
  /// Empty destructor
  virtual ~ACLinearSolver();
  
  /// Set the matrix A
  void set_matrix_A(ACMatrix *const matrix_pt);
  
  /// Clean up for any dynamically stored data
  void clean_up();
  
  /// Virtual function to solve a system of equations with input
  /// A_mat. We specify the right-hand side B and the X matrices where
  /// the results are returned. We assume that the input/output
  /// matrices have the correct dimensions: A_mat.ncolumns() x
  /// A_mat.nrows() for B, and A_mat.nrows() x A_mat.ncolumns() for X.
  virtual void solve(ACMatrix *const A_mat_pt, const ACMatrix *const B_pt, ACMatrix *const X_pt) = 0;
  
  /// Virtual function to solve a system of equations with input
  /// A_mat. We specify the right-hand side b and the x vector where
  /// the result is returned. We assume that the input/output vectors
  /// have the correct dimensions: A_mat.ncolumns() for b, and
  /// A_mat.nrows() for x.
  virtual void solve(ACMatrix *const A_mat_pt, const ACVector *const b_pt, ACVector *const x_pt) = 0;
  
  /// Virtual function to solve a system of equations with the already
  /// stored matrix A. We specify the right-hand side B and the X
  /// matrices where the results are returned. We assume that the
  /// input/output matrices have the correct dimensions: A.ncolumns() x
  /// A.nrows() for B, and A.nrows() x A.ncolumns() for X.
  virtual void solve(const ACMatrix *const B_pt, ACMatrix *const X_pt) = 0;
  
  /// Virtual function to solve a system of equations with the already
  /// stored matrix A. We specify the right-hand side b and the x
  /// vectors where the result is returned. We assume that the
  /// input/output vectors have the correct dimensions: A.ncolumns()
  /// for b, and A.nrows() for x.
  virtual void solve(const ACVector *const b_pt, ACVector *const x_pt) = 0;
  
  /// Virtual function to re-solve a system of equations with the
  /// already stored matrix A (re-use of the LU decomposition or call
  /// the solve method for an iterative solver). BROKEN beacuse
  /// iterative solvers may not implement it. We specify the right-hand
  /// side B and the X matrices where the results are returned. We
  /// assume that the input/output vectors have the correct dimensions:
  /// A.ncolumns() x A.nrows() for B, and A.nrows() x A.ncolumns() for
  /// X.
  virtual void resolve(const ACMatrix *const B_pt, ACMatrix *const X_pt)
  {
   /// Error message
   std::ostringstream error_message;
   error_message << "Virtual function to resolve systems of equations should be\n"
                 << "implemented in derived class" << std::endl;
   throw SciCellxxLibError(error_message.str(),
                          SCICELLXX_CURRENT_FUNCTION,
                          SCICELLXX_EXCEPTION_LOCATION);
  }
  
  /// Virtual function to re-solve a system of equations with the
  /// already stored matrix A (re-use of the LU decomposition or call
  /// the solve method for an iterative solver). BROKEN beacuse
  /// iterative solvers may not implement it. We specify the right-hand
  /// side b and the x vector where the result is returned. We assume
  /// that the input/output vectors have the correct dimensions:
  /// A.ncolumns() for b, and A.nrows() for x.
  virtual void resolve(const ACVector *const b_pt, ACVector *const x_pt)
  {
   // Error message
   std::ostringstream error_message;
   error_message << "Virtual function to resolve systems of equations should be\n"
                 << "implemented in derived class" << std::endl;
   throw SciCellxxLibError(error_message.str(),
                          SCICELLXX_CURRENT_FUNCTION,
                          SCICELLXX_EXCEPTION_LOCATION);
  }
  
 protected:
  
  /// The matrix A
  ACMatrix *A_pt;
  
  /// Flag to indicate whether the matrix A has been set
  bool Matrix_A_has_been_set;
  
 private:
 
  /// Copy constructor (we do not want this class to be copiable). Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  ACLinearSolver(const ACLinearSolver &copy)
   {
    BrokenCopy::broken_copy("ACLinearSolver");
   }
 
  /// Assignment operator (we do not want this class to be
  /// copiable. Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  void operator=(const ACLinearSolver &copy)
   {
    BrokenCopy::broken_assign("ACLinearSolver");
   }
 
 };

}
 
#endif // #ifndef ACLINEARSOLVER_H
