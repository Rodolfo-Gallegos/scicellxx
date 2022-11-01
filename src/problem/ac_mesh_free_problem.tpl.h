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
#ifndef ACMESHFREEPROBLEM_TPL_H
#define ACMESHFREEPROBLEM_TPL_H

#include "../general/common_includes.h"
#include "../general/utilities.h"

#include "../data_structures/cc_node.h"
#include "../time_steppers/ac_time_stepper.h"
#include "ac_problem.h"

// Matrices
#ifdef SCICELLXX_USES_ARMADILLO
#include "../matrices/cc_matrix_armadillo.h"
#else // If armadillo is not used then use the default solver
#include "../matrices/cc_matrix.h"
#endif // #ifdef SCICELLXX_USES_ARMADILLO

// Factory for linear solver
#include "../linear_solvers/cc_factory_linear_solver.h"

namespace scicellxx
{
 
 /// @class ACMeshFreeProblem ac_mesh_free_problem.h
 
 /// This class implements the interfaces for a mesh free problem for
 /// PDEs. It specifies a template for solving a problem thus one
 /// needs to create a class that inherents from this one to solve a
 /// particular Mesh Free Problem for a particular PDE
 template<class NODE_TYPE>
  class ACMeshFreeProblem : public virtual ACProblem
 {
  
 public:
  
  /// Constructor
  ACMeshFreeProblem();
  
  /// Destructor
  virtual ~ACMeshFreeProblem();
  
  /// Clean up, free allocated memory
  void clean_up_linear_solver();
  
  /// Sets default linear solver
  void set_default_linear_solver();
  
  /// Set the Linear solver
  void set_linear_solver(ACLinearSolver<MAT_TYPE, VEC_TYPE> *linear_solver_pt);
  
  // -------------------------------------------------------------------------
  // THESE METHODS MUST BE IMPLEMENTED IN THE CONCRETE PROBLEM CLASS [BEGIN]
  // -------------------------------------------------------------------------
  /// Set initial conditions
  virtual void set_initial_conditions() = 0;
  
  /// Set boundary conditions
  virtual void set_boundary_conditions() = 0;
  
  /// Document solution
  virtual void document_solution() = 0;
  // -------------------------------------------------------------------------
  // THESE METHODS MUST BE IMPLEMENTED IN THE CONCRETE PROBLEM CLASS [END]
  // -------------------------------------------------------------------------
  
  /// A helper function to complete the problem setup
  void complete_problem_setup() { }
  
  /// Problem steady solve (empty)
  void steady_solve();
  
  /// Problem unsteady solve
  void unsteady_solve();
  
  /// Set the time stepper for the time dependent problem
  inline void set_time_stepper(ACTimeStepper *time_stepper_pt) {Time_stepper_pt = time_stepper_pt;}
  
  /// Set radial base function
  inline void set_radial_base_function(ACRadialBaseFunction *rbf_pt) {RBF_pt = rbf_pt;}
  
  /* /// Gets access to the dimension of the problem */
  /* inline const unsigned dim() const {return Dim;} */
  
  /// Grants access to the total numer of nodes of the problem
  inline const unsigned n_nodes() const {return Nodes_pt.size();}
  
  /* /// Grants access to the number of nodes per dimension */
  /* inline const unsigned n_nodes_per_dim() const {return N_nodes_per_dim;} */
  
  /* /// Grants access to the number of nodes per boundary */
  /* inline const unsigned n_nodes_per_boundary() const {return N_nodes_per_boundary;} */
  
  /// Grants access to the number of boundary nodes
  inline const unsigned n_boundary_nodes() const {return Boundary_nodes_pt.size();}
  
  /// In charge of creating the nodes considering the geometry of the
  /// domain
  void create_nodes();
  
  /// Get access to i-th node
  const CCNode<Real> *node_pt(const unsigned i);
  
  /// Get access to i-th boundary node
  const CCBoundaryNode<Real> *boundary_node_pt(const unsigned i);
  
#ifdef SCICELLXX_PANIC_MODE
  /// Output nodes information (for debugging purposes)
  void output_nodes_information();
  
  /// Output boundary nodes information (for debugging purposes)
  void output_boundary_nodes_information();  
#endif // #ifdef SCICELLXX_PANIC_MODE
  
 protected:
  
  /// Copy constructor (we do not want this class to be
  /// copiable). Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
 ACMeshFreeProblem(const ACMeshFreeProblem &copy)
  : ACProblem()
   {
    BrokenCopy::broken_copy("ACMeshFreeProblem");
   }
  
  /// Assignment operator (we do not want this class to be
  /// copiable. Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  void operator=(const ACMeshFreeProblem &copy)
   {
    BrokenCopy::broken_assign("ACMeshFreeProblem");
   } 
  
  /// The set of actions to be performed before a time stepping
  virtual void actions_before_time_stepping() { }
  
  /// The set of actions to be performed after a time stepping
  virtual void actions_after_time_stepping() { } 
  
  /// The set of actions to be performed before newton solve
  virtual void actions_before_newton_solve() { }
  
  /// The set of actions to be performed after newton solve
  virtual void actions_after_newton_solve() { }
  
  /// The Time Stepper to approximate a solution to the ODEs
  ACTimeStepper *Time_stepper_pt;
  
  /// The PDE with the RBF
  ACRadialBaseFunction *RBF_pt;
  
  /// The vector with the nodes
  std::vector<CCNode<Real> *> Nodes_pt;
  
  /// BOUNDARY NODES SHOULD BE TREATED AS COMMON NODES AS WELL
  /// The vector with the boundary nodes
  std::vector<CCBoundaryNode<Real> *> Boundary_nodes_pt;
  
  /* /// The dimension of the problem */
  /* unsigned Dim; */
  
  /* /// The total number of nodes */
  /* unsigned N_nodes; */
  
  /* /// The number of nodes per dimension */
  /* unsigned N_nodes_per_dim; */
  
  /* /// The number of nodes per boundary */
  /* unsigned N_nodes_per_boundary; */
  
  /* /// The number of boundary nodes */
  /* unsigned N_boundary_nodes; */
  
  /// Interpolation matrix
  MAT_TYPE Interpolation_matrix;
  
  // Right-hand side vector
  VEC_TYPE RHS_vector;
  
  /// A pointer to the linear solver
  ACLinearSolver<MAT_TYPE, VEC_TYPE> *Linear_solver_pt;
  
  /// Flag to indicate whether to free the memory allocated to the
  /// linear solver
  bool Free_memory_for_linear_solver;
  
  /// Generate the interpolation matrix and the right-hand side vector
  virtual void generate_interpolation_matrix_and_rhs_vector(std::vector<CCNode<Real> *> &center_nodes) = 0;
  
  /// Generate interpolation matrix
  void generate_interpolation_matrix(std::vector<CCNode<Real> *> &center_nodes);
  
  /// Generate right-hand side vector
  void generate_rhs_vector();
  
 };
 
}

#endif // #ifndef ACMESHFREEPROBLEM_TPL_H

