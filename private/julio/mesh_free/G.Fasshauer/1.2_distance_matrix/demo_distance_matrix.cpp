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
// This demo driver is based on the Program 1.2 DistanceMatrixFit.m
// from the book "Meshfree Approximation Methods with MATLAB, Gregory
// E. Fasshauer", World Scientific Publishing, 2007

// Include SciCell++ libraries
#include "../../../../../src/scicellxx.h"

using namespace scicellxx;

template<class MAT_TYPE, class VEC_TYPE>
class Problem
{
public:
 Problem();
 ~Problem();
 
 void initialise_problem();
 void configure_problem();
 void solve_problem();
 void document_solution();
};

// This fucntion has it maximum value at the center, depending on the
// dimension s. At the boundaries it is zero.
const Real test_function(ACVector *x_pt, const unsigned s)
{
 Real prod=1.0;
 for (unsigned i = 0; i < s; i++)
  {
   Real x_i = x_pt->get_value(i);
   prod*=x_i*(1.0-x_i);
  }
 return pow(4, s)*prod;
}

void compute_distance_matrix(ACMatrix *data_sites_pt, ACMatrix *centers_pt,
                             ACMatrix *distance_matrix_pt)
{
 // Get the number of "vector points" on "data_sites_pt"
 // Get the number of "vector points" on "centers_pt"
 const unsigned n_vector_points_data_sites = data_sites_pt->n_columns();
 const unsigned n_vector_points_centers = centers_pt->n_columns();

 // The dimension of input vector points must be the same, otherwise
 // there is an error
 const unsigned dimension = data_sites_pt->n_rows();
 const unsigned tmp_dimension = centers_pt->n_rows();

 if (dimension != tmp_dimension)
  {
   // Error message
   std::ostringstream error_message;
   error_message << "The dimensions of the data sites vector and the\n"
                 << "centers vector are different\n"
                 << "dim(data_site):" << dimension
                 << "\ndim(centers):" << tmp_dimension
                 << std::endl;
   throw SciCellxxLibError(error_message.str(),
                          SCICELLXX_CURRENT_FUNCTION,
                          SCICELLXX_EXCEPTION_LOCATION);
  }
 
 // A factory to create matrices and vectors
 CCFactoryMatrices factory_matrices_and_vectors;
 ACVector *distance_pt = factory_matrices_and_vectors->create_vector(dimension);
 
 // Loop over all the data points in the first matrix
 for (unsigned m = 0; m < n_vector_points_data_sites; m++)
  {
   // Loop over all the data points in the second matrix
   for (unsigned n = 0; n < n_vector_points_centers; n++)
    {
     // Loop over the elements of both vectors
     for (unsigned k = 0; k < dimension; k++)
      {
       Real dis = data_sites_pt->get_value(k, m) - centers_pt->get_value(k, n);
       distance_pt->set_value(k, dis);
      }
     Real norm2 = distance_pt->norm_2();
     distance_matrix_pt->set_value(m,n, norm2);
    }
  }

 delete distance_pt;
 
}

struct Args {
 argparse::ArgValue<unsigned> dimension;
 argparse::ArgValue<unsigned> degree;
 argparse::ArgValue<unsigned> n_evaluation_points_per_dimension;
};

// ==================================================================
// ==================================================================
// ==================================================================
// Main function
// ==================================================================
// ==================================================================
// ==================================================================
int main(int argc, char *argv[])
{
 // Initialise scicellxx
 initialise_scicellxx();

 // A factory to create matrices and vectors
 CCFactoryMatrices factory_matrices_and_vectors;
 
 // A factory to create the linear solver
 CCFactoryLinearSolver factory_linear_solver;
 
 // Instantiate parser
 Args args;
 auto parser = argparse::ArgumentParser(argv[0], "Distance Matrix");
 
 // Add arguments 
 parser.add_argument<unsigned>(args.dimension, "--dim")
  .help("Dimension")
  .default_value("2")
  .choices({"1", "2", "3"});
 
 parser.add_argument<unsigned>(args.degree, "--degree")
  .help("Degree")
  .default_value("3")
  .choices({"1", "2", "3"});
 
 parser.add_argument<unsigned>(args.n_evaluation_points_per_dimension, "--n_evaluation_points_per_dimension")
  .help("Number of evaluation points per dimension")
  .default_value("1000");
 
 // Parse the input arguments
 parser.parse_args(argc, argv);
 
 // --------------------------------------------------------------
 // Domain specification
 // --------------------------------------------------------------
 // TODO: Create a DOMAIN (mesh?) type class
 
 // Dimension of the problem
 unsigned dim = args.dimension; // (if you change the dimension then
                                // also change
                                // 'n_evaluation_points_per_dimension')
 
 // Interpolant degree
 unsigned degree = args.degree;
 
 // Specify the one-dimensional lenght of the domain
 const unsigned L = 1.0;
 
 // --------------------------------------------------------------
 // Create and give position to nodes
 // --------------------------------------------------------------
 // Nodes per dimension
 const unsigned n_nodes_per_dim = pow(2, degree+1);
 // The number of nodes
 const unsigned n_nodes = pow(n_nodes_per_dim, dim);
 // A vector of nodes
 std::vector<CCNode *> nodes_pt(n_nodes);
 
 // Number of variables stored in the node
 const unsigned n_variables = 1;
 // Number of history values
 const unsigned n_history_values = 1;
 
 // Distance between a pair of nodes
 const Real h = L / (Real)(n_nodes_per_dim - 1);
 std::vector<Real> x(dim, 0.0);
 // Create the nodes
 for (unsigned i = 0; i < n_nodes; i++)
  {
   nodes_pt[i] = new CCNode(dim, n_variables, n_history_values);
  }

 // --------------------------------------------------------------
 // Output supporting nodes
 // --------------------------------------------------------------
 std::ofstream nodes_file("RESLT/nodes.dat");
 
 // Assign positions
 for (unsigned i = 0; i < n_nodes; i++)
  {
   for (unsigned k = 0; k < dim; k++)
    {
     const Real r = rand();
     const Real position = static_cast<Real>(r / RAND_MAX) * L;
     // Generate position and assign it
     //const Real position = x[k];
     nodes_pt[i]->set_position(position, k); 
     //x[k]+=h;
     nodes_file << position;
     if (k + 1 < dim)
      {
       nodes_file << " ";
      }
    }
   nodes_file << std::endl;
  }
 
 // Close support nodes file
 nodes_file.close();
 
 // --------------------------------------------------------------
 // Set initial conditions
 // --------------------------------------------------------------
 for (unsigned i = 0; i < n_nodes; i++)
  {
   for (unsigned j = 0; j < n_variables; j++)
    {
     const Real u = 0.0;
     nodes_pt[i]->set_variable(u, j);
    }
  }
 
 // --------------------------------------------------------------
 // Set boundary conditions
 // --------------------------------------------------------------
 
 // Move the first and the last node to the boundary of the domain
 //nodes_pt[0]->set_position(0.0, 0);
 //nodes_pt[0]->set_variable(0.0, 0);
 //nodes_pt[n_nodes-1]->set_position(1.0, 0);
 //nodes_pt[n_nodes-1]->set_variable(1.0, 0);
 
 // --------------------------------------------------------------
 // Set the problem and solve it
 // --------------------------------------------------------------
 
 // TODO: The distance matrix may be formed while we loop over the
 // nodes to extract their position
 
 // --------------------------------------------------------------
 // Loop over the nodes and extract their position and store them in a
 // matrix
 // --------------------------------------------------------------
 ACMatrix *nodes_position_pt = factory_matrices_and_vectors.create_matrix(dim, n_nodes);
 // Each column stores the vector position of a node
  for (unsigned i = 0; i < n_nodes; i++)
  {
   for (unsigned j = 0; j < dim; j++)
    {
     Real pos = nodes_pt[i]->get_position(j);
     nodes_position_pt->set_value(j, i, pos);
    }
  }
 
 // -------------------------------------------------------------- 
 // Create the distance matrix
 // --------------------------------------------------------------
  ACMatrix *distance_matrix_pt = factory_matrices_and_vectors.create_matrix(n_nodes, n_nodes);
 // --------------------------------------------------------------
 // Generate the distance matrix using the nodes position centers
 // shifted by the same nodes position
 // --------------------------------------------------------------
  
  compute_distance_matrix(nodes_position_pt, nodes_position_pt, distance_matrix_pt);
  
  // --------------------------------------------------------------
  // Set right-hand side
  // --------------------------------------------------------------
  ACVector *rhs_pt = factory_matrices_and_vectors.create_vector(n_nodes);
  ACVector *tmp_v_pt = factory_matrices_and_vectors.create_vector(dim);
 for (unsigned i = 0; i < n_nodes; i++)
  {
   for (unsigned j = 0; j < dim; j++)
    {
     Real pos = nodes_pt[i]->get_position(j);
     tmp_v_pt->set_value(j, pos);
    }
   // --------------------------------------------------------------
   // Evaluate the KNOWN function at the centers positions
   // --------------------------------------------------------------
   Real test_function_value = test_function(tmp_v_pt, dim);
   rhs_pt->set_value(i, test_function_value);
  }
 
 // The solution vector (with the respective number of rows) stores
 // the coefficients for the interpolant polynomials
 ACVector *sol_pt = factory_matrices_and_vectors.create_vector(n_nodes);
 
 // --------------------------------------------------------------
 // Solve
 // -------------------------------------------------------------- 
 // Create the linear solver
 ACLinearSolver *linear_solver_pt = factory_linear_solver.create_linear_solver();
 
 std::cerr << "Distance matrix" << std::endl;
 //distance_matrix.print();
 
 // --------------------------------------------------------------
 // Solve the system of equations
 // --------------------------------------------------------------
 linear_solver_pt->solve(distance_matrix_pt, rhs_pt, sol_pt);
 std::cerr << "Solution vector" << std::endl;
 //sol.print();

 std::cerr << "Nodes positions and values" << std::endl;
 /*
 // Show results
 for (unsigned i = 0; i < n_nodes; i++)
  {
   nodes_pt[i]->print(true);
  }
 */
 
 // --------------------------------------------------------------
 // --------------------------------------------------------------
 // EVALUATION STAGE
 // --------------------------------------------------------------
 // --------------------------------------------------------------
 std::cerr << "\nEVALUATION\n" << std::endl;
 
 // --------------------------------------------------------------
 // Evaluate (compute error RMSE)
 // --------------------------------------------------------------
 unsigned n_evaluation_points_per_dimension = args.n_evaluation_points_per_dimension;
 const unsigned n_data_in_evaluation_points = pow(n_evaluation_points_per_dimension, dim);
 // Distance between a pair of nodes
 const Real h_test = L / (Real)(n_evaluation_points_per_dimension - 1);
 
 // Compute approximated solution at new positions
 ACMatrix *approx_solution_position_pt = factory_matrices_and_vectors->create_matrix(dim, n_evaluation_points_per_dimension);
 //approx_solution_position.allocate_memory();
 // --------------------------------------------------------------
 // Assign positions
 // --------------------------------------------------------------
 std::vector<Real> x_eval(dim, 0.0);
 for (unsigned i = 0; i < n_evaluation_points_per_dimension; i++)
  {
   for (unsigned k = 0; k < dim; k++)
    {
     const Real r = rand();
     const Real position = static_cast<Real>(r / RAND_MAX) * L;
     // Generate position and assign it
     //const Real position = x_eval[k];
     approx_solution_position_pt->set_value(k, i, position);
     //x_eval[k]+=h_test;
    }
  }
 
 // Compute distance matrix with new positions
 ACMatrix *approx_distance_matrix_pt = factory_matrices_and_vectors->create_matrix(n_evaluation_points_per_dimension, n_nodes);
 // --------------------------------------------------------------
 // Generate the distance matrix using the nodes position centers
 // shifted by the new positions
 // --------------------------------------------------------------
 compute_distance_matrix(approx_solution_position_pt, nodes_position_pt, approx_distance_matrix_pt);
 //approx_distance_matrix.print();

 // Approximated solution
 ACVector *approx_sol_pt = factory_matrices_and_vectors->create_vector(n_evaluation_points_per_dimension);

 // HERE HERE HERE
 // Approximate solution at given points
 multiply_matrix_times_vector(approx_distance_matrix_pt, sol_pt, approx_sol_pt);
 (*approx_sol_pt) = (*approx_distance_matrix_pt) * (*sol_pt);
 approx_distance_matrix_pt->multiply by vector
 
 // --------------------------------------------------------------
 // Output data for plotting
 // --------------------------------------------------------------
 std::ofstream output_file("RESLT/output.dat");
 for (unsigned i = 0; i < n_evaluation_points_per_dimension; i++)
  {
   for (unsigned k = 0; k < dim; k++)
    {
     output_file << approx_solution_position(k, i) << " ";
    }
   output_file << approx_sol(i) << std::endl;
  }
 
 // Close output file
 output_file.close();
 
 // --------------------------------------------------------------
 // Get real solution at given points and get the error 
 // --------------------------------------------------------------
#ifdef SCICELLXX_USES_ARMADILLO
 CCVectorArmadillo<Real> real_sol(n_evaluation_points_per_dimension);
#else 
 CCVector<Real> real_sol(n_evaluation_points_per_dimension);
#endif // #ifdef SCICELLXX_USES_ARMADILLO
 //real_sol.allocate_memory();
 for (unsigned i = 0; i < n_evaluation_points_per_dimension; i++)
  {
#ifdef SCICELLXX_USES_ARMADILLO
   CCVectorArmadillo<Real> tmp_v(dim);
#else
   CCVector<Real> tmp_v(dim);
#endif // #ifdef SCICELLXX_USES_ARMADILLO
   //tmp_v.allocate_memory();
   for (unsigned j = 0; j < dim; j++)
    {
     tmp_v(j) = approx_solution_position(j, i);
    }
   // ------------------------
   // Evaluation at approx_solution_position
   real_sol(i) = test_function(tmp_v, dim);
  }
 
 // --------------------------------------------------------------
 // Compute error
 // --------------------------------------------------------------
#ifdef SCICELLXX_USES_ARMADILLO
 CCVectorArmadillo<Real> error(n_evaluation_points_per_dimension);
#else
 CCVector<Real> error(n_evaluation_points_per_dimension);
#endif // #ifdef SCICELLXX_USES_ARMADILLO
 //error.allocate_memory();
 std::cerr << "ERRORS" << std::endl;
 for (unsigned i = 0; i < n_evaluation_points_per_dimension; i++)
  {
   error(i) = real_sol(i) - approx_sol(i);
   //std::cerr << i << ": " << std::fabs(error(i)) << std::endl;
   //std::cerr << i << ": " << real_sol(i) << ":" << approx_sol(i) << std::endl;
  }
 
 const Real rms_error = error.norm_2() / sqrt(n_data_in_evaluation_points);
 
 // --------------------------------------------------------------
 // Output error
 // --------------------------------------------------------------
 std::ofstream error_file("RESLT/error.dat");
 for (unsigned i = 0; i < n_evaluation_points_per_dimension; i++)
  {
   for (unsigned k = 0; k < dim; k++)
    {
     error_file << approx_solution_position(k, i) << " ";
    }
   error_file << error(i) << std::endl;
  }
 
 // Close error file
 error_file.close();
 
 // --------------------------------------------------------------
 // Summary
 // --------------------------------------------------------------
 std::cerr << std::endl;
 std::cerr << "Polynomial degree: " << degree << std::endl;
 std::cerr << "N. nodes per dimension: " << n_nodes_per_dim << std::endl;
 std::cerr << "N. total nodes: " << n_nodes << std::endl; 
 std::cerr << "RMS-error: " << rms_error << std::endl;
 
 // ==============================================================
 // ==============================================================

 // Free memory
 delete linear_solver_pt;
 
 delete sol_pt;
 delete tmp_v_pt;
 delete rhs_pt;
 delete distance_matrix_pt;
 delete nodes_position_pt;
 
 // --------------------------------------------------------------
 // Delete nodes storage
 // --------------------------------------------------------------
 for (unsigned i = 0; i < n_nodes; i++)
  {
   delete nodes_pt[i];
  }
 
 // Finalise chapcom
 finalise_scicellxx();
 
 return 0;
 
}

