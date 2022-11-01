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
/// IN THIS FILE: The implementation of the factory for the creation
/// of linear solvers

#include "cc_factory_linear_solver.h"

namespace scicellxx
{
 
 // ===================================================================
 /// Empty constructor
 // ===================================================================
 CCFactoryLinearSolver::CCFactoryLinearSolver()
 {
  
 }
 
 // ===================================================================
 /// Empty destructor
 // ===================================================================
 CCFactoryLinearSolver::~CCFactoryLinearSolver()
 { 
  
 }

 // ===================================================================
 /// Returns the corresponding linear solver (based on compilation
 /// options)
 //===================================================================
 ACLinearSolver* CCFactoryLinearSolver::create_linear_solver()
 {
  // ------------------------------------------------------
  // Check what linear solver we need to create
  // ------------------------------------------------------
#ifdef SCICELLXX_USES_ARMADILLO
  return new CCSolverArmadillo();
#else
  return new CCLUSolverNumericalRecipes();
#endif // #ifdef SCICELLXX_USES_ARMADILLO
 }
 
 // ===================================================================
 /// Returns the specified linear solver
 // ===================================================================
 ACLinearSolver* CCFactoryLinearSolver::create_linear_solver(std::string linear_solver_name)
 {
  // Get the string and change it to lower case 
  std::transform(linear_solver_name.begin(), linear_solver_name.end(),
                 linear_solver_name.begin(), ::tolower);
  
  // ------------------------------------------------------
  // Check what linear solver we need to create
  // ------------------------------------------------------
  // LU solver from numerical recipes
  if (linear_solver_name.compare("numerical_recipes")==0)
   {
    return new CCLUSolverNumericalRecipes();
   }
#ifdef SCICELLXX_USES_ARMADILLO
  // Linear solver from Armadillo
  else if (linear_solver_name.compare("armadillo")==0)
   {
    return new CCSolverArmadillo();
   }
#endif // #ifdef SCICELLXX_USES_ARMADILLO
  else
   {
    std::ostringstream error_message;
    error_message << "The linear solver you want to use is not implemented yet.\n"
                  << "Please implement it yourself or select another one\n\n"
                  << "Availables ones\n"
                  << "- LU linear solver from Numerical Recipes (numerical_recipes)\n"
                  << "- Armadillo Linear Solver (armadillo) - only if support for armadiilo library is enabled\n"
                  << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
 }
 
}

