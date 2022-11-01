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
#ifndef ACPROBLEM_H
#define ACPROBLEM_H

#include "../general/general.h"

#include "../time_steppers/time_steppers.h"

namespace scicellxx
{

 /// @class ACProblem ac_problem.h

 /// This class implements the interfaces for the problem class. It
 /// specifies a template for solving a problem thus one needs to
 /// create a class that inherents from this one to solve a particular
 /// problem
 class ACProblem
 {
  
 public:

  /// Constructor, in charge of initialising any stuff required for
  /// the framework
  ACProblem();
  
  /// Destructor
  virtual ~ACProblem();

  /// Complete the problem setup (things you could not do in the
  /// construtor as initial conditions/boundary conditions/solver/etc)
  virtual void complete_problem_setup() = 0;
  
  /// Every derived class must implement its own solve method (calling
  /// the corresponding steady_solve() and unsteady_solve() methods)
  virtual void solve() = 0;
    
  /// Document solution
  virtual void document_solution()
  {
   // Error message
   std::ostringstream error_message;
   error_message << "Virtual function in ACProblem class, you should implement\n"
                 << "it to document your solution" << std::endl;
   throw SciCellxxLibError(error_message.str(),
                          SCICELLXX_CURRENT_FUNCTION,
                          SCICELLXX_EXCEPTION_LOCATION);
  }
  
  /// Write access to the current time step
  inline unsigned &output_file_index() {return Output_file_index;}
  
  /// Read-only access to the current time step
  inline unsigned output_file_index() const {return Output_file_index;}
  
 protected:
  
  /// Copy constructor (we do not want this class to be
  /// copiable). Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  ACProblem(const ACProblem &copy)
   {
    BrokenCopy::broken_copy("ACProblem");
   }
  
  /// Assignment operator (we do not want this class to be
  /// copiable. Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  void operator=(const ACProblem &copy)
   {
    BrokenCopy::broken_assign("ACProblem");
   }
  
  /// Problem steady solve
  virtual void steady_solve() = 0;
  
  /// Problem unsteady solve
  virtual void unsteady_solve() = 0;
  
  /// The set of actions to be performed before newton solve
  virtual void actions_before_newton_solve() = 0;
  
  /// The set of actions to be performed after newton solve
  virtual void actions_after_newton_solve() = 0;
  
  /// A counter to store the current output file index
  unsigned Output_file_index;
  
 };
 
}

#endif // #ifndef ACPROBLEM_H
