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
