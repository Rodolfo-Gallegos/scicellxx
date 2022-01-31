#ifndef ACIVPFORODES_H
#define ACIVPFORODES_H

#include "../general/general.h"

#include "../data_structures/data_structures.h"

#include "ac_ibvp.h"

namespace scicellxx
{

 /// @class ACIBVPForODEs ac_ibvp_for_odes.h

 /// This class implements the interfaces for an initial boundary
 /// value problem for ODEs. It specifies a template for solving a
 /// problem. You need to create a class that inherents from this one
 /// to solve a particular initial boundary value problem for ODEs
 class ACIBVPForODEs : public virtual ACIBVP
 {
  
 public:
  
  /// Constructor
  ACIBVPForODEs(ACODEs *odes_pt, ACTimeStepperForODEs *time_stepper_pt);
  
  /// Destructor
  virtual ~ACIBVPForODEs();
    
  /// Complete the problem setup (initial conditions/boundary
  /// conditions/set solver/etc)
  virtual void complete_problem_setup() = 0;
  
  /// Set initial conditions
  virtual void set_initial_conditions() = 0;
  
  /// We perform an unsteady solve by default, if you require a
  /// different solving strategy then override this method
  void solve();
  
  /// Document solution
  virtual void document_solution() = 0;
  
  /// Add an ODE time stepper
  void add_time_stepper_for_odes(ACTimeStepperForODEs *time_stepper_pt,
                                 const Real initial_time = 0,
                                 const Real time_step = 0);
  
  /// Read-only access to the time i-th stepper pointer
  ACTimeStepperForODEs *ode_time_stepper_pt(const unsigned i = 0) const;
  
 protected:
  
  /// Copy constructor (we do not want this class to be
  /// copiable). Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
 ACIBVPForODEs(const ACIBVPForODEs &copy)
  : ACIBVP()
   {
    BrokenCopy::broken_copy("ACIBVP");
   }
  
  /// Assignment operator (we do not want this class to be
  /// copiable. Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  void operator=(const ACIBVPForODEs &copy)
   {
    BrokenCopy::broken_assign("ACIBVPForODEs");
   } 
  
  /// Problem unsteady solve
  void unsteady_solve();
  
  /// An ODEs time steppers vector, possibly to solve a problem with
  /// different time scales
  std::vector<ACTimeStepperForODEs *> Time_stepper_for_odes_pt;
  
  /// The ODEs
  ACODEs *ODEs_pt;
  
 };
 
}

#endif // #ifndef ACIVPFORODES_H

