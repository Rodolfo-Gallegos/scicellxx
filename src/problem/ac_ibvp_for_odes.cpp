#include "ac_ibvp_for_odes.h"

namespace scicellxx
{
 
 // ===================================================================
 /// Constructor, sets the ODEs and the time stepper
 // ===================================================================
 ACIBVPForODEs::ACIBVPForODEs(ACODEs *odes_pt, ACTimeStepperForODEs *time_stepper_pt)
  : ACIBVP(),
    ODEs_pt(odes_pt)
 {
  // Add time stepper
  add_time_stepper_for_odes(time_stepper_pt);
  // Get the number of odes
  const unsigned n_odes = odes_pt->n_odes();
  const unsigned n_history_values = time_stepper_pt->n_history_values();
  U_pt = new CCData(n_odes, n_history_values);
 }
 
 // ===================================================================
 /// Destructor
 // ===================================================================
 ACIBVPForODEs::~ACIBVPForODEs()
 {
  // Free memory
  delete U_pt;
  // Set pointer to null
  U_pt = 0;
 }
 
 // ===================================================================
 /// We perform an unsteady solve by default, if you require a
 /// different solving strategy then override this method
 // ===================================================================
 void ACIBVPForODEs::solve()
 {
  // Solve the ODEs
  unsteady_solve();
 }
 
 // ===================================================================
 /// Problem unsteady solve
 // ===================================================================
 void ACIBVPForODEs::unsteady_solve()
 {
  // Call actions before time stepping
  actions_before_time_stepping();
  
  // Get the number of time steppers
  const unsigned n_time_steppers = this->n_time_steppers();
  
  // If there are more than one time stepper then throw an error since
  // this feature has not been tested and it possibly wont work at
  // all, please carefully implement any necessary changes
  if (n_time_steppers > 1)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "This feature has not been fully tested\n"
                  << "Please create a demo driver that uses more than one\n"
                  << "time stepper and implement the necessary changes to\n"
                  << "the code.\n"
                  << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Loop over all the time steppers
  for (unsigned i = 0; i < n_time_steppers; i++)
   {
    const Real t = time(i);
    const Real h = time_step(i);
    
    // Time step (apply the Time stepper to time integrate the ODEs)
    ode_time_stepper_pt(i)->time_step((*ODEs_pt), h, t, (*U_pt));
   }
  
  // Call actions after time stepping
  actions_after_time_stepping();
  
 }
 
 // ===================================================================
 /// Add a time stepper
 //===================================================================
 void ACIBVPForODEs::add_time_stepper_for_odes(ACTimeStepperForODEs *time_stepper_pt,
                                               const Real initial_time,
                                               const Real time_step)
 {
  // Add the time stepper to the time stepper for odes vector
  Time_stepper_for_odes_pt.push_back(time_stepper_pt);
  // Call the parent method to add the time stepper
  add_time_stepper(time_stepper_pt, initial_time, time_step);
 }
 
 // ===================================================================
 /// Read-only access to the time i-th stepper pointer
 // ===================================================================
 ACTimeStepperForODEs *ACIBVPForODEs::ode_time_stepper_pt(const unsigned i) const
 {
#ifdef SCICELLXX_PANIC_MODE
  // Get the number of time steppers
  const unsigned n_time_steppers = Time_stepper_pt.size();
  if (i >= n_time_steppers)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "You are trying to access to a time stepper that is not\n"
                  << "available at the time stepper container\n"
                  << "Maximum index of time steppers in the container: ["<<n_time_steppers<<"]\n"
                  << "The index of the time stepper you want to access: ["<<i<<"]\n"
                  << std::endl;
    throw SciCellxxLibError(error_message.str(),
                            SCICELLXX_CURRENT_FUNCTION,
                            SCICELLXX_EXCEPTION_LOCATION);
   }
#endif // #ifdef SCICELLXX_PANIC_MODE
  return Time_stepper_for_odes_pt[i];
 }
 
}

 
