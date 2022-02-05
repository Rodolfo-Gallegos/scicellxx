#include "ac_ibvp.tpl.h"

namespace scicellxx
{
 
 // ===================================================================
 /// Constructor, sets the ODEs and the time stepper
 // ===================================================================
 template<class EQUATIONS_TYPE>
 ACIBVP<EQUATIONS_TYPE>::ACIBVP()
  : ACProblem(),
    Time(0.0),
    Time_step(0.0)
 {
  // Create an empty time steppers vector
  Time_stepper_pt.clear();
  // Create an empty time vector
  Time.clear();
  // Create and empty time step vector
  Time_step.clear();
 }
 
 // ===================================================================
 /// Destructor
 // ===================================================================
 template<class EQUATIONS_TYPE>
 ACIBVP<EQUATIONS_TYPE>::~ACIBVP()
 {
  
 }

 // ===================================================================
 /// Add a time stepper
 //===================================================================
 template<class EQUATIONS_TYPE>
 void ACIBVP<EQUATIONS_TYPE>::add_time_stepper(ACTimeStepper<EQUATIONS_TYPE> *time_stepper_pt,
                                               const Real initial_time,
                                               const Real time_step)
 {
  // Add the time stepper to the time stepper vector
  Time_stepper_pt.push_back(time_stepper_pt);
  // Add the initial time step for the time stepper
  Time.push_back(initial_time);
  // Add the time step for the time stepper
  Time_step.push_back(time_step);
 }

 // ===================================================================
 /// Read-only access to the time i-th stepper pointer
 // ===================================================================
 template<class EQUATIONS_TYPE>
 ACTimeStepper<EQUATIONS_TYPE> *ACIBVP<EQUATIONS_TYPE>::time_stepper_pt(const unsigned i) const
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
  return Time_stepper_pt[i];
 }
 
 // ===================================================================
 /// Write access to the current time
 // ===================================================================
 template<class EQUATIONS_TYPE>
 Real &ACIBVP<EQUATIONS_TYPE>::time(const unsigned i)
 {
#ifdef SCICELLXX_PANIC_MODE
  // Get the size of the time container
  const unsigned n_time = Time.size();
  if (i >= n_time)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "You are trying to access to a time that is not\n"
                  << "available at the time container\n"
                  << "Maximum index of the time container: ["<<n_time<<"]\n"
                  << "The index on the container you want to access: ["<<i<<"]\n"
                  << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
#endif // #ifdef SCICELLXX_PANIC_MODE
  return Time[i];
 }
 
 // ===================================================================
 /// Read-only access to the current time
 // ===================================================================
 template<class EQUATIONS_TYPE>
 Real ACIBVP<EQUATIONS_TYPE>::time(const unsigned i) const
 {
#ifdef SCICELLXX_PANIC_MODE
  // Get the size of the time container
  const unsigned n_time = Time.size();
  if (i >= n_time)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "You are trying to access to a time that is not\n"
                  << "available at the time container\n"
                  << "Maximum index of the time container: ["<<n_time<<"]\n"
                  << "The index on the container you want to access: ["<<i<<"]\n"
                  << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
#endif // #ifdef SCICELLXX_PANIC_MODE
  return Time[i];
 }
 
 // ===================================================================
 /// Write access to the current time step
 // ===================================================================
 template<class EQUATIONS_TYPE>
 Real &ACIBVP<EQUATIONS_TYPE>::time_step(const unsigned i)
 {
#ifdef SCICELLXX_PANIC_MODE
  // Get the size of the time step container
  const unsigned n_time_step = Time_step.size();
  if (i >= n_time_step)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "You are trying to access to a time step that is not\n"
                  << "available at the time step container\n"
                  << "Maximum index of the time step container: ["<<n_time_step<<"]\n"
                  << "The index on the container you want to access: ["<<i<<"]\n"
                  << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
#endif // #ifdef SCICELLXX_PANIC_MODE
  return Time_step[i];
 }
 
 // ===================================================================
 /// Read-only access to the current time step
 // ===================================================================
 template<class EQUATIONS_TYPE>
 Real ACIBVP<EQUATIONS_TYPE>::time_step(const unsigned i) const
 {
#ifdef SCICELLXX_PANIC_MODE
  // Get the size of the time step container
  const unsigned n_time_step = Time_step.size();
  if (i >= n_time_step)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "You are trying to access to a time step that is not\n"
                  << "available at the time step container\n"
                  << "Maximum index of the time step container: ["<<n_time_step<<"]\n"
                  << "The index on the container you want to access: ["<<i<<"]\n"
                  << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
#endif // #ifdef SCICELLXX_PANIC_MODE
  return Time_step[i];
 }
 
 // ===================================================================
 /// We perform an unsteady solve by default, if you require a
 /// different solving strategy then override this method
 // ===================================================================
 template<class EQUATIONS_TYPE>
 void ACIBVP<EQUATIONS_TYPE>::solve()
 {
  // Solve the Initial Boundary Value Problem
  unsteady_solve();
 }
 
 // ===================================================================
 /// Problem unsteady solve
 // ===================================================================
 template<class EQUATIONS_TYPE>
 void ACIBVP<EQUATIONS_TYPE>::unsteady_solve()
 {
  // Call actions before time stepping
  actions_before_time_stepping();

  // Perform actions before newton solve
  actions_before_newton_solve();

  // HERE
  // Do newton solve
  
  // Perform actions after newton solve
  actions_after_newton_solve();
  
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
    
    // Time step (apply the Time stepper to time integrate the Equations)
    time_stepper_pt(i)->time_step((*Equations_pt), h, t, (*U_pt));
   }
    
  // Call actions after time stepping
  actions_after_time_stepping();
  
 }
 
}
