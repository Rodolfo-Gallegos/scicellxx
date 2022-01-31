#ifndef ACTIMESTEPPERFRORODES_H
#define ACTIMESTEPPERFRORODES_H

#include "../general/general.h"
#include "../data_structures/data_structures.h"
#include "ac_time_stepper.h"

namespace scicellxx
{ 
 /// @class ACTimeStepperForODEs ac_time_stepper_for_odes.h

 /// This class implements the interfaces for integration methods to
 /// solve ODE's
 class ACTimeStepperForODEs : public virtual ACTimeStepper
 {
 
 public:
 
  /// Empty constructor
  ACTimeStepperForODEs();
 
  /// Empty destructor
  virtual ~ACTimeStepperForODEs();
  
  /// Override this function with an error message to indicate the
  /// user needs to specify the odes to work with
  void time_step(const Real h,
                 const Real t,
                 CCData &u,
                 unsigned k = 0)
  {
   // Error message
   std::ostringstream error_message;
   error_message << "You need to pass the ODEs to solve for. Use the version of\n"
                 << "this method that receives the ODEs as an argument\n"
                 << std::endl;
   throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
  }
  
  /// Performs a time step applying a time integration method to the
  /// given odes from the current time "t" to the time "t+h".
  /// Previous the call of the method, the values of u at time "t"
  /// should be stored at index k (default k = 0). After the call, the
  /// values at time "t+h" will be stored at index k, therefore the
  /// values at time "t" will be at index k+1
  virtual void time_step(ACODEs &odes,
                         const Real h,
                         const Real t,
                         CCData &u,
                         unsigned k = 0) = 0;
  
  /// Resets the time stepper to its initial state.
  
  /// For the BDF 2 method we require to re-enable the computation of
  /// the initial guess for the value (t+2h) only the first time that
  /// the methods is called.
  
  /// For ADAPTIVE time steppers we need to indicate no previous "time
  /// step (h)" has been computed. Thus the given time step should be
  /// considered as the initial time step
  virtual void reset() { }
  
 protected:
 
  /// Copy constructor (we do not want this class to be
  /// copiable). Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  ACTimeStepperForODEs(const ACTimeStepperForODEs &copy)
   {
    BrokenCopy::broken_copy("ACTimeStepperForODEs");
   }

  /// Assignment operator (we do not want this class to be
  /// copiable. Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  void operator=(const ACTimeStepperForODEs &copy)
   {
    BrokenCopy::broken_assign("ACTimeStepperForODEs");
   }
  
 };

}
 
#endif // #ifndef ACTIMESTEPPERFRORODES_H
