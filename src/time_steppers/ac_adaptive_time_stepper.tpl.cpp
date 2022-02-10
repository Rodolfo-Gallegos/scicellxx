#include "ac_adaptive_time_stepper.tpl.h"

namespace scicellxx
{
 // ===================================================================
 // Constructor
 // ===================================================================
 template<class EQUATIONS_TYPE>
 ACAdaptiveTimeStepper<EQUATIONS_TYPE>::ACAdaptiveTimeStepper()
  : ACTimeStepper<EQUATIONS_TYPE>(),
    Free_memory_for_new_time_step_strategy(false),
    New_time_step_strategy_has_been_set(false),
    Maximum_iterations(DEFAULT_ADAPTIVE_TIME_STEPPER_MAXIMUM_ITERATIONS),
    Maximum_step_size(DEFAULT_ADAPTIVE_TIME_STEPPER_MAXIMUM_STEP_SIZE),
    Minimum_step_size(DEFAULT_ADAPTIVE_TIME_STEPPER_MINIMUM_STEP_SIZE),
    Maximum_tolerance(DEFAULT_ADAPTIVE_TIME_STEPPER_MAXIMUM_TOLERANCE),
    Minimum_tolerance(DEFAULT_ADAPTIVE_TIME_STEPPER_MINIMUM_TOLERANCE),
    Output_messages(false)
 {
  // Sets the default new time step strategy
  set_default_configuration();
 }
 
 // ===================================================================
 // Empty destructor
 // ===================================================================
 template<class EQUATIONS_TYPE>
 ACAdaptiveTimeStepper<EQUATIONS_TYPE>::~ACAdaptiveTimeStepper()
 {
  clean_up();
 }
 
 // ===================================================================
 // In charge of free memory (if any given to the strategy to compute
 // the new step size)
 // ===================================================================
 template<class EQUATIONS_TYPE>
 void ACAdaptiveTimeStepper<EQUATIONS_TYPE>::clean_up()
 {
  if (Free_memory_for_new_time_step_strategy)
   {
    delete New_time_step_strategy_pt;
   }
  
  // Set pointer to null
  New_time_step_strategy_pt = 0;
  
  // Indicate the memory can not be free for the strategy to compute
  // the new step size
  Free_memory_for_new_time_step_strategy = false;;
  
  // Indicate there is no strategy to compute the new step size
  New_time_step_strategy_has_been_set = false;
  
 }
 
 // ===================================================================
 // Set the default configuration
 // ===================================================================
 template<class EQUATIONS_TYPE>
 void ACAdaptiveTimeStepper<EQUATIONS_TYPE>::set_default_configuration()
 {
  // Call reset (throw any previously automatically computed step
  // size)
  reset();

  // Call clean up (free memory assigned to new step size computation
  // strategies)
  clean_up();

  // Set default strategy for new step size computation
  set_default_new_step_size_strategy();
  
  set_default_maximum_iterations();
  set_default_maximum_step_size();
  set_default_minimum_step_size();
  set_default_maximum_tolerance(); 
  set_default_minimum_tolerance();

  // By default output messages are disabled
  disable_output_messages();
  
 }
 
 // ===================================================================
 // Set the default to compute the new time step
 // ===================================================================
 template<class EQUATIONS_TYPE>
 void ACAdaptiveTimeStepper<EQUATIONS_TYPE>::set_default_new_step_size_strategy()
 {
  // Clean up
  clean_up();
  
  // Create the strategy to compute the new step size
  New_time_step_strategy_pt = new CCAdaptiveNewStepSizeHalfDouble();
  // Configure new step size strategy tolerance values
  New_time_step_strategy_pt->set_new_maximum_tolerance(Maximum_tolerance);
  New_time_step_strategy_pt->set_new_minimum_tolerance(Minimum_tolerance);
  
  // Allow free memory of the strategy to compute the time step
  Free_memory_for_new_time_step_strategy = true;
  // Set new time step strategy has been set
  New_time_step_strategy_has_been_set = true;
  
 }
 
 // ===================================================================
 // Set the strategy to compute the new time step
 // ===================================================================
 template<class EQUATIONS_TYPE>
 void ACAdaptiveTimeStepper<EQUATIONS_TYPE>::set_new_step_size_strategy(ACAdaptiveNewStepSizeStrategy *new_time_step_strategy_pt)
 {
  // Clean up
  clean_up();
  // Set the strategy to compute the new step size
  New_time_step_strategy_pt = new_time_step_strategy_pt;
  // Set new time step strategy has been set
  New_time_step_strategy_has_been_set = true;
  // Do not allow free memory of the strategy to compute the time step
  Free_memory_for_new_time_step_strategy = false;
 }
 
}
