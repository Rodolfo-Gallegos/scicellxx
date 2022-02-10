#include "ac_time_stepper.tpl.h"

namespace scicellxx
{

 // ===================================================================
 // Empty constructor
 // ===================================================================
 template<class EQUATIONS_TYPE>
 ACTimeStepper<EQUATIONS_TYPE>::ACTimeStepper()
 { 
  // Initialise the number of history values
  N_history_values = 0;
 }

 // ===================================================================
 // Empty destructor
 // ===================================================================
 template<class EQUATIONS_TYPE>
 ACTimeStepper<EQUATIONS_TYPE>::~ACTimeStepper()
 { 

 }

}
