#ifndef CCFACTORYTIMESTEPPERFORODES_H
#define CCFACTORYTIMESTEPPERFORODES_H

// Include the integration methods for odes (time steppers)
#include "ac_time_stepper_for_odes.h"
#include "cc_euler_method.h"
#include "cc_runge_kutta_4_method.h"
#include "cc_backward_euler_predictor_corrector_method.h"
#include "cc_adams_moulton_2_predictor_corrector_method.h"
#include "cc_backward_euler_method.h"
#include "cc_adams_moulton_2_method.h"
#include "cc_bdf_2_method.h"
#include "cc_adaptive_runge_kutta_45F_method.h"
#include "cc_adaptive_runge_kutta_45DP_method.h"

namespace scicellxx
{

 /// @class CCFactoryTimeStepperForODEs cc_factory_time_stepper_for_odes.h

 /// This class implements a factory for the integration methods for
 /// ODEs (time steppers)
 class CCFactoryTimeStepperForODEs
 {
  
 public:
  
  /// Empty constructor
  CCFactoryTimeStepperForODEs();
  
  /// Empty destructor
  virtual ~CCFactoryTimeStepperForODEs();
  
  /// Returns the specified time stepper (integration method)
  ACTimeStepperForODEs* create_time_stepper(std::string time_stepper_name);
  
 protected:
  
  /// Copy constructor (we do not want this class to be
  /// copiable). Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  CCFactoryTimeStepperForODEs(const CCFactoryTimeStepperForODEs &copy)
   {
    BrokenCopy::broken_copy("CCFactoryTimeStepperForODEs");
   }
  
  /// Assignment operator (we do not want this class to be
  /// copiable. Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  void operator=(const CCFactoryTimeStepperForODEs &copy)
   {
    BrokenCopy::broken_assign("CCFactoryTimeStepperForODEs");
   }
 
 };

}
 
#endif // #ifndef CCFACTORYTIMESTEPPERFORODES_H

