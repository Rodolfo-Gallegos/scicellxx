#ifndef CCBACKWARDEULERPREDICTORCORRECTORMETHOD_TPL_H
#define CCBACKWARDEULERPREDICTORCORRECTORMETHOD_TPL_H

#include "ac_predictor_corrector_time_stepper.h"

// Include matrices and vectors
#include "../matrices/matrices.h"

namespace scicellxx
{

 /// @class CCBackwardEulerPCMethod
 /// cc_backward_euler_predictor_corrector_method.h This class
 /// implements Backward Euler method as a predictor corrector to
 /// integrate ODE's.
 template<class EQUATIONS_TYPE>
 class CCBackwardEulerPCMethod : public virtual ACPredictorCorrectorTimeStepper<EQUATIONS_TYPE>
 {
 
 public:

  /// Constructor
  CCBackwardEulerPCMethod();
  
  /// Empty destructor
  virtual ~CCBackwardEulerPCMethod();
  
  /// Applies Backward Eulers method implemented as
  /// Predictor-Corrector to the given odes from the current time "t"
  /// to the time "t+h".  The values of u at time t+h will be stored
  /// at index k (default k = 0).
  void time_step(EQUATIONS_TYPE &odes, const Real h, const Real t,
                 CCData &u, const unsigned k = 0);
  
 protected:
 
  /// Copy constructor (we do not want this class to be
  /// copiable). Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  CCBackwardEulerPCMethod(const CCBackwardEulerPCMethod &copy)
   : ACPredictorCorrectorTimeStepper<EQUATIONS_TYPE>()
  //   : ACTimeStepperForODEs() -- try this HERE
  //   : ACTimeStepper()
   {
    BrokenCopy::broken_copy("CCBackwardEulerPCMethod");
   }
 
  /// Assignment operator (we do not want this class to be
  /// copiable. Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  void operator=(const CCBackwardEulerPCMethod &copy)
   {
    BrokenCopy::broken_assign("CCBackwardEulerPCMethod");
   }
    
 };

}
 
#endif // #ifndef CCBACKWARDEULERPREDICTORCORRECTORMETHOD_TPL_H
