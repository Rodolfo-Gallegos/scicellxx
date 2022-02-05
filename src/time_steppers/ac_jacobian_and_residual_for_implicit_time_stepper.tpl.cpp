#include "ac_jacobian_and_residual_for_implicit_time_stepper.tpl.h"

namespace scicellxx
{
 // ===================================================================
 /// Empty constructor
 // ===================================================================
 template<class EQUATIONS_TYPE>
 ACJacobianAndResidualForImplicitTimeStepper<EQUATIONS_TYPE>::ACJacobianAndResidualForImplicitTimeStepper()
  : ACJacobianAndResidual(),
    ODEs_pt(NULL),
    U_pt(NULL),
    Data_for_jacobian_and_residual_has_been_set(false),
    Jacobian_FY_strategy_pt(NULL)
 {
  
 }
 
 // ===================================================================
 /// Destructor
 // ===================================================================
 template<class EQUATIONS_TYPE>
 ACJacobianAndResidualForImplicitTimeStepper<EQUATIONS_TYPE>::~ACJacobianAndResidualForImplicitTimeStepper()
 {
  // Free the pointer
  Jacobian_FY_strategy_pt = NULL;
 }
 
 // ===================================================================
 /// Set data for Jacobian and residual computation. The odes, the time
 /// step 'h', the current time 't', the values of 'u' and the index
 /// where the values of 'u' at time 't+h' will be stored
 // ===================================================================
 template<class EQUATIONS_TYPE>
 void ACJacobianAndResidualForImplicitTimeStepper<EQUATIONS_TYPE>::
 set_data_for_jacobian_and_residual(EQUATIONS_TYPE *odes_pt, const Real h, const Real t,
                                    CCData *u_pt, const unsigned k)
 {
  // Set the odes
  ODEs_pt = odes_pt;
  
  // Set the time step 
  Time_step = h;
  
  // Set the constant time
  Current_time = t;
  
  // Set the storage of the data
  U_pt = u_pt;
  
  // Set the index of where the values of u at time 't+h' should be
  // stored
  History_index = k;
  
  // Change the flag indicating that the data for the computation of
  // the Jacobian and the residual has been set
  Data_for_jacobian_and_residual_has_been_set = true;
  
 }
 
 // ===================================================================
 /// Set the strategy to compute the ODE's Jacobian
 // ===================================================================
 template<class EQUATIONS_TYPE>
 void ACJacobianAndResidualForImplicitTimeStepper<EQUATIONS_TYPE>::
 set_strategy_for_odes_jacobian(ACJacobianAndResidualForImplicitTimeStepper *jacobian_strategy_for_odes_pt)
 {
  if (jacobian_strategy_for_odes_pt != NULL)
   {
    Jacobian_FY_strategy_pt = jacobian_strategy_for_odes_pt;
   }
  else
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The strategy for the computation of the Jacobian of the ODEs\n"
                  << "that you are setting is not valid (NULL)\n"
                  << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
 }
 
 // =================================================================== 
 /// Get access to the strategy to compute the Jacobian of the ODEs
 // ===================================================================
 template<class EQUATIONS_TYPE>
 ACJacobianAndResidualForImplicitTimeStepper<EQUATIONS_TYPE> *ACJacobianAndResidualForImplicitTimeStepper<EQUATIONS_TYPE>::
 jacobian_FY_strategy_pt()
 {
  if (Jacobian_FY_strategy_pt != NULL)
   {
    return Jacobian_FY_strategy_pt;
   }
  else
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The strategy for the computation of the Jacobian of the ODEs\n"
                  << "is not valid (NULL)\n"
                  << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
 }
 
}

