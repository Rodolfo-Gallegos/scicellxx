#include "ac_odes.h"

namespace scicellxx
{

 /// ===================================================================
 /// Constructor, sets the number of odes
 /// ===================================================================
 ACODEs::ACODEs(const unsigned n_equations)
  : N_equations(n_equations)
 {
 
  // Resize the container storing the number of calls to each
  // ODE or equation. Initialise it to zero
  N_calls_equation.resize(n_equations, 0);
 
 }

 /// ===================================================================
 /// Destructor
 /// ===================================================================
 ACODEs::~ACODEs()
 {
  // Clear the storage of the vector
  N_calls_equation.clear();
 
 }

}
