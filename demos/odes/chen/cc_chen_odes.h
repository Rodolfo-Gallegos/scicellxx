#ifndef CCCHENODES_H
#define CCCHENODES_H

// Include SciCell++ libraries
#include "../../../src/scicellxx.h"

namespace scicellxx
{

 /// \class CCChenODEs cc_chen_odes.h
    
 /// This class implements the Chen ODEs
 ///
 /// \frac{du_{1}}{dt} = a*(u_{2} - u_{1})
 /// \frac{du_{2}}{dt} = (c-a)*u_{1} - u_{1}*u_{3} + c*u_{2}
 /// \frac{du_{3}}{dt} = u_{1}*u_{2} - b*u_{3}
 class CCChenODEs : public virtual ACODEs
 {
 
 public:
  
  /// Constructor
  CCChenODEs(Real _a, Real _b, Real _c);
  
  /// Empty destructor
  virtual ~CCChenODEs();
  
  /// Evaluates the system of odes at time 't', using the history
  /// values of u at index k
  void evaluate_time_derivatives(const Real t, CCData &u, CCData &dudt, const unsigned k = 0);
  
 protected:
  
  // Specific ODEs parameters
  // ------------------------------------
  Real a;
  Real b;
  Real c;
  
 };
 
}

#endif // #ifndef CCCHENODES_H
