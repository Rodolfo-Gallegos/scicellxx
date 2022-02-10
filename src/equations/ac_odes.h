#ifndef ACODES_H
#define ACODES_H

#include "../general/general.h"

#include "../data_structures/data_structures.h"

namespace scicellxx
{

 /// @class ACODEs ac_odes.h
    
 /// This class implements the interface to the odes to be solved by
 /// the integration methods
 class ACODEs
 {
  
 public:
  
  /// Constructor, sets the number of odes or equations
  ACODEs(const unsigned n_equations);
  
  /// Empty destructor
  virtual ~ACODEs();
  
  /// Gets the number of equations
  inline unsigned n_equations() const {return N_equations;}
  
  /// Gets the vector storing the number of calls to the odes
  const std::vector<unsigned> &ncalls_equations()
  {return N_calls_equation;}
  
  /// Gets the number of calls to an specific equation
  unsigned ncalls_equation(const unsigned i) const
  {return N_calls_equation[i];}
  
  /// Evaluates the system of odes at time 't'. The index k states the
  /// history values that are being evaluated (default current values,
  /// k=0). The evaluation produces results in the vector dudt(i). The
  /// values of the i-th function at previous times are accessible via
  /// u(i,1), u(i,2) and so on.
  virtual void evaluate_time_derivatives(const Real t, CCData &u, CCData &dudt, const unsigned k = 0) = 0;
  
 protected:
   
  /// Copy constructor (we do not want this class to be
  /// copiable). Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  ACODEs(const ACODEs &copy)
   {
    BrokenCopy::broken_copy("ACODEs");
   }
  
  /// Assignment operator (we do not want this class to be
  /// copiable. Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  void operator=(const ACODEs &copy)
   {
    BrokenCopy::broken_assign("ACODEs");
   }
  
  /// The number of odes or equations
  unsigned N_equations;
  
  /// The number of calls for each equation
  std::vector<unsigned> N_calls_equation;
  
 };
 
}

#endif // #ifndef ACODES_TPL_H
