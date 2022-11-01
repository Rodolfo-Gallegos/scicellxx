//LIC// ====================================================================
//LIC// This file forms part of SciCell++, an object-oriented, 
//LIC// framework for the the simulation of biological and physical
//LIC// phenomena modelled as continuous or discrete processes.
//LIC// 
//LIC// You can find a copy at https://github.com/tachidok/scicellxx
//LIC// 
//LIC//    Version 0.6.0
//LIC//
//LIC// 31/10/2022
//LIC// 
//LIC// SciCell++ Copyright (C) 2016-2022 Julio César Pérez Sansalvador
//LIC// 
//LIC// This framework is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU GENERAL PUBLIC LICENSE
//LIC// published by the Free Software Foundation; either version 3 of
//LIC// the License, or (at your option) any later version.
//LIC// 
//LIC// This framework is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//LIC// GNU GENERAL PUBLIC LICENSE for more details.
//LIC// 
//LIC// You should have received a copy of the GNU GENERAL PUBLIC LICENSE
//LIC// along with this framework; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The author may be contacted at jcp.sansalvador@inaoep.mx
//LIC// 
//LIC// ====================================================================
#ifndef ACTIMESTEPPER_TPL_H
#define ACTIMESTEPPER_TPL_H

#include "../general/general.h"

#include "../data_structures/data_structures.h"

namespace scicellxx
{ 
 /// @class ACTimeStepper ac_time_stepper.tpl.h

 /// This class implements the interfaces for time
 /// steppers/integration methods
 template<class EQUATIONS_TYPE>
 class ACTimeStepper
 {
 
 public:
 
  /// Empty constructor
  ACTimeStepper();
  
  /// Empty destructor
  virtual ~ACTimeStepper();
  
  /// Performs a time step applying a time integration method to the
  /// given equations from the current time "t" to the time
  /// "t+h". Previous to the call of the method, the values of u at
  /// time "t" should be stored at index k (default k = 0). After the
  /// call, the values at time "t+h" will be stored at index k,
  /// therefore the values at time "t" will be at index k+1
  virtual void time_step(EQUATIONS_TYPE &equations,
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
  
  /// Get the associated number of history values (each method is in
  /// charge of setting this value based on the number of history
  /// values it requires)
  inline unsigned n_history_values() const {return N_history_values;}
  
 protected:
 
  /// Copy constructor (we do not want this class to be
  /// copiable). Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  ACTimeStepper(const ACTimeStepper &copy)
   {
    BrokenCopy::broken_copy("ACTimeStepper");
   }

  /// Assignment operator (we do not want this class to be
  /// copiable. Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  void operator=(const ACTimeStepper &copy)
   {
    BrokenCopy::broken_assign("ACTimeStepper");
   }
  
  /// The number of history values
  unsigned N_history_values;
  
 };

}
 
#endif // #ifndef ACTIMESTEPPER_TPL_H
