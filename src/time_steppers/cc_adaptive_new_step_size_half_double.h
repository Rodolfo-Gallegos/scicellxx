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
#ifndef CCADAPTIVENEWSTEPSIZEHALFDOUBLE_H
#define CCADAPTIVENEWSTEPSIZEHALFDOUBLE_H

#include "ac_adaptive_new_step_size_strategy.h"

namespace scicellxx
{ 
 /// @class CCAdaptiveNewStepSizeHalfDouble cc_adaptive_new_step_size_half_double.h
 
 // ==============================================================
 // @class CCAdaptiveNewStepSizeHalfDouble This class implements a
 // concrete strategy to compute the new step size in adaptive time
 // stepping methods
 // ==============================================================
 class CCAdaptiveNewStepSizeHalfDouble : public virtual ACAdaptiveNewStepSizeStrategy
 {
 public:
  // Constructor (empty)
  CCAdaptiveNewStepSizeHalfDouble();
  
  // Destructor (empty)
  virtual ~CCAdaptiveNewStepSizeHalfDouble();
  
  // The strategy to compute the new step size
  Real new_step_size(const Real local_error, const Real h);
  
 private:
  
  /// Copy constructor (we do not want this class to be
  /// copiable). Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  CCAdaptiveNewStepSizeHalfDouble(const CCAdaptiveNewStepSizeHalfDouble &copy)
   : ACAdaptiveNewStepSizeStrategy()
   {
    BrokenCopy::broken_copy("CCAdaptiveNewStepSizeHalfDouble");
   }
  
  /// Assignment operator (we do not want this class to be
  /// copiable. Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  void operator=(const CCAdaptiveNewStepSizeHalfDouble &copy)
   {
    BrokenCopy::broken_assign("CCAdaptiveNewStepSizeHalfDouble");
   }
  
 };
 
}

#endif // #ifndef CCADAPTIVENEWSTEPSIZEHALFDOUBLE_H
