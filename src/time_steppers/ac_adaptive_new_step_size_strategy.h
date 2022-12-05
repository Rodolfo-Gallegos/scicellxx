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
#ifndef ACADAPTIVENEWSTEPSIZESTRATEGY_H
#define ACADAPTIVENEWSTEPSIZESTRATEGY_H

#include "../general/general.h"

namespace scicellxx
{
#ifdef TYPEDEF_REAL_IS_DOUBLE
#define DEFAULT_ADAPTIVE_NEW_STEP_SIZE_MAXIMUM_TOLERANCE 1.0e-3
#define DEFAULT_ADAPTIVE_NEW_STEP_SIZE_MINIMUM_TOLERANCE 1.0e-8
#else
#define DEFAULT_ADAPTIVE_NEW_STEP_SIZE_MAXIMUM_TOLERANCE 1.0e-3
#define DEFAULT_ADAPTIVE_NEW_STEP_SIZE_MINIMUM_TOLERANCE 1.0e-6
#endif // #ifdef TYPEDEF_REAL_IS_DOUBLE
 
 /// @class ACAdaptiveNewStepSizeStrategy ac_adaptive_new_step_size_strategy.h
 
 // ==============================================================
 /// @class ACAdaptiveNewStepSizeStrategy This class implements the
 /// interface for the strategies to compute the new step size in
 /// adaptive time stepping methods
 // ==============================================================
 class ACAdaptiveNewStepSizeStrategy
 {
 public:
  
  // Constructor (empty)
  ACAdaptiveNewStepSizeStrategy();
  
  // Destructor (empty)
  virtual ~ACAdaptiveNewStepSizeStrategy();
  
  // Set default maximum error tolerance
  inline void set_default_maximum_tolerance()
  {Maximum_tolerance = DEFAULT_ADAPTIVE_NEW_STEP_SIZE_MAXIMUM_TOLERANCE;}
  
  // Set default minimum error tolerance
  inline void set_default_minimum_tolerance()
  {Minimum_tolerance = DEFAULT_ADAPTIVE_NEW_STEP_SIZE_MINIMUM_TOLERANCE;}
  
  // Set new maximum error tolerance
  inline void set_new_maximum_tolerance(const Real new_maximum_tolerance)
  {Maximum_tolerance = new_maximum_tolerance;}
  
  // Set new minimum error tolerance
  inline void set_new_minimum_tolerance(const Real new_minimum_tolerance)
  {Minimum_tolerance = new_minimum_tolerance;}
  
  // Read-only maximum error tolerance
  inline const Real maximum_tolerance() const {return Maximum_tolerance;}
  
  // Read-only minimum error tolerance
  inline const Real minimum_tolerance() const {return Minimum_tolerance;}
   
  // The strategy to compute the new step size
  virtual Real new_step_size(const Real local_error, const Real h) = 0;
  
 private:
  
  /// Copy constructor (we do not want this class to be
  /// copiable). Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  ACAdaptiveNewStepSizeStrategy(const ACAdaptiveNewStepSizeStrategy &copy)
   {
    BrokenCopy::broken_copy("ACAdaptiveNewStepSizeStrategy");
   }
  
  /// Assignment operator (we do not want this class to be
  /// copiable. Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  void operator=(const ACAdaptiveNewStepSizeStrategy &copy)
   {
    BrokenCopy::broken_assign("ACAdaptiveNewStepSizeStrategy");
   }
  
  // Maximum error tolerance
  Real Maximum_tolerance;
  
  // Minimum error tolerance
  Real Minimum_tolerance;
  
 };
 
}

#endif // #ifndef ACADAPTIVENEWSTEPSIZESTRATEGY_H
