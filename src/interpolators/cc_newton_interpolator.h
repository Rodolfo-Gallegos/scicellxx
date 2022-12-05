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
#ifndef CCNEWTONINTERPOLATOR_H
#define CCNEWTONINTERPOLATOR_H

#include "ac_interpolator.h"

namespace scicellxx
{

 /// @class CCNewtonInterpolator cc_newton_interpolator.h
 /// This class implements the Newton form of the interpolation
 /// polynomial
 class CCNewtonInterpolator : public virtual ACInterpolator
 {
 
 public:

  /// Empty constructor
  CCNewtonInterpolator();

  /// Empty destructor
  ~CCNewtonInterpolator();
  
  /// Does the interpolation specifying the set data points, the order
  /// of the interpolation and the desired "x_interpolate" value to
  /// interpolate. We use Newton's polynomial formula to construct a
  /// given order polynomial and interpolate.
  /// N(n) = b0 + b1(x_interpolate-x0) +
  /// b2(x_interpolate-x0)(x_interpolate-x1) +
  /// b3(x_interpolate-x0)(x_interpolate-x1)(x_interpolate-x2) ...
  Real interpolate_1D(std::vector<Real> &x,
                      std::vector<Real> &fx,
                      const Real x_interpolate,
                      const unsigned order);
  
  /// Does the interpolation specifying the set data points, the order
  /// of the interpolation and the desired "x_interpolate" values to
  /// interpolate. We use Newton's polynomial formula to construct a
  /// given order polynomial and interpolate.
  /// N(n) = b0 + b1(x_interpolate-x0) +
  /// b2(x_interpolate-x0)(x_interpolate-x1) +
  /// b3(x_interpolate-x0)(x_interpolate-x1)(x_interpolate-x2) ...
  void interpolate_1D(std::vector<Real> &x,
                      std::vector<Real> &fx,
                      std::vector<Real> &x_interpolate,
                      std::vector<Real> &fx_interpolated,
                      const unsigned order);
  
 protected:
 
  /// Copy constructor (we do not want this class to be
  /// copiable). Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  CCNewtonInterpolator(const CCNewtonInterpolator &copy)
   : ACInterpolator()
   {
    BrokenCopy::broken_copy("CCNewtonInterpolator");
   }
 
  /// Assignment operator (we do not want this class to be
  /// copiable. Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  void operator=(const CCNewtonInterpolator &copy)
   {
    BrokenCopy::broken_assign("CCNewtonInterpolator");
   }

 };

}
 
#endif // #ifndef CCNEWTONINTERPOLATOR_H
