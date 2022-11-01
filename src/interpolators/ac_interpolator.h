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
#ifndef ACINTERPOLATOR_H
#define ACINTERPOLATOR_H

#include "../general/general.h"

namespace scicellxx
{

 /// @class ACInterpolator ac_interpolator.h

 /// This class implements the interfaces for interpolation methods
 class ACInterpolator
 {
 
 public:

  /// Empty constructor
  ACInterpolator();
  
  /// Empty destructor
  virtual ~ACInterpolator();
  
  /// Does 1D interpolation specifying the data points, the order of
  /// the interpolation and the desired "x" value to interpolate
  virtual Real interpolate_1D(std::vector<Real> &x_points,
                              std::vector<Real> &fx_points,
                              const Real,
                              const unsigned order);
  
  /// Does 1D interpolation specifying the data points, the order of
  /// the interpolation and the desired "x" values to interpolate
  virtual void interpolate_1D(std::vector<Real> &x_points,
                              std::vector<Real> &fx_points,
                              std::vector<Real> &x,
                              std::vector<Real> &fx,
                              const unsigned order);
  
  /// Does 2D interpolation specifying the data points, the order of
  /// the interpolation and the desired "x" value to interpolate
  virtual Real interpolate_2D(std::vector<std::vector<Real> > &x_points,
                              std::vector<Real> &fx_points,
                              std::vector<Real> &x,
                              const unsigned order);
  
  /// Does 2D interpolation specifying the data points, the order of the
  /// interpolation and the desired "x" values to interpolate
  virtual void interpolate_2D(std::vector<std::vector<Real> > &x_points,
                              std::vector<Real> &fx_points,
                              std::vector<std::vector<Real> > &x,
                              std::vector<Real> &fx,
                              const unsigned order);
  
  /// Does 3D interpolation specifying the data points, the order of
  /// the interpolation and the desired "x" value to interpolate
  virtual Real interpolate_3D(std::vector<std::vector<Real> > &x_points,
                              std::vector<Real> &fx_points,
                              std::vector<std::vector<Real> > &x,
                              const unsigned order);
  
  /// Does 3D interpolation specifying the data points, the order of the
  /// interpolation and the desired "x" values to interpolate
  virtual void interpolate_3D(std::vector<std::vector<Real> > &x_points,
                              std::vector<Real> &fx_points,
                              std::vector<std::vector<Real> > &x,
                              std::vector<Real> &fx,
                              const unsigned order);
  
 protected:
  
  /// Copy constructor (we do not want this class to be
  /// copiable). Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  ACInterpolator(const ACInterpolator &copy)
   {
    BrokenCopy::broken_copy("ACInterpolator");
   }

  /// Assignment operator (we do not want this class to be
  /// copiable. Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  void operator=(const ACInterpolator &copy)
   {
    BrokenCopy::broken_assign("ACInterpolator");
   }

 };

}
 
#endif // #ifndef ACINTERPOLATOR_H
