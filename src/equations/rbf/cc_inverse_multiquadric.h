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
#ifndef CCINVERSEMULTIQUADRIC_H
#define CCINVERSEMULTIQUADRIC_H

#include "../../general/general.h"

#include "ac_radial_base_function.h"

namespace scicellxx
{
 /// A class than implements the Inverse Muliquadric RBF
 class CCInverseMultiquadric : public virtual ACRadialBaseFunction
 {
  
 public:
  
  /// Constructor
  CCInverseMultiquadric(Real epsilon = 9.0);
  
  /// Destructor
  virtual ~CCInverseMultiquadric();
  
  /// Radial function
  const Real psi(const Real r);
  
  /// Read-only version
  Real epsilon() const {return Epsilon;}

  /// Write-access
  Real &epsilon() {return Epsilon;}
  
 private:
  
  // Copy constructor (we do not want this class to be
  // copiable. Check
  // http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  CCInverseMultiquadric(const CCInverseMultiquadric &copy)
   {
    BrokenCopy::broken_copy("CCInverseMultiquadric");
   }
  
  // Copy constructor (we do not want this class to be
  // copiable. Check
  // http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  void operator=(const CCInverseMultiquadric &copy)
   {
    BrokenCopy::broken_assign("CCInverseMultiquadric");
   }
  
 protected:

  /// Epsilon
  Real Epsilon;
  
 };
 
}

#endif // #ifndef CCINVERSEMULTIQUADRIC_H
