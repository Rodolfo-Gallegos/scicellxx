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
#ifndef CCODESBASIC3BODY_H
#define CCODESBASIC3BODY_H

// Include SciCell++ libraries
#include "../../../src/scicellxx.h"

// The dimension of the problem, the number of coordinates for the
// 3-bodies
#define DIM 3 // This specialised implementation assumes we are
              // working with 3 bodies in 3 dimensions
#define NBODIES 3

namespace scicellxx
{

 /// \class CCODEsBasic3Body cc_odes_basic_3_body.h
    
 /// This class implements a set of odes associated with the n-body
 /// problem
 class CCODEsBasic3Body : public virtual ACODEs
 {
 
 public:
  
  /// Constructor, sets the number of bodies (this specialised
  // implementation assumes we are working with 3 bodies in 3
  // dimensions)
  CCODEsBasic3Body(const unsigned n_bodies = NBODIES);
  
  /// Empty destructor
  virtual ~CCODEsBasic3Body();
  
  /// Evaluates the system of odes at time "t". The values of the i-th
  /// function at previous times are accessible via u(i,1), u(i,2) and
  /// so on. The evaluation produces results in the vector dudt.
  void evaluate_time_derivatives(const Real t, CCData &u, CCData &dudt, const unsigned k = 0);
  
  // Gets access to the masses vector
  inline const Real m(const unsigned i) const {return M[i];}
  
  // Sets the value of the i-th body
  inline Real &m(const unsigned i) {return M[i];}
  
  // Gets access to the gravity vector
  inline const Real g(const unsigned i) const {return G[i];}
  
  // Sets the value of the i-th body
  inline Real &g(const unsigned i) {return G[i];}
  
 protected:
  
  /// Copy constructor (we do not want this class to be
  /// copiable). Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
 CCODEsBasic3Body(const CCODEsBasic3Body &copy)
  : ACODEs(copy), N_bodies(0)
   {
    BrokenCopy::broken_copy("CCODEsBasic3Body");
   }
  
  /// Assignment operator (we do not want this class to be
  /// copiable. Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  void operator=(const CCODEsBasic3Body &copy)
   {
    BrokenCopy::broken_assign("CCODEsBasic3Body");
   }
  
  // The number of bodies
  const unsigned N_bodies;

  // The masses of the bodies
  std::vector<Real> M;
  
  // The gravitational constant
  std::vector<Real> G;
  
 };
 
}

#endif // #ifndef CCODESBASIC3BODY_H
