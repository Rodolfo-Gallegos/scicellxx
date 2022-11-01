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
#ifndef CCFACTORYMATRICES_TPL_H
#define CCFACTORYMATRICES_TPL_H

#include "../general/general.h"

#include "../matrices/ac_vector.h"
#include "../matrices/ac_matrix.h"

#include "../matrices/cc_vector.h"
#include "../matrices/cc_matrix.h"
#ifdef SCICELLXX_USES_ARMADILLO
#include "../matrices/cc_vector_armadillo.h"
#include "../matrices/cc_matrix_armadillo.h"
#endif // #ifdef SCICELLXX_USES_ARMADILLO

namespace scicellxx
{
 
 /// @class CCFactoryMatrices cc_factory_matrices.h
 
 /// This class implements a factory for the instantiation of
 /// matrices. The class is templated so that it will return matrices
 /// and vectors of the templated typeit helps to choose the right
 /// matrices types based on compilation flags
#ifdef SCICELLXX_USES_ARMADILLO
 template<class MAT_TYPE = CCMatrixArmadillo, class VEC_TYPE = CCVectorArmadillo>
#else
  template<class MAT_TYPE = CCMatrix, class VEC_TYPE = CCVector>
#endif
 class CCFactoryMatrices
 {
  
 public:
  
  /// Empty constructor
  CCFactoryMatrices();
  
  /// Empty destructor
  virtual ~CCFactoryMatrices();
  
  /// Returns a matrix pointer (based on compilation options)
  MAT_TYPE* create_matrix();
  
  /// Returns a matrix pointer (based on compilation options)
  MAT_TYPE* create_matrix(const unsigned long m, const unsigned long n);
  
  /// Returns a vector pointer (based on compilation options)
  VEC_TYPE* create_vector();
   
  /// Returns a vector pointer (based on compilation options)
  VEC_TYPE* create_vector(const unsigned long n, bool is_column_vector = true);
  
 protected:
   
  /// Copy constructor (we do not want this class to be
  /// copiable). Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  CCFactoryMatrices(const CCFactoryMatrices &copy)
   {
    BrokenCopy::broken_copy("CCFactoryMatrices");
   }
   
  /// Assignment operator (we do not want this class to be
  /// copiable. Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  void operator=(const CCFactoryMatrices &copy)
   {
    BrokenCopy::broken_assign("CCFactoryMatrices");
   }
   
 };
 
}
 
#endif // #ifndef CCFACTORYMATRICES_TPL_H

