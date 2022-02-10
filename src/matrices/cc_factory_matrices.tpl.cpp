#include "cc_factory_matrices.tpl.h"

namespace scicellxx
{

 // ===================================================================
 /// Empty constructor
 // ===================================================================
 template<class MAT_TYPE, class VEC_TYPE>
 CCFactoryMatrices<MAT_TYPE, VEC_TYPE>::CCFactoryMatrices()
 { 

 }

 // ===================================================================
 /// Empty destructor
 // ===================================================================
 template<class MAT_TYPE, class VEC_TYPE>
 CCFactoryMatrices<MAT_TYPE, VEC_TYPE>::~CCFactoryMatrices()
 { 

 }
 
 // ===================================================================
 /// Returns a matrix pointer (based on templated parameters)
 // ===================================================================
 template<class MAT_TYPE, class VEC_TYPE>
 MAT_TYPE* CCFactoryMatrices<MAT_TYPE, VEC_TYPE>::create_matrix()
 {
  // Return the matrix
  return new MAT_TYPE();
 }

 // ===================================================================
 /// Returns a matrix pointer (based on compilation options)
 // ===================================================================
 template<class MAT_TYPE, class VEC_TYPE>
 MAT_TYPE* CCFactoryMatrices<MAT_TYPE, VEC_TYPE>::create_matrix(const unsigned long m, const unsigned long n)
 {
  // Return the matrix with allocated memory
  return new MAT_TYPE(m, n);
 }
 
 // ===================================================================
 /// Returns a vector pointer (based on compilation options)
 // ===================================================================
 template<class MAT_TYPE, class VEC_TYPE>
 VEC_TYPE* CCFactoryMatrices<MAT_TYPE, VEC_TYPE>::create_vector()
 {
  // Return the vector
  return new VEC_TYPE();
 }

 // ===================================================================
 /// Returns a vector pointer (based on compilation options)
 // ===================================================================
 template<class MAT_TYPE, class VEC_TYPE>
 VEC_TYPE* CCFactoryMatrices<MAT_TYPE, VEC_TYPE>::create_vector(const unsigned long n, bool is_column_vector)
 {
  // Return the vector with allocated memory
  return new VEC_TYPE(n, is_column_vector);
 }
  
}
