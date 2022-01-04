// IN THIS FILE: Implementation of an abstract class to represent
// matrices

#include "ac_matrix.h"

namespace scicellxx
{

 // ===================================================================
 // Empty constructor
 // ===================================================================
 ACMatrix<Real>::ACMatrix() 
  : NRows(0), NColumns(0), Is_own_memory_allocated(false), Delete_matrix(true)
 { }
 
 // ===================================================================
 // Constructor to create an m X n zero matrix
 // ===================================================================
 ACMatrix<Real>::ACMatrix(const unsigned long m, const unsigned long n)
  : NRows(m), NColumns(n), Is_own_memory_allocated(false), Delete_matrix(true)
 { }
 
 // ===================================================================
 // Destructor
 // ===================================================================
 ACMatrix<Real>::~ACMatrix()
 { }
 
}
