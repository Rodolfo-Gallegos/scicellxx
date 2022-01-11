// IN THIS FILE: Implementation of an abstract class to represent
// vectors

#include "ac_vector.h"

namespace scicellxx
{

 // ===================================================================
 // Empty constructor
 // ===================================================================
 ACVector::ACVector() 
  : NValues(0), Is_own_memory_allocated(false), Delete_vector(true), Is_column_vector(true)
 { }
 
 // ===================================================================
 // Constructor to create an n size zero vector
 // ===================================================================
 ACVector::ACVector(const unsigned long n, bool is_column_vector)
  : NValues(n), Is_own_memory_allocated(false), Delete_vector(true),
    Is_column_vector(is_column_vector)
 { }
 
 // ===================================================================
 // Destructor
 // ===================================================================
ACVector::~ACVector()
{ }
 
}
