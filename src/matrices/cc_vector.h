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
// IN THIS FILE: The definition of a concrete class to store and work
// with vectors. This is the simplest implementation

// Check whether the class has been already defined
#ifndef CCVECTOR_H
#define CCVECTOR_H

// The parent class
#include "ac_vector.h"

namespace scicellxx
{
 
 // Concrete class to represent vectors
  class CCVector : public virtual ACVector
  {
   
  public:
   
   // Empty constructor
   CCVector();
   
   // Constructor to create an n size zero vector (we assume vectors
   // are created as column vectors, if you need a row vector then
   // pass "false" as the second parameter)
   CCVector(const unsigned long n, bool is_column_vector = true);
   
   // Constructor where we pass the data for the vector of size n.
   CCVector(Real *vector_pt, const unsigned long n, bool is_column_vector = true);
   
   // Copy constructor (we require to define this if we want to use
   // operators overloading as sum and assignment)
   CCVector(const CCVector &copy);
   
   // Destructor
   virtual ~CCVector();
   
   // Assignment operator
   CCVector& operator=(const CCVector &source_vector);
   
   // += operator
   CCVector& operator+=(const CCVector &vector);
   
   // -= operator
   CCVector& operator-=(const CCVector &vector);
   
   // Add operator
   CCVector operator+(const CCVector &vector);
   
   // Substraction operator
   CCVector operator-(const CCVector &vector);
   
   // Element by element multipliation
   CCVector operator*(const CCVector &vector);
   
   // Allows to create a vector with the given size but with no data 
   void allocate_memory(const unsigned long n);
   
   // Allocates memory to store entries of the vector
   //void allocate_memory();
   
   // Fills the vector with zeroes
   void fill_with_zeroes();
   
   // Performs dot product with the current vector
   Real dot(const CCVector &right_vector);
   
   // Transforms the input vector to a vector class type
   void set_vector(const Real *vector_pt,
                   const unsigned long n,
                   bool is_column_vector = true);
   
   // Clean up for any dynamically stored data
   void clean_up();
   
   // Free allocated memory for vector
   void free_memory_for_vector();
   
   // Performs sum of vectors
   void add_vector(const CCVector &vector, CCVector &solution_vector);
   
   // Performs substraction of vectors
   void substract_vector(const CCVector &vector,
                         CCVector &solution_vector);
   
   // Performs multiplication of vectors (one by one entries)
   void multiply_element_by_element_vector(const CCVector &vector,
                                           CCVector &solution_vector);
   
   // Computes the transpose and store it in the transpose vector
   void transpose(CCVector &transposed_vector);
   
   // Transpose the vector
   inline void transpose()
   {
    this->Is_column_vector=!(this->Is_column_vector);
   }
   
   // Get the specified value from the vector (read-only)
   const Real value(const unsigned long i) const;
   
   // Set values in the vector (write version)
   Real &value(const unsigned long i);
   
   // Output the vector
   void output(bool output_indexes = false) const ;
   
   // Output to file
   void output(std::ofstream &outfile, bool output_indexes = false) const;
   
   // Get access to the Vector_pt
   inline Real *vector_pt() const {return Vector_pt;}
   
   // Computes the norm-1 of the vector
   Real norm_1();
   
   // Computes the norm-2 of the vector
   Real norm_2();

   // Computes the infinite norm
   Real norm_inf();
   
   // Computes the maximum value
   Real max();
   
   // Computes the minimum value
   Real min();
   
  protected:
      
   // The vector
   Real *Vector_pt;
   
  };
  
  // ================================================================
  // Extra methods to work with vectors, we do not need them to be
  // friends of the class since all their operations are performed
  // using the class methods
  // ================================================================
  
  // Dot product of vectors
  Real dot_vectors(const CCVector &left_vector, const CCVector &right_vector);
  
  // Addition of vectors
  void add_vectors(const CCVector &vector_one,
                   const CCVector &vector_two,
                   CCVector &solution_vector);
  
  // Substraction of vectors
  void substract_vectors(const CCVector &vector_one,
                         const CCVector &vector_two,
                         CCVector &solution_vector);
 
  // Performs multiplication of vectors (one by one entries)
  void multiply_element_by_element_vectors(const CCVector &vector_one,
                                           const CCVector &vector_two,
                                           CCVector &solution_vector);
  
}

#endif // #ifndef CCVECTOR_H
