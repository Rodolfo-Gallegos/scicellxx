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
// with matrices. This implementation makes use of Armadillo's
// library, thus this is only a wrap for Armadillo's methods

// Check whether the class has been already defined
#ifndef CCMATRIXARMADILLO_H
#define CCMATRIXARMADILLO_H

// The parent class
#include "ac_matrix.h"
// We include the cc_matrix and cc_vector include files to deal with
// transformations from CCMatrix and CCVector classes to
// CCMatrixArmadillo
#include "cc_vector.h"
#include "cc_matrix.h"

// We iclude the cc_vector_armadillo files to deal with
// transformations from the CCVectorArmadillo to CCMatrixArmadillo
#include "cc_vector_armadillo.h"

// Add Armadillo's includes
#include <armadillo>

namespace scicellxx
{
 
 // Concrete class to represent matrices
  class CCMatrixArmadillo : public virtual ACMatrix
  {
   
  public:
   
   // Empty constructor
   CCMatrixArmadillo();
   
   // Constructor to create an m X n matrix.
   CCMatrixArmadillo(const unsigned long m, const unsigned long n);
  
   // Constructor where we pass the data for the matrix of size m X n
   CCMatrixArmadillo(Real *matrix_pt, const unsigned long m, const unsigned long n);
   
   // Constructor that creates an Armadillo's matrix from a CCMatrix
   CCMatrixArmadillo(const CCMatrix &matrix);
   
   // Constructor that creates an Armadillo's matrix from a CCVector
   CCMatrixArmadillo(const CCVector &vector);
   
   // Constructor that creates an Armadillo's matrix from a CCVectorArmadillo
   CCMatrixArmadillo(const CCVectorArmadillo &vector);
   
   // Copy constructor (we require to define this if we want to use
   // operators overloading as sum and assignment)
   CCMatrixArmadillo(const CCMatrixArmadillo &copy);
   
   // Destructor
   virtual ~CCMatrixArmadillo();
   
   // Assignment operator
   CCMatrixArmadillo &operator=(const CCMatrixArmadillo &source_matrix);
  
   // += operator
   CCMatrixArmadillo &operator+=(const CCMatrixArmadillo &matrix);
  
   // -= operator
   CCMatrixArmadillo &operator-=(const CCMatrixArmadillo &matrix);
  
   // Add operator
   CCMatrixArmadillo operator+(const CCMatrixArmadillo &matrix);
   
   // Substraction operator
   CCMatrixArmadillo operator-(const CCMatrixArmadillo &matrix);
   
   // Multiplication operator with vector
   CCVectorArmadillo operator*(const CCVectorArmadillo &right_vector);
   
   // Multiplication operator with matrix
   CCMatrixArmadillo operator*(const CCMatrixArmadillo &right_matrix);
   
   // Allows to create a matrix with the given size but with no data
   void allocate_memory(const unsigned long m,
                        const unsigned long n);
   
   // Allocates memory to store entries of the vector
   //void allocate_memory();
   
   // Fills the vector with zeroes
   void fill_with_zeroes();
   
   // Transforms the input vector to an Armadillo's matrix class type
   void set_matrix(const Real *matrix_pt,
                   const unsigned long m,
                   const unsigned long n);
   
   // Receives an armadillo type matrix
   void set_matrix(arma::Mat<Real> *arma_matrix_pt,
                   const unsigned long m,
                   const unsigned long n);
   
   // Clean up for any dynamically stored data
   void clean_up();
  
   // Free allocated memory for matrix
   void free_memory_for_matrix();
   
   // Performs sum of matrices
   void add_matrix(const CCMatrixArmadillo &matrix, CCMatrixArmadillo &solution_matrix);
  
   // Performs substraction of matrices
   void substract_matrix(const CCMatrixArmadillo &matrix, CCMatrixArmadillo &solution_matrix);
   
   // Performs multiplication of matrix times vector
   void multiply_by_vector(const CCVectorArmadillo &right_vector, CCVectorArmadillo &solution_vector);
    
   // Performs multiplication of matrices
   void multiply_by_matrix(const CCMatrixArmadillo &right_matrix, CCMatrixArmadillo &solution_matrix);
   
   // Computes the transpose and store it in the transpose matrix
   void transpose(CCMatrixArmadillo &transposed_matrix);
   
   // Transpose the matrix
   void transpose();
   
   // Get the specified value from the matrix (read-only)
   const Real value(const unsigned long i, const unsigned long j) const;
   
   // Set values in the matrix (write version)
   Real &value(const unsigned long i, const unsigned long j);
   
   /// Permute the rows in the list
   void permute_rows(std::vector<std::pair<unsigned long, unsigned long> > &permute_list);
  
   /// Permute the columns in the list
   void permute_columns(std::vector<std::pair<unsigned long, unsigned long> > &permute_list);
  
   /// Permute rows i and j
   void permute_rows(const unsigned long &i, const unsigned long &j);
   
   /// Permute columns i and j
   void permute_columns(const unsigned long &i, const unsigned long &j);
  
   // Output the matrix
   void output(bool output_indexes = false) const;
  
   // Output to file
   void output(std::ofstream &outfile, bool output_indexes = false) const;
   
   // Get access to the Armadillo's matrix
   inline arma::Mat<Real> *arma_matrix_pt() const {return Arma_matrix_pt;}
   
  protected:
   
   // The Aramadillo's type matrix
   arma::Mat<Real> *Arma_matrix_pt;
   
  };
 
 // ================================================================
 // Extra methods to work with matrices, we do not need them to be
 // friends of the class since all their operations are performed
 // using the class methods
 // ================================================================
 
 // Performs sum of matrices
  void add_matrices(const CCMatrixArmadillo &matrix_one,
                    const CCMatrixArmadillo &matrix_two,
                    CCMatrixArmadillo &solution_matrix);
 
 // Performs substraction of matrices
  void substract_matrices(const CCMatrixArmadillo &matrix_one,
                          const CCMatrixArmadillo &matrix_two,
                          CCMatrixArmadillo &solution_matrix);
 
 // Performs multiplication of matrices
  void multiply_matrices(const CCMatrixArmadillo &left_matrix,
                         const CCMatrixArmadillo &right_matrix,
                         CCMatrixArmadillo &solution_matrix);
 
 // Concatenate matrices horizontally
  void concatenate_matrices_horizontally(const CCMatrixArmadillo &left_matrix,
                                         const CCMatrixArmadillo &right_matrix,
                                         CCMatrixArmadillo &concatenated_matrix);
 
  // Concatenate matrices vertically
  void concatenate_matrices_vertically(const CCMatrixArmadillo &upper_matrix,
                                       const CCMatrixArmadillo &lower_matrix,
                                       CCMatrixArmadillo &concatenated_matrix);
 
 // ================================================================
 // Extra methods to work with vector and matrices operations
 // ================================================================
 
 // Multiply vector times vector (if you want to perform dot product
 // use the dot() method defined in the cc_vector_armadillo.tpl.h file
 // instead)
  void multiply_vector_times_vector(const CCVectorArmadillo &left_vector,
                                    const CCVectorArmadillo &right_vector,
                                    CCMatrixArmadillo &solution_matrix);
 
 // Multiply vector times matrix
  void multiply_vector_times_matrix(const CCVectorArmadillo &vector,
                                    const CCMatrixArmadillo &matrix,
                                    CCMatrixArmadillo &solution_matrix);
 
  // Multiply matrix times vector
  void multiply_matrix_times_vector(const CCMatrixArmadillo &matrix,
                                    const CCVectorArmadillo &vector,
                                    CCVectorArmadillo &solution_vector);
  
}

#endif // #ifndef CCMATRIXARMADILLO_H
