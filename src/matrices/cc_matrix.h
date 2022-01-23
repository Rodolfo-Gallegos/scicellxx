// IN THIS FILE: The definition of a concrete class to store and work
// with matrices. This is the simplest implementation

// Check whether the class has been already defined
#ifndef CCMATRIX_H
#define CCMATRIX_H

// The parent class
#include "ac_matrix.h"
// We include the cc_vector class to deal with vector-matrices
// operations
#include "cc_vector.h"

namespace scicellxx
{
 
 // Concrete class to represent matrices
 class CCMatrix : public virtual ACMatrix
 {
  
 public:
  
  // Empty constructor
  CCMatrix();
  
  // Constructor to create an m X n matrix.
  CCMatrix(const unsigned long m, const unsigned long n);
  
  // Constructor where we pass the data for the matrix of size m X n
  CCMatrix(Real *matrix_pt, const unsigned long m, const unsigned long n);
  
  // Constructor that creates a matrix from a vector
  CCMatrix(const CCVector &vector);
  
  // Copy constructor (we require to define this if we want to use
  // operators overloading as sum and assignment)
  CCMatrix(const CCMatrix &copy);
  
  // Destructor
  virtual ~CCMatrix();
  
  // Assignment operator
  CCMatrix &operator=(const CCMatrix &source_matrix);
  
  // += operator
  CCMatrix &operator+=(const CCMatrix &matrix);
  
  // -= operator
  CCMatrix &operator-=(const CCMatrix &matrix);
  
  // Add operator
  CCMatrix operator+(const CCMatrix &matrix);
  
  // Substraction operator
  CCMatrix operator-(const CCMatrix &matrix);
  
  // Multiplication operator with vector
  CCVector operator*(const CCVector &right_vector);
  
  // Multiplication operator with matrix
  CCMatrix operator*(const CCMatrix &right_matrix);
  
  // Allows to create a matrix with the given size but with no data
  void allocate_memory(const unsigned long m,
                       const unsigned long n); 
  
  // Allocates memory to store entries of the matrix
  //void allocate_memory();
  
  // Fills the vector with zeroes
  void fill_with_zeroes();
  
  // Transforms the input vector to a matrix class type (virtual such
  // that each derived class has to implement it)
  void set_matrix(const Real *matrix_pt,
                  const unsigned long m,
                  const unsigned long n);
  
  // Clean up for any dynamically stored data
  void clean_up();
  
  // Free allocated memory for matrix
  void free_memory_for_matrix();
  
  // Performs sum of matrices
  void add_matrix(const CCMatrix &matrix, CCMatrix &solution_matrix);
  
  // Performs substraction of matrices
  void substract_matrix(const CCMatrix &matrix, CCMatrix &solution_matrix);
  
  // Performs multiplication of matrices
  void multiply_by_matrix(const CCMatrix &right_matrix, CCMatrix &solution_matrix);
  
  // Performs multiplication of matrix times vector
  void multiply_by_vector(const CCVector &right_vector, CCVector &solution_vector);
  
  // Computes the transpose and store it in the transpose matrix
  void transpose(CCMatrix &transposed_matrix);
  
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
  
  // Get access to the Matrix_pt
  inline Real *matrix_pt() const {return Matrix_pt;}
    
 protected:
    
  // The matrix
  Real *Matrix_pt;
    
 };
 
 // ================================================================
 // Extra methods to work with matrices, we do not need them to be
 // friends of the class since all their operations are performed
 // using the class methods
 // ================================================================
 
 // Performs sum of matrices
 void add_matrices(const CCMatrix &matrix_one,
                   const CCMatrix &matrix_two,
                   CCMatrix &solution_matrix);
 
 // Performs substraction of matrices
 void substract_matrices(const CCMatrix &matrix_one,
                         const CCMatrix &matrix_two,
                         CCMatrix &solution_matrix);
 
 // Performs multiplication of matrices
 void multiply_matrices(const CCMatrix &left_matrix,
                        const CCMatrix &right_matrix,
                        CCMatrix &solution_matrix);
 
 // Concatenate matrices horizontally
 void concatenate_matrices_horizontally(const CCMatrix &left_matrix,
                                        const CCMatrix &right_matrix,
                                        CCMatrix &concatenated_matrix);
 
 // Concatenate matrices vertically
 void concatenate_matrices_vertically(const CCMatrix &upper_matrix,
                                      const CCMatrix &lower_matrix,
                                      CCMatrix &concatenated_matrix);
 
 // ================================================================
 // Extra methods to work with vector and matrices operations
 // ================================================================
 // Multiply vector times matrix
 void multiply_vector_times_matrix(const CCVector &vector, const CCMatrix &matrix,
                                   CCMatrix &solution_matrix);
 
 // Multiply matrix times vector
 void multiply_matrix_times_vector(const CCMatrix &matrix, const CCVector &vector,
                                   CCVector &solution_vector);
 
}

#endif // #ifndef CCMATRIX_H
