// IN THIS FILE: The definition of an abstract class to store and work
// with matrices. Most of the matrices implemented in this library use
// this class as the base class

// Check whether the class has been already defined
#ifndef ACMATRIX_H
#define ACMATRIX_H

#include "../general/common_includes.h"
#include "../general/utilities.h"

#ifdef SCICELLXX_USES_ARMADILLO
// Add Armadillo's includes (only for the arma_matrix_pt() method)
#include <armadillo>
#endif // #ifdef SCICELLXX_USES_ARMADILLO

namespace scicellxx
{
 
 /// @class ACMatrix ac_matrix.h
 
 // Abstract class to represent matrices
 class ACMatrix
 {
  
 public:
  
  // Empty constructor
  ACMatrix();
  
  // Constructor to create an m X n zero matrix
  ACMatrix(const unsigned long m, const unsigned long n);
   
  // Destructor
  virtual ~ACMatrix();
   
  // Allows to create a matrix with the given size but with no data
  virtual void allocate_memory(const unsigned long m,
                               const unsigned long n) = 0;
  
  // Fills the matrix with zeroes
  virtual void fill_with_zeroes() = 0;
   
  // Transforms the input vector to a matrix class type (virtual such
  // that each derived class has to implement it)
  virtual void set_matrix(const Real *matrix_pt,
                          const unsigned long m,
                          const unsigned long n) = 0;
   
  // Clean up for any dynamically stored data
  virtual void clean_up() = 0;
   
  // Free allocated memory for matrix
  virtual void free_memory_for_matrix() = 0;
   
  // Get the specified value from the matrix (read-only)
  virtual const Real value(const unsigned long i, const unsigned long j) const = 0;
   
  // Set values in the matrix (write version)
  virtual Real &value(const unsigned long i, const unsigned long j) = 0;
   
  // Get the specified value from the matrix
  inline Real get_value(const unsigned long i, const unsigned long j) const
  {return value(i,j);}
   
  // Set values in the matrix
  inline void set_value(const unsigned long i, const unsigned long j, Real v)
  {value(i, j) = v;}
   
  /// Get access using brackets as matrix(i,j). Read-only version
  inline virtual Real operator()(const unsigned long &i, 
                              const unsigned long &j) const
  {return value(i, j);}
   
  /// Get access using brackets as matrix(i,j). Read-write version
  inline virtual Real &operator()(const unsigned long &i, 
                               const unsigned long &j)
  {return value(i,j);}
   
  /// Permute the rows in the list
  virtual void permute_rows(std::vector<std::pair<unsigned long, unsigned long> > &permute_list) = 0;
   
  /// Permute the columns in the list
  virtual void permute_columns(std::vector<std::pair<unsigned long, unsigned long> > &permute_list) = 0;
      
  /// Permute rows i and j
  virtual void permute_rows(const unsigned long &i,
                            const unsigned long &j) = 0;
   
  /// Permute columns i and j
  virtual void permute_columns(const unsigned long &i,
                               const unsigned long &j) = 0;
   
  // Output the matrix
  virtual void output(bool output_indexes = false) const = 0;
   
  // Output to file
  virtual void output(std::ofstream &outfile,
                      bool output_indexes = false) const = 0;
   
  // Output the matrix
  inline void print(bool output_indexes = false) const
  {output(output_indexes);}
   
  // Output to file
  inline void print(std::ofstream &outfile,
                    bool output_indexes = false) const
  {output(outfile, output_indexes);}
   
  // Return the number of rows of the matrix
  inline unsigned long n_rows() const {return NRows;}
   
  // Return the number of columns of the matrix
  inline unsigned long n_columns() const {return NColumns;}
   
  // Checks whether the memory of the matrix has been allocated by
  // this class
  inline bool is_own_memory_allocated() const {return Is_own_memory_allocated;}
   
  // Checks whether the matrix is allowed to be deleted
  inline bool delete_matrix() const {return Delete_matrix;}
   
  // Enables the deletion of the matrix by itself
  inline void enable_delete_matrix() {Delete_matrix=true;}
   
  // Disables the deletion of the matrix by itself
  inline void disable_delete_matrix() {Delete_matrix=false;}
   
  // Get access to the Matrix_pt
  virtual Real *matrix_pt() const
  {
   // Error message
   std::ostringstream error_message;
   error_message << "Virtual function to resolve matrix pointer, should be\n"
                 << "implemented in derived class" << std::endl;
   throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
  }
   
#ifdef SCICELLXX_USES_ARMADILLO
  // Get access to the Armadillo's matrix
  virtual arma::Mat<Real> *arma_matrix_pt() const
  {
   // Error message
   std::ostringstream error_message;
   error_message << "Virtual function to resolve armadillo matrix pointer, should be\n"
                 << "implemented in derived class if you want to use the armadillo solver\n"
                 << std::endl;
   throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
  }
#endif // #ifdef SCICELLXX_USES_ARMADILLO
   
 protected:
   
  // The size of the matrix
  unsigned long NRows;
  unsigned long NColumns;
   
  // Flag to indicate whether the memory of the matrix has been
  // allocated by this class
  bool Is_own_memory_allocated;
   
  // Flag to indicate whether to delete (free) the allocated memory
  // for the matrix. For example when the matrix is transformed to an
  // specific matrix type (Armadillo matrix, SuperLU matrix, Trilinos
  // matrix) we need to deallocate the memory used for THIS matrix to
  // avoid having multiple copies of it. The deletion of the matrix
  // is true by default.
  bool Delete_matrix;
   
 private:
   
  // Copy constructor (we do not want this class to be copiable). Check
  // http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  ACMatrix(const ACMatrix &copy)
   {
    BrokenCopy::broken_copy("ACMatrix");
   }
   
  // Assignment operator (we do not want this class to be
  // copiable. Check
  // http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  //ACMatrix& operator=(const ACMatrix &copy)
  void operator=(const ACMatrix &copy)
   {
    BrokenCopy::broken_assign("ACMatrix");
   }
   
 };
 
}

#endif // #ifndef ACMATRIX_H
