// IN THIS FILE: The implementation of a concrete class to store and
// work with matrices. This implementation makes use of Armadillo's
// library, thus this is only a wrap for Armadillo's methods

#include "cc_matrix_armadillo.h"

namespace scicellxx
{

 // ===================================================================
 // Empty constructor
 // ===================================================================
 CCMatrixArmadillo::CCMatrixArmadillo()
  : ACMatrix()
 {
  // Delete any data in memory
  clean_up();
 }
 
 // ===================================================================
 // Constructor to create an m X n matrix.
 // ===================================================================
 CCMatrixArmadillo::CCMatrixArmadillo(const unsigned long m, const unsigned long n)
  : ACMatrix(m, n)
 {
  // Allocate memory
  allocate_memory(m, n);
  
  // DELETE DELETE DELETE DELETE
  // // Delete any data in memory
  // clean_up();
  
  // // Create an Armadillo's matrix
  // Arma_matrix_pt = new arma::Mat(m, n);
  // DELETE DELETE DELETE DELETE
 }
 
 // ===================================================================
 // Constructor where we pass the data for the matrix of size m X n
 // ===================================================================
 CCMatrixArmadillo::CCMatrixArmadillo(Real *matrix_pt,
                       const unsigned long m,
                       const unsigned long n)
 : ACMatrix(m, n)
 {
  // Copy the data from the input vector to the Matrix_pt vector
  set_matrix(matrix_pt, m, n);
 }
 
 // ===================================================================
 // Constructor that creates an Armadillo's matrix from a CCMatrix
 // ===================================================================
 CCMatrixArmadillo::CCMatrixArmadillo(CCMatrix &matrix)
  : ACMatrix()
 {
  // Get the pointer to the matrix data
  Real *matrix_pt = matrix.matrix_pt();
  // Get the dimension of the new matrix
  unsigned long m = matrix.n_rows();
  unsigned long n = matrix.n_columns();
  
  // Copy the data from the vector to the Matrix_pt vector
  set_matrix(matrix_pt, m, n);
 }
 
 // ===================================================================
 // Constructor that creates an Armadillo's matrix from a CCVector
 // ===================================================================
 CCMatrixArmadillo::CCMatrixArmadillo(CCVector &vector)
  : ACMatrix()
 {
  // Get the pointer to the vector data
  Real *vector_pt = vector.vector_pt();
  // Compute the dimension of the new matrix by checking whether the
  // vector is transposed or not
  unsigned long m = 0;
  unsigned long n = 0;
  if (!vector.is_column_vector()) // a row vector
   {
    m = 1;
    n = vector.n_values();
   }
  else // a column vector
   {
    m = vector.n_values();
    n = 1;
   }
  // Copy the data from the vector to the Matrix_pt vector
  set_matrix(vector_pt, m, n);
  
 }

 // ===================================================================
 // Constructor that creates an Armadillo's matrix from a CCVectorArmadillo
 // ===================================================================
 CCMatrixArmadillo::CCMatrixArmadillo(CCVectorArmadillo &vector)
 {
  // Compute the dimension of the new matrix by checking whether the
  // vector is transposed or not
  unsigned long m = 0;
  unsigned long n = 0;
  if (!vector.is_column_vector()) // a row vector
   {
    m = 1;
    n = vector.n_values();
   }
  else // a column vector
   {
    m = vector.n_values();
    n = 1;
   }
  set_matrix(vector.arma_vector_pt(), m, n);
 }
 
 // ===================================================================
 // Copy constructor
 // ===================================================================
 CCMatrixArmadillo::CCMatrixArmadillo(const CCMatrixArmadillo &copy)
  : ACMatrix(copy.n_rows(), copy.n_columns())
 {
  // Clean any possible previously allocated memory
  clean_up();
  
  // Call the copy constructor of Armadillo
  Arma_matrix_pt = new arma::Mat<Real>(*(copy.arma_matrix_pt()));
  
  // Mark the matrix as having its own memory
  this->Is_own_memory_allocated = true;
  
 }
 
 // ===================================================================
 // Empty destructor
 // ===================================================================
 CCMatrixArmadillo::~CCMatrixArmadillo()
 {
  // Deallocate memory
  clean_up();
 }
 
 // ===================================================================
 // Assignment operator
 // ===================================================================
 CCMatrixArmadillo& CCMatrixArmadillo::operator=(const CCMatrixArmadillo &source_matrix)
 {
  // Clean-up and set values
  set_matrix(source_matrix.arma_matrix_pt(),
             source_matrix.n_rows(),
             source_matrix.n_columns());
  
  // Return this (de-referenced pointer)
  return *this;
 }
 
 // ===================================================================
 // += operator
 // ===================================================================
 CCMatrixArmadillo& CCMatrixArmadillo::operator+=(const CCMatrixArmadillo &matrix)
 {
  // Call the method to perform the addition
  add_matrix(matrix, *this);
  // Return the solution matrix
  return *this;
 }
 
 // ===================================================================
 // -= operator
 // ===================================================================
 CCMatrixArmadillo& CCMatrixArmadillo::operator-=(const CCMatrixArmadillo &matrix)
 {
  // Call the method to perform the operation
  substract_matrix(matrix, *this);
  // Return the solution matrix
  return *this; 
 }

 // ===================================================================
 // Add operator
 // ===================================================================
 CCMatrixArmadillo CCMatrixArmadillo::operator+(const CCMatrixArmadillo &matrix)
 {
  // Create a zero matrix where to store the result
  CCMatrixArmadillo solution(this->NRows, this->NColumns);
  // Call the method to perform the addition
  add_matrix(matrix, solution);
  // Return the solution matrix
  return solution;
 }

 // ===================================================================
 // Substraction operator
 // ===================================================================
 CCMatrixArmadillo CCMatrixArmadillo::operator-(const CCMatrixArmadillo &matrix)
 {
  // Create a zero matrix where to store the result
  CCMatrixArmadillo solution(this->NRows, this->NColumns);
  // Call the method to perform the operation
  substract_matrix(matrix, solution);
  return solution;
 }

 // ===================================================================
 // Multiplication operator
 // ===================================================================
 CCMatrixArmadillo CCMatrixArmadillo::operator*(const CCMatrixArmadillo &right_matrix)
 { 
  // Create a zero matrix where to store the result
  CCMatrixArmadillo solution(this->NRows, right_matrix.n_columns());
  // Perform the multiplication
  multiply_by_matrix(right_matrix, solution);
  // Return the solution matrix
  return solution;
 }
 
 // ===================================================================
 // Transforms the input vector to an armadillo matrix class type
 // (virtual such that each derived class has to implement it)
 // ===================================================================
 void CCMatrixArmadillo::set_matrix(const Real *matrix_pt,
                                       const unsigned long m,
                                       const unsigned long n)
 {
  // Clean any possible previously allocated memory
  clean_up();
  
  // Set the number of rows and columns
  this->NRows = m;
  this->NColumns = n;
  
  // Create an Armadillo's matrix (makes an own copy of the data,
  // therefore 'matrix_pt' may be deleted safely)
  Arma_matrix_pt = new arma::Mat<Real>(matrix_pt, m, n);
  
  // Mark the matrix as having its own memory
  this->Is_own_memory_allocated = true;
  
 }
 
 // ===================================================================
 // Receives an armadillo type matrix
 // ===================================================================
 void CCMatrixArmadillo::set_matrix(arma::Mat<Real> *arma_matrix_pt,
                                       const unsigned long m,
                                       const unsigned long n)
 {
  // Clean any possible previously allocated memory
  clean_up();
  
  // Set the number of rows and columns
  this->NRows = m;
  this->NColumns = n;
  
  // Call the copy constructor of Armadillo
  Arma_matrix_pt = new arma::Mat<Real>(*arma_matrix_pt);
  
  // Mark the matrix as having its own memory
  this->Is_own_memory_allocated = true;
  
 }
 
 // ===================================================================
 // Clean up for any dynamically stored data
 // ===================================================================
 void CCMatrixArmadillo::clean_up()
 {
  // Check whether the Matrix allocated its own memory
  if (this->Is_own_memory_allocated)
   {
    // Mark the matrix as deleteable
    this->Delete_matrix = true;
    // Free the memory allocated for the matrix
    free_memory_for_matrix();
   }
  else // If empty
   {
    // Set the pointer of the matrix to NULL
    Arma_matrix_pt = 0;
   }
  
 }
 
 // ===================================================================
 // Free allocated memory for matrix
 // ===================================================================
 void CCMatrixArmadillo::free_memory_for_matrix()
 {
  // Is the matrix allowed for deletion. If this method is called from
  // an external source we need to check whether the matrix has been
  // marked for deletion
  if (this->Delete_matrix)
   {
    delete Arma_matrix_pt;
    Arma_matrix_pt = 0; 
    
    // Mark the matrix as not having allocated memory
    this->Is_own_memory_allocated=false;
    
   } // if (Delete_matrix)
  else
   {
    // Error message
    std::ostringstream error_message;
    error_message << "You are trying to free the memory of a matrix that is\n"
                  << "not marked as deletable" << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
 }
 
 // ===================================================================
 // Performs sum of matrices
 // ===================================================================
 void CCMatrixArmadillo::add_matrix(const CCMatrixArmadillo &matrix,
                                       CCMatrixArmadillo &solution_matrix)
 {
  // Check that THIS and the other matrix have memory allocated
  if (!this->Is_own_memory_allocated || !matrix.is_own_memory_allocated())
   {
    // Error message
    std::ostringstream error_message;
    error_message << "One of the matrices to operate with has no memory allocated\n"
                  << "this->Is_own_memory_allocated = "
                  << this->Is_own_memory_allocated << "\n"
                  << "matrix.is_own_memory_allocated() = "
                  << matrix.is_own_memory_allocated() << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Check whether the dimensions of the matrices are the same
  const unsigned long n_rows_input_matrix = matrix.n_rows();
  const unsigned long n_columns_input_matrix = matrix.n_columns();
  const unsigned long n_rows = this->NRows;
  const unsigned long n_columns = this->NColumns;
  if (n_rows != n_rows_input_matrix || n_columns != n_columns_input_matrix)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The dimension of the matrices is not the same:\n"
                  << "dim(matrix) = (" << n_rows_input_matrix << ", "
                  << n_columns_input_matrix << ")\n"
                  << "dim(this) = (" << n_rows << ", " << n_columns
                  << ")\n" << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }

  // Check whether the solution matrix has allocated memory, otherwise
  // allocate it here!!!
  if (!solution_matrix.is_own_memory_allocated())
   {
    // Allocate memory for the matrix
    solution_matrix.allocate_memory(n_rows, n_columns);
   }
  else
   {
    // Check whether the dimension of the solution matrix are correct
    const unsigned long n_rows_solution_matrix = solution_matrix.n_rows();
    const unsigned long n_columns_solution_matrix = solution_matrix.n_columns();
    if (n_rows != n_rows_solution_matrix || n_columns != n_columns_solution_matrix)
     {
      // Error message
      std::ostringstream error_message;
      error_message << "The dimension of the matrices is not the same:\n"
                    << "dim(this) = (" << n_rows << ", "
                    << n_columns << ")\n"
                    << "dim(solution_matrix) = (" << n_rows_solution_matrix
                    << ", " << n_columns_solution_matrix << ")\n" << std::endl;
      throw SciCellxxLibError(error_message.str(),
                             SCICELLXX_CURRENT_FUNCTION,
                             SCICELLXX_EXCEPTION_LOCATION);
     }
    
   }
  
  // Get the matrix pointer of the solution matrix
  arma::Mat<Real> *arma_solution_matrix_pt = solution_matrix.arma_matrix_pt();
  
  // Get the matrix pointer of the input matrix
  arma::Mat<Real> *arma_matrix_pt = matrix.arma_matrix_pt();
  
  // Perform the operation
  (*arma_solution_matrix_pt) = (*Arma_matrix_pt) + (*arma_matrix_pt);
  
 }
 
 // ===================================================================
 // Performs substraction of matrices
 // ===================================================================
 void CCMatrixArmadillo::substract_matrix(const CCMatrixArmadillo &matrix,
                                    CCMatrixArmadillo &solution_matrix)
 {
  // Check that THIS and the other matrix have memory allocated
  if (!this->Is_own_memory_allocated || !matrix.is_own_memory_allocated())
   {
    // Error message
    std::ostringstream error_message;
    error_message << "One of the matrices to operate with has no memory allocated\n"
                  << "this->Is_own_memory_allocated = "
                  << this->Is_own_memory_allocated << "\n"
                  << "matrix.is_own_memory_allocated() = "
                  << matrix.is_own_memory_allocated() << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Check whether the dimensions of the matrices are the same
  const unsigned long n_rows_input_matrix = matrix.n_rows();
  const unsigned long n_columns_input_matrix = matrix.n_columns();
  const unsigned long n_rows = this->NRows;
  const unsigned long n_columns = this->NColumns;
  if (n_rows != n_rows_input_matrix || n_columns != n_columns_input_matrix)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The dimension of the matrices is not the same:\n"
                  << "dim(matrix) = (" << n_rows_input_matrix << ", "
                  << n_columns_input_matrix << ")\n"
                  << "dim(this) = (" << n_rows << ", " << n_columns
                  << ")\n" << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }

  // Check whether the solution matrix has allocated memory, otherwise
  // allocate it here!!!
  if (!solution_matrix.is_own_memory_allocated())
   {
    // Allocate memory for the matrix
    solution_matrix.allocate_memory(n_rows, n_columns);
   }
  else
   {
    // Check whether the dimension of the solution matrix are correct
    const unsigned long n_rows_solution_matrix = solution_matrix.n_rows();
    const unsigned long n_columns_solution_matrix = solution_matrix.n_columns();
    if (n_rows != n_rows_solution_matrix || n_columns != n_columns_solution_matrix)
     {
      // Error message
      std::ostringstream error_message;
      error_message << "The dimension of the matrices is not the same:\n"
                    << "dim(this) = (" << n_rows << ", "
                    << n_columns << ")\n"
                    << "dim(solution_matrix) = (" << n_rows_solution_matrix
                    << ", " << n_columns_solution_matrix << ")\n" << std::endl;
      throw SciCellxxLibError(error_message.str(),
                             SCICELLXX_CURRENT_FUNCTION,
                             SCICELLXX_EXCEPTION_LOCATION);
     }
    
   }
  
  // Get the matrix pointer of the solution matrix
  arma::Mat<Real> *arma_solution_matrix_pt = solution_matrix.arma_matrix_pt();
  
  // Get the matrix pointer of the input matrix
  arma::Mat<Real> *arma_matrix_pt = matrix.arma_matrix_pt();
  
  // Perform the operation
  (*arma_solution_matrix_pt) = (*Arma_matrix_pt) - (*arma_matrix_pt);
  
 }

 // ===================================================================
 // Performs multiplication of matrices
 // ===================================================================
 void CCMatrixArmadillo::multiply_by_matrix(const CCMatrixArmadillo &right_matrix,
                                               CCMatrixArmadillo &solution_matrix)
 {
  // Check that THIS and the right matrix have memory allocated
  if (!this->Is_own_memory_allocated || !right_matrix.is_own_memory_allocated())
   {
    // Error message
    std::ostringstream error_message;
    error_message << "One of the matrices to operate with has no memory allocated\n"
                  << "this->Is_own_memory_allocated = "
                  << this->Is_own_memory_allocated << "\n"
                  << "right_matrix.is_own_memory_allocated() = "
                  << right_matrix.is_own_memory_allocated() << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Check whether the dimensions of the matrices allow for
  // multiplication
  const unsigned long n_rows_right_matrix = right_matrix.n_rows();
  const unsigned long n_columns_right_matrix = right_matrix.n_columns();
  const unsigned long n_rows_left_matrix = this->NRows;
  const unsigned long n_columns_left_matrix = this->NColumns;
  if (n_columns_left_matrix != n_rows_right_matrix)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The dimension of the matrices does not allow "
                  << "multiplication:\n"
                  << "dim(left_matrix) = (" << n_rows_left_matrix << ", "
                  << n_columns_left_matrix << ")\n"
                  << "dim(right_matrix) = (" << n_rows_right_matrix << ", "
                  << n_columns_right_matrix << ")\n" << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }

  // Check whether the solution matrix has allocated memory, otherwise
  // allocate it here!!!
  if (!solution_matrix.is_own_memory_allocated())
   {
    // Allocate memory for the matrix
    solution_matrix.allocate_memory(n_rows_left_matrix, n_columns_right_matrix);
   }
  else
   {
    // Check whether the dimension of the solution matrix are correct
    const unsigned long n_rows_solution_matrix = solution_matrix.n_rows();
    const unsigned long n_columns_solution_matrix = solution_matrix.n_columns();
    if (n_rows_left_matrix != n_rows_solution_matrix ||
        n_columns_right_matrix != n_columns_solution_matrix)
     {
      // Error message
      std::ostringstream error_message;
      error_message << "The dimension of the solution matrix is not appropiate for\n"
                    << "the operation:\n"
                    << "dim(left_matrix) = (" << n_rows_left_matrix << ", "
                    << n_columns_left_matrix << ")\n"
                    << "dim(right_matrix) = (" << n_rows_right_matrix << ", "
                    << n_columns_right_matrix << ")\n"
                    << "dim(solution_matrix) = (" << n_rows_solution_matrix
                    << ", " << n_columns_solution_matrix << ")\n" << std::endl;
      throw SciCellxxLibError(error_message.str(),
                             SCICELLXX_CURRENT_FUNCTION,
                             SCICELLXX_EXCEPTION_LOCATION);
     }
    
   }
  
  // Get the matrix pointer of the solution matrix
  arma::Mat<Real> *arma_solution_matrix_pt = solution_matrix.arma_matrix_pt();
  
  // Get the matrix pointer of the right matrix
  arma::Mat<Real> *arma_right_matrix_pt = right_matrix.arma_matrix_pt();
  
  // Perform the operation
  (*arma_solution_matrix_pt) = (*Arma_matrix_pt) * (*arma_right_matrix_pt);
  
 }
 
 // ===================================================================
 // Computes the transpose and store in the solution matrix
 // ===================================================================
 void CCMatrixArmadillo::transpose(CCMatrixArmadillo &transposed_matrix)
 {
  // Compute transpose
  arma::Mat<Real> arma_transposed_matrix_pt = Arma_matrix_pt->t();
  
  transposed_matrix.set_matrix(&arma_transposed_matrix_pt, arma_transposed_matrix_pt.n_rows, arma_transposed_matrix_pt.n_cols);
 }
 
 // ===================================================================
 // Computes the transpose and stores it in itself
 // ===================================================================
 void CCMatrixArmadillo::transpose()
 {
  // Check that THIS matrix has memory allocated
  if (!this->Is_own_memory_allocated)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "THIS matrix has no memory allocated\n"
                  << "this->Is_own_memory_allocated = "
                  << this->Is_own_memory_allocated << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Transpose itself
  arma::inplace_trans(*(this->Arma_matrix_pt));

  // Update nrows and columns of transposed_matrix
  unsigned long NRows_tmp = this->NRows;
  this->NRows = this->NColumns;
  this->NColumns = NRows_tmp;
  
 }
 
 // ===================================================================
 // Get the specified value from the matrix (read-only)
 // ===================================================================
 const Real CCMatrixArmadillo::value(const unsigned long i, const unsigned long j) const
 {
#ifdef SCICELLXX_RANGE_CHECK
  if (!(this->is_own_memory_allocated()))
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The matrix has no memory allocated"
                  << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  if (i > this->n_rows())
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The row you are trying to access is out of range\n"
                  << "Number of rows: " << this->n_rows() << std::endl
                  << "Requested row: " << i << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  if (j > this->n_columns())
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The column you are trying to access is out of range\n"
                  << "Number of columns: " << this->n_columns() << std::endl
                  << "Requested column: " << j << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
#endif // #ifdef SCICELLXX_RANGE_CHECK
  // Return the value at row i and column j
  return (*Arma_matrix_pt)(i, j);
 }

 // ===================================================================
 // Set values in the matrix (write version)
 // ===================================================================
 Real &CCMatrixArmadillo::value(const unsigned long i, const unsigned long j)
 {
#ifdef SCICELLXX_RANGE_CHECK
  if (!(this->is_own_memory_allocated()))
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The matrix has no memory allocated"
                  << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  if (i > this->n_rows())
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The row you are trying to access is out of range\n"
                  << "Number of rows: " << this->n_rows() << std::endl
                  << "Requested row: " << i << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  if (j > this->n_columns())
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The column you are trying to access is out of range\n"
                  << "Number of columns: " << this->n_columns() << std::endl
                  << "Requested column: " << j << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
#endif // #ifdef SCICELLXX_RANGE_CHECK
  // Return the value at row i and column j
  return (*Arma_matrix_pt)(i, j);
 }
 
 // ===================================================================
 /// Permute the rows in the list
 // ===================================================================
 void CCMatrixArmadillo::permute_rows(std::vector<std::pair<unsigned long, unsigned long> > &permute_list)
 {
  // Get the number of elements in the permute list
  const unsigned n_permute_list = permute_list.size();
  // Check that the rows numbers are within the range of the number of
  // rows of the matrix
  for (unsigned i = 0; i < n_permute_list; i++)
   {
    if (permute_list[i].first >= this->NRows || permute_list[i].second >= this->NRows)
     {
      // Error message
      std::ostringstream error_message;
      error_message << "The rows of the permute list in the " << i << "-th position\n"
                    << "are out of the range of the number of rows of the matrix"
                    << "permute_list[i].first: " << permute_list[i].first
                    << "permute_list[i].second: " << permute_list[i].second << std::endl;
      throw SciCellxxLibError(error_message.str(),
                             SCICELLXX_CURRENT_FUNCTION,
                             SCICELLXX_EXCEPTION_LOCATION);
     }
    
   }

  // ------------------------------
  // Do the permutations
  
  // Cache the number of columns
  const unsigned long n_columns = this->NColumns;
  for (unsigned long k = 0; k < n_columns; k++)
   {
    // Loop over the permute list
    for (unsigned ii = 0; ii < n_permute_list; ii++)
     {
      // Get the i-th row and the j-th row to permute
      const unsigned long i = permute_list[ii].first;
      const unsigned long j = permute_list[ii].second;
      
      // Swap values in rows i and j
      
      // Temporarly storage to swap values
      Real tmp = value(j, k);
      // Swap
      value(j, k) = value(i, k);
      value(i, k) = tmp;
     } // for (ii < n_permute_list)
    
   } // for (k < n_columns)
  
 }
 
 // ===================================================================
 /// Permute the columns in the list
 // ===================================================================
 void CCMatrixArmadillo::permute_columns(std::vector<std::pair<unsigned long, unsigned long> > &permute_list)
 {
  // Get the number of elements in the permute list
  const unsigned n_permute_list = permute_list.size();
  // Check that the columns numbers are within the range of the number
  // of columns of the matrix
  for (unsigned i = 0; i < n_permute_list; i++)
   {
    if (permute_list[i].first >= this->NColumns || permute_list[i].second >= this->NColumns)
     {
      // Error message
      std::ostringstream error_message;
      error_message << "The columns of the permute list in the " << i << "-th position\n"
                    << "are out of the range of the number of columns of the matrix"
                    << "permute_list[i].first: " << permute_list[i].first
                    << "permute_list[i].second: " << permute_list[i].second << std::endl;
      throw SciCellxxLibError(error_message.str(),
                             SCICELLXX_CURRENT_FUNCTION,
                             SCICELLXX_EXCEPTION_LOCATION);
     }
    
   }
  
  // ------------------------------
  // Do the permutations
  
  // Cache the number of rows
  const unsigned long n_rows = this->NRows;
  for (unsigned long k = 0; k < n_rows; k++)
   {
    // Loop over the permute list
    for (unsigned ii = 0; ii < n_permute_list; ii++)
     {
      // Get the i-th column and the j-th column to permute
      const unsigned long i = permute_list[ii].first;
      const unsigned long j = permute_list[ii].second;
      
      // Swap values in columns i and j
      
      // Temporarly storage to swap values
      Real tmp = value(k, j);
      // Swap
      value(k, j) = value(k, i);
      value(k, i) = tmp;
     } // for (ii < n_permute_list)
    
   } // for (k < n_rows)
  
 }
 
 // ===================================================================
 /// Permute rows i and j
 // ===================================================================
 void CCMatrixArmadillo::permute_rows(const unsigned long &i, const unsigned long &j)
 {
  // Check that both input columns are within the range
  if (i < this->NRows && j < this->NRows)
   {
    // Cache the number of columns
    const unsigned long n_columns = this->NColumns;
    // Swap values in columns i and j
    for (unsigned long k = 0; k < n_columns; k++)
     {
      // Temporarly storage to swap values
      Real tmp = value(j, k);
      // Swap
      value(j, k) = value(i, k);
      value(i, k) = tmp;
     }
   }
  else
   {
    // Error message
    std::ostringstream error_message;
    error_message << "One of the selected rows to permute is larger than\n"
                  << "the number of rows of the matrix"
                  << "i: " << i << " j: " << j << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
 }

 // ================================================================== 
 // Permute columns i and j
 // ===================================================================
 void CCMatrixArmadillo::permute_columns(const unsigned long &i, const unsigned long &j)
 {
  // Check that both input columns are within the range
  if (i < this->NColumns && j < this->NColumns)
   {
    // Cache the number of rows
    const unsigned long n_rows = this->NRows;
    // Swap values in columns i and j
    for (unsigned long k = 0; k < n_rows; k++)
     {
      // Temporarly storage to swap values
      Real tmp = value(k, j);
      // Swap
      value(k, j) = value(k, i);
      value(k, i) = tmp;
     }
   }
  else
   {
    // Error message
    std::ostringstream error_message;
    error_message << "One of the selected columns to permute is larger than\n"
                  << "the number of columns of the matrix"
                  << "i: " << i << " j: " << j << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
   
 }
 
 // ===================================================================
 // Output the matrix
 // ===================================================================
 void CCMatrixArmadillo::output(bool output_indexes) const
 {
  if (!this->Is_own_memory_allocated)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The matrix has no memory allocated. It is empty" << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  else
   {
    // Check whether we should output the indexes
    if (output_indexes)
     {
      for (unsigned long i = 0; i < this->NRows; i++)
       {
        for (unsigned long j = 0; j < this->NColumns; j++)
         {
          std::cout << "(" << i << ", " << j << "): "
                    << value(i,j) << std::endl; 
         } // for (j < this->NColumns)
       } // for (i < this->NRows)
     } // if (output_indexes)
   else
     {
      Arma_matrix_pt->print();
     } // else if (output_indexes)
   
   }
 
 }

 // ===================================================================
 // Output the matrix
 // ===================================================================
 void CCMatrixArmadillo::output(std::ofstream &outfile,
                                   bool output_indexes) const
 {
  if (!this->Is_own_memory_allocated)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The matrix has no memory allocated. It is empty" << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  else
   {
    // Check whether we should output the indexes
    if (output_indexes)
     {
      for (unsigned long i = 0; i < this->NRows; i++)
       {
        for (unsigned long j = 0; j < this->NColumns; j++)
         {
          outfile << "(" << i << ", " << j << "): "
                  << value(i,j) << std::endl; 
         } // for (j < this->NColumns)
       } // for (i < this->NRows)
     } // if (output_indexes)
    else
     {
      Arma_matrix_pt->print(outfile);
     } // else if (output_indexes)
  }
      
 }
 
 // ===================================================================
 // Allows to create a matrix with the given size but with no data
 // ===================================================================
 void CCMatrixArmadillo::allocate_memory(const unsigned long m,
                                            const unsigned long n)
 {
  // Clean any possibly stored data
  clean_up();
  
  // Set the number of rows and columns of the matrix
  this->NRows = m;
  this->NColumns = n;
  
  // Allocate memory
  //allocate_memory();
  
  // Allocate memory for the matrix
  Arma_matrix_pt = new arma::Mat<Real>(this->NRows, this->NColumns);
  
  // Mark the matrix as allocated its own memory
  this->Is_own_memory_allocated=true;
  
 }
 
 // // ===================================================================
 // // Allocates memory to store entries of the matrix
 // // ===================================================================
 // void CCMatrixArmadillo::allocate_memory()
 // {
 //  // Delete any data in memory
 //  clean_up();
  
 //  // Allocate memory for the matrix
 //  Arma_matrix_pt = new arma::Mat<Real>(this->NRows, this->NColumns);
  
 //  // Mark the matrix as allocated its own memory
 //  this->Is_own_memory_allocated=true;
 // }
 
 // ===================================================================
 // Fills the matrix with zeroes
 // ===================================================================
 void CCMatrixArmadillo::fill_with_zeroes()
 {
  // Check that the matrix has memory allocated
  if (this->Is_own_memory_allocated)
   {
    // Fill the matrix with zeroes
    Arma_matrix_pt->zeros();
   }
  else
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The matrix has no memory allocated\n"
                  << "this->Is_own_memory_allocated = "
                  << this->Is_own_memory_allocated << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);    
   }
  
 }
 
 // ================================================================
 // Extra methods to work with matrices, we do not need them to be
 // friends of the class since all their operations are performed
 // using the class methods
 // ================================================================
 
 // ===================================================================
 // Performs sum of matrices
 // ===================================================================
 void add_matrices(const CCMatrixArmadillo &matrix_one,
                   const CCMatrixArmadillo &matrix_two,
                   CCMatrixArmadillo &solution_matrix)
 {
  // Check that both matrices have memory allocated
  if (!matrix_one.is_own_memory_allocated() || !matrix_two.is_own_memory_allocated())
   {
    // Error message
    std::ostringstream error_message;
    error_message << "One of the matrices to operate with has no memory allocated\n"
                  << "matrix_one.is_own_memory_allocated() = "
                  << matrix_one.is_own_memory_allocated() << "\n"
                  << "matrix_two.is_own_memory_allocated() = "
                  << matrix_two.is_own_memory_allocated() << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Check whether the dimensions of the matrices are the same
  const unsigned long n_rows_matrix_one = matrix_one.n_rows();
  const unsigned long n_columns_matrix_one = matrix_one.n_columns();
  const unsigned long n_rows_matrix_two = matrix_two.n_rows();
  const unsigned long n_columns_matrix_two = matrix_two.n_columns();
  if (n_rows_matrix_one != n_rows_matrix_two ||
      n_columns_matrix_one != n_columns_matrix_two)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The dimension of the matrices is not the same:\n"
                  << "dim(matrix_one) = (" << n_rows_matrix_one << ", "
                  << n_columns_matrix_one << ")\n"
                  << "dim(matrix_two) = (" << n_rows_matrix_two << ", "
                  << n_columns_matrix_two
                  << ")\n" << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }

  // Check whether the solution matrix has allocated memory, otherwise
  // allocate it here!!!
  if (!solution_matrix.is_own_memory_allocated())
   {
    // Allocate memory for the matrix
    solution_matrix.allocate_memory(n_rows_matrix_one, n_columns_matrix_one);
   }
  else
   {
    // Check whether the dimension of the solution matrix are correct
    const unsigned long n_rows_solution_matrix = solution_matrix.n_rows();
    const unsigned long n_columns_solution_matrix = solution_matrix.n_columns();
    if (n_rows_matrix_one != n_rows_solution_matrix ||
        n_columns_matrix_one != n_columns_solution_matrix)
     {
      // Error message
      std::ostringstream error_message;
      error_message << "The dimension of the matrices is not the same:\n"
                    << "dim(matrix_one) = (" << n_rows_matrix_one << ", "
                    << n_columns_matrix_one << ")\n"
                    << "dim(solution_matrix) = (" << n_rows_solution_matrix
                    << ", " << n_columns_solution_matrix << ")\n" << std::endl;
      throw SciCellxxLibError(error_message.str(),
                             SCICELLXX_CURRENT_FUNCTION,
                             SCICELLXX_EXCEPTION_LOCATION);
     }

   }
  
  // Get the matrix pointer of the solution matrix
  arma::Mat<Real> *arma_solution_matrix_pt = solution_matrix.arma_matrix_pt();
  
  // Get the matrix pointer of the input matrices
  arma::Mat<Real> *arma_matrix_one_pt = matrix_one.arma_matrix_pt();
  arma::Mat<Real> *arma_matrix_two_pt = matrix_two.arma_matrix_pt();
  
  // Perform the addition
  (*arma_solution_matrix_pt) = (*arma_matrix_one_pt) + (*arma_matrix_two_pt);
  
 }
 
 // ===================================================================
 // Performs substraction of matrices
 // ===================================================================
 void substract_matrices(const CCMatrixArmadillo &matrix_one,
                         const CCMatrixArmadillo &matrix_two,
                         CCMatrixArmadillo &solution_matrix)
 {
  // Check that both matrices have memory allocated
  if (!matrix_one.is_own_memory_allocated() || !matrix_two.is_own_memory_allocated())
   {
    // Error message
    std::ostringstream error_message;
    error_message << "One of the matrices to operate with has no memory allocated\n"
                  << "matrix_one.is_own_memory_allocated() = "
                  << matrix_one.is_own_memory_allocated() << "\n"
                  << "matrix_two.is_own_memory_allocated() = "
                  << matrix_two.is_own_memory_allocated() << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Check whether the dimensions of the matrices are the same
  const unsigned long n_rows_matrix_one = matrix_one.n_rows();
  const unsigned long n_columns_matrix_one = matrix_one.n_columns();
  const unsigned long n_rows_matrix_two = matrix_two.n_rows();
  const unsigned long n_columns_matrix_two = matrix_two.n_columns();
  if (n_rows_matrix_one != n_rows_matrix_two ||
      n_columns_matrix_one != n_columns_matrix_two)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The dimension of the matrices is not the same:\n"
                  << "dim(matrix_one) = (" << n_rows_matrix_one << ", "
                  << n_columns_matrix_one << ")\n"
                  << "dim(matrix_two) = (" << n_rows_matrix_two << ", "
                  << n_columns_matrix_two
                  << ")\n" << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Check whether the solution matrix has allocated memory, otherwise
  // allocate it here!!!
  if (!solution_matrix.is_own_memory_allocated())
   {
    // Allocate memory for the matrix
    solution_matrix.allocate_memory(n_rows_matrix_one, n_columns_matrix_one);
   }
  else
   {
    // Check whether the dimension of the solution matrix are correct
    const unsigned long n_rows_solution_matrix = solution_matrix.n_rows();
    const unsigned long n_columns_solution_matrix = solution_matrix.n_columns();
    if (n_rows_matrix_one != n_rows_solution_matrix ||
        n_columns_matrix_one != n_columns_solution_matrix)
     {
      // Error message
      std::ostringstream error_message;
      error_message << "The dimension of the matrices is not the same:\n"
                    << "dim(matrix_one) = (" << n_rows_matrix_one << ", "
                    << n_columns_matrix_one << ")\n"
                    << "dim(solution_matrix) = (" << n_rows_solution_matrix
                    << ", " << n_columns_solution_matrix << ")\n" << std::endl;
      throw SciCellxxLibError(error_message.str(),
                             SCICELLXX_CURRENT_FUNCTION,
                             SCICELLXX_EXCEPTION_LOCATION);
     }

   }
  
  // Get the matrix pointer of the solution matrix
  arma::Mat<Real> *arma_solution_matrix_pt = solution_matrix.arma_matrix_pt();
  
  // Get the matrix pointer of the input matrices
  arma::Mat<Real> *arma_matrix_one_pt = matrix_one.arma_matrix_pt();
  arma::Mat<Real> *arma_matrix_two_pt = matrix_two.arma_matrix_pt();
  
  // Perform the addition
  (*arma_solution_matrix_pt) = (*arma_matrix_one_pt) - (*arma_matrix_two_pt);
  
 }
 
 // ===================================================================
 // Performs multiplication of matrices
 // ===================================================================
 void multiply_matrices(const CCMatrixArmadillo &left_matrix,
                        const CCMatrixArmadillo &right_matrix,
                        CCMatrixArmadillo &solution_matrix)
 {
  // Check that both matrices have memory allocated
  if (!left_matrix.is_own_memory_allocated() || !right_matrix.is_own_memory_allocated())
   {
    // Error message
    std::ostringstream error_message;
    error_message << "One of the matrices to operate with has no memory allocated\n"
                  << "left_matrix.is_own_memory_allocated() = "
                  << left_matrix.is_own_memory_allocated() << "\n"
                  << "right_matrix.is_own_memory_allocated() = "
                  << right_matrix.is_own_memory_allocated() << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Check whether the dimensions of the matrices allow for
  // multiplication
  const unsigned long n_rows_left_matrix = left_matrix.n_rows();
  const unsigned long n_columns_left_matrix = left_matrix.n_columns();
  const unsigned long n_rows_right_matrix = right_matrix.n_rows();
  const unsigned long n_columns_right_matrix = right_matrix.n_columns();
  if (n_columns_left_matrix != n_rows_right_matrix)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The dimension of the matrices does not allow "
                  << "multiplication:\n"
                  << "dim(left_matrix) = (" << n_rows_left_matrix << ", "
                  << n_columns_left_matrix << ")\n"
                  << "dim(right_matrix) = (" << n_rows_right_matrix << ", "
                  << n_columns_right_matrix << ")\n" << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Check whether the solution matrix has allocated memory, otherwise
  // allocate it here!!!
  if (!solution_matrix.is_own_memory_allocated())
   {
    // Allocate memory for the matrix
    solution_matrix.allocate_memory(n_rows_left_matrix, n_columns_right_matrix);
   }
  else
   {
    // Check whether the dimension of the solution matrix are correct
    const unsigned long n_rows_solution_matrix = solution_matrix.n_rows();
    const unsigned long n_columns_solution_matrix = solution_matrix.n_columns();
    if (n_rows_left_matrix != n_rows_solution_matrix ||
        n_columns_right_matrix != n_columns_solution_matrix)
     {
      // Error message
      std::ostringstream error_message;
      error_message << "The dimension of the solution matrix is not appropiate for\n"
                    << "the operation:\n"
                    << "dim(left_matrix) = (" << n_rows_left_matrix << ", "
                    << n_columns_left_matrix << ")\n"
                    << "dim(right_matrix) = (" << n_rows_right_matrix << ", "
                    << n_columns_right_matrix << ")\n"
                    << "dim(solution_matrix) = (" << n_rows_solution_matrix
                    << ", " << n_columns_solution_matrix << ")\n" << std::endl;
      throw SciCellxxLibError(error_message.str(),
                             SCICELLXX_CURRENT_FUNCTION,
                             SCICELLXX_EXCEPTION_LOCATION);
     }
    
   }
  
  // Get the matrix pointer of the solution matrix
  arma::Mat<Real> *arma_solution_matrix_pt = solution_matrix.arma_matrix_pt();
  
  // Get the matrix pointer of the input matrices
  arma::Mat<Real> *arma_matrix_one_pt = left_matrix.arma_matrix_pt();
  arma::Mat<Real> *arma_matrix_two_pt = right_matrix.arma_matrix_pt();
  
  // Perform the addition
  (*arma_solution_matrix_pt) = (*arma_matrix_one_pt) * (*arma_matrix_two_pt);
  
 }
 
 // ================================================================
 // Concatenate matrices horizontally
 // ================================================================
  void concatenate_matrices_horizontally(const CCMatrixArmadillo &left_matrix,
                                         const CCMatrixArmadillo &right_matrix,
                                         CCMatrixArmadillo &concatenated_matrix)
 {
  // Check that both matrices have memory allocated
  if (!left_matrix.is_own_memory_allocated() || !right_matrix.is_own_memory_allocated())
   {
    // Error message
    std::ostringstream error_message;
    error_message << "One of the matrices to operate with has no memory allocated\n"
                  << "left_matrix.is_own_memory_allocated() = "
                  << left_matrix.is_own_memory_allocated() << "\n"
                  << "right_matrix.is_own_memory_allocated() = "
                  << right_matrix.is_own_memory_allocated() << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Check whether the number of rows of the matrices are the same
  const unsigned long n_rows_left_matrix = left_matrix.n_rows();
  const unsigned long n_rows_right_matrix = right_matrix.n_rows();
  if (n_rows_left_matrix != n_rows_right_matrix)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The number of rows of the matrices is not the same:\n"
                  << "n_rows(left_matrix) = (" << n_rows_left_matrix << ")\n"
                  << "n_rows(right_matrix) = (" << n_rows_right_matrix << ")\n"
                  << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Get the number of columns of each matrix and compute the new
  // number of columns for the new matrix
  const unsigned long n_columns_left_matrix = left_matrix.n_columns();
  const unsigned long n_columns_right_matrix = right_matrix.n_columns();
  const unsigned long n_columns_new_concatenated_matrix =
   n_columns_left_matrix + n_columns_right_matrix;
  const unsigned long n_rows_new_concatenated_matrix = n_rows_left_matrix;
  
  // Check whether the concatenated matrix has allocated memory,
  // otherwise allocate it here!!!
  if (!concatenated_matrix.is_own_memory_allocated())
   {
    // Allocate memory for the matrix
    concatenated_matrix.allocate_memory(n_rows_new_concatenated_matrix,
                                        n_columns_new_concatenated_matrix);
   }
  
  // Check whether the dimension of the concatenated matrix is correct
  const unsigned long n_rows_concatenated_matrix = concatenated_matrix.n_rows();
  const unsigned long n_columns_concatenated_matrix = concatenated_matrix.n_columns();
  if (n_rows_concatenated_matrix != n_rows_new_concatenated_matrix ||
      n_columns_concatenated_matrix != n_columns_new_concatenated_matrix)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The dimension of the matrix is not as expected:\n"
                  << "dim(concatenated_matrix) = (" << n_rows_concatenated_matrix << ", "
                  << n_columns_concatenated_matrix << ")\n"
                  << "dim(expected_matrix) = (" << n_rows_new_concatenated_matrix
                  << ", " << n_columns_new_concatenated_matrix << ")\n" << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Get the matrix pointer of the concatenated matrix
  arma::Mat<Real> *arma_concatenated_matrix_pt = concatenated_matrix.arma_matrix_pt();
  
  // Get the matrix pointer of the input matrices
  arma::Mat<Real> *arma_left_matrix_pt = left_matrix.arma_matrix_pt();
  arma::Mat<Real> *arma_right_matrix_pt = right_matrix.arma_matrix_pt();
  
  // Perform the concatenation of rows/horizontal_concatenation
  (*arma_concatenated_matrix_pt) = arma::join_rows((*arma_left_matrix_pt), (*arma_right_matrix_pt));
  
 }
 
 // ================================================================
 // Concatenate matrices vertically
 // ================================================================
  void concatenate_matrices_vertically(const CCMatrixArmadillo &upper_matrix,
                                       const CCMatrixArmadillo &lower_matrix,
                                       CCMatrixArmadillo &concatenated_matrix)
 {
  // Check that both matrices have memory allocated
  if (!upper_matrix.is_own_memory_allocated() || !lower_matrix.is_own_memory_allocated())
   {
    // Error message
    std::ostringstream error_message;
    error_message << "One of the matrices to operate with has no memory allocated\n"
                  << "upper_matrix.is_own_memory_allocated() = "
                  << upper_matrix.is_own_memory_allocated() << "\n"
                  << "lower_matrix.is_own_memory_allocated() = "
                  << lower_matrix.is_own_memory_allocated() << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Check whether the number of columns of the matrices are the same
  const unsigned long n_columns_upper_matrix = upper_matrix.n_columns();
  const unsigned long n_columns_lower_matrix = lower_matrix.n_columns();
  if (n_columns_upper_matrix != n_columns_lower_matrix)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The number of columns of the matrices is not the same:\n"
                  << "n_columns(upper_matrix) = (" << n_columns_upper_matrix << ")\n"
                  << "n_columns(lower_matrix) = (" << n_columns_lower_matrix << ")\n"
                  << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Get the number of rows of each matrix and compute the new number
  // of rows for the new matrix
  const unsigned long n_rows_upper_matrix = upper_matrix.n_rows();
  const unsigned long n_rows_lower_matrix = lower_matrix.n_rows();
  const unsigned long n_rows_new_concatenated_matrix =
   n_rows_upper_matrix + n_rows_lower_matrix;
  const unsigned long n_columns_new_concatenated_matrix = n_columns_upper_matrix;
  
  // Check whether the concatenated matrix has allocated memory,
  // otherwise allocate it here!!!
  if (!concatenated_matrix.is_own_memory_allocated())
   {
    // Allocate memory for the matrix
    concatenated_matrix.allocate_memory(n_rows_new_concatenated_matrix,
                                        n_columns_new_concatenated_matrix);
   }
  
  // Check whether the dimension of the concatenated matrix is correct
  const unsigned long n_rows_concatenated_matrix = concatenated_matrix.n_rows();
  const unsigned long n_columns_concatenated_matrix = concatenated_matrix.n_columns();
  if (n_rows_concatenated_matrix != n_rows_new_concatenated_matrix ||
      n_columns_concatenated_matrix != n_columns_new_concatenated_matrix)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The dimension of the new matrix is not as expected:\n"
                  << "dim(concatenated_matrix) = (" << n_rows_concatenated_matrix << ", "
                  << n_columns_concatenated_matrix << ")\n"
                  << "dim(expected_matrix) = (" << n_rows_new_concatenated_matrix
                  << ", " << n_columns_new_concatenated_matrix << ")\n" << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Get the matrix pointer of the concatenated matrix
  arma::Mat<Real> *arma_concatenated_matrix_pt = concatenated_matrix.arma_matrix_pt();
  
  // Get the matrix pointer of the input matrices
  arma::Mat<Real> *arma_upper_matrix_pt = upper_matrix.arma_matrix_pt();
  arma::Mat<Real> *arma_lower_matrix_pt = lower_matrix.arma_matrix_pt();
  
  // Perform the concatenation of rows/horizontal_concatenation
  (*arma_concatenated_matrix_pt) = arma::join_cols((*arma_upper_matrix_pt), (*arma_lower_matrix_pt));
  
 }
 
 // ================================================================
 // Extra methods to work with vector and matrices operations
 // ================================================================

 // ================================================================
 // Multiply vector times vector (if you want to perform dot product
 // use the dot() method defined in the cc_vector_armadillo.tpl.h file
 // instead)
 // ================================================================
 void multiply_vector_times_vector(const CCVectorArmadillo &left_vector,
                                   const CCVectorArmadillo &right_vector,
                                   CCMatrixArmadillo &solution_matrix)
 {
  // Check that the left and the right vectors have memory allocated
  if (!left_vector.is_own_memory_allocated() || !right_vector.is_own_memory_allocated())
   {
    // Error message
    std::ostringstream error_message;
    error_message << "One of the vectors to operate with has no memory allocated\n"
                  << "left_vector.is_own_memory_allocated() = "
                  << left_vector.is_own_memory_allocated() << "\n"
                  << "right_vector.is_own_memory_allocated() = "
                  << right_vector.is_own_memory_allocated() << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Check whether the dimensions of the vectors allow the operation
  unsigned n_rows_left_vector = 0;
  unsigned n_columns_left_vector = 0;
  if (!left_vector.is_column_vector()) // a row vector
   {
    n_rows_left_vector = 1;
    n_columns_left_vector = left_vector.n_values();
   }
  else // a column vector
   {
    n_rows_left_vector = left_vector.n_values();
    n_columns_left_vector = 1;
   }
  
  unsigned n_rows_right_vector = 0;
  unsigned n_columns_right_vector = 0;
  if (!right_vector.is_column_vector()) // a row vector
   {
    n_rows_right_vector = 1;
    n_columns_right_vector = right_vector.n_values();
   }
  else // a column vector
   {
    n_rows_right_vector = right_vector.n_values();
    n_columns_right_vector = 1;
   }
  
  // Check that the dimension of the vectors allow the operation
  if (n_columns_left_vector != n_rows_right_vector)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The dimension of the vectors does not allow "
                  << "multiplication:\n"
                  << "dim(left_vector) = (" << n_rows_left_vector << ", "
                  << n_columns_left_vector << ")\n"
                  << "dim(right_vector) = (" << n_rows_right_vector << ", "
                  << n_columns_right_vector << ")\n" << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }

  // Check whether the solution matrix has allocated memory, otherwise
  // allocate it here!!!
  if (!solution_matrix.is_own_memory_allocated())
   {
    // Allocate memory for the matrix
    solution_matrix.allocate_memory(n_rows_left_vector, n_columns_right_vector);
   }
  else
   {
    // Check whether the dimension of the solution matrix are correct
    const unsigned long n_rows_solution_matrix = solution_matrix.n_rows();
    const unsigned long n_columns_solution_matrix = solution_matrix.n_columns();
    if (n_rows_left_vector != n_rows_solution_matrix ||
        n_columns_right_vector != n_columns_solution_matrix)
     {
      // Error message
      std::ostringstream error_message;
      error_message << "The dimension of the solution matrix is not appropiate for\n"
                    << "the operation:\n"
                    << "dim(left_vector) = (" << n_rows_left_vector << ", "
                    << n_columns_left_vector << ")\n"
                    << "dim(right_vector) = (" << n_rows_right_vector << ", "
                    << n_columns_right_vector << ")\n"
                    << "dim(solution_matrix) = (" << n_rows_solution_matrix
                    << ", " << n_columns_solution_matrix << ")\n" << std::endl;
      throw SciCellxxLibError(error_message.str(),
                             SCICELLXX_CURRENT_FUNCTION,
                             SCICELLXX_EXCEPTION_LOCATION);
     }
    
   }
  
  // Get the vector pointers of the input vectors
  arma::Mat<Real> *arma_vector_left_pt = left_vector.arma_vector_pt();
  arma::Mat<Real> *arma_vector_right_pt = right_vector.arma_vector_pt();
  
  // Get the matrix pointer to the solution matrix
  arma::Mat<Real> *arma_solution_matrix_pt = solution_matrix.arma_matrix_pt();
  
  // Perform the multiplication
  (*arma_solution_matrix_pt) = (*arma_vector_left_pt) * (*arma_vector_right_pt);
  
 }
 
 // ================================================================
 // Multiply vector times matrix
 // ================================================================
 void multiply_vector_times_matrix(const CCVectorArmadillo &vector,
                                   const CCMatrixArmadillo &matrix,
                                   CCMatrixArmadillo &solution_matrix)
 {
  // Check that the vector and the matrix have memory allocated
  if (!vector.is_own_memory_allocated() || !matrix.is_own_memory_allocated())
   {
    // Error message
    std::ostringstream error_message;
    error_message << "Either the vector or the matrix have no memory allocated\n"
                  << "vector.is_own_memory_allocated() = "
                  << vector.is_own_memory_allocated() << "\n"
                  << "matrix.is_own_memory_allocated() = "
                  << matrix.is_own_memory_allocated() << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Check whether the dimensions of the vector and the matrix allow
  // the operation
  unsigned n_rows_vector = 0;
  unsigned n_columns_vector = 0;
  if (!vector.is_column_vector()) // a row vector
   {
    n_rows_vector = 1;
    n_columns_vector = vector.n_values();
   }
  else // a column vector
   {
    n_rows_vector = vector.n_values();
    n_columns_vector = 1;
   }
  
  unsigned n_rows_matrix = matrix.n_rows();
  unsigned n_columns_matrix = matrix.n_columns();
  
  // Check that the dimension of the vector and the matrix allow the
  // operation
  if (n_columns_vector != n_rows_matrix)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The dimension of the vector and the matrix does\n"
                  << "not allow multiplication:\n"
                  << "dim(vector) = (" << n_rows_vector << ", "
                  << n_columns_vector << ")\n"
                  << "dim(matrix) = (" << n_rows_matrix << ", "
                  << n_columns_matrix << ")\n" << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Check whether the solution matrix has allocated memory, otherwise
  // allocate it here!!!
  if (!solution_matrix.is_own_memory_allocated())
   {
    // Allocate memory for the matrix
    solution_matrix.allocate_memory(n_rows_vector, n_columns_matrix);
   }
  else
   {
    // Check whether the dimension of the solution matrix are correct
    const unsigned long n_rows_solution_matrix = solution_matrix.n_rows();
    const unsigned long n_columns_solution_matrix = solution_matrix.n_columns();
    if (n_rows_vector != n_rows_solution_matrix ||
        n_columns_matrix != n_columns_solution_matrix)
     {
      // Error message
      std::ostringstream error_message;
      error_message << "The dimension of the solution matrix is not appropiate for\n"
                    << "the operation:\n"
                    << "dim(vector) = (" << n_rows_vector << ", "
                    << n_columns_vector << ")\n"
                    << "dim(matrix) = (" << n_rows_matrix << ", "
                    << n_columns_matrix << ")\n"
                    << "dim(solution_matrix) = (" << n_rows_solution_matrix
                    << ", " << n_columns_solution_matrix << ")\n" << std::endl;
      throw SciCellxxLibError(error_message.str(),
                             SCICELLXX_CURRENT_FUNCTION,
                             SCICELLXX_EXCEPTION_LOCATION);
     }
    
   }
  
  // Get the vector and the matrix pointer
  arma::Mat<Real> *arma_vector_pt = vector.arma_vector_pt();
  arma::Mat<Real> *arma_matrix_pt = matrix.arma_matrix_pt();
  arma::Mat<Real> *arma_solution_matrix_pt = solution_matrix.arma_matrix_pt();
  
  // Perform the multiplication
  (*arma_solution_matrix_pt) = (*arma_vector_pt) * (*arma_matrix_pt);
  
 }
 
 // ================================================================
 // Multiply matrix times vector
 // ================================================================
 void multiply_matrix_times_vector(const CCMatrixArmadillo &matrix,
                                   const CCVectorArmadillo &vector,
                                   CCVectorArmadillo &solution_vector)
 {
  // Check that the matrix and the vector have memory allocated
  if (!matrix.is_own_memory_allocated() || !vector.is_own_memory_allocated())
   {
    // Error message
    std::ostringstream error_message;
    error_message << "Either the matrix or the vector have no memory allocated\n"
                  << "matrix.is_own_memory_allocated() = "
                  << matrix.is_own_memory_allocated() << "\n"
                  << "vector.is_own_memory_allocated() = "
                  << vector.is_own_memory_allocated() << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }

  // Check whether the vector is a column vector
  if (!vector.is_column_vector())
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The vector to multiply the matrix is not a column vector\n"
                  << "vector.is_column_vector(): " << vector.is_column_vector()
                  << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Check whether the dimensions of the matrix and the vector allow
  // the operation
  const unsigned n_rows_matrix = matrix.n_rows();
  const unsigned n_columns_matrix = matrix.n_columns();
  
  const unsigned n_rows_vector = vector.n_values();
  const unsigned n_columns_vector = 1;
  
  // Check that the dimension of the vector and the matrix allow the
  // operation
  if (n_columns_matrix != n_rows_vector)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The dimension of the matrix and the vector does\n"
                  << "not allow multiplication:\n"
                  << "dim(matrix) = (" << n_rows_matrix << ", "
                  << n_columns_matrix << ")\n"
                  << "dim(vector) = (" << n_rows_vector << ", "
                  << n_columns_vector << ")\n" << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Check whether the solution vector is a column vector
  if (!solution_vector.is_column_vector())
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The solution vector is not a column vector\n"
                  << "solution_vector.is_column_vector(): "
                  << solution_vector.is_column_vector()
                  << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Compute the dimensions for the solution vector
  const unsigned n_rows_solution_vector = n_rows_matrix;
  const unsigned n_columns_solution_vector = n_columns_vector;
  
  // Check whether the solution vector has allocated memory, otherwise
  // allocate it here!!!
  if (!solution_vector.is_own_memory_allocated())
   {
    // Allocate memory for the matrix
    solution_vector.allocate_memory(n_rows_solution_vector);
   }
  else
   {
    if (solution_vector.n_values() != n_rows_solution_vector)
     {
      // Error message
      std::ostringstream error_message;
      error_message << "The dimension of the solution vector is not appropiate for\n"
                    << "the operation:\n"
                    << "dim(matrix) = (" << n_rows_matrix << ", "
                    << n_columns_matrix << ")\n"
                    << "dim(vector) = (" << n_rows_vector << ", "
                    << n_columns_vector << ")\n"
                    << "dim(solution_vector) = (" << n_rows_solution_vector
                    << ", " << n_columns_solution_vector << ")\n" << std::endl;
      throw SciCellxxLibError(error_message.str(),
                             SCICELLXX_CURRENT_FUNCTION,
                             SCICELLXX_EXCEPTION_LOCATION);
     }
    
   }
  
  // Get the vector pointer of the solution vector
  arma::Mat<Real> *arma_solution_vector_pt = solution_vector.arma_vector_pt();
  
  // Get both the vector and the matrix pointer
  arma::Mat<Real> *arma_vector_pt = vector.arma_vector_pt();
  arma::Mat<Real> *arma_matrix_pt = matrix.arma_matrix_pt();
  
  // Perform the multiplication
  (*arma_solution_vector_pt) = (*arma_matrix_pt) * (*arma_vector_pt);  
  
#if 0
  // Perform the multiplication
  for (unsigned long i = 0; i < n_rows_solution_vector; i++)
   {
    const unsigned offset_vector = i * n_columns_vector;
    // We can skip one loop since there is only one column
    
    // Initialise
    solution_vector_pt[offset_vector] = 0;
    for (unsigned long k = 0; k < n_columns_matrix; k++)
     {
      solution_vector_pt[offset_vector]+=
       matrix(i, k) * vector_pt[k];
     }
   }
#endif // #if 0
  
 }
 
}
