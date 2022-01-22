// IN THIS FILE: Implementation of a concrete class to represent
// vectors. This is the simplest implementation

#include "cc_vector.h"

namespace scicellxx
{

 // ===================================================================
 // Empty constructor
 // ===================================================================
 CCVector::CCVector() 
  : ACVector()
 {
  // Delete any data in memory
  clean_up();
 }
 
 // ===================================================================
 // Constructor to create an n size zero vector (we assume vectors
 // are created as column vectors, if you need a row vector then
 // pass "false" as the second parameter)
 // ===================================================================
 CCVector::CCVector(const unsigned long n, bool is_column_vector)
  : ACVector(n, is_column_vector)
 {
  // Allocate memory
  allocate_memory(n);
 }
 
 // ===================================================================
 // Constructor where we pass the data for the vector of size n
 // ===================================================================
 CCVector::CCVector(Real *vector_pt, const unsigned long n,
                    bool is_column_vector)
  : ACVector(n, is_column_vector)
 {
  // Copy the data from the input vector to the Vector_pt vector
  set_vector(vector_pt, n, is_column_vector);
 }
 
 // ===================================================================
 // Copy constructor
 // ===================================================================
 CCVector::CCVector(const CCVector &copy)
  : ACVector(copy.n_values(), copy.is_column_vector())
 {
  // Copy the data from the input vector to the Vector_pt vector
  set_vector(copy.vector_pt(), this->NValues, copy.is_column_vector());
 }
 
 // ===================================================================
 // Empty destructor
 // ===================================================================
 CCVector::~CCVector()
 {
  // Deallocate memory
  clean_up();
 }
 
 // ===================================================================
 // Assignment operator
 // ===================================================================
 CCVector& CCVector::operator=(const CCVector &source_vector)
 {
  // Clean-up and set values
  set_vector(source_vector.vector_pt(),
             source_vector.n_values(),
             source_vector.is_column_vector());
  // Return this (de-referenced pointer)
  return *this;
  
 }
 
 // ===================================================================
 // += operator
 // ===================================================================
 CCVector& CCVector::operator+=(const CCVector &vector)
 {  
  // Call the method to perform the addition
  add_vector(vector, *this);
  // Return the solution vector
  return *this;
 }
 
 // ===================================================================
 // -= operator
 // ===================================================================
 CCVector& CCVector::operator-=(const CCVector &vector)
 {
  // Call the method to perform the operation
  substract_vector(vector, *this);
  // Return the solution vector
  return *this; 
 }
 
 // ===================================================================
 // Add operator
 // ===================================================================
 CCVector CCVector::operator+(const CCVector &vector)
 {
  // Create a zero vector where to store the result
  CCVector solution(this->NValues);
  // Call the method to perform the addition
  add_vector(vector, solution);
  // Return the solution vector
  return solution;
 }
 
 // ===================================================================
 // Substraction operator
 // ===================================================================
 CCVector CCVector::operator-(const CCVector &vector)
 {
  // Create a zero vector where to store the result
  CCVector solution(this->NValues);
  // Call the method to perform the operation
  substract_vector(vector, solution);
  return solution;
 }

 // ===================================================================
 // Element by element multipliation
 // ===================================================================
 CCVector CCVector::operator*(const CCVector &vector)
 {
  // Create a zero vector where to store the result
  CCVector solution_vector(this->NValues);
  multiply_element_by_element_vector(vector, solution_vector);
  return solution_vector;
 }

 // ===================================================================
 // Performs dot product with the current vector
 // ===================================================================
 Real CCVector::dot(const CCVector &right_vector)
 {
  return dot_vectors(*this, right_vector);
 }
  
 // ===================================================================
 // Transforms the input vector to a vector class type (virtual such
 // that each derived class has to implement it)
 // ===================================================================
 void CCVector::set_vector(const Real *vector_pt,
                              const unsigned long n,
                              bool is_column_vector)
 {
  // Clean any possible previously allocated memory
  clean_up();
  
  // Set the number of values
  this->NValues = n;
  
  // Allocate memory for the vector
  Vector_pt = new Real[n];
  
  // Mark the vector as allocated its own memory
  this->Is_own_memory_allocated = true;
  
  // Copy the vector (an element by element copy, uff!!)
  std::memcpy(Vector_pt, vector_pt, n*sizeof(Real));
  
  // Set the transposed status
  this->set_as_column_vector(is_column_vector);
  
 }
 
 // ===================================================================
 // Clean up for any dynamically stored data
 // ===================================================================
 void CCVector::clean_up()
 {
  // Check whether the Vector allocated its own memory
  if (this->Is_own_memory_allocated)
   {
    // Mark the vector as deleteable
    this->Delete_vector = true;
    // Free the memory allocated for the vector
    free_memory_for_vector();
   }
  else // if empty
   {
    // Set the pointer of the vector to NULL
    Vector_pt = 0;
   }
  
 }
 
 // ===================================================================
 // Free allocated memory for vector
 // ===================================================================
 void CCVector::free_memory_for_vector()
 {
  // Is the vector allowed for deletion. If this method is called from
  // an external source we need to check whether the vector has been
  // marked for deletion
  if (this->Delete_vector)
   {
    delete [] Vector_pt;
    Vector_pt = 0; 
    
    // Mark the vector as not having memory allocated
    this->Is_own_memory_allocated=false;
    
   } // if (Delete_vector)
  else
   {
    // Error message
    std::ostringstream error_message;
    error_message << "You are trying to free the memory of a vector that is\n"
                  << "not marked as deletable" << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
 
 }
 
 // ===================================================================
 // Performs sum of vectors
 // ===================================================================
 void CCVector::add_vector(const CCVector &vector,
                              CCVector &solution_vector)
 {
  add_vectors(*this, vector, solution_vector);
 }
 
 // ===================================================================
 // Performs substraction of vectors
 // ===================================================================
 void CCVector::substract_vector(const CCVector &vector,
                                    CCVector &solution_vector)
 {
  substract_vectors(*this, vector, solution_vector);
 }

 // ===================================================================
 // Performs multiplication of vectors (element by element)
 // ===================================================================
 void CCVector::
 multiply_element_by_element_vector(const CCVector &vector,
                                    CCVector &solution_vector)
 {
  multiply_element_by_element_vectors(*this, vector, solution_vector);
 }
 
 // ===================================================================
 // Computes the transpose and store in the transposed vector
 // ===================================================================
 void CCVector::transpose(CCVector &transposed_vector)
 {
  // Check that THIS vector has memory allocated
  if (!this->Is_own_memory_allocated)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "THIS vector has no memory allocated\n"
                  << "this->Is_own_memory_allocated = "
                  << this->Is_own_memory_allocated << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Copy the vector into the tranposed vector
  transposed_vector = (*this);
  // Get the current "column vector" status of the vector and set the
  // transposed status of the new vector
  transposed_vector.set_as_column_vector(!(this->Is_column_vector));
 }
 
 // ===================================================================
 // Get the specified value from the vector (read-only)
 // ===================================================================
 const Real CCVector::value(const unsigned long i) const
 {
#ifdef SCICELLXX_RANGE_CHECK
  if (!(this->is_own_memory_allocated()))
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The vector has no memory allocated"
                  << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  if (i > this->n_values())
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The entry you are trying to access is out of range\n"
                  << "Number of values: " << this->n_values() << std::endl
                  << "Requested entry: " << i << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   } 
#endif // #ifdef SCICELLXX_RANGE_CHECK
  // Return the value at position i
  return Vector_pt[i];
 }
 
 // ===================================================================
 // Set values in the vector (write version)
 // ===================================================================
 Real &CCVector::value(const unsigned long i)
 {
#ifdef SCICELLXX_RANGE_CHECK
  if (!(this->is_own_memory_allocated()))
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The vector has no memory allocated"
                  << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  if (i > this->n_values())
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The entry you are trying to access is out of range\n"
                  << "Number of values: " << this->n_values() << std::endl
                  << "Requested entry: " << i << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   } 
#endif // #ifdef SCICELLXX_RANGE_CHECK
  // Return the value at row i and column j
  return Vector_pt[i];
 }

 // ===================================================================
 // Output the vector
 // ===================================================================
 void CCVector::output(bool output_indexes) const
 {
  if (!this->Is_own_memory_allocated)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The vector has no memory allocated. It is empty" << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  else
   {
    // Check whether we should output the indexes
    if (output_indexes)
     {
      for (unsigned long i = 0; i < this->NValues; i++)
       {
        std::cout << "(" << i << "): " << value(i)
                  << std::endl; 
       } // for (i < this->NValues)
     } // if (output_indexes)
    else
     {
      for (unsigned long i = 0; i < this->NValues; i++)
       {
        std::cout << value(i) << " ";
       } // for (i < this->NValues)
      std::cout << std::endl;
     } // else if (output_indexes)
   }
  
 }
 
 // ===================================================================
 // Output the vector
 // ===================================================================
 void CCVector::output(std::ofstream &outfile, bool output_indexes) const
 {
  if (!this->Is_own_memory_allocated)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The vector has no memory allocated. It is empty" << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  else
   {
    // Check whether we should output the indexes
    if (output_indexes)
     {
      for (unsigned long i = 0; i < this->NValues; i++)
       {
        outfile << "(" << i << "): " << value(i)
                << std::endl; 
       } // for (i < this->NValues)
     } // if (output_indexes)
    else
     {
      for (unsigned long i = 0; i < this->NValues; i++)
       {
        outfile << value(i) << " ";
       } // for (i < this->NValues)
      outfile << std::endl;
     } // else if (output_indexes)
   }
  
 }

 // ===================================================================
 // Computes the norm-1 of the vector
 // ===================================================================
 Real CCVector::norm_1()
 {
  // Sum
  Real sum = 0.0;
  // Check whether the vector has memory allocated
  if (this->Is_own_memory_allocated)
   {
    // Compute the norm
    for (unsigned long i = 0; i < this->NValues; i++)
     {
      sum+= std::fabs(Vector_pt[i]);
     }
   }
  else
   {
    // Error message
    std::ostringstream error_message;
    error_message << "We can not compute the norm of a vector with no memory allocated\n"
                  << "this->Is_own_memory_allocated = "
                  << this->Is_own_memory_allocated << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  return sum;
  
 }
 
 // ===================================================================
 // Computes the norm-2 of the vector
 // ===================================================================
 Real CCVector::norm_2()
 {
  // Sum
  Real sum = 0.0;
  // Check whether the vector has memory allocated
  if (this->Is_own_memory_allocated)
   {
    // Compute the norm 2
    for (unsigned long i = 0; i < this->NValues; i++)
     {
      sum+= (Vector_pt[i]*Vector_pt[i]);
     }
   }
  else
   {
    // Error message
    std::ostringstream error_message;
    error_message << "We can not compute the norm of a vector with no memory allocated\n"
                  << "this->Is_own_memory_allocated = "
                  << this->Is_own_memory_allocated << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  return sqrt(sum);
  
 }

 // ===================================================================
 // Computes the infinite norm
 // ===================================================================
 Real CCVector::norm_inf()
 {
  // Maximum
  Real norm = 0.0;
  // Check whether the vector has memory allocated
  if (this->Is_own_memory_allocated)
   {
    // Initialise the maximum with the first value
    norm = std::fabs(Vector_pt[0]);
    // Compute the norm
    for (unsigned long i = 1; i < this->NValues; i++)
     {
      if (std::fabs(Vector_pt[i]) > norm)
       norm = std::fabs(Vector_pt[i]);
     }
   }
  else
   {
    // Error message
    std::ostringstream error_message;
    error_message << "We can not compute the infinite norm of a vector with no memory allocated\n"
                  << "this->Is_own_memory_allocated = "
                  << this->Is_own_memory_allocated << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  return norm;
  
 }
 
 // ===================================================================
 // Computes the maximum value
 // ===================================================================
 Real CCVector::max()
 {
  // Maximum
  Real max = 0.0;
  // Check whether the vector has memory allocated
  if (this->Is_own_memory_allocated)
   {
    // Initialise the maximum with the first value
    max = Vector_pt[0];
    // Compute the maximum
    for (unsigned long i = 1; i < this->NValues; i++)
     {
      if (Vector_pt[i] > max)
       max = Vector_pt[i];
     }
   }
  else
   {
    // Error message
    std::ostringstream error_message;
    error_message << "We can not compute the maximum of a vector with no memory allocated\n"
                  << "this->Is_own_memory_allocated = "
                  << this->Is_own_memory_allocated << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  return max;
  
 }
 
 // ===================================================================
 // Computes the minimum value
 // ===================================================================
 Real CCVector::min()
 {
  // Minimum
  Real min = 0.0;
  // Check whether the vector has memory allocated
  if (this->Is_own_memory_allocated)
   {
    // Initialise the minimum with the first value
    min = Vector_pt[0];
    // Compute the minimum
    for (unsigned long i = 1; i < this->NValues; i++)
     {
      if (Vector_pt[i] < min)
       min = Vector_pt[i];
     }
   }
  else
   {
    // Error message
    std::ostringstream error_message;
    error_message << "We can not compute the minimum of a vector with no memory allocated\n"
                  << "this->Is_own_memory_allocated = "
                  << this->Is_own_memory_allocated << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  return min;
  
 }
 
 // ===================================================================
 // Allows to create a vector with the given size but with no data
 // ===================================================================
 void CCVector::allocate_memory(const unsigned long n)
 {
  // Clean any possibly stored data
  clean_up();
  
  // Set the number of rows and columns of the matrix
  this->NValues = n;

  //allocate_memory();
  
  // Allocate memory for the vector
  Vector_pt = new Real[this->NValues];
  
  // Mark the vector as allocated its own memory
  this->Is_own_memory_allocated=true;
  
 }
 
 // // ===================================================================
 // // Allocates memory to store entries of the vector
 // // ===================================================================
 // void CCVector::allocate_memory()
 // {
 //  // Delete any data in memory
 //  clean_up();
  
 //  // Allocate memory for the vector
 //  Vector_pt = new Real[this->NValues];
  
 //  // Mark the vector as allocated its own memory
 //  this->Is_own_memory_allocated=true;
 // }
 
 // ===================================================================
 // Fills the vector with zeroes
 // ===================================================================
 void CCVector::fill_with_zeroes()
 {
  // Check that the vector has memory allocated
  if (this->Is_own_memory_allocated)
   {
    // Fill the vector with zeroes
    std::memset(Vector_pt, 0, this->NValues*sizeof(Real));
   }
  else
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The vector has no memory allocated\n"
                  << "this->Is_own_memory_allocated = "
                  << this->Is_own_memory_allocated << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);    
   }
  
 }
 
 // ================================================================
 // Extra methods to work with vectors, we do not need them to be
 // friends of the class since all their operations are performed
 // using the class methods
 // ================================================================
 
 // ================================================================
 // Dot product of vectors
 // ================================================================
 Real dot_vectors(const CCVector &left_vector, const CCVector &right_vector)
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
  const unsigned long n_values_left_vector = left_vector.n_values();
  const unsigned long n_values_right_vector = right_vector.n_values();
  if (n_values_left_vector != n_values_right_vector)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The dimension of the vectors is not the same:\n"
                  << "dim(left_vector) = (" << n_values_left_vector << ")\n"
                  << "dim(right_vector) = (" << n_values_right_vector << ")\n"
                  << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Check that the left vector is a row vector and that the right
  // vector is a column vector
  if (left_vector.is_column_vector())
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The left vector should be a row vector\n"
                  << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  if (!right_vector.is_column_vector())
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The right vector should be a column vector\n"
                  << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Get the vector pointer of the left vector
  Real *left_vector_pt = left_vector.vector_pt();
  // Get the vector pointer of the right vector
  Real *right_vector_pt = right_vector.vector_pt();
  
  // Store the dot product of the vectors
  Real dot_product = 0.0;
  
  // Compute the dot product
  for (unsigned long i = 0; i < n_values_left_vector; i++)
   {
    dot_product+= left_vector_pt[i] * right_vector_pt[i];
   }
  
  return dot_product;
  
 }

 // ================================================================
 // Addition of vectors
 // ================================================================
 void add_vectors(const CCVector &vector_one,
                  const CCVector &vector_two,
                  CCVector &solution_vector)
 {
  // Check that the vectors have memory allocated
  if (!vector_one.is_own_memory_allocated() || !vector_two.is_own_memory_allocated())
   {
    // Error message
    std::ostringstream error_message;
    error_message << "One of the vectors to operate with has no memory allocated\n"
                  << "vector_one.is_own_memory_allocated() = "
                  << vector_one.is_own_memory_allocated() << "\n"
                  << "vector_two.is_own_memory_allocated() = "
                  << vector_two.is_own_memory_allocated() << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Check whether the dimensions of the vectors are the same
  const unsigned long n_values_vector_one = vector_one.n_values();
  const unsigned long n_values_vector_two = vector_two.n_values();
  if (n_values_vector_one != n_values_vector_two)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The dimension of the vectors is not the same:\n"
                  << "dim(vector_one) = (" << n_values_vector_one << ")\n"
                  << "dim(vector_two) = (" << n_values_vector_two << ")\n"
                  << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Check that the three vectors have the same column vector status
  if (vector_one.is_column_vector() != vector_two.is_column_vector() ||
      solution_vector.is_column_vector() != vector_two.is_column_vector())
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The three vectors MUST BE either column or row vectors\n"
                  << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Check whether the solution vector has allocated memory, otherwise
  // allocate it here!!!
  if (!solution_vector.is_own_memory_allocated())
   {
    // Allocate memory for the vector
    solution_vector.allocate_memory(n_values_vector_one);
   }
  else
   {
    // Check that the already allocated memory is correct (n_values)
    if (solution_vector.n_values() != vector_one.n_values())
     {
      // Error message
      std::ostringstream error_message;
      error_message << "The number of elements for the solution vector is\n"
                    << "not as expected\n"
                    << "dim(solution_vector): " << solution_vector.n_values()
                    << "\ndim(vector_one): " << vector_one.n_values()
                    << std::endl;
      throw SciCellxxLibError(error_message.str(),
                             SCICELLXX_CURRENT_FUNCTION,
                             SCICELLXX_EXCEPTION_LOCATION);
     }
    
   }
  
  // Get the vector pointer of the solution vector
  Real *solution_vector_pt = solution_vector.vector_pt();
  
  // Get the vector pointer of the vector one
  Real *vector_one_pt = vector_one.vector_pt();
  // Get the vector pointer of the vector two
  Real *vector_two_pt = vector_two.vector_pt();
  
  // Perform the addition
  for (unsigned long i = 0; i < n_values_vector_one; i++)
   {
    solution_vector_pt[i] = vector_one_pt[i] + vector_two_pt[i];
   }
  
 }

 // ================================================================
 // Substraction of vectors
 // ================================================================
 void substract_vectors(const CCVector &vector_one,
                        const CCVector &vector_two,
                        CCVector &solution_vector)
 {
  // Check that the vectors have no memory allocated
  if (!vector_one.is_own_memory_allocated() || !vector_two.is_own_memory_allocated())
   {
    // Error message
    std::ostringstream error_message;
    error_message << "One of the vectors to operate with has no memory allocated\n"
                  << "vector_one.is_own_memory_allocated() = "
                  << vector_one.is_own_memory_allocated() << "\n"
                  << "vector_two.is_own_memory_allocated() = "
                  << vector_two.is_own_memory_allocated() << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Check whether the dimensions of the vectors are the same
  const unsigned long n_values_vector_one = vector_one.n_values();
  const unsigned long n_values_vector_two = vector_two.n_values();
  if (n_values_vector_one != n_values_vector_two)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The dimension of the vectors is not the same:\n"
                  << "dim(vector_one) = (" << n_values_vector_one << ")\n"
                  << "dim(vector_two) = (" << n_values_vector_two << ")\n"
                  << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }

  // Check that the three vectors have the same column vector status
  if (vector_one.is_column_vector() != vector_two.is_column_vector() ||
      solution_vector.is_column_vector() != vector_two.is_column_vector())
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The three vectors MUST BE either column or row vectors\n"
                  << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Check whether the solution vector has allocated memory, otherwise
  // allocate it here!!!
  if (!solution_vector.is_own_memory_allocated())
   {
    // Allocate memory for the vector
    solution_vector.allocate_memory(n_values_vector_one);
   }
  else
   {
    // Check that the already allocated memory is correct (n_values)
    if (solution_vector.n_values() != vector_one.n_values())
     {
      // Error message
      std::ostringstream error_message;
      error_message << "The number of elements for the solution vector is\n"
                    << "not as expected\n"
                    << "dim(solution_vector): " << solution_vector.n_values()
                    << "\ndim(vector_one): " << vector_one.n_values()
                    << std::endl;
      throw SciCellxxLibError(error_message.str(),
                             SCICELLXX_CURRENT_FUNCTION,
                             SCICELLXX_EXCEPTION_LOCATION);
     }
    
   }
  
  // Get the vector pointer of the solution vector
  Real *solution_vector_pt = solution_vector.vector_pt();
    
  // Get the vector pointer of the vector one
  Real *vector_one_pt = vector_one.vector_pt();
  // Get the vector pointer of the vector two
  Real *vector_two_pt = vector_two.vector_pt();
  
  // Perform the substraction of vectors
  for (unsigned long i = 0; i < n_values_vector_one; i++)
   {
    solution_vector_pt[i] = vector_one_pt[i] - vector_two_pt[i];
   }
  
 }

 // ================================================================
 // Performs multiplication of vectors (one by one entries)
 // ================================================================
 void multiply_element_by_element_vectors(const CCVector &vector_one,
                                          const CCVector &vector_two,
                                          CCVector &solution_vector)
 {
  // Check that the vectors have memory allocated
  if (!vector_one.is_own_memory_allocated() || !vector_two.is_own_memory_allocated())
   {
    // Error message
    std::ostringstream error_message;
    error_message << "One of the vectors to operate with has no memory allocated\n"
                  << "vector_one.is_own_memory_allocated() = "
                  << vector_one.is_own_memory_allocated() << "\n"
                  << "vector_two.is_own_memory_allocated() = "
                  << vector_two.is_own_memory_allocated() << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Check whether the dimensions of the vectors are the same
  const unsigned long n_values_vector_one = vector_one.n_values();
  const unsigned long n_values_vector_two = vector_two.n_values();
  if (n_values_vector_one != n_values_vector_two)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The dimension of the vectors is not the same:\n"
                  << "dim(vector_one) = (" << n_values_vector_one << ")\n"
                  << "dim(vector_two) = (" << n_values_vector_two << ")\n"
                  << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }

  // Check that the three vectors have the same column vector status
  if (vector_one.is_column_vector() != vector_two.is_column_vector() ||
      solution_vector.is_column_vector() != vector_two.is_column_vector())
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The three vectors MUST BE either column or row vectors\n"
                  << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Check whether the solution vector has allocated memory, otherwise
  // allocate it here!!!
  if (!solution_vector.is_own_memory_allocated())
   {
    // Allocate memory for the vector
    solution_vector.allocate_memory(n_values_vector_one);
   }
  else
   {
    // Check that the already allocated memory is correct (n_values)
    if (solution_vector.n_values() != vector_one.n_values())
     {
      // Error message
      std::ostringstream error_message;
      error_message << "The number of elements for the solution vector is\n"
                    << "not as expected\n"
                    << "dim(solution_vector): " << solution_vector.n_values()
                    << "\ndim(vector_one): " << vector_one.n_values()
                    << std::endl;
      throw SciCellxxLibError(error_message.str(),
                             SCICELLXX_CURRENT_FUNCTION,
                             SCICELLXX_EXCEPTION_LOCATION);
     }
    
   }
  
  // Get the vector pointer of the solution vector
  Real *solution_vector_pt = solution_vector.vector_pt();
  
  // Get the vector pointer of the vector one
  Real *vector_one_pt = vector_one.vector_pt();
  // Get the vector pointer of the vector two
  Real *vector_two_pt = vector_two.vector_pt();
  
  // Perform the substraction of vectors
  for (unsigned long i = 0; i < n_values_vector_one; i++)
   {
    solution_vector_pt[i] = vector_one_pt[i] * vector_two_pt[i];
   }
  
 }
 
}
