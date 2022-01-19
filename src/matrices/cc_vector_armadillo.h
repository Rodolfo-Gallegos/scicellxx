// IN THIS FILE: The definition of a concrete class to store and work
// with vectors. This implementation makes use of Armadillo's library,
// thus this is only a wrap for Armadillo's methods

// Check whether the class has been already defined
#ifndef CCVECTORARMADILLO_H
#define CCVECTORARMADILLO_H

// The parent class
#include "ac_vector.h"

// We include the cc_vector include file to deal with transformations
// from CCVector class to CCVectorArmadillo
#include "cc_vector.h"

// Add Armadillo's includes
#include <armadillo>

namespace scicellxx
{
 
 // Concrete class to represent vectors
  class CCVectorArmadillo : public virtual ACVector
  {
   
  public:
   
   // Empty constructor
   CCVectorArmadillo();
   
   // Constructor to create an n size zero vector (we assume vectors
   // are created as column vectors, if you need a row vector then
   // pass "false" as the second parameter)
   CCVectorArmadillo(const unsigned long n, bool is_column_vector = true);
   
   // Constructor where we pass the data for the vector of size n.
   CCVectorArmadillo(Real *vector_pt, const unsigned long n, bool is_column_vector = true);
   
   // Constructor that creates an Armadillo's vector from a CCVector
   CCVectorArmadillo(CCVector &vector);
   
   // Copy constructor (we require to define this if we want to use
   // operators overloading as sum and assignment)
   CCVectorArmadillo(const CCVectorArmadillo &copy);
   
   // Destructor
   virtual ~CCVectorArmadillo();
   
   // Assignment operator
   CCVectorArmadillo& operator=(const CCVectorArmadillo &source_vector);
   
   // += operator
   CCVectorArmadillo& operator+=(const CCVectorArmadillo &vector);
   
   // -= operator
   CCVectorArmadillo& operator-=(const CCVectorArmadillo &vector);
   
   // Add operator
   CCVectorArmadillo operator+(const CCVectorArmadillo &vector);
   
   // Substraction operator
   CCVectorArmadillo operator-(const CCVectorArmadillo &vector);
   
   // Allows to create a vector with the given size but with no data 
   void allocate_memory(const unsigned long n);
   
   // Allocates memory to store entries of the vector
   //void allocate_memory();
   
   // Fills the vector with zeroes
   void fill_with_zeroes();
   
   // Performs dot product with the current vector
   Real dot(const CCVectorArmadillo &right_vector);
   
   // Transforms the input vector to an Armadillo vector class type
   // (virtual such that each derived class has to implement it)
   void set_vector(const Real *vector_pt,
                   const unsigned long n, bool is_column_vector = true);
   
   // Receives an armadillo type Mat
   void set_vector(arma::Mat<Real> *arma_vector_pt,
                   const unsigned long n, bool is_column_vector = true);
   
   // Clean up for any dynamically stored data
   void clean_up();
   
   // Free allocated memory for vector
   void free_memory_for_vector();
   
   // Performs sum of vectors
   void add_vector(const CCVectorArmadillo &vector, CCVectorArmadillo &solution_vector);
   
   // Performs substraction of vectors
   void substract_vector(const CCVectorArmadillo &vector,
                         CCVectorArmadillo &solution_vector);
   
   // Performs multiplication of vectors (one by one entries)
   void multiply_element_by_element_vector(const CCVectorArmadillo &vector,
                                           CCVectorArmadillo &solution_vector);
   
   // Computes the transpose and store it in the transpose vector
   void transpose(CCVectorArmadillo &transposed_vector);
   
   // Transpose the vector
   void transpose();
   
   // Get the specified value from the vector (read-only)
   const Real value(const unsigned long i) const;
   
   // Set values in the vector (write version)
   Real &value(const unsigned long i);
   
   // Output the vector
   void output(bool output_indexes = false) const ;
   
   // Output to file
   void output(std::ofstream &outfile, bool output_indexes = false) const;
   
   // Get access to the Armadillo's vector
   inline arma::Mat<Real> *arma_vector_pt() const {return Arma_vector_pt;}
   
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
   
   // We use an Armadillo type matrix wrapped with specific functions
   // for vectors. We avoided using Col<Real> type because it was not
   // possible to compute the transpose if the same type of object. We
   // did not checked but it may be required to create a Row<Real> type
   arma::Mat<Real> *Arma_vector_pt;
   
  };
  
  // ================================================================
 // Extra methods to work with vectors, we do not need them to be
 // friends of the class since all their operations are performed
 // using the class methods
 // ================================================================
  
 // Dot product of vectors
  Real dot_vectors(const CCVectorArmadillo &left_vector, const CCVectorArmadillo &right_vector);
  
  // Addition of vectors
  void add_vectors(const CCVectorArmadillo &vector_one,
                   const CCVectorArmadillo &vector_two,
                   CCVectorArmadillo &solution_vector);
 
  // Substraction of vectors
  void substract_vectors(const CCVectorArmadillo &vector_one,
                         const CCVectorArmadillo &vector_two,
                         CCVectorArmadillo &solution_vector);
  
  // Performs multiplication of vectors (one by one entries)
  void multiply_element_by_element_vectors(const CCVectorArmadillo &vector_one,
                                           const CCVectorArmadillo &vector_two,
                                           CCVectorArmadillo &solution_vector);
  
  // Concatenate vector horizontally
  void concatenate_vectors_horizontally(const CCVectorArmadillo &left_vector,
                                        const CCVectorArmadillo &right_vector,
                                        CCVectorArmadillo &concatenated_vector);
 
  // Concatenate matrices vertically
  void concatenate_vectors_vertically(const CCVectorArmadillo &upper_vector,
                                      const CCVectorArmadillo &lower_vector,
                                      CCVectorArmadillo &concatenated_vector);
  
}

#endif // #ifndef CCVECTORARMADILLO_H
