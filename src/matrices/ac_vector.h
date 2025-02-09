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
/// IN THIS FILE: The definition of an abstract class to store and work
/// with vectors. Most of the vectors implemented in this library use
/// this class as the base class

// Check whether the class has been already defined
#ifndef ACVECTOR_H
#define ACVECTOR_H

#include "../general/general.h"

#ifdef SCICELLXX_USES_ARMADILLO
// Add Armadillo's includes (only for the arma_vector methods)
#include <armadillo>
#endif // #ifdef SCICELLXX_USES_ARMADILLO

namespace scicellxx
{
 
 /// @class ACVector ac_vector.h
 
 // Abstract class to represent vector
 class ACVector
 {
   
 public:
   
  // Empty constructor
  ACVector();
   
  // Constructor to create an n size zero vector (we assume vectors
  // are created as column vectors, if you need a row vector then
  // pass "false" as the second parameter)
  ACVector(const unsigned long n, bool is_column_vector = true);
   
  // Destructor
  virtual ~ACVector();
   
  // Allows to create a vector with the given size but with no data
  virtual void allocate_memory(const unsigned long n) = 0;
   
  // Allocates memory to store entries of the vector
  //virtual void allocate_memory() = 0;
   
  // Fills the vector with zeroes
  virtual void fill_with_zeroes() = 0;
   
  // Transforms the input vector to a vector class type (virtual such
  // that each derived class has to implement it)
  virtual void set_vector(const Real *vector_pt,
                          const unsigned long n,
                          bool is_column_vector = true) = 0;
   
  // Clean up for any dynamically stored data
  virtual void clean_up() = 0;
   
  // Free allocated memory for vector
  virtual void free_memory_for_vector() = 0;
   
  // Get the specified value from the vector (read-only)
  virtual const Real value(const unsigned long i) const = 0;
   
  // Set values in the vector (write version)
  virtual Real &value(const unsigned long i) = 0;
   
  // Get the specified value from the vector
  inline Real get_value(const unsigned long i) const
  {return value(i);}
   
  // Set values in the vector
  inline void set_value(const unsigned long i, Real v)
  {value(i) = v;}
   
  /// Get access using brackets as vector(i). Read-only version
  inline virtual Real operator()(const unsigned long &i) const
  {return value(i);}
   
  /// Get access using brackets as vector(i). Read-write version
  inline virtual Real &operator()(const unsigned long &i)
  {return value(i);}
   
  // Output the vector (output horizontally without indexes by
  // default, otherwise output vertically with indexes)
  virtual void output(bool output_indexes = false) const = 0;
  // Output to file (output horizontally without indexes by default,
  // otherwise output vertically with indexes)
  virtual void output(std::ofstream &outfile,
                      bool output_indexes = false) const = 0;
   
  // Output the vector
  inline void print(bool output_indexes = false) const
  {output(output_indexes);}
   
  // Output to file
  inline void print(std::ofstream &outfile,
                    bool output_indexes = false) const
  {output(outfile, output_indexes);}
   
  // Return the number of entries of the vector
  inline unsigned long n_values() const {return NValues;}
   
  // Return the number of entries of the vector
  inline unsigned long n_entries() const {return NValues;}
   
  // Return the number of entries of the vector
  inline unsigned long size() const {return NValues;}
   
  // Check whether the vector should be treated as "row vector" (we
  // assume all vectors are created as column vectors by default)
  inline bool is_column_vector() const {return Is_column_vector;}
   
  // Set as column vector status
  inline void set_as_column_vector(bool status = true)
  {Is_column_vector = status;}
   
  // Checks whether the memory of the vector has been allocated by
  // this class
  inline bool is_own_memory_allocated() const {return Is_own_memory_allocated;}
   
  // Checks whether the vector is allowed to be deleted
  inline bool delete_vector() const {return Delete_vector;}
   
  // Enables the deletion of the vector by itself
  inline void enable_delete_vector() {Delete_vector=true;}
   
  // Disables the deletion of the vector by itself
  inline void disable_delete_vector() {Delete_vector=false;}
   
  // Computes the norm-1 of the vector
  virtual Real norm_1() = 0;
   
  // Computes the norm-2 of the vector
  virtual Real norm_2() = 0;
   
  // Computes the infinite norm
  virtual Real norm_inf() = 0;
   
  // Computes the maximum value
  virtual Real max() = 0;
   
  // Computes the minimum value
  virtual Real min() = 0;
   
  // Get access to the Vector_pt
  virtual Real *vector_pt() const
  {
   // Error message
   std::ostringstream error_message;
   error_message << "Virtual function to resolve vector pointer, should be\n"
                 << "implemented in derived class" << std::endl;
   throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
  }
   
#ifdef SCICELLXX_USES_ARMADILLO
  // Get access to the Armadillo's vector
  virtual arma::Mat<Real> *arma_vector_pt() const
  {
   // Error message
   std::ostringstream error_message;
   error_message << "Virtual function to resolve armadillo vector pointer, should be\n"
                 << "implemented in derived class if you want to use the armadillo solver\n"
                 << std::endl;
   throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
  }
#endif // #ifdef SCICELLXX_USES_ARMADILLO
   
 protected:
   
  // The size of the vector
  unsigned long NValues;
   
  // Flag to indicate whether the memory of the vector has been
  // allocated by this class
  bool Is_own_memory_allocated;
   
  // Flag to indicate whether to delete (free) the allocated memory
  // for the vector. For example when the vector is transformed to an
  // specific vector type (Armadillo vector, SuperLU vector, Trilinos
  // vector) we need to deallocate the memory used for THIS vector to
  // avoid having multiple copies of it. The deletion of the vector
  // is true by default.
  bool Delete_vector;
   
  // Is column vector (we assume all vectors are created as column
  // vectors, thus this flag is set to true by default). We check
  // this variable for all its operations such that they are valid
  // according to the vector dimensions
  bool Is_column_vector;
   
 private:
   
  // Copy constructor (we do not want this class to be copiable). Check
  // http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  ACVector(const ACVector &copy)
   {
    BrokenCopy::broken_copy("ACVector");
   }
   
  // Assignment operator (we do not want this class to be
  // copiable. Check
  // http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  //ACVector& operator=(const ACVector &copy)
  void operator=(const ACVector &copy)
   {
    BrokenCopy::broken_assign("ACVector");
   }
   
 };
 
}

#endif // #ifndef ACVECTOR_H
