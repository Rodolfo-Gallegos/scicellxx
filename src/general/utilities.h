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
#ifndef UTILITIES_H
#define UTILITIES_H

#include "common_includes.h"

namespace scicellxx
{

 // ======================================================================
 /// The error messages are based on oomph-lib's implementation to deal
 /// with errors
 //=======================================================================

 namespace BrokenCopy
 {
  /// Issue error message and terminate execution
  extern void broken_assign(const std::string& class_name);
 
  /// Issue error message and terminate execution
  extern void broken_copy(const std::string& class_name);
 }

 //=======================================================================
 /// Helper namespace for set_terminate function -- used to spawn
 /// messages from uncaught errors (their destructor may not be called)
 ///======================================================================
 namespace TerminateHelper
 {
  /// Setup terminate helper
  extern void setup();

  /// \short Suppress error messages (e.g. because error has been caught)
  extern void suppress_exception_error_messages();

  /// Function to spawn messages from uncaught errors
  extern void spawn_errors_from_uncaught_errors();

  /// Stream to output error messages
  extern std::ostream* Error_message_stream_pt;

  /// String stream that records the error message
  extern std::stringstream* Exception_stringstream_pt;

 }

 //=======================================================================
 /// Helper namespace for file system operations -- mainly used to
 /// create the RESLT folder
 ///======================================================================
 namespace SciCellxxFileSystem
 {
  // Check whether a given directory exists
  extern bool directory_exists(std::string &directory_name);

  // Create a directory
  extern bool create_directory(std::string &directory_name);
 }

 //=======================================================================
 /// Helper namespace to generate the catesian product of multiple
 // /vectors
 ///======================================================================
 namespace SciCellxxCartesianProduct
 {
  // Print the cartesian product of the vectors
  extern void print(const std::vector<std::vector<Real > >& v);
  // Compute the cartesian product of a set of vectors
  extern std::vector<std::vector<Real> > product(const std::vector<std::vector<Real> >& lists);
 }

 //=======================================================================
 /// Helper namespace to a linearspace in the range [min_value,
 /// max_value] with the specified number of points
 ///======================================================================
 namespace SciCellxxLinearSpace
 {
  // Create a linear space with Real values
  extern void create_linear_space(std::vector<Real> &linear_space,
                                  const Real min_value, const Real max_value,
                                  const Real step, const unsigned n_points);
  
  // Create a linear space with integer values
  extern void create_linear_space(std::vector<int> &linear_space,
                                  const int min_value, const int max_value,
                                  const int step, const unsigned n_points);
  
  // Create a linear space with unsigned values
  extern void create_linear_space(std::vector<unsigned> &linear_space,
                                  const unsigned min_value,
                                  const unsigned max_value,
                                  const unsigned, const unsigned n_points);
  
  // Print the linear space (unsigned)
  extern void print_linear_space(std::vector<unsigned> &linear_space);
  // Print the linear space (int)
  extern void print_linear_space(std::vector<int> &linear_space);
  // Print the linear space (Real)
  extern void print_linear_space(std::vector<Real> &linear_space);
  
  
 }
 
 //=====================================================================
 /// Run-time exception handling  (error and warning).
 ///
 /// The (protected) constructor combines its string arguments into a
 /// standard format for uniform exception reports which are written to
 /// the specified output stream.
 ///
 //=====================================================================

 /// The class can only be instantiated by the derived classes
 /// SciCellxxLibError and SciCellxxLibWarning.
 class SciCellxxLibException : public std::runtime_error
 {
 
 public:
 
  /// Suppress error message in destructor (useful if error is caught
  /// successfully!)
  void disable_error_message();
 
 protected:
 
  /// Constructor takes the error description, function name and a
  /// location string provided by the SCICELLXX_EXCEPTION_LOCATION macro
  /// and combines them into a standard header. The exception type will
  /// be the string "WARNING" or "ERROR" and the message is written to
  /// the exception_stream, with a specified output_width. Optionally
  /// provide a traceback of the function calls.
  SciCellxxLibException(const std::string &error_description,
                       const std::string &function_name,
                       const char *location,
                       const std::string &exception_type,
                       std::ostream &exception_stream,
                       const unsigned &output_width);
 
  /// The destructor cannot throw an exception (C++ STL standard)
  ~SciCellxxLibException() throw(); 
 
  /// Exception stream to which we write message in destructor        
  std::ostream* Exception_stream_pt;
 
  /// String stream that records the error message
  std::stringstream* Exception_string_stream_pt;
 
  /// Boolean to suppress issuing of the error message in destructor
  /// (useful if error is caught successfully!)
  bool Suppress_error_message;
 
 };

 /// ====================================================================
 /// Throw this object when an run-time / error is encountered. The
 /// error stream and stream width can be specified. The default is
 /// cerr with a width of 70 characters.
 /// ====================================================================
 class SciCellxxLibError : public SciCellxxLibException
 {
  /// Output stream that is used to write the errors
  static std::ostream *Stream_pt;
 
  /// Width in characters of the output report
  static unsigned Output_width;

 public:

  /// Constructor requires the error description and the function in
  /// which the error occured and the location provided by the
  /// SCICELLXX_EXCEPTION_LOCATION macro
 SciCellxxLibError(const std::string &error_description,
		  const std::string &function_name,
		  const char *location) :
  SciCellxxLibException(error_description,function_name,location,"ERROR",
                       *Stream_pt, Output_width) 
   { }
 
  /// Static member function used to specify the error stream, which
  /// must be passed as a pointer because streams cannot be copied.
  static inline void set_stream_pt(std::ostream* const &stream_pt)
  {Stream_pt = stream_pt;}
 
  /// Static member function used to specify the width (in characters)
  /// of the error stream
  static inline void set_output_width(const unsigned &output_width)
  {Output_width = output_width;}
 
 };

 //====================================================================
 /// An SciCellxxLibWarning object which should be created as a temporary
 /// object to issue a warning. The warning stream and stream width can
 /// be specified. The default is cerr with a width of 70 characters.
 //====================================================================
 class SciCellxxLibWarning : public SciCellxxLibException
 {
  /// Output stream that is used to write the errors
  static std::ostream *Stream_pt;

  /// Width of output
  static unsigned Output_width;

 public:

  /// Constructor requires the warning description and the function
  /// in which the warning occurred.
 SciCellxxLibWarning(const std::string &warning_description,
		    const std::string &function_name,
		    const char* location) :
  SciCellxxLibException(warning_description,function_name, location,
                       "WARNING",*Stream_pt,Output_width) { }
 
  /// Static member function used to specify the error stream, which
  /// must be passed as a pointer because streams cannot be copied.
  static inline void set_stream_pt(std::ostream* const &stream_pt)
  {Stream_pt = stream_pt;}
 
  /// \short Static member function used to specify the width (in characters)
  /// of the error stream
  static inline void set_output_width(const unsigned &output_width)
  {Output_width = output_width;}
 
 };

 ////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////
 
 //========================================================================
 // Wrapper to a stream and an output modifier used to control the
 // output from scicellxx. Its instantiation can be used like std::cout.
 // =======================================================================
 class SciCellxxOutput
 {
  
 private:
  
  ///Pointer to the output stream -- defaults to std::cout
  std::ostream *Stream_pt;
  
 public:
  
  ///\short Set default values for the output stream (cout)
  ///and modifier (no modification)
  SciCellxxOutput();
  
#if 0
  ///\short Overload the << operator, writing output to the stream addressed by
  ///Stream_pt and calling the function defined by the object addressed by
  ///Output_modifier_pt
  template<class _Tp>
   std::ostream &operator<<(_Tp argument);
#endif // #if 0
  
#if 1
  template<class _Tp>
   std::ostream &operator<<(_Tp argument)
   {
    *Stream_pt << argument;
    return (*Stream_pt);
   }
#endif // #if 1
  
  ///Access function for the stream pointer
  std::ostream* &stream_pt() {return Stream_pt;}
 
  ///Overload insertor to handle stream modifiers
  std::ostream &operator<<(std::ostream& (*f)(std::ostream &))
   {
    return f(*Stream_pt);
   }
 
 };
 
 //========================================================================
 // Single (global) instantiation of the SciCellxxOutput object -- this
 // is used throughout the library as a "replacement" for std::cout
 //========================================================================
 extern SciCellxxOutput scicellxx_output;
 
 //==================================================================
 // Utility method to time a program
 //==================================================================
 namespace Timing
 {
  // ================================================================
  // Use this method when you want to know the "wall time" a section
  // of code, method or program takes for its execution. Call it at
  // the beggining of the section and then at the end of the section
  // you want to time
  // ================================================================
  inline time_t wall_time() {return time(0);}
  
  // ================================================================
  // Use this method to get the difference (IN SECONDS) between two
  // "wall_time()" method calls
  // ================================================================
  inline double diff_wall_time(time_t &initial_time, time_t &final_time)
  {return difftime(final_time, initial_time);}
  
  // ================================================================
  // Use this method when you want to know the "cpu time" a section of
  // code, method or program takes for its execution. Call it at the
  // beggining of the section and then at the end of the section you
  // want to time
  // ================================================================
  inline clock_t cpu_clock_time() {return clock();}
  
  // ================================================================
  // Use this method to get the difference (IN SECONDS) between two
  // "cpu_clock_time()" method calls
  // ================================================================
  inline double diff_cpu_clock_time(clock_t &initial_cpu_time,
                                    clock_t &final_cpu_time)
  {return
    static_cast<double>(final_cpu_time-initial_cpu_time)/CLOCKS_PER_SEC;}
  
 }
 
} 

#endif // #ifndef UTILITIES_H
