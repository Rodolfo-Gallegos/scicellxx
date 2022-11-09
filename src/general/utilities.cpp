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
#include "utilities.h"

namespace scicellxx
{
 
 // =======================================================================
 /// The error messages are based on oomph-lib's implementation to deal
 /// with errors
 // =======================================================================

 //=======================================================================
 /// The error messages for broken copy constructors and assignment
 /// operators. Check
 /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
 //=======================================================================

 namespace BrokenCopy
 {
  /// Issue error message and terminate execution
  void broken_assign(const std::string& class_name)
  {
   //Write the error message into a string
   std::ostringstream error_message;
   error_message << "No assignment operator allowed for class\n\n";
   error_message << class_name;
   error_message << "\n\n";
   error_message << "Deliberately broken to avoid the accidental \n";
   error_message << "use of the inappropriate C++ default.\n";
   error_message << "If you really need an assignment operator\n";
   error_message << "for this class, write it yourself...\n";
   throw SciCellxxLibError(error_message.str(),
			  SCICELLXX_CURRENT_FUNCTION,
			  SCICELLXX_EXCEPTION_LOCATION);
  }
 
  /// Issue error message and terminate execution
  void broken_copy(const std::string& class_name)
  {
   //Write the error message into a string
   std::ostringstream error_message;
   error_message << "No copy constructor allowed for class\n\n";
   error_message << class_name;
   error_message << "\n\n";
   error_message << "is deliberately broken to avoid the accidental\n";
   error_message << "use of the inappropriate C++ default.\n";
   error_message <<
    "All function arguments should be passed by reference or\n";
   error_message <<
    "constant reference. If you really need a copy constructor\n";
   error_message << "for this class, write it yourself...\n";
   throw SciCellxxLibError(error_message.str(),
			  SCICELLXX_CURRENT_FUNCTION,
			  SCICELLXX_EXCEPTION_LOCATION);
  }

 }

 //=======================================================================
 /// Helper namespace for set_terminate function -- used to spawn
 /// messages from uncaught errors
 ///=======================================================================
 namespace TerminateHelper
 {

  /// Setup terminate helper
  void setup()
  {
   if (Exception_stringstream_pt!=0) delete Exception_stringstream_pt;
   Exception_stringstream_pt=new std::stringstream;
   std::set_terminate(spawn_errors_from_uncaught_errors);
  }
  
  /// \short Flush string stream of error messages (call when error has been
  /// caught)
  void suppress_exception_error_messages()
  {
   delete Exception_stringstream_pt;
   Exception_stringstream_pt=new std::stringstream;
  }

  /// Function to spawn messages from uncaught errors
  void spawn_errors_from_uncaught_errors()
  {
   (*Error_message_stream_pt) << (*Exception_stringstream_pt).str();
  }

  /// Stream to output error messages
  std::ostream* Error_message_stream_pt=&std::cerr;

  /// String stream that records the error message
  std::stringstream* Exception_stringstream_pt=0;
 }

 //=======================================================================
 /// Helper namespace for file system operations -- mainly used to
 /// create the RESLT folder
 ///======================================================================
 namespace SciCellxxFileSystem
 {
  // Check whether a given directory exists
  bool directory_exists(std::string &directory_name)
  { 
   struct stat info;
   
   if( stat( directory_name.c_str(), &info ) != 0 )
    {
     //printf( "cannot access %s\n", directory_name.c_str() );
     return false; // The file or directory does not exist
    }
   else if( info.st_mode & S_IFDIR )  // S_ISDIR() doesn't exist on my windows
    {
     //printf( "%s is a directory\n", directory_name.c_str() );
     return true;
    }
   else
    {
     //printf( "%s is no directory\n", directory_name.c_str() );
     return false; // The name corresponds to a file, not to a directory
    }
   
   return false;

  }

  // Create a directory
  bool create_directory(std::string &directory_name)
  {
   if (directory_exists(directory_name))
    {
     // Error message
     std::ostringstream error_message;
     error_message << "The [" << directory_name << "] folder already exists (or there is a file with the same folder name).\n"
                   << "We will not overwrite the data in that folder.\n"
                   << "Use another folder name.\n"
                   << std::endl;
     throw SciCellxxLibError(error_message.str(),
                             SCICELLXX_CURRENT_FUNCTION,
                             SCICELLXX_EXCEPTION_LOCATION);
    }
   
   // Check whether the directory was created
   int check = mkdir(directory_name.c_str(), 0777);
   
   // check if directory is created or not
   if (!check)
    {
     //printf("Directory created\n");
     return true;
    }
   else
    {
     //printf("Unable to create directory\n");
    }

   // Error message
   std::ostringstream error_message;
   error_message << "The [" << directory_name << "] could not be created\n"
                 << "Check you have the right permissions to work on the filesystem.\n"
                 << std::endl;
   throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   
   return false;
  }
  
 }
 
 //=======================================================================
 /// Helper namespace to generate the catesian product of multiple
 // /vectors
 ///======================================================================
 namespace SciCellxxCartesianProduct
 {
  // Print the cartesian product of the vectors
  void print(const std::vector<std::vector<int> >& v)
  {
   std::cout << "{ ";
   for (const auto& p : v) {
    std::cout << "(";
    for (const auto& e : p) {
     std::cout << e << " ";
    }
    std::cout << ") ";
   }
   std::cout << "}" << std::endl;
  }

  // Compute the cartesian product of a set of vectors
  auto product(const std::vector<std::vector<int> >& lists) -> std::vector<std::vector<int> >
  {
   std::vector<std::vector<int> > result;
   if (std::find_if(std::begin(lists), std::end(lists), 
                    [](auto e) -> bool { return e.size() == 0; }) != std::end(lists))
    {
     return result;
    }
   
   for (auto& e : lists[0])
    {
    result.push_back({ e });
    }
   
   for (size_t i = 1; i < lists.size(); ++i)
    {
     std::vector<std::vector<int>> temp;
     for (auto& e : result)
      {
      for (auto f : lists[i])
       {
        auto e_tmp = e;
        e_tmp.push_back(f);
        temp.push_back(e_tmp);
       }
      }
     result = temp;
    }
   return result;
  }
  
  /*  
  int main() {
   std::vector<std::vector<int>> prods[] = {
    { { 1, 2 }, { 3, 4 } },
    { { 3, 4 }, { 1, 2} },
    { { 1, 2 }, { } },
    { { }, { 1, 2 } },
    { { 1776, 1789 }, { 7, 12 }, { 4, 14, 23 }, { 0, 1 } },
    { { 1, 2, 3 }, { 30 }, { 500, 100 } },
    { { 1, 2, 3 }, { }, { 500, 100 } }
   };
   for (const auto& p : prods) {
    print(product(p));
   }
   std::cin.ignore();
   std::cin.get();
   return 0;
  }
  */
 }
 
 //=======================================================================
 /// Helper namespace to a linearspace in the range [min_value,
 /// max_value] with the specified number of points
 ///======================================================================
 namespace SciCellxxLinearSpace
 {
  // Create a linear space with Real values
  void create_linear_space(std::vector<Real> &linear_space,
                           const Real min_value, const Real max_value,
                           const unsigned n_points)
  {
   const Real step = (max_value - min_value) / n_points;
   Real ival = min_value;
   for (unsigned i = 0; i < n_points; i++)
    {
     linear_space.push_back(ival);
     ival+=step;
    }
   
  }
  
  // Create a linear space with integer values
  void create_linear_space(std::vector<int> &linear_space,
                           const int min_value, const int max_value,
                           const unsigned n_points)
  {
   int ival = min_value;
   for (unsigned i = 0; i < n_points; i++)
    {
     linear_space.push_back(ival);
     ival++;
    }
   
  }
  
  // Create a linear space with unsigned values
  void create_linear_space(std::vector<unsigned> &linear_space,
                           const unsigned min_value, const unsigned max_value,
                           const unsigned n_points)
  {
   unsigned ival = min_value;
   for (unsigned i = 0; i < n_points; i++)
    {
     linear_space.push_back(ival);
     ival++;
    }

  }
  
  // Print the linear space (unsigned)
  void print_linear_space(std::vector<unsigned> &linear_space)
  {
   // Get the number of elements
   const unsigned n_ele = linear_space.size();
   // Print them
   scicellxx_output << "["
   for (unsigned i = 0; i < n_ele-1; i++)
    {
     scicellxx_output << linear_space[i] << ",";
    }
   scicellxx_output linear_space[i] << "]";
   
  }
  
  // Print the linear space (int)
  void print_linear_space(std::vector<int> &linear_space)
  {
   // Get the number of elements
   const unsigned n_ele = linear_space.size();
   // Print them
   scicellxx_output << "["
   for (unsigned i = 0; i < n_ele-1; i++)
    {
     scicellxx_output << linear_space[i] << ",";
    }
   scicellxx_output linear_space[i] << "]";
   
  }
  
  // Print the linear space (Real)
  void print_linear_space(std::vector<Real> &linear_space)
  {
   // Get the number of elements
   const unsigned n_ele = linear_space.size();
   // Print them
   scicellxx_output << "["
   for (unsigned i = 0; i < n_ele-1; i++)
    {
     scicellxx_output << linear_space[i] << ",";
    }
   scicellxx_output linear_space[i] << "]";
   
  }
  
 }
 
 ///////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////

 //==========================================================================
 /// Constructor takes the error description, function name and a
 /// location string provided by the EXCEPTION_LOCATION macro and
 /// combines them into a standard header. The exception type will be
 /// the string "WARNING" or "ERROR" and the message is written to the
 /// exception_stream, with a specified output_width. Optionally
 /// provide a traceback of the function calls.
 //==========================================================================
 SciCellxxLibException::SciCellxxLibException(const std::string &error_description,
                                            const std::string &function_name,
                                            const char *location,
                                            const std::string &exception_type,
                                            std::ostream &exception_stream,
                                            const unsigned &output_width) : 
  std::runtime_error("Exception")
 {
  // By default we shout
  Suppress_error_message=false;
 
  // Store exception stream
  Exception_stream_pt=&exception_stream;
 
  // Create storage for error message
  Exception_string_stream_pt=new std::stringstream;
 
  //Build an exception header string from the information passed
  //Start with a couple of new lines to space things out
  std::string exception_header="\n\n";
 
  //Now add a dividing line
  for(unsigned i=0;i<output_width;i++) {exception_header += "=";}
  exception_header += "\n";
 
  //Write the type of exception
  exception_header += exception_type;

  //Add the function in which it occurs
  exception_header += "\n\n at ";
  exception_header += location;
  exception_header += "\n\n in ";
  exception_header += function_name;

  //Finish with two new lines 
  exception_header +="\n\n";

  //and a closing line
  for(unsigned i=0;i<(unsigned)(0.8*output_width);i++)  
   {exception_header +=  "-";}
 
  //Output the error header to the stream
  (*Exception_string_stream_pt) << exception_header << std::endl;

  //Report the error
  (*Exception_string_stream_pt) << std::endl << error_description << std::endl;
  
  //Finish off with another set of double lines
  for(unsigned i=0;i<output_width;i++) {(*Exception_string_stream_pt) << "=";}
  (*Exception_string_stream_pt) << std::endl << std::endl;
 
  // Copy message to stream in terminate helper in case the message
  // doesn't get caught and/or doesn/t make it to the destructor
  (*TerminateHelper::Exception_stringstream_pt)<<
   (*Exception_string_stream_pt).str();
 
 }

 //========================================================================
 /// The SciCellxxLibException destructor actually spawns the error message
 /// created in the constructor (unless suppresed)
 //==========================================================================
 SciCellxxLibException::~SciCellxxLibException() throw()
 {
  if (!Suppress_error_message)
   {
    (*Exception_stream_pt) << (*Exception_string_stream_pt).str();
   }
  delete Exception_string_stream_pt;
  Exception_string_stream_pt=0;
 }

 //========================================================================
 /// Default output stream for SciCellxxLibErrors (cerr)
 //========================================================================
 std::ostream *SciCellxxLibError::Stream_pt = &std::cerr;

 //=======================================================================
 /// Default output width for SciCellxxLibErrors (70)
 //=======================================================================
 unsigned SciCellxxLibError::Output_width = 70;

 //=======================================================================
 /// Default output stream for SciCellxxLibWarnings (cerr)
 //=======================================================================
 std::ostream *SciCellxxLibWarning::Stream_pt = &std::cerr;

 //=======================================================================
 /// Default output width for SciCellxxLibWarnings (70)
 //=======================================================================
 unsigned SciCellxxLibWarning::Output_width = 70;

 ////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////

 //=======================================================================
 /// Set default values for the output stream (cout) and modifier (no
 /// modification)
 // =======================================================================
 SciCellxxOutput::SciCellxxOutput() : Stream_pt(&std::cout)
 { }

#if 0
 // =======================================================================
 /// Overload the << operator, writing output to the stream addressed by
 /// Stream_pt and calling the function defined by the object addressed
 /// by Output_modifier_pt
 // =======================================================================
 template<class _Tp>
 std::ostream &SciCellxxOutput::operator<<(_Tp argument)
 {
  *Stream_pt << argument;
  return (*Stream_pt);
 }
#endif // #if 0
 
 //========================================================================
 /// Single (global) instantiation of the SciCellxxOutput class -- this
 /// is used throughout the library as a "replacement" for std::cout
 // ========================================================================
 SciCellxxOutput scicellxx_output;
  
}
