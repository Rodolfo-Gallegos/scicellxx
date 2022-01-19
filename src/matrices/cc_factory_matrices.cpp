#include "cc_factory_matrices.h"

namespace scicellxx
{

 // ===================================================================
 /// Empty constructor
 // ===================================================================
 CCFactoryMatrices::CCFactoryMatrices()
 { 

 }

 // ===================================================================
 /// Empty destructor
 // ===================================================================
 CCFactoryMatrices::~CCFactoryMatrices()
 { 

 }
 
 // ===================================================================
 /// Returns a matrix pointer (based on compilation options)
 // ===================================================================
 ACMatrix* CCFactoryMatrices::create_matrix()
 {
  // ------------------------------------------------------
  // Check what matrix type we need to create
  // ------------------------------------------------------
#ifdef SCICELLXX_USES_ARMADILLO
  return new CCMatrixArmadillo();
#else
  return new CCMatrix();
#endif // #ifdef SCICELLXX_USES_ARMADILLO
 }

 // ===================================================================
 /// Returns a matrix pointer (based on compilation options)
 // ===================================================================
 ACMatrix* CCFactoryMatrices::create_matrix(const unsigned long m, const unsigned long n)
 {
    // ------------------------------------------------------
  // Check what matrix type we need to create
  // ------------------------------------------------------
#ifdef SCICELLXX_USES_ARMADILLO
  return new CCMatrixArmadillo(m, n);
#else
  return new CCMatrix(m, n);
#endif // #ifdef SCICELLXX_USES_ARMADILLO
 }
 
 // ===================================================================
 /// Returns the specified matrix pointer
 // ===================================================================
 ACMatrix* CCFactoryMatrices::create_matrix(std::string matrix_type_name)
 {
  // Get the string and change it to lower case 
  std::transform(matrix_type_name.begin(), matrix_type_name.end(),
                 matrix_type_name.begin(), ::tolower);
  
  // ------------------------------------------------------
  // Check what matrix type we need to create
  // ------------------------------------------------------
  // Default type
  if (matrix_type_name.compare("default")==0)
   {
    return new CCMatrix();
   }
#ifdef SCICELLXX_USES_ARMADILLO
  // Armadillo type
  else if (matrix_type_name.compare("armadillo")==0)
   {
    return new CCMatrixArmadillo();
   }
#endif // #ifdef SCICELLXX_USES_ARMADILLO
  else
   {
    std::ostringstream error_message;
    error_message << "The matrix time you want to use is not implemented yet.\n"
                  << "Please implement it yourself or select from the available ones\n\n"
                  << "- Default (default)\n"
                  << "- Armadillo matrices (armadillo) - only supported when armadillo library is enabled\n"
                  << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
 }

 // ===================================================================
 /// Returns the specified matrix pointer
 // ===================================================================
 ACMatrix* CCFactoryMatrices::create_matrix(std::string matrix_type_name, const unsigned long m, const unsigned long n)
 {
    // Get the string and change it to lower case 
  std::transform(matrix_type_name.begin(), matrix_type_name.end(),
                 matrix_type_name.begin(), ::tolower);
  
  // ------------------------------------------------------
  // Check what matrix type we need to create
  // ------------------------------------------------------
  // Default type
  if (matrix_type_name.compare("default")==0)
   {
    return new CCMatrix(m, n);
   }
#ifdef SCICELLXX_USES_ARMADILLO
  // Armadillo type
  else if (matrix_type_name.compare("armadillo")==0)
   {
    return new CCMatrixArmadillo(m, n);
   }
#endif // #ifdef SCICELLXX_USES_ARMADILLO
  else
   {
    std::ostringstream error_message;
    error_message << "The matrix time you want to use is not implemented yet.\n"
                  << "Please implement it yourself or select from the availables ones\n\n"
                  << "- Default (default)\n"
                  << "- Armadillo matrices (armadillo) - only supported when armadillo library is enabled\n"
                  << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
 }
 
 // ===================================================================
 /// Returns a vector pointer (based on compilation options)
 // ===================================================================
 ACVector* CCFactoryMatrices::create_vector()
 {
#ifdef SCICELLXX_USES_ARMADILLO
  return new CCVectorArmadillo();
#else
  return new CCVector();
#endif // #ifdef SCICELLXX_USES_ARMADILLO
 }

 // ===================================================================
 /// Returns a vector pointer (based on compilation options)
 // ===================================================================
 ACVector* CCFactoryMatrices::create_vector(const unsigned long n, bool is_column_vector)
 {
#ifdef SCICELLXX_USES_ARMADILLO
  return new CCVectorArmadillo(n, is_column_vector);
#else
  return new CCVector(n, is_column_vector);
#endif // #ifdef SCICELLXX_USES_ARMADILLO
 }
 
 // ===================================================================
 /// Returns the specified vector pointer
 // ===================================================================
 ACVector* CCFactoryMatrices::create_vector(std::string vector_type_name)
 {
  // Get the string and change it to lower case 
  std::transform(vector_type_name.begin(), vector_type_name.end(),
                 vector_type_name.begin(), ::tolower);
  
  // ------------------------------------------------------
  // Check what vector type we need to create
  // ------------------------------------------------------
  // Default type
  if (vector_type_name.compare("default")==0)
   {
    return new CCVector();
   }
#ifdef SCICELLXX_USES_ARMADILLO
  // Armadillo type
  else if (vector_type_name.compare("armadillo")==0)
   {
    return new CCVectorArmadillo();
   }
#endif // #ifdef SCICELLXX_USES_ARMADILLO
  else
   {
    std::ostringstream error_message;
    error_message << "The vector type you want to use is not implemented yet.\n"
                  << "Please implement it yourself or select from the availables ones\n\n"
                  << "- Default (default)\n"
                  << "- Armadillo vector (armadillo) - only supported when armadillo library is enabled\n"
                  << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
 }
 
 // ===================================================================
 /// Returns the specified vector pointer
 // ===================================================================
 ACVector* create_vector(std::string vector_type_name, const unsigned long n, bool is_column_vector)
 {
    // Get the string and change it to lower case 
  std::transform(vector_type_name.begin(), vector_type_name.end(),
                 vector_type_name.begin(), ::tolower);
  
  // ------------------------------------------------------
  // Check what vector type we need to create
  // ------------------------------------------------------
  // Default type
  if (vector_type_name.compare("default")==0)
   {
    return new CCVector(n, is_column_vector);
   }
#ifdef SCICELLXX_USES_ARMADILLO
  // Armadillo type
  else if (vector_type_name.compare("armadillo")==0)
   {
    return new CCVectorArmadillo(n, is_column_vector);
   }
#endif // #ifdef SCICELLXX_USES_ARMADILLO
  else
   {
    std::ostringstream error_message;
    error_message << "The vector type you want to use is not implemented yet.\n"
                  << "Please implement it yourself or select from the availables ones\n\n"
                  << "- Default (default)\n"
                  << "- Armadillo vector (armadillo) - only supported when armadillo library is enabled\n"
                  << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
 }
 
}
