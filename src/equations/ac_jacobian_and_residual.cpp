#include "ac_jacobian_and_residual.h"

namespace scicellxx
{
 
 // ===================================================================
 /// Constructor
 // ===================================================================
 ACJacobianAndResidual::ACJacobianAndResidual()
 {
  // Create an instance of the factory for matrices and vectors
#ifdef SCICELLXX_USES_ARMADILLO
  CCFactoryMatrices<CCMatrixArmadillo, CCVectorArmadillo> factory_matrices_and_vectors;
#else
  CCFactoryMatrices<CCMatrix, CCVector> factory_matrices_and_vectors;
#endif // #ifdef SCICELLXX_USES_ARMADILLO
  
  /// Instantiate Jacobian and Residual based on the type of libraries
  /// we are using.
  Jacobian_pt = factory_matrices_and_vectors.create_matrix();
  Residual_pt = factory_matrices_and_vectors.create_vector();
  
 }

 // ===================================================================
 /// Destructor
 // ===================================================================
 ACJacobianAndResidual::~ACJacobianAndResidual()
 {
  delete Jacobian_pt;
  Jacobian_pt = NULL;
  
  delete Residual_pt;
  Residual_pt = NULL;
  
 }
 
}

