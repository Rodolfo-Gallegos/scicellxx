
#ifndef ACJACOBIANANDRESIDUAL_H
#define ACJACOBIANANDRESIDUAL_H

#include "../general/general.h"
#include "../data_structures/data_structures.h"
#include "../matrices/matrices.h"

namespace scicellxx
{
 /// An abstract class used as a template for the algorithms to
 /// compute the Jacobian matrix
 class ACJacobianAndResidual
 {
  
 public:
  
  /// Constructor
  ACJacobianAndResidual();
  
  /// Destructor
  virtual ~ACJacobianAndResidual();
  
  /// In charge of computing the Jacobian based on the particular
  /// strategy implemented in the derived class
  virtual void compute_jacobian() = 0;
  
  /// Get access to the Jacobian
  inline ACMatrix *jacobian_pt() {return Jacobian_pt;}
  
  /// In charge of computing the residual vector based on the
  /// particular strategy implemented in the derived class
  virtual void compute_residual() = 0;
   
  /// Get access to the residual
  inline ACVector *residual_pt() {return Residual_pt;}
   
 protected:
  
  /// Storage for the Jacobian matrix
  ACMatrix *Jacobian_pt;
  
  /// Storage for the residual vector
  ACVector *Residual_pt;
  
 private:
   
  /// Copy constructor (we do not want this class to be
  /// copiable. Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  ACJacobianAndResidual(const ACJacobianAndResidual &copy)
   {
    BrokenCopy::broken_copy("ACJacobianAndResidual");
   }
 
  /// Copy constructor (we do not want this class to be
  /// copiable. Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  void operator=(const ACJacobianAndResidual &copy)
   {
    BrokenCopy::broken_assign("ACJacobianAndResidual");
   }
  
 };
 
}

#endif // #ifndef ACJACOBIANANDRESIDUAL_H

