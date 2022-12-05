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
#include "cc_jacobian_by_fd_from_odes_with_shampine_h_computation.h"

namespace scicellxx
{
 // ===================================================================
 /// Constructor
 // ===================================================================
 CCJacobianByFDFromODEs::CCJacobianByFDFromODEs()
  : Set_of_ODEs_have_been_set(false)
 {
 
 }
 
 // ===================================================================
 /// Destructor
 // ===================================================================
 CCJacobianByFDFromODEs::~CCJacobianByFDFromODEs()
 {
 
 }

 // ===================================================================
 /// In charge of computing the Jacobian using Finite Differences
 // ===================================================================
 void CCJacobianByFDFromODEs::compute_jacobian(const double time)
 {
  // Check whether the ODEs have been set
  if (!Set_of_ODEs_have_been_set)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "You have not established the ODEs that should be\n"
                  << "used to compute the Jacobian matrix\n."
                  << "Call the method set_ODEs() to specify the set of ODEs\n"
                  << "used to compute the Jacobian matrix." << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  int i, j, nEq;
  
  const double raiz_uno_entre_cinco = 4.47213595499957e-01;
  const double epsilon = 8.881784e-16; // >> 2.220446 * 4
  
  const double hmin = 1.0e-12;
  const double hmax = 1.0e-01;
  
  double h;
  
  double hmaxesc;
  //double hminesc;
  
  nEq = odes->getnODEs();
  
  double *fy = new double[nEq];
  double *fyh = new double[nEq];
  double *yh = new double[nEq];
  
  double absy;
  
  //memcpy(fy, _y, sizeof(double) * nEq);
  
  // >> Evaluate and save in fy
  //for (i = idxI; i <= idxF; i++)
  for (i = 0; i < nEq; i++)
   (*odes)[i](_t, _y, fy);
  
  //modelo->evalua(_t, _y, fy);
  
  for (j = 0; j < nEq; j++) {
   
   absy = fabs(_y[j]);
   
   // >> Get the maximum scale h
   hmaxesc = hmax * max(1.0, absy);
   
   // >> Get the h
   if (absy <= 1.0e-10)
    h = hmin;
   else
    h = raiz_uno_entre_cinco*sqrt(epsilon*absy);
   
   // >> Verifying low h quote
   if (h < hmin)
    h = hmin;
   
   // >> Scaling h value
   h = h * max(1.0, absy);
   
   // >> Verifying up h quote
   if (h > hmaxesc)
    h = hmaxesc;
   
   // >> Major estimation of h according with Shampine
   h = (_y[j] + h) - _y[j];
   
   // >> Forward h
   memcpy(yh, _y, sizeof(double) * nEq);
   yh[j]+= h;
   
   // >> Evaluate for (y + h)
   for (i = 0; i < nEq; i++)
    (*odes)[i](_t, yh, fyh);
   // modelo->evalua(_t, yh, fyh);
   
   for (i = 0; i < nEq; i++)
    _J[i][j] = (fyh[i] - fy[i]) / h;
   
  }
  
  delete [] fy;
  delete [] fyh;
  delete [] yh;
  
  fy = fyh = yh = NULL;
  
 }

 // ===================================================================
 // Set the ODEs to compute the Jacobian
 // ===================================================================
 void CCJacobianByFDFromODEs::set_ODEs(const ACODEs &odes)
 {
  // Set a pointer to the odes
  ODEs = odes;

  // Indicate that the ODEs have been set
  Set_of_ODEs_have_been_set = true;
  
 }

}
