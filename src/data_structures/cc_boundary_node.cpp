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
/// IN THIS FILE: Implementation of a concrete class to represent
/// boundary nodes. This is the simplest implementation

#include "cc_boundary_node.h"

namespace scicellxx
{

 /// ===================================================================
 /// Constructor
 /// ===================================================================
 CCBoundaryNode::CCBoundaryNode(const unsigned boundary,
                                   const unsigned dimension,
                                   const unsigned n_variables,
                                   const unsigned n_history_values)
  : CCNode(dimension, n_variables, n_history_values)
 {
  // Add the node to the boundary
  add_to_boundary(boundary);
 }
 
 /// ===================================================================
 /// Empty destructor
 /// ===================================================================
 CCBoundaryNode::~CCBoundaryNode()
 {
  
 }
 
 /// ===================================================================
 /// Add the node to the spefified boundary
 /// ===================================================================
 void CCBoundaryNode::add_to_boundary(const unsigned b)
 {
  Boundaries.insert(b);
 }
 
 /// ===================================================================
 /// Remove the node from an specific boundary
 /// ===================================================================
 void CCBoundaryNode::remove_from_boundary(const unsigned b)
 {
#ifdef SCICELLXX_PANIC_MODE
  // First check whether the node lies on the given boundary
  if(!is_on_boundary(b))
   {
    std::ostringstream error_message;
    error_message << "The node is not on boundary " << b << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
#endif // #ifdef SCICELLXX_PANIC_MODE
  
  //Remove the boundary from the set 
  Boundaries.erase(b);
  
  //Need to delete the equivalent entry in the Boundary coordinate
  //map, if the storage has actually been allocated
  std::map<unsigned, std::vector<Real> >::iterator it = Boundary_coordinates.find(b);
  if (it != Boundary_coordinates.end())
   {
    Boundary_coordinates.erase(it);
   }
  else
   {
    std::ostringstream error_message;
    error_message << "The node is on boundary " << b << " but has no"
                  << "boundary coordinates assigned\n" << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
 }
 
 /// ===================================================================
 /// Test whether the node is on any boundary
 /// ===================================================================
 bool CCBoundaryNode::is_on_boundary()
 {
  if (Boundaries.size() != 0)
   {
    return true;
   }
  
  return false;
  
 }
 
 /// ===================================================================
 /// Test whether the node is on a given boundary
 /// ===================================================================
 bool CCBoundaryNode::is_on_boundary(const unsigned b)
 {
  std::set<unsigned>::iterator it = Boundaries.find(b);
  if (it != Boundaries.end())
   {
    return true;
   }
  
  return false;
  
 }
 
 /// ===================================================================
 /// Set the boundary coordinates zeta of the node on boundary b
 /// ===================================================================
 void CCBoundaryNode::set_boundary_coordinates(const unsigned b, const std::vector<Real> &zeta)
 {
#ifdef SCICELLXX_PANIC_MODE
  // First check whether the node lies on the given boundary
  if(!is_on_boundary(b))
   {
    std::ostringstream error_message;
    error_message << "The node is not on boundary " << b << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
#endif // #ifdef SCICELLXX_PANIC_MODE
  
  // Get the number of boundary coordinates and copy them
  const unsigned n_boundary_coordinates = zeta.size();
  Boundary_coordinates[b].resize(n_boundary_coordinates);
  for (unsigned i = 0; i < n_boundary_coordinates; i++)
   {
    Boundary_coordinates[b][i] = zeta[i];
   }
  
 }
 
 /// ===================================================================
 /// Get the boundary coordinates zeta of the node on boundary b
 /// ===================================================================
 void CCBoundaryNode::get_boundary_coordinates(const unsigned b, std::vector<Real> &zeta)
 {
#ifdef SCICELLXX_PANIC_MODE
  // First check whether the node lies on the given boundary
  if(!is_on_boundary(b))
   {
    std::ostringstream error_message;
    error_message << "The node is not on boundary " << b << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
  
  // Check whether there are boundary coordinates for this boundary
  std::map<unsigned, std::vector<Real> >::iterator it = Boundary_coordinates.find(b);
  if (it == Boundary_coordinates.end())
   {
    std::ostringstream error_message;
    error_message << "The node is on boundary " << b << " but has no"
                  << "boundary coordinates assigned\n" << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
#endif // #ifdef SCICELLXX_PANIC_MODE
  
  // Check whether the number of boundary coordinates matches that of
  // the output vector
  const unsigned n_boundary_coordinates = Boundary_coordinates[b].size();
#ifdef SCICELLXX_PANIC_MODE
  const unsigned n_output_boundary_coordinates = zeta.size();
  if (n_boundary_coordinates != n_output_boundary_coordinates)
   {
    std::ostringstream error_message;
    error_message << "The number of boundary coordinates of the output vector is\n"
                  << "not the same as the number of boundary coordinates for\n"
                  << "boundary " << b << ".\n" << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);
   }
#endif // #ifdef SCICELLXX_PANIC_MODE
  
  // Copy the boundary coordinates
  for (unsigned i = 0; i < n_boundary_coordinates; i++)
   {
    zeta[i] = Boundary_coordinates[b][i];
   }
  
 }
 
 /// ===================================================================
 /// Output the data stored at the node (boundary ids, position and
 /// values at time t)
 /// ===================================================================
 void CCBoundaryNode::output_boundary_position_and_value(const unsigned t)
 {
  std::set<unsigned>::iterator it;
  // Output boundaries
  std::cout << "B: ";
  for (it = Boundaries.begin(); it != Boundaries.end(); it++)
   {
    std::cout << (*it) << "\t";
   }
  
  // Output positions
  std::cout << "X: ";
  for (unsigned i = 0; i < this->Dimension; i++)
   {
    std::cout << this->X.value(i,t) << "\t";
   } // for (i < this->Dimension)

  // Output data
  std::cout << "U: ";
  for (unsigned i = 0; i < this->N_variables; i++)
   {
    std::cout << this->U.value(i,t) << "\t";
   } // for (i < this->N_variables)
  
  std::cout << std::endl;
  
 }
 
 /// ===================================================================
 /// Output the data stored at the node (boundary ids, position and
 /// values at time t)
 /// ===================================================================
 void CCBoundaryNode::output_boundary_position_and_value(std::ofstream &outfile, const unsigned t) const
 {
  // Output boundaries
  std::cout << "B: ";
  std::set<unsigned>::iterator it;
  for (it = Boundaries.begin(); it != Boundaries.end(); it++)
   {
    outfile << (*it) << "\t";
   }
  
  // Output positions
  std::cout << "X: ";
  for (unsigned i = 0; i < this->Dimension; i++)
   {
    outfile << this->X.value(i,t) << "\t";
   } // for (i < this->Dimension)
  
  // Output data
  std::cout << "U: ";
  for (unsigned i = 0; i < this->N_variables; i++)
   {
    outfile << this->U.value(i,t) << "\t";
   } // for (i < this->N_variables)
  
  outfile << std::endl;  
  
 }
 
}
