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
/// IN THIS FILE: The definition of a class to represent boundary
/// nodes

// Check whether the class has been already defined
#ifndef CCBOUNDARYNODE_H
#define CCBOUNDARYNODE_H

#include "./cc_node.h"

namespace scicellxx
{
 
 /// @class CCBoundaryNode cc_boundary_node.h
 
 /// Abstract class to represent boundary nodes
 class CCBoundaryNode : public virtual CCNode
 {
  
 public:
  
  /// Constructor
  CCBoundaryNode(const unsigned boundary,
                 const unsigned dimension, const unsigned n_variables, const unsigned n_history_values=1);
  
  /// Empty destructor
  virtual ~CCBoundaryNode();
  
  /// Add the node to the spefified boundary
  void add_to_boundary(const unsigned b);
  
  /// Remove the node from an specific boundary
  void remove_from_boundary(const unsigned b);
  
  /// Test whether the node is on any boundary
  bool is_on_boundary();
  
  /// Test whether the node is on a given boundary
  bool is_on_boundary(const unsigned b);
  
  /// Set the boundary coordinates zeta of the node on boundary b
  void set_boundary_coordinates(const unsigned b, const std::vector<Real> &zeta);
  
  /// Get the boundary coordinates zeta of the node on boundary b
  void get_boundary_coordinates(const unsigned b, std::vector<Real> &zeta);
  
  /// Get the boundaries the node is on (read only)
  inline const std::set<unsigned> &boundaries() const {return Boundaries;}
  
  /// Output the data stored at the node (boundary ids, position and
  /// values at time t)
  void output_boundary_position_and_value(const unsigned t = 0);
  
  /// Output the data stored at the node (boundary ids, position and
  /// values at time t)
  void output_boundary_position_and_value(std::ofstream &outfile, const unsigned t = 0) const;
  
 protected:
  
  // Copy constructor (we do not want this class to be copiable because
  // it contains dynamically allocated variables, A in this
  // case). Check
  // http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  CCBoundaryNode(const CCBoundaryNode& boundary_node)
   : CCNode(boundary_node.dimension(), boundary_node.n_variables(),
            boundary_node.n_history_values())
   {
    BrokenCopy::broken_copy("CCBoundaryNode");
   }
  
  // Copy constructor (we do not want this class to be copiable because
  // it contains dynamically allocated variables, A in this
  // case). Check
  // http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  void operator=(const CCBoundaryNode&) 
   {
    BrokenCopy::broken_assign("CCBoundaryNode");
   }
  
  /// A set to store the boundaries where the node lies on
  std::set<unsigned> Boundaries;
  
  /// Intrinsic boundary coordiantes for a given boundary
  std::map<unsigned, std::vector<Real> > Boundary_coordinates;
  
 };
 
}

#endif // #ifndef CCBOUNDARYNODE_H
