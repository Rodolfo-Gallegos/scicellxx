#include "ac_problem.h"

namespace scicellxx
{

 // ===================================================================
 /// Constructor, in charge of initialising any stuff required for the
 /// framework
 // ===================================================================
 ACProblem::ACProblem(const unsigned dim)
  : Output_file_index(0),
    Dim(dim),
    Allow_free_memory_for_U(false)
 { 
  
 }
 
 // ===================================================================
 /// Destructor
 // ===================================================================
 ACProblem::~ACProblem()
 {
  if (Allow_free_memory_for_U)
   {
    // Free memory
    delete U_pt;
    // Set pointer to null
    U_pt = 0;
    // Disabled free memory for U
    Allow_free_memory_for_U = false;
   }
  
 }
 
 // ===================================================================
 /// Document nodes positions
 // ===================================================================
 void ACProblem::document_nodes_positions(std::string &filename)
 {
  // Output supporting nodes
  std::ofstream nodes_file(filename.c_str());
  
  // Get the number of nodes
  const unsigned long nnodes = n_nodes();
  
  // Loop over the nodes and output their position to a file
  for (unsigned long i = 0; i < nnodes; i++)
   {
    // Cache the node
    CCNode *inode_pt = node_pt(i);
    // Output the node index
    nodes_file << i;
    // Loop over the dimensions
    for (unsigned k = 0; k < Dim; k++)
     {
      // Get the node position
      Real position = inode_pt->get_position(k);
      
      // Output nodes positions to file
      nodes_file << "," << position;
     }
    nodes_file << std::endl;
   }
  
  // Close support nodes file
  nodes_file.close();
 }
 
 // ===================================================================
 /// Set/get the i-th node
 // ===================================================================
 CCNode* ACProblem::node_pt(const unsigned long i)
 {
#ifdef SCICELLXX_RANGE_CHECK
  const unsigned long nnodes = n_nodes();
  
  // Check whether the requested node index is in the range [0, n_nodes - 1]
  if (i >= nnodes)
   {
    // Error message
    std::ostringstream error_message;
    error_message << "The node you are trying to access is out of range\n"
                  << "Number of nodes: " << nnodes << std::endl
                  << "Requested node: " << i << std::endl;
    throw SciCellxxLibError(error_message.str(),
                           SCICELLXX_CURRENT_FUNCTION,
                           SCICELLXX_EXCEPTION_LOCATION);    
   }
#endif // #ifdef SCICELLXX_RANGE_CHECK
  
  return Nodes_pt[i];
  
 }

 // ===================================================================
 /// Initialise the u vector (solution)
 // ===================================================================
 void ACProblem::initialise_u(const unsigned n_equations, const unsigned n_history_values)
 {
  U_pt = new CCData(n_equations, n_history_values);
  Allow_free_memory_for_U = true;
 }

 // ===================================================================
 /// Assign equations number
 // ===================================================================
 void ACProblem::assign_equations_number()
 {
  // Loop over the nodes and assign an equation number
  unsigned long index = 0;
  const unsigned long nnodes = n_nodes();
  
  while (index < n_nodes)
   {
    // Get the i-th node
    CNode *inode_pt = node_pt(index);
    
    // Get the number of variables in the node
    const unsigned nvariables = inode_pt->n_variables();

    // Assign an equation number for each variable
    for (unsigned i = 0; i < nvariables; i++)
     {
      // Assign the equation number
      inode_pt->equation_number(i) = index;
      // Increase the equations number
      index++;
     }

   }
  // Create a map that the key is the number of equation, and the
  // value is the pair (node, variable number)

  // When creating the vector U we need to loop over the equations
  // number and get the data in the node and variable number
  
  
 }
 
}
