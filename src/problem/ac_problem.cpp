#include "ac_problem.h"

namespace scicellxx
{

 // ===================================================================
 /// Constructor, in charge of initialising any stuff required for the
 /// framework
 // ===================================================================
 ACProblem::ACProblem(const unsigned dim)
  : Output_file_index(0),
    Dim(dim)
 { 
  
 }
 
 // ===================================================================
 /// Destructor
 // ===================================================================
 ACProblem::~ACProblem()
 {
  
 }
 
 // ===================================================================
 /// Document nodes positions
 // ===================================================================
 void ACProblem::document_nodes_positions(const char *filename)
 {
  // Output supporting nodes
  std::ofstream nodes_file(filename);
  
  // Get the number of nodes
  const unsigned long = nnodes = n_nodes();
  
  // Loop over the nodes and output their position to a file
  for (unsigned long i = 0; i < nnodes; i++)
   {
    // Cache the node
    CCNode *node_pt = nodes_pt(i);
    // Output the node index
    nodes_file << i;
    // Loop over the dimensions
    for (unsigned k = 0; k < Dim; k++)
     {
      // Get the node position
      Real position = node_pt->get_position(k);
      
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
 const CCNode* ACProblem::node_pt(const unsigned long i)
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
 
}
