#ifndef ACPROBLEM_H
#define ACPROBLEM_H

#include "../general/general.h"

#include "../time_steppers/time_steppers.h"

namespace scicellxx
{

 /// @class ACProblem ac_problem.h

 /// This class implements the interfaces for the problem class. It
 /// specifies a template for solving a problem thus one needs to
 /// create a class that inherents from this one to solve a particular
 /// problem
 class ACProblem
 {
  
 public:
  
  /// Constructor
  ACProblem(const unsigned dim = 1);
  
  /// Destructor
  virtual ~ACProblem();

  /// Complete the problem setup (things you could not do in the
  /// construtor as initial conditions/boundary conditions/solver/etc)
  virtual void complete_problem_setup() = 0;
  
  /// Every derived class must implement its own solve method (calling
  /// the corresponding steady_solve() and unsteady_solve() methods)
  virtual void solve() = 0;

  /// Document solution
  virtual void document_solution() = 0;
  
  /// Write access to the current time step
  inline unsigned &output_file_index() {return Output_file_index;}
  
  /// Read-only access to the current time step
  inline unsigned output_file_index() const {return Output_file_index;}
  
  /// Document nodes positions
  void document_nodes_positions(std::string &filename);
  
  /// The dimension of the problem
  inline unsigned dim() const {return Dim;}
  
  /// Get the number of nodes
  inline unsigned long n_nodes() const {return Nodes_pt.size();}
  
  /// Set/get the i-th node
  CCNode* node_pt(const unsigned long i);

  /// Get access to the U vector
  CCData *u_pt() const {return U_pt;}
  
  /// Read-only access to the vector U values
  inline const Real u(const unsigned i, const unsigned t = 0) const {return U_pt->value(i,t);}
  
  /// Write access to the vector U values
  inline Real &u(const unsigned i, const unsigned t = 0) {return U_pt->value(i,t);}
  
  /// Initialise the u vector (solution)
  void initialise_u(const unsigned n_equations, const unsigned n_history_values);
  
  /// Assign equations number
  void assign_equations_number();
  
 protected:
  
  /// Copy constructor (we do not want this class to be
  /// copiable). Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
 ACProblem(const ACProblem &copy)
  : Output_file_index(0),
   Dim(0)
    {
     BrokenCopy::broken_copy("ACProblem");
    }
  
  /// Assignment operator (we do not want this class to be
  /// copiable. Check
  /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
  void operator=(const ACProblem &copy)
   {
    BrokenCopy::broken_assign("ACProblem");
   }
  
  /// Problem steady solve
  virtual void steady_solve() = 0;
  
  /// Problem unsteady solve
  virtual void unsteady_solve() = 0;
  
  /// The set of actions to be performed before newton solve
  virtual void actions_before_newton_solve() = 0;
  
  /// The set of actions to be performed after newton solve
  virtual void actions_after_newton_solve() = 0;
  
  /// A counter to store the current output file index
  unsigned Output_file_index;
    
  /// Dimension
  const unsigned Dim;
  
  // Total number of nodes
  unsigned N_nodes;
  
  // The nodes
  std::vector<CCNode *> Nodes_pt;
  
  /// The storage for the computed solution
  CCData *U_pt;
  
  /// Flag to allow release of memory by the class
  bool Allow_free_memory_for_U;
  
  /// Store the number of equations of the problem
  unsigned N_equations;
  
 };
 
}

#endif // #ifndef ACPROBLEM_H
