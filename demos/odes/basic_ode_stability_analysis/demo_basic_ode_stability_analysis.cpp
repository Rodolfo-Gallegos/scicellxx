// Include SciCell++ libraries
#include "../../../src/scicellxx.h"

using namespace scicellxx;
// =================================================================
// =================================================================
// =================================================================
/// This class implements the system of ODEs to be solved
///
/// \frac{du}{dt} = -u^{2}, with initial values u(0) = 1
// =================================================================
// =================================================================
// =================================================================
class CCBasicODEs : public virtual ACODEs
{
 
public:
 
 /// Constructor
 CCBasicODEs()
  :ACODEs(1) // The number of equations
 {
  
 }
 
 /// Empty destructor
 virtual ~CCBasicODEs()
 {
  
 }
 
 /// Evaluates the system of odes at time 't', using the history
 /// values of u at index k
 void evaluate_time_derivatives(const Real t, CCData &u, CCData &dudt, const unsigned k = 0)
 {
  // \frac{du}{dt} = -u^{2}
  dudt(0) = -(u(0,k)*u(0,k));
 }
 
protected:
 
 /// Copy constructor (we do not want this class to be
 /// copiable). Check
 /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
 CCBasicODEs(const CCBasicODEs &copy)
  : ACODEs(copy)
 {
  BrokenCopy::broken_copy("CCBasicODEs");
 }
 
 /// Assignment operator (we do not want this class to be
 /// copiable. Check
 /// http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/
 void operator=(const CCBasicODEs &copy)
 {
  BrokenCopy::broken_assign("CCBasicODEs");
 }
 
};

// =================================================================
// =================================================================
// =================================================================
// This class inherits from the ACIBVP class and solves the system of
// ODEs from above
// =================================================================
// =================================================================
// =================================================================
template<class EQUATIONS_TYPE>
class CCStabilityAnalysisProblem : public virtual ACIBVP<EQUATIONS_TYPE>
{
 
public:
 
 /// Constructor
 CCStabilityAnalysisProblem(EQUATIONS_TYPE *odes_pt,
                            ACTimeStepper<EQUATIONS_TYPE> *time_stepper_pt,
                            std::ostringstream &output_filename_prefix,
                            std::ostringstream &output_filename_error_prefix,
                            std::ostringstream &output_filename_stability)
  : ACIBVP<EQUATIONS_TYPE>(odes_pt, time_stepper_pt),
    Output_filename_prefix(output_filename_prefix.str()),
    Output_error_filename_prefix(output_filename_error_prefix.str()),
    Output_stability_filename(output_filename_stability.str())
 {
  std::ostringstream output_stability_filename;
  output_stability_filename << Output_stability_filename.str() << ".dat";
  Output_stability_file.open((output_stability_filename.str()).c_str());
 }
 
 /// Destructor
 ~CCStabilityAnalysisProblem()
 {
  // Close open files
  Output_file.close();
  Output_error_file.close();
  Output_stability_file.close();
 }

 // Complete problem configuration
 void complete_problem_setup()
 {
   // Prepare time integration data
   Initial_time = 0.0;
   Final_time = 20.0;
   
   // Open files
   prepare_files_for_output();
   
   // Reset problem current time
   this->time() = 0.0;
   
   // Reset initial conditions
   set_initial_conditions();
   
 }
 
 // Prepare the filenames for new output
 void prepare_files_for_output()
 {
  // Close previouly open files
  Output_file.close();
  Output_error_file.close();
  
  // Generate the new files names
  std::ostringstream output_filename;
  std::ostringstream output_filename_error;
  
  output_filename << Output_filename_prefix.str() << "_h_" << this->time_step() << ".dat";
  output_filename_error << Output_error_filename_prefix.str() << "_h_" << this->time_step() << ".dat";
  
  Output_file.open((output_filename.str()).c_str());
  Output_error_file.open((output_filename_error.str()).c_str());
 }
 
 // Solve the problem with a new time step
 void solve()
 {  
  // Document initial solution
  this->document_solution();
     
  // Flag to indicate whether to continue processing
  bool LOOP = true;
   
  // Main LOOP (loop until reaching final time)
  while(LOOP)
   {
    // Solve (unsteady solve) - PARENT VERSION
    ACIBVP<EQUATIONS_TYPE>::solve();
    
    // Update time of the problem
    this->time()+=this->time_step();
    
    // Check whether we have reached the final time
    if (this->time() >= Final_time)
     {
      LOOP = false;
     }
     
    this->document_solution();
     
   } // while(LOOP)
  
  // Document error at $u = Final_time$
  const Real u_analytical = 1.0/(1.0+Final_time);
  Output_stability_file << this->time_step() << "\t" << fabs(this->u(0)- u_analytical) << std::endl;
  
 }
 
 // Set initial conditions
 void set_initial_conditions()
 {
  // Initial conditions
  this->u(0) = 1.0; 
 }
 
 // Document the solution
 void document_solution()
 {
  // Initial problem configuration
  Output_file << this->time() << "\t" << this->u(0) << std::endl;
  output_error();
 }

 // Output error
 void output_error()
 {
  // Compute the error 
  const Real t = this->time();
  const Real u_analytical = 1.0/(1.0+t);
  const Real error = std::fabs(this->u(0)-u_analytical);
  Output_error_file << t << "\t" << error << std::endl;
 }
 
protected:
 
 // The output file prefix
 std::ostringstream Output_filename_prefix;
 // The error output file prefix
 std::ostringstream Output_error_filename_prefix; 
 // The stability output file prefix
 std::ostringstream Output_stability_filename;
 
 // The initial time
 Real Initial_time; 
 // The final time
 Real Final_time;
  
 // The output file
 std::ofstream Output_file;
 // The error output file
 std::ofstream Output_error_file;
 // The stability output file
 std::ofstream Output_stability_file;
  
}; // class CCStabilityAnalysisProblem

// ==================================================================
// ==================================================================
// ==================================================================
// Main function
// ==================================================================
// ==================================================================
// ==================================================================
int main(int argc, char *argv[])
{
 // Initialise scicellxx
 initialise_scicellxx();
 
 // Create the factory for the time steppers (integration methods)
 CCFactoryTimeStepper<CCBasicODEs> factory_time_stepper;
 
 // Euler stability analysis
 {
  std::cout << "--------------------------" << std::endl;
  std::cout << "Euler - Stability analysis" << std::endl;
  // -----------------------------------------------------------------
  // Instantiation of the ODEs
  // -----------------------------------------------------------------
  CCBasicODEs odes;
  
  // ----------------------------------------------------------------
  // Time stepper
  // ----------------------------------------------------------------
  // Create an instance of the integration method
  ACTimeStepper<CCBasicODEs> *time_stepper_pt =
   factory_time_stepper.create_time_stepper("Euler");
  
  // ----------------------------------------------------------------
  // Prepare the output file name prefix
  // ----------------------------------------------------------------
  std::ostringstream output_filename_prefix;
  output_filename_prefix << "RESLT/euler";
  
  // ----------------------------------------------------------------
  // Prepare the output error file name
  // ----------------------------------------------------------------
  std::ostringstream output_error_filename_prefix;
  output_error_filename_prefix << "RESLT/euler_error";
  
  // ----------------------------------------------------------------
  // Prepare the output stability file name
  // ----------------------------------------------------------------
  std::ostringstream output_stability_filename_prefix;
  output_stability_filename_prefix << "RESLT/euler_stability";
  
  // Create an instance of the problem
  CCStabilityAnalysisProblem<CCBasicODEs> stability_analysis_problem(&odes,
                                                                     time_stepper_pt,
                                                                     output_filename_prefix,
                                                                     output_error_filename_prefix,
                                                                     output_stability_filename_prefix);
  
  // ----------------------------------------------------------------
  // Configure the problem with different time steps to perform
  // stability analisys
  // ----------------------------------------------------------------
  Real time_step = 1.0;
  stability_analysis_problem.time_step() = time_step;
  // Complete problem setup
  stability_analysis_problem.complete_problem_setup();
  // Solve
  stability_analysis_problem.solve();
  
  time_step = 0.1;
  stability_analysis_problem.time_step() = time_step;
  // Complete problem setup
  stability_analysis_problem.complete_problem_setup();
  // Solve
  stability_analysis_problem.solve();
  
  time_step = 0.01;
  stability_analysis_problem.time_step() = time_step;
  // Complete problem setup
  stability_analysis_problem.complete_problem_setup();
  // Solve
  stability_analysis_problem.solve();
  
  time_step = 0.001;
  stability_analysis_problem.time_step() = time_step;
  // Complete problem setup
  stability_analysis_problem.complete_problem_setup();
  // Solve
  stability_analysis_problem.solve();
    
  std::cout << "[FINISHING UP] ... " << std::endl;
  
  // Free memory
  delete time_stepper_pt;
  time_stepper_pt = 0;
  
 }
 
 {
  std::cout << "--------------------------" << std::endl;
  std::cout << "Runge-Kutta 4 - Stability analysis" << std::endl;
  // -----------------------------------------------------------------
  // Instantiation of the ODEs
  // -----------------------------------------------------------------
  CCBasicODEs odes;
  
  // ----------------------------------------------------------------
  // Time stepper
  // ----------------------------------------------------------------
  ACTimeStepper<CCBasicODEs> *time_stepper_pt =
   factory_time_stepper.create_time_stepper("RK4");
  
  // ----------------------------------------------------------------
  // Prepare the output file name prefix
  // ----------------------------------------------------------------
  std::ostringstream output_filename_prefix;
  output_filename_prefix << "RESLT/rk4";
  
  // ----------------------------------------------------------------
  // Prepare the output error file name
  // ----------------------------------------------------------------
  std::ostringstream output_error_filename_prefix;
  output_error_filename_prefix << "RESLT/rk4_error";
  
  // ----------------------------------------------------------------
  // Prepare the output stability file name
  // ----------------------------------------------------------------
  std::ostringstream output_stability_filename_prefix;
  output_stability_filename_prefix << "RESLT/rk4_stability";
  
  // Create an instance of the problem
  CCStabilityAnalysisProblem<CCBasicODEs> stability_analysis_problem(&odes,
                                                                     time_stepper_pt,
                                                                     output_filename_prefix,
                                                                     output_error_filename_prefix,
                                                                     output_stability_filename_prefix);
  
  // ----------------------------------------------------------------
  // Configure the problem with different time steps to perform
  // stability analisys
  // ----------------------------------------------------------------
  Real time_step = 1.0;
  stability_analysis_problem.time_step() = time_step;
  // Complete problem setup
  stability_analysis_problem.complete_problem_setup();
  // Solve
  stability_analysis_problem.solve();
  
  time_step = 0.1;
  stability_analysis_problem.time_step() = time_step;
  // Complete problem setup
  stability_analysis_problem.complete_problem_setup();
  // Solve
  stability_analysis_problem.solve();
  
  time_step = 0.01;
  stability_analysis_problem.time_step() = time_step;
  // Complete problem setup
  stability_analysis_problem.complete_problem_setup();
  // Solve
  stability_analysis_problem.solve();
  
  time_step = 0.001;
  stability_analysis_problem.time_step() = time_step;
  // Complete problem setup
  stability_analysis_problem.complete_problem_setup();
  // Solve
  stability_analysis_problem.solve();
      
  std::cout << "[FINISHING UP] ... " << std::endl;
  
  // Free memory
  delete time_stepper_pt;
  time_stepper_pt = 0;
  
 }
 
 {
  std::cout << "--------------------------" << std::endl;
  std::cout << "Backward Euler - Predictor-Corrector - Stability analysis" << std::endl;
  // -----------------------------------------------------------------
  // Instantiation of the ODEs
  // -----------------------------------------------------------------
  CCBasicODEs odes;
  
  // ----------------------------------------------------------------
  // Time stepper
  // ----------------------------------------------------------------
  ACTimeStepper<CCBasicODEs> *time_stepper_pt =
   factory_time_stepper.create_time_stepper("BEPC");
  
  // ----------------------------------------------------------------
  // Prepare the output file name prefix
  // ----------------------------------------------------------------
  std::ostringstream output_filename_prefix;
  output_filename_prefix << "RESLT/bepc";
  
  // ----------------------------------------------------------------
  // Prepare the output error file name
  // ----------------------------------------------------------------
  std::ostringstream output_error_filename_prefix;
  output_error_filename_prefix << "RESLT/bepc_error"; 
  
  // ----------------------------------------------------------------
  // Prepare the output stability file name
  // ----------------------------------------------------------------
  std::ostringstream output_stability_filename_prefix;
  output_stability_filename_prefix << "RESLT/bepc_stability";
  
  // Create an instance of the problem
  CCStabilityAnalysisProblem<CCBasicODEs> stability_analysis_problem(&odes,
                                                                     time_stepper_pt,
                                                                     output_filename_prefix,
                                                                     output_error_filename_prefix,
                                                                     output_stability_filename_prefix);
  
  // ----------------------------------------------------------------
  // Configure the problem with different time steps to perform
  // stability analisys
  // ----------------------------------------------------------------
  Real time_step = 1.0;
  stability_analysis_problem.time_step() = time_step;
  // Complete problem setup
  stability_analysis_problem.complete_problem_setup();
  // Solve
  stability_analysis_problem.solve();

  time_step = 0.1;
  stability_analysis_problem.time_step() = time_step;
  // Complete problem setup
  stability_analysis_problem.complete_problem_setup();
  // Solve
  stability_analysis_problem.solve();
  
  time_step = 0.01;
  stability_analysis_problem.time_step() = time_step;
  // Complete problem setup
  stability_analysis_problem.complete_problem_setup();
  // Solve
  stability_analysis_problem.solve();
  
  time_step = 0.001;
  stability_analysis_problem.time_step() = time_step;
  // Complete problem setup
  stability_analysis_problem.complete_problem_setup();
  // Solve
  stability_analysis_problem.solve();
  
  std::cout << "[FINISHING UP] ... " << std::endl;
  
  // Free memory
  delete time_stepper_pt;
  time_stepper_pt = 0;
  
 }
 
 {
  std::cout << "--------------------------" << std::endl;
  std::cout << "Adams-Moulton 2 or Trapezoidal Rule - Predictor-Corrector - Stability analysis" << std::endl;
  // -----------------------------------------------------------------
  // Instantiation of the ODEs
  // -----------------------------------------------------------------
  CCBasicODEs odes;
  
  // ----------------------------------------------------------------
  // Time stepper
  // ----------------------------------------------------------------
  ACTimeStepper<CCBasicODEs> *time_stepper_pt =
   factory_time_stepper.create_time_stepper("AM2PC");

  // ----------------------------------------------------------------
  // Prepare the output file name prefix
  // ----------------------------------------------------------------
  std::ostringstream output_filename_prefix;
  output_filename_prefix << "RESLT/am2pc";
   
  // ----------------------------------------------------------------
  // Prepare the output error file name
  // ----------------------------------------------------------------
  std::ostringstream output_error_filename_prefix;
  output_error_filename_prefix << "RESLT/am2pc_error"; 
  
  // ----------------------------------------------------------------
  // Prepare the output stability file name
  // ----------------------------------------------------------------
  std::ostringstream output_stability_filename_prefix;
  output_stability_filename_prefix << "RESLT/am2pc_stability";
  
  // Create an instance of the problem
  CCStabilityAnalysisProblem<CCBasicODEs> stability_analysis_problem(&odes,
                                                                     time_stepper_pt,
                                                                     output_filename_prefix,
                                                                     output_error_filename_prefix,
                                                                     output_stability_filename_prefix);
  
  // ----------------------------------------------------------------
  // Configure the problem with different time steps to perform
  // stability analisys
  // ----------------------------------------------------------------
  Real time_step = 1.0;
  stability_analysis_problem.time_step() = time_step;
  // Complete problem setup
  stability_analysis_problem.complete_problem_setup();
  // Solve
  stability_analysis_problem.solve();

  time_step = 0.1;
  stability_analysis_problem.time_step() = time_step;
  // Complete problem setup
  stability_analysis_problem.complete_problem_setup();
  // Solve
  stability_analysis_problem.solve();
  
  time_step = 0.01;
  stability_analysis_problem.time_step() = time_step;
  // Complete problem setup
  stability_analysis_problem.complete_problem_setup();
  // Solve
  stability_analysis_problem.solve();
  
  time_step = 0.001;
  stability_analysis_problem.time_step() = time_step;
  // Complete problem setup
  stability_analysis_problem.complete_problem_setup();
  // Solve
  stability_analysis_problem.solve();
    
  std::cout << "[FINISHING UP] ... " << std::endl;
  
  // Free memory
  delete time_stepper_pt;
  time_stepper_pt = 0;
  
 }
  
 {
  std::cout << "--------------------------" << std::endl;
  std::cout << "Backward Differentiation Formula 1 - Fully Implicit - Stability analysis" << std::endl;
  // -----------------------------------------------------------------
  // Instantiation of the ODEs
  // -----------------------------------------------------------------
  CCBasicODEs odes;
   
  // ----------------------------------------------------------------
  // Time stepper
  // ----------------------------------------------------------------
  ACTimeStepper<CCBasicODEs> *time_stepper_pt =
   factory_time_stepper.create_time_stepper("BDF1");
  
  // ----------------------------------------------------------------
  // Prepare the output file name prefix
  // ----------------------------------------------------------------
  std::ostringstream output_filename_prefix;
  output_filename_prefix << "RESLT/bdf1";
  
  // ----------------------------------------------------------------
  // Prepare the output error file name
  // ----------------------------------------------------------------
  std::ostringstream output_error_filename_prefix;
  output_error_filename_prefix << "RESLT/bdf1_error";
  
  // ----------------------------------------------------------------
  // Prepare the output stability file name
  // ----------------------------------------------------------------
  std::ostringstream output_stability_filename_prefix;
  output_stability_filename_prefix << "RESLT/bdf1_stability";
  
  // Create an instance of the problem
  CCStabilityAnalysisProblem<CCBasicODEs> stability_analysis_problem(&odes,
                                                                     time_stepper_pt,
                                                                     output_filename_prefix,
                                                                     output_error_filename_prefix,
                                                                     output_stability_filename_prefix);
  
  // ----------------------------------------------------------------
  // Configure the problem with different time steps to perform
  // stability analisys
  // ----------------------------------------------------------------
  Real time_step = 1.0;
  stability_analysis_problem.time_step() = time_step;
  // Complete problem setup
  stability_analysis_problem.complete_problem_setup();
  // Solve
  stability_analysis_problem.solve();
  
  time_step = 0.1;
  stability_analysis_problem.time_step() = time_step;
  // Complete problem setup
  stability_analysis_problem.complete_problem_setup();
  // Solve
  stability_analysis_problem.solve();
  
  time_step = 0.01;
  stability_analysis_problem.time_step() = time_step;
  // Complete problem setup
  stability_analysis_problem.complete_problem_setup();
  // Solve
  stability_analysis_problem.solve();
  
  time_step = 0.001;
  stability_analysis_problem.time_step() = time_step;
  // Complete problem setup
  stability_analysis_problem.complete_problem_setup();
  // Solve
  stability_analysis_problem.solve();
    
  std::cout << "[FINISHING UP] ... " << std::endl;
  
  // Free memory
  delete time_stepper_pt;
  time_stepper_pt = 0;
  
 }
  
 {
  std::cout << "--------------------------" << std::endl;
  std::cout << "Adams-Moulton 2 or Trapezoidal Rule - Fully Implicit - Stability analysis" << std::endl;
  // -----------------------------------------------------------------
  // Instantiation of the ODEs
  // -----------------------------------------------------------------
  CCBasicODEs odes;
  
  // ----------------------------------------------------------------
  // Time stepper
  // ----------------------------------------------------------------
  ACTimeStepper<CCBasicODEs> *time_stepper_pt =
   factory_time_stepper.create_time_stepper("AM2");

  // ----------------------------------------------------------------
  // Prepare the output file name prefix
  // ----------------------------------------------------------------
  std::ostringstream output_filename_prefix;
  output_filename_prefix << "RESLT/am2";
  
  // ----------------------------------------------------------------
  // Prepare the output error file name
  // ----------------------------------------------------------------
  std::ostringstream output_error_filename_prefix;
  output_error_filename_prefix << "RESLT/am2_error";
  
  // ----------------------------------------------------------------
  // Prepare the output stability file name
  // ----------------------------------------------------------------
  std::ostringstream output_stability_filename_prefix;
  output_stability_filename_prefix << "RESLT/am2_stability";
  
  // Create an instance of the problem
  CCStabilityAnalysisProblem<CCBasicODEs> stability_analysis_problem(&odes,
                                                                     time_stepper_pt,
                                                                     output_filename_prefix,
                                                                     output_error_filename_prefix,
                                                                     output_stability_filename_prefix);
  
  // ----------------------------------------------------------------
  // Configure the problem with different time steps to perform
  // stability analisys
  // ----------------------------------------------------------------
  Real time_step = 1.0;
  stability_analysis_problem.time_step() = time_step;
  // Complete problem setup
  stability_analysis_problem.complete_problem_setup();
  // Solve
  stability_analysis_problem.solve();
  
  time_step = 0.1;
  stability_analysis_problem.time_step() = time_step;
  // Complete problem setup
  stability_analysis_problem.complete_problem_setup();
  // Solve
  stability_analysis_problem.solve();
  
  time_step = 0.01;
  stability_analysis_problem.time_step() = time_step;
  // Complete problem setup
  stability_analysis_problem.complete_problem_setup();
  // Solve
  stability_analysis_problem.solve();
  
  time_step = 0.001;
  stability_analysis_problem.time_step() = time_step;
  // Complete problem setup
  stability_analysis_problem.complete_problem_setup();
  // Solve
  stability_analysis_problem.solve();
  
  std::cout << "[FINISHING UP] ... " << std::endl;
  
  // Free memory
  delete time_stepper_pt;
  time_stepper_pt = 0;
  
 } 

 {
  std::cout << "--------------------------" << std::endl;
  std::cout << "Backward Differentiation Formula 2 - Fully Implicit - Stability analysis" << std::endl;
  // -----------------------------------------------------------------
  // Instantiation of the ODEs
  // -----------------------------------------------------------------
  CCBasicODEs odes;
  
  // ----------------------------------------------------------------
  // Time stepper
  // ----------------------------------------------------------------
  ACTimeStepper<CCBasicODEs> *time_stepper_pt =
   factory_time_stepper.create_time_stepper("BDF2");
  
  // ----------------------------------------------------------------
  // Prepare the output file name prefix
  // ----------------------------------------------------------------
  std::ostringstream output_filename_prefix;
  output_filename_prefix << "RESLT/bdf2";
  
  // ----------------------------------------------------------------
  // Prepare the output error file name
  // ----------------------------------------------------------------
  std::ostringstream output_error_filename_prefix;
  output_error_filename_prefix << "RESLT/bdf2_error";

  // ----------------------------------------------------------------
  // Prepare the output stability file name
  // ----------------------------------------------------------------
  std::ostringstream output_stability_filename_prefix;
  output_stability_filename_prefix << "RESLT/bdf2_stability";
  
  // Create an instance of the problem
  CCStabilityAnalysisProblem<CCBasicODEs> stability_analysis_problem(&odes,
                                                                     time_stepper_pt,
                                                                     output_filename_prefix,
                                                                     output_error_filename_prefix,
                                                                     output_stability_filename_prefix);
  
  // ----------------------------------------------------------------
  // Configure the problem with different time steps to perform
  // stability analisys
  // ----------------------------------------------------------------
  Real time_step = 1.0;
  // Resets time stepper
  //time_stepper_pt->reset();
  stability_analysis_problem.time_step() = time_step;
  // Complete problem setup
  stability_analysis_problem.complete_problem_setup();
  // Solve
  stability_analysis_problem.solve();
  
  time_step = 0.1;
  // Resets time stepper
  time_stepper_pt->reset();
  stability_analysis_problem.time_step() = time_step;
  // Complete problem setup
  stability_analysis_problem.complete_problem_setup();
  // Solve
  stability_analysis_problem.solve();
  
  time_step = 0.01;
  // Resets time stepper
  time_stepper_pt->reset();
  stability_analysis_problem.time_step() = time_step;
  // Complete problem setup
  stability_analysis_problem.complete_problem_setup();
  // Solve
  stability_analysis_problem.solve();
  
  time_step = 0.001;
    // Resets time stepper
  time_stepper_pt->reset();
  stability_analysis_problem.time_step() = time_step;
  // Complete problem setup
  stability_analysis_problem.complete_problem_setup();
  // Solve
  stability_analysis_problem.solve();
  
  std::cout << "[FINISHING UP] ... " << std::endl;
  
  // Free memory
  delete time_stepper_pt;
  time_stepper_pt = 0;
  
 }

 // Finalise scicellxx
 finalise_scicellxx();
 
 return 0;
 
}
