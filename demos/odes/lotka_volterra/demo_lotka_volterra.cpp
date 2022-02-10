// Include SciCell++ libraries
#include "../../../src/scicellxx.h"

// Odes for Lotka-Volkaterra problem
#include "cc_lotka_volterra_odes.h"

using namespace scicellxx;

/// This class implements inherits from the ACIBVP class, we
/// implement specific functions to solve the Lotka-Volterra equations
template<class EQUATIONS_TYPE>
class CCLotkaVolterraProblem : public virtual ACIBVP<EQUATIONS_TYPE>
{
  
public:
 
 /// Constructor
 CCLotkaVolterraProblem(EQUATIONS_TYPE *odes_pt,
                        ACTimeStepper<EQUATIONS_TYPE> *time_stepper_pt,
                        std::ostringstream &output_filename)
  : ACIBVP<EQUATIONS_TYPE>(odes_pt, time_stepper_pt)
 {
  Output_file.open((output_filename.str()).c_str());
 }
 
 /// Destructor
 ~CCLotkaVolterraProblem()
 {
  Output_file.close();
 }
 
 // Complete problem configuration
 void complete_problem_setup()
 {
   // Prepare time integration data
   const Real initial_time = 0.0;
   const Real final_time = 40.0;
   
   // Initial time
   this->time() = initial_time;
   
   // Final time
   this->Final_time = final_time;
   
   // Set initial conditions
   set_initial_conditions();
 }
 
 // Set initial conditions
 void set_initial_conditions()
 {
  // Initial conditions
  this->u(0) = 2.0; // Initial number of prey
  this->u(1) = 1.0; // Initial number of predators
  //u(0) = 0.9; // Initial number of prey
  //u(1) = 0.9; // Initial number of predators
 }
 
 // Document the solution
 void document_solution()
 {
  // Initial problem configuration
  Output_file << this->time() << "\t" << this->u(0) << "\t" << this->u(1) << std::endl;
 }
 
 // Return the final time
 inline Real final_time() const {return Final_time;}
 
protected:

 // The output file
 std::ofstream Output_file;

 // Final time
 Real Final_time;
 
}; // class CCLotkaVolterraProblem

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
 CCFactoryTimeStepper<CCLotkaVolterraODEs> factory_time_stepper;
 
 // Euler method test
 {
  std::cout << "Euler test" << std::endl;
  // -----------------------------------------------------------------
  // Instantiation of the ODEs
  // -----------------------------------------------------------------
  CCLotkaVolterraODEs odes(1.2, 0.6, 0.8, 0.3);
  //CCLotkaVolterraODEs odes(2.0/3.0, 4.0/3.0, 1.0, 1.0);
  
  // ----------------------------------------------------------------
  // Time stepper
  // ----------------------------------------------------------------
  // Create an instance of the integration method
  ACTimeStepper<CCLotkaVolterraODEs> *time_stepper_pt =
   factory_time_stepper.create_time_stepper("Euler");
  
  // ----------------------------------------------------------------
  // Prepare the output file name
  // ----------------------------------------------------------------
   std::ostringstream output_filename;
   output_filename << "RESLT/euler.dat";
   output_filename.precision(8);
  
   // Create an instance of the problem
   CCLotkaVolterraProblem<CCLotkaVolterraODEs> lotka_volterra_problem(&odes,
                                                                      time_stepper_pt,
                                                                      output_filename);
   
   // Prepare time integration data
   const Real time_step = 0.0625;
   
   // Time step
   lotka_volterra_problem.time_step() = time_step;
   
   // Complete problem configuration
   lotka_volterra_problem.complete_problem_setup();
   
   // Document initial solution
   lotka_volterra_problem.document_solution();
   
   // Flag to indicate whether to continue processing
   bool LOOP = true;
 
   // Main LOOP (loop until reaching final time)
   while(LOOP)
    {
     // Solve (unsteady solve)
     lotka_volterra_problem.solve();
    
     // Update time of the problem
     lotka_volterra_problem.time()+=lotka_volterra_problem.time_step();
    
     // Check whether we have reached the final time
     if (lotka_volterra_problem.time() >= lotka_volterra_problem.final_time())
      {
       LOOP = false;
      }
    
     lotka_volterra_problem.document_solution();
    
    } // while(LOOP)
  
   // Free memory
   delete time_stepper_pt;
   time_stepper_pt = 0;
  
  }
 
 {
  std::cout << "Runge-Kutta 4 test" << std::endl;
  // -----------------------------------------------------------------
  // Instantiation of the ODEs
  // -----------------------------------------------------------------
  CCLotkaVolterraODEs odes(1.2, 0.6, 0.8, 0.3);
  //CCLotkaVolterraODEs odes(2.0/3.0, 4.0/3.0, 1.0, 1.0);
  
  // ----------------------------------------------------------------
  // Time stepper
  // ----------------------------------------------------------------
  ACTimeStepper<CCLotkaVolterraODEs> *time_stepper_pt =
   factory_time_stepper.create_time_stepper("RK4");
  
  // ----------------------------------------------------------------
  // Prepare the output file name
  // ----------------------------------------------------------------
  std::ostringstream output_filename;
  output_filename << "RESLT/rk4.dat";
  output_filename.precision(8);
  
  // Create an instance of the problem
  CCLotkaVolterraProblem<CCLotkaVolterraODEs> lotka_volterra_problem(&odes,
                                                                     time_stepper_pt,
                                                                     output_filename);
  
  // Prepare time integration data
  const Real time_step = 0.0625;
  
  // Time step
  lotka_volterra_problem.time_step() = time_step;
  
  // Complete problem configuration
  lotka_volterra_problem.complete_problem_setup();
  
  // Document initial solution
  lotka_volterra_problem.document_solution();
  
  // Flag to indicate whether to continue processing
  bool LOOP = true;
  
  // Main LOOP (loop until reaching final time)
  while(LOOP)
   {
    // Solve (unsteady solve)
    lotka_volterra_problem.solve();
    
    // Update time of the problem
    lotka_volterra_problem.time()+=lotka_volterra_problem.time_step();
    
    // Check whether we have reached the final time
    if (lotka_volterra_problem.time() >= lotka_volterra_problem.final_time())
     {
      LOOP = false;
     }
    
    lotka_volterra_problem.document_solution();
    
   } // while(LOOP)
  
  // Free memory
  delete time_stepper_pt;
  time_stepper_pt = 0;
  
 }
 
 {
  std::cout << "Adams-Moulton 2 or Trapezoidal Rule - Predictor-Corrector test" << std::endl;
  // -----------------------------------------------------------------
  // Instantiation of the ODEs
  // -----------------------------------------------------------------
  CCLotkaVolterraODEs odes(1.2, 0.6, 0.8, 0.3);
  //CCLotkaVolterraODEs odes(2.0/3.0, 4.0/3.0, 1.0, 1.0);
  
  // ----------------------------------------------------------------
  // Time stepper
  // ----------------------------------------------------------------
  ACTimeStepper<CCLotkaVolterraODEs> *time_stepper_pt =
   factory_time_stepper.create_time_stepper("AM2PC");
  
  // ----------------------------------------------------------------
  // Prepare the output file name
  // ----------------------------------------------------------------
  std::ostringstream output_filename;
  output_filename << "RESLT/am2pc.dat";
  output_filename.precision(8);
  
  // Create an instance of the problem
  CCLotkaVolterraProblem<CCLotkaVolterraODEs> lotka_volterra_problem(&odes,
                                                                     time_stepper_pt,
                                                                     output_filename);
  
  // Prepare time integration data
  const Real time_step = 0.1;
  
  // Initial time step
  lotka_volterra_problem.time_step() = time_step;
  
  // Complete problem configuration
  lotka_volterra_problem.complete_problem_setup();
  
  // Document initial solution
  lotka_volterra_problem.document_solution();
  
  // Flag to indicate whether to continue processing
  bool LOOP = true;
  
  // Main LOOP (loop until reaching final time)
  while(LOOP)
   {
    // Solve (unsteady solve)
    lotka_volterra_problem.solve();
    
    // Update time of the problem
    lotka_volterra_problem.time()+=lotka_volterra_problem.time_step();
    
    // Check whether we have reached the final time
    if (lotka_volterra_problem.time() >= lotka_volterra_problem.final_time())
     {
      LOOP = false;
     }
    
    lotka_volterra_problem.document_solution();
    
   } // while(LOOP)
  
  // Free memory
  delete time_stepper_pt;
  time_stepper_pt = 0;
  
 }
 
 {
  std::cout << "Backward Differentiation Formula 1 - Fully Implicit test" << std::endl;
  // -----------------------------------------------------------------
  // Instantiation of the ODEs
  // -----------------------------------------------------------------
  CCLotkaVolterraODEs odes(1.2, 0.6, 0.8, 0.3);
  //CCLotkaVolterraODEs odes(2.0/3.0, 4.0/3.0, 1.0, 1.0);
  
  // ----------------------------------------------------------------
  // Time stepper
  // ----------------------------------------------------------------
  ACTimeStepper<CCLotkaVolterraODEs> *time_stepper_pt =
   factory_time_stepper.create_time_stepper("BDF1");
  
  // ----------------------------------------------------------------
  // Prepare the output file name
  // ----------------------------------------------------------------
  std::ostringstream output_filename;
  output_filename << "RESLT/bdf1.dat";
  output_filename.precision(8);
  
  // Create an instance of the problem
  CCLotkaVolterraProblem<CCLotkaVolterraODEs> lotka_volterra_problem(&odes,
                                                                     time_stepper_pt,
                                                                     output_filename);
  
  // Prepare time integration data
  const Real time_step = 0.1;
  
  // Time step
  lotka_volterra_problem.time_step() = time_step;
    
  // Complete problem configuration
  lotka_volterra_problem.complete_problem_setup();
  
  // Document initial solution
  lotka_volterra_problem.document_solution();
  
  // Flag to indicate whether to continue processing
  bool LOOP = true;
  
  // Main LOOP (loop until reaching final time)
  while(LOOP)
   {
    // Solve (unsteady solve)
    lotka_volterra_problem.solve();
    
    // Update time of the problem
    lotka_volterra_problem.time()+=lotka_volterra_problem.time_step();
    
    // Check whether we have reached the final time
    if (lotka_volterra_problem.time() >= lotka_volterra_problem.final_time())
     {
      LOOP = false;
     }
    
    lotka_volterra_problem.document_solution();
    
   } // while(LOOP)
  
  // Free memory
  delete time_stepper_pt;
  time_stepper_pt = 0;

 }
 
 {
  std::cout << "Adams-Moulton 2 or Trapezoidal Rule - Fully Implicit test" << std::endl;
  // -----------------------------------------------------------------
  // Instantiation of the ODEs
  // -----------------------------------------------------------------
  CCLotkaVolterraODEs odes(1.2, 0.6, 0.8, 0.3);
  //CCLotkaVolterraODEs odes(2.0/3.0, 4.0/3.0, 1.0, 1.0);
  
  // ----------------------------------------------------------------
  // Time stepper
  // ----------------------------------------------------------------
  ACTimeStepper<CCLotkaVolterraODEs> *time_stepper_pt =
   factory_time_stepper.create_time_stepper("AM2");
  
  // ----------------------------------------------------------------
  // Prepare the output file name
  // ----------------------------------------------------------------
  std::ostringstream output_filename;
  output_filename << "RESLT/am2.dat";
  output_filename.precision(8);
  
  // Create an instance of the problem
  CCLotkaVolterraProblem<CCLotkaVolterraODEs> lotka_volterra_problem(&odes,
                                                                     time_stepper_pt,
                                                                     output_filename);
  
  // Prepare time integration data
  const Real time_step = 0.1;
  
  // Time step
  lotka_volterra_problem.time_step() = time_step;
  
  // Complete problem configuration
  lotka_volterra_problem.complete_problem_setup();
  
  // Document initial solution
  lotka_volterra_problem.document_solution();
  
  // Flag to indicate whether to continue processing
  bool LOOP = true;
  
  // Main LOOP (loop until reaching final time)
  while(LOOP)
   {
    // Solve (unsteady solve)
    lotka_volterra_problem.solve();
    
    // Update time of the problem
    lotka_volterra_problem.time()+=lotka_volterra_problem.time_step();
    
    // Check whether we have reached the final time
    if (lotka_volterra_problem.time() >= lotka_volterra_problem.final_time())
     {
      LOOP = false;
     }
    
    lotka_volterra_problem.document_solution();
    
   } // while(LOOP)
  
  // Free memory
  delete time_stepper_pt;
  time_stepper_pt = 0;

 }
 
 {
  std::cout << "Backward Differentiation Formula 2 - Fully Implicit test" << std::endl;
  // -----------------------------------------------------------------
  // Instantiation of the ODEs
  // -----------------------------------------------------------------
  CCLotkaVolterraODEs odes(1.2, 0.6, 0.8, 0.3);
  //CCLotkaVolterraODEs odes(2.0/3.0, 4.0/3.0, 1.0, 1.0);
  
  // ----------------------------------------------------------------
  // Time stepper
  // ----------------------------------------------------------------
  ACTimeStepper<CCLotkaVolterraODEs> *time_stepper_pt =
   factory_time_stepper.create_time_stepper("BDF2");
  
  // ----------------------------------------------------------------
  // Prepare the output file name
  // ----------------------------------------------------------------
  std::ostringstream output_filename;
  output_filename << "RESLT/bdf2.dat";
  output_filename.precision(8);
  
  // Create an instance of the problem
  CCLotkaVolterraProblem<CCLotkaVolterraODEs> lotka_volterra_problem(&odes,
                                                                     time_stepper_pt,
                                                                     output_filename);
  
  // Prepare time integration data
  const Real time_step = 0.1;
  
  // Time step
  lotka_volterra_problem.time_step() = time_step;
  
  // Complete problem configuration
  lotka_volterra_problem.complete_problem_setup();
  
  // Document initial solution
  lotka_volterra_problem.document_solution();
  
  // Flag to indicate whether to continue processing
  bool LOOP = true;
  
  // Main LOOP (loop until reaching final time)
  while(LOOP)
   {
    // Solve (unsteady solve)
    lotka_volterra_problem.solve();
    
    // Update time of the problem
    lotka_volterra_problem.time()+=lotka_volterra_problem.time_step();
    
    // Check whether we have reached the final time
    if (lotka_volterra_problem.time() >= lotka_volterra_problem.final_time())
     {
      LOOP = false;
     }
    
    lotka_volterra_problem.document_solution();
    
   } // while(LOOP)
  
  // Free memory
  delete time_stepper_pt;
  time_stepper_pt = 0;

 }
 
 // Finalise scicellxx
 finalise_scicellxx();
 
 return 0;
 
}
