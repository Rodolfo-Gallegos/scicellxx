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
/// This demo uses the floor field based on the paper Kirchner, Ansgar
/// and Schadschneider, Andreas, Simulation of evacuation processes
/// using a bionics-inspired cellular automaton model for pedestrian
/// dynamics, Physica A, Elsevier, 2002.

// Include SciCell++ libraries
#include "../../../src/scicellxx.h"

// Include mTASEP algorithm
//#include "cc_mTASEP.h"

// Use the namespace of the framework
using namespace scicellxx;

// Used to define arguments
struct Args {
 argparse::ArgValue<unsigned> L;
 argparse::ArgValue<unsigned> N;
 argparse::ArgValue<Real> alpha_min;
 argparse::ArgValue<Real> alpha_max;
 argparse::ArgValue<unsigned> alpha_n_points;
 argparse::ArgValue<Real> beta_min;
 argparse::ArgValue<Real> beta_max;
 argparse::ArgValue<unsigned> beta_n_points;
 argparse::ArgValue<Real> rho_min;
 argparse::ArgValue<Real> rho_max;
 argparse::ArgValue<unsigned> rho_n_points;
 argparse::ArgValue<Real> omega_in_min;
 argparse::ArgValue<Real> omega_in_max;
 argparse::ArgValue<unsigned> omega_in_n_points;
 argparse::ArgValue<Real> omega_out_min;
 argparse::ArgValue<Real> omega_out_max;
 argparse::ArgValue<unsigned> omega_out_n_points;
 argparse::ArgValue<unsigned> lateral_movement;
 argparse::ArgValue<unsigned> max_experiments;
 argparse::ArgValue<unsigned> max_simulations_per_experiment;
 argparse::ArgValue<unsigned> simulation_step_to_start_gathering_data;
 argparse::ArgValue<std::string> root_output_folder;
 argparse::ArgValue<unsigned> output_experiment_state;
 argparse::ArgValue<unsigned> output_microtubule_state;
};

// Output parameters to a file
void output_parameters_to_file(std::string &filename, const int argc, const char **argv,  struct Args &args)
{
  // Output file
  std::ofstream output_parameters(filename, std::ios_base::out);
  
  // Write the command line as a comment into the parameters file
  output_parameters << "#";
  for (int i = 0; i < argc-1; i++)
   {
    output_parameters << argv[i] << " ";
   }
  output_parameters << argv[argc-1] << std::endl;

  const unsigned precision_real_values = 4;
  
  output_parameters << "L:" << args.L << std::endl;
  output_parameters << "N:" << args.N << std::endl;
  output_parameters << "alpha_min:" << setprecision(precision_real_values) << args.alpha_min << std::endl;
  output_parameters << "alpha_max:" << setprecision(precision_real_values) << args.alpha_max << std::endl;
  output_parameters << "alpha_n_points:" << args.alpha_n_points << std::endl;
  output_parameters << "beta_min:" << setprecision(precision_real_values) << args.beta_min << std::endl;
  output_parameters << "beta_max:" << setprecision(precision_real_values) << args.beta_max << std::endl;
  output_parameters << "beta_n_points:" << args.beta_n_points << std::endl;
  output_parameters << "rho_min:" << setprecision(precision_real_values) << args.rho_min << std::endl;
  output_parameters << "rho_max:" << setprecision(precision_real_values) << args.rho_max << std::endl;
  output_parameters << "rho_n_points:" << args.rho_n_points << std::endl;
  output_parameters << "omega_in_min:" << setprecision(precision_real_values) << args.omega_in_min << std::endl;
  output_parameters << "omega_in_max:" << setprecision(precision_real_values) << args.omega_in_max << std::endl;
  output_parameters << "omega_in_n_points:" << args.omega_in_n_points << std::endl;
  output_parameters << "omega_out_min:" << setprecision(precision_real_values) << args.omega_out_min << std::endl;
  output_parameters << "omega_out_max:" << setprecision(precision_real_values) << args.omega_out_max << std::endl;
  output_parameters << "omega_out_n_points:" << args.omega_out_n_points << std::endl;
  output_parameters << "lateral_movement:" << args.lateral_movement << std::endl;
  output_parameters << "max_experiments:" << args.max_experiments << std::endl;
  output_parameters << "max_simulations_per_experiment:" << args.max_simulations_per_experiment << std::endl;
  output_parameters << "simulation_step_to_start_gathering_data:" << args.simulation_step_to_start_gathering_data << std::endl;
  output_parameters << "root_output_folder:" << args.root_output_folder << std::endl;
  output_parameters << "output_experiment_state:" << args.output_experiment_state << std::endl;
  output_parameters << "output_microtubule_state:" << args.output_microtubule_state << std::endl;
  
  // Close the parameters file
  output_parameters.close(); 
}

/// Tries to perform a lateral movement on the particule at position
/// (k, i) in the microtubule
void try_lateral_movement(bool **m, const unsigned N, const unsigned L, const unsigned k, const unsigned i)
{
 /*
   Performs lateral movement on a TASEP model, this check for up/down
   possibilities, does not accounts for diagonal movement
 */
 
 // Double-check there is a particle at the current position
 if (m[k][i] == 0)
  {
   return;
  }

 // Store the indices for the microtubule above and below
 unsigned index_microtubule_above = k - 1;
 unsigned index_microtubule_below = k + 1;
 
 // Correct the indexes for the microtubules above and below
 // (periodic-boundary conditions)
 if (k == 0)
  {
   index_microtubule_above = N - 1;
  }
 
 if (k == N-1)
  {
   index_microtubule_below = 0;
  }
 
 // Used to get a seed for the random number engine
 std::random_device rd;
 // Standard mersenne_twister_engine seeded with rd()
 std::mt19937 gen(rd());
 
 // Use dist to generate a random number into a Real in the range
 // [0,1)
 std::uniform_real_distribution<> dis(0.0, 1.0);
 
 // Generate a random number to choose between the above or below
 // microtubule
 const Real r = dis(gen);
 
 // First choose the microtubule from above (preferred due to probability)
 if (r <= 0.5)
  {
   // If the space at the above microtubule is free then move there!
   if (m[index_microtubule_above][i] == false)
    {
     m[k][i] = false;
     m[index_microtubule_above][i] = true;
    }
   // If the space at the below microtubule is free then move there!
   else if (m[index_microtubule_below][i] == false)
    {
     m[k][i] = false;
     m[index_microtubule_below][i] = true;
    }
  } // if (r <= 0.5)
 else // First choose the microtubule from below (preferred due to
      // probability)
  {
   if (m[index_microtubule_below][i] == false)
    {
     m[k][i] = false;
     m[index_microtubule_below][i] = true;
    }
   else if (m[index_microtubule_above][i] == false)
    {
     m[k][i] = false;
     m[index_microtubule_above][i] = true;
    }
  } // else if (r <= 0.5)
 
}

/// Perform mTASEP method
void mTASEP(bool **m, const unsigned N, const unsigned L,
            const Real alpha, const Real beta, const Real rho,
            const Real omega_in, const Real omega_out,
            bool lateral_movement)
{
 /*
    Applies TASEP algorithm to a multichannel microtubule

    Return the updated microtubule as an numpy matrix

    This version implements the "(d) parallel_update" strategy as
    described in section 2.1 of "The Asymmetric Exclusion Process:
    Comparison of Update Procedures, N. Rajewsky et. al., Journal of
    Statistical Physics, Vol. 92, 1998"
 */
 
 // Used to get a seed for the random number engine
 std::random_device rd;
 // Standard mersenne_twister_engine seeded with rd()
 std::mt19937 gen(rd());
 
 // Use dist to generate a random number into a Real in the range
 // [0,1)
 std::uniform_real_distribution<> dis(0.0, 1.0);
 
 // Used to generate a random position in the microtubule (avoid first
 // and last positions)
 std::uniform_int_distribution<> dis_microtubule_size(1, L-2);
 
 // Perform the method for each channel
 for (unsigned k = 0; k < N; k++)
  {
   // *******************************************
   // Update the boundaries (left and right)
   // *******************************************
   
   // -------------------------------------------
   /// LEFT END
   // -------------------------------------------
   
   // Compute a probability to add a molecule at the beginning of the
   // microtubule
   const Real a = dis(gen);
   // A flag indicating whether a molecule was added at the beginning
   // of the microtubule
   bool added_to_start = false;
   // Is a <= alpha and the first space is free?
   if (a <= alpha && m[k][0] == 0)
    {
     // Add a molecule to the start
     m[k][0] = 1;
     // Indicate we added molecule to the beginning of the microtubule
     // so there is no need to update its position
     added_to_start = true;
    }
   
   // -------------------------------------------
   /// RIGHT END
   // -------------------------------------------
   
   // Compute a probability to remove the last molecule of the
   // microtubule
   const Real b = dis(gen);
   // A flag indicating whether a molecule was removed from the last
   // cell of the microtubule
   bool removed_from_end = false;
   // Is a <= alpha and the first space is free?
   if (b <= beta && m[k][L-1] == 1)
    {
     // Remove the molecule from the microtubule
     m[k][L-1] = 0;
     // Indicate we removed a molecule from the last position of the
     // microtubule so other molecule should not step in this cell
     removed_from_end = true;
    }
   
   // *******************************************
   // Update particles in the microtubule
   // *******************************************
   
   // Update particles in the microtubule using a parallel update
   // strategy (delay movement, but keep track in a list of those
   // particles to update position)
   
   // -----------------------------------------------------------
   
   // Update the first and last indexes based on the added_to_start
   // and removed_from_end flags
   unsigned end_index = L - 2;
   if (removed_from_end)
    {
     end_index = L - 3;     
    }
   
   unsigned start_index = 0;
   if (added_to_start)
    {
     start_index = 1;
    }
   
   // Vector with particles' indexes that will move (step-forward)
   std::vector<unsigned> step_forward_particles_list;
   step_forward_particles_list.reserve(L);
      
   // Vector with particles' indexes that will try lateral movement
   std::vector<unsigned> step_lateral_particles_list;
   step_lateral_particles_list.reserve(L);
   
   // Update the microtubule from left to right
   for (unsigned i = start_index; i <= end_index; i++)
    {
     // Check whether there is a particle at the current cell space
     if (m[k][i] == 1)
      {
       // Check whether there is a free space at the next cell space
       if (m[k][i+1] == 0)
        {
         // Compute a probability to move to the next space
         const Real p = dis(gen);
         if (p <= rho)
          {
           // Add the particle index to the step forward list
           step_forward_particles_list.push_back(i);
          }

         // Skip the cell at the next space since there is no particle
         // there
         ++i;
         
        } // if (m[k][i+1] == 0)
       else
        {
         // Try lateral movement since there is a particle on the next
         // cell
         step_lateral_particles_list.push_back(i);
         
        } // else if (m[k][i+1] == 0)
       
      } // if (m[k][i] == 1)

    } // for (i <= end_index)
        
   // Perform movement of particles based on the step forward list
   // Cache the size of the list
   const unsigned step_forward_particle_list_size = step_forward_particles_list.size();
   for (unsigned j = 0; j < step_forward_particle_list_size; j++)
    {
     // Get the index of the particle
     const unsigned i = step_forward_particles_list[j];
     m[k][i] = 0;
     m[k][i+1] = 1;
    }
   
   // *******************************************
   // Apply lateral movement
   // *******************************************
   if (lateral_movement)
    {
     const unsigned step_lateral_particle_list_size = step_lateral_particles_list.size();
     for (unsigned j = 0; j < step_lateral_particle_list_size; j++)
      {
       // Get the index of the particle
       const unsigned i = step_lateral_particles_list[j];
       try_lateral_movement(m, N, L, k, i);
      }
     
    }
   
   // *******************************************
   // Apply TASEP-LK rules
   // *******************************************

   // -------------------------------------------
   /// Attach (omega in)
   // -------------------------------------------
   if (omega_in > 0.0)
    {
     // Choose a random position at the microtubule
     const unsigned r_pos = dis_microtubule_size(gen);

     // Generate a random number
     const Real r = dis(gen);

     // if r <= omega_in and m[k][r_pos] == 0 then add a molecule to
     // the microtubule in that position
     if (r <= omega_in && m[k][r_pos] == 0)
      {
       // Attach a molecule
       m[k][r_pos] = 1;
      }
    }
   
   // -------------------------------------------
   /// Detach (omega out)
   // -------------------------------------------
   if (omega_out > 0.0)
    {
     // Choose a random position at the microtubule
     const unsigned r_pos = dis_microtubule_size(gen);

     // Generate a random number
     const Real r = dis(gen);

     // if r <= omega_out and m[k][r_pos] == 1 then remove the
     // molecule from that position of the microtubule
     if (r <= omega_out && m[k][r_pos] == 1)
      {
       // Detach the molecule
       m[k][r_pos] = 0;
      }
    }
   
  } // for (k < N)
 
}

/// Output matrix into a csv file
void matrix_to_csv_file(bool **m, const unsigned nrows, const unsigned ncolumns, std::string &file_name)
{
 // Create file
 std::ofstream output_file(file_name, std::ios_base::out);
 
 for (unsigned i = 0; i < nrows; i++)
  {
   for (unsigned j = 0; j < ncolumns-1; j++)
    {
     output_file << m[i][j] << ",";
    }
   // The last element without the ','
   output_file << m[i][ncolumns-1] << std::endl;
  }
 
 // Close the file
 output_file.close();
 
} // matrix_to_csv_file

//#define OUTPUT_ONE_SINGLE_RUN

int main(int argc, const char** argv)
{
 // Initialise scicellxx
 initialise_scicellxx();
 
 // Instantiate parser
 Args args;
 auto parser = argparse::ArgumentParser(argv[0], "mtasep-LK model for a multichannel microtubule with parallel update");
 
 // Add arguments
 
 // Optional
 parser.add_argument<unsigned>(args.L, "--L")
  .help("Size of the microtubule")
  .default_value("101");
 
 parser.add_argument<unsigned>(args.N, "--N")
  .help("Number of channels on the microtubule")
  .default_value("4");
 
 parser.add_argument<Real>(args.alpha_min, "--alpha_min")
  .help("The minimum probability to introduce a particle to the microtubule (at the first --left-- position of the microtubule)")
  .default_value("0.1");
 
 parser.add_argument<Real>(args.alpha_max, "--alpha_max")
  .help("The maximum probability to introduce a particle to the microtubule (at the first --left-- position of the microtubule)")
  .default_value("1.0");
 
 parser.add_argument<unsigned>(args.alpha_n_points, "--alpha_n_points")
  .help("The number of points in the range alpha_min to alpha_max")
  .default_value("10");

 parser.add_argument<Real>(args.beta_min, "--beta_min")
  .help("The minimum probability to remove a particle from the microtubule (from the last --right-- position of the microtubule)")
  .default_value("0.1");

 parser.add_argument<Real>(args.beta_max, "--beta_max")
  .help("The maximum probability to remove a particle from the microtubule (from the last --right-- position of the microtubule)")
  .default_value("1.0");
 
 parser.add_argument<unsigned>(args.beta_n_points, "--beta_n_points")
  .help("The number of points in the range beta_min to beta_max")
  .default_value("10");
 
 parser.add_argument<Real>(args.rho_min, "--rho_min")
  .help("The minimum probability for a particle to step-forward")
  .default_value("0.8");
 
 parser.add_argument<Real>(args.rho_max, "--rho_max")
  .help("The maximum probability for a particle to step-forward")
  .default_value("1.0");
 
 parser.add_argument<unsigned>(args.rho_n_points, "--rho_n_points")
  .help("The number of points in the range rho_min to rho_max")
  .default_value("3");
 
 parser.add_argument<Real>(args.omega_in_min, "--omega_in_min")
  .help("The minimum probability for a particle to attach to a site in the microtubule")
  .default_value("0.0");
 
 parser.add_argument<Real>(args.omega_in_max, "--omega_in_max")
  .help("The maximum probability for a particle to attach to a site in the microtubule")
  .default_value("1.0");
 
 parser.add_argument<unsigned>(args.omega_in_n_points, "--omega_in_n_points")
  .help("The number of points in the range omega_in_min to omega_in_max")
  .default_value("11");
 
 parser.add_argument<Real>(args.omega_out_min, "--omega_out_min")
  .help("The minimum probability for a particle to detach from a site in the microtubule")
  .default_value("0.0");
 
 parser.add_argument<Real>(args.omega_out_max, "--omega_out_max")
  .help("The maximum probability for a particle to detach from a site in the microtubule")
  .default_value("1.0");
 
 parser.add_argument<unsigned>(args.omega_out_n_points, "--omega_out_n_points")
  .help("The number of points in the range omega_out_min to omega_out_max")
  .default_value("11");
 
 parser.add_argument<unsigned>(args.lateral_movement, "--lateral_movement")
  .help("Enables/disables lateral movement")
  .default_value("0");

 parser.add_argument<unsigned>(args.max_experiments, "--max_experiments")
  .help("Set the maximum number of experiments to perform")
  .default_value("1");
 
 parser.add_argument<unsigned>(args.max_simulations_per_experiment, "--max_simulations_per_experiment")
  .help("Set the maximum number of simulations per experiments")
  .default_value("200");
 
 parser.add_argument<unsigned>(args.simulation_step_to_start_gathering_data, "--simulation_step_to_start_gathering_data")
  .help("Set the simulation step to start gathering data (must be smaller than max_simulations_per_experiment)")
  .default_value("1");

 parser.add_argument<std::string>(args.root_output_folder, "--root_output_folder")
  .help("The root output folder")
  .default_value("RESLT");
 
 parser.add_argument<unsigned>(args.output_experiment_state, "--output_experiment_state")
  .help("Enables/disables output of experiments state (averaged space-time diagrams for all channels). Disable when performing large simulations")
  .default_value("0");
 
  parser.add_argument<unsigned>(args.output_microtubule_state, "--output_microtubule_state")
  .help("Enables/disables output of microtubule state (space-time diagrams for all channels). Disable when performing large simulations")
  .default_value("0");
 
  // Parse the input arguments
  parser.parse_args(argc, argv);

  // Configure problem
  const unsigned L = args.L; // Size of the microtubule
  const unsigned N = args.N; // Number of channels of the microtubule
  const Real alpha_min = args.alpha_min;
  const Real alpha_max = args.alpha_max;
  const unsigned alpha_n_points = args.alpha_n_points;
  const Real beta_min = args.beta_min;
  const Real beta_max = args.beta_max;
  const unsigned beta_n_points = args.beta_n_points;
  const Real rho_min = args.rho_min;
  const Real rho_max = args.rho_max;
  const unsigned rho_n_points = args.rho_n_points;
  const Real omega_in_min = args.omega_in_min;
  const Real omega_in_max = args.omega_in_max;
  const unsigned omega_in_n_points = args.omega_in_n_points;
  const Real omega_out_min = args.omega_out_min;
  const Real omega_out_max = args.omega_out_max;
  const unsigned omega_out_n_points = args.omega_out_n_points;
  bool lateral_movement = false;
  if (args.lateral_movement)
   {
    lateral_movement = true;
   }
  const unsigned max_experiments = args.max_experiments;
  const unsigned max_simulations_per_experiment = args.max_simulations_per_experiment;
  const unsigned simulation_step_to_start_gathering_data = args.simulation_step_to_start_gathering_data;
  std::string root_output_folder(args.root_output_folder);
  bool output_experiment_state = false;
  if (args.output_experiment_state)
   {
    output_experiment_state = true;
   }
  bool output_microtubule_state = false;
  if (args.output_microtubule_state)
   {
    output_microtubule_state = true;
   }
  
  // Create output directory
  SciCellxxFileSystem::create_directory(root_output_folder);
  
  // Output parameters to a file
  std::string parameters_filename(root_output_folder + "/parameters.txt");
  output_parameters_to_file(parameters_filename, argc, argv, args);

  // ------------------------------------------------------------
  // Generate all configurations as the cartesian product of the
  // ranges of parameters
  // ------------------------------------------------------------
  
  // Alphas
  std::vector<Real> alphas;
  SciCellxxLinearSpace::create_linear_space(alphas, alpha_min, alpha_max, alpha_n_points);
  // Betas
  std::vector<Real> betas;
  SciCellxxLinearSpace::create_linear_space(betas, beta_min, beta_max, beta_n_points);
  
  // Rhos
  std::vector<Real> rhos;
  SciCellxxLinearSpace::create_linear_space(rhos, rho_min, rho_max, rho_n_points);
  
  // Omegas_in
  std::vector<Real> omegas_in;
  SciCellxxLinearSpace::create_linear_space(omegas_in, omega_in_min, omega_in_max, omega_in_n_points);
  
  // Omegas_out
  std::vector<Real> omegas_out;
  SciCellxxLinearSpace::create_linear_space(omegas_out, omega_out_min, omega_out_max, omega_out_n_points);
  
  SciCellxxLinearSpace::print_linear_space(alphas);
  SciCellxxLinearSpace::print_linear_space(betas);
  SciCellxxLinearSpace::print_linear_space(rhos);
  SciCellxxLinearSpace::print_linear_space(omegas_in);
  SciCellxxLinearSpace::print_linear_space(omegas_out);

  return 0;
  
  //
  // HERE HERE HERE

  
  const unsigned all_configurations = 1;
  
  // Run all configurations
  for (unsigned i_configuration = 0; i_configuration < all_configurations; i_configuration++)
   {
    
    // Get values for a simulation configuration
    Real alpha = alpha_min;
    Real beta = beta_min;
    Real rho = rho_min;
    Real omega_in = omega_in_min;
    Real omega_out = omega_out_min;
    
    // Run all experiments for the current configuration
    for (unsigned i_experiment = 0; i_experiment < max_experiments; i_experiment++)
     {
      // Construct the microtubule with N channels and L cells on each channel
      bool **m = new bool*[N];
      for (unsigned i_channel = 0; i_channel < N; i_channel++)
       {
        m[i_channel] = new bool[L];
       }

      // Initialize microtubule with zeroes
      for (unsigned i_channel = 0; i_channel < N; i_channel++)
       {
        for (unsigned i_cell = 0; i_cell < L; i_cell++)
         {
          m[i_channel][i_cell] = 0;
         }
       }
      
      // Start simulation
      for (unsigned i_simulation_step = 0; i_simulation_step < max_simulations_per_experiment; i_simulation_step++)
       {
        // Apply mTASEP
        mTASEP(m, N, L, alpha, beta, rho, omega_in, omega_out, lateral_movement);
        
        // Store csv file with the microtubule state
        if (output_microtubule_state)
         {
          // HERE HERE HERE HERE
          // Consider zeroes for filenames
          // Check for outputfolder
          //                     experiment_folder_name = "experiment_{:05d}_a_{:.4f}_b_{:.4f}_r_{:.4f}".format(iexperiment, alpha, beta, rho)
          //                      path_experiment_folder = os.path.join(current_path, experiment_folder_name)
          //          microtubule_filename = path_experiment_folder + "/microtubule_{:05d}.csv".format(simulation)
          const unsigned width_number = 5;
          const char fill_char = '0';
          std::ostringstream ss;
          ss << std::setw(width_number) << std::setfill(fill_char) << std::to_string(i_simulation_step);
          
          //std::string csv_filename(root_output_folder + std::string("/example_") + std::to_string(i_simulation_step) + std::string(".csv"));
          std::string csv_filename(root_output_folder + std::string("/example_") + ss.str() + std::string(".csv"));
          matrix_to_csv_file(m, N, L, csv_filename);
         }
       }

      // Free memory for multichannel-microtubule
      for (unsigned i_channel = 0; i_channel < N; i_channel++)
       {
        delete [] m[i_channel];
       }
      delete [] m;
      
     } // for (i_experiment < max_experiments)
    
   } // for (i_configuration < all_configurations)
  
  // Finalise chapcom
  finalise_scicellxx();
  
  return 0;
  
}

