// cuda_functions.h

#ifndef CUDA_FUNCTIONS_H
#define CUDA_FUNCTIONS_H

#ifdef TYPEDEF_REAL_IS_DOUBLE
 typedef double Real;
#else
 typedef float Real;
#endif // #ifdef TYPEDEF_REAL_IS_DOUBLE

#include "../../../src/scicellxx.h"

// __device__ void d_try_lateral_movement(bool* m, const unsigned N, const unsigned L, const unsigned k, const unsigned i, const unsigned e_m);

// __device__ void d_compute_mean_channels_density(const bool* m, const unsigned N, const unsigned L, unsigned e_m, Real* density);

// __device__ void d_mTASEP(bool* d_m, const unsigned e_m, const unsigned i_simulation, 
//             const unsigned d_N, const unsigned d_L,
//             const Real alpha, const Real beta, const Real rho,
//             const Real omega_in, const Real omega_out,
//             bool lateral_movement, Real &mean_current,
//             Real* mean_current_per_channel,                                      
//             unsigned* step_forward_particles_list,
//             unsigned* step_lateral_particles_list);

__device__ __constant__ unsigned d_N;
__device__ __constant__ unsigned d_L;
__device__ __constant__ unsigned d_n_all_configurations;
__device__ __constant__ unsigned d_max_experiments;
__device__ __constant__ unsigned d_n_data_to_gather;
__device__ __constant__ unsigned d_max_simulations_per_experiment;
__device__ __constant__ unsigned d_simulation_step_to_start_gathering_data;
__device__ __constant__ unsigned d_tam_experiment;
__device__ __constant__ unsigned d_tam_simulation;
__device__ __constant__ bool d_lateral_movement;


#endif // CUDA_FUNCTIONS_H
