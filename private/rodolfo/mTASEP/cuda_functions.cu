#include "cuda_functions.h"
#include "../../../src/scicellxx.h"
#include <curand_kernel.h>


__device__ void d_compute_mean_channels_density(const bool* m, const unsigned N, const unsigned L, unsigned e_m, Real* density)
{
    // Moverse a lo largo del microtúbulo y calcular la densidad promedio
    for (unsigned i = 0; i < L; i++)
    {
        Real sum_occupation = 0.0;

        for (unsigned k = 0; k < N; k++)
        {
            // Calcular el índice del elemento actual en el arreglo plano
            unsigned idx = e_m + k * L + i;

            sum_occupation += static_cast<Real>(m[idx]);
        }
        // Calcular la densidad promedio para esta celda
        density[e_m + i * L] = sum_occupation / static_cast<Real>(N);
    }
}


// CUDA function for lateral movement
__device__ void d_try_lateral_movement(bool *m, const unsigned N, const unsigned L, const unsigned k, const unsigned i, const unsigned e_m)
{
    // Calculate the linear index
    unsigned index = e_m + k * L + i;

    // Double-check there is a particle at the current position
    if (!m[index])
    {
        return;
    }

    // Store the indices for the microtubule above and below
    unsigned index_microtubule_above = (k == 0) ? N - 1 : k - 1;
    unsigned index_microtubule_below = (k == N - 1) ? 0 : k + 1;

    // Generate a random number to choose between the above or below microtubule
    unsigned seed = threadIdx.x + blockIdx.x * blockDim.x;
    curandState_t state;
    curand_init(seed, 0, 0, &state);
    const Real r = curand_uniform(&state);

    // First choose the microtubule from above (preferred due to probability)
    if (r <= 0.5f)
    {
        // If the space at the above microtubule is free then move there!
        if (!m[e_m + index_microtubule_above * L + i])
        {
            m[index] = false;
            m[e_m + index_microtubule_above * L + i] = true;
        }
        // If the space at the below microtubule is free then move there!
        else if (!m[e_m + index_microtubule_below * L + i])
        {
            m[index] = false;
            m[e_m + index_microtubule_below * L + i] = true;
        }
    }
    // First choose the microtubule from below (preferred due to probability)
    else
    {
        if (!m[e_m + index_microtubule_below * L + i])
        {
            m[index] = false;
            m[e_m + index_microtubule_below * L + i] = true;
        }
        else if (!m[e_m + index_microtubule_above * L + i])
        {
            m[index] = false;
            m[e_m + index_microtubule_above * L + i] = true;
        }
    }
}


__device__ void d_mTASEP(bool* d_m, const unsigned e_m, const unsigned i_simulation, 
            const unsigned d_N, const unsigned d_L,
            const Real alpha, const Real beta, const Real rho,
            const Real omega_in, const Real omega_out,
            bool lateral_movement, Real &mean_current,
            Real* mean_current_per_channel,                                      
            unsigned* step_forward_particles_list,
            unsigned* step_lateral_particles_list)
{
 /*
    Applies TASEP algorithm to a multichannel microtubule

    Return the updated microtubule as an numpy matrix

    This version implements the "(d) parallel_update" strategy as
    described in section 2.1 of "The Asymmetric Exclusion Process:
    Comparison of Update Procedures, N. Rajewsky et. al., Journal of
    Statistical Physics, Vol. 92, 1998"
 */
 
    // Inicializar el estado del generador de números aleatorios
    curandState state;
    unsigned long long seed = threadIdx.x + blockIdx.x * blockDim.x;
    curand_init(seed, 0, 0, &state);

    // Generar un número real aleatorio uniformemente distribuido en el rango [0,1)
    Real r = curand_uniform(&state);

    // Generar un número entero aleatorio uniformemente distribuido en el rango [0, L-1]
    unsigned random_index = static_cast<unsigned>(curand_uniform(&state) * d_L);
    
    // Compute the current for each channel on the microtubule

    // Perform the method for each channel
    for (unsigned k = 0; k < d_N; k++)
    {
        // *******************************************
        // Apply TASEP-LK rules
        // *******************************************
        
        // -------------------------------------------
        /// Attach (omega in)
        // -------------------------------------------
        
        // Keep track of the position of the particle attached by omega_in
        // because the just attached particle cannot move
        bool omega_in_attached_a_particle = false;
        unsigned r_pos_omega_in = 0;
        
        if (omega_in > 0.0)
        {
            // Choose a random position at the microtubule
            const unsigned r_pos = static_cast<unsigned>(curand_uniform(&state) * d_L);
            
            // Generate a random number
            const Real r = curand_uniform(&state);
            
            // if r <= omega_in and d_m[k][r_pos] == 0 then add a particle to
            // the microtubule in that position
            if (r <= omega_in && d_m[e_m +  k * d_L + r_pos] == 0)
            {
                // Attach a particle
                d_m[e_m +  k * d_L + r_pos] = 1;
                
                // Keep track of the position of the particle attached by omega
                // in
                r_pos_omega_in = r_pos;
                omega_in_attached_a_particle = true;
            }
        }
        
        // -------------------------------------------
        /// Detach (omega out)
        // -------------------------------------------
        
        // Keep track of the position free by omega_out becasue this
        // position cannot be occupied by other particles in the current
        // simulation step
        bool omega_out_dettached_a_particle = false;
        unsigned r_pos_omega_out = 0;
        
        if (omega_out > 0.0)
        {
            // Choose a random position at the microtubule
            const unsigned r_pos = static_cast<unsigned>(curand_uniform(&state) * d_L);
            
            // Generate a random number
            const Real r = curand_uniform(&state);
            
            // if r <= omega_out and d_m[k][r_pos] == 1 then remove the
            // particle from that position of the microtubule ALSO check that
            // the particle is not there becasue it was just attached by the
            // omega_in process
            if (r <= omega_out && d_m[e_m +  k * d_L + r_pos] == 1 && !(omega_in_attached_a_particle && r_pos == r_pos_omega_in))
            {
            // Detach the particle
            d_m[e_m +  k * d_L + r_pos] = 0;
            
            r_pos_omega_out = r_pos;
            omega_out_dettached_a_particle = true;
            }    
        }
        
        // *******************************************
        // Update the boundaries (left and right)
        // *******************************************
        
        // -------------------------------------------
        /// LEFT END
        // -------------------------------------------
        
        // Compute a probability to add a particle at the beginning of the
        // microtubule
        const Real a = curand_uniform(&state);
        // A flag indicating whether a particle was added at the beginning
        // of the microtubule
        bool added_to_start = false;
        // Is a <= alpha and the first space is free? ALSO check whether
        // the space is free due to the dettaching process by omega_out, if
        // that is the case then we cannot add a particle to the start
        if (a <= alpha && d_m[e_m + k * d_L] == 0 && !(omega_out_dettached_a_particle && r_pos_omega_out == 0))
            {
            // Add a particle to the start
            d_m[e_m + k * d_L] = 1;
            // Indicate we added a particle at the beginning of the
            // microtubule so there is no need to update its position
            added_to_start = true;
            
            }
        
        // -------------------------------------------
        /// RIGHT END
        // -------------------------------------------
        
        // Compute a probability to remove the last particle of the
        // microtubule
        const Real b = curand_uniform(&state);
        // A flag indicating whether a particle was removed from the last
        // cell of the microtubule
        bool removed_from_end = false;
        // Is b <= beta and the last space occupied, ALSO, was the particle
        // not attached by omega_in?
        if (b <= beta && d_m[e_m + k * d_L + d_L-1] == 1 && !(omega_in_attached_a_particle && r_pos_omega_in == d_L-1))
        {
            // Remove the particle from the microtubule
            d_m[e_m + k * d_L + d_L-1] = 0;
            // Indicate we removed a particle from the last position of the
            // microtubule so other particle should not step in this cell
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
        unsigned end_index = d_L - 2;
        if (removed_from_end)
        {
            end_index = d_L - 3;     
        }
        
        unsigned start_index = 0;
        if (added_to_start)
        {
            start_index = 1;
        }

        unsigned forward_count = 0;
        unsigned lateral_count = 0;
        
        // // Vector with particles' indexes that will move (step-forward)
        // std::vector<unsigned> step_forward_particles_list;
        // step_forward_particles_list.reserve(d_L);
            
        // // Vector with particles' indexes that will try lateral movement
        // std::vector<unsigned> step_lateral_particles_list;
        // step_lateral_particles_list.reserve(d_L);
        
        // Update the microtubule from left to right
        for (unsigned i = start_index; i <= end_index; i++)
        {
            // Check whether there is a particle at the current cell space
            // and ensure that particle was not introduced by the omega_in
            // process
            if (d_m[e_m + k * d_L + i] == 1 && !(omega_in_attached_a_particle && r_pos_omega_in == i))
            {
                // Check whether there is a free space at the next cell space
                // and ensure that space is not there becase the omega_out
                // process
                if (d_m[e_m + k * d_L + i + 1] == 0 && !(omega_out_dettached_a_particle && r_pos_omega_out == i+1))
                {
                    // Compute a probability to move to the next space
                    const Real p = curand_uniform(&state);
                    if (p <= rho)
                    {
                        // Add the particle index to the step forward list
                        // step_forward_particles_list.push_back(i);
                        step_forward_particles_list[e_m + k * d_L + forward_count] = i;
                        forward_count++;
                    }

                    // Skip the cell at the next space since there is no particle
                    // there
                    ++i;
                
                } // if (d_m[k][i+1] == 0)
                else
                {
                    // Try lateral movement since there is a particle on the next
                    // cell OR the particle on the next cell was just removed by
                    // the omega_out process
                    // step_lateral_particles_list.push_back(i);
                    step_lateral_particles_list[e_m + k * d_L + lateral_count] = i;
                    lateral_count++;
                    
                } // else if (d_m[k][i+1] == 0 && !(omega_out_dettached_a_particle && r_pos_omega_out == i+1))
            
            } // if (d_m[k][i] == 1 && !(omega_in_attached_a_particle && r_pos_omega_in == i))

        } // for (i <= end_index)
                
        // Perform movement of particles based on the step forward list
        // Cache the size of the list
        const unsigned step_forward_particle_list_size = forward_count;

        for (unsigned j = 0; j < step_forward_particle_list_size; j++)
        {
        // Get the index of the particle
        const unsigned i = step_forward_particles_list[e_m + k * d_L + j];
        d_m[e_m + k * d_L + i] = 0;
        d_m[e_m + k * d_L + i + 1] = 1;
        // Add the movement to the current per channel vector
        mean_current_per_channel[e_m + k * d_L]+=1;
        }
        
        // *******************************************
        // Apply lateral movement
        // *******************************************
        if (d_lateral_movement)
        {
            const unsigned step_lateral_particle_list_size = lateral_count;

            for (unsigned j = 0; j < step_lateral_particle_list_size; j++)
            {
                // Get the index of the particle
                const unsigned i = step_lateral_particles_list[e_m + k * d_L + j];

                d_try_lateral_movement(d_m, d_N, d_L, k, i, e_m);
            }
        }
    
    } // for (k < d_N)
    
    // Compute the averaged current on all the microtubules
    mean_current = 0.0;
    for (unsigned i = 0; i < d_N; i++)
    {
        mean_current+=mean_current_per_channel[e_m + i * d_L];
    }
    
    const Real factor = 1.0/static_cast<Real>(d_L*d_N);
    
    mean_current=mean_current*factor;
}
