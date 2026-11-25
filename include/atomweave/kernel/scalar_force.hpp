#pragma once

#include <algorithm>
#include <cmath>
#include "domain/periodic_cube.hpp"
#include "physics/lennard_jones.hpp"
#include "state/particle_state.hpp"


namespace atomweave::kernel
{

template <typename TReal> struct ScalarForceKernel
{
    TReal evaluate(state::ParticleState<TReal> &state,
                   const domain::PeriodicCube<TReal> &periodic_cube,
                   const physics::LennardJones<TReal> &lennard_jones)
    {
        TReal total_energy = 0.0;

        std::fill(state.force.x.begin(), state.force.x.end(), 0);
        std::fill(state.force.y.begin(), state.force.y.end(), 0);
        std::fill(state.force.z.begin(), state.force.z.end(), 0);


        for (std::size_t j = 1; j < state.size(); ++j)
        {
            for (std::size_t i = 0; i < j; ++i)
            {
                TReal dx = periodic_cube.minimum_image(state.position.x[i] -
                                                       state.position.x[j]);
                TReal dy = periodic_cube.minimum_image(state.position.y[i] -
                                                       state.position.y[j]);
                TReal dz = periodic_cube.minimum_image(state.position.z[i] -
                                                       state.position.z[j]);

                TReal r_squared = dx * dx + dy * dy + dz * dz;

                if (r_squared > lennard_jones.cutoff_squared())
                {
                    continue;
                }

                TReal force_coeff = lennard_jones.force_coefficient(r_squared);

                TReal F_x = force_coeff * dx;
                TReal F_y = force_coeff * dy;
                TReal F_z = force_coeff * dz;

                state.force.x[i] += F_x;
                state.force.y[i] += F_y;
                state.force.z[i] += F_z;

                state.force.x[j] -= F_x;
                state.force.y[j] -= F_y;
                state.force.z[j] -= F_z;

                total_energy += lennard_jones.potential(std::sqrt(r_squared));
            }
        }
        return total_energy;
    }
};


} // namespace atomweave::kernel
