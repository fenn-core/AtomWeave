#pragma once

#include "domain/periodic_cube.hpp"
#include "physics/lennard_jones.hpp"
#include "state/particle_state.hpp"


namespace atomweave::dynamics
{

template <typename TReal, typename TKernel> struct VelocityVerlet
{
    TReal step(state::ParticleState<TReal> &state,
               const domain::PeriodicCube<TReal> &periodic_cube,
               const physics::LennardJones<TReal> &lennard_jones, TKernel &kernel,
               TReal dt, TReal mass)
    {
        TReal half_dt_inv_mass = 0.5 * dt * (1.0 / mass);

        for (std::size_t i = 0; i < state.size(); ++i)
        {
            state.velocity.x[i] += half_dt_inv_mass * state.force.x[i];
            state.velocity.y[i] += half_dt_inv_mass * state.force.y[i];
            state.velocity.z[i] += half_dt_inv_mass * state.force.z[i];

            state.position.x[i] = periodic_cube.wrap_coordinate(state.position.x[i] +
                                                                dt * state.velocity.x[i]);
            state.position.y[i] = periodic_cube.wrap_coordinate(state.position.y[i] +
                                                                dt * state.velocity.y[i]);
            state.position.z[i] = periodic_cube.wrap_coordinate(state.position.z[i] +
                                                                dt * state.velocity.z[i]);
        }

        TReal energy = kernel.evaluate(state, periodic_cube, lennard_jones);

        for (std::size_t i = 0; i < state.size(); ++i)
        {
            state.velocity.x[i] += half_dt_inv_mass * state.force.x[i];
            state.velocity.y[i] += half_dt_inv_mass * state.force.y[i];
            state.velocity.z[i] += half_dt_inv_mass * state.force.z[i];
        }

        return energy;
    }
};

} // namespace atomweave::dynamics