#pragma once

#include <array>
#include "state/particle_state.hpp"


namespace atomweave::diagnostics
{

template <typename TReal>
TReal kinetic_energy(const state::ParticleState<TReal> &state, TReal mass)
{
    TReal total = 0.0;
    for (std::size_t i = 0; i < state.size(); ++i)
    {
        total += (state.velocity.x[i] * state.velocity.x[i] +
                  state.velocity.y[i] * state.velocity.y[i] +
                  state.velocity.z[i] * state.velocity.z[i]);
    }
    return 0.5 * mass * total;
}

template <typename TReal>
std::array<TReal, 3> total_momentum(const state::ParticleState<TReal> &state, TReal mass)
{
    TReal total_x = 0.0;
    TReal total_y = 0.0;
    TReal total_z = 0.0;

    for (std::size_t i = 0; i < state.size(); ++i)
    {
        total_x += state.velocity.x[i];
        total_y += state.velocity.y[i];
        total_z += state.velocity.z[i];
    }

    std::array<TReal, 3> total_momentum = {mass * total_x, mass * total_y,
                                           mass * total_z};

    return total_momentum;
}

} // namespace atomweave::diagnostics