#pragma once

#include <cmath>
#include <random>
#include <stdexcept>
#include "diagnostics/observables.hpp"
#include "state/particle_state.hpp"


namespace atomweave::init
{

template <typename TReal>
void cubic_lattice(state::ParticleState<TReal> &state, TReal box_length,
                   std::size_t particles_per_axis)
{
    if (particles_per_axis == 0 ||
        particles_per_axis > state.size() / particles_per_axis / particles_per_axis ||
        particles_per_axis * particles_per_axis * particles_per_axis != state.size())
    {
        throw std::invalid_argument("Invalid cubic lattice dimensions");
    }

    TReal dx = box_length / particles_per_axis;

    std::size_t idx = 0;
    for (std::size_t i = 0; i < particles_per_axis; ++i)
    {
        const auto x = (i + 0.5) * dx;
        for (std::size_t j = 0; j < particles_per_axis; ++j)
        {
            const auto y = (j + 0.5) * dx;
            for (std::size_t k = 0; k < particles_per_axis; ++k, ++idx)
            {
                state.position.x[idx] = x;
                state.position.y[idx] = y;
                state.position.z[idx] = (k + 0.5) * dx;
            }
        }
    }
}

template <typename TReal>
void gaussian_velocities(state::ParticleState<TReal> &state, TReal temperature,
                         TReal mass, unsigned int seed)
{
    if (state.size() < 2 || !std::isfinite(mass) || mass <= 0 ||
        !std::isfinite(temperature) || temperature <= 0)
    {
        throw std::invalid_argument("Invalid gaussian state");
    }

    TReal stddev = std::sqrt(temperature / mass);

    std::mt19937 gen(seed);
    std::normal_distribution<TReal> dist(0.0, stddev);

    TReal total_x = 0.0;
    TReal total_y = 0.0;
    TReal total_z = 0.0;

    for (std::size_t i = 0; i < state.size(); ++i)
    {
        state.velocity.x[i] = dist(gen);
        total_x += state.velocity.x[i];

        state.velocity.y[i] = dist(gen);
        total_y += state.velocity.y[i];

        state.velocity.z[i] = dist(gen);
        total_z += state.velocity.z[i];
    }

    TReal mean_x = total_x / state.size();
    TReal mean_y = total_y / state.size();
    TReal mean_z = total_z / state.size();

    for (std::size_t i = 0; i < state.size(); ++i)
    {
        state.velocity.x[i] -= mean_x;
        state.velocity.y[i] -= mean_y;
        state.velocity.z[i] -= mean_z;
    }

    TReal kinetic_energy = diagnostics::kinetic_energy(state, mass);

    if (!std::isfinite(kinetic_energy) || kinetic_energy <= 0)
    {
        throw std::runtime_error("Non-finite or non-positive kinetic energy");
    }

    TReal scaling_factor =
        std::sqrt(temperature / (2 * kinetic_energy / (3 * state.size() - 3)));

    for (std::size_t i = 0; i < state.size(); ++i)
    {
        state.velocity.x[i] *= scaling_factor;
        state.velocity.y[i] *= scaling_factor;
        state.velocity.z[i] *= scaling_factor;
    }
}


} // namespace atomweave::init
