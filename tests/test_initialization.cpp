
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include "state/particle_state.hpp"
#include "initialization/init.hpp"
#include "diagnostics/observables.hpp"
#include "tests/test.hpp"


int main()
{
    using namespace atomweave;

    // temperature and momentum

    {
        double mass = 1.0;
        double temperature = 1.0;

        state::ParticleState<double> state(64);
        init::gaussian_velocities(state, temperature, mass, 42);

        double kinetic_energy = diagnostics::kinetic_energy(state, mass);
        std::array total_momentum = diagnostics::total_momentum(state, mass);

        test::assert_near(kinetic_energy, 94.5,
                          "Gaussian velocities: kinetic energy");

        test::assert_near(total_momentum[0], 0.0,
                          "Gaussian velocities: x momentum");
        test::assert_near(total_momentum[1], 0.0,
                          "Gaussian velocities: y momentum");
        test::assert_near(total_momentum[2], 0.0,
                          "Gaussian velocities: z momentum");

        for (std::size_t i = 0; i < state.size(); ++i)
        {
            assert(std::isfinite(state.velocity.x[i]) &&
                   "Gaussian velocities: non-finite x velocity");
            assert(std::isfinite(state.velocity.y[i]) &&
                   "Gaussian velocities: non-finite y velocity");
            assert(std::isfinite(state.velocity.z[i]) &&
                   "Gaussian velocities: non-finite z velocity");
        }

        std::cout << "Gaussian velocity initialization test passed.\n";
    }

    // non-unit mass and temperature

    {
        double mass = 2.0;
        double temperature = 1.5;

        state::ParticleState<double> state(32);
        init::gaussian_velocities(state, temperature, mass, 123);

        double kinetic_energy = diagnostics::kinetic_energy(state, mass);
        std::array total_momentum = diagnostics::total_momentum(state, mass);

        test::assert_near(kinetic_energy, 69.75,
                          "Non-unit parameters: kinetic energy");

        test::assert_near(total_momentum[0], 0.0,
                          "Non-unit parameters: x momentum");
        test::assert_near(total_momentum[1], 0.0,
                          "Non-unit parameters: y momentum");
        test::assert_near(total_momentum[2], 0.0,
                          "Non-unit parameters: z momentum");

        std::cout << "Non-unit mass and temperature test passed.\n";
    }

    // reproducibility

    {
        state::ParticleState<double> state_a(64);
        state::ParticleState<double> state_b(64);

        init::gaussian_velocities(state_a, 1.0, 1.0, 42);
        init::gaussian_velocities(state_b, 1.0, 1.0, 42);

        for (std::size_t i = 0; i < state_a.size(); ++i)
        {
            test::assert_near(state_a.velocity.x[i], state_b.velocity.x[i],
                              "Reproducibility: x velocity", 0.0, 0.0);
            test::assert_near(state_a.velocity.y[i], state_b.velocity.y[i],
                              "Reproducibility: y velocity", 0.0, 0.0);
            test::assert_near(state_a.velocity.z[i], state_b.velocity.z[i],
                              "Reproducibility: z velocity", 0.0, 0.0);
        }

        std::cout << "Velocity reproducibility test passed.\n";
    }

    // invalid inputs

    {
        auto throws_invalid_argument = [](std::size_t count,
                                          double temperature, double mass)
        {
            state::ParticleState<double> state(count);

            try
            {
                init::gaussian_velocities(state, temperature, mass, 42);
            }
            catch (const std::invalid_argument &)
            {
                return true;
            }

            return false;
        };

        double inf = std::numeric_limits<double>::infinity();
        double nan = std::numeric_limits<double>::quiet_NaN();

        assert(throws_invalid_argument(0, 1.0, 1.0) &&
               "Invalid input: zero particles");
        assert(throws_invalid_argument(1, 1.0, 1.0) &&
               "Invalid input: single particle");

        assert(throws_invalid_argument(64, 1.0, 0.0) &&
               "Invalid input: zero mass");
        assert(throws_invalid_argument(64, 1.0, -1.0) &&
               "Invalid input: negative mass");
        assert(throws_invalid_argument(64, 1.0, inf) &&
               "Invalid input: infinite mass");
        assert(throws_invalid_argument(64, 1.0, nan) &&
               "Invalid input: NaN mass");

        assert(throws_invalid_argument(64, 0.0, 1.0) &&
               "Invalid input: zero temperature");
        assert(throws_invalid_argument(64, -1.0, 1.0) &&
               "Invalid input: negative temperature");
        assert(throws_invalid_argument(64, inf, 1.0) &&
               "Invalid input: infinite temperature");
        assert(throws_invalid_argument(64, nan, 1.0) &&
               "Invalid input: NaN temperature");

        std::cout << "Invalid initialization input tests passed.\n";
    }

    std::cout << "All Gaussian initialization tests passed.\n";

    return 0;
}
