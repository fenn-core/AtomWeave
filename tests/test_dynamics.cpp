
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include "domain/periodic_cube.hpp"
#include "dynamics/velocity_verlet.hpp"
#include "kernel/scalar_force.hpp"
#include "physics/lennard_jones.hpp"
#include "state/particle_state.hpp"
#include "tests/test.hpp"


int main()
{
    using namespace atomweave;

    // free particle test

    {
        state::ParticleState<double> state(1);
        physics::LennardJones lennard_jones(1.0, 1.0, 2.5);
        domain::PeriodicCube periodic_cube(10.0);
        kernel::ScalarForceKernel<double> kernel;
        dynamics::VelocityVerlet<double, kernel::ScalarForceKernel<double>>
            velocity_verlet;

        state.set_position(0, 9.9, 0.1, 5.0);
        state.set_velocity(0, 1, -1, 0.5);

        kernel.evaluate(state, periodic_cube, lennard_jones);

        double energy =
            velocity_verlet.step(state, periodic_cube, lennard_jones, kernel, 0.25, 1.0);

        // position after one timestep

        test::assert_near(state.position.x[0], 0.15, "Free particle: x position");
        test::assert_near(state.position.y[0], 9.85, "Free particle: y position");
        test::assert_near(state.position.z[0], 5.125, "Free particle: z position");

        // velocity remains constant

        test::assert_near(state.velocity.x[0], 1.0, "Free particle: x velocity");
        test::assert_near(state.velocity.y[0], -1.0, "Free particle: y velocity");
        test::assert_near(state.velocity.z[0], 0.5, "Free particle: z velocity");

        // no interactions for a single particle

        test::assert_near(state.force.x[0], 0.0, "Free particle: x force");
        test::assert_near(state.force.y[0], 0.0, "Free particle: y force");
        test::assert_near(state.force.z[0], 0.0, "Free particle: z force");

        test::assert_near(energy, 0.0, "Free particle: potential energy");

        std::cout << "Free particle dynamics test passed.\n";
    }

    // two particle interaction test

    {
        state::ParticleState<double> state(2);
        physics::LennardJones lennard_jones(1.0, 1.0, 2.5);
        domain::PeriodicCube periodic_cube(10.0);
        kernel::ScalarForceKernel<double> kernel;
        dynamics::VelocityVerlet<double, kernel::ScalarForceKernel<double>>
            velocity_verlet;

        state.set_position(0, 0.0, 0.0, 0.0);
        state.set_position(1, 1.5, 0.0, 0.0);

        state.set_velocity(0, 0.0, 0.0, 0.0);
        state.set_velocity(1, 0.0, 0.0, 0.0);

        kernel.evaluate(state, periodic_cube, lennard_jones);

        double initial_force = state.force.x[0];
        double dt = 0.001;

        double energy =
            velocity_verlet.step(state, periodic_cube, lennard_jones, kernel, dt, 1.0);

        // positions after one timestep

        test::assert_near(state.position.x[0], 5.595146768e-7,
                          "Two particles: particle 0 position");
        test::assert_near(state.position.x[1], 1.5 - 5.595146768e-7,
                          "Two particles: particle 1 position");

        // velocity Verlet half-kicks

        double expected_velocity = 0.5 * dt * (initial_force + state.force.x[0]);

        test::assert_near(state.velocity.x[0], expected_velocity,
                          "Two particles: particle 0 velocity");
        test::assert_near(state.velocity.x[1], -expected_velocity,
                          "Two particles: particle 1 velocity");

        // momentum and force balance

        test::assert_near(state.velocity.x[0] + state.velocity.x[1], 0.0,
                          "Two particles: x momentum");
        test::assert_near(state.force.x[0] + state.force.x[1], 0.0,
                          "Two particles: x force balance");

        test::assert_near(state.velocity.y[0], 0.0,
                          "Two particles: particle 0 y velocity");
        test::assert_near(state.velocity.y[1], 0.0,
                          "Two particles: particle 1 y velocity");

        test::assert_near(state.velocity.z[0], 0.0,
                          "Two particles: particle 0 z velocity");
        test::assert_near(state.velocity.z[1], 0.0,
                          "Two particles: particle 1 z velocity");

        // returned potential energy corresponds to updated positions

        double separation = std::abs(state.position.x[1] - state.position.x[0]);

        test::assert_near(energy, lennard_jones.potential(separation),
                          "Two particles: updated potential energy");

        std::cout << "Two-particle dynamics test passed.\n";
    }

    // energy conservation test

    {
        state::ParticleState<double> state(2);
        physics::LennardJones lennard_jones(1.0, 1.0, 2.5);
        domain::PeriodicCube periodic_cube(10.0);
        kernel::ScalarForceKernel<double> kernel;
        dynamics::VelocityVerlet<double, kernel::ScalarForceKernel<double>>
            velocity_verlet;

        state.set_position(0, 4.0, 0.0, 0.0);
        state.set_position(1, 5.5, 0.0, 0.0);

        state.set_velocity(0, 0.0, 0.0, 0.0);
        state.set_velocity(1, 0.0, 0.0, 0.0);

        double initial_energy = kernel.evaluate(state, periodic_cube, lennard_jones);

        double dt = 0.001;
        double mass = 1.0;

        // total mechanical energy E = K + U

        for (int step = 0; step < 1000; ++step)
        {
            double potential_energy = velocity_verlet.step(
                state, periodic_cube, lennard_jones, kernel, dt, mass);

            double kinetic_energy = 0.0;

            for (std::size_t i = 0; i < state.size(); ++i)
            {
                double vx = state.velocity.x[i];
                double vy = state.velocity.y[i];
                double vz = state.velocity.z[i];

                kinetic_energy += 0.5 * mass * (vx * vx + vy * vy + vz * vz);
            }

            double total_energy = kinetic_energy + potential_energy;

            test::assert_near(total_energy, initial_energy,
                              "Energy conservation: step " + std::to_string(step + 1),
                              5e-5, 0.0);
        }

        // total momentum remains zero

        test::assert_near(state.velocity.x[0] + state.velocity.x[1], 0.0,
                          "Energy test: x momentum");

        std::cout << "Energy conservation test passed.\n";
    }

    std::cout << "All dynamics tests passed." << std::endl;

    return 0;
}
