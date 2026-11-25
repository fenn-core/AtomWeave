
#include <cmath>
#include <iostream>
#include <string>
#include "domain/periodic_cube.hpp"
#include "kernel/scalar_force.hpp"
#include "physics/lennard_jones.hpp"
#include "state/particle_state.hpp"
#include "tests/test.hpp"


int main()
{
    using namespace atomweave;

    state::ParticleState<double> state(2);
    physics::LennardJones lennard_jones(1.0, 1.0, 2.5);
    domain::PeriodicCube periodic_cube(10.0);

    state.set_position(0, 0, 0, 0);
    state.set_position(1, 0, 0, 0);

    double r_vals[] = {0.9, 1.0, 1.5, 2.4, 2.5, 2.6};

    double expected_force[] = {
        -138.6986234717, -24.0389994775, +1.1190293536, +0.0127806915, 0, 0};

    double total_potential[] = {
        6.7148350083, 0.0748161073, -0.2650202257, -0.0006047567, 0, 0};

    kernel::ScalarForceKernel<double> kernel;
    double energy;

    for (std::size_t i = 0; i < std::size(r_vals); ++i)
    {
        state.position.x[1] = r_vals[i];
        energy = kernel.evaluate(state, periodic_cube, lennard_jones);

        test::assert_near(energy, total_potential[i],
                          "Pair r=" + std::to_string(r_vals[i]) + ": potential energy");

        test::assert_near(state.force.x[0], expected_force[i],
                          "Pair r=" + std::to_string(r_vals[i]) + ": x force");
    }

    // periodic boundary test

    state.set_position(0, 0.2, 0, 0);
    state.set_position(1, 9.3, 0, 0);

    double expected_dx = 0.9;
    double expected_force_0 = 138.6986234717;
    double expected_force_1 = -138.6986234717;
    double expected_total_energy = 6.7148350083;

    double dx = periodic_cube.minimum_image(state.position.x[0] - state.position.x[1]);

    energy = kernel.evaluate(state, periodic_cube, lennard_jones);

    test::assert_near(dx, expected_dx, "Periodic boundary: minimum-image displacement");
    test::assert_near(state.force.x[0], expected_force_0,
                      "Periodic boundary: particle 0 x force");
    test::assert_near(state.force.x[1], expected_force_1,
                      "Periodic boundary: particle 1 x force");
    test::assert_near(energy, expected_total_energy,
                      "Periodic boundary: potential energy");

    // newton's third law and repeated evaluation

    state::ParticleState<double> multi_state(4);

    multi_state.set_position(0, 0.2, 0.3, 0.4);
    multi_state.set_position(1, 1.7, 0.8, 0.4);
    multi_state.set_position(2, 0.8, 2.0, 0.8);
    multi_state.set_position(3, 1.3, 1.7, 2.0);

    double initial_energy = kernel.evaluate(multi_state, periodic_cube, lennard_jones);

    double sum_x = 0;
    double sum_y = 0;
    double sum_z = 0;

    for (std::size_t i = 0; i < multi_state.size(); ++i)
    {
        sum_x += multi_state.force.x[i];
        sum_y += multi_state.force.y[i];
        sum_z += multi_state.force.z[i];
    }

    test::assert_near(sum_x, 0.0, "Newton's third law: net x force");
    test::assert_near(sum_y, 0.0, "Newton's third law: net y force");
    test::assert_near(sum_z, 0.0, "Newton's third law: net z force");

    const auto force_x = multi_state.force.x;
    const auto force_y = multi_state.force.y;
    const auto force_z = multi_state.force.z;

    double repeated_energy = kernel.evaluate(multi_state, periodic_cube, lennard_jones);

    test::assert_near(repeated_energy, initial_energy,
                      "Repeated evaluation: potential energy");

    for (std::size_t i = 0; i < multi_state.size(); ++i)
    {
        test::assert_near(multi_state.force.x[i], force_x[i],
                          "Repeated evaluation: particle " + std::to_string(i) +
                              " x force");
        test::assert_near(multi_state.force.y[i], force_y[i],
                          "Repeated evaluation: particle " + std::to_string(i) +
                              " y force");
        test::assert_near(multi_state.force.z[i], force_z[i],
                          "Repeated evaluation: particle " + std::to_string(i) +
                              " z force");
    }

    // force potential derivative consistency

    double r = 1.5;
    double h = 1e-5;

    double numerical_force =
        -(lennard_jones.potential(r + h) - lennard_jones.potential(r - h)) / (2 * h);

    double analytical_force = r * lennard_jones.force_coefficient(r * r);

    test::assert_near(numerical_force, analytical_force,
                      "Force-potential derivative consistency", 1e-8, 0.0);

    std::cout << "All force validation tests passed.\n";

    return 0;
}
