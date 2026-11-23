#include <iomanip>
#include <iostream>
#include "domain/periodic_cube.hpp"
#include "kernel/scalar_force.hpp"
#include "physics/lennard_jones.hpp"
#include "state/particle_state.hpp"


int main()
{
    std::setprecision(15);
    using namespace atomweave;

    state::ParticleState<double> state(2);
    physics::LennardJones<double> lennard_jones(1.0, 1.0, 2.5);
    domain::PeriodicCube<double> periodic_cube(10.0);


    state.position.x[0] = 0;
    state.position.y[0] = 0;
    state.position.z[0] = 0;

    state.position.x[1] = 1.5;
    state.position.y[1] = 0;
    state.position.z[1] = 0;

    kernel::ScalarForceKernel<double> kernel;
    double energy;

    energy = kernel.evaluate(state, periodic_cube, lennard_jones);
    std::cout << state.force.x[0] << " | " << energy;

}