#include <array>
#include <iostream>
#include "diagnostics/observables.hpp"
#include "state/particle_state.hpp"
#include "tests/test.hpp"


int main()
{
    using namespace atomweave;

    state::ParticleState<double> state(2);

    state.set_velocity(0, 1, 0, 0);
    state.set_velocity(1, 0, 2, 0);

    std::array<double, 3> total_momentum = diagnostics::total_momentum(state, 2.0);

    test::assert_near(total_momentum[0], 2.0, "Total momentum: x");
    test::assert_near(total_momentum[1], 4.0, "Total momentum: y");
    test::assert_near(total_momentum[2], 0.0, "Total momentum: z");

    test::assert_near(diagnostics::kinetic_energy(state, 2.0), 5.0,
                      "Kinetic energy: particle");

    std::cout << "All observable validation tests passed." << std::endl;

    return 0;
}