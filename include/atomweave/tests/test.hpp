#pragma once

#include <cassert>
#include <cmath>


namespace atomweave::test
{

template <typename TReal>
void assert_near(TReal actual, TReal expected, TReal atol = 1e-10, TReal rtol = 1e-10)
{
    assert(std::abs(actual - expected) <= atol + std::abs(expected) * rtol);
}

} // namespace atomweave::test
