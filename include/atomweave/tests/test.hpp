#pragma once

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>


namespace atomweave::test
{

template <typename TReal>
void assert_near(TReal actual, TReal expected, const std::string &err_msg,
                 TReal atol = 1e-10, TReal rtol = 1e-10)
{
#ifndef NDEBUG

    std::cerr << std::setprecision(std::numeric_limits<TReal>::max_digits10);

    if ((!std::isfinite(actual)) || (!std::isfinite(expected)))
    {
        std::cerr << "[FAIL] " << err_msg << "\n"
                  << "Expected: " << expected << "\n"
                  << "Actual: " << actual << "\n"
                  << "Non-finite value present " << std::endl;

        std::abort();
    }

    if (std::abs(actual - expected) <= atol + std::abs(expected) * rtol)
    {
        return;
    }

    std::cerr << "[FAIL] " << err_msg << "\n"
              << "Expected: " << expected << "\n"
              << "Actual: " << actual << "\n"
              << "Absolute Error: " << std::abs(actual - expected) << "\n"
              << std::endl;

    std::abort();

#endif
}

} // namespace atomweave::test
