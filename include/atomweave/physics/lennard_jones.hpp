#pragma once

#include <cmath>


namespace atomweave::physics
{

template <typename TReal> struct LennardJones
{
    TReal sigma;
    TReal epsilon;
    TReal r_c;

  private:
    TReal cutoff_potential_;
    TReal cutoff_radial_force_;

    static TReal pow3(const TReal x)
    {
        return x * x * x;
    }

    static TReal pow6(const TReal x)
    {
        TReal x2 = x * x;
        return x2 * x2 * x2;
    }

    TReal unshifted_potential(const TReal r, const TReal q) const
    {
        return 4 * epsilon * (q * q - q);
    }

    TReal unshifted_radial_force(const TReal r, const TReal q) const
    {
        return (24 * epsilon / r) * (2 * q * q - q);
    }

  public:
    LennardJones(TReal sigma, TReal epsilon, TReal r_c)
        : sigma(sigma), epsilon(epsilon), r_c(r_c)
    {
        TReal q_c = pow6(sigma / r_c);
        cutoff_potential_ = unshifted_potential(r_c, q_c);
        cutoff_radial_force_ = unshifted_radial_force(r_c, q_c);
    }

    TReal potential(const TReal r) const
    {
        TReal q = pow6(sigma / r);

        return unshifted_potential(r, q) - cutoff_potential_ +
            (r - r_c) * cutoff_radial_force_;
    }

    TReal force_coefficient(TReal r_squared) const
    {
        if (r_squared >= cutoff_squared())
        {
            return 0;
        }
        TReal r = std::sqrt(r_squared);
        TReal q = pow6(sigma / r);

        return (unshifted_radial_force(r, q) - cutoff_radial_force_) / r;
    }

    TReal cutoff_squared() const
    {
        return r_c * r_c;
    }
};

} // namespace atomweave::physics