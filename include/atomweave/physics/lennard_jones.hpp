#pragma once


namespace atomweave::physics
{

template <typename TReal> struct LennardJones
{
    TReal sigma;
    TReal epsilon;
    TReal r_c;

  private:
    static TReal pow3(const TReal x)
    {
        return x * x * x;
    }

    static TReal pow6(const TReal x)
    {
        TReal x2 = x * x;
        return x2 * x2 * x2;
    }

  public:
    LennardJones(TReal sigma, TReal epsilon, TReal r_c)
        : sigma(sigma), epsilon(epsilon), r_c(r_c)
    {
    }

    TReal potential(const TReal r) const
    {
        TReal sigma_r = sigma / r;
        TReal sr6 = pow6(sigma_r);

        return (4 * epsilon * (sr6 * sr6 - sr6));
    }

    TReal force_coefficient(TReal r_squared) const
    {
        TReal ratio_squared = (sigma * sigma) / (r_squared);

        return (24 * epsilon / r_squared) * ((2 * pow6(ratio_squared)) - (pow3(ratio_squared)));
    }

    TReal cutoff_squared() const
    {
        return r_c * r_c;
    }

};

} // namespace atomweave::physics