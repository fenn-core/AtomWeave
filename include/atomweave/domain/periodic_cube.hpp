#pragma once

#include <cmath>
#include <stdexcept>


namespace atomweave::domain
{

template <typename TReal> struct PeriodicCube
{

  private:
    TReal L;

  public:
    explicit PeriodicCube(TReal L) : L(L)
    {
        if (this->L <= 0)
        {
            throw std::invalid_argument("L must be greater than 0");
        }
    }

    PeriodicCube(const PeriodicCube &other) = default;


    TReal wrap_coordinate(TReal x) const
    {
        return (x - L * std::floor(x / L));
    }

    TReal minimum_image(TReal dx) const
    {
        return (dx - L * std::round(dx / L));
    }

    TReal side_length() const
    {
        return L;
    }

    TReal volume() const
    {
        return (L * L * L);
    }
};


} // namespace atomweave::domain