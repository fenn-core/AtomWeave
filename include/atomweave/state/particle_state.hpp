#pragma once

#include <vector>


namespace atomweave::state
{

template <typename TReal> struct VectorSoA
{
    std::vector<TReal> x;
    std::vector<TReal> y;
    std::vector<TReal> z;
};

template <typename TReal> struct ParticleState
{
    VectorSoA<TReal> position;
    VectorSoA<TReal> velocity;
    VectorSoA<TReal> force;

    explicit ParticleState(std::size_t count)
    {
        position.x.resize(count);
        position.y.resize(count);
        position.z.resize(count);

        velocity.x.resize(count);
        velocity.y.resize(count);
        velocity.z.resize(count);

        force.x.resize(count);
        force.y.resize(count);
        force.z.resize(count);
    }

    std::size_t size() const
    {
        return position.x.size();
    }
};

} // namespace atomweave::state
