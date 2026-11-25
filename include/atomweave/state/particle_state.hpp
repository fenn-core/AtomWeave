#pragma once

#include <cassert>
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

  private:
    void assert_consistent_size() const
    {
#ifndef NDEBUG
        const auto n = position.x.size();

        assert(position.y.size() == n);
        assert(position.z.size() == n);

        assert(velocity.x.size() == n);
        assert(velocity.y.size() == n);
        assert(velocity.z.size() == n);

        assert(force.x.size() == n);
        assert(force.y.size() == n);
        assert(force.z.size() == n);
#endif
    }

  public:
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

        assert_consistent_size();
    }

    std::size_t size() const
    {


        return position.x.size();
    }

    void set_position(std::size_t index, TReal x, TReal y, TReal z)
    {
        assert(index < position.x.size());
        assert_consistent_size();

        position.x[index] = x;
        position.y[index] = y;
        position.z[index] = z;
    }
};

} // namespace atomweave::state
