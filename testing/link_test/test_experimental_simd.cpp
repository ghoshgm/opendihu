
#include <experimental/simd>
#include <iostream>
#include <cstdint>
#include <Vc/Vc>

#define SIZE 10

namespace stdx = std::experimental;
 
int main()
{
    const stdx::native_simd<std::int64_t> a = 3;
    for (std::size_t i = 0; i != a.size(); ++i)
        std::cout << a[i] << ' ';
    std::cout << '\n';

    /*
    The following block of code is just to check if
    libvc and libstdsimd interfere with each other.
    */
    Vc::Memory<Vc::float_v, SIZE> x_mem;

    for(size_t ii = 0; ii < x_mem.vectorsCount(); ++ii)
    {
      x_mem.vector(ii) = Vc::float_v::Random() * 2.f - 1.f;
    }

    for(size_t ii = 0; ii < x_mem.vectorsCount(); ++ii)
    {
      std::cout << x_mem[ii] << std::endl;
    }
}