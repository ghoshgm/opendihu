
#include <Vc/Vc>
#include <iostream>

#define SIZE 10

int main()
{
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