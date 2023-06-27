
#include <omp.h>
#include <iostream>

int main(int argc, char** argv)
{
  omp_set_num_threads(2);
  int tnum = omp_get_thread_num();
  std::cout << tnum << std::endl;
  return EXIT_SUCCESS;
}