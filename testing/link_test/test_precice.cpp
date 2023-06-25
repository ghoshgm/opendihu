#include <iostream>
#include <precice/SolverInterface.hpp>

int main(int argc, char** argv)
{
  std::cout << precice::getVersionInformation() << std::endl;
  
  return EXIT_SUCCESS;
}
