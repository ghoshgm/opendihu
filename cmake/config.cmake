
include(CheckIncludeFile)
include(CheckSymbolExists)
include(ProcessorCount)

# Generic routine to set appropriate cores for all tests.
# Typically the tests are executed with maximum 2 cores.
cmake_host_system_information(RESULT Ncpu QUERY NUMBER_OF_PHYSICAL_CORES)
if(Ncpu LESS 2)
  ProcessorCount(n)
  if(n GREATER Ncpu)
    set(Ncpu ${n})
  endif()
  set(MPIEXEC_NUMPROC_MAX 1)
else()
  set(MPIEXEC_NUMPROC_MAX 2)
endif()