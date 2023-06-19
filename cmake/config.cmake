
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

if(${MPI_FOUND})
  set(CMAKE_REQUIRED_LIBRARIES MPI::MPI_C)
  check_symbol_exists(MPI_Put mpi.h USE_MPI_PUT)
  check_symbol_exists(MPI_Win_lock mpi.h USE_MPI_WIN_LOCK)
  check_symbol_exists(MPI_Win_free mpi.h USE_MPI_WIN_FREE)
  check_symbol_exists(MPI_Win_unlock mpi.h USE_MPI_WIN_UNLOCK)
  check_symbol_exists(MPI_Win_create mpi.h USE_MPI_WIN_CREATE)
  check_symbol_exists(MPI_Win_allocate_shared mpi.h USE_MPI_ALLOC)
endif()