#pragma once

#include <Python.h>  // has to be the first included header
#include "solver/solver.h"

#include <petscksp.h>
#include <memory>

namespace Solver
{

/**
 * A linear solver configuration that contains all parameters for PETSc and the KSP object.
 */
class Linear : public Solver
{
public:
  //! construct solver from python settings
  Linear(PythonConfig specificSettings, MPI_Comm mpiCommunicator, std::string name);

  //! return the KSP object that is used for solving
  std::shared_ptr<KSP> ksp();

protected:

  //! parse the solver and preconditioner type from settings
  void parseSolverTypes(KSPType &kspType, PCType &pcType);

  std::shared_ptr<KSP> ksp_;   ///< the PETSc KSP (Krylov subspace) object
  double relativeTolerance_;    ///< relative solver tolerance
  long int maxIterations_;     ///< maximum number of iterations
};

class LinearMG : public Solver
{
public:
  //! construct solver from python settings
  LinearMG(PythonConfig specificSettings, MPI_Comm mpiCommunicator, std::string name);

  //! return the KSP object that is used for solving
  std::shared_ptr<KSP> ksp();
  std::shared_ptr<KSP> ksp2();
  std::shared_ptr<int> numCycles();
protected:

  //! parse the solver and preconditioner type from settings
  void parseSolverTypes(KSPType &kspType, PCType &pcType);
  
  std::shared_ptr<KSP> ksp_;   ///< the PETSc KSP (Krylov subspace) object
  std::shared_ptr<KSP> ksp2_;   ///< the PETSc KSP (Krylov subspace) object
  double relativeTolerance_;    ///< relative solver tolerance
  long int maxIterations_;     ///< maximum number of iterations
  std::shared_ptr<int> numCycles_;	///<number of multigrid cycles
};

}  // namespace
