#include "solver/linear.h"

#include "utility/python_utility.h"
#include "control/diagnostic_tool/performance_measurement.h"
#include "partition/partitioned_petsc_mat/partitioned_petsc_mat.h"
#include "control/diagnostic_tool/memory_leak_finder.h"

namespace Solver
{

Linear::Linear(PythonConfig specificSettings, MPI_Comm mpiCommunicator, std::string name) :
  Solver(specificSettings, name)
{
  if (VLOG_IS_ON(1))
  {
    int size;
    PetscErrorCode ierr;
    ierr = MPI_Comm_size(mpiCommunicator, &size); CHKERRV(ierr);
    VLOG(1) << "Create linear solver on " << size << (size == 1? " rank." : " ranks.");
  }

  mpiCommunicator_ = mpiCommunicator;

  parseOptions();
}

void Linear::parseOptions()
{
  // parse options
  relativeTolerance_ = this->specificSettings_.getOptionDouble("relativeTolerance", 1e-5, PythonUtility::Positive);
  absoluteTolerance_ = this->specificSettings_.getOptionDouble("absoluteTolerance", 0, PythonUtility::NonNegative);  // 0 means disabled
  maxIterations_ = this->specificSettings_.getOptionDouble("maxIterations", 10000, PythonUtility::Positive);

  //parse information to use for dumping matrices and vectors
  dumpFormat_ = this->specificSettings_.getOptionString("dumpFormat", "default");
  dumpFilename_ = this->specificSettings_.getOptionString("dumpFilename", "");

  // prepare log keys to log number of iterations and residual norm
  std::stringstream nIterationsLogKey;
  nIterationsLogKey << "nIterations_" << name_;
  nIterationsLogKey_ = nIterationsLogKey.str();
  
  std::stringstream residualNormLogKey;
  residualNormLogKey << "residualNorm_" << name_;
  residualNormLogKey_ = residualNormLogKey.str();
  
  std::stringstream nIterationsTotalLogKey;
  nIterationsTotalLogKey << "nIterationsTotal_" << name_;
  nIterationsTotalLogKey_ = nIterationsTotalLogKey.str();
}

void Linear::initialize()
{
  // set up KSP object
  ksp_ = std::make_shared<KSP>();
  PetscErrorCode ierr = KSPCreate(mpiCommunicator_, ksp_.get()); CHKERRV(ierr);
  KSP ksp = *ksp_;

  setupKsp(ksp);
}

void Linear::setupKsp(KSP ksp)
{
  // parse the solver and preconditioner types from settings
  parseSolverTypes();

  // set solver type
  PetscErrorCode ierr;
  ierr = KSPSetType(ksp, kspType_); CHKERRV(ierr);

  // set options from command line, this overrides the python config
  ierr = KSPSetFromOptions(ksp); CHKERRV(ierr);

  // extract preconditioner context
  PC pc;
  ierr = KSPGetPC(ksp, &pc); CHKERRV(ierr);

  // set type of preconditioner
  ierr = PCSetType(pc, pcType_); CHKERRV(ierr);

  // for multigrid set number of levels and cycle type
  if (pcType_ == std::string(PCGAMG))
  {
    int nLevels = this->specificSettings_.getOptionInt("nLevels", 25, PythonUtility::Positive);
    ierr = PCMGSetLevels(pc, nLevels, NULL); CHKERRV(ierr);
    
    std::string mgType = this->specificSettings_.getOptionString("gamgType", "agg");
    
    PCGAMGType gamgType = PCGAMGCLASSICAL;
    
    if (mgType == "agg")
    {
      gamgType = PCGAMGAGG;
    }
    else if (mgType == "geo")
    {
      gamgType = PCGAMGGEO;
    }
    
    ierr = PCGAMGSetType(pc,gamgType); CHKERRV(ierr);
    
    mgType = this->specificSettings_.getOptionString("cycleType", "cycleV");
    
    PCMGCycleType cycleType = PC_MG_CYCLE_V;
    
    if (mgType == "cycleW")
    {
      cycleType = PC_MG_CYCLE_W;
    }
    
    ierr = PCMGSetCycleType(pc, cycleType); CHKERRV(ierr);    
  }
  // set Hypre Options from Python config
  else if (pcType_ == std::string(PCHYPRE))
  {
    std::string hypreOptions = this->specificSettings_.getOptionString("hypreOptions", "");
    PetscOptionsInsertString(NULL, hypreOptions.c_str());
    
    // if one of the hypre preconditioners is in preconditionerType_, pcType_ was set to HYPRE, now set the chosen preconditioner as -pc_hypre_type
    if (preconditionerType_ == "euclid" || preconditionerType_ == "pilut" || preconditionerType_ == "parasails" 
      || preconditionerType_ == "boomeramg" || preconditionerType_ == "ams" || preconditionerType_ == "ads")
    {
#if defined(PETSC_HAVE_HYPRE)
      ierr = PCHYPRESetType(pc, preconditionerType_.c_str()); CHKERRV(ierr);
#else
      LOG(ERROR) << "Petsc is not compiled with HYPRE!";
#endif
      
      LOG(DEBUG) << "set pc_hypre_type to " << preconditionerType_;
    }
  }

  // set options from command line, this overrides the python config
  ierr = PCSetFromOptions(pc); CHKERRV(ierr);

  //                                    relative tol,      absolute tol,  diverg tol.,   max_iterations
  ierr = KSPSetTolerances (ksp, relativeTolerance_, absoluteTolerance_, PETSC_DEFAULT, maxIterations_); CHKERRV(ierr);

}

void Linear::parseSolverTypes()
{
  // parse solver type
  solverType_ = this->specificSettings_.getOptionString("solverType", "gmres");

  // parse preconditioner type
  preconditionerType_ = this->specificSettings_.getOptionString("preconditionerType", "none");
  

  // all pc types: https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html
  pcType_ = PCNONE;
  if (preconditionerType_ == "jacobi")
  {
    pcType_ = PCJACOBI;
  }
  else if (preconditionerType_ == "sor")
  {
    pcType_ = PCSOR;
  }
  else if (preconditionerType_ == "lu")
  {
    pcType_ = PCLU;
  }
  else if (preconditionerType_ == "ilu")
  {
    pcType_ = PCILU;
  }
  else if (preconditionerType_ == "gamg")
  {
    pcType_ = PCGAMG;
  }
  else if (preconditionerType_ == "pcmg")
  {
    pcType_ = PCMG;
  }
  // the hypre boomeramg as the only solver does not provide the correct solution 
  else if (preconditionerType_ == "pchypre")
  {
    pcType_ = PCHYPRE;
  }
  // if one of the hypre preconditioners is in preconditionerType_, set pcType_ to HYPRE and set the chosen preconditioner as -pc_hypre_type
  else if (preconditionerType_ == "euclid" || preconditionerType_ == "pilut" || preconditionerType_ == "parasails" 
    || preconditionerType_ == "boomeramg" || preconditionerType_ == "ams" || preconditionerType_ == "ads")
  {
    pcType_ = PCHYPRE;
  }
  else if (preconditionerType_ == "none")
  {
    pcType_ = PCNONE;
  }
  else if (preconditionerType_ != "none" && preconditionerType_ != "")
  {
    pcType_ = preconditionerType_.c_str();
  }
  
  // all ksp types: https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPType.html#KSPType
  kspType_ = KSPGMRES;
  if (solverType_ == "richardson")
  {
    kspType_ = KSPRICHARDSON ;
  }
  else if (solverType_ == "chebyshev")
  {
    kspType_ = KSPCHEBYSHEV;
  }
  else if (solverType_ == "cg")
  {
    kspType_ = KSPCG;
  }
  else if (solverType_ == "bcgs")
  {
    kspType_ = KSPBCGS;
  }
  else if (solverType_ == "preonly")
  {
    kspType_ = KSPPREONLY;
  }
  else if (solverType_ == "lu")
  {
    kspType_ = KSPPREONLY;
    pcType_ = PCLU;
  }
  else if (solverType_ == "cholesky")
  {
    kspType_ = KSPPREONLY;
    pcType_ = PCCHOLESKY;
  }
  else if (solverType_ == "gamg")
  {
    kspType_ = KSPPREONLY;
    pcType_ = PCGAMG;
  }
  else if (solverType_ == "jacobi")
  {
    kspType_ = KSPPREONLY;
    pcType_ = PCJACOBI;
  }
  else if (solverType_ == "sor")
  {
    kspType_ = KSPPREONLY;
    pcType_ = PCSOR;
  }
  else if (solverType_ == "gmres")
  {
    kspType_ = KSPGMRES;
  }
  else if (solverType_ != "")
  {
    kspType_ = solverType_.c_str();
  }

  std::stringstream optionKey;
  optionKey << this->name_ << "_solverType";
  Control::PerformanceMeasurement::setParameter(optionKey.str(), solverType_);

  optionKey.str("");
  optionKey << this->name_ << "_preconditionerType";
  Control::PerformanceMeasurement::setParameter(optionKey.str(), preconditionerType_);

  LOG(DEBUG) << "linear solver type: " << solverType_ << " (" << kspType_ << "), preconditionerType_: " << preconditionerType_ << " (" << pcType_ << ")";
}

std::shared_ptr<KSP> Linear::ksp()
{
  return ksp_;
}

void Linear::dumpMatrixRightHandSideSolution(Vec rightHandSide, Vec solution)
{
  // dump files containing rhs and system matrix
  if (!dumpFilename_.empty())
  {
    PetscUtility::dumpVector(dumpFilename_+std::string("_rhs"), dumpFormat_, rightHandSide, mpiCommunicator_);
    PetscUtility::dumpVector(dumpFilename_+std::string("_solution"), dumpFormat_, solution, mpiCommunicator_);

    // get matrix
    Mat matrix;
    Mat preconditionerMatrix;
    KSPGetOperators(*ksp_, &matrix, &preconditionerMatrix);
    
    PetscUtility::dumpMatrix(dumpFilename_+std::string("_matrix"), dumpFormat_, matrix, mpiCommunicator_);
    if (matrix != preconditionerMatrix)
    {
      PetscUtility::dumpMatrix(dumpFilename_+std::string("_preconditioner_matrix"), dumpFormat_, matrix, mpiCommunicator_);
    }
  }
}

bool Linear::solve(Vec rightHandSide, Vec solution, std::string message)
{
  PetscErrorCode ierr;

  Control::PerformanceMeasurement::start(this->durationLogKey_);

  // reset memory count in MemoryLeakFinder
  Control::MemoryLeakFinder::nBytesIncreaseSinceLastCheck();

  // solve the system
  ierr = KSPSolve(*ksp_, rightHandSide, solution); CHKERRQ(ierr);

  // output a warning if the memory increased by over 1 MB
  Control::MemoryLeakFinder::warnIfMemoryConsumptionIncreases("In Linear::solve, after KSPSolve");
    
  Control::PerformanceMeasurement::stop(this->durationLogKey_);

  // dump files of rhs, solution and system matrix for debugging
  dumpMatrixRightHandSideSolution(rightHandSide, solution);

  // determine meta data
  PetscReal residualNorm = 0.0;
  PetscInt nDofsGlobal = 0;

  ierr = KSPGetIterationNumber(*ksp_, &lastNumberOfIterations_); CHKERRQ(ierr);
  ierr = KSPGetResidualNorm(*ksp_, &residualNorm); CHKERRQ(ierr);

  KSPConvergedReason convergedReason;
  ierr = KSPGetConvergedReason(*ksp_, &convergedReason); CHKERRQ(ierr);
  ierr = VecGetSize(rightHandSide, &nDofsGlobal); CHKERRQ(ierr);

  // compute residual norm
  if (kspType_ == KSPPREONLY && (pcType_ == PCLU || pcType_ == PCILU))
  {
    if (!residual_)
    {
      temporaryVectorLeft_ = std::make_shared<Vec>();     //< temporary vector for computation of residual for direct solvers
      temporaryVectorRight_ = std::make_shared<Vec>();    //< temporary vector for computation of residual for direct solvers
      residual_ = std::make_shared<Vec>();    //< residual vector for direct solvers

      Mat systemMatrix;
      ierr = KSPGetOperators(*ksp_, &systemMatrix, NULL); CHKERRQ(ierr);
      ierr = MatCreateVecs(systemMatrix, &(*temporaryVectorRight_), &(*temporaryVectorLeft_)); CHKERRQ(ierr);

      LOG(DEBUG) << "create temporary vectors";
    }
    LOG(DEBUG) << "compute residual";

    // compute residual
    ierr = KSPBuildResidual(*ksp_, *temporaryVectorLeft_, *temporaryVectorRight_, &(*residual_)); CHKERRQ(ierr);

    // compute norm of residual
    ierr = VecNorm(*residual_, NORM_2, &residualNorm); CHKERRQ(ierr);
  }
  
  // output message
  if (message != "")
  {
    // example for output: "Linear system of multidomain problem solved in 373 iterations, 3633 dofs, residual norm 9.471e-11: KSP_CONVERGED_ATOL: residual 2-norm less than abstol"
    LOG(INFO) << message << " in " << lastNumberOfIterations_ << " iterations, " << nDofsGlobal << " dofs, residual norm " << residualNorm
      << ": " << PetscUtility::getStringLinearConvergedReason(convergedReason);
  }

  // store parameter values to be logged
  Control::PerformanceMeasurement::setParameter(nIterationsLogKey_, lastNumberOfIterations_);
  Control::PerformanceMeasurement::setParameter(residualNormLogKey_, residualNorm);
  Control::PerformanceMeasurement::countNumber(nIterationsTotalLogKey_, lastNumberOfIterations_);

  // if convergedReason > 0 then it converged
  return convergedReason > 0;
}

int Linear::lastNumberOfIterations()
{
  return lastNumberOfIterations_;
}

}   //namespace
