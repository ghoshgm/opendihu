#include "cellml/01_rhs_routine_handler.h"

#include <Python.h>  // has to be the first included header

#include <list>
#include <sstream>
#include <sys/stat.h>  // stat() to check if file exists

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "utility/string_utility.h"
#include "mesh/mesh_manager/mesh_manager.h"

#include <unistd.h>  //dlopen
#include <dlfcn.h>
#include <ctime>

// forward declaration
template <int nStates,int nIntermediates_,typename FunctionSpaceType>
class CellmlAdapter;

template<int nStates, int nIntermediates_, typename FunctionSpaceType>
void RhsRoutineHandler<nStates,nIntermediates_,FunctionSpaceType>::
initializeRhsRoutine()
{
  // if "libraryFilename" is specified: use the given library
  // else: compile library from source

  // output deprecation warnings for old settings keys
  if (this->specificSettings_.hasKey("forceRecompileRhs"))
  {
    LOG(WARNING) << "Option \"forceRecompileRhs\" is no longer used, simply specify \"sourceFilename\" to compile the function.";
  }
  if (this->specificSettings_.hasKey("useGivenLibrary"))
  {
    LOG(WARNING) << "Option \"useGivenLibrary\" is no longer used, simply specify \"libraryFilename\" to use an existing library.";
  }
  if (this->specificSettings_.hasKey("simdSourceFilename"))
  {
    LOG(WARNING) << "Option \"simdSourceFilename\" is no longer used, simply specify \"sourceFilename\" to create a simd version and compile the function.";
  }

  // determine library filename, create library if necessary
  std::string libraryFilename;
  if (this->specificSettings_.hasKey("libraryFilename"))
  {
    libraryFilename = this->specificSettings_.getOptionString("libraryFilename", "lib.so");
  }
  else
  {
    // load type
    optimizationType_ = this->specificSettings_.getOptionString("optimizationType", "vc");

    if (optimizationType_ != "simd" && optimizationType_ != "vc" && optimizationType_ != "openmp")
    {
      LOG(ERROR) << "Option \"optimizationType\" is \"" << optimizationType_ << "\" but valid values are \"simd\", \"vc\" or \"openmp\"."
       << " Now setting to \"vc\".";
      optimizationType_ = "vc";
    }

    // load compiler flags
    std::string compilerFlags = this->specificSettings_.getOptionString("compilerFlags", "-fPIC -finstrument-functions -ftree-vectorize -fopt-info-vec-optimized=vectorizer_optimized.log -shared ");

    // for GPU: -ta=host,tesla,time

    // compose compile command
    std::stringstream s;
#ifdef NDEBUG
    // release mode
    // other possible options
    // -fopt-info-vec-missed=vectorizer_missed.log
    // -fopt-info-vec-all=vectorizer_all.log
    s << C_COMPILER_COMMAND << " -O3 " << compilerFlags << " ";
#else
    // debug mode
    s << C_COMPILER_COMMAND << " -O0 -ggdb " << compilerFlags << " ";
#endif
    std::string compileCommandOptions = s.str();

    // compile source file to a library
    s.str();
    s << "lib/"+StringUtility::extractBasename(this->cellmlSourceCodeGenerator_->sourceFilename()) << "_" << optimizationType_
      << "_" << this->nInstances_ << ".so";
    libraryFilename = s.str();

    int rankNoWorldCommunicator = DihuContext::ownRankNoCommWorld();
    s.str("");
    s << "src/"+StringUtility::extractBasename(this->cellmlSourceCodeGenerator_->sourceFilename()) << "_" << optimizationType_
      << "_" << this->nInstances_ << "." << rankNoWorldCommunicator << ".c";
    sourceToCompileFilename_ = s.str();

    // create path of library filename if it does not exist
    if (libraryFilename.find("/") != std::string::npos)
    {
      std::string path = libraryFilename.substr(0, libraryFilename.rfind("/"));
      // if directory does not yet exist, create it
      struct stat info;
      if (stat(path.c_str(), &info) != 0)
      {
        int ret = system((std::string("mkdir -p ")+path).c_str());

        if (ret != 0)
        {
          LOG(ERROR) << "Could not create path \"" << path << "\".";
        }
      }
    }

    // gather what number of instances all ranks have
    int nRanksCommunicator = this->functionSpace_->meshPartition()->nRanks();
    int ownRankNoCommunicator = this->functionSpace_->meshPartition()->ownRankNo();
    std::vector<int> nInstancesRanks(nRanksCommunicator);
    nInstancesRanks[ownRankNoCommunicator] = this->nInstances_;

    //std::cout << "ownRankNoCommunicator: " << ownRankNoCommunicator << ", Communicator has " << nRanksCommunicator << " ranks, nInstancesRanks: " << nInstancesRanks << std::endl;

    MPIUtility::handleReturnValue(MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, nInstancesRanks.data(),
                                                1, MPI_INT, this->functionSpace_->meshPartition()->mpiCommunicator()), "MPI_Allgather");

    // check if the library already exists by a previous compilation
    struct stat buffer;
    if (stat(libraryFilename.c_str(), &buffer) == 0)
    {
      LOG(DEBUG) << "Library \"" << libraryFilename << "\" already exists.";
    }
    else
    {
      // compile the library on only one rank
      createLibraryOnOneRank(libraryFilename, compileCommandOptions, nInstancesRanks);
    }

    // barrier to wait until the one rank that compiles the library has finished
    MPIUtility::handleReturnValue(MPI_Barrier(this->functionSpace_->meshPartition()->mpiCommunicator()), "MPI_Barrier");
  }

  loadRhsLibrary(libraryFilename);
}

template<int nStates, int nIntermediates_, typename FunctionSpaceType>
bool RhsRoutineHandler<nStates,nIntermediates_,FunctionSpaceType>::
loadRhsLibrary(std::string libraryFilename)
{
  // load dynamic library
  //
  std::string currentWorkingDirectory = getcwd(NULL,0);

  // append slash if there is none at the end
  if (currentWorkingDirectory[currentWorkingDirectory.length()-1] != '/')
    currentWorkingDirectory += "/";

  void *handle = NULL;
  for (int i = 0; handle == NULL && i < 50; i++)  // wait maximum 2.5 ms for rank 0 to finish
  {
    // load dynamic library
    handle = dlopen((currentWorkingDirectory+libraryFilename).c_str(), RTLD_LOCAL | RTLD_LAZY);
    if (i > 30)
    {
      std::this_thread::yield();
    }
    else if (i > 35)
    {
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    else if (i > 40)
    {
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    else if (i > 45)
    {
      std::this_thread::sleep_for(std::chrono::milliseconds(500));
    }
  }

  if (handle)
  {
    // load rhs method

    // try to load several routine names
    rhsRoutine_ = (void (*)(void *,double,double*,double*,double*,double*)) dlsym(handle, "computeCellMLRightHandSide");
    rhsRoutineGPU_ = (void (*)(void *,double,double*,double*,double*,double*)) dlsym(handle, "computeGPUCellMLRightHandSide");
    initConstsOpenCOR_ = (void(*)(double*, double*, double*)) dlsym(handle, "initConsts");
    computeRatesOpenCOR_ = (void(*)(double, double*, double*, double*, double*)) dlsym(handle, "computeRates");
    computeVariablesOpenCOR_  = (void(*)(double, double*, double*, double*, double*)) dlsym(handle, "computeVariables");

    LOG(DEBUG) << "Library \"" << libraryFilename << "\" loaded. "
      << "rhsRoutine_: " << (rhsRoutine_==NULL? "NULL" : "yes")
      << ", rhsRoutineGPU_: " << (rhsRoutineGPU_==NULL? "NULL" : "yes")
      << ", initConstsOpenCOR_: " << (initConstsOpenCOR_==NULL? "NULL" : "yes")
      << ", computeRatesOpenCOR_: " << (computeRatesOpenCOR_==NULL? "NULL" : "yes")
      << ", computeVariablesOpenCOR_: " << (computeVariablesOpenCOR_==NULL? "NULL" : "yes");

    if (rhsRoutine_)
    {
      // if the opendihu-generated rhs function (with openmp pragmas) is present in the library, we can directly use it and we are done in this method.
      return true;
    }
    else if (rhsRoutineGPU_)
    {
      // if the gpu-processed rhs function (with openmp 4.5 pragmas) is present in the library, we can directly use it and we are done in this method.
      rhsRoutine_ = rhsRoutineGPU_;
      return true;
    }
    else if (initConstsOpenCOR_)
    {
      LOG(FATAL) << "Could not load rhs routine from dynamic library \"" << libraryFilename << "\", but a computeRates method from OpenCOR was detected. "
        << "Please provide the sourceFilename of the source instead of the library.";
      return false;
    }
    else
    {
      LOG(FATAL) << "Could not load rhs routine from dynamic library \"" << libraryFilename << "\".";
      return false;
    }
  }
  else
  {
    LOG(FATAL) << "Could not load dynamic library \"" << libraryFilename << "\". Reason: " << dlerror();
  }
  return false;
}

template<int nStates, int nIntermediates_, typename FunctionSpaceType>
void RhsRoutineHandler<nStates,nIntermediates_,FunctionSpaceType>::
createLibraryOnOneRank(std::string libraryFilename, std::string compileCommandOptions, const std::vector<int> &nInstancesRanks)
{
  // get the global rank no, needed for the output filenames
  int rankNoWorldCommunicator = DihuContext::ownRankNoCommWorld();

  // determine if this rank should do compilation, such that each nInstances is compiled only once, by the rank with lowest number
  int i = 0;
  int rankWhichCompilesLibrary = 0;
  for (std::vector<int>::const_iterator iter = nInstancesRanks.begin(); iter != nInstancesRanks.end(); iter++, i++)
  {
    if (*iter == this->nInstances_)
    {
      rankWhichCompilesLibrary = i;
      break;
    }
  }

  LOG(DEBUG) << "Library will be compiled on rank " << rankWhichCompilesLibrary;

  int ownRankNoCommunicator = this->functionSpace_->meshPartition()->ownRankNo();
  if (rankWhichCompilesLibrary == ownRankNoCommunicator)
  {
    LOG(DEBUG) << "compile on this rank";

    // create source file
    this->cellmlSourceCodeGenerator_->generateSourceFile(sourceToCompileFilename_, optimizationType_);

    // create library file
    if (libraryFilename.find("/") != std::string::npos)
    {
      std::string path = libraryFilename.substr(0, libraryFilename.rfind("/"));
      int ret = system((std::string("mkdir -p ")+path).c_str());

      if (ret != 0)
      {
        LOG(ERROR) << "Could not create path \"" << path << "\" for library file.";
      }
    }

    std::stringstream compileCommand;

    // compile library to filename with "*.rankNoWorldCommunicator", then wait (different wait times for ranks), then rename file to without "*.rankNoWorldCommunicator"
    compileCommand << compileCommandOptions
      << " -o " << libraryFilename << "." << rankNoWorldCommunicator << " " << sourceToCompileFilename_
      //<< " && sleep " << int((rankNoWorldCommunicator%100)/10+1)
      << " && mv " << libraryFilename << "." << rankNoWorldCommunicator << " " << libraryFilename;

    int ret = system(compileCommand.str().c_str());
    if (ret != 0)
    {
      LOG(ERROR) << "Compilation failed. Command: \"" << compileCommand.str() << "\".";
      libraryFilename = "";
    }
    else
    {
      LOG(DEBUG) << "Compilation successful. Command: \"" << compileCommand.str() << "\".";
    }
  }
  else
  {
    LOG(DEBUG) << "we are the wrong rank, do not compile library "
      << "wait until library has been compiled";
  }
}
