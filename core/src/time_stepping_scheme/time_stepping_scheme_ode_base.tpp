#include "time_stepping_scheme/time_stepping_scheme_ode.h"

#include <Python.h>  // has to be the first included header
#include <vector>

#include "utility/python_utility.h"

namespace TimeSteppingScheme
{

template<typename FunctionSpaceType, int nComponents>
TimeSteppingSchemeOdeBase<FunctionSpaceType,nComponents>::
TimeSteppingSchemeOdeBase(DihuContext context, std::string name) :
TimeSteppingScheme(context[name]), initialized_(false)
{
  // get python config
  this->specificSettings_ = this->context_.getPythonConfig();

  // initialize output writers
  this->outputWriterManager_.initialize(context_, this->specificSettings_);
}

template<typename FunctionSpaceType, int nComponents>
Data::TimeStepping<FunctionSpaceType, nComponents> &TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::
data()
{
  return *data_;
}

template<typename FunctionSpaceType, int nComponents>
void TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::
setInitialValues()
{
  // set initial values as given in settings, or set to zero if not given
  std::vector<double> localValues;

  bool inputMeshIsGlobal = this->specificSettings_.getOptionBool("inputMeshIsGlobal", true);
  if (inputMeshIsGlobal)
  {
    assert(this->data_);
    assert(this->data_->functionSpace());
    const int nDofsGlobal = this->data_->functionSpace()->nDofsGlobal();
    LOG(DEBUG) << "setInitialValues, nDofsGlobal = " << nDofsGlobal;

    this->specificSettings_.getOptionVector("initialValues", nDofsGlobal, localValues);

    this->data_->functionSpace()->meshPartition()->extractLocalDofsWithoutGhosts(localValues);
  }
  else
  {
    const int nDofsLocal = this->data_->functionSpace()->nDofsLocalWithoutGhosts();
    this->specificSettings_.getOptionVector("initialValues", nDofsLocal, localValues);
  }
  VLOG(1) << "set initial values to " << localValues;

  // set the first component of the solution variable by the given values
  data_->solution()->setValuesWithoutGhosts(0, localValues);

  VLOG(1) << data_->solution();
}

template<typename FunctionSpaceType, int nComponents>
typename TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::TransferableSolutionDataType TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::
getSolutionForTransfer()
{
  return data_->getSolutionForTransfer();
}

template<typename FunctionSpaceType, int nComponents>
void TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::
setRankSubset(Partition::RankSubset rankSubset)
{
  data_->setRankSubset(rankSubset);
}

template<typename FunctionSpaceType, int nComponents>
void TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::
reset()
{
  TimeSteppingScheme::reset();
  initialized_ = false;
}

template<typename FunctionSpaceType, int nComponents>
void TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::
initialize()
{
  if (initialized_)
    return;

  TimeSteppingScheme::initialize();
  LOG(TRACE) << "TimeSteppingSchemeOdeBase::initialize";

  initialized_ = true;
}

template<typename FunctionSpaceType, int nComponents>
void TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::
run()
{
  // initialize
  this->initialize();

  // do simulations
  this->advanceTimeSpan();
}

//! output the given data for debugging
template<typename FunctionSpaceType, int nComponents>
std::string TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::
getString(typename TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::TransferableSolutionDataType &data)
{
  return data_->getString(data);
}

} // namespace
