#include "time_stepping_scheme.h"

#include "utility/python_utility.h"

namespace TimeSteppingScheme
{

TimeSteppingScheme::TimeSteppingScheme(DihuContext context) :
  context_(context), initialized_(false)
{
  specificSettings_ = NULL;   // needs to be set by deriving class
  isTimeStepWidthSignificant_ = false;
}

void TimeSteppingScheme::setTimeStepWidth(double timeStepWidth)
{
  double epsilon = 1e-15;
  numberTimeSteps_ = int(std::ceil((endTime_ - startTime_) / timeStepWidth - epsilon));
}

void TimeSteppingScheme::setNumberTimeSteps(int numberTimeSteps)
{
  //numberTimeSteps_ = numberTimeSteps;
  timeStepWidth_=(endTime_ - startTime_) /numberTimeSteps;
  //LOG(DEBUG) << "timeStepWidth_ in setNumberTimeSteps: " << timeStepWidth_;
}

void TimeSteppingScheme::setTimeSpan(double startTime, double endTime)
{
  startTime_ = startTime;
  endTime_ = endTime;

  if(isTimeStepWidthSignificant_)
  {
    setTimeStepWidth(timeStepWidth_);
    LOG(DEBUG) << "set number of time steps to "<<numberTimeSteps_<<" from timeStepWidth "<<timeStepWidth_;
  }
}

void TimeSteppingScheme::initialize()
{
  if (initialized_)
    return;
  
  // initialize time stepping values
  startTime_ = 0.0;
  endTime_ = 1.0;
  if (PythonUtility::hasKey(specificSettings_, "endTime"))
    endTime_ = PythonUtility::getOptionDouble(specificSettings_, "endTime", 1.0, PythonUtility::Positive);

  LOG(DEBUG) << "  TimeSteppingScheme::initialize read endTime="<<endTime_;

  if (PythonUtility::hasKey(specificSettings_, "timeStepWidth"))
  {
    timeStepWidth_ = PythonUtility::getOptionDouble(specificSettings_, "timeStepWidth", 0.001, PythonUtility::Positive);
    setTimeStepWidth(timeStepWidth_);

    LOG(DEBUG) << "  TimeSteppingScheme::initialize, timeStepWidth="
      <<PythonUtility::getOptionDouble(specificSettings_, "timeStepWidth", 0.001, PythonUtility::Positive)
      <<", compute numberTimeSteps="<<numberTimeSteps_;

    if (PythonUtility::hasKey(specificSettings_, "numberTimeSteps"))
    {
      numberTimeSteps_ = PythonUtility::getOptionInt(specificSettings_, "numberTimeSteps", 10, PythonUtility::Positive);      
      isTimeStepWidthSignificant_ = false;
      LOG(WARNING) << "Time step width will be overridden by number of time steps (" << numberTimeSteps_ << ")";
    }
    else
    {
      isTimeStepWidthSignificant_ = true;
    }
  }
  else
  {
    numberTimeSteps_ = PythonUtility::getOptionInt(specificSettings_, "numberTimeSteps", 10, PythonUtility::Positive);
    LOG(DEBUG) << "  TimeSteppingScheme::initialize, timeStepWidth not specified, read numberTimeSteps_="<<numberTimeSteps_;
    setNumberTimeSteps(numberTimeSteps_);
  }

  LOG(INFO) << "Time span: [" << startTime_ << "," << endTime_ << "], Number of time steps: " << numberTimeSteps_
    << ", time step width: " << timeStepWidth_;
    
  initialized_ = true;
}

};  // namespace