#pragma once

#include <vector>

#include "interfaces/runnable.h"
#include "interfaces/discretizable_in_time.h"
#include "control/dihu_context.h"
#include "output_writer/manager.h"
#include "function_space/function_space.h"
#include "basis_function/lagrange.h"
#include "cellml/01_rhs_routine_handler.h"

/** The is a class that contains cellml equations and can be used with a time stepping scheme.
 *  The nStates template parameter specifies the number of state variables that should be used with the integrator.
 *  It is necessary that this value is fixed at compile time because the timestepping scheme needs to know which field variable types is has to construct.
 *  This class can also be computed easily in multiple instances along the nodes of a mesh.
 * 
 *  Naming:
 *   Intermediate (opendihu) = KNOWN (OpenCMISS) = Algebraic (OpenCOR)
 *   Parameter (opendihu, OpenCMISS) = KNOWN (OpenCMISS), in OpenCOR also algebraic
 *   Constant - these are constants that are only present in the source files
 *   State: state variable
 *   Rate: the time derivative of the state variable, i.e. the increment value in an explicit Euler stepping
 */
template <int nStates, int nIntermediates_, typename FunctionSpaceType>
class CallbackHandler :
  public RhsRoutineHandler<nStates,nIntermediates_,FunctionSpaceType>,
  public DiscretizableInTime
{
public:

  //! constructor
  CallbackHandler(DihuContext context, bool initializeOutputWriter);

  //! constructor
  CallbackHandler(DihuContext context);

  //! destructor
  virtual ~CallbackHandler();

  //! register a callback function setParameters that can set parameter values before each computation
  void registerSetParameters(void (*setParameters) (void *context, int nInstances, int timeStepNo, double currentTime,
                                                    std::vector<double> &parameters));

  //! register a callback function setSpecificParameters that can set parameter values before each computation
  void registerSetSpecificParameters(void (*setSpecificParameters) (void *context, int nInstances, int timeStepNo, double currentTime,
                                                                    std::vector<double> &localParameters));

  //! register a callback function setSpecificStates that can set state values before each computation
  void registerSetSpecificStates(void (*setSpecificStatesParameters) (void *context, int nInstances, int timeStepNo, double currentTime,
                                                                      double *states));

  //! register a callbackfunction handleResult that gets called after each new values are available
  void registerHandleResult(void (*handleResult) (void *context, int nInstances, int timeStepNo, double currentTime,
                                                  double *states, double intermediates[]));
  
  //! directly call the python callback if it exists
  void callPythonSetParametersFunction(int nInstances, int timeStepNo, double currentTime, double *parameterValues, int nParameters);

  //! directly call the python callback if it exists
  void callPythonSetSpecificParametersFunction(int nInstances, int timeStepNo, double currentTime, double *localParameterValues, int nLocalParameters);

  //! directly call the python callback if it exists
  void callPythonSetSpecificStatesFunction(int nInstances, int timeStepNo, double currentTime, double *states);

  //! directly call the python callback if it exists
  void callPythonHandleResultFunction(int nInstances, int timeStepNo, double currentTime, double *states, double *intermediates);

  //! get the values of this->lastCallSpecificStatesTime
  double lastCallSpecificStatesTime();

  //! set the value of this->lastCallSpecificStatesTime
  void setLastCallSpecificStatesTime(double lastCallSpecificStatesTime);

protected:
 
  //! construct the python call back functions from config
  virtual void initializeCallbackFunctions();

  //! call Py_CLEAR on all python objects
  void clearPyObjects();

  void (*setParameters_) (void *context, int nInstances, int timeStepNo, double currentTime, double *parameterValues, int nParameterValues);  //< callback function that will be called before new states are computed. It can set new parameters ("known" variables) for the computation.
  void (*setSpecificParameters_) (void *context, int nInstances, int timeStepNo, double currentTime, double *localParameterValues, int nLocalParameters);  //< callback function that will be called before new states are computed. It can set values for global parameters ("known" variables) for the computation.
  void (*setSpecificStates_) (void *context, int nInstances, int timeStepNo, double currentTime, double *states);  //< callback function that will be called before new states are computed. It can set values for states.
  void (*handleResult_) (void *context, int nInstances, int timeStepNo, double currentTime, double *states, double *intermediates);   //< callback function that will be called after new states and intermediates were computed

  int setParametersCallInterval_;                 //< setParameters_ will be called every callInterval_ time steps
  int setSpecificParametersCallInterval_;         //< setSpecificParameters_ will be called every callInterval_ time steps
  int setSpecificStatesCallInterval_;             //< setSpecificStates_ will be called every callInterval_ time steps
  int handleResultCallInterval_;                  //< handleResult will be called every callInterval_ time steps
  double setSpecificStatesCallFrequency_;         //< frequency, after which the setSpecificStates callback function will be called, either this condition or the condition with setSpecificStatesCallInterval_ is used
  std::vector<double> setSpecificStatesFrequencyJitter_;   //< relative jitter values: factors of setSpecificStatesCallFrequency_, random jitter to add or substract from frequency
  double currentJitter_;                          //< the absolute value of the current jitter
  int jitterIndex_;                               //< which of the stored jitter values in setSpecificStatesFrequencyJitter_ to use

  double lastCallSpecificStatesTime_;             //< last time the setSpecificStates_ method was called
  double setSpecificStatesRepeatAfterFirstCall_;  //< duration of continuation of calling the setSpecificStates callback after it was triggered
  double setSpecificStatesCallEnableBegin_;       //< first time when setSpecificStates should be called
 
  PyObject *pythonSetParametersFunction_;         //< Python function handle that is called to set parameters to the CellML problem from the python config
  PyObject *pythonSetSpecificParametersFunction_; //< Python function handle that is called to set parameters to the CellML problem from the python config
  PyObject *pythonSetSpecificStatesFunction_;     //< Python function handle that is called to set states to the CellML problem from the python config
  PyObject *pythonHandleResultFunction_;          //< Python function handle that is called to process results from CellML problem from the python config

  PyObject *pySetFunctionAdditionalParameter_;    //< an additional python object that will be passed as last argument to the setParameters, setSpecificParameters and setSpecificStates callback function
  PyObject *pyHandleResultFunctionAdditionalParameter_;   //< an additional python object that will be passed as last argument to the handleResult callback function

  PyObject *pyGlobalNaturalDofsList_;             //< python list of global dof nos
};

#include "cellml/02_callback_handler.tpp"
