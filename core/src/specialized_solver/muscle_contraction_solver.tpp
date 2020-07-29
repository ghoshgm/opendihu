#include "specialized_solver/muscle_contraction_solver.h"

#include <omp.h>
#include <sstream>

#include "utility/math_utility.h"
#include "control/diagnostic_tool/solver_structure_visualizer.h"

template<typename MeshType>
MuscleContractionSolver<MeshType>::
MuscleContractionSolver(DihuContext context) :
  Runnable(),
  ::TimeSteppingScheme::TimeSteppingScheme(context["MuscleContractionSolver"]),
  data_(this->context_), initialized_(false)
{
  // get python settings object from context
  this->specificSettings_ = this->context_.getPythonConfig();

  // determine if the dynamic or the quasi-static formulation is used
  isDynamic_ = this->specificSettings_.getOptionBool("dynamic", true);

  if (isDynamic_)
  {
    LOG(DEBUG) << "use dynamic hyperelasticity solver";
    dynamicHyperelasticitySolver_ = std::make_shared<DynamicHyperelasticitySolverType>(this->context_);
  }
  else
  {
    LOG(DEBUG) << "use static hyperelasticity solver";
    staticHyperelasticitySolver_ = std::make_shared<StaticHyperelasticitySolverType>(this->context_);
  }

  // initialize output writers
  this->outputWriterManager_.initialize(this->context_, this->specificSettings_);

  // parse options
  pmax_ = this->specificSettings_.getOptionDouble("Pmax", 1.0, PythonUtility::Positive);

  this->specificSettings_.template getOptionVector<std::string>("mapGeometryToMeshes", meshNamesOfGeometryToMapTo_);
}

template<typename MeshType>
void MuscleContractionSolver<MeshType>::
advanceTimeSpan()
{
  // This method computes some time steps of the simulation by running a for loop over the time steps.
  // The number of steps, timestep width and current time are all set by the parent class, TimeSteppingScheme.

  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  // compute time span of this method
  double timeSpan = this->endTime_ - this->startTime_;

  LOG_N_TIMES(3,INFO) << "durationComputeMappingBetweenMeshes: " << Control::PerformanceMeasurement::getDuration("durationComputeMappingBetweenMeshes");

  // output for debugging
  LOG(DEBUG) << "MuscleContractionSolver<MeshType>::advanceTimeSpan, timeSpan=" << timeSpan<< ", timeStepWidth=" << this->timeStepWidth_
    << " n steps: " << this->numberTimeSteps_;

  // loop over time steps
  double currentTime = this->startTime_;
  for (int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    // in defined intervals (settings "timeStepOutputInterval") print out the current timestep
    if (timeStepNo % this->timeStepOutputInterval_ == 0 && (this->timeStepOutputInterval_ <= 10 || timeStepNo > 0))  // show first timestep only if timeStepOutputInterval is <= 10
    {
      LOG(INFO) << "MuscleContractionSolver, timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;
    }

    // compute the current active stress
    computeActiveStress();

    if (isDynamic_)
    {
      this->dynamicHyperelasticitySolver_->setTimeSpan(currentTime, currentTime+this->timeStepWidth_);

      LOG(DEBUG) << "call dynamic hyperelasticitySolver";

      // advance the simulation by the specified time span
      dynamicHyperelasticitySolver_->advanceTimeSpan();
    }
    else
    {
      staticHyperelasticitySolver_->run();
    }

    // compute new values of λ and λ_dot, to be transferred to the CellML solvers
    computeLambda();

    // advance simulation time
    timeStepNo++;

    // compute new current simulation time
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;

    // stop duration measurement
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::stop(this->durationLogKey_);

    // write current output values using the output writers
    this->outputWriterManager_.writeOutput(this->data_, timeStepNo, currentTime);

    // start duration measurement
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::start(this->durationLogKey_);
  }

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);

  mapGeometryToGivenMeshes();
}

template<typename MeshType>
void MuscleContractionSolver<MeshType>::
initialize()
{
  // only initialize once
  if (initialized_)
    return;

  // initialize() will be called before the simulation starts.

  // call initialize of the parent class, this parses the timestepping settings from the settings file
  TimeSteppingScheme::TimeSteppingScheme::initialize();

  // add this solver to the solvers diagram
  DihuContext::solverStructureVisualizer()->addSolver("MuscleContractionSolver");

  // indicate in solverStructureVisualizer that now a child solver will be initialized
  DihuContext::solverStructureVisualizer()->beginChild();

  // call initialize of the nested timestepping solver
  if (isDynamic_)
  {
    dynamicHyperelasticitySolver_->initialize();
  }
  else
  {
    staticHyperelasticitySolver_->initialize();
  }

  // indicate in solverStructureVisualizer that the child solver initialization is done
  DihuContext::solverStructureVisualizer()->endChild();

  // In order to initialize the data object and actuall create all variables, we first need to assign a function space to the data object.
  // A function space object of type FunctionSpace<MeshType,BasisFunctionType> (see "function_space/function_space.h")
  // is an object that stores the mesh (e.g., all nodes and elements) as well as the basis function (e.g. linear Lagrange basis functions).
  // The dynamicHyperelasticitySolver_ solver already created a function space that we should use. We already have a typedef "FunctionSpace" that is the class of dynamicHyperelasticitySolver_'s function space type.
  std::shared_ptr<FunctionSpace> functionSpace;
  if (isDynamic_)
    functionSpace = dynamicHyperelasticitySolver_->data().functionSpace();
  else
    functionSpace = staticHyperelasticitySolver_->data().functionSpace();

  // Pass the function space to the data object. data_ stores field variables.
  // It needs to know the number of nodes and degrees of freedom (dof) from the function space in order to create the vectors with the right size.
  data_.setFunctionSpace(functionSpace);

  // now call initialize, data will then create all variables (Petsc Vec's)
  data_.initialize();

  bool setGeometryFieldForTransfer = meshNamesOfGeometryToMapTo_.empty();

  if (isDynamic_)
  {
    typename DynamicHyperelasticitySolverType::HyperelasticitySolverType &hyperelasticitySolver = dynamicHyperelasticitySolver_->hyperelasticitySolver();

    // set field variables from dynamicHyperelasticitySolver in data_ such that they can be output by the output writer
    data_.setFieldVariables(dynamicHyperelasticitySolver_->data().displacements(),
                            dynamicHyperelasticitySolver_->data().velocities(),
                            hyperelasticitySolver.data().activePK2Stress(),
                            hyperelasticitySolver.data().pK2Stress(),
                            hyperelasticitySolver.data().fiberDirection(),
                            hyperelasticitySolver.data().materialTraction(),
                            setGeometryFieldForTransfer);
  }
  else
  {
    // set field variables from dynamicHyperelasticitySolver in data_ such that they can be output by the output writer
    data_.setFieldVariables(staticHyperelasticitySolver_->data().displacements(),
                            staticHyperelasticitySolver_->data().velocities(),
                            staticHyperelasticitySolver_->data().activePK2Stress(),
                            staticHyperelasticitySolver_->data().pK2Stress(),
                            staticHyperelasticitySolver_->data().fiberDirection(),
                            staticHyperelasticitySolver_->data().materialTraction(),
                            setGeometryFieldForTransfer);
  }

  // set the slotConnectorData for the solverStructureVisualizer to appear in the solver diagram
  DihuContext::solverStructureVisualizer()->setSlotConnectorData(getSlotConnectorData());

  initialized_ = true;
}

template<typename MeshType>
void MuscleContractionSolver<MeshType>::
run()
{
  initialize();

  advanceTimeSpan();
}

template<typename MeshType>
void MuscleContractionSolver<MeshType>::
reset()
{
  if (isDynamic_)
    dynamicHyperelasticitySolver_->reset();
  else
    staticHyperelasticitySolver_->reset();

  // "uninitialize" everything
}

template<typename MeshType>
void MuscleContractionSolver<MeshType>::
computeLambda()
{
  typedef typename DynamicHyperelasticitySolverType::HyperelasticitySolverType::DisplacementsFieldVariableType DisplacementsFieldVariableType;
  typedef typename DynamicHyperelasticitySolverType::HyperelasticitySolverType::Data::DeformationGradientFieldVariableType DeformationGradientFieldVariableType;
  typedef typename Data::ScalarFieldVariableType FieldVariableType;

  std::shared_ptr<DisplacementsFieldVariableType> fiberDirectionVariable;
  std::shared_ptr<DeformationGradientFieldVariableType> deformationGradientVariable;
  std::shared_ptr<DeformationGradientFieldVariableType> fDotVariable;
  std::shared_ptr<DisplacementsFieldVariableType> velocitiesVariable;

  if (isDynamic_)
  {
    fiberDirectionVariable = dynamicHyperelasticitySolver_->hyperelasticitySolver().data().fiberDirection();
    deformationGradientVariable = dynamicHyperelasticitySolver_->hyperelasticitySolver().data().deformationGradient();
    fDotVariable = dynamicHyperelasticitySolver_->hyperelasticitySolver().data().deformationGradientTimeDerivative();
    velocitiesVariable = dynamicHyperelasticitySolver_->data().velocities();
  }
  else
  {
    fiberDirectionVariable = staticHyperelasticitySolver_->data().fiberDirection();
    deformationGradientVariable = staticHyperelasticitySolver_->data().deformationGradient();
    fDotVariable = staticHyperelasticitySolver_->data().deformationGradientTimeDerivative();
    velocitiesVariable = staticHyperelasticitySolver_->data().velocities();
  }

  // compute lambda and \dot{lambda} (contraction velocity)
  //std::shared_ptr<DisplacementsFieldVariableType> displacementsVariable = dynamicHyperelasticitySolver_->data().displacements();

  std::shared_ptr<FieldVariableType> lambdaVariable = data_.lambda();
  std::shared_ptr<FieldVariableType> lambdaDotVariable = data_.lambdaDot();

  // loop over local degrees of freedom
  for (dof_no_t dofNoLocal = 0; dofNoLocal < data_.functionSpace()->nDofsLocalWithoutGhosts(); dofNoLocal++)
  {
    const Vec3 fiberDirection = fiberDirectionVariable->getValue(dofNoLocal);
    //const Vec3 displacement = displacementsVariable->getValue(dofNoLocal);
    const VecD<9> deformationGradientValues = deformationGradientVariable->getValue(dofNoLocal);
    const VecD<9> fDotValues = fDotVariable->getValue(dofNoLocal);

    // create matrix, deformationGradientValues are in row-major order
    MathUtility::Matrix<3,3> deformationGradient(deformationGradientValues);
    MathUtility::Matrix<3,3> fDot(fDotValues);

    // fiberDirection is normalized
    assert (MathUtility::norm<3>(fiberDirection) - 1.0 < 1e-10);

    // get deformation gradient, project lambda and lambda dot
    // dx = F dX, dx^2 = C dX^2
    // lambda = ||dx•a0||
    // project displacements on normalized fiberDirection, a0
    //const double lambda = displacement[0] * fiberDirection[0] + displacement[1] * fiberDirection[1] + displacement[2] * fiberDirection[2];

    Vec3 fiberDirectionCurrentConfiguration = deformationGradient * fiberDirection;
    const double lambda = MathUtility::norm<3>(fiberDirectionCurrentConfiguration);

    // compute lambda lambda dot
    // d/dt dx = d/dt F
    // d/dt lambda = d/dt ||dx•a0|| = 1 / ||Fa|| (Fa0 • dot{F}a0) = 1/lambda (Fa • Fdot a0)
    // d/dt lambda = d/dt
    const double lambdaDot = 0.0;

    //

    lambdaVariable->setValue(dofNoLocal, lambda);
    lambdaDotVariable->setValue(dofNoLocal, lambdaDot);
  }

  lambdaVariable->zeroGhostBuffer();
  lambdaVariable->finishGhostManipulation();
  lambdaVariable->startGhostManipulation();

  lambdaDotVariable->zeroGhostBuffer();
  lambdaDotVariable->finishGhostManipulation();
  lambdaDotVariable->startGhostManipulation();
}

template<typename MeshType>
void MuscleContractionSolver<MeshType>::
computeActiveStress()
{
  LOG(DEBUG) << "computeActiveStress";

  typedef typename DynamicHyperelasticitySolverType::HyperelasticitySolverType::StressFieldVariableType StressFieldVariableType;
  typedef typename DynamicHyperelasticitySolverType::HyperelasticitySolverType::DisplacementsFieldVariableType DisplacementsFieldVariableType;
  typedef typename Data::ScalarFieldVariableType FieldVariableType;

  std::shared_ptr<StressFieldVariableType> activePK2StressVariable;
  std::shared_ptr<DisplacementsFieldVariableType> fiberDirectionVariable;

  if (isDynamic_)
  {
    activePK2StressVariable = dynamicHyperelasticitySolver_->hyperelasticitySolver().data().activePK2Stress();
    fiberDirectionVariable = dynamicHyperelasticitySolver_->hyperelasticitySolver().data().fiberDirection();
  }
  else
  {
    activePK2StressVariable = staticHyperelasticitySolver_->data().activePK2Stress();
    fiberDirectionVariable = staticHyperelasticitySolver_->data().fiberDirection();
  }

  std::shared_ptr<FieldVariableType> lambdaVariable = data_.lambda();
  std::shared_ptr<FieldVariableType> gammaVariable = data_.gamma();

  // Heidlauf 2013: "Modeling the Chemoelectromechanical Behavior of Skeletal Muscle Using the Parallel Open-Source Software Library OpenCMISS", p.4, Eq. (11)

  const double lambdaOpt = 1.2;

  // loop over local degrees of freedom
  for (dof_no_t dofNoLocal = 0; dofNoLocal < data_.functionSpace()->nDofsLocalWithoutGhosts(); dofNoLocal++)
  {
    const Vec3 fiberDirection = fiberDirectionVariable->getValue(dofNoLocal);

    const double lambda = lambdaVariable->getValue(dofNoLocal);
    const double gamma = gammaVariable->getValue(dofNoLocal);
    const double lambdaRelative = lambda / lambdaOpt;

    // compute f function
    double f = 0;
    if (0.6 <= lambdaRelative && lambdaRelative <= 1.4)
    {
      f = -25./4 * lambdaRelative*lambdaRelative + 25./2 * lambdaRelative - 5.25;
    }

    const double factor = 1./lambda * pmax_ * f * gamma;

    // Voigt notation:
    // [0][0] -> [0];
    // [1][1] -> [1];
    // [2][2] -> [2];
    // [0][1] -> [3];
    // [1][0] -> [3];
    // [1][2] -> [4];
    // [2][1] -> [4];
    // [0][2] -> [5];
    // [2][0] -> [5];

    VecD<6> activeStress;
    activeStress[0] = factor * fiberDirection[0] * fiberDirection[0];
    activeStress[1] = factor * fiberDirection[1] * fiberDirection[1];
    activeStress[2] = factor * fiberDirection[2] * fiberDirection[2];
    activeStress[3] = factor * fiberDirection[0] * fiberDirection[1];
    activeStress[4] = factor * fiberDirection[1] * fiberDirection[2];
    activeStress[5] = factor * fiberDirection[0] * fiberDirection[2];

    LOG(DEBUG) << "dof " << dofNoLocal << ", lambda: " << lambda << ", lambdaRelative: " << lambdaRelative
      << ", pmax_: " << pmax_ << ", f: " << f << ", gamma: " << gamma << ", => factor: " << factor << ", fiberDirection: " << fiberDirection;

    // if lambda is not yet computed (before first computation), set active stress to zero
    if (fabs(lambda)  < 1e-12)
    {
      activeStress = VecD<6>{0.0};
    }

    LOG(DEBUG) << "set active stress to " << activeStress;
    activePK2StressVariable->setValue(dofNoLocal, activeStress, INSERT_VALUES);
  }

  activePK2StressVariable->zeroGhostBuffer();
  activePK2StressVariable->finishGhostManipulation();
  activePK2StressVariable->startGhostManipulation();
}

template<typename MeshType>
void MuscleContractionSolver<MeshType>::
mapGeometryToGivenMeshes()
{
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_+std::string("_map_geometry"));

  LOG(DEBUG) << "mapGeometryToGivenMeshes: meshNamesOfGeometryToMapTo: " << meshNamesOfGeometryToMapTo_;
  if (!meshNamesOfGeometryToMapTo_.empty())
  {
    using SourceFunctionSpaceType = typename StaticHyperelasticitySolverType::DisplacementsFunctionSpace;
    using SourceFieldVariableType = FieldVariable::FieldVariable<SourceFunctionSpaceType,3>;
    
    assert(data_.functionSpace());

    // get source field variable
    std::shared_ptr<SourceFieldVariableType> geometryFieldSource = std::make_shared<SourceFieldVariableType>(data_.functionSpace()->geometryField());

    std::vector<Vec3> geometryValuesSource;
    geometryFieldSource->getValuesWithoutGhosts(geometryValuesSource);

    LOG(DEBUG) << "geometryValuesSource: " << geometryValuesSource;

    // loop over all given mesh names to which we should transfer the geometry
    for (std::string meshName : meshNamesOfGeometryToMapTo_)
    {
      // for first order meshes
      using TargetFunctionSpaceType1 = ::FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunction::LagrangeOfOrder<1>>;
      using TargetFieldVariableType1 = FieldVariable::FieldVariable<TargetFunctionSpaceType1,3>;
      
      // if the mesh name corresponds to a linear mesh
      if (this->context_.meshManager()->hasFunctionSpaceOfType<TargetFunctionSpaceType1>(meshName))
      {
        // get target geometry field variable
        std::shared_ptr<TargetFieldVariableType1> geometryFieldTarget = std::make_shared<TargetFieldVariableType1>(
          this->context_.meshManager()->functionSpace<TargetFunctionSpaceType1>(meshName)->geometryField());

        LOG(DEBUG) << "transfer geometry field to linear mesh, " << geometryFieldSource->functionSpace()->meshName() << " -> "
          << geometryFieldTarget->functionSpace()->meshName();
        LOG(DEBUG) << StringUtility::demangle(typeid(SourceFunctionSpaceType).name()) << " -> " << StringUtility::demangle(typeid(TargetFunctionSpaceType1).name());

        // perform the mapping
        DihuContext::mappingBetweenMeshesManager()->template prepareMapping<SourceFieldVariableType,TargetFieldVariableType1>(geometryFieldSource, geometryFieldTarget, -1);

        // map the whole geometry field (all components), do not avoid copy
        DihuContext::mappingBetweenMeshesManager()->template map<SourceFieldVariableType,TargetFieldVariableType1>(geometryFieldSource, geometryFieldTarget, -1, -1, false);
        DihuContext::mappingBetweenMeshesManager()->template finalizeMapping<SourceFieldVariableType,TargetFieldVariableType1>(geometryFieldSource, geometryFieldTarget, -1, -1, false);
      }

      // for second order meshes
      using TargetFunctionSpaceType2 = ::FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunction::LagrangeOfOrder<2>>;
      using TargetFieldVariableType2 = FieldVariable::FieldVariable<TargetFunctionSpaceType2,3>;
      
      // if the mesh name corresponds to a quadratic mesh
      if (this->context_.meshManager()->hasFunctionSpaceOfType<TargetFunctionSpaceType2>(meshName))
      {
        // get target geometry field variable
        std::shared_ptr<TargetFieldVariableType2> geometryFieldTarget = std::make_shared<TargetFieldVariableType2>(
          this->context_.meshManager()->functionSpace<TargetFunctionSpaceType2>(meshName)->geometryField());

        LOG(DEBUG) << "transfer geometry field to quadratic mesh, " << geometryFieldSource->functionSpace()->meshName() << " -> "
          << geometryFieldTarget->functionSpace()->meshName();
        LOG(DEBUG) << StringUtility::demangle(typeid(SourceFunctionSpaceType).name()) << " -> " << StringUtility::demangle(typeid(TargetFunctionSpaceType2).name());

        // perform the mapping
        DihuContext::mappingBetweenMeshesManager()->template prepareMapping<SourceFieldVariableType,TargetFieldVariableType2>(geometryFieldSource, geometryFieldTarget, -1);

        // map the whole geometry field (all components), do not avoid copy
        DihuContext::mappingBetweenMeshesManager()->template map<SourceFieldVariableType,TargetFieldVariableType2>(geometryFieldSource, geometryFieldTarget, -1, -1, false);
        DihuContext::mappingBetweenMeshesManager()->template finalizeMapping<SourceFieldVariableType,TargetFieldVariableType2>(geometryFieldSource, geometryFieldTarget, -1, -1, false);
      }
    }
  }

  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_+std::string("_map_geometry"));
}


//! get a reference to the DynamicHyperelasticitySolverType
template<typename MeshType>
std::shared_ptr<typename MuscleContractionSolver<MeshType>::DynamicHyperelasticitySolverType> MuscleContractionSolver<MeshType>::
dynamicHyperelasticitySolver()
{
  return dynamicHyperelasticitySolver_;
}

template<typename MeshType>
typename MuscleContractionSolver<MeshType>::Data &MuscleContractionSolver<MeshType>::
data()
{
  return data_;
}

//! get the data that will be transferred in the operator splitting to the other term of the splitting
//! the transfer is done by the slot_connector_data_transfer class
template<typename MeshType>
std::shared_ptr<typename MuscleContractionSolver<MeshType>::SlotConnectorDataType> MuscleContractionSolver<MeshType>::
getSlotConnectorData()
{
  return data_.getSlotConnectorData();
}
