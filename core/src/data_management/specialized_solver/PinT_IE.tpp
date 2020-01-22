#include "data_management/specialized_solver/PinT_IE.h"

#include "easylogging++.h"

#include "utility/python_utility.h"
#include "control/dihu_context.h"
#include "utility/petsc_utility.h"

namespace Data
{

template<typename FunctionSpaceType>
PinT<FunctionSpaceType>::
PinT(DihuContext context) :
  Data<FunctionSpaceType>::Data(context)
{
}

template<typename FunctionSpaceType>
void PinT<FunctionSpaceType>::
initialize()
{
  // call initialize of base class
  Data<FunctionSpaceType>::initialize();

  // create th output connector data object
  outputConnectorData_ = std::make_shared<OutputConnectorDataType>();

  // add all needed field variables to be transferred

  // add component 0 of fieldvariableA_
  outputConnectorData_->addFieldVariable(this->solution_, 0);

  // There is addFieldVariable(...) and addFieldVariable2(...) for the two different field variable types,
  // Refer to "data_management/output_connector_data.h" for details.

  // you can also access settings from the python config here:
  std::string option1 = this->context_.getPythonConfig().getOptionString("option1", "default string");
  LOG(DEBUG) << "In data object, parsed option1: [" << option1 << "].";
}

template<typename FunctionSpaceType>
void PinT<FunctionSpaceType>::
createPetscObjects()
{
  assert(this->functionSpace_);

  // Here, the actual field variables will be created.
  // Make sure, the number of components matches. The string is the name of the field variable. It will also be used in the VTK output files.
  this->fieldVariableB_ = this->functionSpace_->template createFieldVariable<1>("b");
}

// ... add a "getter" method for each fieldvariable with the same name as the field variable (but without underscore)
template<typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> PinT<FunctionSpaceType>::
solution()
{
  return this->solution_;
}

template<typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> PinT<FunctionSpaceType>::
fieldVariableB()
{
  return this->fieldVariableB_;
}

template<typename FunctionSpaceType>
void PinT<FunctionSpaceType>::
setSolutionVariable(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> solution)
{
  // Store the solution field variable that was given as argument.
  // Because this is a pointer, is has access to the original solution variable and no data copy is performed.
  solution_ = solution;
}


template<typename FunctionSpaceType>
std::shared_ptr<typename PinT<FunctionSpaceType>::OutputConnectorDataType> PinT<FunctionSpaceType>::
getOutputConnectorData()
{
  // return the output connector data object
  return this->outputConnectorData_;
}

template<typename FunctionSpaceType>
typename PinT<FunctionSpaceType>::FieldVariablesForOutputWriter PinT<FunctionSpaceType>::
getFieldVariablesForOutputWriter()
{
  // these field variables will be written to output files by the output writer

  // get the geometry field, which is always needed, from the function space
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> geometryField
    = std::make_shared<FieldVariable::FieldVariable<FunctionSpaceType,3>>(this->functionSpace_->geometryField());

  return std::make_tuple(
    geometryField,
    this->solution_,
    this->fieldVariableB_   // add all field variables that should appear in the output file. Of course, this list has to match the type in the header file.
  );
}

} // namespace
