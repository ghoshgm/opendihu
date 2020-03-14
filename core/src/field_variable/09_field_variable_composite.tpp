#include "field_variable/09_field_variable_composite.h"

namespace FieldVariable
{

//! split this field variables into field variables for the sub function space and set the corresponding values
template<int D,typename BasisFunctionType,int nComponents>
void FieldVariableComposite<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,nComponents>::
getSubFieldVariables(std::vector<std::shared_ptr<FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>, BasisFunctionType>, nComponents>>> &subFieldVariables)
{
  subFieldVariables.clear();

  updateSubFieldVariables();
  subFieldVariables = subFieldVariables_;
}

//! get the sub field variable no i
template<int D,typename BasisFunctionType,int nComponents>
std::shared_ptr<FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>, BasisFunctionType>, nComponents>> FieldVariableComposite<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,nComponents>::
subFieldVariable(int i)
{
  updateSubFieldVariables();

  if (i == -1)
    i = subFieldVariables_.size() - 1;
  assert(i < subFieldVariables_.size());
  return subFieldVariables_[i];
}

//! get the sub field variable no i
template<int D,typename BasisFunctionType,int nComponents>
std::shared_ptr<FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>, BasisFunctionType>, nComponents>> FieldVariableComposite<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,nComponents>::
subFieldVariableWithoutUpdate(int i)
{
  // Do not update the sub field variables because this would need a collective call by all ranks of this field variable. This is not the case for face computations.

  if (subFieldVariables_.empty())
  {
    LOG(FATAL) << "subFieldVariableWithoutUpdate(" << i << ") called when sub field variables have not yet been initialized.";
  }

  if (i == -1)
    i = subFieldVariables_.size() - 1;

  if (i >= subFieldVariables_.size())
  {
    LOG(FATAL) << "Index of subFieldVariable is wrong. i: " << i << ", size: " << subFieldVariables_.size();
  }
  assert(i < subFieldVariables_.size());
  return subFieldVariables_[i];
}

template<int D,typename BasisFunctionType,int nComponents>
void FieldVariableComposite<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,nComponents>::
updateSubFieldVariables()
{
  // update sub field variables

  // get subFunctionSpaces
  std::vector<std::shared_ptr<SubFunctionSpaceType>> subFunctionSpaces = this->functionSpace_->subFunctionSpaces();

  subFieldVariables_.resize(subFunctionSpaces.size());

  // get own values
  std::vector<VecD<nComponents>> ownValues;
  this->getValuesWithoutGhosts(ownValues);
  assert(ownValues.size() == this->functionSpace_->meshPartition()->nDofsLocalWithoutGhosts());

  VLOG(1) << "getSubFieldVariables for \"" << this->name_ << "\", " << subFunctionSpaces.size() << " sub meshes";
  VLOG(2) << ownValues.size() << " ownValues: " << ownValues;

  // loop over sub function spaces
  int subMeshNo = 0;
  for (std::shared_ptr<SubFunctionSpaceType> &subFunctionSpace : subFunctionSpaces)
  {
    VLOG(1) << "subMeshNo " << subMeshNo;

    // create sub field variable if it does not yet exist
    if (!subFieldVariables_[subMeshNo])
    {
      // determine name of sub field variable
      std::stringstream fieldVariableName;
      fieldVariableName << this->name_;  // << "_" << subFunctionSpace->meshName();

      // create sub field variable
      std::vector<std::string> componentNames(this->componentNames_.begin(), this->componentNames_.end());

      subFieldVariables_[subMeshNo] = subFunctionSpace->template createFieldVariable<nComponents>(fieldVariableName.str(), componentNames);
      subFieldVariables_[subMeshNo]->setIsGeometryField(this->isGeometryField());
    }

    // determine values for the sub field variable
    subFieldVariableValues_.resize(subFunctionSpace->nDofsLocalWithoutGhosts());

    // loop over dofs
    for (node_no_t nodeNoLocal = 0; nodeNoLocal < subFunctionSpace->nNodesLocalWithoutGhosts(); nodeNoLocal++)
    {
      bool nodeIsShared = false;
      node_no_t compositeNodeNo = this->functionSpace_->meshPartition()->getNodeNoLocalFromSubmesh(subMeshNo, nodeNoLocal, nodeIsShared);

      const int nDofsPerNode = subFunctionSpace->nDofsPerNode();
      for (int nodalDofNo = 0; nodalDofNo < nDofsPerNode; nodalDofNo++)
      {
        dof_no_t compositeDofNo = compositeNodeNo*nDofsPerNode + nodalDofNo;
        dof_no_t subFunctionSpaceDofNo = nodeNoLocal*nDofsPerNode + nodalDofNo;

        assert (subFunctionSpaceDofNo < subFunctionSpace->nDofsLocalWithoutGhosts());
        assert (compositeDofNo < this->functionSpace_->meshPartition()->nDofsLocalWithoutGhosts());

        VLOG(2) << "subFieldVariableValues[" << subFunctionSpaceDofNo << "] = values[ " << compositeDofNo << "], nodeIsShared: " << nodeIsShared;
        subFieldVariableValues_[subFunctionSpaceDofNo] = ownValues[compositeDofNo];
      }
    }

    // set values in the sub field variable
    //subFieldVariables_[subMeshNo]->startGhostManipulation();
    subFieldVariables_[subMeshNo]->setValuesWithoutGhosts(subFieldVariableValues_);
    subFieldVariables_[subMeshNo]->finishGhostManipulation();
    subFieldVariables_[subMeshNo]->startGhostManipulation();
    subFieldVariables_[subMeshNo]->zeroGhostBuffer();
    subFieldVariables_[subMeshNo]->finishGhostManipulation();

    VLOG(1) << "new field variable : " << *subFieldVariables_[subMeshNo];

    subMeshNo++;
  }
}

} // namespace
