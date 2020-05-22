#include "function_space/function_space.h"

#include <cmath>
#include <array>
#include <sstream>

#include "easylogging++.h"

namespace FunctionSpace
{

template<typename MeshType, typename BasisFunctionType>
std::array<dof_no_t,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()> FunctionSpace<MeshType,BasisFunctionType>::
getElementDofNosLocal(element_no_t elementNo) const
{
  std::array<dof_no_t,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()> dofNosLocal;
  for (int dofIndex = 0; dofIndex < FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement(); dofIndex++)
  {
    dofNosLocal[dofIndex] = this->getDofNo(elementNo, dofIndex);
  }
  return dofNosLocal;
}

//! vectorized version of getElementDofNosLocal
template<typename MeshType, typename BasisFunctionType>
std::array<Vc::int_v,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()> FunctionSpace<MeshType,BasisFunctionType>::
getElementDofNosLocal(Vc::int_v elementNo) const
{
  std::array<Vc::int_v,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()> dofNosLocal;
  for (int dofIndex = 0; dofIndex < FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement(); dofIndex++)
  {
    for (int vcComponent = 0; vcComponent < Vc::double_v::size(); vcComponent++)
    {
      if (elementNo[vcComponent] == -1)
        dofNosLocal[dofIndex][vcComponent] = -1;
      else
        dofNosLocal[dofIndex][vcComponent] = this->getDofNo(elementNo[vcComponent], dofIndex);
    }
  }
  return dofNosLocal;
}

template<typename MeshType, typename BasisFunctionType>
void FunctionSpace<MeshType,BasisFunctionType>::
getElementDofNosLocal(element_no_t elementNo, std::vector<dof_no_t> &dofNosLocal) const
{
  dofNosLocal.resize(FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement());
  for (int dofIndex = 0; dofIndex < FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement(); dofIndex++)
  {
    dofNosLocal[dofIndex] = this->getDofNo(elementNo, dofIndex);
  }
}

template<typename MeshType, typename BasisFunctionType>
void FunctionSpace<MeshType,BasisFunctionType>::
getElementDofNosLocalWithoutGhosts(element_no_t elementNo, std::vector<dof_no_t> &dofNosLocal) const
{
  dofNosLocal.reserve(FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement());
  for (int dofIndex = 0; dofIndex < FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement(); dofIndex++)
  {
    dof_no_t dofNoLocal = this->getDofNo(elementNo, dofIndex);
    if (dofNoLocal < this->meshPartition_->nDofsLocalWithoutGhosts())
    {
      dofNosLocal.push_back(dofNoLocal);
    }
  }
}


} // namespace
