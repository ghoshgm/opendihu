#include "spatial_discretization/finite_element_method/02_boundary_conditions.h"

#include <Python.h>
#include <memory>
#include <vector>
#include <petscsys.h>

#include "quadrature/tensor_product.h"
#include "utility/vector_operators.h"

namespace SpatialDiscretization
{

template<typename FunctionSpaceType,typename QuadratureType,typename Term,typename Dummy>
void BoundaryConditions<FunctionSpaceType,QuadratureType,Term,Dummy>::
setBoundaryConditionHandlingEnabled(bool boundaryConditionHandlingEnabled)
{
  boundaryConditionHandlingEnabled_ = boundaryConditionHandlingEnabled;
}

template<typename FunctionSpaceType,typename QuadratureType,typename Term,typename Dummy>
void BoundaryConditions<FunctionSpaceType,QuadratureType,Term,Dummy>::
setDirichletBoundaryConditions(std::shared_ptr<DirichletBoundaryConditions<FunctionSpaceType,1>> dirichletBoundaryConditions)
{
  this->dirichletBoundaryConditions_ = dirichletBoundaryConditions;
}

template<typename FunctionSpaceType,typename QuadratureType,typename Term,typename Dummy>
void BoundaryConditions<FunctionSpaceType,QuadratureType,Term,Dummy>::
applyBoundaryConditions()
{
  if (!boundaryConditionHandlingEnabled_)
  {
    if (this->specificSettings_.hasKey("dirichletBoundaryConditions"))
    {
      LOG(WARNING) << "You have specified dirichlet boundary conditions in FiniteElementMethod via the key \"dirichletBoundaryConditions\". "
        << "They are not used here, e.g. because the FiniteElementMethod is wrapped by a time stepping scheme. Consider only setting dirichlet boundary conditions in the time stepping scheme.";
    }

    VLOG(1) << "do not handle boundary conditions in finite element method, because boundaryConditionHandlingEnabled=false";
    return;
  }

  LOG(TRACE) << "FiniteElementMethod::applyBoundaryConditions";

  if (VLOG_IS_ON(4))
  {
    VLOG(4) << "Finite element data before applyBoundaryConditions";
    this->data_.print();
  }

  if (dirichletBoundaryConditions_ == nullptr)
  {
    dirichletBoundaryConditions_ = std::make_shared<DirichletBoundaryConditions<FunctionSpaceType,1>>();
    dirichletBoundaryConditions_->initialize(this->specificSettings_, this->data_.functionSpace());
  }

  // get abbreviations
  std::shared_ptr<FunctionSpaceType> functionSpace = this->data_.functionSpace();
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> rightHandSide = this->data_.rightHandSide();
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> stiffnessMatrix = this->data_.stiffnessMatrix();

  // apply the boundary conditions in stiffness matrix (set bc rows and columns of matrix to 0 and diagonal to 1), also add terms with matrix entries to rhs
  dirichletBoundaryConditions_->applyInSystemMatrix(stiffnessMatrix, rightHandSide);

  // set prescribed values in rhs
  dirichletBoundaryConditions_->applyInRightHandSide(rightHandSide, rightHandSide);

  if (VLOG_IS_ON(4))
  {
    VLOG(4) << "Finite element data after applyBoundaryConditions";
    this->data_.print();
  }
}

} // namespace
