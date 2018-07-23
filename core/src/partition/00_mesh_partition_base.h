#pragma once

#include <memory>
#include <petscdmda.h>

#include "control/types.h"
#include "partition/rank_subset.h"

namespace Partition
{

/** base class for mesh partition */
class MeshPartitionBase
{
public:
 
  //! constructor, store the rankSubset
  MeshPartitionBase(std::shared_ptr<RankSubset> rankSubset);
  
  //! virtual destructor
  virtual ~MeshPartitionBase();
 
  //! number of ranks
  int nRanks();
  
  //! number of entries in the current partition (this usually refers to the elements)
  virtual element_no_t localSize() = 0;
  
  //! number of nodes in total
  virtual global_no_t globalSize() = 0;
  
  //! return reference to a vector containing all local dofs, i.e. a vector with {0,1,2,...,localSize()-1}
  std::vector<PetscInt> &localDofs();
  
  //! get the MPI communicator that is needed for the work portion
  MPI_Comm mpiCommunicator();
  
  //! get an AO object
  virtual AO &applicationOrdering() = 0;
  
  //! from a vector of global numbers remove all that are non-local
  template <typename T>
  virtual void extractLocalNumbers(std::vector<T> &vector) = 0;
  
protected:
 
  //! initialize the localDofs_ vector to values {0,1,...,localSize()-1}, this needs localSize set previously
  void initializeLocalDofs();
 
  std::vector<PetscInt> localDofs_;     ///< list of local dofs, which is {0,1,...,localSize()-1}. This is needed for calls to Petsc functions that access all local data, e.g. within fieldVariable->getValues
  std::shared_ptr<RankSubset> rankSubset_;  ///< the set of ranks that compute something where this partition is a part of, also holds the MPI communciator
};

}  // namespace
#include "partition/00_mesh_partition_base.tpp"
