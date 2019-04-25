#include "output_writer/loop_count_n_field_variables_of_mesh.h"

#include <cstdlib>

namespace OutputWriter
{

namespace LoopOverTuple
{
 
/** Static recursive loop from 0 to number of entries in the tuple
 * Loop body
 */
template<typename OutputFieldVariablesType, int i>
inline typename std::enable_if<i < std::tuple_size<OutputFieldVariablesType>::value, void>::type
loopCountNFieldVariablesOfMesh(const OutputFieldVariablesType &fieldVariables, std::string meshName,
                               int &nFieldVariablesOfMesh)
{
  // call what to do in the loop body
  if (countNFieldVariablesOfMesh<typename std::tuple_element<i,OutputFieldVariablesType>::type>(std::get<i>(fieldVariables), meshName, nFieldVariablesOfMesh))
    return;
  
  // advance iteration to next tuple element
  loopCountNFieldVariablesOfMesh<OutputFieldVariablesType, i+1>(fieldVariables, meshName, nFieldVariablesOfMesh);
}
 
// current element is of pointer type (not vector)
template<typename CurrentFieldVariableType>
typename std::enable_if<!TypeUtility::isTuple<CurrentFieldVariableType>::value && !TypeUtility::isVector<CurrentFieldVariableType>::value, bool>::type
countNFieldVariablesOfMesh(CurrentFieldVariableType currentFieldVariable, std::string meshName,
                           int &nFieldVariablesOfMesh)
{
  // if mesh name is the specified meshName, count this mesh
  if (currentFieldVariable->functionSpace()->meshName() == meshName)
  {
    nFieldVariablesOfMesh++;
  }
  
  return false;  // do not break iteration 
}

// element i is of vector type
template<typename VectorType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
countNFieldVariablesOfMesh(VectorType currentFieldVariableVector, std::string meshName,
                           int &nFieldVariablesOfMesh)
{
  for (auto& currentFieldVariable : currentFieldVariableVector)
  {
    // call function on all vector entries
    if (countNFieldVariablesOfMesh<typename VectorType::value_type>(currentFieldVariable, meshName, nFieldVariablesOfMesh))
      return true;
  }
  
  return false;  // do not break iteration
}

// element i is of tuple type
template<typename TupleType>
typename std::enable_if<TypeUtility::isTuple<TupleType>::value, bool>::type
countNFieldVariablesOfMesh(TupleType currentFieldVariableTuple, std::string meshName,
                           int &nFieldVariablesOfMesh)
{
  // call for tuple element
  loopCountNFieldVariablesOfMesh<TupleType>(currentFieldVariableTuple, meshName, nFieldVariablesOfMesh);
 
  return false;  // do not break iteration
}

}  // namespace LoopOverTuple
}  // namespace OutputWriter
