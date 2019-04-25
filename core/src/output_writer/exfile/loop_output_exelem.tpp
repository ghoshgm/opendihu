#include "output_writer/exfile/loop_output_exelem.h"

#include "output_writer/loop_count_n_field_variables_of_mesh.h"

#include <cstdlib>

namespace OutputWriter
{

namespace ExfileLoopOverTuple
{
 
 /** Static recursive loop from 0 to number of entries in the tuple
 * Loop body
 */
template<typename OutputFieldVariablesType, typename AllOutputFieldVariablesType, int i>
inline typename std::enable_if<i < std::tuple_size<OutputFieldVariablesType>::value, void>::type
loopOutputExelem(const OutputFieldVariablesType &fieldVariables, const AllOutputFieldVariablesType &allFieldVariables,
                 std::string meshName, std::ofstream &file, std::shared_ptr<Mesh::Mesh> &mesh
)
{
  // call what to do in the loop body
  if (outputExelem<typename std::tuple_element<i,OutputFieldVariablesType>::type, AllOutputFieldVariablesType>(
        std::get<i>(fieldVariables), allFieldVariables, meshName, file, mesh))
    return;
  
  // advance iteration to next tuple element
  loopOutputExelem<OutputFieldVariablesType, AllOutputFieldVariablesType, i+1>(fieldVariables, allFieldVariables, meshName, file, mesh);
}
 
// current element is of pointer type (not vector)
template<typename CurrentFieldVariableType, typename OutputFieldVariablesType>
typename std::enable_if<!TypeUtility::isTuple<CurrentFieldVariableType>::value && !TypeUtility::isVector<CurrentFieldVariableType>::value, bool>::type
outputExelem(CurrentFieldVariableType currentFieldVariable, const OutputFieldVariablesType &fieldVariables, std::string meshName, 
             std::ofstream &file, std::shared_ptr<Mesh::Mesh> &mesh)
{
  // if mesh name is the specified meshName
  if (currentFieldVariable->functionSpace()->meshName() == meshName)
  {
    // here we have the type of the mesh with meshName (which is typedef to FunctionSpace)
    typedef typename CurrentFieldVariableType::element_type::FunctionSpace FunctionSpace;
   
    int nFieldVariablesInMesh = 0;
    LoopOverTuple::loopCountNFieldVariablesOfMesh(fieldVariables, meshName, nFieldVariablesInMesh);
    // get mesh
    mesh = currentFieldVariable->functionSpace();
    
    // call exfile writer to output all field variables with the meshName
    ExfileWriter<FunctionSpace, OutputFieldVariablesType>::outputExelem(file, fieldVariables, meshName, currentFieldVariable->functionSpace(), nFieldVariablesInMesh);
   
    return true;  // break iteration
  }
  
  return false;  // do not break iteration 
}

// Elementent i is of vector type
template<typename VectorType, typename OutputFieldVariablesType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
outputExelem(VectorType currentFieldVariableVector, const OutputFieldVariablesType &fieldVariables, std::string meshName, 
             std::ofstream &file, std::shared_ptr<Mesh::Mesh> &mesh)
{
  for (auto& currentFieldVariable : currentFieldVariableVector)
  {
    // call function on all vector entries
    if (outputExelem<typename VectorType::value_type,OutputFieldVariablesType>(currentFieldVariable, fieldVariables, meshName, file, mesh))
      return true; // break iteration
  }
  return false;  // do not break iteration 
}

// element i is of tuple type
template<typename TupleType, typename AllOutputFieldVariablesType>
typename std::enable_if<TypeUtility::isTuple<TupleType>::value, bool>::type
outputExelem(TupleType currentFieldVariableTuple, const AllOutputFieldVariablesType &fieldVariables, std::string meshName, 
             std::ofstream &file, std::shared_ptr<Mesh::Mesh> &mesh)
{
  // call for tuple element
  loopOutputExelem<TupleType, AllOutputFieldVariablesType>(currentFieldVariableTuple, fieldVariables, meshName, file, mesh);
  
  return false;  // do not break iteration 
}

}  // namespace ExfileLoopOverTuple
}  // namespace OutputWriter
