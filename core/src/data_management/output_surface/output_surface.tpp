#include "data_management/output_surface/output_surface.h"

namespace Data
{

template<typename Data3D>
OutputSurface<Data3D>::
OutputSurface(DihuContext context) :
  Data<typename Data3D::FunctionSpace>(context)
{
  std::string faceStr = this->context_.getPythonConfig().getOptionString("face","2-");
  face_ = Mesh::parseFace(faceStr);

  if (face_ != Mesh::face_t::face2Minus && face_ != Mesh::face_t::face2Plus)
  {
    LOG(ERROR) << this->context_.getPythonConfig() << ": only \"2-\" and \"2+\" are supported! Now using \"2-\".";
    face_ = Mesh::face_t::face2Minus;
  }

  LOG(DEBUG) << "OutputSurface: parsed face \"" << Mesh::getString(face_) << "\".";
}

template<typename Data3D>
void OutputSurface<Data3D>::
initialize()
{
}

template<typename Data3D>
void OutputSurface<Data3D>::
setData(Data3D &data3d)
{
  this->data3D_ = std::make_shared<Data3D>(data3d);
}

template<typename Data3D>
void OutputSurface<Data3D>::
createPetscObjects()
{
}

template<typename Data3D>
void OutputSurface<Data3D>::
print()
{

}

template<typename Data3D>
typename OutputSurface<Data3D>::OutputFieldVariables OutputSurface<Data3D>::
getOutputFieldVariables()
{
  typename Data3D::OutputFieldVariables outputFieldVariables3D = this->data3d_->getOutputFieldVariables();
  OutputFieldVariables outputFieldVariables2D;

  // create surface field variables for 3D field variables (the template magic happens in convert_output_field_variables.h)
  ConvertOutputFieldVariables<typename Data3D::OutputFieldVariables>::convert(outputFieldVariables3D, outputFieldVariables2D, face_);
  return outputFieldVariables2D;
}

} // namespace
