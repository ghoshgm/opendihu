#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <vector>
#include <atomic>

#include "control/types.h"
#include "output_writer/generic.h"

namespace OutputWriter
{

struct BoundingBox
{
  Vec3 min, max;   // minimum,maximum value of bounding box
  BoundingBox();
};

struct MegaMolWriterContext
{
  std::vector<Vec3> geometryFieldValues;
  std::vector<double> scalarFieldVariableValues;
  double approximateDistanceBetweenFibers;
};

class MegaMol : public Generic
{
public:

  //! constructor
  MegaMol(DihuContext context, PythonConfig specificSettings);

  //! write out solution to given filename, if timeStepNo is not -1, this value will be part of the filename
  template<typename DataType>
  void write(DataType &data, int timeStepNo = -1, double currentTime = -1);

#ifdef HAVE_ADIOS
  struct adios_writer_t
  {
    MegaMolWriterContext megaMolWriterContext;   ///< all variables that are combined and collected during the calls until combineNInstances_ calls have been made

    std::shared_ptr<adios2::Engine> engine;    //< adios writer
    std::atomic_int nOpenWriters;            ///< this writer could be in multiple instances, called at the same time. Then just the first open should really open the file and the last close should close it. This variables counts the number of opens.
  };

private:

#if defined(HAVE_MEGAMOL) && defined(HAVE_ADIOS)
  void notifyMegaMol();
#endif

  static std::array<std::shared_ptr<adios_writer_t>, 2> adiosWriters_;   ///< the writers for front and back buffer
  int currentOpenWriterIndex_ = 0;   ///< which writer the last opened is, 0 or 1

  std::shared_ptr<adios2::Variable<double>> adiosFieldVariableGeometry_;    ///< the adios field variable for the geometry, with name "xyz"
  std::shared_ptr<adios2::Variable<double>> adiosFieldVariableScalar_;    ///< the adios scalar field variable, with name "i"
  std::shared_ptr<adios2::Variable<double>> boxVariable_;                   ///< the adios variable containing the bounding box information, with name "box"
  std::shared_ptr<adios2::Variable<double>> globalRadiusVariable_;         ///< the adios variable for the radius for visualization of spheres, with name "global_radius"
  std::shared_ptr<adios2::Variable<int>> globalNumberOfNodesVariable_;     ///< the adios variable for the global number of nodes, with name "p_count" / particle count

  double globalRadius_;      ///< the global radius to be used for the visualization
  int nNodesGlobal_ = 0;         ///< the global number of nodes or the particle count
  std::vector<double> boundingBoxValues_;    ///< the global bounding box that is needed for the visualization

  std::string currentFilename_;   ///< the file to which is currently being written, the file is not yet ready
  std::string lastFilename_;     ///< the last used filename for output

  int combineNInstances_;       ///< number of calls to MegaMol::write that should be combined into one timestep, this is needed when the output writer is given in a MultipleInstances environment
  bool useFrontBackBuffer_;      ///< if the two buffers, *_0 and *_1 should be used, otherwise is uses a standard file naming scheme with increasing counter number
#endif

};

} // namespace

#include "output_writer/megamol/megamol.tpp"
