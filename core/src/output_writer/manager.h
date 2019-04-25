#pragma once

#include <Python.h>  // has to be the first included header
#include <list>
#include <memory>

#include "control/types.h"
#include "data_management/data.h"
#include "output_writer/generic.h"

namespace OutputWriter
{

class Manager
{
public:
  //! call all output writers to write output, timeStepNo of -1 means no time step number in output filename
  template<typename DataType>
  void writeOutput(DataType &problemData, int timeStepNo = -1, double currentTime = 0.0) const;

  //! parse settings and create output writers from specification in "OutputWriter" list
  void initialize(DihuContext context, PythonConfig settings);

  //! if this manager contains any output writers
  bool hasOutputWriters();

protected:

  //! helper function that creates an outputWriter
  void createOutputWriterFromSettings(DihuContext context, PythonConfig settings);

  std::list<std::shared_ptr<Generic>> outputWriter_;    ///< list of active output writers
};

}  // namespace OutputWriter

#include "output_writer/manager.tpp"
