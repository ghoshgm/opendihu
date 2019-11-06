#include "output_writer/megamol/megamol.h"

#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>

#include "easylogging++.h"
#include <utility/python_utility.h>

#ifdef HAVE_MEGAMOL
#include <libzmq/zmq.hpp>
#endif
namespace OutputWriter
{

#ifdef HAVE_ADIOS
std::map<std::string,std::array<std::shared_ptr<MegaMol::adios_writer_t>, 2>> MegaMol::adiosWriters_;
std::map<std::string,std::shared_ptr<MegaMol::adios_writer_t>> MegaMol::adiosWriter_;
#endif

BoundingBox::BoundingBox():
  min(Vec3({0.0,0.0,0.0})),
  max(Vec3({0.0,0.0,0.0}))
{

}

#ifdef HAVE_ADIOS
MegaMol::MegaMol(DihuContext context, PythonConfig settings, std::shared_ptr<Partition::RankSubset> rankSubset) :
  Generic(context, settings, rankSubset), currentOpenWriterIndex_(0)
{
  combineNInstances_ = specificSettings_.getOptionInt("combineNInstances", 1);
  useFrontBackBuffer_ = specificSettings_.getOptionBool("useFrontBackBuffer", true);
}
#else

MegaMol::MegaMol(DihuContext context, PythonConfig settings, std::shared_ptr<Partition::RankSubset> rankSubset) :
Generic(context, settings, rankSubset)
{
}
#endif

#ifdef HAVE_ADIOS

void MegaMol::writeAdiosVariables()
{
  std::shared_ptr<Partition::RankSubset> rankSubset = DihuContext::partitionManager()->rankSubsetForCollectiveOperations();

  int ownRankNo = rankSubset->ownRankNo();
  int nRanks = rankSubset->size();


  std::shared_ptr<adios_writer_t> &currentWriter = adiosWriter_[this->filenameBase_];

  std::vector<Vec3> &geometryFieldValues = currentWriter->megaMolWriterContext.geometryFieldValues;
  double approximateDistanceBetweenFibers = currentWriter->megaMolWriterContext.approximateDistanceBetweenFibers;

  // compute local bounding box
  BoundingBox localBoundingBox;
  localBoundingBox.min = geometryFieldValues[0];
  localBoundingBox.max = geometryFieldValues[0];

  // loop over geometry field entries
  for(std::vector<Vec3>::iterator iter = geometryFieldValues.begin(); iter != geometryFieldValues.end(); iter++)
  {
    for (int i = 0; i < 3; i++)
    {
      if ((*iter)[i] < localBoundingBox.min[i])
      {
        localBoundingBox.min[i] = (*iter)[i];
      }
      if ((*iter)[i] > localBoundingBox.max[i])
      {
        localBoundingBox.max[i] = (*iter)[i];
      }
    }
  }

  // reduce the bounding box values
  BoundingBox globalBoundingBox;
  MPI_Reduce(localBoundingBox.min.data(), globalBoundingBox.min.data(), 3, MPI_DOUBLE, MPI_MIN, 0, rankSubset->mpiCommunicator());
  MPI_Reduce(localBoundingBox.max.data(), globalBoundingBox.max.data(), 3, MPI_DOUBLE, MPI_MAX, 0, rankSubset->mpiCommunicator());

  LOG(DEBUG) << "reduced bounding box: " << globalBoundingBox.min << ", " << globalBoundingBox.max;

  std::array<double,6> globalBoundingBoxValues;
  globalBoundingBoxValues[0] = globalBoundingBox.min[0];
  globalBoundingBoxValues[1] = globalBoundingBox.min[1];
  globalBoundingBoxValues[2] = globalBoundingBox.max[2];
  globalBoundingBoxValues[3] = globalBoundingBox.max[0];
  globalBoundingBoxValues[4] = globalBoundingBox.max[1];
  globalBoundingBoxValues[5] = globalBoundingBox.min[2];

  std::array<double,6> localBoundingBoxValues;
  localBoundingBoxValues[0] = localBoundingBox.min[0];
  localBoundingBoxValues[1] = localBoundingBox.min[1];
  localBoundingBoxValues[2] = localBoundingBox.max[2];
  localBoundingBoxValues[3] = localBoundingBox.max[0];
  localBoundingBoxValues[4] = localBoundingBox.max[1];
  localBoundingBoxValues[5] = localBoundingBox.min[2];

  // reduce global number of nodes
  int nNodesLocal = geometryFieldValues.size();
  nNodesGlobal_ = 0;
  MPI_Reduce(&nNodesLocal, &nNodesGlobal_, 1, MPI_INT, MPI_SUM, 0, rankSubset->mpiCommunicator());

  // write geometry field

  // convert data to be send to ADIOS from vector<Vec3> to vector<double>

  std::vector<double> geometryFieldScalarValues(3*geometryFieldValues.size());
  for (int i = 0; i < geometryFieldValues.size(); i++)
  {
    for (int j = 0; j != 3; j++)
    {
      geometryFieldScalarValues[3*i + j] = geometryFieldValues[i][j];
    }
  }

  // define variable
  if (!adiosFieldVariableGeometry_)
  {
    std::string variableName = "xyz";

    // communicate offset into the global values array
    long long localSize = geometryFieldScalarValues.size();
    long long offset = 0;
    long long globalSize = 0;

    MPI_Exscan(&localSize, &offset, 1, MPI_LONG_LONG, MPI_SUM, rankSubset->mpiCommunicator());
    MPI_Allreduce(&localSize, &globalSize, 1, MPI_LONG_LONG, MPI_SUM, rankSubset->mpiCommunicator());

    LOG(DEBUG) << ownRankNo << "/" << nRanks << ": \"" << this->filenameBase_ << "\" "
      << "define variable \"" << variableName
      << "\", localSize: " << localSize << ", offset: " << offset << ", globalSize: " << globalSize;

    // name, global size, offset, local size
    adiosFieldVariableGeometry_ = std::make_shared<adios2::Variable<double>>(currentWriter->adiosIo->DefineVariable<double>(
      variableName, {(long unsigned int)globalSize}, {(long unsigned int)offset}, {(long unsigned int)localSize}, adios2::ConstantDims
    ));
  }

  // write data to file
  currentWriter->engine->Put<double>(*adiosFieldVariableGeometry_.get(), geometryFieldScalarValues.data());

  // write local bounding box
  if (!localBoundingBoxVariable_)
  {
    localBoundingBoxVariable_ = std::make_shared<adios2::Variable<double>>(
      currentWriter->adiosIo->DefineVariable<double>("localBoundingBox", {(long unsigned int)nRanks,6}, {(long unsigned int)ownRankNo,6}, {1,6}, adios2::ConstantDims));
  }
  currentWriter->engine->Put<double>(*localBoundingBoxVariable_.get(), localBoundingBoxValues.data());

  // write variables that have to be written by only one rank (rank 0)
  if (ownRankNo == 0)
  {
    // write global bounding box
    if (!globalBoundingBoxVariable_)
    {
      globalBoundingBoxVariable_ = std::make_shared<adios2::Variable<double>>(
        currentWriter->adiosIo->DefineVariable<double>("globalBoundingBox", {6},{0},{6}, adios2::ConstantDims));
    }
    currentWriter->engine->Put<double>(*globalBoundingBoxVariable_.get(), globalBoundingBoxValues.data());


    // write global nPointsPerCoordinateDirection
    if (!adiosNPointsPerCoordinateDirection_)
    {
      adiosNPointsPerCoordinateDirection_ = std::make_shared<adios2::Variable<int>>(
        currentWriter->adiosIo->DefineVariable<int>("nPointsPerCoordinateDirection", {3},{0},{3}, adios2::ConstantDims));
    }

    std::array<int,3> &nPointsPerCoordinateDirection = currentWriter->megaMolWriterContext.nPointsPerCoordinateDirection;
    currentWriter->engine->Put<int>(*adiosNPointsPerCoordinateDirection_.get(), nPointsPerCoordinateDirection.data());


    // write radius
    globalRadius_ = 0.1;
    if (approximateDistanceBetweenFibers > 1e-3)
    {
      globalRadius_ = approximateDistanceBetweenFibers*0.1;
    }
    if (!globalRadiusVariable_)
    {
      globalRadiusVariable_ = std::make_shared<adios2::Variable<double>>(currentWriter->adiosIo->DefineVariable<double>("global_radius"));
    }
    currentWriter->engine->Put<double>(*globalRadiusVariable_.get(), &globalRadius_);

    // write global number of nodes
    if (!globalNumberOfNodesVariable_)
    {
      globalNumberOfNodesVariable_ = std::make_shared<adios2::Variable<int>>(currentWriter->adiosIo->DefineVariable<int>("node_count"));
    }
    int nNodesGlobal = nNodesGlobal_;
    currentWriter->engine->Put<int>(*globalNumberOfNodesVariable_.get(), &nNodesGlobal);
  }


  // write scalar fields

  std::vector<double> &vmValues = currentWriter->megaMolWriterContext.vmValues;
  std::vector<double> &emgValues = currentWriter->megaMolWriterContext.emgValues;
  std::vector<double> &transmembraneFlowValues = currentWriter->megaMolWriterContext.transmembraneFlowValues;


  // communicate offset and global size
  int localSize = vmValues.size();
  int offset = 0;
  int globalSize = 0;

  MPI_Exscan(&localSize, &offset, 1, MPI_INT, MPI_SUM, rankSubset->mpiCommunicator());
  MPI_Allreduce(&localSize, &globalSize, 1, MPI_INT, MPI_SUM, rankSubset->mpiCommunicator());

  // define variable
  if (!adiosFieldVariableVm_)
  {
    // name, global size, offset, local size
    adiosFieldVariableVm_ = std::make_shared<adios2::Variable<double>>(currentWriter->adiosIo->DefineVariable<double>(
      "vm", {(long unsigned int)globalSize}, {(long unsigned int)offset}, {(long unsigned int)localSize}, adios2::ConstantDims
    ));

    LOG(DEBUG) << "define variable \"vm\", localSize: " << localSize << ", offset: " << offset << ", globalSize: " << globalSize;

    if (!emgValues.empty())
    {
      assert(vmValues.size() == emgValues.size());
      assert(vmValues.size() == transmembraneFlowValues.size());

      // name, global size, offset, local size
      adiosFieldVariableEmg_ = std::make_shared<adios2::Variable<double>>(currentWriter->adiosIo->DefineVariable<double>(
        "emg", {(long unsigned int)globalSize}, {(long unsigned int)offset}, {(long unsigned int)localSize}, adios2::ConstantDims
      ));

      // name, global size, offset, local size
      adiosFieldVariableTransmembraneFlow_ = std::make_shared<adios2::Variable<double>>(currentWriter->adiosIo->DefineVariable<double>(
        "transmembraneFlow", {(long unsigned int)globalSize}, {(long unsigned int)offset}, {(long unsigned int)localSize}, adios2::ConstantDims
      ));
    }

    offsetsVariable_ = std::make_shared<adios2::Variable<int>>(currentWriter->adiosIo->DefineVariable<int>(
      "offsets", {(long unsigned int)nRanks}, {(long unsigned int)ownRankNo}, {1}, adios2::ConstantDims
    ));
  }

  // write data to file
  int offsetUint64 = offset;
  LOG(DEBUG) << "write offset: " << offsetUint64 << ", nRanks: " << nRanks << ", ownRankNo: " << ownRankNo;
  currentWriter->engine->Put<int>(*offsetsVariable_.get(), &offsetUint64);

  // write data to file
  currentWriter->engine->Put<double>(*adiosFieldVariableVm_.get(), vmValues.data());

  if (!emgValues.empty())
  {
    // write data to file
    currentWriter->engine->Put<double>(*adiosFieldVariableEmg_.get(), emgValues.data());

    // write data to file
    currentWriter->engine->Put<double>(*adiosFieldVariableTransmembraneFlow_.get(), transmembraneFlowValues.data());
  }

}

#endif

#if defined(HAVE_MEGAMOL) && defined(HAVE_ADIOS)

void MegaMol::notifyMegaMol()
{
  if (DihuContext::ownRankNoCommWorld() == 0 && context_.zmqSocket())
  {
    std::stringstream message;
    //message << "return mmHelp()";

    std::stringstream adiosOutputFilename;
    adiosOutputFilename << lastFilename_ << ".bp";
    message << "mmSetParamValue(\"::dat::filename\", \"" << adiosOutputFilename.str() << "\")";

    std::string messageStr = message.str();
    zmq::message_t messageZmq(messageStr.length());

    memcpy(messageZmq.data(), messageStr.data(), messageStr.length());

    bool MegaMolReady = false;
    while(!MegaMolReady)
    {

      LOG(DEBUG) << "send message \"" << messageStr << "\" to MegaMol";
      assert(context_.zmqSocket());
      context_.zmqSocket()->send(messageZmq);

      LOG(DEBUG) << "recv message from MegaMol";
      zmq::message_t receivedMessage(1000);
      int replySize = context_.zmqSocket()->recv(&receivedMessage);

      if (replySize != 0)
      {
        std::string receivedMessageStr;
        receivedMessageStr = (char *)receivedMessage.data();

        LOG(DEBUG) << "received message: \"" << receivedMessageStr << "\", replySize: " << replySize;
        if (receivedMessageStr.find("READY") == std::string::npos)
        {
          LOG(DEBUG) << "MegaMol accepted request.";
          MegaMolReady = true;
        }
        else
        {
          std::this_thread::sleep_for(std::chrono::milliseconds(200));
        }
      }
      else
      {
        LOG(DEBUG) << "replySize: " << replySize;
      }
    }
  }
}
/*
MegaMol::~MegaMol()
{
  std::shared_ptr<adios_writer_t> &currentWriter = adiosWriter_[this->filenameBase_];

  currentWriter->engine->Close();
  currentWriter->engine = nullptr;

  LOG(DEBUG) << "close ADIOS engine";
}*/

#endif

} // namespace
