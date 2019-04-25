#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
createSeedPoints(const std::array<bool,4> &subdomainIsAtBorder, int seedPointsZIndex, const std::vector<Vec3> &nodePositions, std::vector<Vec3> &seedPoints)
{
  // nodePositions contains all node positions in the current 3D mesh
  // naming of borders: (0-,1-,0+,1+)
  //     ___1+__
  //    |   |   |
  // 0- |___|___| 0+
  // ^  |   |   |
  // |  |___|___|
  // +-->   1-

  LOG(DEBUG) << "createSeedPoints, seedPointsZIndex: " << seedPointsZIndex << ", subdomainIsAtBorder: " << std::boolalpha << subdomainIsAtBorder;

  int subdomainNNodesX = nBorderPointsXNew_;
  int subdomainNNodesY = nBorderPointsXNew_;

  // boundary indices for face0Minus and face0Plus (vertical direction)
  int iBeginVertical = 0;
  int iEndVertical = nBorderPointsXNew_;

  if (subdomainIsAtBorder[(int)Mesh::face_t::face1Minus])
    iBeginVertical += 1;

  if (subdomainIsAtBorder[(int)Mesh::face_t::face1Plus])
    iEndVertical -= 1;

  // boundary indices for face1Minus and face1Plus (horizontal direction)
  int iBeginHorizontal = 0;
  int iEndHorizontal = nBorderPointsXNew_;

  if (subdomainIsAtBorder[(int)Mesh::face_t::face0Minus])
    iBeginHorizontal += 1;

  if (subdomainIsAtBorder[(int)Mesh::face_t::face0Plus])
    iEndHorizontal -= 1;

  // face0Minus
  if (!subdomainIsAtBorder[(int)Mesh::face_t::face0Minus])
  {
    for (int i = iBeginVertical; i < iEndVertical; i++)
    {
      seedPoints.push_back(nodePositions[seedPointsZIndex*subdomainNNodesX*subdomainNNodesY + i*subdomainNNodesX + 0]);
    }
  }

  // face0Plus
  if (!subdomainIsAtBorder[(int)Mesh::face_t::face0Plus])
  {
    for (int i = iBeginVertical; i < iEndVertical; i++)
    {
      seedPoints.push_back(nodePositions[seedPointsZIndex*subdomainNNodesX*subdomainNNodesY + i*subdomainNNodesX + (subdomainNNodesX-1)]);
    }
  }

  // face1Minus (with corner points)
  if (!subdomainIsAtBorder[(int)Mesh::face_t::face1Minus])
  {
    for (int i = iBeginHorizontal; i < iEndHorizontal; i++)
    {
      seedPoints.push_back(nodePositions[seedPointsZIndex*subdomainNNodesX*subdomainNNodesY + i]);
    }
  }

  // face1Plus (with corner points)
  if (!subdomainIsAtBorder[(int)Mesh::face_t::face1Plus])
  {
    for (int i = iBeginHorizontal; i < iEndHorizontal; i++)
    {
      seedPoints.push_back(nodePositions[seedPointsZIndex*subdomainNNodesX*subdomainNNodesY + (subdomainNNodesY-1)*subdomainNNodesX + i]);
    }
  }


  LOG(DEBUG) << "seedPoints: starting with horizontal center line, streamlineIndex = " << seedPoints.size();

  // horizontal center line (with corner points)
  for (int i = iBeginHorizontal; i < iEndHorizontal; i++)
  {
    seedPoints.push_back(nodePositions[seedPointsZIndex*subdomainNNodesX*subdomainNNodesY + int(subdomainNNodesY/2)*subdomainNNodesX + i]);
  }

  LOG(DEBUG) << "seedPoints: starting with vertical center line, streamlineIndex = " << seedPoints.size();

  // vertical center line (with corner points and center point)
  for (int i = iBeginVertical; i < iEndVertical; i++)
  {
    seedPoints.push_back(nodePositions[seedPointsZIndex*subdomainNNodesX*subdomainNNodesY + i*subdomainNNodesX + int(subdomainNNodesX/2)]);
  }

  LOG(DEBUG) << "seedPoints: end, streamlineIndex = " << seedPoints.size();


  // create seed points at corners for border points
  //     ___1+__
  //    |   |   |
  // 0- |___|___| 0+
  // ^  |   |   |
  // |  |___|___|
  // +-->   1-

  std::vector<std::pair<int,int>> coordinates =
  {
    std::pair<int,int>{1,1},    // bottom left
    std::pair<int,int>{subdomainNNodesX-2, 1},    // bottom right
    std::pair<int,int>{1, subdomainNNodesY-2},    // top left
    std::pair<int,int>{subdomainNNodesX-2, subdomainNNodesY-2}    // top right
  };

  for (std::vector<std::pair<int,int>>::iterator iter = coordinates.begin(); iter != coordinates.end(); iter++)
  {
    int i = iter->first;
    int j = iter->second;
    LOG(DEBUG) << "(i,j) = " << i << "," << j << ", index: " << seedPointsZIndex*subdomainNNodesX*subdomainNNodesY + j*subdomainNNodesX + i
      << ", size of nodePositions: " << nodePositions.size();
    seedPoints.push_back(nodePositions[seedPointsZIndex*subdomainNNodesX*subdomainNNodesY + j*subdomainNNodesX + i]);
  }
  LOG(DEBUG) << "seedPoints: after cornerStreamlines, streamlineIndex = " << seedPoints.size();

}

} // namespace
