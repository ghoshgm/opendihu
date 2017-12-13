#include <iostream>
#include <cstdlib>
#include <fstream>

#include "gtest/gtest.h"
#include "opendihu.h"
#include "utility/petsc_utility.h"
#include "arg.h"
#include "stiffness_matrix_tester.h"

namespace SpatialDiscretization
{
  
TEST(PoissonTest, MatrixIsCorrect1DSmall)
{
  std::string pythonConfig = R"(
# Laplace 1D
n = 5
    
# boundary conditions
bc = {}
bc[0] = 1.0
bc[n] = 0.0

config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "FiniteElementMethod" : {
    "nElements": n,
    "physicalExtend": 4.0,
    "DirichletBoundaryCondition": bc,
    "relativeTolerance": 1e-15,
  },
  "OutputWriter" : [
    {"format": "Paraview", "interval": 1, "filename": "out", "binaryOutput": "false", "fixedFormat": False},
    #{"format": "Python", "filename": "p"}
  ]
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  FiniteElementMethod<
    Mesh::RegularFixed<1>,
    BasisFunction::Lagrange,
    Integrator::None,
    Equation::Static::Poisson
  > equationDiscretized(settings);
  
  Computation computation(settings, equationDiscretized);
  computation.run();

  std::vector<double> referenceMatrix = {
    1, 0, 0, 0, 0, 0,
    0, -2.5, 1.25, 0, 0, 0,
    0, 1.25, -2.5, 1.25, 0, 0,
    0, 0, 1.25, -2.5, 1.25, 0, 
    0, 0, 0, 1.25, -2.5, 0,
    0, 0, 0, 0, 0, 1
  };
  std::vector<double> referenceRhs = {
    1, -1.25, 0, 0, 0, 0
  };
  std::map<int, double> dirichletBC = {{0, 1.0}, {5,0.0}};
  StiffnessMatrixTester::compareMatrix(equationDiscretized, referenceMatrix);
  StiffnessMatrixTester::compareRhs(equationDiscretized, referenceRhs);
  StiffnessMatrixTester::checkDirichletBCInSolution(equationDiscretized, dirichletBC);
}

TEST(PoissonTest, MatrixIsCorrect1DBig)
{
  std::string pythonConfig = R"(
# Laplace 1D
n = 10
    
# boundary conditions
bc = {}
bc[0] = 1.0
bc[n] = 2.0

config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "FiniteElementMethod" : {
    "nElements": n,
    "physicalExtend": 0.5,
    "DirichletBoundaryCondition": bc,
    "relativeTolerance": 1e-15,
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  FiniteElementMethod<
    Mesh::RegularFixed<1>,
    BasisFunction::Lagrange,
    Integrator::None,
    Equation::Static::Poisson
  > equationDiscretized(settings);
  
  Computation computation(settings, equationDiscretized);
  computation.run();

  std::vector<double> referenceMatrix = {
    1,    0,    0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, -40, 20, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 20, -40, 20, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 20, -40, 20, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 20, -40, 20, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 20, -40, 20, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 20, -40, 20, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 20, -40, 20, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 20, -40, 20, 0,
    0, 0, 0, 0, 0, 0, 0, 0,    20, -40, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
  };
  std::vector<double> referenceRhs = {
    1, -20, 0, 0, 0, 0, 0, 0, 0, -40, 2
  };
  std::map<int, double> dirichletBC = {{0, 1.0}, {10,2.0}};
  StiffnessMatrixTester::compareMatrix(equationDiscretized, referenceMatrix);
  StiffnessMatrixTester::compareRhs(equationDiscretized, referenceRhs);
  StiffnessMatrixTester::checkDirichletBCInSolution(equationDiscretized, dirichletBC);
}

TEST(PoissonTest, MatrixIsCorrect1DStencils)
{
  std::string pythonConfig = R"(
# Laplace 1D
n = 10
    
# no boundary conditions
config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "FiniteElementMethod" : {
    "nElements": n,
    "physicalExtend": n,
    "relativeTolerance": 1e-15,
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  FiniteElementMethod<
    Mesh::RegularFixed<1>,
    BasisFunction::Lagrange,
    Integrator::None,
    Equation::Static::Poisson
  > equationDiscretized(settings);
  
  Computation computation(settings, equationDiscretized);
  computation.run();

  std::vector<double> referenceMatrix(121, 0.0);
  std::array<double, 2> stencil = {-1.0, 1.0};
    // stencil for -Δu in 1D: [1 _-2_ 1] (element contribution: [_-1_ 1])
  
  // fill with stencil values
  auto matrixIndex = [](int x, int y){return y*11 + x;};
  
  // loop over elements
  for (int i=0; i<10; i++)
  {
    // nodes:
    // 0--1
    int node0 = i;
    int node1 = i+1;
    
    // add contribution from node 0 to all nodes
    referenceMatrix[matrixIndex(node0, node0)] += stencil[0];
    referenceMatrix[matrixIndex(node0, node1)] += stencil[1];
    
    // add contribution from node 1 to all nodes
    referenceMatrix[matrixIndex(node1, node0)] += stencil[1];
    referenceMatrix[matrixIndex(node1, node1)] += stencil[0];
  }
  
  StiffnessMatrixTester::compareMatrix(equationDiscretized, referenceMatrix);
}

TEST(PoissonTest, MatrixIsCorrect2DSmall)
{
  std::string pythonConfig = R"(
# Laplace 2D
print "start"
n = 2.
m = 3.

# boundary conditions
bc = {}
for i in range(int(n+1)):
  x = i/(n+1.)
  bc[i] = x
  i2 = (n+1)*m + i
  bc[int(i2)] = x
  
  print "bc[{}] = {}, bc[{}] = {}".format(i, bc[i], i2, bc[i2])

config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "FiniteElementMethod" : {
    "nElements": [n, m],
    "physicalExtend": [6.0, 9.0],
    "DirichletBoundaryCondition": bc,
    "relativeTolerance": 1e-15,
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  FiniteElementMethod<
    Mesh::RegularFixed<2>,
    BasisFunction::Lagrange,
    Integrator::None,
    Equation::Static::Poisson
  > equationDiscretized(settings);
  
  Computation computation(settings, equationDiscretized);
  computation.run();

  std::vector<double> referenceMatrix = {
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0,-12/9.,3/9., 0,1.5/9.,3/9., 0, 0, 0, 0,
    0, 0, 0,3/9.,-24/9.,3/9.,3/9.,3/9., 3/9., 0, 0, 0,
    0, 0, 0, 0,3/9.,-12/9.,0,3/9.,1.5/9., 0, 0, 0,
    0, 0, 0,1.5/9.,3/9.,0,-12/9.,3/9., 0, 0, 0, 0,
    0, 0, 0,3/9.,3/9.,3/9.,3/9.,-24/9., 3/9., 0, 0, 0,
    0, 0, 0,0,3/9.,1.5/9.,0,3/9.,-12/9., 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
  };
  std::vector<double> referenceRhs = {
    0, 1./3, 2./3, -1./9, -3./9, -2./9, -1./9, -3/9., -2./9, 0, 1./3, 2./3
  };
  std::map<int, double> dirichletBC = {{0, 0.0}, {1,1./3}, {2,2./3}, {9,0}, {10,1./3}, {11,2./3}};
  
  StiffnessMatrixTester::compareMatrix(equationDiscretized, referenceMatrix);
  StiffnessMatrixTester::compareRhs(equationDiscretized, referenceRhs);
  StiffnessMatrixTester::checkDirichletBCInSolution(equationDiscretized, dirichletBC);
}

TEST(PoissonTest, MatrixIsCorrect2DStencils)
{
  std::string pythonConfig = R"(
# Laplace 2D
config = {
  "disablePrinting": False,
  "disableMatrixPrinting": True,
  "FiniteElementMethod" : {
    "nElements": [4, 4],
    "physicalExtend": [4.0, 4.0],
    "relativeTolerance": 1e-15,
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  FiniteElementMethod<
    Mesh::RegularFixed<2>,
    BasisFunction::Lagrange,
    Integrator::None,
    Equation::Static::Poisson
  > equationDiscretized(settings);
  
  Computation computation(settings, equationDiscretized);
  computation.run();

  std::vector<double> referenceMatrix(625, 0.0);
  
  // stencil for -Δu in 2D:     [1  1   1] (element contribution: [  1/6  1/3])
  //                        1/3*[1 _-8_ 1]                        [_-2/3_ 1/6]
  //                            [1  1   1]
  double stencil[2][2] = {{-2./3, 1./6}, {1./6, 1./3}};
  
  // fill with stencil values
  auto matrixIndex = [](int x, int y){return y*25 + x;};
  
  // loop over elements
  for (int x=0; x<4; x++)
  {
    for (int y=0; y<4; y++)
    {
      // nodes:
      // 2--3
      // 0--1
      int node0 = y*5 + x;
      int node1 = y*5 + x+1;
      int node2 = (y+1)*5 + x;
      int node3 = (y+1)*5 + x+1;
      
      // add contribution from node 0 to all nodes
      referenceMatrix[matrixIndex(node0, node0)] += stencil[0][0];
      referenceMatrix[matrixIndex(node0, node1)] += stencil[0][1];
      referenceMatrix[matrixIndex(node0, node2)] += stencil[1][0];
      referenceMatrix[matrixIndex(node0, node3)] += stencil[1][1];
      
      // add contribution from node 1 to all nodes
      referenceMatrix[matrixIndex(node1, node0)] += stencil[0][1];
      referenceMatrix[matrixIndex(node1, node1)] += stencil[0][0];
      referenceMatrix[matrixIndex(node1, node2)] += stencil[1][1];
      referenceMatrix[matrixIndex(node1, node3)] += stencil[1][0];
      
      // add contribution from node 2 to all nodes
      referenceMatrix[matrixIndex(node2, node0)] += stencil[1][0];
      referenceMatrix[matrixIndex(node2, node1)] += stencil[1][1];
      referenceMatrix[matrixIndex(node2, node2)] += stencil[0][0];
      referenceMatrix[matrixIndex(node2, node3)] += stencil[0][1];
      
      // add contribution from node 3 to all nodes
      referenceMatrix[matrixIndex(node3, node0)] += stencil[1][1];
      referenceMatrix[matrixIndex(node3, node1)] += stencil[1][0];
      referenceMatrix[matrixIndex(node3, node2)] += stencil[0][1];
      referenceMatrix[matrixIndex(node3, node3)] += stencil[0][0];
    }
  }
    
  StiffnessMatrixTester::compareMatrix(equationDiscretized, referenceMatrix);
}

TEST(PoissonTest, MatrixIsCorrect3DStencils)
{
  std::string pythonConfig = R"(
# Laplace 3D
config = {
  "disablePrinting": False,
  "disableMatrixPrinting": True,
  "FiniteElementMethod" : {
    "nElements": [4, 4, 4],
    "physicalExtend": [4.0, 4.0, 4.0],
    "relativeTolerance": 1e-15,
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  FiniteElementMethod<
    Mesh::RegularFixed<3>,
    BasisFunction::Lagrange,
    Integrator::None,
    Equation::Static::Poisson
  > equationDiscretized(settings);
  
  Computation computation(settings, equationDiscretized);
  computation.run();

  std::vector<double> referenceMatrix(15625, 0.0);
  
  const double stencil[2][2][2] = {
    {{-4./12, 0./12},
    {0./12, 1./12}},    //center
    {{0./12, 1./12},
    {1./12, 1./12}},    //bottom
  };
  
  // fill with stencil values
  auto matrixIndex = [](int i, int j){return j*125+i;};
  
  // loop over elements
  for (int x=0; x<4; x++)
  {
    for (int y=0; y<4; y++)
    {
      for (int z=0; z<4; z++)
      {
        // nodes:
        // 6--7
        // 4--5
        //
        // 2--3
        // 0--1
        int node0 = z*25 + y*5 + x;
        int node1 = z*25 + y*5 + x+1;
        int node2 = z*25 + (y+1)*5 + x;
        int node3 = z*25 + (y+1)*5 + x+1;
        int node4 = (z+1)*25 + y*5 + x;
        int node5 = (z+1)*25 + y*5 + x+1;
        int node6 = (z+1)*25 + (y+1)*5 + x;
        int node7 = (z+1)*25 + (y+1)*5 + x+1;
        
        // add contribution from node 0 to all nodes          x  y  z
        referenceMatrix[matrixIndex(node0, node0)] += stencil[0][0][0];
        referenceMatrix[matrixIndex(node0, node1)] += stencil[0][1][0];
        referenceMatrix[matrixIndex(node0, node2)] += stencil[1][0][0];
        referenceMatrix[matrixIndex(node0, node3)] += stencil[1][1][0];
        referenceMatrix[matrixIndex(node0, node4)] += stencil[0][0][1];
        referenceMatrix[matrixIndex(node0, node5)] += stencil[0][1][1];
        referenceMatrix[matrixIndex(node0, node6)] += stencil[1][0][1];
        referenceMatrix[matrixIndex(node0, node7)] += stencil[1][1][1];
        
        // add contribution from node 1 to all nodes
        referenceMatrix[matrixIndex(node1, node0)] += stencil[0][1][0];
        referenceMatrix[matrixIndex(node1, node1)] += stencil[0][0][0];
        referenceMatrix[matrixIndex(node1, node2)] += stencil[1][1][0];
        referenceMatrix[matrixIndex(node1, node3)] += stencil[1][0][0];
        referenceMatrix[matrixIndex(node1, node4)] += stencil[0][1][1];
        referenceMatrix[matrixIndex(node1, node5)] += stencil[0][0][1];
        referenceMatrix[matrixIndex(node1, node6)] += stencil[1][1][1];
        referenceMatrix[matrixIndex(node1, node7)] += stencil[1][0][1];
        
        // add contribution from node 2 to all nodes
        referenceMatrix[matrixIndex(node2, node0)] += stencil[1][0][0];
        referenceMatrix[matrixIndex(node2, node1)] += stencil[1][1][0];
        referenceMatrix[matrixIndex(node2, node2)] += stencil[0][0][0];
        referenceMatrix[matrixIndex(node2, node3)] += stencil[0][1][0];
        referenceMatrix[matrixIndex(node2, node4)] += stencil[1][0][1];
        referenceMatrix[matrixIndex(node2, node5)] += stencil[1][1][1];
        referenceMatrix[matrixIndex(node2, node6)] += stencil[0][0][1];
        referenceMatrix[matrixIndex(node2, node7)] += stencil[0][1][1];
        
        // add contribution from node 3 to all nodes
        referenceMatrix[matrixIndex(node3, node0)] += stencil[1][1][0];
        referenceMatrix[matrixIndex(node3, node1)] += stencil[1][0][0];
        referenceMatrix[matrixIndex(node3, node2)] += stencil[0][1][0];
        referenceMatrix[matrixIndex(node3, node3)] += stencil[0][0][0];
        referenceMatrix[matrixIndex(node3, node4)] += stencil[1][1][1];
        referenceMatrix[matrixIndex(node3, node5)] += stencil[1][0][1];
        referenceMatrix[matrixIndex(node3, node6)] += stencil[0][1][1];
        referenceMatrix[matrixIndex(node3, node7)] += stencil[0][0][1];
        
        // add contribution from node 4 to all nodes          x  y  z
        referenceMatrix[matrixIndex(node4, node0)] += stencil[0][0][1];
        referenceMatrix[matrixIndex(node4, node1)] += stencil[0][1][1];
        referenceMatrix[matrixIndex(node4, node2)] += stencil[1][0][1];
        referenceMatrix[matrixIndex(node4, node3)] += stencil[1][1][1];
        referenceMatrix[matrixIndex(node4, node4)] += stencil[0][0][0];
        referenceMatrix[matrixIndex(node4, node5)] += stencil[0][1][0];
        referenceMatrix[matrixIndex(node4, node6)] += stencil[1][0][0];
        referenceMatrix[matrixIndex(node4, node7)] += stencil[1][1][0];
        
        // add contribution from node 5 to all nodes
        referenceMatrix[matrixIndex(node5, node0)] += stencil[0][1][1];
        referenceMatrix[matrixIndex(node5, node1)] += stencil[0][0][1];
        referenceMatrix[matrixIndex(node5, node2)] += stencil[1][1][1];
        referenceMatrix[matrixIndex(node5, node3)] += stencil[1][0][1];
        referenceMatrix[matrixIndex(node5, node4)] += stencil[0][1][0];
        referenceMatrix[matrixIndex(node5, node5)] += stencil[0][0][0];
        referenceMatrix[matrixIndex(node5, node6)] += stencil[1][1][0];
        referenceMatrix[matrixIndex(node5, node7)] += stencil[1][0][0];
        
        // add contribution from node 6 to all nodes
        referenceMatrix[matrixIndex(node6, node0)] += stencil[1][0][1];
        referenceMatrix[matrixIndex(node6, node1)] += stencil[1][1][1];
        referenceMatrix[matrixIndex(node6, node2)] += stencil[0][0][1];
        referenceMatrix[matrixIndex(node6, node3)] += stencil[0][1][1];
        referenceMatrix[matrixIndex(node6, node4)] += stencil[1][0][0];
        referenceMatrix[matrixIndex(node6, node5)] += stencil[1][1][0];
        referenceMatrix[matrixIndex(node6, node6)] += stencil[0][0][0];
        referenceMatrix[matrixIndex(node6, node7)] += stencil[0][1][0];
        
        // add contribution from node 7 to all nodes
        referenceMatrix[matrixIndex(node7, node0)] += stencil[1][1][1];
        referenceMatrix[matrixIndex(node7, node1)] += stencil[1][0][1];
        referenceMatrix[matrixIndex(node7, node2)] += stencil[0][1][1];
        referenceMatrix[matrixIndex(node7, node3)] += stencil[0][0][1];
        referenceMatrix[matrixIndex(node7, node4)] += stencil[1][1][0];
        referenceMatrix[matrixIndex(node7, node5)] += stencil[1][0][0];
        referenceMatrix[matrixIndex(node7, node6)] += stencil[0][1][0];
        referenceMatrix[matrixIndex(node7, node7)] += stencil[0][0][0];
      }
    }
  }
    
  StiffnessMatrixTester::compareMatrix(equationDiscretized, referenceMatrix);
}

TEST(PoissonTest, RhsIsParsedCorrectly)
{
  std::string pythonConfig1 = R"(
# Laplace 2D
config = {
  "disablePrinting": False,
  "disableMatrixPrinting": True,
  "FiniteElementMethod" : {
    "nElements": [4, 5],
    "physicalExtend": [1.0, 2.0, 3.0],
    "rightHandSide": [1, 4.0, 5, 9.0, 0.0, 0, 5, 7, 3],
    "DirichletBoundaryCondition": {0:0},
    "relativeTolerance": 1e-15,
  },
}
)";

  std::string pythonConfig2 = R"(
# Laplace 2D
config = {
  "disablePrinting": False,
  "disableMatrixPrinting": True,
  "FiniteElementMethod" : {
    "nElements": [4, 5],
    "physicalExtend": [1.0, 2.0],
    "rightHandSide": {2:5, 7:7.0, 3:9, "8":3, 0:1, 1:4, 6.0:5, "10":0},
    "DirichletBoundaryCondition": {0:0},
    "relativeTolerance": 1e-15,
  },
}
)";

  DihuContext settings1(argc, argv, pythonConfig1);
  FiniteElementMethod<
    Mesh::RegularFixed<2>,
    BasisFunction::Lagrange,
    Integrator::None,
    Equation::Static::Poisson
  > equationDiscretized1(settings1);
  
  Computation computation1(settings1, equationDiscretized1);
  computation1.run();
  
  
  DihuContext settings2(argc, argv, pythonConfig2);
  FiniteElementMethod<
    Mesh::RegularFixed<2>,
    BasisFunction::Lagrange,
    Integrator::None,
    Equation::Static::Poisson
  > equationDiscretized2(settings2);
  
  Computation computation2(settings2, equationDiscretized2);
  computation2.run();
  
  StiffnessMatrixTester::checkEqual(equationDiscretized1, equationDiscretized2);
}

TEST(PoissonTest, RhsDiscretizationMatrix1DIsCorrect)
{
  
  std::string pythonConfig = R"(
# Laplace 1D
config = {
  "disablePrinting": False,
  "disableMatrixPrinting": True,
  "FiniteElementMethod" : {
    "nElements": [10],
    "physicalExtend": [4.0],
    "rightHandSide": [1,8,3,4,5,6,7,8,9,10,11,12],
    "relativeTolerance": 1e-15,
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  FiniteElementMethod<
    Mesh::RegularFixed<1>,
    BasisFunction::Lagrange,
    Integrator::None,
    Equation::Static::Poisson
  > equationDiscretized(settings);
  
  Computation computation(settings, equationDiscretized);
  computation.run();
  
  FiniteElementMethod<
    Mesh::RegularFixed<1>,
    BasisFunction::Lagrange,
    Integrator::None,
    Equation::Dynamic::Diffusion
  > equationDiscretized2(settings);
  
  Computation computation2(settings, equationDiscretized2);
  computation2.run();

  std::vector<double> initialRhs = {
    1,8,3,4,5,6,7,8,9,10,11,12
  };
  
  StiffnessMatrixTester::testRhsDiscretizationMatrix(equationDiscretized, equationDiscretized2, initialRhs);

}

TEST(PoissonTest, DeformableRhsDiscretizationMatrix1DIsCorrect)
{
  
  std::string pythonConfig = R"(
# Laplace 1D
config = {
  "disablePrinting": False,
  "disableMatrixPrinting": True,
  "FiniteElementMethod" : {
    "nElements": [10],
    "physicalExtend": [4.0],
    "rightHandSide": [1,8,3,4,5,6,7,8,9,10,11,12],
    "relativeTolerance": 1e-15,
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  FiniteElementMethod<
    Mesh::Deformable<1>,
    BasisFunction::Lagrange,
    Integrator::Gauss<2>,
    Equation::Static::Poisson
  > equationDiscretized(settings);
  
  Computation computation(settings, equationDiscretized);
  computation.run();
  
  FiniteElementMethod<
    Mesh::RegularFixed<1>,
    BasisFunction::Lagrange,
    Integrator::None,
    Equation::Dynamic::Diffusion
  > equationDiscretized2(settings);
  
  Computation computation2(settings, equationDiscretized2);
  computation2.run();

  std::vector<double> initialRhs = {
    1,8,3,4,5,6,7,8,9,10,11,12
  };
  
  StiffnessMatrixTester::testRhsDiscretizationMatrix(equationDiscretized, equationDiscretized2, initialRhs);

}

TEST(PoissonTest, RhsDiscretizationMatrix2DIsCorrect)
{
  
  std::string pythonConfig = R"(
# Laplace 2D
config = {
  "disablePrinting": False,
  "disableMatrixPrinting": True,
  "FiniteElementMethod" : {
    "nElements": [3, 2],
    "physicalExtend": [4.0, 3.0],
    "rightHandSide": [1,2,3,4,5,6,7,8,9,10,11,12],
    "relativeTolerance": 1e-15,
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  FiniteElementMethod<
    Mesh::RegularFixed<2>,
    BasisFunction::Lagrange,
    Integrator::None,
    Equation::Static::Poisson
  > equationDiscretized(settings);
  
  Computation computation(settings, equationDiscretized);
  computation.run();
  
  FiniteElementMethod<
    Mesh::RegularFixed<2>,
    BasisFunction::Lagrange,
    Integrator::None,
    Equation::Dynamic::Diffusion
  > equationDiscretized2(settings);
  
  Computation computation2(settings, equationDiscretized2);
  computation2.run();

  std::vector<double> initialRhs = {
    1,2,3,4,5,6,7,8,9,10,11,12
  };
  
  StiffnessMatrixTester::testRhsDiscretizationMatrix(equationDiscretized, equationDiscretized2, initialRhs);

}

TEST(PoissonTest, DeformableRhsDiscretizationMatrix2DIsCorrect)
{
  
  std::string pythonConfig = R"(
# Laplace 2D
config = {
  "disablePrinting": False,
  "disableMatrixPrinting": True,
  "FiniteElementMethod" : {
    "nElements": [3, 2],
    "physicalExtend": [4.0, 3.0],
    "rightHandSide": [1,2,3,4,5,6,7,8,9,10,11,12],
    "relativeTolerance": 1e-15,
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  FiniteElementMethod<
    Mesh::RegularFixed<2>,
    BasisFunction::Lagrange,
    Integrator::None,
    Equation::Static::Poisson
  > equationDiscretized(settings);
  
  Computation computation(settings, equationDiscretized);
  computation.run();
  
  FiniteElementMethod<
    Mesh::Deformable<2>,
    BasisFunction::Lagrange,
    Integrator::Gauss<2>,
    Equation::Dynamic::Diffusion
  > equationDiscretized2(settings);
  
  Computation computation2(settings, equationDiscretized2);
  computation2.run();

  std::vector<double> initialRhs = {
    1,2,3,4,5,6,7,8,9,10,11,12
  };
  
  StiffnessMatrixTester::testRhsDiscretizationMatrix(equationDiscretized, equationDiscretized2, initialRhs);

}

TEST(PoissonTest, RhsDiscretizationMatrix3DIsCorrect)
{
  
  std::string pythonConfig = R"(
# Laplace 3D
config = {
  "disablePrinting": False,
  "disableMatrixPrinting": True,
  "FiniteElementMethod" : {
    "nElements": [3, 2, 2],
    "physicalExtend": [4.0, 3.0, 2.0],
    "rightHandSide": [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24],
    "relativeTolerance": 1e-15,
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  FiniteElementMethod<
    Mesh::RegularFixed<3>,
    BasisFunction::Lagrange,
    Integrator::None,
    Equation::Static::Poisson
  > equationDiscretized(settings);
  
  std::cout<<"Note: This should print an error message about mesh being non-uniform.";
  Computation computation(settings, equationDiscretized);
  computation.run();
  
  FiniteElementMethod<
    Mesh::RegularFixed<3>,
    BasisFunction::Lagrange,
    Integrator::None,
    Equation::Dynamic::Diffusion
  > equationDiscretized2(settings);
  
  Computation computation2(settings, equationDiscretized2);
  computation2.run();

  std::vector<double> initialRhs = {
    1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24
  };
  
  StiffnessMatrixTester::testRhsDiscretizationMatrix(equationDiscretized, equationDiscretized2, initialRhs);

}

TEST(PoissonTest, DeformableRhsDiscretizationMatrix3DIsCorrect)
{
  
  std::string pythonConfig = R"(
# Laplace 3D
config = {
  "disablePrinting": False,
  "disableMatrixPrinting": True,
  "FiniteElementMethod" : {
    "nElements": [3, 2, 2],
    "physicalExtend": [4.0, 3.0, 2.0],
    "rightHandSide": [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24],
    "relativeTolerance": 1e-15,
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  FiniteElementMethod<
    Mesh::RegularFixed<3>,
    BasisFunction::Lagrange,
    Integrator::None,
    Equation::Static::Poisson
  > equationDiscretized(settings);
  
  std::cout<<"Note: This should print an error message about mesh being non-uniform.";
  Computation computation(settings, equationDiscretized);
  computation.run();
  
  FiniteElementMethod<
    Mesh::Deformable<3>,
    BasisFunction::Lagrange,
    Integrator::Gauss<2>,
    Equation::Dynamic::Diffusion
  > equationDiscretized2(settings);
  
  Computation computation2(settings, equationDiscretized2);
  computation2.run();

  std::vector<double> initialRhs = {
    1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24
  };
  
  StiffnessMatrixTester::testRhsDiscretizationMatrix(equationDiscretized, equationDiscretized2, initialRhs);

}

};

