import numpy as np
import sys

# 3D laplace problem

nx = 1
ny = 1
nz = 4

local = False

# Dirichlet boundary conditions
bc = {}

if not local:
  for j in range(int(ny+1)):
    for i in range(int(nx+1)):
      # x = 0 plane
      x = i/nx
      y = j/ny
      index = int(j*(nx+1) + i)
    
      # z- plane (bottom)
      bc[index] = np.sin(x*np.pi)
      bc[index] = 1.0
      
      # z+ plane (top)
      index += int(nz*(nx+1)*(ny+1))
      bc[index] = np.sin(y*np.pi) + 2.0
      bc[index] = 2.0

rank_no = (int)(sys.argv[-2])

n_elements = [1,1,4]
if local:
  n_elements = [1,1,1]
  if rank_no == 0:
    n_elements = [1,1,2]

  # boundary conditions
  bc = {}
  if rank_no == 0:
    bc = {dof:1.0 for dof in range(4)}
  elif rank_no == 2:
    bc = {-1-dof:2.0 for dof in range(4)}

config = {
  "OutputSurface": {
    "face": "2-",
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/surface", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},      
      {"format": "PythonFile", "filename": "out/surface", "outputInterval": 1, "binary":False, "onlyNodalValues":True}
    ],
    "FiniteElementMethod" : {
      "nElements": n_elements,
      "inputMeshIsGlobal": not local,
      "physicalExtent": [1.0, 1.0, 3.0],
      "outputInterval": 1.0,
      "prefactor": 1,
      "dirichletBoundaryConditions": bc,
      "relativeTolerance": 1e-15,
      "solverType": "gmres",
      "preconditionerType": "none",
      "maxIterations": 10000,
      "OutputWriter" : [
        {"format": "Paraview", "outputInterval": 1, "filename": "out/laplace", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},      
        {"format": "PythonFile", "filename": "out/laplace", "outputInterval": 1, "binary":False, "onlyNodalValues":True}
      ]
    }
  }
}
