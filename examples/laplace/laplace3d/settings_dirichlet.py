import numpy as np
import sys

# 3D laplace problem

nx = 1
ny = 1
nz = 4

local = True

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
  "logFormat":                      "csv",     # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "FiniteElementMethod" : {
    "nElements": n_elements,
    "nRanks": [1,1,1],
    "inputMeshIsGlobal": not local,
    "physicalExtent": [1.0, 1.0, 3.0],
    "outputInterval": 1.0,
    
    "dirichletBoundaryConditions": bc,
    "neumannBoundaryConditions": [],
    "prefactor": 1,
    
    "solverType": "gmres",
    "preconditionerType": "none",
    "relativeTolerance": 1e-15,
    "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual        
    "maxIterations": 10000,
    "dumpFormat": "default",
    "dumpFilename": "",
    
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/laplace", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},      
      {"format": "PythonFile", "filename": "out/laplace", "outputInterval": 1, "binary":False, "onlyNodalValues":True, "fileNumbering": "incremental"}
    ]
  },
}
