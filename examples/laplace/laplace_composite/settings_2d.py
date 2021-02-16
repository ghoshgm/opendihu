import numpy as np

nx1 = 3
nx2 = 2
ny = 2

# boundary conditions (for quadratic elements)
bc = {}
n_nodes_x1 = 2*nx1+1
for i in range(int(n_nodes_x1)):
  x = i/n_nodes_x1
  bc[i] = np.sin(x*np.pi)
  
n_nodes_x2 = 2*nx2+1
for i in range(int(n_nodes_x2)):
  x = i/n_nodes_x2
  i2 = (2*nx1+1)*(2*ny+1) - 3 + n_nodes_x2*(2*ny) + i
  bc[i2] = np.sin(x*np.pi)
  
  print("bc[{}] = {}, bc[{}] = {}".format(i, bc[i], i2, bc[i2]))

config = {
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "logFormat":                      "csv",                      # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "scenarioName":                   "laplace_composite",        # scenario name to find the run in the log file
  "mappingsBetweenMeshesLogFile":   "",                         # a log file about mappings between meshes, here we do not want that because there are no mappings
  "Meshes": {
    "submesh0": {
      "nElements": [nx1, ny],
      "inputMeshIsGlobal": True,
      "physicalExtent": [2.0, 1.0],
      "physicalOffset": [0.0, 0.0],
    },
    "submesh1": {
      "nElements": [nx2, ny],
      "inputMeshIsGlobal": True,
      "physicalExtent": [2.0, 1.0],
      "physicalOffset": [2.0, 0.5],
    },
  },
  "FiniteElementMethod": {
    "inputMeshIsGlobal": True,
    "outputInterval": 1.0,
    
    "dirichletBoundaryConditions": bc,
    "dirichletOutputFilename":     None,                                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    "neumannBoundaryConditions": [],
    "prefactor": [1,2],
    
    "meshName": ["submesh0", "submesh1"],
    
    # solver parameters
    "solverType": "gmres",
    "preconditionerType": "none",
    "relativeTolerance": 1e-15,
    "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual        
    "maxIterations": 10000,
    "dumpFormat": "default",
    "dumpFilename": "",
    "slotName": "",
    
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/2d", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      {"format": "PythonFile", "filename": "out/2d", "outputInterval": 1, "binary":False, "onlyNodalValues":True, "fileNumbering": "incremental"},
    ]
  },
}
