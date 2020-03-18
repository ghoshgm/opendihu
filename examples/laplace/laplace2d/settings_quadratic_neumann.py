# Laplace 2D, quadratic, Neumann BC

nx = 3
ny = nx

neumann_bc = []

# bottom, top
if True:
  for i in range(nx):
    neumann_bc += [
        {"element": i, "constantValue": -1.0, "face": "1-"},               # bottom
        {"element": (ny-1)*nx + i, "constantValue": 1.0, "face": "1+"},   # top
      ]

# left, right
if False:
  for j in range(ny):
    neumann_bc += [
        {"element": j*nx, "constantValue": -1, "face": "0-"},            # left
        {"element": j*nx + (nx-1), "constantValue": 1, "face": "0+"},   # top
      ]
      
config = {
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "FiniteElementMethod" : {
    "nElements": [nx, ny],
    "inputMeshIsGlobal": True,
    "physicalExtent": [1.0, 1.0],
    "nodePositions": [[float(i)/(2*nx), float(j)/(2*ny)] for j in range(2*ny+1) for i in range(2*nx+1)],
    "elements": [[2*j*(2*nx+1)+2*i, 2*j*(2*nx+1)+2*i+1, 2*j*(2*nx+1)+2*i+2, (2*j+1)*(2*nx+1)+2*i, (2*j+1)*(2*nx+1)+2*i+1, (2*j+1)*(2*nx+1)+2*i+2, (2*j+2)*(2*nx+1)+2*i, (2*j+2)*(2*nx+1)+2*i+1, (2*j+2)*(2*nx+1)+2*i+2] for j in range(ny) for i in range(nx)],
    
    #"dirichletBoundaryConditions": {
    #  0:0, 1:0.2, 2:0.4, 3:0.6, 4:0.8, 5:1,   
    #  30:0, 31:0.2, 32:0.4, 33:0.6, 34:0.8, 35:1
    #},
    "dirichletBoundaryConditions": {0:0},
    "neumannBoundaryConditions": neumann_bc,
    "prefactor": 1,
    
    "solverType": "gmres",
    "preconditionerType": "none",
    "relativeTolerance": 1e-15,
    "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual        
    "maxIterations": 10000,
    "dumpFormat": "default",
    "dumpFilename": "",
    
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/p", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":False},
      {"format": "PythonFile", "filename": "out/p", "outputInterval": 1, "binary":False, "onlyNodalValues":True},
    ]
  }
}

print(config)
print(len(config["FiniteElementMethod"]["nodePositions"]))
