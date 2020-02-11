# isotropic Mooney Rivlin
import numpy as np
import sys, os

# number of elements
nx = 2    # 2
ny = 2    # 2
nz = 5    # 5

# boundary conditions (for quadratic elements)
dirichlet_bc = {}

xpos = 0.0
ypos = 0.0
zpos = 0.0

# bottom plane
for j in range(0,2*ny+1):
  for i in range(0,2*nx+1):
    dirichlet_bc[j*(2*nx+1) + i] = [None,None,zpos]

if False:
  # left plane
  for k in range(0,2*nz+1):
    for j in range(0,2*ny+1):
      dirichlet_bc[k*(2*nx+1)*(2*ny+1) + j*(2*nx+1)] = [xpos,None,None]

  # front plane
  for k in range(0,2*nz+1):
    for i in range(0,2*nx+1):
      dirichlet_bc[k*(2*nx+1)*(2*ny+1) + i] = [None,ypos,None]

  # vertical edge
  for k in range(0,2*nz+1):
    dirichlet_bc[k*(2*nx+1)*(2*ny+1)] = [xpos,ypos,None]

if True:
  # horizontal edge y=0
  for i in range(0,2*nx+1):
    dirichlet_bc[i] = [None,ypos,zpos]

if True:
  # horizontal edge, x=0
  for j in range(0,2*ny+1):
    dirichlet_bc[j*(2*nx+1)] = [xpos,None,zpos]

# corner
dirichlet_bc[0] = [xpos,ypos,zpos]
dirichlet_bc[1] = [None,ypos,zpos]

neumann_bc = [{"element": (nz-1)*nx*ny + j*nx + i, "constantVector": [0,1e-1,5e-1], "face": "2+"} for j in range(ny) for i in range(nx)]

#dirichlet_bc = {}
neumann_bc = []

config = {
  "scenarioName": "3d_box",
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "HyperelasticitySolver": {
    "durationLogKey": "nonlinear",
    
    #"materialParameters": [1.5,2.0],
    "materialParameters": [0.0,1.0],
    "displacementsScalingFactor": 1.0,   # scaling factor for displacements
    "constantBodyForce": [5e-2, 0.0, 0.0],   # body force in whole body region
    "residualNormLogFilename": "log_residual_norm.txt",
    "useAnalyticJacobian": True,
    "useNumericJacobian": False,   # Only works in parallel execution. If both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
      
    "dumpDenseMatlabVariables": False,   # extra output of matlab vectors, x,r, jacobian matrix
    # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables are all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
    
    # mesh
    "nElements": [nx, ny, nz],
    "inputMeshIsGlobal": True,
    "physicalExtent": [nx, ny, nz],
    
    # solver
    "relativeTolerance": 1e-10,
    "solverType": "preonly",          # cg groppcg pipecg pipecgrr cgne nash stcg gltr richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs fbcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly qcg bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres tsirm cgls
    "preconditionerType": "lu",
    "maxIterations": 1e4,
    "dumpFilename": "out/m",
    "dumpFormat": "matlab",   # default, ascii, matlab
    
    "dirichletBoundaryConditions": dirichlet_bc,
    "neumannBoundaryConditions": neumann_bc,
    
    "OutputWriter" : [   # output files for displacements function space (quadratic elements)
      {"format": "Paraview", "outputInterval": 1, "filename": "out/u", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},
      {"format": "PythonFile", "filename": "out/u", "outputInterval": 1, "binary":False, "onlyNodalValues":True},
    ],
    "pressure": {   # output files for pressure function space (linear elements)
      "OutputWriter" : [
        {"format": "Paraview", "outputInterval": 1, "filename": "out/p", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},
        {"format": "PythonFile", "filename": "out/p", "outputInterval": 1, "binary":False, "onlyNodalValues":True},
      ]
    }
  },
}
