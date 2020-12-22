# Mooney-Rivlin transversely isotropic with parameters from Heidlauf 2013

import numpy as np
import sys, os

# number of elements
nx = 2    # 2
ny = 1    # 2
nz = 5    # 5

# number of nodes
mx = 2*nx + 1
my = 2*ny + 1
mz = 2*nz + 1

# boundary conditions (for quadratic elements)
dirichlet_bc = {}

xpos = 0.0
ypos = 0.0
zpos = 0.0

own_rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

# bottom plane
for j in range(0,my):
  for i in range(0,mx):
    dirichlet_bc[j*(mx) + i] = [None,None,zpos]

if False:
  # left plane
  for k in range(0,mz):
    for j in range(0,my):
      dirichlet_bc[k*mx*my + j*mx + 0] = [xpos,None,None]

  # front plane
  for k in range(0,mz):
    for i in range(0,mx):
      dirichlet_bc[k*mx*my + i] = [None,ypos,None]

  # vertical edge
  for k in range(0,mz):
    dirichlet_bc[k*mx*my] = [xpos,ypos,None]

if True:
  # horizontal edge y=0
  for i in range(0,mx):
    dirichlet_bc[i] = [None,ypos,zpos]

if True:
  # horizontal edge, x=0
  for j in range(0,my):
    dirichlet_bc[j*(mx)] = [xpos,None,zpos]

# corner
dirichlet_bc[0] = [xpos,ypos,zpos]
dirichlet_bc[1] = [None,ypos,zpos]

# force to the top
neumann_bc = [{"element": (nz-1)*nx*ny + j*nx + i, "constantVector": [0,0,10], "face": "2+"} for j in range(ny) for i in range(nx)]

#dirichlet_bc = {}
#neumann_bc = []

# fiber directions
fiber_meshes = {}
fiber_mesh_names = []

# each fiber will have mz nodes
n_elements_fiber = mz-1

for j in range(my):
  for i in range(mx):
    fiber_no = j*nx + i
    
    # determine start position of fiber in (x,y)-plane
    x = 0.5 + i
    y = -2.0 + 0.2 + j
    angle = 40./180.*np.pi
    
    # loop over points of a single fiber
    node_positions = []
    for z in range(mz):
        
      h = 5.0*z/nz
      x_pos = x
      y_pos = y + np.sin(angle)*h
      z_pos = 0.0 + np.cos(angle)*h
      node_positions.append([x_pos,y_pos,z_pos])
    
    mesh_name = "fiber{}".format(fiber_no)
    fiber_mesh_names.append(mesh_name)
    
    fiber_meshes[mesh_name] = {
      "nodePositions": node_positions,
      "nElements": [n_elements_fiber],
      "inputMeshIsGlobal": True,
      "nRanks": [n_ranks],
    }
    
config = {
  "scenarioName": "3d_box",
  "logFormat":    "csv",     # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile":   "mappings_between_meshes.txt",   # log file for mappings between meshes
  
  "Meshes": fiber_meshes,
  "HyperelasticitySolver": {
    "durationLogKey": "nonlinear",
    
    "materialParameters": [2, 3, 4, 5],  # c1, c2, b1, d1
    "displacementsScalingFactor": 10.0,   # scaling factor for displacements
    "constantBodyForce":          [0.0, 0.0, 0.0],   # body force in whole body region
    "residualNormLogFilename": "log_residual_norm.txt",
    "useAnalyticJacobian": True,
    "useNumericJacobian": False,   # Only works in parallel execution. If both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
      
    "dumpDenseMatlabVariables": False,   # extra output of matlab vectors, x,r, jacobian matrix
    # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables are all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
    
    "fiberMeshNames": fiber_mesh_names,   # fiber meshes that will be used to determine the fiber direction
    
    # mesh
    "nElements": [nx, ny, nz],
    "inputMeshIsGlobal": True,
    "physicalExtent": [2.0, 2.0, 5.0],
    "physicalOffset": [0, 0, 0],        # offset/translation where the whole mesh begins
    
    # nonlinear solver
    "relativeTolerance": 1e-10,         # 1e-10 relative tolerance of the linear solver
    "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual in the linear solver
    "solverType": "preonly",            # type of the linear solver: cg groppcg pipecg pipecgrr cgne nash stcg gltr richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs fbcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly qcg bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres tsirm cgls
    "preconditionerType": "lu",         # type of the preconditioner
    "maxIterations": 1e4,               # maximum number of iterations in the linear solver
    "dumpFilename": "",#"out/m",            # filename for output of solver matrix
    "dumpFormat": "matlab",             # default, ascii, matlab
    "snesMaxFunctionEvaluations": 1e8,  # maximum number of function iterations
    "snesMaxIterations": 15,             # maximum number of iterations in the nonlinear solver
    "snesRebuildJacobianFrequency": 5,  # frequency with which the jacobian is newly computed
    "snesRelativeTolerance": 1e-5,      # relative tolerance of the nonlinear solver
    "snesLineSearchType": "l2",         # type of linesearch, possible values: "bt" "nleqerr" "basic" "l2" "cp" "ncglinear"
    "snesAbsoluteTolerance": 1e-5,      # absolute tolerance of the nonlinear solver
    "loadFactorGiveUpThreshold": 0.1,   # if the adaptive time stepping produces a load factor smaller than this value, the solution will be accepted for the current timestep, even if it did not converge fully to the tolerance
    
    #"loadFactors":  [0.1, 0.2, 0.35, 0.5, 1.0],   # load factors for every timestep
    "loadFactors": [],                  # no load factors, solve problem directly
    "nNonlinearSolveCalls": 1,          # how often the nonlinear solve should be repeated
    
    # boundary conditions
    "dirichletBoundaryConditions": dirichlet_bc,
    "neumannBoundaryConditions": neumann_bc,
    "divideNeumannBoundaryConditionValuesByTotalArea": True,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
    #"updateDirichletBoundaryConditionsFunction": update_dirichlet_boundary_conditions,
    "updateDirichletBoundaryConditionsFunction": None,
    "updateDirichletBoundaryConditionsFunctionCallInterval": 1,
    
    "dirichletOutputFilename":  "out/dirichlet_boundary_conditions",    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values
    
    "OutputWriter" : [   # output files for displacements function space (quadratic elements)
      {"format": "Paraview", "outputInterval": 1, "filename": "out/u", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      {"format": "PythonFile", "filename": "out/u", "outputInterval": 1, "binary":False, "onlyNodalValues":True, "fileNumbering": "incremental"},
    ],
    "pressure": {   # output files for pressure function space (linear elements)
      "OutputWriter" : [
        {"format": "Paraview", "outputInterval": 1, "filename": "out/p", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
        {"format": "PythonFile", "filename": "out/p", "outputInterval": 1, "binary":False, "onlyNodalValues":True, "fileNumbering": "incremental"},
      ]
    },
    # output writer for debugging, outputs files after each load increment, the geometry is not changed but u and v are written
    "LoadIncrements": {   
      "OutputWriter" : [
        {"format": "Paraview", "outputInterval": 1, "filename": "out/load_increments", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      ]
    },
  },
}
