# # Modified diffusion 1D
# # Compute
# import numpy as np
#
# nx = 10
# ny = 10
#
# # boundary conditions (for linear elements)
# # set bottom and top border to 1
# bc = {}
# for i in range(int(nx+1)):
#   x = i/(nx+1.)
#   #bc[i] = np.sin(x*np.pi)
#   bc[i] = 1
#   i2 = (nx+1)*ny + i
#   #bc[i2] = np.sin(x*np.pi)
#   bc[i2] = 1
#
#   print("bc[{}] = {}, bc[{}] = {}".format(i, bc[i], i2, bc[i2]))
#
# config = {
#   "PinT": {        # this is the name of the solver, as given in the constructor to the timestepping object
#     "myOption": 42,                   # example option that is parsed in the constructor
#     "option1": "blabla",              # another example option that is parsed in the data object
#
#     # settings for the nested solver
#     "FiniteElementMethod" : {
#       # mesh
#       "nElements": [nx, ny],               # number of elements in x and y direction
#       "physicalExtent": [1.0, 1.0],        # the physical size of the domain
#
#       # solver
#       "solverType": "gmres",          # the solver type, refer to PETSc documentation about implemented solvers and preconditioners (https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPType.html)
#       "preconditionerType": "none",   # the preconditioner type (https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html#PCType)
#       "relativeTolerance": 1e-15,     # the relative tolerance of the linear solver
#       "maxIterations": 1e4,           # maximum number of iterations of the linear solver
#       "dumpFilename": "out/",             # a filename to dump the system matrix and rhs in every timestep, "" to disable
#       "dumpFormat": "matlab",         # the file format of the dump file, "matlab", "ascii" or "default"
#
#       # problem definition
#       "prefactor": 1,                      # prefactor of the laplace operator
#       "dirichletBoundaryConditions": bc,   # Dirichlet boundary conditions as dict
#       "neumannBoundaryConditions": [{"element": j*nx, "constantVector": 1.0, "face": "0-"} for j in range(ny)],   # Neumann boundary conditions: at border
#       "inputMeshIsGlobal": True,           # boundary condition values are given for all dofs, even if executed in parallel
#     },
#
#     "OutputWriter" : [
#       {"format": "Paraview", "outputInterval": 1, "filename": "out/paraview", "binary": True, "fixedFormat": False, "onlyNodalValues": True, "combineFiles": False},
#       {"format": "PythonFile", "filename": "out/python", "outputInterval": 1, "binary":False, "onlyNodalValues": True}
#     ]
#   }
# }

# Modified diffusion 1D
# Compute
import numpy as np

n = 5   # number of elements

config = {
 "Solvers": {
   "linearSolver": {
     "solverType": "gmres",          # the solver type, refer to PETSc documentation about implemented solvers and preconditioners (https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPType.html)
     "preconditionerType": "none",   # the preconditioner type (https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html#PCType)
     "relativeTolerance": 1e-15,     # the relative tolerance of the linear solver
     "maxIterations": 1e4,           # maximum number of iterations of the linear solver
     "dumpFilename": "",             # a filename to dump the system matrix and rhs in every timestep, "" to disable
     "dumpFormat": "matlab",         # the file format of the dump file, "matlab", "ascii" or "default"
   },
},
# for k in range (n):
"Meshes":{
    "mesh_{}".format(k): {
      "nElements": (2 ** k)-1,                 # number of elements
      # "nElements": 5,                 # number of elements
      "physicalExtent": 4.0,          # the physical size of the domain
      "inputMeshIsGlobal": True,      # boundary conditions are given as global indices
    }
  for k in range(n)
  },
"PinT": {        # this is the name of the solver, as given in the constructor to the timestepping object
    "myOption": 42,                   # example option that is parsed in the constructor
    "option1": "blabla",              # another example option that is parsed in the data object

    "TimeSteppingScheme": [
      {"ImplicitEuler": {
            "numberTimeSteps": 5,
            "endTime": 0.1,
            "initialValues": [2,2,4,5,2,2],    # the initial values
            "dirichletBoundaryConditions": {}, # Dirichlet boundary conditions as dict
            "inputMeshIsGlobal": True,         # initial values and BC's are given for all dofs, even if executed in parallel
            "timeStepOutputInterval": 1,       # how often to print the current timestep to console
            "nAdditionalFieldVariables": 0,    # for more complex nested solvers, the number of additional field variables that will be transferred without being touched
            "solverName": "linearSolver",      # the solver to use, referes to what was defined under "Solvers"

            "FiniteElementMethod" : {
               # mesh
               "meshName": "mesh_{}".format(k),
               #"nElements": n,                 # number of elements
               #"physicalExtent": 4.0,          # the physical size of the domain
               "solverName": "linearSolver",   # the solver to use, referes to what was defined under "Solvers"
               "prefactor": 5.0,               # the prefactor 'c' of 'du/dt = c du^2/dx^2'
               "inputMeshIsGlobal": True,      # boundary conditions are given as global indices
               "nAdditionalFieldVariables": 0, # for more complex nested solvers, the number of additional field variables that will be transferred without being touched
            },
            "OutputWriter" : [
              #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binaryOutput": "false", "fixedFormat": False, "onlyNodalValues": True},
              {"format": "PythonFile", "filename": "out/diffusion1d_implicit", "outputInterval": 1, "binary":False, "onlyNodalValues": True}
            ]
          }
        } for k in range(n)],
 # for j in range(n)
 }
}
