import sys, os
import timeit
import argparse
import importlib
import distutils.util

# parse rank arguments
rank_no = int(sys.argv[-2])
n_ranks = int(sys.argv[-1])

# add variables subfolder to python path where the variables script is located
script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_path)
sys.path.insert(0, os.path.join(script_path,'variables'))

import variables              # file variables.py, defined default values for all parameters, you can set the parameters there


# initialize some dependend variables
if variables.n_subdomains is not None:
  variables.n_subdomains_x = variables.n_subdomains[0]
  variables.n_subdomains_y = variables.n_subdomains[1]
  variables.n_subdomains_z = variables.n_subdomains[2]

variables.n_subdomains = variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z

# output information of run
if rank_no == 0:
  print("scenario_name: {},  n_subdomains: {} {} {},  n_ranks: {},  end_time: {}".format(variables.scenario_name, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, n_ranks, variables.end_time))
  print("dt_0D:         {:0.1e}, dt_1D:    {:0.1e}".format(variables.dt_0D, variables.dt_1D))
  print("dt_elasticity:           {:0.1e}".format(variables.dt_elasticity))
  print("cellml_file:             {}".format(variables.cellml_file))
  print("fiber_distribution_file: {}".format(variables.fiber_distribution_file))
  print("firing_times_file:       {}".format(variables.firing_times_file))
  print("********************************************************************************")

  print("prefactor: sigma_eff/(Am*Cm) = {} = {} / ({}*{})".format(variables.Conductivity/(variables.Am*variables.Cm), variables.Conductivity, variables.Am, variables.Cm))

# initialize all helper variables
from helper import *

variables.n_subdomains_xy = variables.n_subdomains_x * variables.n_subdomains_y
variables.n_fibers_total = variables.n_fibers_x * variables.n_fibers_y



#### set Dirichlet BC and Neumann BC for the free side of the muscle

[nx, ny, nz] = [elem + 1 for elem in variables.n_elements_muscle1]
[mx, my, mz] = [elem // 2 for elem in variables.n_elements_muscle1] # quadratic elements consist of 2 linear elements along each axis

variables.elasticity_dirichlet_bc = {}

k = 0 #free side of the muscle

for j in range(ny):
    for i in range(nx):
      variables.elasticity_dirichlet_bc[k*nx*ny + j*nx + i] = [None, None, 0.0, None,None,None] # displacement ux uy uz, velocity vx vy vz

for k in range(nz):
    variables.elasticity_dirichlet_bc[k*nx*ny] = [0.0, 0.0, None, None,None,None] # displacement ux uy uz, velocity vx vy vz

variables.elasticity_dirichlet_bc[0] = [0.0, 0.0, 0.0, None,None,None] # displacement ux uy uz, velocity vx vy vz

# meshes

# add neuron meshes
neuron_meshes = {
  "motoneuronMesh": {
    # each muscle has the same number of muscle spindels
    "nElements" :         (2*variables.n_motoneurons)-1 if n_ranks == 1 else (2*variables.n_motoneurons)*n_ranks,  # the last dof is empty in parallel
    "physicalExtent":     1, # has no special meaning. only to seperate data points in paraview
    "physicalOffset":     0,
    "logKey":             "motoneuron",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks
  },
  "muscleSpindleMesh": {
    # each muscle has the same number of muscle spindels
    "nElements" :         (2*variables.n_muscle_spindles)-1 if n_ranks == 1 else (2*variables.n_muscle_spindles)*n_ranks,
    "physicalExtent":     1, # has no special meaning. only to seperate data points in paraview
    "physicalOffset":     0,
    "logKey":             "muscle_spindle",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks
  },
  "muscle1Mesh": {
    "nElements" :         variables.n_elements_muscle1,
    "physicalExtent":     variables.muscle1_extent,
    "physicalOffset":     [0,0,0],
    "logKey":             "muscle1",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks
  },
  # needed for mechanics solver
  "muscle1Mesh_quadratic": {
    "nElements" :         [elems // 2 for elems in variables.n_elements_muscle1],
    "physicalExtent":     variables.muscle1_extent,
    "physicalOffset":     [0,0,0],
    "logKey":             "muscle1_quadratic",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks,
  }
}
variables.meshes.update(neuron_meshes)
variables.meshes.update(fiber_meshes)


# define the config dict
config = {
  "scenarioName":                   variables.scenario_name,    # scenario name which will appear in the log file
  "logFormat":                      "csv",                      # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile":   "out/" + variables.scenario_name + "/mappings_between_meshes.txt",
  "meta": {                 # additional fields that will appear in the log
    "partitioning": [variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z]
  },
  "Meshes":                variables.meshes,
  # meshes know their coordinates and the mapping happens automatically. We could also define other parameters such as a mapping tolerance
  "MappingsBetweenMeshes": {"muscle1_fiber{}".format(f) : ["muscle1Mesh", "muscle1Mesh_quadratic"] for f in range(variables.n_fibers_total)},
  "Solvers": {
    "diffusionTermSolver": {# solver for the implicit timestepping scheme of the diffusion time step
      "relativeTolerance":  variables.diffusion_solver_reltol,
      "absoluteTolerance":  1e-10,         # 1e-10 absolute tolerance of the residual
      "maxIterations":      variables.diffusion_solver_maxit,
      "solverType":         variables.diffusion_solver_type,
      "preconditionerType": variables.diffusion_preconditioner_type,
      "dumpFilename":       "",   # "out/dump_"
      "dumpFormat":         "matlab",
    },

    "mechanicsSolver": {   # solver for the dynamic mechanics problem
      "relativeTolerance":   variables.linear_relative_tolerance,           # 1e-10 relative tolerance of the linear solver
      "absoluteTolerance":   variables.linear_absolute_tolerance,           # 1e-10 absolute tolerance of the residual of the linear solver
      "solverType":          variables.elasticity_solver_type,            # type of the linear solver
      "preconditionerType":  variables.elasticity_preconditioner_type,    # type of the preconditioner
      "maxIterations":       1e4,                                         # maximum number of iterations in the linear solver
      "snesMaxFunctionEvaluations": 1e8,                                  # maximum number of function iterations
      "snesMaxIterations":   variables.snes_max_iterations,               # maximum number of iterations in the nonlinear solver
      "snesRelativeTolerance": variables.snes_relative_tolerance,         # relative tolerance of the nonlinear solver
      "snesAbsoluteTolerance": variables.snes_absolute_tolerance,         # absolute tolerance of the nonlinear solver
      "snesLineSearchType": "l2",                                         # type of linesearch, possible values: "bt" "nleqerr" "basic" "l2" "cp" "ncglinear"
      "snesRebuildJacobianFrequency": variables.snes_rebuild_jacobian_frequency,    # how often the jacobian should be recomputed, -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again
      "hypreOptions":        "",                                          # additional options for the hypre solvers could be given here
      "dumpFilename":        "",                                          # dump system matrix and right hand side after every solve
      "dumpFormat":          "matlab",                                    # default, ascii, matlab
    }
  },


  # connections of the slots, identified by slot name
  "connectedSlots": [
    # global slots only support named slots (connectedSlotsTerm1To2 also allows indices)

    # use global slot, because automatic connection of "Razumova/activestress" does not work for some reason
    # "Razumova/activestress" from CellML to Muscle contaction solver
    ("m1gout", "m1g_in"),
    ("m2gout", "m2g_in"),

    # lambda and derived values (by MapDofs) -> input of muscle splindel simulation
    ("ms0",    "ms_in0"),
    ("ms1",    "ms_in1"),
    ("ms2",    "ms_in2"),
    ("ms3",    "ms_in3"),
    ("ms4",    "ms_in4"),

    ("in_g",   "in_in"),
  ],


  "Coupling": {
    'description':            "everything",
    "timeStepWidth":          variables.dt_elasticity,
    "logTimeStepWidthAsKey":  "dt_elasticity",
    "durationLogKey":         "duration_coupling",
    "timeStepOutputInterval": 1,
    "endTime":                variables.end_time,
    "connectedSlotsTerm1To2": None,       # transfer nothing. only use numbers here!
    "connectedSlotsTerm2To1": None,       # transfer nothing back

    "Term1": {
      "MultipleCoupling": {
        "description":            "small time steps",
        "timeStepWidth":          variables.dt_neuron_system,
        "logTimeStepWidthAsKey":  "dt_multiple_coupling",
        "durationLogKey":         "duration_multiple_coupling",
        "timeStepOutputInterval": 500,
        "deferInitialization":    True,       # initialize nested solvers only right before computation of first timestep
        "connectedSlotsTerm1To2": None,       # connect lambda to slot 0 and gamma to slot 2
        "connectedSlotsTerm2To1": None,       # transfer nothing back

        # muscle spindles
        "Term1": {

          # mapping muscle spindles output -> motor neuron input
          "MapDofs": {
            "description":                "muscle_spindles_to_motoneurons",   # description that will be shown in solver structure visualization
            "nAdditionalFieldVariables":  1,                                  # number of additional field variables that are defined by this object. They have 1 component, use the templated function space and mesh given by meshName.
            "additionalSlotNames":        ["ms&in"],
            "meshName":                   "muscleSpindleMesh",                   # the mesh on which the additional field variables will be defined
            "beforeComputation": None,
            "afterComputation": [                                         # transfer/mapping of dofs that will be performed after the computation of the nested solver
              # each motoneuron gets input from all muscle spindels in the muscle
              # ms_out = primary_afferent is max(0, value) an can therefore be zero
              {
                "fromConnectorSlot":                "ms_out",
                "toConnectorSlots":                 ["ms&in"],
                "fromSlotConnectorArrayIndex":      0,                    # which fiber/compartment
                "toSlotConnectorArrayIndex":        0,
                "mode":                             "callback",           # "copyLocal", "copyLocalIfPositive", "localSetIfAboveThreshold" or "communicate"
                "fromDofNosNumbering":              "local",
                "toDofNosNumbering":                "local",
                "dofsMapping":                      None,
                "inputDofs":                        muscle1_spindle_indices,   # [0,1,...,n_muscle_spindles], this is for mesh "muscleSpindleMesh"
                "outputDofs":                       [in_ms_mesh_muscle1_motoneuron_indices],   # [0,1,...,n_motor_neurons], this is for mesh "motoneuronMesh"
                "callback":                         variables.callback_muscle_spindles_to_motoneurons,
              }
            ],

            "Heun": {
              "description":                  "muscle spindle",
              "timeStepWidth":                variables.dt_muscle_spindles,
              "logTimeStepWidthAsKey":        "dt_muscle_spindles",
              "durationLogKey":               "duration_muscle_spindles",
              "initialValues":                [],
              "timeStepOutputInterval":       500,
              "inputMeshIsGlobal":            True,
              "dirichletBoundaryConditions":  {},
              "dirichletOutputFilename":      None,                                 # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
              "checkForNanInf":               True,             # check if the solution vector contains nan or +/-inf values, if yes, an error is printed. This is a time-consuming check.
              "nAdditionalFieldVariables":    0,
              "additionalSlotNames":          [],

              # cellml model of muscle spindle
              "CellML" : {
                "modelFilename":                          variables.muscle_spindle_cellml_file,           # input C++ source file or cellml XML file
                "initializeStatesToEquilibrium":          False,                                          # if the equilibrium values of the states should be computed before the simulation starts
                "initializeStatesToEquilibriumTimestepWidth": 1e-4,                                       # if initializeStatesToEquilibrium is enable, the timestep width to use to solve the equilibrium equation

                # optimization parameters
                "optimizationType":                       "vc",                                           # "vc", "simd", "openmp" type of generated optimizated source file
                "approximateExponentialFunction":         True,                                           # if optimizationType is "vc", whether the exponential function exp(x) should be approximate by (1+x/n)^n with n=1024
                "compilerFlags":                          "-fPIC -O3 -march=native -shared ",             # compiler flags used to compile the optimized model code
                "maximumNumberOfThreads":                 0,                                              # if optimizationType is "openmp", the maximum number of threads to use. Default value 0 means no restriction.

                # stimulation callbacks, motor neuron is not stimulated by a callback function, but has a constant stimulation current
                "setSpecificStatesFunction":              None,                                           # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
                #"setSpecificStatesCallInterval":         0,                                              # set_specific_states should be called variables.stimulation_frequency times per ms, the factor 2 is needed because every Heun step includes two calls to rhs
                "setSpecificStatesCallInterval":          0,                                              # 0 means disabled
                "setSpecificStatesCallFrequency":         0,                                              # set_specific_states should be called variables.stimulation_frequency times per ms
                "setSpecificStatesFrequencyJitter":       0,                                              # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
                "setSpecificStatesRepeatAfterFirstCall":  0.01,                                           # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
                "setSpecificStatesCallEnableBegin":       0,                                              # [ms] first time when to call setSpecificStates
                "additionalArgument":                     None,

                "handleResultFunction":                   None, #handle_result,              # callback function that gets all current values and can do something with them
                "handleResultCallInterval":               1,                          # interval in which handle_result will be called
                "handleResultFunctionAdditionalParameter": None,                      # additional last argument for handle_result

                "mappings":                               variables.muscle_spindle_mappings,              # mappings between parameters and algebraics/constants and between outputConnectorSlots and states, algebraics or parameters, they are defined in helper.py
                "parametersInitialValues":                variables.muscle_spindle_parameters_initial_values,  # initial values for the parameters, either once for all instances, or for all instances in array of struct ordering with nParameters_ parameters per instance: [inst0p0, inst0p1, ... inst0pn, inst1p0, inst1p1, ...]

                "meshName":                               "muscleSpindleMesh",
                "stimulationLogFilename":                 "out/stimulation.log",

                # output writer for states, algebraics and parameters
                "OutputWriter" : [
                  {"format": "Paraview",   "outputInterval": int(2./variables.dt_muscle_spindles*variables.output_timestep_spindles), "filename": "out/" + variables.scenario_name + "/muscle_spindles", "binary": True, "fixedFormat": False, "combineFiles": True, "onlyNodalValues": True, "fileNumbering": "incremental"},
                  {"format": "PythonFile", "outputInterval": int(2./variables.dt_muscle_spindles*variables.output_timestep_spindles), "filename": "out/" + variables.scenario_name + "/muscle_spindles", "binary": True, "fixedFormat": False, "combineFiles": True, "onlyNodalValues": True, "fileNumbering": "incremental"},
                ]
              }
            }
          }
        },


        # motor neurons
        "Term2": {

          # mapping signals (from spindles and interneurons) to motor neuron + cortical input to actual inputs
          "MapDofs": {
            "description":                "motoneurons_input",   # description that will be shown in solver structure visualization
            "nAdditionalFieldVariables":  1,                              # number of additional field variables that are defined by this object. They have 1 component, use the templated function space and mesh given by meshName.
            "additionalSlotNames":        ["ms&in"],
            "meshName":                   "muscleSpindleMesh",               # the mesh on which the additional field variables will be defined
            "beforeComputation": [                                        # transfer/mapping of dofs that will be performed before the computation of the nested solver
              {
                "fromConnectorSlot":                "ms&in",
                "toConnectorSlots":                 "mn_in",
                "fromSlotConnectorArrayIndex":      0,
                "toSlotConnectorArrayIndex":        0,
                "mode":                             "callback",          # "copyLocal", "copyLocalIfPositive", "localSetIfAboveThreshold" or "communicate"
                "fromDofNosNumbering":              "local",
                "toDofNosNumbering":                "local",
                "dofsMapping":                      None,
                "inputDofs":                        in_ms_mesh_muscle1_motoneuron_indices,   # [0,1,...,n_muscle_spindles,...,n_muscle_spindles+n_interneurons]
                "outputDofs":                       [muscle1_motoneuron_indices],   # [0,1,...,n_motoneurons]
                "callback":                         variables.callback_motoneurons_input,
              }
            ],
            "afterComputation": None,

            "Heun" : {
              "description":                  "motoneurons",
              "timeStepWidth":                variables.dt_motoneuron,
              "logTimeStepWidthAsKey":        "dt_motoneuron",
              "durationLogKey":               "duration_motoneuron",
              "initialValues":                [],
              "timeStepOutputInterval":       500,
              "inputMeshIsGlobal":            True,
              "dirichletBoundaryConditions":  {},
              "dirichletOutputFilename":      None,                                 # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
              "checkForNanInf":               True,             # check if the solution vector contains nan or +/-inf values, if yes, an error is printed. This is a time-consuming check.
              "nAdditionalFieldVariables":    0,
              "additionalSlotNames":          [],

              # cellml model of motorneuron
              "CellML" : {
                "modelFilename":                          variables.motoneuron_cellml_file,               # input C++ source file or cellml XML file
                "initializeStatesToEquilibrium":          False,                                          # if the equilibrium values of the states should be computed before the simulation starts
                "initializeStatesToEquilibriumTimestepWidth": 1e-4,                                       # if initializeStatesToEquilibrium is enable, the timestep width to use to solve the equilibrium equation

                # optimization parameters
                "optimizationType":                       "vc",                                           # "vc", "simd", "openmp" type of generated optimizated source file
                "approximateExponentialFunction":         True,                                           # if optimizationType is "vc", whether the exponential function exp(x) should be approximate by (1+x/n)^n with n=1024
                "compilerFlags":                          "-fPIC -O3 -march=native -shared ",             # compiler flags used to compile the optimized model code
                "maximumNumberOfThreads":                 0,                                              # if optimizationType is "openmp", the maximum number of threads to use. Default value 0 means no restriction.

                # stimulation callbacks, motor neuron is not stimulated by a callback function, but has a constant stimulation current
                "setSpecificStatesFunction":              None,                                           # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
                #"setSpecificStatesCallInterval":         0,                                              # set_specific_states should be called variables.stimulation_frequency times per ms, the factor 2 is needed because every Heun step includes two calls to rhs
                "setSpecificStatesCallInterval":          0,                                              # 0 means disabled
                "setSpecificStatesCallFrequency":         0,                                              # set_specific_states should be called variables.stimulation_frequency times per ms
                "setSpecificStatesFrequencyJitter":       0,                                              # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
                "setSpecificStatesRepeatAfterFirstCall":  0.01,                                           # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
                "setSpecificStatesCallEnableBegin":       0,                                              # [ms] first time when to call setSpecificStates
                "additionalArgument":                     None,

                "mappings":                               variables.motoneuron_mappings,                  # mappings between parameters and algebraics/constants and between outputConnectorSlots and states, algebraics or parameters, they are defined in helper.py
                "parametersInitialValues":                variables.motoneuron_parameters_initial_values * 2, # initial values for the parameters, either once for all instances, or for all instances in array of struct ordering with nParameters_ parameters per instance: [inst0p0, inst0p1, ... inst0pn, inst1p0, inst1p1, ...]

                "meshName":                               "motoneuronMesh",                               # use the linear mesh, it was partitioned by the helper.py script which called opendihu/scripts/create_partitioned_meshes_for_settings.py
                "stimulationLogFilename":                 "out/stimulation.log",

              },
              # output writer for states, algebraics and parameters
              "OutputWriter" : [
                {"format": "Paraview",   "outputInterval": 100, "filename": "out/" + variables.scenario_name + "/motoneurons", "binary": True, "fixedFormat": False, "combineFiles": True, "onlyNodalValues": True, "fileNumbering": "incremental"},
                {"format": "PythonFile", "outputInterval": 100, "filename": "out/" + variables.scenario_name + "/motoneurons", "binary": True, "fixedFormat": False, "combineFiles": True, "onlyNodalValues": True, "fileNumbering": "incremental"},
                # {"format": "PythonCallback", "outputInterval": 1+0*int(2./variables.dt_motoneuron*variables.output_timestep_motoneuron), "callback": lambda *args, **kwargs: print('output callback')},
              ]
            }
          }
        },

        # muscle1: bidoamin + 1D monodomain + 0D
        "Term3": {

          # map from motoneuronMesh to stimulated nodes
          "MapDofs": {
            # TODO move both mn->V mappings around the Heun solver?
            "description":                "motoneurons->stimulated nodes",  # description that will be shown in solver structure visualization
            "nAdditionalFieldVariables":  1,                              # number of additional field variables that are defined by this object. They have 1 component, use the templated function space and mesh given by meshName.
            "additionalSlotNames":        ["mn_out"],
            "meshName":                   "motoneuronMesh",               # the mesh on which the additional field variables will be defined
            "afterComputation":          None,

            # mapping from motoneuronMesh which contains on every rank as many nodes as there are motoneurons to the 3D domain
            # map from motoneuronMesh (algebraics) to the fiber meshes (solution)
            "beforeComputation": [                                        # transfer/mapping of dofs that will be performed before the computation of the nested solver
              {
                "fromConnectorSlot":                "mn_out",                # source slot of the dofs mapping
                "toConnectorSlots":                 "m1vm",                # target slot of the dofs mapping
                "fromSlotConnectorArrayIndex":      0,
                "toSlotConnectorArrayIndex":        get_fiber_index_in_motor_unit(fiber_index, motor_unit_no),      # which fiber in this motor unit
                "mode":                             "localSetIfAboveThreshold",          # "copyLocal", "copyLocalIfPositive", "localSetIfAboveThreshold" or "communicate"
                "fromDofNosNumbering":              "local",
                "toDofNosNumbering":                "global",
                "dofsMapping":                      {muscle1_motoneuron_indices[motor_unit_no]: stimulation_node_nos},   # map from the motor unit to the stimulated node of the fiber mesh
                "inputDofs":                        None,                # this option is only needed in mode "callback"
                "outputDofs":                       None,                # this option is only needed in mode "callback"
                "callback":                         None,                # this option is only needed in mode "callback"
                "thresholdValue":                   20,                  # if mode is "localSetIfAboveThreshold", this is the threshold, if the value is above it, set the value `valueToSet`
                "valueToSet":                       variables.vm_value_stimulated,       # if mode is "localSetIfAboveThreshold", this is the value to set the target dof to, if the source dof is above thresholdValue.
              }
              # iterate over all motor units and all fibers in every motor unit
              for motor_unit_no in range(variables.n_motor_units)
                for fiber_index in range(get_n_fibers_in_motor_unit(motor_unit_no))
            ],

            "MultipleInstances": {
              "logKey":                     "duration_subdomains_xy_muscle1",
              "ranksAllComputedInstances":  list(range(n_ranks)),
              "nInstances":                 variables.n_subdomains_xy,
              "instances": [
                {
                  "ranks": list(range(subdomain_coordinate_y*variables.n_subdomains_x + subdomain_coordinate_x, n_ranks, variables.n_subdomains_x*variables.n_subdomains_y)),
                  "StrangSplitting": {
                    #"numberTimeSteps": 1,
                    "timeStepWidth":          variables.dt_splitting_0D1D,  # 1e-1
                    "logTimeStepWidthAsKey":  "dt_splitting",
                    "durationLogKey":         "duration_monodomain_muscle1",
                    "timeStepOutputInterval": 100,
                    "connectedSlotsTerm1To2": [0],   # transfer slot 0 = state Vm from Term1 (CellML) to Term2 (Diffusion), for elasticity also transfer gamma
                    "connectedSlotsTerm2To1": [0],   # transfer the same back, this avoids data copy

                    # CellML, i.e. reaction term of Monodomain equation
                    "Term1": {
                      "MultipleInstances": {
                        "logKey":             "duration_subdomains_z_muscle1",
                        "nInstances":         n_fibers_in_subdomain_x(subdomain_coordinate_x)*n_fibers_in_subdomain_y(subdomain_coordinate_y),
                        "instances": [
                          {
                            "ranks":                          list(range(variables.n_subdomains_z)),    # these rank nos are local nos to the outer instance of MultipleInstances, i.e. from 0 to number of ranks in z direction
                            "Heun" : {
                              "timeStepWidth":                variables.dt_0D,                         # timestep width of 0D problem
                              "logTimeStepWidthAsKey":        "dt_0D",                                 # key under which the time step width will be written to the log file
                              "durationLogKey":               "duration_0D_muscle1",                           # log key of duration for this solver
                              "timeStepOutputInterval":       1e4,                                     # how often to print the current timestep
                              "initialValues":                [],                                      # no initial values are specified
                              "dirichletBoundaryConditions":  {},                                      # no Dirichlet boundary conditions are specified
                              "dirichletOutputFilename":      None,                                    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable

                              "inputMeshIsGlobal":            True,                                    # the boundary conditions and initial values would be given as global numbers
                              "checkForNanInf":               False,                                   # abort execution if the solution contains nan or inf values
                              "nAdditionalFieldVariables":    0,                                       # number of additional field variables
                              "additionalSlotNames":          [],                                      # names for the additional slots

                              "CellML" : {
                                "modelFilename":                          variables.cellml_file,                          # input C++ source file or cellml XML file
                                #"statesInitialValues":                   [],                                             # if given, the initial values for the the states of one instance
                                "initializeStatesToEquilibrium":          False,                                          # if the equilibrium values of the states should be computed before the simulation starts
                                "initializeStatesToEquilibriumTimestepWidth": 1e-4,                                       # if initializeStatesToEquilibrium is enable, the timestep width to use to solve the equilibrium equation

                                # optimization parameters
                                "optimizationType":                       variables.optimization_type,                    # "vc", "simd", "openmp" type of generated optimizated source file
                                "approximateExponentialFunction":         variables.approximate_exponential_function,     # if optimizationType is "vc", whether the exponential function exp(x) should be approximate by (1+x/n)^n with n=1024
                                "compilerFlags":                          "-fPIC -O3 -march=native -shared ",             # compiler flags used to compile the optimized model code
                                "maximumNumberOfThreads":                 variables.maximum_number_of_threads,            # if optimizationType is "openmp", the maximum number of threads to use. Default value 0 means no restriction.
                                "useAoVSMemoryLayout":                    variables.use_aovs_memory_layout,               # if optimizationType is "vc", whether to use the Array-of-Vectorized-Struct (AoVS) memory layout instead of the Struct-of-Vectorized-Array (SoVA) memory layout. Setting to True is faster.

                                # stimulation callbacks
                                #"libraryFilename":                       "cellml_simd_lib.so",                           # compiled library
                                #"setSpecificParametersFunction":         set_specific_parameters,                        # callback function that sets parameters like stimulation current
                                #"setSpecificParametersCallInterval":     int(1./variables.stimulation_frequency/variables.dt_0D),         # set_specific_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
                                "setSpecificStatesFunction":              None,                                             # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
                                #"setSpecificStatesCallInterval":          2*int(1./variables.stimulation_frequency/variables.dt_0D),       # set_specific_states should be called variables.stimulation_frequency times per ms, the factor 2 is needed because every Heun step includes two calls to rhs
                                "setSpecificStatesCallInterval":          0,                                                               # 0 means disabled
                                "setSpecificStatesCallFrequency":         variables.get_specific_states_call_frequency(fiber_no, motor_unit_no),   # set_specific_states should be called variables.stimulation_frequency times per ms
                                "setSpecificStatesFrequencyJitter":       variables.get_specific_states_frequency_jitter(fiber_no, motor_unit_no), # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
                                "setSpecificStatesRepeatAfterFirstCall":  0.01,                                                            # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
                                "setSpecificStatesCallEnableBegin":       variables.get_specific_states_call_enable_begin(fiber_no, motor_unit_no),# [ms] first time when to call setSpecificStates
                                "additionalArgument":                     fiber_no,                                       # last argument that will be passed to the callback functions set_specific_states, set_specific_parameters, etc.

                                # parameters to the cellml model
                                "parametersInitialValues":                variables.parameters_initial_values,            #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
                                "mappings":                               variables.muscle1_mappings,                             # mappings between parameters and algebraics/constants and between outputConnectorSlots and states, algebraics or parameters, they are defined in helper.py

                                "meshName":                               "muscle1_fiber{}".format(fiber_no),                # reference to the fiber mesh
                                "stimulationLogFilename":                 "out/" + variables.scenario_name + "/stimulation_muscle1.log",                          # a file that will contain the times of stimulations
                              },
                              "OutputWriter" : [
                                {"format": "Paraview", "outputInterval": 1, "filename": "out/" + variables.scenario_name + "/muscle1_0D_states({},{})".format(fiber_in_subdomain_coordinate_x,fiber_in_subdomain_coordinate_y), "binary": True, "fixedFormat": False, "combineFiles": True}
                              ] if variables.states_output else []
                            },

                          } for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(subdomain_coordinate_y)) \
                              for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(subdomain_coordinate_x)) \
                                for fiber_no in [get_fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)] \
                                  for motor_unit_no in [get_motor_unit_no(fiber_no)]
                        ],
                      }
                    },

                    # Diffusion
                    "Term2": {
                      "MultipleInstances": {
                        "nInstances": n_fibers_in_subdomain_x(subdomain_coordinate_x)*n_fibers_in_subdomain_y(subdomain_coordinate_y),
                        "instances": [
                          {
                            "ranks":                         list(range(variables.n_subdomains_z)),    # these rank nos are local nos to the outer instance of MultipleInstances, i.e. from 0 to number of ranks in z direction
                            "ImplicitEuler": {
                              "initialValues":               [],                                      # initial values to be set in the solution vector prior to the first timestep
                              #"numberTimeSteps":            1,
                              "timeStepWidth":               variables.dt_1D,                         # timestep width for the diffusion problem
                              "timeStepWidthRelativeTolerance": 1e-10,                                # tolerance for the time step width, when to rebuild the system matrix
                              "logTimeStepWidthAsKey":       "dt_1D",                                 # key under which the time step width will be written to the log file
                              "durationLogKey":              "duration_1D_muscle1",                           # log key of duration for this solver
                              "timeStepOutputInterval":      1e4,                                     # how often to print the current timestep to console
                              "dirichletBoundaryConditions": {},                                      # old Dirichlet BC that are not used in FastMonodomainSolver: {0: -75.0036, -1: -75.0036},
                              "dirichletOutputFilename":     None,                                    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
                              "inputMeshIsGlobal":           True,                                    # initial values would be given as global numbers
                              "solverName":                  "diffusionTermSolver",                   # reference to the linear solver
                              "checkForNanInf":              False,                                   # if the solution should be checked for NaN and inf values, this requires a lot of runtimes
                              "nAdditionalFieldVariables":   2,    # number of additional field variables that should be added and potentially written to output files, these field variables can be used for receiving data from other solvers
                              "additionalSlotNames":         [],                                      # slot names for the additional field variables
                              "FiniteElementMethod" : {
                                "inputMeshIsGlobal":         True,
                                "meshName":                  "muscle1_fiber{}".format(fiber_no),
                                "prefactor":                 get_diffusion_prefactor(fiber_no, motor_unit_no),  # resolves to Conductivity / (Am * Cm)
                                "solverName":                "diffusionTermSolver",
                                "slotName":                  "",
                              },
                              "OutputWriter" : [
                              ]
                            },
                          } for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(subdomain_coordinate_y)) \
                              for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(subdomain_coordinate_x)) \
                                for fiber_no in [get_fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)] \
                                  for motor_unit_no in [get_motor_unit_no(fiber_no)]
                        ],
                      },
                    },
                  }
                } if (subdomain_coordinate_x,subdomain_coordinate_y) == (variables.own_subdomain_coordinate_x,variables.own_subdomain_coordinate_y) else None
                for subdomain_coordinate_y in range(variables.n_subdomains_y)
                    for subdomain_coordinate_x in range(variables.n_subdomains_x)
              ]
                },
                "fiberDistributionFile":    variables.fiber_distribution_file,   # for FastMonodomainSolver, e.g. MU_fibre_distribution_3780.txt
                "firingTimesFile":          variables.firing_times_file,         # for FastMonodomainSolver, e.g. MU_firing_times_real.txt
                "onlyComputeIfHasBeenStimulated": variables.fast_monodomain_solver_optimizations,                          # only compute fibers after they have been stimulated for the first time
                "disableComputationWhenStatesAreCloseToEquilibrium": variables.fast_monodomain_solver_optimizations,       # optimization where states that are close to their equilibrium will not be computed again
                "valueForStimulatedPoint":  variables.vm_value_stimulated,       # to which value of Vm the stimulated node should be set
                "neuromuscularJunctionRelativeSize": 0.1,                        # range where the neuromuscular junction is located around the center, relative to fiber length. The actual position is draws randomly from the interval [0.5-s/2, 0.5+s/2) with s being this option. 0 means sharply at the center, 0.1 means located approximately at the center, but it can vary 10% in total between all fibers.
                "generateGPUSource":        True,                                # (set to True) only effective if optimizationType=="gpu", whether the source code for the GPU should be generated. If False, an existing source code file (which has to have the correct name) is used and compiled, i.e. the code generator is bypassed. This is useful for debugging, such that you can adjust the source code yourself. (You can also add "-g -save-temps " to compilerFlags under CellMLAdapter)
                "useSinglePrecision":       False,                               # only effective if optimizationType=="gpu", whether single precision computation should be used on the GPU. Some GPUs have poor double precision performance. Note, this drastically increases the error and, in consequence, the timestep widths should be reduced.
                #"preCompileCommand":        "bash -c 'module load argon-tesla/gcc/11-20210110-openmp; module list; gcc --version",     # only effective if optimizationType=="gpu", system command to be executed right before the compilation
                #"postCompileCommand":       "'",   # only effective if optimizationType=="gpu", system command to be executed right after the compilation
              },
            }
        },
      },
    "Term2": {

      # map from λ in the 3D mesh to muscle spindles input
      "MapDofs": {
        "description":                "muscle_spindles_input",        # description that will be shown in solver structure visualization
        "nAdditionalFieldVariables":  5,                              # number of additional field variables that are defined by this object. They have 1 component, use the templated function space and mesh given by meshName.
        "additionalSlotNames":        ["ms0","ms1","ms2","ms3","ms4"],
        "meshName":                   "muscleSpindleMesh",            # the mesh on which the additional field variables will be defined
        "beforeComputation": None,
        "afterComputation": [       # transfer/mapping of dofs that will be performed before the computation of the nested solver
          # read spindle stretch (slot m_lda) and communicate to all processes
          # {
          #   "fromConnectorSlot":                "m1lda",
          #   "toConnectorSlots":                 "m1ms0",
          #   "fromSlotConnectorArrayIndex":      0,                    # which fiber/compartment
          #   "toSlotConnectorArrayIndex":        0,
          #   "mode":                             "communicate",        # "copyLocal", "copyLocalIfPositive", "localSetIfAboveThreshold" or "communicate"
          #   "fromDofNosNumbering":              "global",
          #   "toDofNosNumbering":                "global",
          #   "dofsMapping":
          #     dbg({muscle_spindle_dof : [rank_no*variables.n_muscle_spindles + i for rank_no in range(n_ranks)]
          #      for i,muscle_spindle_dof in enumerate(muscle_spindle_node_nos)}),
          #   "inputDofs":                        None,
          #   "outputDofs":                       None,
          #   "callback":                         None,
          #   #"thresholdValue":                   20,                  # if mode is "localSetIfAboveThreshold", this is the threshold, if the value is above it, set the value `valueToSet`
          #   #"valueToSet":                       20,                  # if mode is "localSetIfAboveThreshold", this is the value to set the target dof to, if the source dof is above thresholdValue.
          # },
          # call callback_muscle_spindles_input
          {
            "fromConnectorSlot":                "m1lda",
            "toConnectorSlots":                 ["ms0","ms1","ms2","ms3","ms4"],
            "fromSlotConnectorArrayIndex":      0,                    # which fiber/compartment
            "toSlotConnectorArrayIndex":        0,
            "mode":                             "callback",           # "copyLocal", "copyLocalIfPositive", "localSetIfAboveThreshold" or "communicate"
            "fromDofNosNumbering":              "local",
            "toDofNosNumbering":                "local",
            "dofsMapping":                      None,
            "inputDofs":                        muscle1_spindle_node_nos,
            "outputDofs":                       [muscle1_spindle_indices] * 5,   # dofs for each slot
            "callback":                         variables.callback_muscle_spindles_input,
          }
        ],

        "MuscleContractionSolver": {
            "numberTimeSteps":              1,                         # only use 1 timestep per interval
            "timeStepOutputInterval":       1,
            "Pmax":                         variables.Pmax,            # maximum PK2 active stress
            "enableForceLengthRelation":    True,                      # if the factor f_l(λ_f) modeling the force-length relation (as in Heidlauf2013) should be multiplied. Set to false if this relation is already considered in the CellML model.
            "lambdaDotScalingFactor":       1.0,                       # scaling factor for the output of the lambda dot slot, i.e. the contraction velocity. Use this to scale the unit-less quantity to, e.g., micrometers per millisecond for the subcellular model.
            "slotNames":                    ["m1lda", "m1ldot", "m1g_in", "m1T", "m1ux", "m1uy", "m1uz"],  # slot names of the data connector slots: lambda, lambdaDot, gamma, traction
            "OutputWriter" : [
              {"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/" + variables.scenario_name + "/muscle1_contraction", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles": True, "fileNumbering": "incremental"},
            ],
            "mapGeometryToMeshes":          ["muscle1Mesh"] + [key for key in fiber_meshes.keys() if "muscle1_fiber" in key],    # the mesh names of the meshes that will get the geometry transferred
            "reverseMappingOrder":          True,                      # if the mapping target->own mesh should be used instead of own->target mesh. This gives better results in some cases.
            "dynamic":                      variables.dynamic,                      # if the dynamic solid mechanics solver should be used, else it computes the quasi-static problem

            # the actual solid mechanics solver, this is either "DynamicHyperelasticitySolver" or "HyperelasticitySolver", depending on the value of "dynamic"
            "DynamicHyperelasticitySolver": {
              "timeStepWidth":              variables.dt_elasticity,           # time step width
              "durationLogKey":             "muscle1_duration_mechanics",               # key to find duration of this solver in the log file
              "timeStepOutputInterval":     1,                         # how often the current time step should be printed to console

              "materialParameters":         variables.material_parameters,  # material parameters of the Mooney-Rivlin material
              "density":                    variables.rho,             # density of the material
              "dampingFactor":              variables.damping_factor,  # factor for velocity dependent damping
              "displacementsScalingFactor": 1.0,                       # scaling factor for displacements, only set to sth. other than 1 only to increase visual appearance for very small displacements
              "residualNormLogFilename":    "out/"+variables.scenario_name+"/muscle1_log_residual_norm.txt",   # log file where residual norm values of the nonlinear solver will be written
              "useAnalyticJacobian":        True,                      # whether to use the analytically computed jacobian matrix in the nonlinear solver (fast)
              "useNumericJacobian":         False,                     # whether to use the numerically computed jacobian matrix in the nonlinear solver (slow), only works with non-nested matrices, if both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian

              "dumpDenseMatlabVariables":   False,                     # whether to have extra output of matlab vectors, x,r, jacobian matrix (very slow)
              # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables all all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian

              # mesh
              "inputMeshIsGlobal":          True,                     # boundary conditions and initial values are given as global numbers (every process has all information)
              "meshName":                   "muscle1Mesh_quadratic",       # name of the 3D mesh, it is defined under "Meshes" at the beginning of this config
              "fiberMeshNames":             [],                       # fiber meshes that will be used to determine the fiber direction
              "fiberDirection":             [0,0,1],                  # if fiberMeshNames is empty, directly set the constant fiber direction, in element coordinate system

              # solving
              "solverName":                 "mechanicsSolver",         # name of the nonlinear solver configuration, it is defined under "Solvers" at the beginning of this config
              #"loadFactors":                [0.5, 1.0],                # load factors for every timestep
              "loadFactors":                [],                        # no load factors, solve problem directly
              "loadFactorGiveUpThreshold":  0.25,                       # a threshold for the load factor, when to abort the solve of the current time step. The load factors are adjusted automatically if the nonlinear solver diverged. If the load factors get too small, it aborts the solve.
              "scaleInitialGuess":          False,                     # when load stepping is used, scale initial guess between load steps a and b by sqrt(a*b)/a. This potentially reduces the number of iterations per load step (but not always).
              "nNonlinearSolveCalls":       1,                         # how often the nonlinear solve should be repeated

              # boundary and initial conditions
              "dirichletBoundaryConditions": variables.elasticity_dirichlet_bc,   # the initial Dirichlet boundary conditions that define values for displacements u and velocity v
              "neumannBoundaryConditions":   variables.elasticity_neumann_bc,     # Neumann boundary conditions that define traction forces on surfaces of elements
              "divideNeumannBoundaryConditionValuesByTotalArea": True,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
              "updateDirichletBoundaryConditionsFunction": None,                  # muscle1_update_dirichlet_boundary_conditions_helper, function that updates the dirichlet BCs while the simulation is running
              "updateDirichletBoundaryConditionsFunctionCallInterval": 1,         # every which step the update function should be called, 1 means every time step
              "updateNeumannBoundaryConditionsFunction":   None,                    # function that updates the Neumann BCs while the simulation is running
              "updateNeumannBoundaryConditionsFunctionCallInterval": 1,           # every which step the update function should be called, 1 means every time step


              "initialValuesDisplacements":  [[0.0,0.0,0.0] for _ in range(variables.n_points_global)],     # the initial values for the displacements, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
              "initialValuesVelocities":     [[0.0,0.0,0.0] for _ in range(variables.n_points_global)],     # the initial values for the velocities, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
              "extrapolateInitialGuess":     True,                                # if the initial values for the dynamic nonlinear problem should be computed by extrapolating the previous displacements and velocities
              "constantBodyForce":           [0,0,0],       # a constant force that acts on the whole body, e.g. for gravity

              "dirichletOutputFilename":     "out/"+variables.scenario_name+"/muscle1_dirichlet_boundary_conditions",     # output filename for the dirichlet boundary conditions, set to "" to have no output
              "totalForceLogFilename":       "out/"+variables.scenario_name+"/muscle1_tendon_force.csv",              # filename of a log file that will contain the total (bearing) forces and moments at the top and bottom of the volume
              "totalForceLogOutputInterval":       10,  
              "OutputWriter" : [
                #{"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/" + variables.scenario_name + "/muscle1_contraction", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles": True, "fileNumbering": "incremental"},
              ],                                # output interval when to write the totalForceLog file
              "pressure": {   # output files for pressure function space (linear elements), contains pressure values, as well as displacements and velocities
                # "OutputWriter" : [
                #   {"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/"+variables.scenario_name+"/muscle1_pressure", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
                # ]
              },
              # 3. additional output writer that writes virtual work terms
              "dynamic": {    # output of the dynamic solver, has additional virtual work values
                # "OutputWriter" : [   # output files for displacements function space (quadratic elements)
                #   {"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/muscle1_dynamic", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
                #   {"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/"+variables.scenario_name+"/muscle1_virtual_work", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
                # ],
              },
              # 4. output writer for debugging, outputs files after each load increment, the geometry is not changed but u and v are written
              "LoadIncrements": {
                # "OutputWriter" : [
                #   {"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/"+variables.scenario_name+"/muscle1_load_increments", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
                # ]
              }  
            }   
        },          
      }
    }
  }
  
}
