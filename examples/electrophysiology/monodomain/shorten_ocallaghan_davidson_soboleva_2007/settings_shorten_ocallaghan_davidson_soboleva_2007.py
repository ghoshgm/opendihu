# Electrophysiology
# Monodomain with Shorten model as the rhs. 
#
# parameters: [<scenario_name> [<dt_0D> [<dt_1D>]]]
#
# Plot in the "out" folder:
# plot strang_000*

import sys
import numpy as np

end_time = 10   # [ms] end time of simulation
n_elements = 200
element_size = 1./100   # [cm]
#element_size = 1./10

# global parameters
Conductivity = 3.828    # sigma, conductivity [mS/cm]
Am = 500.0              # surface area to volume ratio [cm^-1]
Cm = 0.58               # membrane capacitance [uF/cm^2]
innervation_zone_width = 0.  # [cm]
solver_type = "gmres"

# timing parameters
stimulation_frequency = 1000*1e-3    # [ms^-1] sampling frequency of stimuli in firing_times_file, in stimulations per ms, number before 1e-3 factor is in Hertz.
dt_0D = 1e-4                     # timestep width of ODEs
dt_1D = 1e-4                     # timestep width of diffusion
dt_splitting = max(dt_0D,dt_1D)  # overall timestep width of splitting

output_timestep = 1e0            # timestep for output files

# input files
#cellml_file = "../../../input/shorten_ocallaghan_davidson_soboleva_2007.c"
#cellml_file = "../../../input/shorten_ocallaghan_davidson_soboleva_2007.cellml"
cellml_file = "../../../input/shorten_ocallaghan_davidson_soboleva_2007_no_stim.cellml"

fiber_distribution_file = "../../../input/MU_fibre_distribution_3780.txt"
#firing_times_file = "../../../input/MU_firing_times_real.txt"
firing_times_file = "../../../input/MU_firing_times_immediately.txt"
firing_times_file = "../../../input/MU_firing_times_once.txt"
#firing_times_file = "../../../input/MU_firing_times_always.txt"

# parse command line options (scenario name)
scenario_name = ""
if len(sys.argv) <= 2:
  scenario_name = ""

if len(sys.argv) >= 2+1:
  scenario_name = sys.argv[0]

if len(sys.argv) >= 2+2:
  dt_0D = float(sys.argv[1])

if len(sys.argv) >= 2+3:
  dt_1D = float(sys.argv[2])
  dt_splitting = max(dt_0D,dt_1D)  # overall timestep width of splitting


rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

if rank_no == 0:
  print("scenario_name: {}".format(scenario_name))

  print("n elements: {}, end time: {}".format(n_elements,end_time))
  print("prefactor: ",Conductivity/(Am*Cm))
  print("dt_0D: {}, dt_1D: {}, dt_splitting: {}".format(dt_0D, dt_1D, dt_splitting))

# set values for cellml model
if "shorten_ocallaghan_davidson_soboleva_2007_no_stim" in cellml_file:
  # parameters: stimulation current I_stim, fiber stretch λ
  mappings = {
    ("parameter", 0):           "wal_environment/I_HH",    # parameter 0 is stimulation current
    ("parameter", 1):           "razumova/l_hs",           # parameter 1 is fiber stretch λ
    ("connectorSlot", 0): "wal_environment/vS",      # expose state 0 = Vm to the operator splitting
    ("connectorSlot", 1): "razumova/activestress",   
  }
  parameters_initial_values = [0.0, 1.0]
  nodal_stimulation_current = 400.
  
  
# load MU distribution and firing times
fiber_distribution = np.genfromtxt(fiber_distribution_file, delimiter=" ")
firing_times = np.genfromtxt(firing_times_file)

def get_motor_unit_no(fiber_no):
  """
  get the no. of the motor unit which fiber fiber_no is part of
  """
  return int(fiber_distribution[fiber_no % len(fiber_distribution)]-1)

def fiber_gets_stimulated(fiber_no, frequency, current_time):
  """
  determine if fiber fiber_no gets stimulated at simulation time current_time
  """

  # determine motor unit
  alpha = 1.0   # 0.8
  mu_no = (int)(get_motor_unit_no(fiber_no)*alpha)

  # determine if fiber fires now
  index = int(np.round(current_time * frequency))
  n_firing_times = np.size(firing_times,0)

  #if firing_times[index % n_firing_times, mu_no] == 1:
  #  print("{}: fiber {} is mu {}, t = {}, row: {}, stimulated: {} {}".format(rank_no, fiber_no, mu_no, current_time, (index % n_firing_times), firing_times[index % n_firing_times, mu_no], "true" if firing_times[index % n_firing_times, mu_no] == 1 else "false"))

  return firing_times[index % n_firing_times, mu_no] == 1

# callback function that can set states, i.e. prescribed values for stimulation
def set_specific_states(n_nodes_global, time_step_no, current_time, states, fiber_no):

  #print("call set_specific_states at time {}".format(current_time))

  # determine if fiber gets stimulated at the current time
  is_fiber_gets_stimulated = fiber_gets_stimulated(fiber_no, stimulation_frequency, current_time)

  if is_fiber_gets_stimulated:
    # determine nodes to stimulate (center node, left and right neighbour)
    innervation_zone_width_n_nodes = innervation_zone_width*100  # 100 nodes per cm
    innervation_node_global = int(n_nodes_global / 2)  # + np.random.randint(-innervation_zone_width_n_nodes/2,innervation_zone_width_n_nodes/2+1)
    nodes_to_stimulate_global = [innervation_node_global]
    if innervation_node_global > 0:
      nodes_to_stimulate_global.insert(0, innervation_node_global-1)
    if innervation_node_global < n_nodes_global-1:
      nodes_to_stimulate_global.append(innervation_node_global+1)
    if rank_no == 0:
      print("t: {}, stimulate fiber {} at nodes {}".format(current_time, fiber_no, nodes_to_stimulate_global))

    for node_no_global in nodes_to_stimulate_global:
      states[(node_no_global,0,0)] = 20.0   # key: ((x,y,z),nodal_dof_index,state_no)

# callback function that can set parameters, i.e. stimulation current
def set_specific_parameters(n_nodes_global, time_step_no, current_time, parameters, fiber_no):
  
  # determine if fibre gets stimulated at the current time
  is_fiber_gets_stimulated = fiber_gets_stimulated(fiber_no, stimulation_frequency, current_time)
  
  # determine nodes to stimulate (center node, left and right neighbour)
  innervation_zone_width_n_nodes = innervation_zone_width*100  # 100 nodes per cm
  innervation_node_global = int(n_nodes_global / 2)  # + np.random.randint(-innervation_zone_width_n_nodes/2,innervation_zone_width_n_nodes/2+1)
  
  nodes_to_stimulate_global = [innervation_node_global]
  if innervation_node_global > 0:
    nodes_to_stimulate_global.insert(0, innervation_node_global-1)
  if innervation_node_global < n_nodes_global-1:
    nodes_to_stimulate_global.append(innervation_node_global+1)

  # stimulation value
  if is_fiber_gets_stimulated:
    stimulation_current = 40.
    
    # output
    if rank_no == 0:
      print("t: {}, stimulate fiber {} at nodes {}".format(current_time, fiber_no, nodes_to_stimulate_global))
  else:
    stimulation_current = 0.

  for node_no_global in nodes_to_stimulate_global:
    parameters[(node_no_global,0,0)] = stimulation_current   # key: ((x,y,z),nodal_dof_index,parameter_no)

# callback function from output writer
def callback(data, shape, nEntries, dim, timeStepNo, currentTime, null):
  pass
    
config = {
  "scenarioName":                 scenario_name,
  "logFormat":                    "csv",     # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":   "solver_structure.txt",     # filename of file that will contain a visualization of the solver structure and data mapping
  "mappingsBetweenMeshesLogFile": "mappings_between_meshes_log.txt",    # log file for mappings 
  
  "Meshes": {
    "MeshFiber": {
      "nElements":          n_elements,
      "physicalExtent":     n_elements*element_size,      # 100 elements per cm
      "physicalOffset":     [0,0,0],
      "logKey":             "Fiber",
      "inputMeshIsGlobal":  True,
    },
  },
  "Solvers": {
    "implicitSolver": {
      "maxIterations":      1e4,
      "relativeTolerance":  1e-5,
      "absoluteTolerance":  1e-10,         # 1e-10 absolute tolerance of the residual          
      "solverType":         solver_type,
      "preconditionerType": "none",
      "dumpFormat":         "default",
      "dumpFilename":       "",   # dump of rhs and system matrix disabled (no filename specified)
    }
  },
  
  "StrangSplitting": {
    "timeStepWidth":              dt_splitting,  # 1e-1
    "endTime":                    end_time,
    "logTimeStepWidthAsKey":      "dt_splitting",
    "durationLogKey":             "duration_total",
    "timeStepOutputInterval":     1000,
    "connectedSlotsTerm1To2":     {0:0, 1:1, 2:2},   # transfer slot 0 = state Vm from Term1 (CellML) to Term2 (Diffusion), slot 1 = stress for output writer in diffusion
    "connectedSlotsTerm2To1":     {0:0, 1:1, 2:2},   # transfer the same back, in order to reuse field variables
    
    "Term1": {      # CellML
      "Heun" : {
        "timeStepWidth":                dt_0D,  # 5e-5
        "initialValues":                [],
        "timeStepOutputInterval":       1e4,
        "logTimeStepWidthAsKey":        "dt_0D",
        "durationLogKey":               "duration_0D",
        "inputMeshIsGlobal":            True,
        "dirichletBoundaryConditions":  {},
        "nAdditionalFieldVariables":    0,
        "checkForNanInf":               True,                                             # check if the solution vector contains nan or +/-inf values, if yes, an error is printed. This is a time-consuming check.                                
        
        "CellML" : {
          "modelFilename":                          cellml_file,                                    # input C++ source file or cellml XML file
          "statesInitialValues":                    [],                                             # if given, the initial values for the the states of one instance
          "initializeStatesToEquilibrium":          False,                                          # if the equilibrium values of the states should be computed before the simulation starts
          "initializeStatesToEquilibriumTimestepWidth": 1e-4,                                       # if initializeStatesToEquilibrium is enable, the timestep width to use to solve the equilibrium equation
          
          # optimization parameters
          "optimizationType":                       "vc",                                           # "vc", "simd", "openmp" type of generated optimizated source file
          "approximateExponentialFunction":         True,                                           # if optimizationType is "vc", whether the exponential function exp(x) should be approximate by (1+x/n)^n with n=1024
          "compilerFlags":                          "-fPIC -O3 -march=native -shared ",             # compiler flags used to compile the optimized model code
          "maximumNumberOfThreads":                 0,                                              # if optimizationType is "openmp", the maximum number of threads to use. Default value 0 means no restriction.
          
          # stimulation callbacks
          #"libraryFilename":                       "cellml_simd_lib.so",                           # compiled library
          
          # stimulation by setting I_Stim
          "setSpecificParametersFunction":          None, #set_specific_parameters,                        # callback function that sets parameters like stimulation current
          "setSpecificParametersCallInterval":      int(1./(5*stimulation_frequency)/dt_0D),        # set_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
          
          # stimulation by setting Vm to a prescribed value
          "setSpecificStatesFunction":              set_specific_states,                            # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
          "setSpecificStatesCallInterval":          0,                                              # int(1./stimulation_frequency/dt_0D),     # set_specific_states should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
          "setSpecificStatesCallFrequency":         stimulation_frequency,                          # set_specific_states should be called stimulation_frequency times per ms, the factor 2 is needed because every Heun step includes two calls to rhs
          "setSpecificStatesFrequencyJitter":       0,                                              # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
          "setSpecificStatesRepeatAfterFirstCall":  0.01,                                           # simulation time span for which the setSpecificStates callback will be called after a call was triggered
          "setSpecificStatesCallEnableBegin":       0,                                              # [ms] first time when to call setSpecificStates
          
          "additionalArgument":                     0,                                              # last argument that will be passed to the callback functions set_specific_states, set_specific_parameters, etc.
          #"handleResultFunction":                   handle_result,                                  # callback function that prints some current values
          #"handleResultCallInterval":               int(1./dt_0D),                                  # call interval for handle_result callback
          "handleResultFunctionAdditionalParameter": None,                                          # last parameter for handle_result callback
          
          # parameters to the cellml model
          "parametersInitialValues":                parameters_initial_values,                      #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
          "mappings":                               mappings,                                       # mappings between parameters and algebraics/constants and between outputConnectorSlots and states, algebraics or parameters, they are defined in helper.py
          
          "meshName":                               "MeshFiber",
          "stimulationLogFilename":                 "out/stimulation.log",                          # a file that will contain the times of stimulations
          
          # output writer for states, algebraics and parameters
          "OutputWriter" : [
#            {"format": "Paraview", "outputInterval": int(1./dt_1D*output_timestep), "filename": "out/cellml", "binary": True, "onlyNodalValues": True, "fixedFormat": True, "combineFiles": True, "fileNumbering": "incremental"},
            {"format": "PythonFile", "outputInterval": int(1./dt_1D*output_timestep), "filename": "out/cellml", "binary": True, "onlyNodalValues": True, "fixedFormat": True, "combineFiles": True, "fileNumbering": "incremental"},
          ],
        },
        
        # output writer only for states
        "OutputWriter" : [
          #{"format": "PythonFile", "outputInterval": int(1./dt_0D*output_timestep), "filename": "out/states", "binary": True, "onlyNodalValues": True, "fileNumbering": "incremental"},
        ],
      },
    },
    "Term2": {     # Diffusion
      "CrankNicolson" : {
        "initialValues": [],
        #"numberTimeSteps": 1,
        "timeStepWidth":                dt_1D,
        "timeStepOutputInterval":       1e4,
        "logTimeStepWidthAsKey":        "dt_1D",
        "durationLogKey":               "duration_1D",
        "inputMeshIsGlobal":            True,
        "dirichletBoundaryConditions":  {},
        "solverName":                   "implicitSolver",
        "checkForNanInf":               True,             # check if the solution vector contains nan or +/-inf values, if yes, an error is printed. This is a time-consuming check.
        "nAdditionalFieldVariables":    1,
        
        "FiniteElementMethod" : {
          "meshName":               "MeshFiber",
          "prefactor":              Conductivity/(Am*Cm),
          "solverName":             "implicitSolver",
          "inputMeshIsGlobal":      True,
        },
        
        # output writer only for the diffusion variable (i.e. state "Vm")
        "OutputWriter" : [
          {"format": "PythonFile", "outputInterval": int(1./dt_1D*output_timestep), "filename": "out/vm", "binary": True, "onlyNodalValues": False, "fileNumbering": "incremental"},
          {"format": "Paraview",   "outputInterval": int(1./dt_1D*output_timestep), "filename": "out/vm", "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"},
          #{"format": "ExFile", "filename": "out/fiber", "outputInterval": 1e5, "sphereSize": "0.02*0.02*0.02"},
        ],
      },
    },
  }
}
