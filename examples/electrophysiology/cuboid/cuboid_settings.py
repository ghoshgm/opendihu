# Cuboid
#
# arguments: <n_processes_per_fiber> <n_fibers> <n_nodes_per_fiber> <scenario_name>

end_time = 0.2

import numpy as np
import pickle
import sys

# global parameters
PMax = 7.3              # maximum stress [N/cm^2]
Conductivity = 3.828    # sigma, conductivity [mS/cm]
Am = 500.0              # surface area to volume ratio [cm^-1]
Cm = 0.58               # membrane capacitance [uF/cm^2]
innervation_zone_width = 1.  # cm
innervation_zone_width = 0.  # cm
  
cellml_file = "../../input/shorten_ocallaghan_davidson_soboleva_2007.c"
cellml_file = "../../input/shorten.cpp"
#cellml_file = "../../input/hodgkin_huxley_1952.c"

fiber_distribution_file = "../../input/MU_fibre_distribution_3780.txt"
firing_times_file = "../../input/MU_firing_times_real.txt"

# timing parameters
stimulation_frequency = 10.0      # stimulations per ms
dt_1D = 1e-3                      # timestep width of diffusion
dt_0D = 3e-3                      # timestep width of ODEs
dt_3D = 3e-3                      # overall timestep width of splitting
output_timestep = 1e0             # timestep for output files
n_nodes_per_fiber = 1000             # number of nodes per fiber
n_fibers = 10

#print("prefactor: ",Conductivity/(Am*Cm))
#print("numpy path: ",np.__path__)

# parse arguments
try:
  n_processes_per_fiber = (int)(sys.argv[0])
  n_fibers = (int)(sys.argv[1])
  n_nodes_per_fiber = (int)(sys.argv[2])
  scenario_name = sys.argv[3]
except:
  print("arguments: <n_processes_per_fiber> <n_fibers> <n_nodes_per_fiber> <scenario_name>")
  sys.exit(0)


solver_type = "cg"
if scenario_name == "Strong_scaling":
  end_time = 1.0
elif scenario_name == "Strong_scalig_LU":
  solver_type = "lu"
  end_time = 1.0
elif scenario_name == "solver_scaling":
  solver_type = sys.argv[4]
  end_time = 1.0

rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

# adjust n_processes_per_fiber if there are too many ranks
surplus = n_processes_per_fiber*n_fibers - n_ranks
while surplus > n_fibers:
  n_processes_per_fiber += 1
  surplus = n_processes_per_fiber*n_fibers - n_ranks


if rank_no == 0:
  print("scenario_name: {}".format(scenario_name))
  print("n_processes_per_fiber: {}, n_fibers: {}, n_nodes_per_fiber: {}".format(n_processes_per_fiber, n_fibers, n_nodes_per_fiber))
  print("solver_type: {}".format(solver_type))

#print("rank: {}/{}".format(rank_no,n_ranks))
   
# load MU distribution and firing times
fiber_distribution = np.genfromtxt(fiber_distribution_file, delimiter=" ")
firing_times = np.genfromtxt(firing_times_file)

# set values for cellml model
if "shorten" in cellml_file:
  parameters_used_as_intermediate = [32]
  parameters_used_as_constant = [65]
  parameters_initial_values = [0.0, 1.0]
  nodal_stimulation_current = 400.
  
elif "hodgkin_huxley" in cellml_file:
  parameters_used_as_intermediate = []
  parameters_used_as_constant = [2]
  parameters_initial_values = [0.0]
  nodal_stimulation_current = 40.

def get_motor_unit_no(fiber_no):
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
    #print("{}: fiber {} is mu {}, t = {}, row: {}, stimulated: {} {}".format(rank_no, fiber_no, mu_no, current_time, (index % n_firing_times), firing_times[index % n_firing_times, mu_no], "true" if firing_times[index % n_firing_times, mu_no] == 1 else "false"))
  
  return firing_times[index % n_firing_times, mu_no] == 1
  
# callback function that can set states, i.e. prescribed values for stimulation
def set_specific_states(n_nodes_global, time_step_no, current_time, states, fiber_no):
  
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

def set_parameters_null(n_nodes_global, time_step_no, current_time, parameters, dof_nos_global, fibre_no):
  pass
  
def set_parameters(n_nodes_global, time_step_no, current_time, parameters, dof_nos_global, fibre_no):
  
  # determine if fibre gets stimulated at the current time
  fibre_gets_stimulated = fibreGetsStimulated(fibre_no, stimulation_frequency, current_time)
  
  # determine nodes to stimulate (center node, left and right neighbour)
  innervation_zone_width_n_nodes = innervation_zone_width*100  # 100 nodes per cm
  innervation_node_global = int(n_nodes_global / 2)  # + np.random.randint(-innervation_zone_width_n_nodes/2,innervation_zone_width_n_nodes/2+1)
  nodes_to_stimulate_global = [innervation_node_global]
  if innervation_node_global > 0:
    nodes_to_stimulate_global.insert(0, innervation_node_global-1)
  if innervation_node_global < n_nodes_global-1:
    nodes_to_stimulate_global.append(innervation_node_global+1)
  
  # stimulation value
  if fibre_gets_stimulated:
    stimulation_current = 400.
  else:
    stimulation_current = 0.
  
  first_dof_global = dof_nos_global[0]
  last_dof_global = dof_nos_global[-1]
    
  for node_no_global in nodes_to_stimulate_global:
    if first_dof_global <= node_no_global <= last_dof_global:
      # get local no for global no (1D)
      dof_no_local = node_no_global - first_dof_global
      parameters[dof_no_local] = stimulation_current
 
      #print("       {}: set stimulation for local dof {}".format(rank_no, dof_no_local))
  
  #print("       {}: setParameters at timestep {}, t={}, n_nodes_global={}, range: [{},{}], fibre no {}, MU {}, stimulated: {}".\
        #format(rank_no, time_step_no, current_time, n_nodes_global, first_dof_global, last_dof_global, fibre_no, getMotorUnitNo(fibre_no), fibre_gets_stimulated))
    
  #wait = input("Press any key to continue...")
    
# callback function that can set parameters, i.e. stimulation current
def set_specific_parameters(n_nodes_global, time_step_no, current_time, parameters, fiber_no):
  
  # determine if fiber gets stimulated at the current time
  is_fiber_gets_stimulated = fiber_gets_stimulated(fiber_no, stimulation_frequency, current_time)
  
  # determine nodes to stimulate (center node, left and right neighbour)
  innervation_zone_width_n_nodes = innervation_zone_width*100  # 100 nodes per cm
  innervation_node_global = int(n_nodes_global / 2)  # + np.random.randint(-innervation_zone_width_n_nodes/2,innervation_zone_width_n_nodes/2+1)
  nodes_to_stimulate_global = [innervation_node_global]
  
  for k in range(0):
    if innervation_node_global-k >= 0:
      nodes_to_stimulate_global.insert(0, innervation_node_global-k)
    if innervation_node_global+k <= n_nodes_global-1:
      nodes_to_stimulate_global.append(innervation_node_global+k)
  
  # stimulation value
  if is_fiber_gets_stimulated:
    stimulation_current = 10*nodal_stimulation_current
    if rank_no == 0:
      print("t: {}, stimulate fiber {} at nodes {}".format(current_time, fiber_no, nodes_to_stimulate_global))

  else:
    stimulation_current = 0.

  for node_no_global in nodes_to_stimulate_global:
    parameters[(node_no_global,0,0)] = stimulation_current   # key: ((x,y,z),nodal_dof_index)

def callback(data, shape, nEntries, dim, timeStepNo, currentTime):
  pass
    
def get_instance_config(i):

  # set ranks list containing the rank nos for fiber i 
  if n_processes_per_fiber > 0:
    ranks = []
    n_normal_fibers = n_fibers - surplus
    if i < n_normal_fibers:
      n_previous_ranks = n_processes_per_fiber*i
    else:
      n_previous_ranks = n_processes_per_fiber*n_normal_fibers + (n_processes_per_fiber+1)*(i-n_normal_fibers)

    n_processes_on_this_fiber = n_processes_per_fiber
    if i >= n_normal_fibers:
      n_processes_on_this_fiber = n_processes_per_fiber+1

    for j in range(n_processes_on_this_fiber):
      rank_id = n_previous_ranks + j
      ranks.append(rank_id)
  else:
    ranks = [int(i/-n_processes_per_fiber)]

  bc = {0: -82.747, -1: -82.747}
  bc = {}
  instance_config = {
    "ranks": [0],
    "StrangSplitting": {
      #"numberTimeSteps": 1,
      "timeStepWidth": dt_3D,  # 1e-1
      "logTimeStepWidthAsKey": "dt_3D",
      "durationLogKey": "duration_total",
      "timeStepOutputInterval" : 100,
      "endTime": end_time,
      "outputData1": False,
      "outputData2": True,

      "Term1": {      # CellML
        "Heun" : {
          "timeStepWidth": dt_0D,  # 5e-5
          "logTimeStepWidthAsKey": "dt_0D",
          "durationLogKey": "duration_0D",
          "initialValues": [],
          "timeStepOutputInterval": 1e4,
          "inputMeshIsGlobal": True,
          "dirichletBoundaryConditions": {},
          
          "CellML" : {
            "sourceFilename": cellml_file,             # input C++ source file, can be either generated by OpenCMISS or OpenCOR from cellml model
            "compilerFlags": "-fPIC -O3 -shared ",
            #"simdSourceFilename" : "simdcode.cpp",     # transformed C++ source file that gets generated from sourceFilename and is ready for multiple instances
            #"libraryFilename": "cellml_simd_lib.so",   # compiled library
            "useGivenLibrary": False,
            #"statesInitialValues": [],
            "setSpecificParametersFunction": set_specific_parameters,    # callback function that sets parameters like stimulation current
            "setSpecificParametersCallInterval": 2*int(1./stimulation_frequency/dt_0D),     # set_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
            #"setSpecificStatesFunction": set_specific_states,    # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
            #"setSpecificStatesCallInterval": 2*int(1./stimulation_frequency/dt_0D),     # set_specific_states should be called stimulation_frequency times per ms, the factor 2 is needed because every Heun step includes two calls to rhs
            "additionalArgument": i,
            
            "outputStateIndex": 0,     # state 0 = Vm, rate 28 = gamma
            "parametersUsedAsIntermediate": parameters_used_as_intermediate,  #[32],       # list of intermediate value indices, that will be set by parameters. Explicitely defined parameters that will be copied to intermediates, this vector contains the indices of the algebraic array. This is ignored if the input is generated from OpenCMISS generated c code.
            "parametersUsedAsConstant": parameters_used_as_constant,          #[65],           # list of constant value indices, that will be set by parameters. This is ignored if the input is generated from OpenCMISS generated c code.
            "parametersInitialValues": parameters_initial_values,            #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
            "meshName": "MeshFiber_{}".format(i),
            "prefactor": 1.0,
          },
        },
      },
      "Term2": {     # Diffusion
        "CrankNicolson" : {
          "initialValues": [],
          #"numberTimeSteps": 1,
          "timeStepWidth": dt_1D,  # 1e-5
          "logTimeStepWidthAsKey": "dt_1D",
          "durationLogKey": "duration_1D",
          "timeStepOutputInterval": 1e4,
          "dirichletBoundaryConditions": bc,
          "inputMeshIsGlobal": True,
          "solverName": "implicitSolver",
          "FiniteElementMethod" : {
            "solverName": "implicitSolver",
            "inputMeshIsGlobal": True,
            "meshName": "MeshFiber_{}".format(i),
            "prefactor": Conductivity/(Am*Cm),
          },
          "OutputWriter" : [
            {"format": "Paraview", "outputInterval": (int)(1./dt_1D*output_timestep), "filename": "out/fibre_"+str(i), "binary": True, "fixedFormat": False, "combineFiles":False},
            #{"format": "Paraview", "outputInterval": 1./dt_1D*output_timestep, "filename": "out/fibre_"+str(i)+"_txt", "binary": False, "fixedFormat": False},
            #{"format": "ExFile", "filename": "out/fibre_"+str(i), "outputInterval": 1./dt_1D*output_timestep, "sphereSize": "0.02*0.02*0.02"},
            #{"format": "PythonFile", "filename": "out/fibre_"+str(i), "outputInterval": int(1./dt_1D*output_timestep), "binary":True, "onlyNodalValues":True, "combineFiles":False},
          ]
        },
      },
    }
  }
  return instance_config
    
config = {
  "scenarioName": scenario_name,
  "Meshes": {
    "MeshFiber_{}".format(i): {
      "nElements": n_nodes_per_fiber-1,
      "nodePositions": [[x,i,0] for x in np.linspace(0,(n_nodes_per_fiber-1)/100.,n_nodes_per_fiber)],
      "inputMeshIsGlobal": True,
      "setHermiteDerivatives": False,
      "logKey": "1D"
    }
    for i in range(n_fibers)
  },
  "Solvers": {
    "implicitSolver": {
      "maxIterations": 1e4,
      "relativeTolerance": 1e-10,
      "solverType": solver_type,
      "preconditionerType": "none",
    }
  },
  "MultipleInstances": {
    "nInstances": n_fibers,
    "instances": [get_instance_config(i) for i in range(n_fibers)],
  }
}
