# Electrophysiology
# Monodomain with either Shorten or Hodgkin-Huxley model as rhs

end_time = 30.0   # [ms] end time of simulation
n_elements = 50

# global parameters
PMax = 7.3              # maximum stress [N/cm^2]
Conductivity = 3.828    # sigma, conductivity [mS/cm]
Am = 500.0              # surface area to volume ratio [cm^-1]
Cm = 0.58               # membrane capacitance [uF/cm^2]
innervation_zone_width = 1.  # cm
innervation_zone_width = 0.  # cm
solver_type = "gmres"

cellml_file = "../input/shorten_ocallaghan_davidson_soboleva_2007.c"
cellml_file = "../input/shorten.cpp"

# timing parameters
stimulation_frequency = 10.0      # stimulations per ms
dt_1D = 1e-3                      # timestep width of diffusion
dt_0D = 3e-3                      # timestep width of ODEs
dt_3D = 3e-3                      # overall timestep width of splitting
output_timestep = 1e0             # timestep for output files

# import needed packages
import numpy as np
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt

print("prefactor: ",Conductivity/(Am*Cm))

# determine if fibre gets stimulation at given time
def fibre_gets_stimulated(current_time):
  a = current_time
  
  if a - int(a) < 0.1 and a < 5:
    return True
  else:
    return False
 
def set_parameters(n_nodes_global, time_step_no, current_time, parameters, dof_nos_global, fibre_no):
  
  # determine if fibre gets stimulated at the current time
  is_fibre_gets_stimulated = fibre_gets_stimulated(current_time)
  
  # determine nodes to stimulate (center node, left and right neighbour)
  innervation_zone_width_n_nodes = innervation_zone_width*100  # 100 nodes per cm
  innervation_node_global = int(n_nodes_global / 2)  # + np.random.randint(-innervation_zone_width_n_nodes/2,innervation_zone_width_n_nodes/2+1)
  nodes_to_stimulate_global = [innervation_node_global]
  
  for k in range(10):
    if innervation_node_global-k >= 0:
      nodes_to_stimulate_global.insert(0, innervation_node_global-k)
    if innervation_node_global+k <= n_nodes_global-1:
      nodes_to_stimulate_global.append(innervation_node_global+k)
  
  # stimulation value
  if fibre_gets_stimulated:
    stimulation_current = 40.
  else:
    stimulation_current = 0.
  
  first_dof_global = dof_nos_global[0]
  last_dof_global = dof_nos_global[-1]
    
  for node_no_global in nodes_to_stimulate_global:
    if first_dof_global <= node_no_global <= last_dof_global:
      # get local no for global no (1D)
      dof_no_local = node_no_global - first_dof_global
      parameters[dof_no_local] = stimulation_current
 
  
def set_specific_parameters(n_nodes_global, time_step_no, current_time, parameters, fibre_no):
  
  # determine if fibre gets stimulated at the current time
  is_fibre_gets_stimulated = fibre_gets_stimulated(current_time)
  
  # determine nodes to stimulate (center node, left and right neighbour)
  innervation_zone_width_n_nodes = innervation_zone_width*100  # 100 nodes per cm
  innervation_node_global = int(n_nodes_global / 2)  # + np.random.randint(-innervation_zone_width_n_nodes/2,innervation_zone_width_n_nodes/2+1)
  nodes_to_stimulate_global = [innervation_node_global]
  
  for k in range(10):
    if innervation_node_global-k >= 0:
      nodes_to_stimulate_global.insert(0, innervation_node_global-k)
    if innervation_node_global+k <= n_nodes_global-1:
      nodes_to_stimulate_global.append(innervation_node_global+k)
 
  nodes_to_stimulate_global += [2,5,10]

  # stimulation value
  if fibre_gets_stimulated:
    stimulation_current = 40.
  else:
    stimulation_current = 0.
  
  for node_no_global in nodes_to_stimulate_global:
    parameters[(node_no_global,0)] = stimulation_current
    
fig = plt.figure(1)
#plt.ion()

# callback function that is called after integration of rhs, generates plots
def handleResult(n_instances, time_step_no, current_time, states, intermediates, null):
  
  # collect data for every instance
  xdata = []
  vm_data = []
  gamma_data = []
  for i in range(nInstances):
    xdata.append(i)
    vm_data.append(states[i])
    gamma_data.append(intermediates[i])
  
  # generate plot of Vm and gamma
  # prepare figure
  plt.figure(1)
  plt.clf()
  plt.xlabel('position $x$')
  ax1 = plt.gca()
  ax1.plot(xdata, vm_data, 'go-', label='$V_m$')
  plt.ylim(-80, 80)
  plt.xlabel('t')
  plt.ylabel('$V_m$')
  ax2 = ax1.twinx()
  
  # plot data
  ax2.plot(xdata, gamma_data, 'ro-', label='$\gamma$')
  plt.ylabel('$\gamma$')    
  plt.ylim(0, 1)
  
  # ask matplotlib for the plotted objects and their labels
  lines, labels = ax1.get_legend_handles_labels()
  lines2, labels2 = ax2.get_legend_handles_labels()
  ax2.legend(lines + lines2, labels + labels2, loc=0)
  
  # save to png file
  filename = "out_{:06.1f}.png".format(currentTime)
  plt.savefig(filename)
  #print "   saved ""{}""".format(filename)
  #plt.draw()
    
# callback function from output writer
def callback(data, shape, nEntries, dim, timeStepNo, currentTime, null):
  pass
    
bc = {0: -75, -1: -75}
config = {
  "Meshes": {
    "MeshFibre": {
      "nElements": n_elements,
      "physicalExtent": n_elements/100.,
      "inputMeshIsGlobal": True,
      "setHermiteDerivatives": False
    },
  },
  "Solvers": {
    "implicitSolver": {
      "maxIterations": 1e4,
      "relativeTolerance": 1e-10,
      "solverType": solver_type,
      "preconditionerType": "none",
    }
  },
  "StrangSplitting": {
    #"numberTimeSteps": 1,
    "timeStepWidth": dt_3D,  # 1e-1
    "endTime": end_time,
    "logTimeStepWidthAsKey": "dt_3D",
    "durationLogKey": "duration_total",
    "timeStepOutputInterval" : 200,
    
    "Term1": {      # CellML
      "Heun" : {
        "timeStepWidth": dt_0D,  # 5e-5
        "initialValues": [],
        "timeStepOutputInterval": 1e4,
        "inputMeshIsGlobal": True,
        "dirichletBoundaryConditions": {},
  
        "CellML" : {
          "sourceFilename": cellml_file,             # input C++ source file, can be either generated by OpenCMISS or OpenCOR from cellml model
          #"simdSourceFilename" : "simdcode.cpp",     # transformed C++ source file that gets generated from sourceFilename and is ready for multiple instances
          #"libraryFilename": "cellml_simd_lib.so",   # compiled library
          "useGivenLibrary": False,
          #"statesInitialValues": [],
#          "setParametersFunction": set_parameters,    # callback function that sets parameters like stimulation current
#          "setParametersCallInterval": int(3./stimulation_frequency/dt_0D),     # set_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
          "setSpecificParametersFunction": set_specific_parameters,    # callback function that sets parameters like stimulation current
          "setSpecificParametersCallInterval": int(3./stimulation_frequency/dt_0D),     # set_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
           
          #"handleResultFunction": handleResult,
          #"handleResultCallInterval": 2e3,
          
          "outputStateIndex": 0,     # state 0 = Vm, rate 28 = gamma
          "parametersUsedAsIntermediate": [32],       # list of intermediate value indices, that will be set by parameters. Explicitely defined parameters that will be copied to intermediates, this vector contains the indices of the algebraic array. This is ignored if the input is generated from OpenCMISS generated c code.
          "parametersUsedAsConstant": [65],           # list of constant value indices, that will be set by parameters. This is ignored if the input is generated from OpenCMISS generated c code.
          "parametersInitialValues": [0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
          "meshName": "MeshFibre",
          "prefactor": 1.0,
        },
        
        "OutputWriter" : [
          #{"format": "PythonFile", "outputInterval": 1e4, "filename": "out/states", "binary": True, "onlyNodalValues":True},
        ],
      },
    },
    "Term2": {     # Diffusion
      "ImplicitEuler" : {
        "initialValues": [],
        #"numberTimeSteps": 1,
        "timeStepWidth": dt_1D,
        "logTimeStepWidthAsKey": "dt_1D",
        "durationLogKey": "duration_1D",
        "timeStepOutputInterval": 1e4,
        "dirichletBoundaryConditions": bc,
        "inputMeshIsGlobal": True,
        "solverName": "implicitSolver",
        "FiniteElementMethod" : {
          "maxIterations": 1e4,
          "relativeTolerance": 1e-10,
          "meshName": "MeshFibre",
          "prefactor": Conductivity/(Am*Cm),
          "solverName": "implicitSolver",
          "inputMeshIsGlobal": True,
        },
        "OutputWriter" : [
          {"format": "Paraview", "outputInterval": int(1./dt_1D*output_timestep), "filename": "out/fibre", "binary": True, "fixedFormat": False, "combineFiles": True, "onlyNodalValues":True},
          #{"format": "Paraview", "outputInterval": 1./dt_1D*output_timestep, "filename": "out/fibre_"+str(i)+"_txt", "binary": False, "fixedFormat": False},
          #{"format": "ExFile", "filename": "out/fibre_"+str(i), "outputInterval": 1./dt_1D*output_timestep, "sphereSize": "0.02*0.02*0.02"},
          #{"format": "PythonFile", "filename": "out/fibre", "outputInterval": int(1./dt_1D*output_timestep), "binary":True, "onlyNodalValues":True},
        ]
      },
    },
  }
}
