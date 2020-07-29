# Electrophysiology debug
#
#  parameters: <fibre no> [<cellml_file> [<end_time>]]

end_time = 10.0

import numpy as np
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle
import sys

# global parameters
PMax = 7.3              # maximum stress [N/cm^2]
Conductivity = 3.828    # sigma, conductivity [mS/cm]
Am = 500.0              # surface area to volume ratio [cm^-1]
Cm = 0.58               # membrane capacitance [uF/cm^2]
innervation_zone_width = 1.  # cm
innervation_zone_width = 0.  # cm
  
cellml_file = "../input/shorten_ocallaghan_davidson_soboleva_2007.c"
cellml_file = "../input/shorten_opencmiss.cpp"

# timing parameters
stimulation_frequency = 10.0      # stimulations per ms
dt_1D = 1e-5                      # timestep width of diffusion
dt_0D = 5e-5                      # timestep width of ODEs
dt_3D = 1e-1                      # overall timestep width of splitting
output_timestep = 1e-1            # timestep for output files

# input files
#fibre_file = "../input/laplace3d_structured_quadratic"
fibre_file = "../input/laplace3d_structured_linear"
fibre_distribution_file = "../input/MU_fibre_distribution_3780.txt"
#firing_times_file = "../input/MU_firing_times_real.txt"
firing_times_file = "../input/MU_firing_times_immediately.txt"

#print("prefactor: ",Conductivity/(Am*Cm))
#print("numpy path: ",np.__path__)

# parse fibre number as first command line parameter
if len(sys.argv) > 0:
  if "check_results.py" not in sys.argv[0]:
    try:
      fibre_no = int(sys.argv[0])
      print("fibre_no: {}".format(fibre_no))
    except:
      print("could not parse fibre_no");
      
    if len(sys.argv) > 1:
      cellml_file = sys.argv[1]
    if len(sys.argv) > 2:
      end_time = float(sys.argv[2])
    
    print("end_time: \"{}\", cellml_file: {}".format(end_time, cellml_file))

subset_length = None
if len(sys.argv) > 1:
  try:
    subset_length = int(sys.argv[1])
    print("subset_length: {}".format(subset_length))
  except:
    print("could not parse subset_length");

if "shorten" in cellml_file:
  parametersUsedAsIntermediate = [32]
  parametersUsedAsConstant = [65]
  parametersInitialValues = [0.0, 1.0]
  
elif "hodgkin_huxley" in cellml_file:
  parametersUsedAsIntermediate = []
  parametersUsedAsConstant = [2]
  parametersInitialValues = [0.0]
  

def getMotorUnitNo(fibre_no):
  return int(fibre_distribution[fibre_no % len(fibre_distribution)]-1)

def fibreGetsStimulated(fibre_no, frequency, current_time):

  # determine motor unit
  mu_no = getMotorUnitNo(fibre_no)
  
  # determine if fibre fires now
  index = int(current_time * frequency)
  n_firing_times = np.size(firing_times,0)
  return firing_times[index % n_firing_times, mu_no] == 1
  
def set_parameters(n_nodes, time_step_no, current_time, parameters, fibre_no):
  
  # determine if fibre gets stimulated at the current time
  fibre_gets_stimulated = fibreGetsStimulated(fibre_no, stimulation_frequency, current_time)
  
  # determine nodes to stimulate (center node, left and right neighbour)
  innervation_zone_width_n_nodes = innervation_zone_width*100  # 100 nodes per cm
  innervation_node = int(n_nodes / 2) + np.random.randint(-innervation_zone_width_n_nodes/2,innervation_zone_width_n_nodes/2+1)
  nodes_to_stimulate = [innervation_node]
  if innervation_node > 0:
    nodes_to_stimulate.insert(0, innervation_node-1)
  if innervation_node < n_nodes-1:
    nodes_to_stimulate.append(innervation_node+1)
  
  # stimulation value
  if fibre_gets_stimulated:
    stimulation_current = 400.
  else:
    stimulation_current = 0.
  
  for node_no in nodes_to_stimulate:
    parameters[node_no] = stimulation_current

  print("       setParameters at timestep {}, t={}, n_nodes={}, fibre no {}, MU {}, stimulated: {}".\
    format(time_step_no, current_time, n_nodes, fibre_no, getMotorUnitNo(fibre_no), fibre_gets_stimulated))
    
 
  print("       set stimulation for nodes {}".format(nodes_to_stimulate))
  
  #wait = input("Press any key to continue...")
    
fig = plt.figure(1)
#plt.ion()

def handleResult(n_nodes, time_step_no, current_time, states, intermediates):
  #print "handleResult: time step {}, t={}, n_nodes: {}, n entries states: {}".format(time_step_no, current_time, n_nodes, len(states))
  #print "states:", states[0:n_nodes]
  
  xdata = []
  vm_data = []
  gamma_data = []
  for i in range(n_nodes):
    xdata.append(i)
    vm_data.append(states[i])
    gamma_data.append(intermediates[i])
  
  # write out a png file
  plt.figure(1)
  plt.clf()
  plt.xlabel('position $x$')
  ax1 = plt.gca()
  ax1.plot(xdata, vm_data, 'go-', label='$V_m$')
  plt.ylim(-80, 80)
  plt.xlabel('t')
  plt.ylabel('$V_m$')
  ax2 = ax1.twinx()
  ax2.plot(xdata, gamma_data, 'ro-', label='$\gamma$')
  plt.ylabel('$\gamma$')    
  plt.ylim(0, 1)
  
  # ask matplotlib for the plotted objects and their labels
  lines, labels = ax1.get_legend_handles_labels()
  lines2, labels2 = ax2.get_legend_handles_labels()
  ax2.legend(lines + lines2, labels + labels2, loc=0)
  
  filename = "out_{:06.1f}.png".format(currentTime)
  plt.savefig(filename)
  #print "   saved ""{}""".format(filename)
  #plt.draw()
    
def callback(data, shape, nEntries, dim, timeStepNo, currentTime):
  pass
    
def get_instance_config(i):
  
  bc = {0: -80.25667, -1: -80.25667}
  instance_config = {
    "StrangSplitting": {
      #"numberTimeSteps": 1,
      "timeStepWidth": dt_3D,  # 1e-1
      "endTime": end_time,
      "outputData1": False,
      "outputData2": True,

      "Term1": {      # CellML
        "Heun" : {
          "timeStepWidth": dt_0D,  # 5e-5
          #"initialValues": [],
          "timeStepOutputInterval": 1e4,
            
          "CellML" : {
            "sourceFilename": cellml_file,             # input C++ source file, can be either generated by OpenCMISS or OpenCOR from cellml model
            #"simdSourceFilename" : "simdcode.cpp",     # transformed C++ source file that gets generated from sourceFilename and is ready for multiple instances
            #"libraryFilename": "cellml_simd_lib.so",   # compiled library
            "useGivenLibrary": False,
            #"statesInitialValues": [],
            "setParametersFunction": set_parameters,    # callback function that sets parameters like stimulation current
            "setParametersCallInterval": 1./stimulation_frequency/dt_0D,     # set_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
            "setParametersFunctionAdditionalParameter": i,
            
            "outputStateIndex": 0,     # state 0 = Vm, rate 28 = gamma
            "parametersUsedAsIntermediate": parametersUsedAsIntermediate,  #[32],       # list of intermediate value indices, that will be set by parameters. Explicitely defined parameters that will be copied to intermediates, this vector contains the indices of the algebraic array. This is ignored if the input is generated from OpenCMISS generated c code.
            "parametersUsedAsConstant": parametersUsedAsConstant,          #[65],           # list of constant value indices, that will be set by parameters. This is ignored if the input is generated from OpenCMISS generated c code.
            "parametersInitialValues": parametersInitialValues,            #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
            "meshName": "MeshFibre"+str(i),
            "prefactor": 1.0,
          },
        },
      },
      "Term2": {     # Diffusion
        "Heun" : {
          #"initialValues": [2,2,4,5,2,2],
          #"numberTimeSteps": 1,
          "timeStepWidth": dt_1D,  # 1e-5
          "timeStepOutputInterval": 1e4,
          "FiniteElementMethod" : {
            "relativeTolerance": 1e-10,
            "meshName": "MeshFibre"+str(i),
            "prefactor": Conductivity/(Am*Cm),
            "dirichletBoundaryConditions": bc,
          },
          "OutputWriter" : [
            {"format": "Paraview", "outputInterval": 1./dt_1D*output_timestep, "filename": "out/fibre_"+str(i)+"_bin", "binary": True, "fixedFormat": False},
            {"format": "Paraview", "outputInterval": 1./dt_1D*output_timestep, "filename": "out/fibre_"+str(i)+"_txt", "binary": False, "fixedFormat": False},
            #{"format": "ExFile", "filename": "out/fibre_"+str(i), "outputInterval": 1./dt_1D*output_timestep, "sphereSize": "0.02*0.02*0.02"},
            {"format": "PythonFile", "filename": "out/fibre_"+str(i), "outputInterval": 1./dt_1D*output_timestep, "binary":True, "onlyNodalValues":True},
          ]
        },
      },
    }
  }
  return instance_config
    
    
# create fibre meshes
meshes = {}

with open(fibre_file, "rb") as f:
  streamlines = pickle.load(f)
    
nInstances = len(streamlines)
print("nInstances: {}".format(nInstances))

#nInstances = 1
    
for i,streamline in enumerate(streamlines):
  
  # extract only part of streamline according to given subset_length
  if subset_length is not None:
    streamline = streamline[1:min(subset_length,len(streamline))]
  
    if i == fibre_no:
      print("element lengths: ")
      distances = []
      for j in range(len(streamline)-1):
        distances.append(np.linalg.norm(np.array(streamline[j]) - np.array(streamline[j+1])))
      print(distances)
  
  # define mesh
  meshes["MeshFibre{}".format(i)] = {
    "nElements": len(streamline)-1,
    "nodePositions": streamline
  }
    
# load MU distribution and firing times
fibre_distribution = np.genfromtxt(fibre_distribution_file, delimiter=" ")
firing_times = np.genfromtxt(firing_times_file)

# determine when the fibres will fire, for debugging output
print("Debugging output about fibre firing: Taking input from file \"{}\"".format(firing_times_file))
for fibre_no_index in range(nInstances):
  first_stimulation = None
  for current_time in np.linspace(0,1./stimulation_frequency*nInstances,nInstances+1):
    if fibreGetsStimulated(fibre_no_index, stimulation_frequency, current_time):
      first_stimulation = current_time
      break

  print("   Fibre {} is of MU {} and will be stimulated for the first time at {}".format(fibre_no_index, getMotorUnitNo(fibre_no_index), first_stimulation))


i = fibre_no
config = {
  "logFormat": "csv",
  "Meshes": meshes,
  "instances": [get_instance_config(i) for i in range(nInstances)],
  "OutputWriter" : [
    #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binary": "false", "fixedFormat": False},
    #{"format": "ExFile", "filename": "out/single_fibre_"+str(i), "outputInterval": 1},
    #{"format": "PythonFile", "filename": "out/single_fibre"+str(i), "binary":True, "onlyNodalValues":True},
  ]
}

# merge config dict with get_instance_config
config.update(get_instance_config(fibre_no));
    
