# This file contains all global variables for the fibers_emg example and their default values. These are the parameters and other internal variables.
# These values will be used by all involved scripts: helper.py, create_partitioned_meshes_for_settings.py and settings_fibers_emg.py
# settings_fibers_emg.py handles setting the parameter values. Those can be overridden on the command line and by specifying a custom variables.py script
# To run the simulation use the settings_fibers_emg.py file, which imports this file, e.g. ./fibers_emg ../settings_fibers_emg.py custom_variables.py

# material parameters
# --------------------
PMax = 7.3                          # maximum stress [N/cm^2]
Conductivity = 3.828                # sigma, conductivity [mS/cm]
Am = 500.0                          # surface area to volume ratio [cm^-1]
Cm = 0.58                           # membrane capacitance [uF/cm^2]
damping_factor = 0                  # velocity dependent damping factor

innervation_zone_width = 0.         # not used [cm], this will later be used to specify a variance of positions of the innervation point at the fibers

# solvers
# -------
diffusion_solver_type = "cg"        # solver and preconditioner for the diffusion part of the Monodomain equation
diffusion_preconditioner_type = "none"      # preconditioner
diffusion_solver_maxit = 1e4
diffusion_solver_reltol = 1e-10
potential_flow_solver_type = "gmres"        # solver and preconditioner for an initial Laplace flow on the domain, from which fiber directions are determined
potential_flow_preconditioner_type = "gamg" # preconditioner
potential_flow_solver_maxit = 1e4
potential_flow_solver_reltol = 1e-10
emg_solver_type = "cg"              # solver and preconditioner for the 3D static Bidomain equation that solves the intra-muscular EMG signal
emg_preconditioner_type = "none"    # preconditioner
emg_initial_guess_nonzero = False   # If the initial guess for the emg linear system should be set to the previous solution
emg_solver_maxit = 1e4
emg_solver_abstol = 1e-5
emg_solver_reltol = 1e-5

# elasticity
elasticity_solver_type = "preonly"
elasticity_preconditioner_type = "lu"
snes_max_iterations = 10                  # maximum number of iterations in the nonlinear solver
snes_rebuild_jacobian_frequency = 2       # how often the jacobian should be recomputed, -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again 
snes_relative_tolerance = 1e-5      # relative tolerance of the nonlinear solver
snes_absolute_tolerance = 1e-5      # absolute tolerance of the nonlinear solver
linear_relative_tolerance = 1e-5           # relative tolerance of the residual of the linear solver
linear_absolute_tolerance = 1e-10          # absolute tolerance of the residual of the linear solver


# timing parameters
# -----------------
end_time = 20.0                   # [ms] end time of the simulation
stimulation_frequency = 100*1e-3    # [ms^-1] sampling frequency of stimuli in firing_times_file, in stimulations per ms, number before 1e-3 factor is in Hertz. This is not used here.
dt_muscle_spindles          = 1e-3  # [ms]
dt_golgi_tendon_organs      = 1e-3  # [ms]
dt_motoneuron               = 1e-3  # [ms]
dt_0D = 0.5e-3                        # [ms] timestep width of ODEs
dt_1D = 1e-3                      # [ms] timestep width of diffusion
dt_bidomain = 1e-2                  # [ms] timestep width of multidomain
dt_splitting_0D1D = 1e-3            # [ms] overall timestep width of strang splitting
dt_elasticity = 1e0                 # [ms] time step width of elasticity solver
output_timestep = 1e0               # [ms] timestep for output files
activation_start_time = 0           # [ms] time when to start checking for stimulation
output_timestep_fibers = 2   # [ms] timestep for multidomain output files
output_timestep_elasticity = 1    # [ms] timestep for elasticity output files
output_timestep_emg = 20    # [ms] timestep for emg output files

output_timestep_golgi_tendon_organs = 20
output_timestep_neurons = 1         # [ms] timestep for output of files for all sensor organs and neurons
output_timestep_motoneuron = 1    # [ms] timestep for output of files for motoneuron
output_timestep_surface = 20



# input files
# -----------
# CellML model, Shorten or Hodgkin-Huxley
#cellml_file = "../../../input/shorten_ocallaghan_davidson_soboleva_2007.c"
#cellml_file = "../../../input/shorten.cpp"
#cellml_file = "../../../input/hodgkin_huxley_1952.c"

debug_output = False                # verbose output in this python script, for debugging the domain decomposition
disable_firing_output = True        # Disables the initial list of fiber firings on the console to save some console space
paraview_output = True             # If the paraview output writer should be enabled
adios_output = False                # If the MegaMol/ADIOS output writer should be enabled
python_output = False               # If the Python output writer should be enabled
exfile_output = False               # If the Exfile output writer should be enabled
initial_guess_nonzero = True        # if the initial guess of the multidomain solver should be set to the previous values, this is only possible if an iterative solver is used
theta = 0.5                         # weighting factor of implicit term in Crank-Nicolson scheme, 0.5 gives the classic, 2nd-order Crank-Nicolson scheme, 1.0 gives implicit euler
use_symmetric_preconditioner_matrix = True    # if the diagonal blocks of the system matrix should be used as preconditioner matrix
use_lumped_mass_matrix = False      # which formulation to use, the formulation with lumped mass matrix (True) is more stable but approximative, the other formulation (False) is exact but needs more iterations
show_linear_solver_output = True    # if convergence information of the linear solver in every timestep should be printed, this is a lot of output for fast computations
optimization_type = "vc"            # the optimization_type used in the cellml adapter, "vc" uses explicit vectorization
approximate_exponential_function = False   # if the exponential function should be approximated by a Taylor series with only 11 FLOPS
dynamic = True                      # if the dynamic hyperelasticity solver should be used

# motor unit stimulation times
firing_times_file = "../../../input/MU_firing_times_real.txt"
#firing_times_file = "../../../input/MU_firing_times_immediately.txt"

# partitioning
# ------------
# this has to match the total number of processes
n_subdomains_x = 1
n_subdomains_y = 1
n_subdomains_z = 1

# stride for sampling the 3D elements from the fiber data
# here any number is possible
sampling_stride_x = 1
sampling_stride_y = 1
sampling_stride_z = 50
sampling_stride_fat = 1

# how much of the multidomain mesh is used for elasticity
sampling_factor_elasticity_x = 0.5    
sampling_factor_elasticity_y = 0.5
sampling_factor_elasticity_z = 0.5
sampling_factor_elasticity_fat_y = 0.5

# scenario name for log file
scenario_name = ""

# functions, here, Am, Cm and Conductivity are constant for all fibers and MU's
# These functions can be redefined differently in a custom variables script
def get_am(fiber_no, mu_no):
  return Am

def get_cm(fiber_no, mu_no):
  return Cm
  
def get_conductivity(fiber_no, mu_no):
  return Conductivity

def get_specific_states_call_frequency(fiber_no, mu_no):
  return stimulation_frequency

def get_specific_states_frequency_jitter(fiber_no, mu_no):
  return [0]

def get_specific_states_call_enable_begin(fiber_no, mu_no):
  return activation_start_time


muscle1_extent = [3.0, 3.0, 14.8] # [cm, cm, cm]
tendon_length = 1.2 # cm
muscle2_extent = [3.0, 3.0, 14.8] # [cm, cm, cm]

n_elements_muscle1 = [2, 2, 4] # linear elements. each qudaratic element uses the combined nodes of 8 linear elements
n_elements_muscle2 = [2, 2, 4]
n_points_whole_fiber = 40
n_fibers_x = 4
n_fibers_y = 4



## currently undefined 
maximum_number_of_threads = 1
use_aovs_memory_layout = True
enable_surface_emg = False
fast_monodomain_solver_optimizations = True


# REMOVE ? 
fiber_file = None # create_partitioned_meshes_for_settings not needed? reads fibers
load_fiber_data = False
include_global_node_positions = False

# further internal variables that will be set by the helper.py script and used in the config in settings_fibers_emg.py
n_fibers_total = None
n_subdomains_xy = None
own_subdomain_coordinate_x = 0 # TODO fix this for parallelization
own_subdomain_coordinate_y = 0 # TODO fix this for parallelization
own_subdomain_coordinate_z = 0 # TODO fix this for parallelization
n_points_3D_mesh_global_x = None
n_points_3D_mesh_global_y = None
n_points_3D_mesh_global_z = None
output_writer_fibers = None
output_writer_emg = None
output_writer_0D_states = None
states_output = False
parameters_used_as_algebraic = None
parameters_used_as_constant = None
parameters_initial_values = None
output_algebraic_index = None
output_state_index = None
nodal_stimulation_current = None
fiber_file_handle = None
fibers = None
firing_times = None
n_fibers_per_subdomain_x = None
n_fibers_per_subdomain_y = None
n_points_per_subdomain_z = None
z_point_index_start = None
z_point_index_end = None
meshes = {}
potential_flow_dirichlet_bc = None
elasticity_dirichlet_bc = None
elasticity_neumann_bc = None
fibers_on_own_rank = None
n_fiber_nodes_on_subdomain = None
fiber_start_node_no = None
generate_linear_3d_mesh = True
generate_quadratic_3d_mesh = True
fat_mesh_n_points = None
fat_mesh_n_points_global = None
local_range_i = None
local_range_k = None
relative_factors = None
n_compartments = None
nx = None
ny = None
nz = None
constant_body_force = None
pmax = None
bottom_traction = None
states_initial_values = []
fix_bottom = False




# general parameters
# -----------------------------
rho = 10                    # [1e-4 kg/cm^3] density of the muscle (density of water)
# Mooney-Rivlin parameters [c1,c2,b,d] of c1*(Ibar1 - 3) + c2*(Ibar2 - 3) + b/d (λ - 1) - b*ln(λ)
# Heidlauf13: [6.352e-10 kPa, 3.627 kPa, 2.756e-5 kPa, 43.373] = [6.352e-11 N/cm^2, 3.627e-1 N/cm^2, 2.756e-6 N/cm^2, 43.373], pmax = 73 kPa = 7.3 N/cm^2
# Heidlauf16: [3.176e-10 N/cm^2, 1.813 N/cm^2, 1.075e-2 N/cm^2, 9.1733], pmax = 7.3 N/cm^2
c1 = 3.176e-10              # [N/cm^2]
c2 = 1.813                  # [N/cm^2]
b  = 1.075e-2               # [N/cm^2] anisotropy parameter
d  = 9.1733                 # [-] anisotropy parameter
# for debugging, b = 0 leads to normal Mooney-Rivlin

material_parameters = [c1, c2, b, d]   # material parameters


#--------------------------------

import os
import numpy as np
input_directory = os.path.join(os.environ.get('OPENDIHU_HOME', '../../../../../'), "examples/electrophysiology/input")

# cellml_file = input_directory+"/hodgkin_huxley_1952.c"
cellml_file = input_directory+"/hodgkin_huxley-razumova.cellml"
fiber_distribution_file = input_directory+"/MU_fibre_distribution_multidomain_67x67_100.txt"
 

# neurons and sensors
# -------------------

# muscle spindles
n_muscle_spindles = 3
muscle_spindle_cellml_file = input_directory+"/mileusenic_muscle_spindle_equilibrium.cellml"
muscle_spindle_mappings = {
  # first: don't define any mappings. then fix all warnings and errors

  # define the order of the parameters for initial values and callbacks
  # opendihu                    cellml
  ("parameter", 0):            "modell/L",                         # Spinndle Stretch (= fiber stretch)
  ("parameter", 1):            "modell/L_dot",                     # Spinndle velocity (= dL/dt)
  ("parameter", 2):            "modell/L_ddot",                    # Spinndle Acceleration (= d^2L/dt^2)
  ("parameter", 3):            "modell/gamma_sta",                 # Static fusimotor drive
  ("parameter", 4):            "modell/gamma_dyn",                 # Dynamic fusimotor drive
  # we first have to define the cellml constants as parameters. primary_afferent ist schon ein state. 
  # opendihu                    cellml
  ("connectorSlot", "ms_out"): "modell/primary_afferent",          # Ia Afferent firing frequency
  ("connectorSlot", "ms_in0"): "modell/L",                         # Spinndle Stretch (= fiber stretch)
  ("connectorSlot", "ms_in1"): "modell/L_dot",                     # Spinndle velocity (= dL/dt)
  ("connectorSlot", "ms_in2"): "modell/L_ddot",                    # Spinndle Acceleration (= d^2L/dt^2)
  ("connectorSlot", "ms_in3"): "modell/gamma_sta",                 # Static fusimotor drive
  ("connectorSlot", "ms_in4"): "modell/gamma_dyn",                 # Dynamic fusimotor drive
}
muscle_spindle_parameters_initial_values = [0, 0, 0, 5, 0]    # [L, L_dot, L_ddot, gamma_sta, gamma_dyn] (see above)
muscle_spindle_delay = 30             # [ms] signal delay between muscle spindle model and motoneuron model

# load cortical input values
cortical_input_file = input_directory+"/cortical_input_realistic.txt"
cortical_input = np.genfromtxt(cortical_input_file, delimiter=",")

# motor neurons
n_motoneurons = 10
motoneuron_cellml_file = input_directory+"/WSBM_1457_MN_Cisi_Kohn_2008.cellml"
motoneuron_mappings = {
  ("parameter", 0):            "motor_neuron/drive",   # stimulation
  ("parameter", 1):            "lumped_geometry_parameters/C_m",
  ("parameter", 2):            "lumped_geometry_parameters/R_i",
  ("parameter", 3):            "lumped_geometry_parameters/R_md",
  ("parameter", 4):            "lumped_geometry_parameters/R_ms",
  ("parameter", 5):            "lumped_geometry_parameters/l_d",
  ("parameter", 6):            "lumped_geometry_parameters/l_s",
  ("parameter", 7):            "lumped_geometry_parameters/r_d",
  ("parameter", 8):            "lumped_geometry_parameters/r_s",
  ("connectorSlot", "mn_out"): "motor_neuron/V_s",     # voltage
  ("connectorSlot", "mn_in"):  "motor_neuron/drive",   # stimulation
}
# initial values for the parameters, either once for all instances, or for all instances in array of struct ordering with nParameters_ parameters per instance: [inst0p0, inst0p1, ... inst0pn, inst1p0, inst1p1, ...]
motoneuron_parameters_initial_values = []
# Parameter ranges for MN Paramters [smallest MN value, largest MN value]
C_m = [1, 1]
Ri = [0.07, 0.07]
r_d = [20.75e-4, 46.25e-4]
l_d = [0.55, 1.06]
r_s = [38.75e-4, 56.5e-4]
l_s = [77.5e-4, 113e-4]
R_ms = [1.15, 0.65]
R_md = [14.4, 6.05]
gg = np.linspace(0,1,num = n_motoneurons)
#
for mu_idx in range(n_motoneurons):
  # compute a scaling factor that distributes variables exponentially between a lower and an upper value
  factor = np.exp(np.log(100)*gg[mu_idx-1])/100  # 100^(gg-1) in [0,1]
  # add parameter values for motoneuron mu_idx
  motoneuron_parameters_initial_values += [0.0, C_m[0]+factor*(C_m[1]-C_m[0]), Ri[0]+factor*(Ri[1]-Ri[0]), 
                                          R_md[0]+factor*(R_md[1]-R_md[0]), R_ms[0]+factor*(R_ms[1]-R_ms[0]), 
                                          l_d[0]+factor*(l_d[1]-l_d[0]), l_s[0]+factor*(l_s[1]-l_s[0]), 
                                          r_d[0]+factor*(r_d[1]-r_d[0]), r_s[0]+factor*(r_s[1]-r_s[0])]


# golgi tendon organs
n_golgi_tendon_organs = 3
golgi_tendon_organ_cellml_file = input_directory+"/hodgkin_huxley_1952.cellml"
golgi_tendon_organ_mappings = {
  ("parameter", 0):            "membrane/i_Stim",   # stimulation
  ("connectorSlot", "gt_out"): "membrane/V",        # voltage
  ("connectorSlot", "gt_in"):  "membrane/i_Stim",   # stimulation
}
golgi_tendon_organ_parameters_initial_values = [0]    # [i_Stim]
golgi_tendon_organ_delay = 300




tendon_force = 1 # TODO
main_constant_body_force = [0, 0, 0]

### mechanics
def muscle1_update_neumann_boundary_conditions(t,m_list):
  [mx,my,mz] = m_list # elements

  # Neumann BC at bottom nodes, force to the right (will be divided by area by config flag)
  # muscle mesh
  muscle1_elasticity_neumann_bc = [{"element": (mz-1)*mx*my + j*mx + i, "constantVector": (0, 0, tendon_force), "face": "2+"} for j in range(my) for i in range(mx)]

  config = {
    "inputMeshIsGlobal": True,
    "divideNeumannBoundaryConditionValuesByTotalArea": True,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
    "neumannBoundaryConditions": muscle1_elasticity_neumann_bc
  }
  print("t: {}, tendon force: {}".format(t, tendon_force))

  return config


def muscle2_update_neumann_boundary_conditions(t,m_list):
  [mx,my,mz] = m_list # elements

  # Neumann BC at bottom nodes, force to the left (will be divided by area by config flag)
  # muscle mesh
  muscle2_elasticity_neumann_bc = [{"element": 0*mx*my + j*mx + i, "constantVector": (0,0,-tendon_force), "face": "2-"} for j in range(my) for i in range(mx)]

  config = {
    "inputMeshIsGlobal": True,
    "divideNeumannBoundaryConditionValuesByTotalArea": True,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
    "neumannBoundaryConditions": muscle2_elasticity_neumann_bc
  }
  print("t: {}, tendon force: {}".format(t, tendon_force))

  return config



def muscle1_postprocess(data):
    print('callback 1')
def muscle2_postprocess(data):
    print('callback 2')