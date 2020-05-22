# debugging scenario that produces a fast simulation, not really physiological conditions (Am is too high therefore more diffusion than normal)
# scenario name for log file
scenario_name = "debug"

# material parameters
# --------------------
sigma_f = 8.93              # [mS/cm] conductivity in fiber direction (f)
sigma_xf = 0                # [mS/cm] conductivity in cross-fiber direction (xf)
sigma_e_f = 6.7             # [mS/cm] conductivity in extracellular space, fiber direction (f)
sigma_e_xf = 3.35           # [mS/cm] conductivity in extracellular space, cross-fiber direction (xf) / transverse

Am = 500.0                  # [cm^-1] surface area to volume ratio, actual values will be set by motor_units, not by this variable
Cm = 0.58                   # [uF/cm^2] membrane capacitance, (1 = fast twitch, 0.58 = slow twitch), actual values will be set by motor_units, not by this variable

Am = 1.0  # set Am lower such that there is more diffusion
# diffusion prefactor = Conductivity/(Am*Cm)

# timing and activation parameters
# -----------------
import random
random.seed(0)  # ensure that random numbers are the same on every rank
#   fiber_no: center MU around this fiber
#   standard_deviation [-]: relative to muscle diameter, 
#   maximum [-]: create f_r as gaussian with standard_deviation and maximum around the fiber given in fiber_no
#   radius: [μm], activation_start_time: [s], stimulation frequency [Hz], jitter [-]
# exponential distribution: low number of fibers per MU, slow twitch (type I), activated first --> high number of fibers per MU, fast twitch (type II), activated last -->
motor_units = [
  {"fiber_no": 10, "standard_deviation": 0.2, "maximum": 0.2, "radius": 40.00, "cm": 0.58, "activation_start_time": 0.0, "stimulation_frequency": 23.92, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},    # low number of fibers
  {"fiber_no": 20, "standard_deviation": 0.2, "maximum": 0.2, "radius": 42.35, "cm": 0.58, "activation_start_time": 0.2, "stimulation_frequency": 23.36, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"fiber_no": 30, "standard_deviation": 0.2, "maximum": 0.2, "radius": 45.00, "cm": 0.58, "activation_start_time": 0.4, "stimulation_frequency": 23.32, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"fiber_no": 40, "standard_deviation": 0.2, "maximum": 0.2, "radius": 48.00, "cm": 0.58, "activation_start_time": 0.6, "stimulation_frequency": 22.46, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"fiber_no": 55, "standard_deviation": 0.2, "maximum": 0.2, "radius": 51.42, "cm": 0.58, "activation_start_time": 0.8, "stimulation_frequency": 20.28, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"fiber_no": 60, "standard_deviation": 0.2, "maximum": 0.2, "radius": 55.38, "cm": 0.58, "activation_start_time": 1.0, "stimulation_frequency": 16.32, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"fiber_no": 70, "standard_deviation": 0.2, "maximum": 0.2, "radius": 60.00, "cm": 0.58, "activation_start_time": 1.2, "stimulation_frequency": 12.05, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"fiber_no": 80, "standard_deviation": 0.2, "maximum": 0.2, "radius": 65.45, "cm": 1.00, "activation_start_time": 1.4, "stimulation_frequency": 10.03, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"fiber_no": 50, "standard_deviation": 0.2, "maximum": 0.2, "radius": 72.00, "cm": 1.00, "activation_start_time": 1.6, "stimulation_frequency": 8.32,  "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"fiber_no": 25, "standard_deviation": 0.2, "maximum": 0.2, "radius": 80.00, "cm": 1.00, "activation_start_time": 1.8, "stimulation_frequency": 7.66,  "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},    # high number of fibers
]
motor_units = motor_units[0:2]  # only some motor units

# solvers
# -------
multidomain_solver_type = "gmres"          # solver for the multidomain problem
#multidomain_preconditioner_type = "bjacobi"   # preconditioner
#multidomain_preconditioner_type = "boomeramg"   # preconditioner
multidomain_preconditioner_type = "euclid"   # preconditioner

multidomain_alternative_solver_type = "gmres"            # alternative solver, used when normal solver diverges
multidomain_alternative_preconditioner_type = "none"    # preconditioner of the alternative solver

# set initial guess to zero for direct solver
initial_guess_nonzero = "lu" not in multidomain_solver_type 

# if using boomeramg, tolerances cannot be as low as 1e-10, otherwise it becomes unstable
multidomain_absolute_tolerance = 1e-15    # absolute residual tolerance for the multidomain solver
multidomain_relative_tolerance = 1e-15    # relative residual tolerance for the multidomain solver
theta = 1.0                               # weighting factor of implicit term in Crank-Nicolson scheme, 0.5 gives the classic, 2nd-order Crank-Nicolson scheme, 1.0 gives implicit euler
use_symmetric_preconditioner_matrix = True   # if the diagonal blocks of the system matrix should be used as preconditioner matrix
use_lumped_mass_matrix = False            # which formulation to use, the formulation with lumped mass matrix (True) is more stable but approximative, the other formulation (False) is exact but needs more iterations
multidomain_max_iterations = 1e3                         # maximum number of iterations
multidomain_alternative_solver_max_iterations = 1e4      # maximum number of iterations of the alternative solver

# timing parameters
# -----------------
end_time = 4000.0                   # [ms] end time of the simulation
stimulation_frequency = 100*1e-3    # [ms^-1] sampling frequency of stimuli in firing_times_file, in stimulations per ms, number before 1e-3 factor is in Hertz.
dt_0D = 1e-3                        # [ms] timestep width of ODEs (1e-3)
dt_multidomain = 1e-3               # [ms] timestep width of the multidomain solver
dt_splitting = 1e-3                 # [ms] overall timestep width of strang splitting (3e-3)
output_timestep_multidomain = dt_multidomain  # [ms] timestep for multidomain output
output_timestep_0D_states = dt_0D    # [ms] timestep for output files of 0D subcellular model states
end_time = 40*dt_0D

scenario_name = "{}_{}_dt{}_atol{}_rtol{}_theta{}_sym{}_lump{}".format(multidomain_solver_type, multidomain_preconditioner_type, dt_splitting, multidomain_absolute_tolerance, multidomain_relative_tolerance, theta, use_symmetric_preconditioner_matrix, use_lumped_mass_matrix)

# input files
#cellml_file = "../../input/hodgkin_huxley_1952.c"
cellml_file = "../../input/new_slow_TK_2014_12_08.c"
#fiber_file = "../../input/left_biceps_brachii_7x7fibers.bin"
#fiber_file = "../../input/small_5x5x11.bin"
fiber_file = "../../input/box_3x3x17.bin"
fat_mesh_file = fiber_file + "_fat.bin"
firing_times_file = "../../input/MU_firing_times_immediately.txt"


# stride for sampling the 3D elements from the fiber data
# a higher number leads to less 3D elements
sampling_stride_x = 1
sampling_stride_y = 1
sampling_stride_z = 1 #74   # faster, but stimulus does not propagate
sampling_stride_fat = 1

# other options
paraview_output = True
adios_output = False
exfile_output = False
python_output = False
states_output = True
disable_firing_output = False

# functions, here, Am, Cm and Conductivity are constant for all fibers and MU's
def get_am(mu_no):
  # get radius in cm, 1 μm = 1e-6 m = 1e-4*1e-2 m = 1e-4 cm
  r = motor_units[mu_no]["radius"]*1e-4
  # cylinder surface: A = 2*π*r*l, V = cylinder volume: π*r^2*l, Am = A/V = 2*π*r*l / (π*r^2*l) = 2/r
  return 2./r
  #return Am

def get_cm(mu_no):
  return motor_units[mu_no % len(motor_units)]["cm"]
  #return Cm
  
def get_conductivity(mu_no):
  return Conductivity

def get_specific_states_call_frequency(mu_no):
  stimulation_frequency = motor_units[mu_no % len(motor_units)]["stimulation_frequency"]
  return stimulation_frequency*1e-3

def get_specific_states_frequency_jitter(mu_no):
  return 0
  #return motor_units[mu_no % len(motor_units)]["jitter"]

def get_specific_states_call_enable_begin(mu_no):
  return motor_units[mu_no % len(motor_units)]["activation_start_time"]*1e3
