
# scenario name for log file
scenario_name = "neurons_multidomain"

# Fixed units in cellMl models:
# These define the unit system.
# 1 cm = 1e-2 m
# 1 ms = 1e-3 s
# 1 uA = 1e-6 A
# 1 uF = 1e-6 F
# 
# Fixed units in mechanics system
# 1 cm = 1e-2 m
# 1 ms = 1e-3 s
# 1 N
# 1 N/cm^2 = (kg*m*s^-2) / (1e-2 m)^2 = 1e4 kg*m^-1*s^-2 = 10 kPa
# (kg = N*s^2*m^-1) => N*ms^2*cm^-1 = N*(1e-3 s)^2 * (1e-2 m)^-1 = 1e-4 N*s^2*m^-1 = 1e-4 kg
# (kg/m^3) => 1 * 1e-4 kg * (1e-2 m)^-3 = 1e2 kg/m^3
# (m/s^2) => 1 cm/ms^2 = 1e-2 m * (1e-3 s)^-2 = 1e4 m*s^-2

# material parameters
# --------------------
# quantities in mechanics unit system

# parameters for precontraction
# -----------------------------
# load
precontraction_constant_body_force = (0,0,20*9.81e-4)   # [cm/ms^2], gravity constant for the body force
precontraction_bottom_traction = [0,0,0]        # [N]
constant_gamma = 0.3    # 0.3 works, the active stress will be pmax*constant_gamma

# parameters for prestretch
# -----------------------------
# load
prestretch_constant_body_force = (0,0,-9.81e-4)   # [cm/ms^2], gravity constant for the body force
prestretch_bottom_traction = [0,0,-10]        # [N]  (-30 also works)

# parameters for multidomain simulation
# load
multidomain_constant_body_force = (0,0,-9.81e-4)   # [cm/ms^2], gravity constant for the body force
multidomain_bottom_traction = [0,0,-10]        # [N]  (-30 works)

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
#b = 0

material_parameters = [c1, c2, b, d]   # material parameters
pmax = 7.3                  # [N/cm^2] maximum isometric active stress

# timing and activation parameters
# -----------------
# motor unit parameters similar to paper Klotz2019 "Modelling the electrical activity of skeletal muscle tissue using a multi‐domain approach"
# however, values from paper fail for mu >= 6, then stimulus gets reflected at the ends of the muscle, therefore fiber radius is set to <= 55

import random
random.seed(0)  # ensure that random numbers are the same on every rank
import scipy
import numpy as np

n_fibers_in_fiber_file = 81
n_motor_units = 3   # number of motor units

motor_units = []
for mu_no in range(n_motor_units):

  # capacitance of the membrane
  if mu_no <= 0.7*n_motor_units:
    cm = 0.58    # slow twitch (type I)
  else:
    cm = 1.0     # fast twitch (type II)

  # fiber radius between 40 and 55 [μm]
  min_value = 40
  max_value = 55

  # ansatz value(i) = c1 + c2*exp(i),
  # value(0) = min = c1 + c2  =>  c1 = min - c2
  # value(n-1) = max = min - c2 + c2*exp(n-1)  =>  max = min + c2*(exp(n-1) - 1)  =>  c2 = (max - min) / (exp(n-1) - 1)
  c2 = (max_value - min_value) / (1.02**(n_motor_units-1) - 1)
  c1 = min_value - c2
  radius = c1 + c2*1.02**(mu_no)

  # standard_deviation
  min_value = 0.1
  max_value = 0.6
  c2 = (max_value - min_value) / (1.02**(n_motor_units-1) - 1)
  c1 = min_value - c2
  standard_deviation = c1 + c2*1.02**mu_no
  maximum = 10.0/n_motor_units*standard_deviation

  # exponential distribution: low number of fibers per MU, slow twitch (type I), activated first --> high number of fibers per MU, fast twitch (type II), activated last
  motor_units.append(
  {
    "fiber_no":              random.randint(0,n_fibers_in_fiber_file),  # [-] fiber from input files that is the center of the motor unit domain
    "maximum":               maximum,                # [-] maximum value of f_r, create f_r as gaussian with standard_deviation and maximum around the fiber 
    "standard_deviation":    standard_deviation,     # [-] standard deviation of f_r
    "radius":                radius,                 # [μm] parameter for motor unit: radius of the fiber, used to compute Am
    "cm":                    cm,                     # [uF/cm^2] parameter Cm
  })

#motor_units = motor_units[0:1]  # for debugging, only 1 motor unit

# solvers
# -------
# potential flow
potential_flow_solver_type = "gmres"        # solver and preconditioner for an initial Laplace flow on the domain, from which fiber directions are determined
potential_flow_preconditioner_type = "none" # preconditioner

# multidomain
multidomain_solver_type = "gmres"          # solver for the multidomain problem
multidomain_preconditioner_type = "euclid"   # preconditioner
multidomain_max_iterations = 1e3                         # maximum number of iterations

multidomain_alternative_solver_type = "gmres"            # alternative solver, used when normal solver diverges
multidomain_alternative_preconditioner_type = "euclid"    # preconditioner of the alternative solver
multidomain_alternative_solver_max_iterations = 1e4      # maximum number of iterations of the alternative solver

multidomain_absolute_tolerance = 1e-15 # absolute residual tolerance for the multidomain solver
multidomain_relative_tolerance = 1e-15 # absolute residual tolerance for the multidomain solver

initial_guess_nonzero = "lu" not in multidomain_solver_type   # set initial guess to zero for direct solver
theta = 1.0                               # weighting factor of implicit term in Crank-Nicolson scheme, 0.5 gives the classic, 2nd-order Crank-Nicolson scheme, 1.0 gives implicit euler
use_symmetric_preconditioner_matrix = True   # if the diagonal blocks of the system matrix should be used as preconditioner matrix
use_lumped_mass_matrix = False            # which formulation to use, the formulation with lumped mass matrix (True) is more stable but approximative, the other formulation (False) is exact but needs more iterations

# elasticity
elasticity_solver_type = "lu"
elasticity_preconditioner_type = "none"
snes_max_iterations = 34                  # maximum number of iterations in the nonlinear solver
snes_rebuild_jacobian_frequency = 1       # how often the jacobian should be recomputed, -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again 
snes_relative_tolerance = 1e-5      # relative tolerance of the nonlinear solver
snes_absolute_tolerance = 1e-4      # absolute tolerance of the nonlinear solver
relative_tolerance = 1e-10           # relative tolerance of the residual of the linear solver
absolute_tolerance = 1e-10          # absolute tolerance of the residual of the linear solver

theta = 1.0                               # weighting factor of implicit term in Crank-Nicolson scheme, 0.5 gives the classic, 2nd-order Crank-Nicolson scheme, 1.0 gives implicit euler
use_symmetric_preconditioner_matrix = True   # if the diagonal blocks of the system matrix should be used as preconditioner matrix
use_lumped_mass_matrix = False            # which formulation to use, the formulation with lumped mass matrix (True) is more stable but approximative, the other formulation (False) is exact but needs more iterations

# timing parameters
# -----------------
end_time = 5_000.0                   # [ms] end time of the simulation
dt_0D = 1e-3                        # [ms] timestep width of ODEs (1e-3)
dt_multidomain = 1e-3               # [ms] timestep width of the multidomain solver, i.e. the diffusion
dt_splitting = dt_multidomain       # [ms] timestep width of strang splitting between 0D and multidomain, this is the same as the dt_multidomain, because we do not want to subcycle for the diffusion part
dt_elasticity = 1e-1                # [ms] time step width of elasticity solver
#dt_elasticity = 1e-2                # [ms] time step width of elasticity solver

dt_neurons = 1e-2                   # [ms] same timestep width for all neuron solvers
dt_golgi_tendon_organs = dt_neurons # [ms] timestep width of cellml solver of golgi tendon organs
dt_muscle_spindles     = 1e-3       # [ms] timestep width of cellml solver of muscle spindles
dt_interneuron         = dt_neurons # [ms] timestep width of the cellml solver for interneurons
dt_motoneuron          = dt_neurons # [ms] timestep width of the cellml solver for motoneurons

dt_neuron_transfer     = dt_elasticity  # [ms] interval when to call callback functions and transfer values between CellML models, increase this to speed up the simulation
#dt_neuron_transfer     = dt_neurons  # [ms] interval when to call callback functions and transfer values between CellML models, increase this to speed up the simulation

output_timestep_multidomain = 2     # [ms] timestep for multidomain solver output
output_timestep_elasticity = 1      # [ms] timestep for elasticity output files
output_timestep_neurons = 1         # [ms] timestep for output of files for all sensor organs and neurons
output_timestep_motoneuron = 0.2    # [ms] timestep for output of files for motoneuron

#output_timestep_multidomain = dt_elasticity
#output_timestep_elasticity = dt_elasticity

# input files
# -----------
import os
input_directory   = os.path.join(os.environ["OPENDIHU_HOME"], "examples/electrophysiology/input")
#cellml_file       = input_directory+"/new_slow_TK_2014_12_08.c"
cellml_file       = input_directory+"/hodgkin_huxley-razumova.cellml"

fiber_file        = input_directory+"/left_biceps_brachii_9x9fibers_b.bin"  # this is a variant of 9x9fibers with a slightly different mesh that somehow works better
#fiber_file        = input_directory+"/left_biceps_brachii_13x13fibers.bin"
fat_mesh_file     = fiber_file + "_fat.bin"
firing_times_file = input_directory+"/MU_firing_times_always.txt"    # use setSpecificStatesCallEnableBegin and setSpecificStatesCallFrequency
firing_times_file = input_directory+"/MU_firing_times_once.txt"    # use setSpecificStatesCallEnableBegin and setSpecificStatesCallFrequency
fiber_distribution_file = input_directory+"/MU_fibre_distribution_10MUs.txt"

# stride for meshes
# -----------------

# stride for sampling the 3D elements from the fiber data
# a higher number leads to less 3D elements
# If you change this, delete the compartment_relative_factors.* files, they have to be generated again.
sampling_stride_x = 1 
sampling_stride_y = 1 
sampling_stride_z = 50
local_sampling_stride_z = 1 
sampling_stride_fat = 1 

# how much of the multidomain mesh is used for elasticity
sampling_factor_elasticity_x = 0.7 
sampling_factor_elasticity_y = 0.7 
sampling_factor_elasticity_z = 0.3 
sampling_factor_elasticity_fat_y = 0.5 

# other options
paraview_output = True
adios_output = False
exfile_output = False
python_output = False
states_output = False    # if also the subcellular states should be output, this produces large files, set output_timestep_0D_states
show_linear_solver_output = False    # if every solve of multidomain diffusion should be printed
disable_firing_output = True   # if information about firing of MUs should be printed

# neurons and sensors
# -------------------

# muscle spindles
n_muscle_spindles = 3
muscle_spindle_cellml_file = input_directory+"/mileusenic_muscle_spindle_equilibrium.cellml"
muscle_spindle_mappings = {
  ("parameter", 0):            "modell/L",                         # Spinndle Stretch (= fiber stretch)
  ("parameter", 1):            "modell/L_dot",                     # Spinndle velocity (= dL/dt)
  ("parameter", 2):            "modell/L_ddot",                    # Spinndle Acceleration (= d^2L/dt^2)
  ("parameter", 3):            "modell/gamma_sta",                 # Static fusimotor drive
  ("parameter", 4):            "modell/gamma_dyn",                 # Dynamic fusimotor drive
  ("connectorSlot", "ms_out"): "modell/primary_afferent",          # Ia Afferent firing frequency
  ("connectorSlot", "ms_in0"): "modell/L",                         # Spinndle Stretch (= fiber stretch)
  ("connectorSlot", "ms_in1"): "modell/L_dot",                     # Spinndle velocity (= dL/dt)
  ("connectorSlot", "ms_in2"): "modell/L_ddot",                    # Spinndle Acceleration (= d^2L/dt^2)
  ("connectorSlot", "ms_in3"): "modell/gamma_sta",                 # Static fusimotor drive
  ("connectorSlot", "ms_in4"): "modell/gamma_dyn",                 # Dynamic fusimotor drive
}
muscle_spindle_parameters_initial_values = [0, 0, 0, 0, 0]    # [L, L_dot, L_ddot, gamma_sta, gamma_dyn]
muscle_spindle_delay = 30             # [ms] signal delay between muscle spindle model and motoneuron model

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

# interneurons
n_interneurons = 3
interneuron_cellml_file = input_directory+"/hodgkin_huxley_1952.cellml"
interneuron_mappings = {
  ("parameter", 0):            "membrane/i_Stim",   # stimulation
  ("parameter", 1):            "membrane/Cm",       # stimulation
  ("connectorSlot", "in_out"): "membrane/V",        # voltage
  ("connectorSlot", "in_in"):  "membrane/i_Stim",   # stimulation
}
# initial values for the parameters, either once for all instances, or for all instances in array of struct ordering with nParameters_ parameters per instance: [inst0p0, inst0p1, ... inst0pn, inst1p0, inst1p1, ...]
interneuron_parameters_initial_values = []    # [i_Stim, Cm]
for i in range(n_interneurons):
  
  # compute a scaling factor that runs exponentially from min_factor to max_factor
  min_factor = 0.5
  max_factor = 2.0
  
  # ansatz scaling_factor(i) = c1 + c2*exp(i),
  # scaling_factor(0) = min = c1 + c2  =>  c1 = min - c2
  # scaling_factor(n-1) = max = min - c2 + c2*exp(n-1)  =>  max = min + c2*(exp(n-1) - 1)  =>  c2 = (max - min) / (exp(n-1) - 1)
  c2 = (max_factor - min_factor) / (np.exp(n_interneurons-1) - 1)
  c1 = min_factor - c2
  scaling_factor = c1 + c2*np.exp(i)
  
  # add parameter values for motoneuron i
  interneuron_parameters_initial_values += [0.0, 1*scaling_factor]
  
#print("interneuron_parameters_initial_values: {}".format(interneuron_parameters_initial_values))
  
# motor neurons
n_motoneurons = 3
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

# neuron and mapping callbacks
# ------------------------------

def callback_muscle_spindles_input(input_values, output_values, current_time, slot_nos, buffer):
  """
  Callback function that transform a number of input_values to a number of output_values.
  This function gets called by a MapDofs object.
  :param input_values: (list of float values) The input values from the slot as defined in the MapDofs settings.
  :param output_values: (list of list of float values) output_values[slotIndex][valueIndex]
                        The output values buffer, potentially for multiple slots.
                        Initially, this is a list of the form [[None, None, ..., None]] with the size matching 
                        the number of required output values. The function should set some of the entries to a computed value.
                        The entries that are not None will be set in the output slot at the dofs defined by MapDofs.
  :param current_time:  Current simulation time.
  :param slot_nos:      List of [fromSlotNo, toSlotNo, fromArrayIndex, toArrayIndex].
  :param buffer:        A persistent helper buffer. This variable can be set to anything and will be provided back to 
                        this function every time. Using this buffer, it is possible to implement a time delay of signals.
  """
  # map from λ in the 3D mesh to muscle spindles model input
  # input_values contains a list of λ values from the muscle spindle nodes
  
  # get number of input and output values
  n_input_values = len(input_values)      # = n_muscle_spindles
  n_output_values = len(output_values[0]) # = n_muscle_spindles (per output slot if there are multiple)
  
  # initialize buffer, buffer is needed to compute velocity and acceleration
  if "stretch" not in buffer:
    buffer["stretch"]      = [1 for _ in range(n_input_values)]
    buffer["velocity"]     = [0 for _ in range(n_input_values)]
    buffer["acceleration"] = [0 for _ in range(n_input_values)]
    buffer["current_time"] = 0
    delta_t = 1
  else:  
    delta_t  = buffer["current_time"] - current_time
    
  # normal stimulation: parse input value of λ
  for i in range(n_input_values):
    stretch = input_values[i]
    
    # the first time when the stretch is not yet computed it has a value of 0, set to 1
    if stretch == 0:
      stretch = 1
      
    # compute velocity and acceleration by difference quotient
    velocity     = (buffer["stretch"][i] - stretch) / delta_t
    acceleration = (buffer["velocity"][i] - velocity) / delta_t
    buffer["stretch"][i]      = stretch
    buffer["velocity"][i]     = velocity
    buffer["acceleration"][i] = acceleration
    print("spindle no. {}, stretch: {}, v: {}, a: {}".format(i, stretch, velocity, acceleration))
      
    # scale value for muscle spindle model 
    output_values[0][i] = stretch                 # fiber stretch, output units: [L0] with L0=cm
    output_values[1][i] = velocity*1e-3           # velocity, output units: [cm*s^-1], 1e-3 cm*ms^-1 = 1e-3 cm*(1e-3s)^-1 = 1e-3 cm*1e3*s^-1 = 1 cm*s^-1
    output_values[2][i] = acceleration*1e-6       # acceleration, output units: [cm*s^-2], 1e-6 cm*ms^-2 = 1e-6 cm*(1e-3s)^-2 = 1e-6 cm*(1e6)*s^-2 = 1 cm*s^-2
    output_values[3][i] = 0                       # static fusimotor drive
    output_values[4][i] = 0                       # dynamic fusimotor drive
  
  # store current time
  buffer["current_time"] = current_time
  
  # artifical muscle spindle input for debugging
  if False:
    for i in range(n_output_values):
      T = 1000   # [ms] cycle duration of stimulus 
      output_values[0][i] = np.sin(current_time / T * np.pi)**2 * 0.001
      output_values[0][i] = 1
      
      # "ms_in0" -> "modell/L",                         # Spinndle Stretch (= fiber stretch)
      # "ms_in1" -> "modell/L_dot",                     # Spinndle velocity (= dL/dt)
      # "ms_in2" -> "modell/L_ddot",                    # Spinndle Acceleration (= d^2L/dt^2)
      # "ms_in3" -> "modell/gamma_sta",                 # Static fusimotor drive
      # "ms_in4" -> "modell/gamma_dyn",                 # Dynamic fusimotor drive
      output_values[0][i] = 1 + np.sin(current_time / T * np.pi) * 0.1                  # Stretch    
      output_values[1][i] = np.cos(current_time / T * np.pi) * 0.1 * (np.pi/T)          # Velocity
      output_values[2][i] = np.sin(current_time / T * np.pi) * (-0.1) * (np.pi/T)**2    # Acceleration
      output_values[3][i] = 0     # static fusimotor drive
      output_values[4][i] = 0     # dynamic fusimotor drive
  
def callback_muscle_spindles_to_motoneurons(input_values, output_values, current_time, slot_nos, buffer):
  """
  Callback function that transform a number of input_values to a number of output_values.
  This function gets called by a MapDofs object.
  :param input_values: (list of float values) The input values from the slot as defined in the MapDofs settings.
  :param output_values: (list of list of float values) output_values[slotIndex][valueIndex]
                        The output values buffer, potentially for multiple slots.
                        Initially, this is a list of the form [[None, None, ..., None]] with the size matching 
                        the number of required output values. The function should set some of the entries to a computed value.
                        The entries that are not None will be set in the output slot at the dofs defined by MapDofs.
  :param current_time:  Current simulation time.
  :param slot_nos:      List of [fromSlotNo, toSlotNo, fromArrayIndex, toArrayIndex].
  :param buffer:        A persistent helper buffer. This variable can be set to anything and will be provided back to 
                        this function every time. Using this buffer, it is possible to implement a time delay of signals.
  """
  # mapping muscle spindles output -> motor neuron signals, delay signals by muscle_spindle_delay
  
  # get number of input and output values
  n_input_values = len(input_values)      # = n_muscle_spindles
  n_output_values = len(output_values[0]) # = n_muscle_spindles (per output slot if there are multiple)
  
  # initialize buffer the first time, buffer later stores time of last activation, to model signal delay
  if 0 not in buffer:
    for muscle_spindle_index in range(n_input_values):
      buffer[muscle_spindle_index] = None
  
  # loop over muscle spindles
  for muscle_spindle_index in range(n_input_values):
    
    # determine spike by threshold
    if input_values[muscle_spindle_index] > 0:
      buffer[muscle_spindle_index] = current_time    # store time of last activation in buffer
      
    # if there has been a stimulation so far
    if buffer[muscle_spindle_index] is not None:
      
      # convolute Dirac delta, kernel is a shifted and scaled gaussian
      t_delay = muscle_spindle_delay             # [ms] delay of the signal
      gaussian_std_dev = 10                      # [ms] width of the gaussian curve
      convolution_kernel = lambda t: scipy.stats.norm.pdf(t, loc=t_delay, scale=gaussian_std_dev)*np.sqrt(2*np.pi)*gaussian_std_dev
      delayed_signal = convolution_kernel(current_time - buffer[muscle_spindle_index]) * 5
        
      # loop over output values and set all to the computed signal, cut off at 1e-5
      if delayed_signal > 1e-5:
        #print("muscle spindle t: {}, last_activation: {}, computed delayed_signal: {}".format(current_time, buffer[muscle_spindle_index], delayed_signal))
        output_values[0][muscle_spindle_index] = delayed_signal
      else:
        output_values[0][muscle_spindle_index] = None     # signal is below 1e-5, do not set any values
        
  print("muscle_spindles_to_motoneurons: {} -> {}".format(input_values, output_values))

def callback_golgi_tendon_organs_input(input_values, output_values, current_time, slot_nos, buffer):
  """
  Callback function that transform a number of input_values to a number of output_values.
  This function gets called by a MapDofs object.
  :param input_values: (list of float values) The input values from the slot as defined in the MapDofs settings.
  :param output_values: (list of list of float values) output_values[slotIndex][valueIndex]
                        The output values buffer, potentially for multiple slots.
                        Initially, this is a list of the form [[None, None, ..., None]] with the size matching 
                        the number of required output values. The function should set some of the entries to a computed value.
                        The entries that are not None will be set in the output slot at the dofs defined by MapDofs.
  :param current_time:  Current simulation time.
  :param slot_nos:      List of [fromSlotNo, toSlotNo, fromArrayIndex, toArrayIndex].
  :param buffer:        A persistent helper buffer. This variable can be set to anything and will be provided back to 
                        this function every time. Using this buffer, it is possible to implement a time delay of signals.
  """
  # map from T in the 3D mesh to golgi tendon organs
  
  # get number of input and output values
  n_input_values = len(input_values)      # = n_golgi_tendon_organs
  n_output_values = len(output_values[0]) # = n_golgi_tendon_organs (per output slot if there are multiple)
  
  # loop over golgi tendon organs
  for i in range(n_input_values):
    stress = input_values[i]
    
    output_values[0][i] = abs(stress)
    
  # artifical muscle spindle input for debugging
  if False:
    for i in range(n_output_values):
      T = 100   # [ms] cycle duration of stimulus 
      
      if 2000 + i*200 < current_time  <= 2100 + i*200:
        output_values[0][i] = np.sin(current_time / T * np.pi)**2 * 20
  
  print("traction at Golgi tendon organs: {}, output: {}".format(input_values, output_values))

def callback_golgi_tendon_organs_to_interneurons(input_values, output_values, current_time, slot_nos, buffer):
  """
  Callback function that transform a number of input_values to a number of output_values.
  This function gets called by a MapDofs object.
  :param input_values: (list of float values) The input values from the slot as defined in the MapDofs settings.
  :param output_values: (list of list of float values) output_values[slotIndex][valueIndex]
                        The output values buffer, potentially for multiple slots.
                        Initially, this is a list of the form [[None, None, ..., None]] with the size matching 
                        the number of required output values. The function should set some of the entries to a computed value.
                        The entries that are not None will be set in the output slot at the dofs defined by MapDofs.
  :param current_time:  Current simulation time.
  :param slot_nos:      List of [fromSlotNo, toSlotNo, fromArrayIndex, toArrayIndex].
  :param buffer:        A persistent helper buffer. This variable can be set to anything and will be provided back to 
                        this function every time. Using this buffer, it is possible to implement a time delay of signals.
  """
  # mapping Golgi tendon organs -> interneurons
  
  # get number of input and output values
  n_input_values = len(input_values)      # = n_golgi_tendon_organs
  n_output_values = len(output_values[0]) # = n_interneurons (per output slot if there are multiple)
  
  # sum up all input signals and set all outputs to this sum
  
  # collect sum of all Golgi tendon organs
  total_signal = 0
  
  # loop over Golgi tendon organs
  for golgi_tendon_organ_index in range(n_input_values):
    
    # if input is active
    if input_values[golgi_tendon_organ_index] > 20:
      total_signal += input_values[golgi_tendon_organ_index]
    
  # set same value to all connected interneurons
  for interneuron_index in range(n_output_values):
    output_values[0][interneuron_index] = total_signal * 0.01
    
  # initialize buffer the first time
  if 0 not in buffer:
    for golgi_tendon_organ_index in range(n_input_values):
      buffer[golgi_tendon_organ_index] = None
  
  # loop over golgi tendon organ inputs
  for golgi_tendon_organ_index in range(n_input_values):
    
    # determine spike by threshold
    if input_values[golgi_tendon_organ_index] > 20:
      buffer[golgi_tendon_organ_index] = current_time    # store time of last activation in buffer
      
    # if there has been a stimulation so far
    if buffer[golgi_tendon_organ_index] is not None:
      
      # convolute Dirac delta, kernel is a shifted and scaled gaussian
      t_delay = 0                   # [ms] delay of the signal
      gaussian_std_dev = 10         # [ms] width of the gaussian curve
      convolution_kernel = lambda t: scipy.stats.norm.pdf(t, loc=t_delay, scale=gaussian_std_dev)*np.sqrt(2*np.pi)*gaussian_std_dev
      delayed_signal = convolution_kernel(current_time - buffer[golgi_tendon_organ_index]) * 5
      # hodgkin-huxley fires from i_Stim(t) > 4 
        
      output_values[0][golgi_tendon_organ_index] = delayed_signal
  
  print("golgi_tendon_organs_to_interneurons input: {}, output: {}".format(input_values, output_values))

def callback_interneurons_to_motoneurons(input_values, output_values, current_time, slot_nos, buffer):
  """
  Callback function that transform a number of input_values to a number of output_values.
  This function gets called by a MapDofs object.
  :param input_values: (list of float values) The input values from the slot as defined in the MapDofs settings.
  :param output_values: (list of list of float values) output_values[slotIndex][valueIndex]
                        The output values buffer, potentially for multiple slots.
                        Initially, this is a list of the form [[None, None, ..., None]] with the size matching 
                        the number of required output values. The function should set some of the entries to a computed value.
                        The entries that are not None will be set in the output slot at the dofs defined by MapDofs.
  :param current_time:  Current simulation time.
  :param slot_nos:      List of [fromSlotNo, toSlotNo, fromArrayIndex, toArrayIndex].
  :param buffer:        A persistent helper buffer. This variable can be set to anything and will be provided back to 
                        this function every time. Using this buffer, it is possible to implement a time delay of signals.
  """
  # mapping interneurons -> input for motor neurons, i.e. signal delay from interneurons to motoneurons
  # the actual N->M mapping from interneurons to motoneurons is done by callback_motoneurons_input
  
  # get number of input and output values
  n_input_values = len(input_values)      # = n_interneurons
  n_output_values = len(output_values[0]) # = n_interneurons (per output slot if there are multiple)
  
  # initialize buffer the first time
  if 0 not in buffer:
    for interneuron_index in range(n_input_values):
      buffer[interneuron_index] = None
  
  # loop over interneurons
  for interneuron_index in range(n_input_values):
    
    # determine spike by threshold
    if input_values[interneuron_index] > 0:
      buffer[interneuron_index] = current_time    # store time of last activation in buffer
      
    # if there has been a stimulation so far
    if buffer[interneuron_index] is not None:
      
      # convolute Dirac delta, kernel is a shifted and scaled gaussian
      t_delay = golgi_tendon_organ_delay          # [ms] delay of the signal
      gaussian_std_dev = 10                       # [ms] width of the gaussian curve
      convolution_kernel = lambda t: scipy.stats.norm.pdf(t, loc=t_delay, scale=gaussian_std_dev)*np.sqrt(2*np.pi)*gaussian_std_dev
      delayed_signal = convolution_kernel(current_time - buffer[interneuron_index])   # motor neuron input should be around 1
        
      # loop over output values and set all to the computed signal, cut off at 1e-5
      if delayed_signal > 1e-5:
        #print("interneuron t: {}, last_activation: {}, computed delayed_signal: {}".format(current_time, buffer[interneuron_index], delayed_signal))
        output_values[0][interneuron_index] = delayed_signal
      else:
        output_values[0][interneuron_index] = None     # signal is below 1e-5, do not set any values
  
  print("interneurons_to_motoneurons input: {}, output: {}".format(input_values, output_values))
  
def callback_motoneurons_input(input_values, output_values, current_time, slot_nos, buffer):
  """
  Callback function that transform a number of input_values to a number of output_values.
  This function gets called by a MapDofs object.
  :param input_values: (list of float values) The input values from the slot as defined in the MapDofs settings.
  :param output_values: (list of list of float values) output_values[slotIndex][valueIndex]
                        The output values buffer, potentially for multiple slots.
                        Initially, this is a list of the form [[None, None, ..., None]] with the size matching 
                        the number of required output values. The function should set some of the entries to a computed value.
                        The entries that are not None will be set in the output slot at the dofs defined by MapDofs.
  :param current_time:  Current simulation time.
  :param slot_nos:      List of [fromSlotNo, toSlotNo, fromArrayIndex, toArrayIndex].
  :param buffer:        A persistent helper buffer. This variable can be set to anything and will be provided back to 
                        this function every time. Using this buffer, it is possible to implement a time delay of signals.
  """
  
  # map from delayed muscle spindle model outputs and delayed interneuron outputs to motoneuron inputs
  
  # get number of input and output values
  n_input_values = len(input_values)      # = n_muscle_spindles + n_interneurons
  n_output_values = len(output_values[0]) # = n_motoneurons (per output slot if there are multiple)
  
  # sum up all input signals and set all outputs to this sum
  
  # collect sum of all Golgi tendon organs
  total_signal = 0
  
  # loop over input values
  for input_index in range(n_input_values):
    total_signal += input_values[input_index] * 1e-3
    
  # add cortical input
  total_signal += 5e-3
  # motor neuron fires with ~14Hz if drive(t) = 5e-3
  
  # set same value to all connected motoneurons
  for motoneuron_index in range(n_output_values):
    output_values[0][motoneuron_index] = total_signal
    
  print("motoneurons input from spindles and interneurons: {}, resulting drive: {}".format(input_values, output_values))
  

# multidomain callbacks
# ----------------------
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

# the following callback functions are not needed as the activation is driven by the motor neurons
if False:  
  def get_specific_states_call_frequency(mu_no):
    stimulation_frequency = motor_units[mu_no % len(motor_units)]["stimulation_frequency"]
    return stimulation_frequency*1e-3

  def get_specific_states_frequency_jitter(mu_no):
    #return 0
    return motor_units[mu_no % len(motor_units)]["jitter"]

  def get_specific_states_call_enable_begin(mu_no):
    #return 1000  # start directly
    #return 0  # start directly
    return motor_units[mu_no % len(motor_units)]["activation_start_time"]*1e3

# callback function for artifical stress values in precontraction computation
def set_gamma_values(n_dofs_global, n_nodes_global_per_coordinate_direction, time_step_no, current_time, values, global_natural_dofs, custom_argument):
    # n_dofs_global:       (int) global number of dofs in the mesh where to set the values
    # n_nodes_global_per_coordinate_direction (list of ints)   [mx, my, mz] number of global nodes in each coordinate direction. 
    #                       For composite meshes, the values are only for the first submesh, for other meshes sum(...) equals n_dofs_global
    # time_step_no:        (int)   current time step number
    # current_time:        (float) the current simulation time
    # values:              (list of floats) all current local values of the field variable, if there are multiple components, they are stored in struct-of-array memory layout 
    #                       i.e. [point0_component0, point0_component1, ... point0_componentN, point1_component0, point1_component1, ...]
    #                       After the call, these values will be assigned to the field variable.
    # global_natural_dofs  (list of ints) for every local dof no. the dof no. in global natural ordering
    # additional_argument: The value of the option "additionalArgument", can be any Python object.
    
    # set all values to 1
    for i in range(len(values)):
      values[i] = constant_gamma

# callback function for artifical lambda values in precontraction computation
def set_lambda_values(n_dofs_global, n_nodes_global_per_coordinate_direction, time_step_no, current_time, values, global_natural_dofs, custom_argument):
    # set all values to 1
    for i in range(len(values)):
      values[i] = 1.0
