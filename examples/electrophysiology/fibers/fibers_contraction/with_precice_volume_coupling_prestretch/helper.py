# Multiple 1D fibers (monodomain) with 3D contraction, biceps geometry
# This is a helper script that sets a lot of the internal variables which are all defined in variables.py

import numpy as np
import scipy
import pickle
import sys,os
import struct
import argparse
import random
import time
sys.path.insert(0, '..')
import variables    # file variables.py
from create_partitioned_meshes_for_settings import *   # file create_partitioned_meshes_for_settings

# parse arguments
rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

variables.n_subdomains = variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z
variables.load_fiber_data = True   # load all local node positions from fiber_file, in order to infer partitioning for fat_layer mesh

# create the partitioning using the script in create_partitioned_meshes_for_settings.py
result = create_partitioned_meshes_for_settings(
    variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, 
    variables.fiber_file, variables.load_fiber_data,
    variables.sampling_stride_x, variables.sampling_stride_y, variables.sampling_stride_z, variables.generate_linear_3d_mesh, variables.generate_quadratic_3d_mesh,
    fiber_set_rank_nos=False, have_fibers=True)
[variables.meshes, variables.own_subdomain_coordinate_x, variables.own_subdomain_coordinate_y, variables.own_subdomain_coordinate_z, variables.n_fibers_x, variables.n_fibers_y, variables.n_points_whole_fiber] = result
  
variables.n_subdomains_xy = variables.n_subdomains_x * variables.n_subdomains_y
variables.n_fibers_total = variables.n_fibers_x * variables.n_fibers_y

# create mappings between meshes
#variables.mappings_between_meshes = {"MeshFiber_{}".format(i) : "3Dmesh" for i in range(variables.n_fibers_total)}
variables.mappings_between_meshes = {"MeshFiber_{}".format(i) : {"name": "3Dmesh", "xiTolerance": 3e-1, "defaultValue": 0} for i in range(variables.n_fibers_total)}

# set output writer    
variables.output_writer_fibers = []
variables.output_writer_elasticity = []
variables.output_writer_emg = []
variables.output_writer_0D_states = []



subfolder = ""
if variables.paraview_output:
  if variables.adios_output:
    subfolder = "paraview/"
  variables.output_writer_emg.append({"format": "Paraview", "outputInterval": int(1./variables.dt_3D*variables.output_timestep_big), "filename": "out/" + subfolder + variables.scenario_name + "/hd_emg", "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"})
  variables.output_writer_elasticity.append({"format": "Paraview", "outputInterval": int(1./variables.dt_3D*variables.output_timestep_big), "filename": "out/" + subfolder + variables.scenario_name + "/elasticity", "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"})
  variables.output_writer_fibers.append({"format": "Paraview", "outputInterval": int(1./variables.dt_splitting*variables.output_timestep_fibers), "filename": "out/" + subfolder + variables.scenario_name + "/fibers", "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"})
  if variables.states_output:
    variables.output_writer_0D_states.append({"format": "Paraview", "outputInterval": 1, "filename": "out/" + subfolder + variables.scenario_name + "/0D_states", "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"})

if variables.adios_output:
  if variables.paraview_output:
    subfolder = "adios/"
  variables.output_writer_emg.append({"format": "MegaMol", "outputInterval": int(1./variables.dt_3D*variables.output_timestep), "filename": "out/" + subfolder + variables.scenario_name + "/hd_emg", "useFrontBackBuffer": False, "combineNInstances": 1, "fileNumbering": "incremental"})
  variables.output_writer_elasticity.append({"format": "MegaMol", "outputInterval": int(1./variables.dt_3D*variables.output_timestep), "filename": "out/" + subfolder + variables.scenario_name + "/elasticity", "useFrontBackBuffer": False, "fileNumbering": "incremental"})
  variables.output_writer_fibers.append({"format": "MegaMol", "outputInterval": int(1./variables.dt_splitting*variables.output_timestep_fibers), "filename": "out/" + subfolder + variables.scenario_name + "/fibers", "combineNInstances": variables.n_subdomains_xy, "useFrontBackBuffer": False, "fileNumbering": "incremental"})
  #variables.output_writer_fibers.append({"format": "MegaMol", "outputInterval": int(1./variables.dt_splitting*variables.output_timestep), "filename": "out/" + variables.scenario_name + "/fibers", "combineNInstances": 1, "useFrontBackBuffer": False, "fileNumbering": "incremental"}

if variables.python_output:
  if variables.adios_output:
    subfolder = "python/"
  variables.output_writer_emg.append({"format": "PythonFile", "outputInterval": int(1./variables.dt_3D*variables.output_timestep), "filename": "out/" + subfolder + variables.scenario_name + "/hd_emg", "binary": True, "fileNumbering": "incremental"})
  variables.output_writer_elasticity.append({"format": "PythonFile", "outputInterval": int(1./variables.dt_3D*variables.output_timestep), "filename": "out/" + subfolder + variables.scenario_name + "/elasticity", "binary": True, "fileNumbering": "incremental"})
  variables.output_writer_fibers.append({"format": "PythonFile", "outputInterval": int(1./variables.dt_splitting*variables.output_timestep_fibers), "filename": "out/" + subfolder + variables.scenario_name + "/fibers", "binary": True, "fileNumbering": "incremental"})

if variables.exfile_output:
  if variables.adios_output:
    subfolder = "exfile/"
  variables.output_writer_emg.append({"format": "Exfile", "outputInterval": int(1./variables.dt_3D*variables.output_timestep), "filename": "out/" + subfolder + variables.scenario_name + "/hd_emg", "fileNumbering": "incremental"})
  variables.output_writer_elasticity.append({"format": "Exfile", "outputInterval": int(1./variables.dt_3D*variables.output_timestep), "filename": "out/" + subfolder + variables.scenario_name + "/elasticity", "fileNumbering": "incremental"})
  variables.output_writer_fibers.append({"format": "Exfile", "outputInterval": int(1./variables.dt_splitting*variables.output_timestep_fibers), "filename": "out/" + subfolder + variables.scenario_name + "/fibers", "fileNumbering": "incremental"})

# set variable mappings for cellml model
if "hodgkin_huxley" in variables.cellml_file and "hodgkin_huxley-razumova" not in variables.cellml_file:
  # parameters: I_stim
  variables.mappings = {
    ("parameter", 0):           ("constant", "membrane/i_Stim"),      # parameter 0 is constant 2 = I_stim
    ("connectorSlot", 0): ("state", "membrane/V"),              # expose state 0 = Vm to the operator splitting
  }
  variables.parameters_initial_values = [0.0]                         # initial value for stimulation current
  variables.nodal_stimulation_current = 40.                           # not used
  variables.vm_value_stimulated = 20.                                 # to which value of Vm the stimulated node should be set (option "valueForStimulatedPoint" of FastMonodomainSolver)

elif "shorten" in variables.cellml_file:
  # parameters: stimulation current I_stim, fiber stretch λ
  variables.mappings = {
    ("parameter", 0):           ("algebraic", "wal_environment/I_HH"), # parameter is algebraic 32
    ("parameter", 1):           ("constant", "razumova/L_x"),             # parameter is constant 65, fiber stretch λ, this indicates how much the fiber has stretched, 1 means no extension
    ("connectorSlot", 0): ("state", "wal_environment/vS"),          # expose state 0 = Vm to the operator splitting
  }
  variables.parameters_initial_values = [0.0, 1.0]                        # stimulation current I_stim, fiber stretch λ
  variables.nodal_stimulation_current = 1200.                             # not used
  variables.vm_value_stimulated = 40.                                 # to which value of Vm the stimulated node should be set (option "valueForStimulatedPoint" of FastMonodomainSolver)
  
elif "slow_TK_2014" in variables.cellml_file:   # this is (3a, "MultiPhysStrain", old tomo mechanics) in OpenCMISS
  # parameters: I_stim, fiber stretch λ
  variables.mappings = {
    ("parameter", 0):           ("constant", "wal_environment/I_HH"), # parameter 0 is constant 54 = I_stim
    ("parameter", 1):           ("constant", "razumova/L_S"),         # parameter 1 is constant 67 = fiber stretch λ
    ("connectorSlot", 0): ("state", "wal_environment/vS"),      # expose state 0 = Vm to the operator splitting
    ("connectorSlot", 1, "stress"): ("algebraic", "razumova/stress"),  # expose algebraic 12 = γ to the operator splitting
  }
  variables.parameters_initial_values = [0.0, 1.0]                    # wal_environment/I_HH = I_stim, razumova/L_S = λ
  variables.nodal_stimulation_current = 40.                           # not used
  variables.vm_value_stimulated = 40.                                 # to which value of Vm the stimulated node should be set (option "valueForStimulatedPoint" of FastMonodomainSolver)
  
elif "Aliev_Panfilov_Razumova_2016_08_22" in variables.cellml_file :   # this is (3, "MultiPhysStrain", numerically more stable) in OpenCMISS, this only computes A1,A2,x1,x2 not the stress
  # parameters: I_stim, fiber stretch λ, fiber contraction velocity \dot{λ}
  variables.mappings = {
    ("parameter", 0):           ("constant", "Aliev_Panfilov/I_HH"),  # parameter 0 is constant 0 = I_stim
    ("parameter", 1):           ("constant", "Razumova/l_hs"),        # parameter 1 is constant 8 = fiber stretch λ
    ("parameter", 2):           ("constant", "Razumova/velo"),        # parameter 2 is constant 9 = fiber contraction velocity \dot{λ}
    ("connectorSlot", 0): ("state", "Aliev_Panfilov/V_m"),      # expose state 0 = Vm to the operator splitting
    ("connectorSlot", 1): ("algebraic", "Razumova/sigma"),   # expose algebraic 0 = γ to the operator splitting
  }
  variables.parameters_initial_values = [0, 1, 0]                     # Aliev_Panfilov/I_HH = I_stim, Razumova/l_hs = λ, Razumova/velo = \dot{λ}
  variables.nodal_stimulation_current = 40.                           # not used
  variables.vm_value_stimulated = 40.                                 # to which value of Vm the stimulated node should be set (option "valueForStimulatedPoint" of FastMonodomainSolver)
  
elif "Aliev_Panfilov_Razumova_Titin" in variables.cellml_file:   # this is (4, "Titin") in OpenCMISS
  # parameters: I_stim, fiber stretch λ, fiber contraction velocity \dot{λ}
  variables.mappings = {
    ("parameter", 0):           ("constant", "Aliev_Panfilov/I_HH"),  # parameter 0 is constant 0 = I_stim
    ("parameter", 1):           ("constant", "Razumova/l_hs"),        # parameter 1 is constant 11 = fiber stretch λ
    ("parameter", 2):           ("constant", "Razumova/rel_velo"),    # parameter 2 is constant 12 = fiber contraction velocity \dot{λ}
    ("connectorSlot", 0): ("state", "Aliev_Panfilov/V_m"),      # expose state 0 = Vm to the operator splitting
    ("connectorSlot", 1): ("algebraic", "Razumova/ActiveStress"),   # expose algebraic 4 = γ to the operator splitting
    ("connectorSlot", 2): ("algebraic", "Razumova/Activation"),     # expose algebraic 5 = α to the operator splitting
  }
  variables.parameters_initial_values = [0, 1, 0]                     # Aliev_Panfilov/I_HH = I_stim, Razumova/l_hs = λ, Razumova/rel_velo = \dot{λ}
  variables.nodal_stimulation_current = 40.                           # not used
  variables.vm_value_stimulated = 40.                                 # to which value of Vm the stimulated node should be set (option "valueForStimulatedPoint" of FastMonodomainSolver)
  
elif "hodgkin_huxley-razumova" in variables.cellml_file:   # this is (4, "Titin") in OpenCMISS
  # parameters: I_stim, fiber stretch λ, fiber contraction velocity \dot{λ}
  variables.mappings = {
    ("parameter", 0):           "membrane/i_Stim",          # parameter 0 is I_stim
    ("parameter", 1):           "Razumova/l_hs",            # parameter 1 is fiber stretch λ
    ("connectorSlot", "vmf"):    "membrane/V",               # expose Vm to the operator splitting
    ("connectorSlot", "stress"):"Razumova/activestress",
    ("connectorSlot", "alpha"): "Razumova/activation",      # expose activation .
  }
  variables.parameters_initial_values = [0, 1]
  variables.nodal_stimulation_current = 40.                           # not used
  variables.vm_value_stimulated = 20.                                 # to which value of Vm the stimulated node should be set (option "valueForStimulatedPoint" of FastMonodomainSolver)

else:
  print("\033[0;31mCellML file {} has no mappings implemented in helper.py\033[0m".format(variables.cellml_file))
  quit()

# load MU distribution and firing times
variables.firing_times = np.genfromtxt(variables.firing_times_file)

# ---------------------
# callback functions
def get_motor_unit_no(fiber_no):
  return int(variables.fiber_distribution[fiber_no % len(variables.fiber_distribution)]-1)

def get_diffusion_prefactor(fiber_no, mu_no):
  return variables.get_conductivity(fiber_no, mu_no) / (variables.get_am(fiber_no, mu_no) * variables.get_cm(fiber_no, mu_no))

def fiber_gets_stimulated(fiber_no, frequency, current_time):
  """
  determine if fiber fiber_no gets stimulated at simulation time current_time
  """

  # determine motor unit
  alpha = 1.0   # 0.8
  mu_no = (int)(get_motor_unit_no(fiber_no)*alpha)
  
  # determine if fiber fires now
  index = int(np.round(current_time * frequency))
  n_firing_times = np.size(variables.firing_times,0)
  
  #if variables.firing_times[index % n_firing_times, mu_no] == 1:
    #print("{}: fiber {} is mu {}, t = {}, row: {}, stimulated: {} {}".format(rank_no, fiber_no, mu_no, current_time, (index % n_firing_times), variables.firing_times[index % n_firing_times, mu_no], "true" if variables.firing_times[index % n_firing_times, mu_no] == 1 else "false"))
  
  return variables.firing_times[index % n_firing_times, mu_no] == 1
  
# callback function that can set states, i.e. prescribed values for stimulation
def set_specific_states(n_nodes_global, time_step_no, current_time, states, fiber_no):

  #print("call set_specific_states at time {}".format(current_time)) 
  
  # determine if fiber gets stimulated at the current time
  is_fiber_gets_stimulated = fiber_gets_stimulated(fiber_no, variables.stimulation_frequency, current_time)

  if is_fiber_gets_stimulated:  
    # determine nodes to stimulate (center node, left and right neighbour)
    innervation_zone_width_n_nodes = variables.innervation_zone_width*100  # 100 nodes per cm
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

# load MU distribution and firing times
variables.fiber_distribution = np.genfromtxt(variables.fiber_distribution_file, delimiter=" ")
variables.firing_times = np.genfromtxt(variables.firing_times_file)

# for debugging output show when the first 20 fibers will fire
if rank_no == 0 and not variables.disable_firing_output:
  print("Debugging output about fiber firing: Taking input from file \"{}\"".format(variables.firing_times_file))
  import timeit
  t_start = timeit.default_timer()
  
  first_stimulation_info = []
  
  n_firing_times = np.size(variables.firing_times,0)
  for fiber_no_index in range(variables.n_fibers_total):
    if fiber_no_index % 100 == 0:
      t_algebraic = timeit.default_timer()
      if t_algebraic - t_start > 100:
        print("Note: break after {}/{} fibers ({:.0f}%) because it already took {:.3f}s".format(fiber_no_index,variables.n_fibers_total,100.0*fiber_no_index/(variables.n_fibers_total-1.),t_algebraic - t_start))
        break
    
    first_stimulation = None
    for current_time in np.linspace(0,1./variables.stimulation_frequency*n_firing_times,n_firing_times):
      if fiber_gets_stimulated(fiber_no_index, variables.stimulation_frequency, current_time):
        first_stimulation = current_time
        break
    mu_no = get_motor_unit_no(fiber_no_index)
    first_stimulation_info.append([fiber_no_index,mu_no,first_stimulation])
  
  first_stimulation_info.sort(key=lambda x: 1e6+1e-6*x[1]+1e-12*x[0] if x[2] is None else x[2]+1e-6*x[1]+1e-12*x[0])
  
  print("First stimulation times")
  print("    Time  MU fibers")
  n_stimulated_mus = 0
  n_not_stimulated_mus = 0
  stimulated_fibers = []
  last_time = 0
  last_mu_no = first_stimulation_info[0][1]
  for stimulation_info in first_stimulation_info:
    mu_no = stimulation_info[1]
    fiber_no = stimulation_info[0]
    if mu_no == last_mu_no:
      stimulated_fibers.append(fiber_no)
    else:
      if last_time is not None:
        if len(stimulated_fibers) > 10:
          print("{:8.2f} {:3} {} (only showing first 10, {} total)".format(last_time,last_mu_no,str(stimulated_fibers[0:10]),len(stimulated_fibers)))
        else:
          print("{:8.2f} {:3} {}".format(last_time,last_mu_no,str(stimulated_fibers)))
        n_stimulated_mus += 1
      else:
        if len(stimulated_fibers) > 10:
          print("  never stimulated: MU {:3}, fibers {} (only showing first 10, {} total)".format(last_mu_no,str(stimulated_fibers[0:10]),len(stimulated_fibers)))
        else:
          print("  never stimulated: MU {:3}, fibers {}".format(last_mu_no,str(stimulated_fibers)))
        n_not_stimulated_mus += 1
      stimulated_fibers = [fiber_no]

    last_time = stimulation_info[2]
    last_mu_no = mu_no
    
  print("stimulated MUs: {}, not stimulated MUs: {}".format(n_stimulated_mus,n_not_stimulated_mus))

  t_end = timeit.default_timer()
  print("duration of assembling this list: {:.3f} s\n".format(t_end-t_start))  
  
# compute partitioning
if rank_no == 0:
  if n_ranks != variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z:
    print("\n\nError! Number of ranks {} does not match given partitioning {} x {} x {} = {}.\n\n".format(n_ranks, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z))
    quit()
  
####################################
# set Dirichlet BC for the flow problem

n_points_3D_mesh_linear_global_x = sum([n_sampled_points_in_subdomain_x(subdomain_coordinate_x) for subdomain_coordinate_x in range(variables.n_subdomains_x)])
n_points_3D_mesh_linear_global_y = sum([n_sampled_points_in_subdomain_y(subdomain_coordinate_y) for subdomain_coordinate_y in range(variables.n_subdomains_y)])
n_points_3D_mesh_linear_global_z = sum([n_sampled_points_in_subdomain_z(subdomain_coordinate_z) for subdomain_coordinate_z in range(variables.n_subdomains_z)])
n_points_3D_mesh_linear_global = n_points_3D_mesh_linear_global_x*n_points_3D_mesh_linear_global_y*n_points_3D_mesh_linear_global_z

n_points_3D_mesh_quadratic_global_x = 2*n_points_3D_mesh_linear_global_x - 1
n_points_3D_mesh_quadratic_global_y = 2*n_points_3D_mesh_linear_global_y - 1
n_points_3D_mesh_quadratic_global_z = 2*n_points_3D_mesh_linear_global_z - 1
 
variables.n_points_global_mesh = n_points_3D_mesh_quadratic_global_x*n_points_3D_mesh_quadratic_global_y*n_points_3D_mesh_quadratic_global_z
 
# set Dirichlet BC values for bottom nodes to 0 and for top nodes to 1
variables.potential_flow_dirichlet_bc = {}
for i in range(n_points_3D_mesh_linear_global_x*n_points_3D_mesh_linear_global_y):
  variables.potential_flow_dirichlet_bc[i] = 0.0
  variables.potential_flow_dirichlet_bc[(n_points_3D_mesh_linear_global_z-1)*n_points_3D_mesh_linear_global_x*n_points_3D_mesh_linear_global_y + i] = 1.0

variables.n_subdomains_xy = variables.n_subdomains_x * variables.n_subdomains_y
variables.n_fibers_total = variables.n_fibers_x * variables.n_fibers_y
  
# set boundary conditions for the elasticity
[mx, my, mz] = variables.meshes["3Dmesh_quadratic"]["nPointsGlobal"]
nx = (mx-1)//2
ny = (my-1)//2
nz = (mz-1)//2

# copy meshes
import copy
variables.meshes["3Dmesh_quadratic_precontraction"] = copy.deepcopy(variables.meshes["3Dmesh_quadratic"])

# set Dirichlet BC at top nodes for linear elasticity problem, fix muscle at top

# parameters for precontraction
# -----------------------------
variables.precontraction_elasticity_dirichlet_bc = {}
# muscle mesh
for j in range(my):
  for i in range(mx):
    variables.precontraction_elasticity_dirichlet_bc[(mz-1)*mx*my + j*mx + i] = [None,None,0.0,None,None,None]

# fix edge
for i in range(mx):
  variables.precontraction_elasticity_dirichlet_bc[(mz-1)*mx*my + 0*mx + i] = [0.0,0.0,0.0,None,None,None]
  
# fix corner completely
variables.precontraction_elasticity_dirichlet_bc[(mz-1)*mx*my + 0] = [0.0,0.0,0.0,None,None,None]

# guide lower end of muscle along z axis
# muscle mesh
for j in range(my):
  for i in range(mx):
    variables.precontraction_elasticity_dirichlet_bc[0*mx*my + j*mx + i] = [0.0,0.0,None,None,None,None]

# Neumann BC at bottom nodes, traction downwards
# muscle mesh
variables.precontraction_elasticity_neumann_bc = [{"element": 0*nx*ny + j*nx + i, "constantVector": variables.precontraction_bottom_traction, "face": "2-"} for j in range(ny) for i in range(nx)]

# parameters for prestretch
# -----------------------------
if hasattr(variables, "prestretch_bottom_traction"):
  variables.prestretch_elasticity_dirichlet_bc = {}
  # muscle mesh
  for j in range(my):
    for i in range(mx):
      variables.prestretch_elasticity_dirichlet_bc[(mz-1)*mx*my + j*mx + i] = [None,None,0.0,None,None,None]

  # fix edge, note: the prestretch simulation does not work without this (linear solver finds no solution)
  for i in range(mx):
    variables.prestretch_elasticity_dirichlet_bc[(mz-1)*mx*my + 0*mx + i] = [0.0,0.0,0.0,None,None,None]
    
  # fix corner completely
  variables.prestretch_elasticity_dirichlet_bc[(mz-1)*mx*my + 0] = [0.0,0.0,0.0,None,None,None]

  # guide lower end of muscle along z axis
  # muscle mesh
  for j in range(my):
    for i in range(mx):
      variables.prestretch_elasticity_dirichlet_bc[0*mx*my + j*mx + i] = [0.0,0.0,None,None,None,None]

  # Neumann BC at bottom nodes, traction downwards
  # muscle mesh
  variables.prestretch_elasticity_neumann_bc = [{"element": 0*nx*ny + j*nx + i, "constantVector": variables.prestretch_bottom_traction, "face": "2-"} for j in range(ny) for i in range(nx)]

# parameters for the main simulation
# ---------------------------------------------
if hasattr(variables, "main_bottom_traction"):
  variables.main_elasticity_dirichlet_bc = {}
  # muscle mesh
  for j in range(my):
    for i in range(mx):
      variables.main_elasticity_dirichlet_bc[(mz-1)*mx*my + j*mx + i] = [None,None,0.0,None,None,None]

# set boundary conditions for the elasticity
[mx, my, mz] = variables.meshes["3Dmesh_quadratic"]["nPointsGlobal"]
nx = (mx-1)//2
ny = (my-1)//2
nz = (mz-1)//2

# set Dirichlet BC at top nodes for elasticity problem, fix muscle at top
variables.elasticity_dirichlet_bc = {}
for j in range(my):
  for i in range(mx):
    variables.elasticity_dirichlet_bc[(mz-1)*mx*my + j*mx + i] = [None,None,0.0,None,None,None]
  
# fix edge
for i in range(mx):
  #variables.elasticity_dirichlet_bc[(mz-1)*mx*my + 0*mx + i] = [0.0,None,0.0,None,None,None]
  variables.elasticity_dirichlet_bc[(mz-1)*mx*my + 0*mx + i] = [0.0,0.0,0.0,None,None,None]
  
# fix corner completely
variables.elasticity_dirichlet_bc[(mz-1)*mx*my + 0] = [0.0,0.0,0.0,None,None,None]


# guide lower end of muscle along z axis
# muscle mesh
for j in range(my):
  for i in range(mx):
    variables.main_elasticity_dirichlet_bc[0*mx*my + j*mx + i] = [0.0,0.0,None,None,None,None]

# also fix z displacements for one point at the center
#if variables.fix_bottom:
#  variables.main_elasticity_dirichlet_bc[0*mx*my + (my//2)*mx + mx//2][2] = 0.0


# Neumann BC at bottom nodes, traction downwards
variables.elasticity_neumann_bc = [{"element": 0*nx*ny + j*nx + i, "constantVector": variables.bottom_traction, "face": "2-"} for j in range(ny) for i in range(nx)]
#variables.elasticity_neumann_bc = []

#print("bottom_traction={}\n elasticity_neumann_bc={}".format(variables.bottom_traction,variables.elasticity_neumann_bc))

#with open("mesh","w") as f:
#  f.write(str(variables.meshes["3Dmesh_quadratic"]))
    
