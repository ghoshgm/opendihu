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
rank_no = int(sys.argv[-2])
n_ranks = int(sys.argv[-1])

variables.n_subdomains = variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z

if variables.n_subdomains != n_ranks:
  print("\n\n\033[0;31mError! Number of ranks {} does not match given partitioning {} x {} x {} = {}.\033[0m\n\n".format(n_ranks, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z))
  quit()

variables.relative_factors_muscle1_file = "compartments_relative_factors.{}.{}_mus_partitioning_{}x{}x{}".\
  format(variables.meshes,len(variables.motor_units),variables.n_subdomains_x,variables.n_subdomains_y,variables.n_subdomains_z)

include_global_node_positions = False
if not os.path.exists(variables.relative_factors_muscle1_file) and rank_no == 0:
  include_global_node_positions = True

variables.relative_factors_muscle2_file = "compartments_relative_factors.{}.{}_mus_partitioning_{}x{}x{}".\
  format(variables.meshes,len(variables.motor_units),variables.n_subdomains_x,variables.n_subdomains_y,variables.n_subdomains_z)

include_global_node_positions = False
if not os.path.exists(variables.relative_factors_muscle2_file) and rank_no == 0:
  include_global_node_positions = True


def n_fibers_in_subdomain_x(_):
    return variables.n_fibers_x
def n_fibers_in_subdomain_y(_):
    return variables.n_fibers_y
def get_fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y):
    return fiber_in_subdomain_coordinate_x + variables.n_fibers_x * fiber_in_subdomain_coordinate_y

variables.n_subdomains_xy = variables.n_subdomains_x * variables.n_subdomains_y
variables.n_fibers_total = variables.n_fibers_x * variables.n_fibers_y

#############################
# create the partitioning using the script in create_partitioned_meshes_for_settings.py
# 3D meshes are created internally
# however, we have to create the fibers on our own

# fiber directions
fiber_meshes = {}
fiber_mesh_names = []

### muscle 1
for j in range(variables.n_fibers_y):
  for i in range(variables.n_fibers_x):
    fiber_no = j*variables.n_fibers_x + i
    
    # determine start position of fiber in (x,y)-plane
    x = 0 + i / (variables.n_fibers_x - 1) * variables.muscle1_extent[0]
    y = 0 + j / (variables.n_fibers_y - 1) * variables.muscle1_extent[1]

    # loop over points of a sifiber_file)ngle fiber
    node_positions = []
    for k in range(variables.n_points_whole_fiber):
      x_pos = x
      y_pos = y
      z_pos = 0.0 + k / (variables.n_points_whole_fiber - 1) * variables.muscle1_extent[2]
      node_positions.append([x_pos,y_pos,z_pos])
    
    mesh_name = "muscle1_fiber{}".format(fiber_no)
    fiber_mesh_names.append(mesh_name)
    
    fiber_meshes[mesh_name] = {
      "nodePositions": node_positions,
      "nElements": [variables.n_points_whole_fiber - 1],
      "inputMeshIsGlobal": True,
      "nRanks": [n_ranks],
    }

### muscle 2: same discretization but can have different extent
for j in range(variables.n_fibers_y):
  for i in range(variables.n_fibers_x):
    fiber_no = j*variables.n_fibers_x + i
    
    # determine start position of fiber in (x,y)-plane
    x = 0 + i / (variables.n_fibers_x - 1) * variables.muscle2_extent[0]
    y = 0 + j / (variables.n_fibers_y - 1) * variables.muscle2_extent[1]

    # loop over points of a single fiber
    node_positions = []
    for k in range(variables.n_points_whole_fiber):
      x_pos = x
      y_pos = y
      z_pos = (variables.muscle1_extent[2] + variables.tendon_length) + k / (variables.n_points_whole_fiber - 1) * variables.muscle2_extent[2]
      node_positions.append([x_pos,y_pos,z_pos])
    
    mesh_name = "muscle2_fiber{}".format(fiber_no)
    fiber_mesh_names.append(mesh_name)
    
    fiber_meshes[mesh_name] = {
      "nodePositions": node_positions,
      "nElements": [variables.n_points_whole_fiber - 1],
      "inputMeshIsGlobal": True,
      "nRanks": [n_ranks],
    }



#### set output writer
variables.output_writer_fibers_muscle1 = []
variables.output_writer_fibers_muscle2 = []
variables.output_writer_emg_muscle1 = []
variables.output_writer_emg_muscle2 = []

subfolder = ""
if variables.paraview_output:
  if variables.adios_output:
    subfolder = "paraview/"
  variables.output_writer_fibers_muscle1.append({"format": "Paraview", "outputInterval": int(1./variables.dt_splitting_0D1D*variables.output_timestep_fibers), "filename": "out/" + subfolder + variables.scenario_name + "/muscle1_fibers", "binary": False, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"})
  variables.output_writer_fibers_muscle2.append({"format": "Paraview", "outputInterval": int(1./variables.dt_splitting_0D1D*variables.output_timestep_fibers), "filename": "out/" + subfolder + variables.scenario_name + "/muscle2_fibers", "binary": False, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"})
  variables.output_writer_emg_muscle1.append({"format": "Paraview", "outputInterval": int(1./variables.dt_bidomain*variables.output_timestep_emg), "filename": "out/" + subfolder + variables.scenario_name + "/muscle1_emg", "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"})
  variables.output_writer_emg_muscle2.append({"format": "Paraview", "outputInterval": int(1./variables.dt_bidomain*variables.output_timestep_emg), "filename": "out/" + subfolder + variables.scenario_name + "/muscle2_emg", "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"})
  
if variables.adios_output:
  if variables.paraview_output:
    subfolder = "adios/"
  variables.output_writer_fibers_muscle1.append({"format": "MegaMol", "outputInterval": int(1./variables.dt_splitting_0D1D*variables.output_timestep_fibers), "filename": "out/" + subfolder + variables.scenario_name + "/muscle1_fibers", "useFrontBackBuffer": False, "combineNInstances": 1, "fileNumbering": "incremental"})
  variables.output_writer_fibers_muscle2.append({"format": "MegaMol", "outputInterval": int(1./variables.dt_splitting_0D1D*variables.output_timestep_fibers), "filename": "out/" + subfolder + variables.scenario_name + "/muscle2_fibers", "useFrontBackBuffer": False, "combineNInstances": 1, "fileNumbering": "incremental"})
  variables.output_writer_emg_muscle1.append({"format": "MegaMol", "outputInterval": int(1./variables.dt_bidomain*variables.output_timestep_emg), "filename": "out/" + subfolder + variables.scenario_name + "/muscle1_emg", "useFrontBackBuffer": False, "combineNInstances": 1, "fileNumbering": "incremental"})
  variables.output_writer_emg_muscle2.append({"format": "MegaMol", "outputInterval": int(1./variables.dt_bidomain*variables.output_timestep_emg), "filename": "out/" + subfolder + variables.scenario_name + "/muscle2_emg", "useFrontBackBuffer": False, "combineNInstances": 1, "fileNumbering": "incremental"})
  
if variables.python_output:
  if variables.adios_output:
    subfolder = "python/"
  variables.output_writer_fibers_muscle1.append({"format": "PythonFile", "outputInterval": int(1./variables.dt_splitting_0D1D*variables.output_timestep_fibers), "filename": "out/" + subfolder + variables.scenario_name + "/muscle1_fibers", "binary": True, "fileNumbering": "incremental"})
  variables.output_writer_fibers_muscle2.append({"format": "PythonFile", "outputInterval": int(1./variables.dt_splitting_0D1D*variables.output_timestep_fibers), "filename": "out/" + subfolder + variables.scenario_name + "/muscle2_fibers", "binary": True, "fileNumbering": "incremental"})
  variables.output_writer_emg_muscle1.append({"format": "PythonFile", "outputInterval": int(1./variables.dt_bidomain*variables.output_timestep_emg), "filename": "out/" + subfolder + variables.scenario_name + "/muscle1_emg", "binary": True, "fileNumbering": "incremental"})
  variables.output_writer_emg_muscle2.append({"format": "PythonFile", "outputInterval": int(1./variables.dt_bidomain*variables.output_timestep_emg), "filename": "out/" + subfolder + variables.scenario_name + "/muscle2_emg", "binary": True, "fileNumbering": "incremental"})
  
if variables.exfile_output:
  if variables.adios_output:
    subfolder = "exfile/"
  variables.output_writer_fibers_muscle1.append({"format": "Exfile", "outputInterval": int(1./variables.dt_splitting_0D1D*variables.output_timestep_fibers), "filename": "out/" + subfolder + variables.scenario_name + "/muscle1_fibers", "fileNumbering": "incremental"})
  variables.output_writer_fibers_muscle2.append({"format": "Exfile", "outputInterval": int(1./variables.dt_splitting_0D1D*variables.output_timestep_fibers), "filename": "out/" + subfolder + variables.scenario_name + "/muscle2_fibers", "fileNumbering": "incremental"})
  variables.output_writer_emg_muscle1.append({"format": "Exfile", "outputInterval": int(1./variables.dt_bidomain*variables.output_timestep_emg), "filename": "out/" + subfolder + variables.scenario_name + "/muscle1_emg", "fileNumbering": "incremental"})
  variables.output_writer_emg_muscle2.append({"format": "Exfile", "outputInterval": int(1./variables.dt_bidomain*variables.output_timestep_emg), "filename": "out/" + subfolder + variables.scenario_name + "/muscle2_emg", "fileNumbering": "incremental"})
  


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
    ("connectorSlot","vm"):     "wal_environment/vS",                 # expose state 0 = Vm to the operator splitting
    ("connectorSlot", "stress"):"razumova/stress",                    # expose algebraic 12 = γ to the operator splitting
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
  variables.muscle1_mappings = {
    ("parameter", 0):           "membrane/i_Stim",          # parameter 0 is I_stim
    ("parameter", 1):           "Razumova/l_hs",            # parameter 1 is fiber stretch λ
    ("connectorSlot", "m1vm"):  "membrane/V",               # expose Vm to the operator splitting
    ("connectorSlot", "m1gout"):"Razumova/activestress",
    ("connectorSlot", "m1alp"): "Razumova/activation",      # expose activation .
    ("connectorSlot", "m1lda"): "Razumova/l_hs",            # fiber stretch λ
  }
  variables.muscle2_mappings = {
    ("parameter", 0):           "membrane/i_Stim",          # parameter 0 is I_stim
    ("parameter", 1):           "Razumova/l_hs",            # parameter 1 is fiber stretch λ
    ("connectorSlot", "m2vm"):  "membrane/V",               # expose Vm to the operator splitting
    ("connectorSlot", "m2gout"):"Razumova/activestress",
    ("connectorSlot", "m2alp"): "Razumova/activation",      # expose activation .
    ("connectorSlot", "m2lda"): "Razumova/l_hs",            # fiber stretch λ
  }
  variables.parameters_initial_values = [0, 1]
  variables.nodal_stimulation_current = 40.                           # not used
  variables.vm_value_stimulated = 20.                                 # to which value of Vm the stimulated node should be set (option "valueForStimulatedPoint" of FastMonodomainSolver)

else:
  print("\033[0;31mCellML file {} has no mappings implemented in helper.py\033[0m".format(variables.cellml_file))
  quit()





# load MU distribution and firing times
variables.fiber_distribution = np.genfromtxt(variables.fiber_distribution_file, delimiter=" ", dtype=int)
variables.firing_times = np.genfromtxt(variables.firing_times_file)

# MU indices in file are 1-based.
needed_motorneurons = np.max(variables.fiber_distribution[:variables.n_fibers_total])
# we cannot reindex the MUs so unused MUs are not computed because things like cortical input depend on the MU numbering
if variables.n_motoneurons < needed_motorneurons:
    raise Exception(f"Not enough motoneurons: {variables.n_fibers_total} fibers with MUs in [0, {needed_motorneurons}], but n_motoneurons={variables.n_motoneurons}")


# ---------------------
# callback functions
def get_motor_unit_no(compartment_no):
  return compartment_no

def compartment_gets_stimulated(compartment_no, frequency, current_time):
  """
  determine if compartment compartment_no gets stimulated at simulation time current_time
  """

  # determine motor unit
  alpha = 1.0   # 0.8
  mu_no = int(get_motor_unit_no(compartment_no)*alpha)
  
  # determine if fiber fires now
  index = int(np.round(current_time * frequency))
  n_firing_times = np.size(variables.firing_times,0)
  
  #if variables.firing_times[index % n_firing_times, mu_no] == 1:
    #print("{}: fiber {} is mu {}, t = {}, row: {}, stimulated: {} {}".format(rank_no, fiber_no, mu_no, current_time, (index % n_firing_times), variables.firing_times[index % n_firing_times, mu_no], "true" if variables.firing_times[index % n_firing_times, mu_no] == 1 else "false"))
  
  return variables.firing_times[index % n_firing_times, mu_no] == 1


def set_dummy_stress(n_dofs_global, time_step_no, current_time, values, global_natural_dofs, compartment_no):
  pass

# for debugging output show when the first 20 fibers will fire
if rank_no == 0 and not variables.disable_firing_output:
  print("\nDebugging output about compartment firing: Taking input from file \"{}\"".format(variables.firing_times_file))
  import timeit
  t_start = timeit.default_timer()
  
  first_stimulation_info = []
  
  n_firing_times = np.size(variables.firing_times,0)
  for compartment_no_index in range(variables.n_compartments):
    if compartment_no_index % 100 == 0:
      t_algebraic = timeit.default_timer()
      if t_algebraic - t_start > 100:
        print("Note: break after {}/{} compartments ({:.0f}%) because it already took {:.3f}s".format(compartment_no_index,variables.n_compartments,100.0*compartment_no_index/(variables.n_compartments-1.),t_algebraic - t_start))
        break
    
    first_stimulation = None
    for current_time in np.linspace(0,1./variables.stimulation_frequency*n_firing_times,n_firing_times):
      if compartment_gets_stimulated(compartment_no_index, variables.stimulation_frequency, current_time):
        first_stimulation = current_time
        break
    mu_no = get_motor_unit_no(compartment_no_index)
    first_stimulation_info.append([compartment_no_index,mu_no,first_stimulation])
  
  first_stimulation_info.sort(key=lambda x: 1e6+1e-6*x[1]+1e-12*x[0] if x[2] is None else x[2]+1e-6*x[1]+1e-12*x[0])
  
  print("First stimulation times")
  print("    Time  MU compartments")
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


####################################
# set Dirichlet BC for the flow problem

[nx, ny, nz] = [elem + 1 for elem in variables.n_elements_muscle1]
[mx, my, mz] = [elem // 2 for elem in variables.n_elements_muscle1] # quadratic elements consist of 2 linear elements along each axis

variables.n_points_global = nx * ny * nz
variables.n_elements_global = mx * my * mz
 
# set Dirichlet BC values for bottom nodes to 0 and for top nodes to 1
variables.potential_flow_dirichlet_bc = {}
for i in range(nx*ny):
  variables.potential_flow_dirichlet_bc[i] = 0.0
  variables.potential_flow_dirichlet_bc[(nz-1)*nx*ny + i] = 1.0
  
# set boundary conditions for the elasticity
# Note, we have a composite mesh, consisting of 3Dmesh_elasticity_quadratic and 3DFatMesh_elasticity_quadratic and this composite mesh has a numbering that goes over all dofs.
# The following works because we index the first sub mesh and the first mesh of a composite mesh always has all own dofs with their normal no.s.
# The 2nd mesh has the shared dofs to the first mesh removed in the internal numbering, i.e. internally, they are not counted twice. 
# However, here in the settings, the numbering is concatenated from both meshes, i.e., first all nodes of mesh 0, then all nodes of mesh 1, etc.

# print all mesh information
if False:
  for mesh_name in ["3Dmesh_elasticity_quadratic","3DFatMesh_elasticity_quadratic"]:
    print("")
    print("{}: mesh {}".format(rank_no, mesh_name))
    print("    nPointsLocal: {}".format(variables.meshes[mesh_name]["nPointsLocal"]))
    print("    nPointsGlobal: {}".format(variables.meshes[mesh_name]["nPointsGlobal"]))
    print("    nElementsGlobal: {}".format(variables.meshes[mesh_name]["nElementsGlobal"]))
    print("    inputMeshIsGlobal: {}".format(variables.meshes[mesh_name]["inputMeshIsGlobal"]))
    print("    nElements: {}".format(variables.meshes[mesh_name]["nElements"]))
    print("    nRanks: {}".format(variables.meshes[mesh_name]["nRanks"]))
    print("    len(nodePositions): {}".format(len(variables.meshes[mesh_name]["nodePositions"][1])))


# set Dirichlet BC at top nodes for linear elasticity problem, fix muscle at top

# parameters for the main simulation
# ---------------------------------------------
# https://uni-stuttgart.conceptboard.com/board/i45e-9bz9-qzb0-ppek-n9s2
#          ----> z
#          +----------+-.  tendon  .-+----------+
#  u_z=0 & | muscle 1 |~:~~~~~~~~~~:~| muscle 2 | u_z=0 &
# one      +----------+-'    ^     '-+----------+  one
# edge                       |                     edge
# u=0                    u_x=u_y=0                 u=0
#
#### muscle 1
variables.muscle1_elasticity_dirichlet_bc = {}
# muscle mesh
for j in range(ny):
    for i in range(nx):
      variables.muscle1_elasticity_dirichlet_bc[0*nx*ny + j*nx + i] = [None,None,0.0, None,None,None] # displacement ux uy uz, velocity vx vy vz

# fix edge, note: the multidomain simulation does not work without this (linear solver finds no solution)
for i in range(nx):
    variables.muscle1_elasticity_dirichlet_bc[0*nx*ny + 0*nx + i] = [0.0,0.0,0.0, None,None,None]
    
# fix corner completely
variables.muscle1_elasticity_dirichlet_bc[0*nx*ny + 0] = [0.0,0.0,0.0, None,None,None]

# guide right end of muscle along z axis
# muscle mesh
for j in range(ny):
    for i in range(nx):
      variables.muscle1_elasticity_dirichlet_bc[(nz-1)*nx*ny + j*nx + i] = [0.0,0.0,None, None,None,None]

# initial Neumann BC at bottom nodes, traction along z axis
# will be set by tendon
variables.muscle1_elasticity_neumann_bc = [{"element": (mz-1)*mx*my + j*mx + i, "constantVector": (0,0,0), "face": "2+"} for j in range(my) for i in range(mx)]


#### muscle 2
variables.muscle2_elasticity_dirichlet_bc = {}
# muscle mesh
for j in range(ny):
    for i in range(nx):
      variables.muscle2_elasticity_dirichlet_bc[(nz-1)*nx*ny + j*nx + i] = [None,None,0.0, None,None,None]

# fix edge, note: the multidomain simulation does not work without this (linear solver finds no solution)
for i in range(nx):
    variables.muscle2_elasticity_dirichlet_bc[(nz-1)*nx*ny + 0*nx + i] = [0.0,0.0,0.0, None,None,None]
    
# fix corner completely
variables.muscle2_elasticity_dirichlet_bc[(nz-1)*nx*ny + 0] = [0.0,0.0,0.0, None,None,None]

# guide left end of muscle along z axis
# muscle mesh
for j in range(ny):
    for i in range(nx):
      variables.muscle2_elasticity_dirichlet_bc[0*nx*ny + j*nx + i] = [0.0,0.0,None, None,None,None]

# initial Neumann BC at bottom nodes, traction along z axis
# will be set by tendon
variables.muscle2_elasticity_neumann_bc = [{"element": 0*mx*my + j*mx + i, "constantVector": (0,0,0), "face": "2-"} for j in range(my) for i in range(mx)]



# callback for dirichlet bc
# Function to update dirichelt boundary conditions over time, t.
# This function returns "neumann_bc". Only those entries can be updated that were also initially set.
def muscle1_update_neumann_boundary_conditions_helper(t):
  return variables.muscle1_update_neumann_boundary_conditions(t, [mx,my,mz])
def muscle2_update_neumann_boundary_conditions_helper(t):
  return variables.muscle2_update_neumann_boundary_conditions(t, [mx,my,mz])

####################################
# compute relative factors fr for compartments
def compute_compartment_relative_factors(mesh_node_positions, n_mesh_points_xy, n_mesh_points_z, fiber_data, motor_units):
  """
  Compute the relative factors, f_r, that are needed in the multidomain formulation as a weighting for compartments.
  Result is relative_factors[motor_unit_no][node_no] for the 3D mesh.
  :param mesh_node_positions:  list of (x,y,z) values, global node positions of the 3D mesh
  :param fiber_data: list of fibers, each fiber is a list of points, i.e. point = fiber_data[xy_index][z_index]
  :param motor_units: a list of dicts, settings for the motor units, [{"fiber_no": 0, "standard_deviation": 0.5, "maximum": 1}]
  """
  
  # list of fibers, fiber = list of points, point = list with 3 coordinate entries
  n_compartments = len(motor_units)
  n_points_fiber = len(fiber_data[0])

  # create relative factors for compartments
  #if rank_no == 0:
  #  print("determine relative factors for {} motor units:\n{}".format(n_compartments, motor_units))

  # determine approximate diameter of muscle at every point is z direction
  diameters = []
    
  # loop over points in z direction
  for z_index_mesh in range(n_mesh_points_z):
    
    z_index_fiber = int(z_index_mesh / (float)(n_mesh_points_z) * n_points_fiber)
    
    # get point on first and last fiber
    point0 = np.array(fiber_data[0][z_index_fiber])
    point4 = np.array(fiber_data[(variables.n_fibers_x-1)//2][z_index_fiber])
    point1 = np.array(fiber_data[variables.n_fibers_x-1][z_index_fiber])
    point2 = np.array(fiber_data[-variables.n_fibers_x][z_index_fiber])
    point5 = np.array(fiber_data[(-variables.n_fibers_x)//2][z_index_fiber])
    point3 = np.array(fiber_data[-1][z_index_fiber])
    
    # their distance is an approximation for the diameter
    distance01 = np.linalg.norm(point0 - point1)
    distance02 = np.linalg.norm(point0 - point2)
    distance03 = np.linalg.norm(point0 - point3)
    distance04 = np.linalg.norm(point0 - point4)
    distance05 = np.linalg.norm(point0 - point5)
    distance12 = np.linalg.norm(point1 - point2)
    distance13 = np.linalg.norm(point1 - point3)
    distance14 = np.linalg.norm(point1 - point4)
    distance15 = np.linalg.norm(point1 - point5)
    distance23 = np.linalg.norm(point2 - point3)
    distance24 = np.linalg.norm(point2 - point4)
    distance25 = np.linalg.norm(point2 - point5)
    distance34 = np.linalg.norm(point3 - point4)
    distance35 = np.linalg.norm(point3 - point5)
    distance45 = np.linalg.norm(point4 - point5)
    distance = max(distance01,distance02,distance03,distance04,distance05,distance12,distance13,distance14,distance15,distance23,distance24,distance25,distance34,distance35,distance45)
    diameters.append(distance)

  #print("diameters: {}".format(diameters))

  # create data structure with 0
  relative_factors = np.zeros((n_compartments, len(mesh_node_positions)))   # each row is one compartment

  # loop over nodes of mesh
  for node_no,node_position in enumerate(mesh_node_positions):
    node_position = np.array(node_position)
    
    z_index_mesh = int((float)(node_no) / n_mesh_points_xy)
    z_index_fiber = int(z_index_mesh / (float)(n_mesh_points_z) * n_points_fiber)
    
    # loop over motor units
    for motor_unit_no,motor_unit in enumerate(motor_units):
      
      # find point on fiber that is closest to current node
      fiber_no = motor_unit["fiber_no"]
      if fiber_no >= len(fiber_data):
        new_fiber_no = fiber_no % len(fiber_data)
        if node_no == 0:
          print("\033[0;31mError with motor unit {} around fiber {}, only {} fibers available, now using fiber {} % {} = {} instead.\033[0m".format(motor_unit_no, fiber_no, len(fiber_data), fiber_no, len(fiber_data), new_fiber_no))
        fiber_no = new_fiber_no
      
      min_distance = None
      search_range = int(1 / (float)(n_mesh_points_z) * n_points_fiber)
      search_range = max(10,search_range)
      z_start = max(0,z_index_fiber - search_range)
      z_end = min(n_points_fiber, z_index_fiber + search_range)
      
      #print("node_position: {}, z_index_fiber: {}, fiber at z index: {}, fiber: {}".format(node_position, z_index_fiber, fiber_data[fiber_no][z_index_fiber], fiber_data[fiber_no][z_start:z_end]))
      #print("search_range: {}".format(search_range))
      
      for k,fiber_point in enumerate(fiber_data[fiber_no][z_start:z_end]):
        d = np.array(fiber_point) - node_position
        distance = np.inner(d,d)
        if min_distance is None or distance < min_distance:
          min_distance = distance
          #print("node_position {}, fiber_point {}, d={}, |d|={}".format(node_position, fiber_point, d, np.sqrt(distance)))
      
      distance = np.sqrt(min_distance)
      
      # compute value as gaussian with given standard_deviation and maximum
      standard_deviation = motor_unit["standard_deviation"]*diameters[z_index_mesh]
      gaussian = scipy.stats.norm(loc = 0., scale = standard_deviation)
      value = gaussian.pdf(distance)*standard_deviation*np.sqrt(2*np.pi)*motor_unit["maximum"]
      relative_factors[motor_unit_no][node_no] += value
      
      #print("motor unit {}, fiber {}, distance {}, value {}".format(motor_unit_no, fiber_no, distance, value))

  return relative_factors

####################################
# load relative factors for motor units

# determine relative factor fields fr(x) for compartments for muscle 1
if not os.path.exists(variables.relative_factors_muscle1_file):

  # the file does not yet exist, create it on rank 0
  if rank_no == 0: 
    
    mesh_node_positions_muscle1 = fiber_meshes["muscle1_fiber{}".format(fiber_no)]["nodePositions"]
    n_points_global_muscle1 = [nx, ny, nz]
    n_mesh_points_xy_muscle1 = n_points_global_muscle1[0]*n_points_global_muscle1[1]
    n_mesh_points_z_muscle1 = n_points_global_muscle1[2]
    
    #print("Computing the relative MU factors, f_r, for {} motor units and {} mesh nodes, {} fibers. This may take a while ...".format(len(variables.motor_units), len(mesh_node_positions_muscle1), len(variables.fibers)))
    variables.relative_factors_muscle1 = compute_compartment_relative_factors(mesh_node_positions_muscle1, n_mesh_points_xy_muscle1, n_mesh_points_z_muscle1, variables.fibers, variables.motor_units)
    if rank_no == 0:
      print("Save relative factors to file \"{}\".".format(variables.relative_factors_muscle1_file))
      with open(variables.relative_factors_muscle1_file, "wb") as f:
        pickle.dump(variables.relative_factors_muscle1, f)
  else:
    # wait until file is created on rank 0
    while not os.path.exists(variables.relative_factors_muscle1_file):
      time.sleep(1)

if os.path.exists(variables.relative_factors_muscle1_file):
  with open(variables.relative_factors_muscle1_file, "rb") as f:
    if rank_no == 0:
      print("Load relative factors, f_r, from file \"{}\"".format(variables.relative_factors_muscle1_file))
    variables.relative_factors_muscle1 = pickle.load(f, encoding='latin1')
else:
  print("\033[0;31mError: Could not load relative factors file \"{}\"\033[0m".format(variables.relative_factors_muscle1_file))
  quit()

# debugging output
if rank_no == 0 and not variables.disable_firing_output:
  for i,factors_list in enumerate(variables.relative_factors_muscle1.tolist()):
    print("MU {}, maximum fr: {}".format(i,max(factors_list)))

####################################
# load relative factors for motor units

# determine relative factor fields fr(x) for compartments for muscle 2
if not os.path.exists(variables.relative_factors_muscle2_file):

  # the file does not yet exist, create it on rank 0
  if rank_no == 0: 
    
    mesh_node_positions_muscle2 = variables.fiber_meshes["muscle2_fiber{}".format(fiber_no)]["nodePositions"]
    n_points_global_muscle2 = [nx, ny, nz]
    n_mesh_points_xy_muscle2 = n_points_global_muscle2[0]*n_points_global_muscle2[1]
    n_mesh_points_z_muscle2 = n_points_global_muscle2[2]
    
    #print("Computing the relative MU factors, f_r, for {} motor units and {} mesh nodes, {} fibers. This may take a while ...".format(len(variables.motor_units), len(mesh_node_positions_muscle2), len(variables.fibers)))
    variables.relative_factors_muscle2 = compute_compartment_relative_factors(mesh_node_positions_muscle2, n_mesh_points_xy_muscle2, n_mesh_points_z_muscle2, variables.fibers, variables.motor_units)
    if rank_no == 0:
      print("Save relative factors to file \"{}\".".format(variables.relative_factors_muscle2_file))
      with open(variables.relative_factors_muscle2_file, "wb") as f:
        pickle.dump(variables.relative_factors_muscle2, f)
  else:
    # wait until file is created on rank 0
    while not os.path.exists(variables.relative_factors_muscle2_file):
      time.sleep(1)

if os.path.exists(variables.relative_factors_muscle2_file):
  with open(variables.relative_factors_muscle2_file, "rb") as f:
    if rank_no == 0:
      print("Load relative factors, f_r, from file \"{}\"".format(variables.relative_factors_muscle2_file))
    variables.relative_factors_muscle2 = pickle.load(f, encoding='latin1')
else:
  print("\033[0;31mError: Could not load relative factors file \"{}\"\033[0m".format(variables.relative_factors_muscle2_file))
  quit()

# debugging output
if rank_no == 0 and not variables.disable_firing_output:
  for i,factors_list in enumerate(variables.relative_factors_muscle2.tolist()):
    print("MU {}, maximum fr: {}".format(i,max(factors_list)))

#######################################
# position sensor organs in the 3D mesh

# determine (random) positions of muscle spindles in elasticity mesh
_muscle_spindle_node_nos = []
for muscle_spindle_no in range(variables.n_muscle_spindles):
  random.seed(1)  
  i = random.randrange(0,nx)
  random.seed(2)  
  j = random.randrange(0,ny)
  random.seed(3)  
  k = random.randrange(0,nz)

  dof_no_global = k*nx*ny + j*nx + i
  _muscle_spindle_node_nos.append(dof_no_global)
muscle1_spindle_node_nos = _muscle_spindle_node_nos
muscle2_spindle_node_nos = _muscle_spindle_node_nos
# the muscle spindle mesh holds muscle spdindels of both muscles
muscle1_spindle_indices = list(range(variables.n_muscle_spindles))
muscle2_spindle_indices = list(range(variables.n_muscle_spindles, 2*variables.n_muscle_spindles))

#######################################
# position Golgi tendon organs in the 3D mesh

# determine (random) positions of Golgi organs in elasticity mesh close to tendons
_golgi_tendon_organ_node_nos = []
for golgi_tendon_organ_no in range(variables.n_golgi_tendon_organs):
  random.seed(4)      
  i = random.randrange(0,nx)
  random.seed(5)    
  j = random.randrange(0,ny)
  # position on left and right tendon
  if golgi_tendon_organ_no % 2 == 0:
    k = int(0.1*nz)
  else:
    k = int(0.9*nz)
 
  dof_no_global = k*nx*ny + j*nx + i
  _golgi_tendon_organ_node_nos.append(dof_no_global)
muscle1_golgi_tendon_organ_node_nos = _golgi_tendon_organ_node_nos
muscle2_golgi_tendon_organ_node_nos = _golgi_tendon_organ_node_nos
# the muscle spindle mesh holds muscle spdindels of both muscles
muscle1_golgi_tendon_organ_indices = list(range(variables.n_golgi_tendon_organs))
muscle2_golgi_tendon_organ_indices = list(range(variables.n_golgi_tendon_organs, 2*variables.n_golgi_tendon_organs))

#######################################
# determine positions of neuromuscular junctions in fiber mesh
# stimulate every fiber at its center
center_index = variables.n_points_whole_fiber // 2
stimulation_node_nos = [center_index-1, center_index, center_index+1]

#######################################
# motoneurons
muscle1_motoneuron_indices = list(range(variables.n_motoneurons))
muscle2_motoneuron_indices = list(range(variables.n_motoneurons, 2*variables.n_motoneurons))

#######################################
# combined interneuron + muscle spindle mesh
# (input for motoneurons)
"""
splindles1 --mean--> mn1 \
                          }--mean--> mn1
golgi      --------> in  X
                          }--mean--> mn2
splindles2 --mean--> mn2 /            ^
                      ^               |
                      |               '- motoneuron input: muscle..._motoneuron_indices
                      '- combinded mesh: in_ms_..._indices
"""
in_ms_mesh_muscle1_motoneuron_indices = list(range(variables.n_motoneurons))
in_ms_mesh_muscle2_motoneuron_indices = list(range(variables.n_motoneurons, 2*variables.n_motoneurons))
in_ms_mesh_interneuron_indices = list(range(2*variables.n_motoneurons, 2*variables.n_motoneurons + variables.n_interneurons))

print(f"Muscle 1 indices in spindle    mesh: {muscle1_spindle_indices[0]:3}...{muscle1_spindle_indices[-1]:3}")
print(f"Muscle 2 indices in spindle    mesh: {muscle2_spindle_indices[0]:3}...{muscle2_spindle_indices[-1]:3}")
print(f"Muscle 1 indices in golgi t.o. mesh: {muscle1_golgi_tendon_organ_indices[0]:3}...{muscle1_golgi_tendon_organ_indices[-1]:3}")
print(f"Muscle 2 indices in golgi t.o. mesh: {muscle1_golgi_tendon_organ_indices[0]:3}...{muscle1_golgi_tendon_organ_indices[-1]:3}")
print(f"Muscle 1 indices in motoneuron mesh: {muscle1_motoneuron_indices[0]:3}...{muscle1_motoneuron_indices[-1]:3}")
print(f"Muscle 2 indices in motoneuron mesh: {muscle2_motoneuron_indices[0]:3}...{muscle2_motoneuron_indices[-1]:3}")
print(f"Muscle 1 indices in sp.+itern. mesh: {in_ms_mesh_muscle1_motoneuron_indices[0]:3}...{in_ms_mesh_muscle1_motoneuron_indices[-1]:3}")
print(f"Muscle 2 indices in sp.+itern. mesh: {in_ms_mesh_muscle2_motoneuron_indices[0]:3}...{in_ms_mesh_muscle2_motoneuron_indices[-1]:3}")
print(f"Intern.  indices in sp.+itern. mesh: {in_ms_mesh_interneuron_indices[0]:3}...{in_ms_mesh_interneuron_indices[-1]:3}")
