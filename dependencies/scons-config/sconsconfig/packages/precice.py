import sys, os, multiprocessing
from .Package import Package
import subprocess

class precice(Package):

  def __init__(self, **kwargs):
    defaults = {
        'download_url': 'https://github.com/precice/precice/archive/v2.0.0.zip',
    }
    defaults.update(kwargs)
    super(precice, self).__init__(**defaults)
    self.ext = '.cpp'
    #self.sub_dirs = [
    #    ('include/mysql', 'lib'),
    #    ('include/mysql', 'lib64'),
    #]
    #self.headers = ['mysql.h']
    #self.libs = ['adios2']
    #self.extra_libs = [['boost_atomic', 'boost_chrono', 'boost_date_time', 'boost_filesystem', 'boost_log', 'boost_log_setup', 'boost_prg_exec_monitor', 'boost_program_options', 'boost_regex', 'boost_system', 'boost_test_exec_monitor', 'boost_thread', 'boost_timer', 'boost_unit_test_framework']]
    self.set_rpath = True
    self.check_text = r'''
    #include <iostream>
    #include <cstdlib>
    #include <fstream>
    #include <precice/SolverInterface.hpp>
    
    int main()
    {
      std::ofstream file("install_precice-config.xml");
      file << R"(<?xml version="1.0"?>
<precice-configuration>
  <solver-interface dimensions="3">
    
    <!-- Data fields that are exchanged between the solvers -->
    <data:scalar name="Data"/>

    <!-- A common mesh that uses these data fields -->
    <mesh name="Mesh">
      <use-data name="Data"/>
    </mesh>
    
    <participant name="Participant1">
      <!-- Makes the named mesh available to the participant. Mesh is provided by the solver directly. -->
      <use-mesh name="Mesh" provide="yes"/>
    </participant>
    <participant name="Participant2">
      <!-- Makes the named mesh available to the participant. Mesh is provided by the solver directly. -->
      <use-mesh name="Mesh" provide="yes"/>
    </participant>
    <m2n:sockets from="Participant1" to="Participant2" network="lo" />

    <coupling-scheme:serial-explicit>
      <participants first="Participant1" second="Participant2"/>
      <time-window-size value="0.01"/>
      <max-time value="0.05"/>
      <exchange data="Data" mesh="Mesh" from="Participant1" to="Participant2"/>
    </coupling-scheme:serial-explicit>
  </solver-interface>
</precice-configuration>
)";
      file.close();
    
      precice::SolverInterface solverInterface("Participant1","install_precice-config.xml",0,1);
      //solverInterface.initialize();
      //solverInterface.finalize();
      
      int ret = system("rm -f install_precice-config.xml");
    
      return EXIT_SUCCESS;
    }
'''

    self.number_output_lines = 15536
      
    self.libs = ['precice']
    self.extra_libs = [[], ['boost_filesystem' 'boost_log' 'boost_log_setup' 'boost_program_options' 'boost_system' 'boost_thread' 'boost_unit_test_framework'],
    									 ['boost_atomic', 'boost_chrono', 'boost_date_time', 'boost_filesystem', 'boost_log', 'boost_log_setup', 'boost_prg_exec_monitor', 'boost_program_options', 'boost_regex', 'boost_system', 'boost_test_exec_monitor', 'boost_thread', 'boost_timer', 'boost_unit_test_framework']]
    self.headers = ["precice/SolverInterface.hpp"]

  def check(self, ctx):
    env = ctx.env
    ctx.Message('Checking for preCICE ...       ')
    self.check_options(env)

    self.set_build_handler([
      'mkdir -p ${PREFIX}/include',
      'cd ${SOURCE_DIR} && wget http://bitbucket.org/eigen/eigen/get/3.3.7.tar.bz2 && tar xf 3.3.7.tar.bz2 && ln -s ${SOURCE_DIR}/eigen-eigen-323c052e1731/Eigen ${PREFIX}/include/Eigen',   # eigen
#      'cd ${SOURCE_DIR} && wget https://dl.bintray.com/boostorg/release/1.65.1/source/boost_1_65_1.tar.gz && tar xf boost_1_65_1.tar.gz && \
#      cd boost_1_65_1 && ./bootstrap.sh --with-libraries=log,thread,system,filesystem,program_options,test --prefix=${PREFIX} && \
#      ./b2 install',
      'cd ${SOURCE_DIR} && mkdir -p build && cd build && '+ctx.env["cmake"]+' -DCMAKE_INSTALL_PREFIX=${PREFIX} \
      -DCMAKE_BUILD_TYPE=RELEASE -DPYTHON_EXECUTABLE=${DEPENDENCIES_DIR}/python/install/bin/python3 \
      -DPETSC_DIR=${PETSC_DIR} -DPETSC_EXECUTABLE_RUNS=TRUE \
      -DLIBXML2_INCLUDE_DIR=${LIBXML2_DIR}/include -DLIBXML2_LIBRARY=${LIBXML2_DIR}/lib/libxml2.so -DPRECICE_ENABLE_FORTRAN=OFF -DMPI_CXX_COMPILER='+ctx.env["mpiCC"]+' .. && \
      cd ${SOURCE_DIR}/build && make all install'
    ])
    
    res = super(precice, self).check(ctx)

    self.check_required(res[0], ctx)
    ctx.Result(res[0])
    return res[0]
