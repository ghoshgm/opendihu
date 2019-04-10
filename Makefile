.PHONY: clean all

all: debug release

#ifneq ("$(wildcard ./dependencies/python/install/bin/python)","")
#python := ./dependencies/python/install/bin/python
#else
python := python2.7
#endif

debug:
	$(python) dependencies/scons/scons.py BUILD_TYPE=DEBUG

release:
	$(python) dependencies/scons/scons.py BUILD_TYPE=RELEASE

clean:
	rm -rf .sconf_temp
	rm .sconsign.dblite	

purge: clean
	rm -rf core/build_debug
	rm -rf core/build_release

purge_dependencies:
	cd dependencies; rm -rf base64/ bzip2/ cython/ easyloggingpp/ googletest/ lapack/ matplotlib/ numpyc/ petsc/ python/ scipy/ semt/ pythonpackages; cd -

rebuild: purge_dependencies purge clean debug release

# on hazel hen rebuild everying including dependencies
rebuild_hazelhen:
	rm -rf dependencies/easyloggingpp/install dependencies/easyloggingpp/src dependencies/python/install dependencies/python/src  dependencies/base64/install dependencies/base64/src dependencies/numpyc/install dependencies/numpyc/src dependencies/semt/src dependencies/semt/install && rm -rf core/build_release && $(python) dependencies/scons/scons.py BUILD_TYPE=RELEASE; cd dependencies/matplotlib && ../python/install/bin/pip3 install *.whl

doc:
	cd doc/doxygen; doxygen
	
# the following targets are just for convenience and could also be deleted
release_without_tests:
	$(python) dependencies/scons/scons.py BUILD_TYPE=RELEASE no_tests=True

system_testing:
	cd testing/system_testing && ./run.sh

solid_mechanics:
	cd testing/system_testing/tests/solid_mechanics && python ../../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

multiple_fibers_system_testing:
	cd testing/system_testing/tests/multiple_fibers && python ../../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

streamline_tracer:
	cd testing/system_testing/tests/fibers && python ../../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

diffusion:
	cd testing/system_testing/tests/diffusion &&  python ../../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

quadrature:
	cd examples/quadrature/own && python ../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

fibers:
	cd testing/system_testing/tests/fibers &&  python ../../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

hodgkin_huxley:
	cd examples/electrophysiology/hodgkin_huxley && python ../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

shorten:
	cd examples/electrophysiology/shorten && python ../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

cellml:
	cd examples/electrophysiology/cellml && python ../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

pod:
	cd examples/diffusion1d && python ../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

multidomain:
	cd examples/multidomain3d && python ../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

parallel_fiber_estimation:
	cd examples/parallel_fiber_estimation && python ../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

load_balancing:
	cd examples/load_balancing && python ../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

multiple_fibers:
	cd examples/electrophysiology/multiple_fibers && python ../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG
	
multiple_fibers_cubes_partitioning:
	cd examples/electrophysiology/multiple_fibers_cubes_partitioning && python ../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG
	
fibers_emg:
	cd examples/electrophysiology/fibers_emg && python ../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

laplace2d:
	cd examples/laplace/laplace2d && python ../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG
