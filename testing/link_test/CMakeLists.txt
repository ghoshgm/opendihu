
# Test for MPI installation.
add_executable(test_mpi test_mpi.cpp)
target_link_libraries(test_mpi ${MPI_CXX_LIBRARIES})
add_test(NAME test_mpi COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_NUMPROC_MAX} ${CMAKE_CURRENT_BINARY_DIR}/test_mpi)

# Test for zlib installation.
add_executable(test_zlib test_zlib.cpp)
target_link_libraries(test_zlib ${ZLIB_LIBRARIES})
add_test(NAME test_zlib COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_zlib)

# Test for Python3 installation.
add_executable(test_python test_python.cpp)
target_link_libraries(test_python ${Python3_LIBRARIES})
add_test(NAME test_python COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_python)

# Test for Eigen3 installation.
add_executable(test_eigen test_eigen.cpp)
target_link_libraries(test_eigen Eigen3::Eigen)
add_test(NAME test_eigen COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_eigen)

# Test for Boost installation.
add_executable(test_boost test_boost.cpp)
target_link_libraries(test_boost ${Boost_LIBRARIES})
add_test(NAME test_boost COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_boost)

# Test for Boost unit_test_framework installation.
add_executable(test_boost_unit_test_framework test_boost_unit_test_framework.cpp)
target_link_libraries(test_boost_unit_test_framework ${Boost_LIBRARIES})
add_test(NAME test_boost_unit_test_framework COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_boost_unit_test_framework)

# Test for Boost installation.
add_executable(test_xml2 test_xml2.cpp)
target_link_libraries(test_xml2 ${LIBXML2_LIBRARIES})
add_test(NAME test_xml2 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_xml2)

# Test for easyloggingpp installation.
add_executable(test_easyloggingpp test_easyloggingpp.cpp)
target_link_libraries(test_easyloggingpp ${EASYLOGGINGPP_LIBRARIES})
add_test(NAME test_easyloggingpp COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_easyloggingpp)

# Test for base64 installtion.
add_executable(test_base64 test_base64.cpp)
target_link_libraries(test_base64 ${BASE64_LIBRARIES})
add_test(NAME test_base64 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_base64)

# Test for Vc installation.
add_executable(test_vc test_vc.cpp)
target_link_libraries(test_vc ${Vc_LIBRARIES})
add_test(NAME test_vc COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_vc)

# Test for PETSc installation.
add_executable(test_petsc test_petsc.cpp)
target_link_libraries(test_petsc ${PETSC_LINK_LIBRARIES} ${MPI_CXX_LIBRARIES})
add_test(NAME test_petsc COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_NUMPROC_MAX} ${CMAKE_CURRENT_BINARY_DIR}/test_petsc)

# Test for ADIOS2 installation.
add_executable(test_adios2 test_adios2.cpp)
target_link_libraries(test_adios2 ${adios2_LIBRARIES} ${MPI_CXX_LIBRARIES})
add_test(NAME test_adios2 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_NUMPROC_MAX} ${CMAKE_CURRENT_BINARY_DIR}/test_adios2)

# Test for xbraid installation.
add_executable(test_xbraid test_xbraid.cpp)
target_link_libraries(test_xbraid ${XBRAID_LIBRARIES} ${MPI_CXX_LIBRARIES})
add_test(NAME test_xbraid COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_NUMPROC_MAX} ${CMAKE_CURRENT_BINARY_DIR}/test_xbraid)

# Test for SEMT installation.
add_executable(test_semt test_semt.cpp)
target_link_libraries(test_semt ${SEMT_LIBRARIES})
add_test(NAME test_semt COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_semt)

# Test for googletest installation.
add_executable(test_gtest test_gtest.cpp)
target_link_libraries(test_gtest ${GTEST_BOTH_LIBRARIES})
add_test(NAME test_gtest COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_gtest)

# Test for preCICE installation.
add_executable(test_precice test_precice.cpp)
target_link_libraries(test_precice precice::precice)
add_test(NAME test_precice COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_precice)

# Test for LAPACK installation.
add_executable(test_lapacke test_lapacke.cpp)
target_link_libraries(test_lapacke ${OpenBLAS_LIBRARIES})
add_test(NAME test_lapacke COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_lapacke)

# Test for BLAS installation.
add_executable(test_blas test_blas.cpp)
target_link_libraries(test_blas ${OpenBLAS_LIBRARIES})
add_test(NAME test_blas COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_blas)