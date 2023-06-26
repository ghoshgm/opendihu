
# Test for MPI installation.
add_executable(test_mpi test_mpi.cpp)
target_include_directories(test_mpi PUBLIC ${MPI_CXX_INCLUDE_DIRS})
target_link_libraries(test_mpi PUBLIC ${MPI_CXX_LIBRARIES})
add_test(NAME test_mpi COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_NUMPROC_MAX} ${CMAKE_CURRENT_BINARY_DIR}/test_mpi)

# Test for zlib installation.
add_executable(test_zlib test_zlib.cpp)
target_include_directories(test_zlib PUBLIC ${ZLIB_INCLUDE_DIRS})
target_link_libraries(test_zlib PUBLIC ${ZLIB_LIBRARIES})
add_test(NAME test_zlib COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_zlib)

# Test for Python3 installation.
add_executable(test_python test_python.cpp)
target_include_directories(test_python PUBLIC ${Python3_INCLUDE_DIRS})
target_link_libraries(test_python PUBLIC ${Python3_LIBRARIES})
add_test(NAME test_python COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_python)

# Test for Eigen3 installation.
add_executable(test_eigen test_eigen.cpp)
target_link_libraries(test_eigen PUBLIC Eigen3::Eigen)
add_test(NAME test_eigen COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_eigen)

# Test for Boost installation.
add_executable(test_boost test_boost.cpp)
target_link_libraries(test_boost PUBLIC ${Boost_LIBRARIES})
add_test(NAME test_boost COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_boost)

# Test for Boost unit_test_framework installation.
add_executable(test_boost_unit_test_framework test_boost_unit_test_framework.cpp)
target_link_libraries(test_boost_unit_test_framework PUBLIC ${Boost_LIBRARIES})
add_test(NAME test_boost_unit_test_framework COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_boost_unit_test_framework)

# Test for LibXml2 installation.
add_executable(test_xml2 test_xml2.cpp)
target_include_directories(test_xml2 PUBLIC ${LIBXML2_INCLUDE_DIR})
target_link_libraries(test_xml2 PUBLIC ${LIBXML2_LIBRARIES})
add_test(NAME test_xml2 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_xml2)

# Test for easyloggingpp installation.
add_executable(test_easyloggingpp test_easyloggingpp.cpp)
target_include_directories(test_easyloggingpp PUBLIC ${EASYLOGGINGPP_INCLUDE_DIRS})
target_link_libraries(test_easyloggingpp PUBLIC ${EASYLOGGINGPP_LIBRARIES})
add_test(NAME test_easyloggingpp COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_easyloggingpp)

# Test for base64 installtion.
add_executable(test_base64 test_base64.cpp)
target_include_directories(test_base64 PUBLIC ${BASE64_INCLUDE_DIR})
target_link_libraries(test_base64 PUBLIC ${BASE64_LIBRARIES})
add_test(NAME test_base64 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_base64)

# Test for Vc installation.
add_executable(test_vc test_vc.cpp)
target_include_directories(test_vc PUBLIC ${Vc_INCLUDE_DIR})
target_link_libraries(test_vc PUBLIC ${Vc_LIBRARIES})
add_test(NAME test_vc COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_vc)

# Test for PETSc installation.
add_executable(test_petsc test_petsc.cpp)
target_include_directories(test_petsc PUBLIC ${PETSC_INCLUDE_DIRS} ${MPI_CXX_INCLUDE_DIRS})
target_link_libraries(test_petsc PUBLIC ${PETSC_LINK_LIBRARIES} ${MPI_CXX_LIBRARIES})
add_test(NAME test_petsc COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_NUMPROC_MAX} ${CMAKE_CURRENT_BINARY_DIR}/test_petsc)

# Test for ADIOS2 installation.
#add_executable(test_adios2 test_adios2.cpp)
#target_link_libraries(test_adios2 ${adios2_LIBRARIES} ${MPI_CXX_LIBRARIES})
#add_test(NAME test_adios2 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_NUMPROC_MAX} ${CMAKE_CURRENT_BINARY_DIR}/test_adios2)

# Test for xbraid installation.
#add_executable(test_xbraid test_xbraid.cpp)
#target_link_libraries(test_xbraid ${XBRAID_LIBRARIES} ${MPI_CXX_LIBRARIES})
#add_test(NAME test_xbraid COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_NUMPROC_MAX} ${CMAKE_CURRENT_BINARY_DIR}/test_xbraid)

# Test for SEMT installation.
add_executable(test_semt test_semt.cpp)
target_include_directories(test_semt PUBLIC ${SEMT_INCLUDE_DIR})
target_link_libraries(test_semt PUBLIC ${SEMT_LIBRARIES})
add_test(NAME test_semt COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_semt)

# Test for googletest installation.
add_executable(test_gtest test_gtest.cpp)
target_include_directories(test_gtest PUBLIC ${GTEST_INCLUDE_DIRS})
target_link_libraries(test_gtest PUBLIC ${GTEST_BOTH_LIBRARIES})
add_test(NAME test_gtest COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_gtest)

# Test for preCICE installation.
add_executable(test_precice test_precice.cpp)
target_link_libraries(test_precice PUBLIC precice::precice)
add_test(NAME test_precice COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_precice)

# Test for LAPACK installation.
add_executable(test_lapacke test_lapacke.cpp)
target_include_directories(test_lapacke PUBLIC ${OpenBLAS_INCLUDE_DIRS})
target_link_libraries(test_lapacke PUBLIC ${OpenBLAS_LIBRARIES})
add_test(NAME test_lapacke COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_lapacke)

# Test for BLAS installation.
add_executable(test_blas test_blas.cpp)
target_include_directories(test_blas PUBLIC ${OpenBLAS_INCLUDE_DIRS})
target_link_libraries(test_blas PUBLIC ${OpenBLAS_LIBRARIES})
add_test(NAME test_blas COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_blas)
