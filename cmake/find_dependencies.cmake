
# Search for MPI installation.
find_package(MPI REQUIRED)

# Search for Vc installation.
find_package(Vc REQUIRED)
if(${Vc_FOUND})
 message(STATUS "Found Vc: ${Vc_DIR} (found version ${Vc_VERSION})")
endif()

# Search for zlib installation.
find_package(ZLIB REQUIRED)

# Search for python3 installation.
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

# Search for Eigen installation.
find_package(Eigen3 REQUIRED NO_MODULE)
if(${Eigen3_FOUND})
  message(STATUS "Found Eigen3: ${Eigen3_DIR} (found version ${Eigen3_VERSION})")
endif()

# Search for Boost installation.
find_package(Boost REQUIRED COMPONENTS log log_setup thread system filesystem program_options unit_test_framework)

# Search for xml2 installation.
find_package(LibXml2 REQUIRED)

# Search for preCICE installation.
find_package(precice REQUIRED)
if(${precice_FOUND})
 message(STATUS "Found preCICE: ${precice_DIR} (found version ${precice_VERSION})")
endif()

# Search for googletest installation.
find_package(GTest REQUIRED)

# Search for ADIOS2 installation.
find_package(adios2 REQUIRED)

# Search for OpenBLAS installation.
# It comes with LAPACK, LAPACKE, BLAS packed together.
find_package(OpenBLAS REQUIRED)
if(${OpenBLAS_FOUND})
  message(STATUS "Found OpenBLAS: ${OpenBLAS_DIR} (found version ${OpenBLAS_VERSION})")
endif()

# Search for Easylogging++ installation.
pkg_check_modules(EASYLOGGINGPP REQUIRED easyloggingpp)
if(${EASYLOGGINGPP_FOUND})
 message(STATUS "Found easyloggingpp at ${EASYLOGGINGPP_PREFIX}")
 add_library(elpp SHARED ${EASYLOGGINGPP_INCLUDE_DIRS}/easylogging++.h ${EASYLOGGINGPP_INCLUDE_DIRS}/easylogging++.cc)
 set(EASYLOGGINGPP_LIBRARIES elpp)
endif()

# Search for PETSc installation.
pkg_check_modules(PETSC REQUIRED PETSc)
if(${PETSC_FOUND})
 message(STATUS "Found PETSc at ${PETSC_PREFIX}")
endif()

# Search for xbraid installation.
# The library does not generate *.pc or *-config.cmake files.
# The includes and libs are pointed directly.
find_path(XBRAID_INCLUDE_DIR NAMES braid.h HINTS ${xbraid_DIR}/include)
find_path(XBRAID_LIBRARIES NAMES libbraid.a HINTS ${xbraid_DIR}/lib)
if(EXISTS ${XBRAID_INCLUDE_DIR} AND EXISTS ${XBRAID_LIBRARIES})
  set(XBRAID_FOUND 1)
  set(XBRAID_LIBRARIES "${xbraid_DIR}/lib/libbraid.a")
  message(STATUS "Found XBraid: ${xbraid_DIR}")
else()
  message(FATAL_ERROR "CMake could not find braid.h and/or static library in ${xbraid_DIR}.
                       Please make sure the path to the spack installation directory is correct.")
endif()

# Search for Base64 installation.
# The library has no support for CMake or Autotools and it does not generate any static/shared libs.
# A shared library is generated manually.
find_path(BASE64_INCLUDE_DIR NAMES base64.h HINTS ${base64_DIR})
if(EXISTS ${BASE64_INCLUDE_DIR})
  message(STATUS "Found base64: ${base64_DIR}")
  add_library(base64 SHARED ${BASE64_INCLUDE_DIR}/base64.h)
  set_target_properties(base64 PROPERTIES LINKER_LANGUAGE CXX)
  set(BASE64_LIBRARIES "base64")
else()
  message(FATAL_ERROR "CMake could not find base64.h in ${base64_DIR}.
                       Please make sure the path to the source code is correct.")
endif()

# Search for SEMT installation.
# The library has no support for CMake or Autotools and it does not generate any static/shared libs.
# A shared library is generated manually.
find_path(SEMT_INCLUDE_DIR NAMES Semt.h HINTS ${semt_DIR}/semt)
if(EXISTS ${SEMT_INCLUDE_DIR})
  message(STATUS "Found SEMT: ${semt_DIR}")
  add_library(semt SHARED
            ${semt_DIR}/semt/BinaryOperators.h
            ${semt_DIR}/semt/BinaryTypes.h
            ${semt_DIR}/semt/Common.h
            ${semt_DIR}/semt/Conditions.h
            ${semt_DIR}/semt/Constexpr.h
            ${semt_DIR}/semt/DifferentiableVectorExpr.h
            ${semt_DIR}/semt/Expression.h
            ${semt_DIR}/semt/Forwards.h
            ${semt_DIR}/semt/Macros.h
            ${semt_DIR}/semt/Parameter.h
            ${semt_DIR}/semt/Semtfwd.h
            ${semt_DIR}/semt/Semt.h
            ${semt_DIR}/semt/Sequence.h
            ${semt_DIR}/semt/Shortcuts.h
            ${semt_DIR}/semt/Simplifications.h
            ${semt_DIR}/semt/Traits.h
            ${semt_DIR}/semt/UnaryOperator.h
            ${semt_DIR}/semt/UnaryTypes.h
            ${semt_DIR}/semt/Variable.h
            ${semt_DIR}/semt/VectorExpr.cpp
            ${semt_DIR}/semt/VectorExpr.h
            ${semt_DIR}/loki/EmptyType.h
            ${semt_DIR}/loki/NullType.h
            ${semt_DIR}/loki/static_check.h
            ${semt_DIR}/loki/Typelist.h
            ${semt_DIR}/loki/TypelistMacros.h
            ${semt_DIR}/loki/TypeManip.h
          )
  set(SEMT_LIBRARIES "semt")
else()
  message(FATAL_ERROR "CMake could not find Semt.h in ${semt_DIR}.
                       Please make sure the path to the source code is correct.")
endif()

include_directories(${MPI_CXX_INCLUDE_DIRS}
                    ${Python3_INCLUDE_DIRS}
                    ${LIBXML2_INCLUDE_DIR}
                    ${ZLIB_INCLUDE_DIRS}
                    ${Vc_INCLUDE_DIR}
                    ${PETSC_INCLUDE_DIRS}
                    ${EASYLOGGINGPP_INCLUDE_DIRS}
                    ${BASE64_INCLUDE_DIR}
                    ${XBRAID_INCLUDE_DIR}
		                ${SEMT_INCLUDE_DIR}
		                ${GTEST_INCLUDE_DIRS}
                    ${OpenBLAS_INCLUDE_DIRS}
                   )
