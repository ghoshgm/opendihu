cmake_minimum_required(VERSION 3.13...3.24)

# package versioning
include(cmake/git.cmake)

# Print package version to console.
message(STATUS "OpenDiHu ${PROJECT_VERSION}")

# Set CXX standard
set(CMAKE_CXX_STANDARD 17)

# Project description
project(OPENDIHU
LANGUAGES C CXX
HOMEPAGE_URL https://opendihu.readthedocs.io/en/latest/introduction.html#overview
DESCRIPTION "OpenDiHu is a scalable Finite Element code for skeletal muscle simulations."
VERSION ${PROJECT_VERSION}
)

message(STATUS "")
message(STATUS "****************************************************")
message(STATUS "           Checking for Modern CXX support          ")
message(STATUS "****************************************************")

# Chech CXX compiler support for modern CXX standards (11 or higher)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  if(CMAKE_CXX_COMPILER_VERSION LESS "5.0")
    message(FATAL_ERROR "GNU gcc v5.0 or higher required to enable modern cxx features.")
  else()
    message(STATUS "Detecting CXX compiler support for modern cxx standards - done")
    message(STATUS "Detected CXX standards:")
    foreach(feature ${CMAKE_CXX_COMPILE_FEATURES})
      if(feature MATCHES "^(cxx_std)")
        message(STATUS "${feature}")
      endif()
    endforeach()
    message(STATUS "Default CXX standard - ${CMAKE_CXX_STANDARD}")
  endif()
else()
  message(WARNING "Currently OpenDiHu only supports the GNU C/CXX compiler.CLang and Intel compilers will be supported in future releases.")
endif()

# Find the packages of all dependencies.
# The configure is aborted if a package cannot be found.
include(cmake/find_dependencies.cmake)

# Configure opendihu
include(cmake/config.cmake)

# Unit testing for dependecies and source code.
include(CTest)
add_subdirectory(testing/link_test)
