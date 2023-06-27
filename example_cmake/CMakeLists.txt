
add_executable(laplace_linear ${PROJECT_SOURCE_DIR}/example_cmake/laplace_1d/laplace_linear.cpp)
target_include_directories(laplace_linear PUBLIC  ${PROJECT_SOURCE_DIR}/core/src
                                                  ${PROJECT_BINARY_DIR}/include
                                                  ${MPI_CXX_INCLUDE_DIRS}
                                                  ${ZLIB_INCLUDE_DIRS}
                                                  ${Python3_INCLUDE_DIRS}
                                                  ${LIBXML2_INCLUDE_DIR}
                                                  ${EASYLOGGINGPP_INCLUDE_DIRS}
                                                  ${BASE64_INCLUDE_DIR}
                                                  ${Vc_INCLUDE_DIR}
                                                  ${PETSC_INCLUDE_DIRS}
                                                  ${semt_DIR}
                                                  ${GTEST_INCLUDE_DIRS}
                                                  ${OpenBLAS_INCLUDE_DIRS})
target_link_libraries(laplace_linear PUBLIC ${PROJECT_BINARY_DIR}/core/libopendihu.so
                                            ${MPI_CXX_LIBRARIES}
                                            ${ZLIB_LIBRARIES}
                                            ${Python3_LIBRARIES}
                                            Eigen3::Eigen
                                            ${Boost_LIBRARIES}
                                            ${LIBXML2_LIBRARIES}
                                            ${EASYLOGGINGPP_LIBRARIES}
                                            ${BASE64_LIBRARIES}
                                            ${Vc_LIBRARIES}
                                            ${PETSC_LINK_LIBRARIES}
                                            ${SEMT_LIBRARIES}
                                            ${GTEST_BOTH_LIBRARIES}
                                            precice::precice
                                            ${OpenBLAS_LIBRARIES}
                                            ${OpenMP_CXX_LIBRARIES}
)