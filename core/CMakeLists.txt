
# Gather all header files.
file(GLOB_RECURSE OPENDIHU_HEADERS src *.h *.tpp)

# Gather all source files.
file(GLOB_RECURSE OPENDIHU_SOURCES src *.cpp)

# Build libopendihu
add_library(opendihu SHARED ${OPENDIHU_HEADERS} ${OPENDIHU_SOURCES})
target_include_directories(opendihu PUBLIC ${PROJECT_SOURCE_DIR}/core/src
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
                                          ${OpenBLAS_INCLUDE_DIRS}
)
target_link_libraries(opendihu PUBLIC ${MPI_CXX_LIBRARIES}
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