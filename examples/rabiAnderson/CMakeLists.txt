cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project(rgflowKondo LANGUAGES CXX)

##########Dependencies########
find_package(HDF5 COMPONENTS CXX HL)
include_directories(${HDF5_INCLUDE_DIR})
set(LINK_LIBS ${LINK_LIBS} ${HDF5_LIBRARIES})
#lapack
########################################################################

################
file(GLOB sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
file(GLOB headers CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
######################################################

add_executable(${PROJECT_NAME}.out ${sources} ${headers})



target_link_libraries(
  ${PROJECT_NAME}.out NRG
  #openblas
  #cblas
  #lapack
  hdf5
  hdf5_cpp
  MKL::MKL
  mkl_intel_lp64 mkl_intel_thread mkl_core
)


