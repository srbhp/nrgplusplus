cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project(twoChannelSiam LANGUAGES CXX)


##########Dependencies########
find_package(HDF5 COMPONENTS CXX HL)
include_directories(${HDF5_INCLUDE_DIR})
set(LINK_LIBS ${LINK_LIBS} ${HDF5_LIBRARIES})
#lapack
find_package(LAPACK)
find_package(BLAS)
set(BLA_VENDOR Intel10_64lp)
################
file(GLOB sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
file(GLOB headers CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(${PROJECT_NAME}.output ${sources} ${headers})
target_link_libraries(
  ${PROJECT_NAME}.output NRG
  hdf5
  hdf5_cpp
	mkl_intel_lp64 mkl_core mkl_intel_thread iomp5

	)
