cmake_minimum_required(VERSION 3.8)
project(
    NRGCore
  VERSION 0.0.1
  LANGUAGES CXX
)
# set standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CLANG_ENABLE_STATIC_ANALYZER ON)
set(LLVM_ENABLE_PROJECTS "clang;clang-tools-extra;compiler-rt;llvm-spirv;sycl;libclc")
set(LLVM_ENABLE_RUNTIMES "libcxx;libcxxabi")

# Add the cmake folder so the FindSphinx module is found
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# add sources and headers
file(GLOB_RECURSE headers CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/include/*/*.hpp")
file(GLOB_RECURSE sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.cpp")
#



include_directories(SYSTEM "${CMAKE_CURRENT_SOURCE_DIR}/include")
#Add dependencies
find_package(Threads REQUIRED)  # Threads may or may not be required. For the profiler, the target that you link against
set(THREADS_PREFER_PTHREAD_FLAG ON) # This is just for us, we do prefer pthreads on our system, because linux etc. We
#MKL
#  The following 2 lines are very important if you wish to use the user instrumentation API detailed here:
#  https://software.intel.com/content/www/us/en/develop/documentation/vtune-help/top/api-support/instrumentation-and-tracing-technology-apis/instrumentation-and-tracing-technology-api-reference.html
#  This is so you can include and then link against their libs. For very basic example usage, see the thread_profiling_example.cpp
# ---- Create library ----

##########Dependencies########
find_package(HDF5 COMPONENTS CXX HL)
include_directories(${HDF5_INCLUDE_DIR})
set(LINK_LIBS ${LINK_LIBS} ${HDF5_LIBRARIES})

# add_dependencies(NRG ${headers} ${sources})

#----------------------------------------------------------------
add_library(NRG "${headers};${sources}" )
# add_library(NRG ${headers} ${sources})
# target_include_directories(NRG PUBLIC .)
# fix rebuild if header or src file changes 
target_include_directories(NRG  PRIVATE
     $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
     $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
)
#

target_include_directories(NRG INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include" )
target_include_directories(NRG PRIVATE /opt/intel/oneapi/vtune/latest/sdk/include)
target_link_directories(NRG PRIVATE /opt/intel/oneapi/vtune/latest/sdk/lib64)
