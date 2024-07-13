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
set(CMAKE_CXX_FLAGS " -D_GLIBCXX_USE_TBB_PAR_BACKEND=0 -fsanitize=address,leak,undefined -Wno-narrowing -fsanitize-recover=all -g -O0" )
# Add the cmake folder so the FindSphinx module is found
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

#Add dependencies
find_package(Threads REQUIRED)  # Threads may or may not be required. For the profiler, the target that you link against
set(THREADS_PREFER_PTHREAD_FLAG ON) # This is just for us, we do prefer pthreads on our system, because linux etc. We
#MKL
#  The following 2 lines are very important if you wish to use the user instrumentation API detailed here:
#  https://software.intel.com/content/www/us/en/develop/documentation/vtune-help/top/api-support/instrumentation-and-tracing-technology-apis/instrumentation-and-tracing-technology-api-reference.html
#  This is so you can include and then link against their libs. For very basic example usage, see the thread_profiling_example.cpp
# ---- Create library ----

#----------------------------------------------------------------
# Create Documentation 
add_subdirectory("nrgcore")
add_subdirectory("docs")
# add_subdirectory("examples/rgflowSIAM")
#add_subdirectory("examples/entropySIAM")
# add all sub directories .. 
file(GLOB V_GLOB LIST_DIRECTORIES true "examples/*")
foreach(item ${V_GLOB})
	if(IS_DIRECTORY ${item})
		add_subdirectory(${item})
	endif()
endforeach()
