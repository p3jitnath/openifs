cmake_minimum_required(VERSION 3.22)
project(oifs_cuda_radiation LANGUAGES C CXX CUDA)

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES 90)
endif()

find_package(Threads REQUIRED)

add_library(oifs_cuda_radiation_loader SHARED oifs_cuda_radiation_loader.c)
target_link_libraries(oifs_cuda_radiation_loader PRIVATE ${CMAKE_DL_LIBS} Threads::Threads)

add_library(oifs_cuda_radiation SHARED oifs_cuda_radiation.cu)
target_include_directories(oifs_cuda_radiation PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(oifs_cuda_radiation PRIVATE cxx_std_17)
set_target_properties(oifs_cuda_radiation PROPERTIES
  CUDA_SEPARABLE_COMPILATION ON
  POSITION_INDEPENDENT_CODE ON)

target_compile_options(oifs_cuda_radiation PRIVATE
  $<$<COMPILE_LANGUAGE:CUDA>:--fmad=false --prec-div=true --prec-sqrt=true>)

add_executable(oifs_cuda_radiation_test test_cuda_radiation.cu)
target_link_libraries(oifs_cuda_radiation_test PRIVATE oifs_cuda_radiation)
add_executable(oifs_cuda_longwave_test test_cuda_longwave.cu)
target_link_libraries(oifs_cuda_longwave_test PRIVATE oifs_cuda_radiation)
enable_testing()
add_test(NAME oifs_cuda_radiation_test COMMAND oifs_cuda_radiation_test)
add_test(NAME oifs_cuda_longwave_test COMMAND oifs_cuda_longwave_test)
