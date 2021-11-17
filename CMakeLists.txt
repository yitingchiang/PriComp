cmake_minimum_required(VERSION 2.8)
project(PriComp_Demo)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/modules/)
find_package(GMP REQUIRED)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib/)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/)
link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
#find_library(libsmc_two NAMES smc_two PATHS ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
#find_library(libpricomp NAMES pricomp PATHS ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -mcmodel=large -m64 -lm")

if(EVAL_TIME STREQUAL "Y")
  message("PriComp programs will generate the time and the number of scalar product protocols performed.\n")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D EVAL_TIME")
endif()

if(DEBUG STREQUAL "Y")
  message("PriComp programs will generate debug messages.\n")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D DEBUG")
endif()

include_directories(./src/ ./src/two_parties/ ./src/two_parties/protocol/ ${GMP_INCLUDE_DIR})

add_subdirectory(src)

add_executable(stat_demo stat_demo.cpp)
add_executable(arith_demo arith_demo.cpp)
add_executable(sp_demo sp_demo.cpp)
add_executable(f_cond_demo f_condition_demo.cpp)
add_executable(millionair millionair.cpp)

target_link_libraries(stat_demo smc_two pricomp ${GMP_LIBRARIES})
target_link_libraries(arith_demo smc_two pricomp ${GMP_LIBRARIES})
target_link_libraries(sp_demo smc_two pricomp ${GMP_LIBRARIES})
target_link_libraries(f_cond_demo smc_two pricomp ${GMP_LIBRARIES})
target_link_libraries(millionair smc_two pricomp ${GMP_LIBRARIES})
