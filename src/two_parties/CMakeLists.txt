cmake_minimum_required(VERSION 3.0)

project(smc_two)

file(GLOB SRCS
     "*.h"
     "*.cpp"
     )

file(GLOB PROT_SRC
     "protocol/*.cpp"
     "protocol/*.h")
include_directories(${GMP_INCLUDE_DIR} ./protocol/)

add_library(smc_two SHARED ${SRCS} tables.cpp client.cpp ${PROT_SRC})
