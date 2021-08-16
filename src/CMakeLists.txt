cmake_minimum_required(VERSION 3.0)

project(pricomp)

include_directories(${GMP_INCLUDE_DIR})
add_subdirectory(two_parties)
file(GLOB SRC
     "*.cpp"
     "*.h")
     
add_library(pricomp SHARED ${SRC})
