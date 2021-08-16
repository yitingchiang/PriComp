cmake_minimum_required(VERSION 3.0)

project(ml_protocol)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../src/ ${GMP_INCLUDE_DIR})

file(GLOB SRC
     "*.cpp"
     "*.h")

add_library(ml_protocol SHARED ${SRC})

target_link_libraries(ml_protocol PUBLIC ${SMC_TWO} ${PRICOMP} ${GMP_LIBRARIES})
